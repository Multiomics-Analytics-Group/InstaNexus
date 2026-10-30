#!/usr/bin/env python3

"""Mass-aware mapping to the reference (issue #70): scaffolds with an indel, an isobaric block or a swap map,
and the ungapped default is unchanged."""

import json
import subprocess
import sys
from pathlib import Path

import pandas as pd
import pytest

from instanexus import assembly, helpers, mass_alignment
from instanexus import visualization as viz

# normalized, so I is L
REFERENCE = "MKTAYLAKQRQLSFVKSHFNRQLEERLGLLEVQAPLLSRVGDGTQDNLSGAEK"
SEGMENT = REFERENCE[4:24]  # YLAKQRQLSFVKSHFNRQLE


def mass_aware(seq: str, max_mismatches: int = 10, min_identity: float = 0.8) -> tuple:
    """Map seq to REFERENCE with the mass-aware alignment.

    Args:
        seq: Sequence to map.
        max_mismatches: Maximum mismatches allowed.
        min_identity: Minimum identity required.

    Returns:
        The mapping, or None.
    """
    return viz.map_to_protein(seq, REFERENCE, max_mismatches, min_identity, alignment="mass-aware")


def ungapped(seq: str, max_mismatches: int = 10, min_identity: float = 0.8) -> tuple:
    """Map seq to REFERENCE with the default ungapped comparison.

    Args:
        seq: Sequence to map.
        max_mismatches: Maximum mismatches allowed.
        min_identity: Minimum identity required.

    Returns:
        The mapping, or None.
    """
    return viz.map_to_protein(seq, REFERENCE, max_mismatches, min_identity)


def differences(mapping: tuple) -> dict:
    """Return the non-zero difference counts of a mass-aware mapping."""
    return {kind: count for kind, count in mapping.alignment.differences.items() if count}


@pytest.mark.parametrize(
    "seq",
    [
        pytest.param(SEGMENT, id="exact"),
        pytest.param(SEGMENT[:6] + "W" + SEGMENT[7:], id="one-substitution"),
        pytest.param("D" + SEGMENT[1:12] + "W" + SEGMENT[13:], id="two-substitutions"),
    ],
)
def test_substitutions_map_as_they_do_ungapped(seq: str) -> None:
    assert tuple(mass_aware(seq)) == ungapped(seq)


def test_one_deleted_residue() -> None:
    seq = SEGMENT[:10] + SEGMENT[11:]

    assert ungapped(seq) is None
    mapping = mass_aware(seq)
    start, end, mismatches, identity = mapping
    # the reference region keeps the residue the scaffold lacks
    assert (start, end) == (4, 24)
    assert mismatches == [10]
    assert identity == pytest.approx(19 / 20)
    assert differences(mapping) == {"deletions": 1}


def test_one_inserted_residue() -> None:
    seq = SEGMENT[:10] + "W" + SEGMENT[10:]

    assert ungapped(seq) is None
    mapping = mass_aware(seq)
    start, end, mismatches, identity = mapping
    assert (start, end) == (4, 24)
    # an inserted residue counts at the reference residue after it
    assert mismatches == [10]
    assert identity == pytest.approx(20 / 21)
    assert differences(mapping) == {"insertions": 1}
    assert mapping.alignment.residue_positions[10] is None
    assert mapping.alignment.residue_positions[11] == 4 + 10


@pytest.mark.parametrize(
    ("residue", "block"),
    [pytest.param("N", "GG", id="N-as-GG"), pytest.param("Q", "AG", id="Q-as-AG")],
)
def test_isobaric_block_read_for_one_residue(residue: str, block: str) -> None:
    offset = SEGMENT.index(residue)
    seq = SEGMENT[:offset] + block + SEGMENT[offset + 1 :]

    assert ungapped(seq) is None
    mapping = mass_aware(seq)
    start, end, mismatches, identity = mapping
    assert (start, end) == (4, 24)
    # two columns, both at the reference residue
    assert mismatches == [offset, offset]
    assert identity == pytest.approx(19 / 21)
    assert differences(mapping) == {"isobaric_blocks": 1}


def test_isobaric_block_in_the_reference() -> None:
    reference = REFERENCE[:30] + "GG" + REFERENCE[30:]
    seq = reference[22:30] + "N" + reference[32:42]

    mapping = viz.map_to_protein(seq, reference, 10, 0.8, alignment="mass-aware")

    assert (mapping[0], mapping[1]) == (22, 42)
    assert mapping[2] == [8, 9]
    assert differences(mapping) == {"isobaric_blocks": 1}


def test_r_and_gv_differ_by_more_than_the_tolerance() -> None:
    # R is 156.1011 Da and GV 156.0899 Da: 0.0112 Da apart
    offset = SEGMENT.index("R")
    seq = SEGMENT[:offset] + "GV" + SEGMENT[offset + 1 :]

    assert mass_alignment.align(seq, REFERENCE).differences["isobaric_blocks"] == 0
    assert mass_alignment.align(seq, REFERENCE, isobaric_tolerance=0.02).differences["isobaric_blocks"] == 1


@pytest.mark.parametrize("size", [2, 3, 4])
def test_shuffled_block(size: int) -> None:
    # a swap at size 2; a permutation has the mass of the block, but is a shuffle, not an isobaric block
    offset = SEGMENT.index("SFVK")
    block = SEGMENT[offset : offset + size]
    seq = SEGMENT[:offset] + block[1:] + block[0] + SEGMENT[offset + size :]

    mapping = mass_aware(seq)

    assert (mapping[0], mapping[1]) == (4, 24)
    assert mapping[2] == list(range(offset, offset + size))
    assert differences(mapping) == {"shuffled_blocks": 1}


def test_overhang_before_the_reference_start() -> None:
    seq = "WW" + REFERENCE[:18]

    mapping = mass_aware(seq)

    assert (mapping[0], mapping[1]) == (0, 18)
    assert mapping[2] == [0, 0]
    assert differences(mapping) == {"insertions": 2}


def test_thresholds_apply_to_the_alignment() -> None:
    seq = SEGMENT[:5] + SEGMENT[6:10] + SEGMENT[11:]

    assert mass_aware(seq, max_mismatches=2, min_identity=0.9) is not None
    assert mass_aware(seq, max_mismatches=1) is None
    assert mass_aware(seq, min_identity=0.95) is None


def test_unrelated_and_empty_sequences_do_not_map() -> None:
    assert mass_aware("WWWWWWWWWW") is None
    assert mass_alignment.align("", REFERENCE) is None


def test_unknown_alignment_raises() -> None:
    with pytest.raises(ValueError, match="alignment must be one of"):
        viz.map_to_protein(SEGMENT, REFERENCE, 10, 0.8, alignment="blosum62")


def test_statistics_count_each_kind_of_difference(tmp_path: Path) -> None:
    deletion = SEGMENT[:10] + SEGMENT[11:]
    isobaric = REFERENCE[30:46] + "GG" + REFERENCE[47:]
    mapped = viz.process_protein_contigs_scaffold([deletion, isobaric], REFERENCE, 10, 0.8, alignment="mass-aware")
    df = viz.create_dataframe_from_mapped_sequences(mapped)

    statistics = helpers.compute_assembly_statistics(df, "scaffolds", str(tmp_path), REFERENCE)

    assert statistics["coverage"] == pytest.approx((20 + len(REFERENCE) - 30) / len(REFERENCE))
    # the reference regions: 20 residues, and the 23 from 30 to the end
    assert (statistics["min_length"], statistics["max_length"]) == (20, len(REFERENCE) - 30)
    assert statistics["total_mismatches"] == 3
    assert statistics["mismatched_positions"] == 2
    assert statistics["perfect_matches"] == 0
    assert {kind: statistics[kind] for kind in mass_alignment.DIFFERENCE_KINDS} == {
        "substitutions": 0,
        "insertions": 0,
        "deletions": 1,
        "isobaric_blocks": 1,
        "shuffled_blocks": 0,
    }


def test_ungapped_statistics_have_no_difference_kinds(tmp_path: Path) -> None:
    df = viz.create_dataframe_from_mapped_sequences(viz.process_protein_contigs_scaffold([SEGMENT], REFERENCE, 10, 0.8))

    statistics = helpers.compute_assembly_statistics(df, "scaffolds", str(tmp_path), REFERENCE)

    assert not set(mass_alignment.DIFFERENCE_KINDS) & set(statistics)


def test_isoleucine_statistics_follow_the_alignment_past_an_insertion() -> None:
    # isoleucine at reference position 11 (normalized: the L of SEGMENT[7])
    reference = REFERENCE[:11] + "I" + REFERENCE[12:]
    normalized = SEGMENT[:3] + "W" + SEGMENT[3:]
    reported = normalized[:8] + "I" + normalized[9:]
    mapped = viz.process_protein_contigs_scaffold([normalized], REFERENCE, 10, 0.8, alignment="mass-aware")

    il = helpers.compute_isoleucine_statistics(mapped, {normalized: reported}, reference)

    # placed at start + offset, every I/L past the insertion would sit one position off
    assert (il["il_positions_covered"], il["il_correct"]) == (3, 3)


@pytest.mark.parametrize("module", ["instanexus.main", "instanexus.assembly"])
def test_cli_offers_the_option(module: str) -> None:
    proc = subprocess.run([sys.executable, "-m", module, "--help"], capture_output=True, text=True, timeout=120)

    assert "--reference-alignment {ungapped,mass-aware}" in proc.stdout


def test_assembly_records_the_alignment_in_its_statistics(tmp_path: Path) -> None:
    metadata = tmp_path / "metadata.json"
    metadata.write_text(json.dumps({"sample": [{"protein": REFERENCE, "proteases": []}]}))
    # overlapping fragments with one residue deleted near the end: ungapped, the scaffold maps with the residues
    # after the deletion as mismatches, and does not reach the end of the reference
    with_deletion = REFERENCE[:46] + REFERENCE[47:]
    fragments = [with_deletion[start : start + 12] for start in range(0, len(with_deletion) - 11, 3)]
    fragments.append(with_deletion[-12:])
    input_csv = tmp_path / "sample_cleaned.csv"
    pd.DataFrame({"cleaned_preds": fragments, "read_residues": fragments}).to_csv(input_csv, index=False)

    for alignment in viz.REFERENCE_ALIGNMENTS:
        output_fasta = tmp_path / alignment / "scaffolds.fasta"
        assembly.main(
            input_csv_path=str(input_csv),
            output_scaffolds_path=str(output_fasta),
            metadata_json_path=str(metadata),
            assembly_mode="greedy",
            kmer_size=7,
            min_overlap=3,
            size_threshold=10,
            reference=True,
            chain="",
            min_identity=0.8,
            max_mismatches=10,
            reference_alignment=alignment,
        )

    ungapped_stats = json.loads((tmp_path / "ungapped" / "statistics" / "scaffolds_stats.json").read_text())
    aligned_stats = json.loads((tmp_path / "mass-aware" / "statistics" / "scaffolds_stats.json").read_text())
    assert ungapped_stats["reference_alignment"] == "ungapped"
    assert aligned_stats["reference_alignment"] == "mass-aware"
    assert aligned_stats["deletions"] == 1
    assert aligned_stats["coverage"] == 1.0
    assert ungapped_stats["coverage"] < 1.0
    assert aligned_stats["total_mismatches"] < ungapped_stats["total_mismatches"]


@pytest.mark.parametrize(
    ("scaffold", "region", "kinds"),
    [
        pytest.param(
            "SNFMLTQPRSVSESPGKTVTLSATRSSGSLGSDYVHWYQQRPGSSPTML",
            (0, 49),
            ["insertion", "substitution", "deletion", "substitution"],
            id="overhang-and-deletion",
        ),
        pytest.param(
            "GEKSSPVKAGVETTTPSKQSNNKAYASS", (155, 182), ["isobaric", "insertion", "shuffle"], id="GE-as-AD-and-a-swap"
        ),
    ],
)
def test_published_ma2_light_chain_scaffolds(scaffold: str, region: tuple, kinds: list) -> None:
    # scaffolds of the published ma2 greedy_light_fdr0.1 assembly (Zenodo 10.5281/zenodo.18414969)
    metadata = json.loads((Path(__file__).parents[1] / "json" / "sample_metadata.json").read_text())
    light_chain = next(entry for entry in metadata["ma2"] if entry["chain"] == "light")["protein"].replace("I", "L")

    mapping = viz.map_to_protein(scaffold, light_chain, 10, 0.8, alignment="mass-aware")

    assert (mapping[0], mapping[1]) == region
    assert [op.kind for op in mapping.alignment.operations if op.kind != "match"] == kinds
