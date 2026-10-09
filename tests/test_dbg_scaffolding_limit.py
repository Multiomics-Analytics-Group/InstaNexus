#!/usr/bin/env python3

"""dbg scaffolding on noisy input must stop with a clear error instead of running indefinitely (issue #58)."""

import logging
import random
import subprocess
import sys
from pathlib import Path

import pandas as pd
import pytest

from instanexus import assembly as asm

# Exit code documented for ScaffoldingLimitExceeded; hard-coded so this test also runs (and fails) without the fix
EXIT_SCAFFOLDING_LIMIT = 3


def noisy_peptides(n: int = 40, seed: int = 0) -> list:
    """Random peptides that start and end with one of a few shared tripeptides.

    Every contig then overlaps many others by 3 residues, like contigs from low-confidence de novo
    peptides, so the number of overlaps explodes from one scaffolding round to the next.

    Args:
        n: Number of peptides.
        seed: Random seed.

    Returns:
        Peptide sequences.
    """
    rng = random.Random(seed)
    amino_acids = "ACDEFGHIKLMNPQRSTVWY"
    motifs = ["GPK", "LSR", "DEK"]

    return [
        rng.choice(motifs) + "".join(rng.choice(amino_acids) for _ in range(rng.randint(6, 9))) + rng.choice(motifs)
        for _ in range(n)
    ]


def clean_peptides(n: int = 150, seed: int = 7) -> list:
    """Overlapping fragments of one protein, as in tests/test_assembly_determinism.py.

    Args:
        n: Number of peptides.
        seed: Random seed.

    Returns:
        Peptide sequences.
    """
    protein = (
        "MKWVTFISLLLLFSSAYSRGVFRRDTHKSEIAHRFKDLGEEHFKGLVLIAFSQYLQQCPFDEHVKLVNELTEFAKTCVADESHAGCEKSLHTLFGDELCKVASL"
        "RETYGDMADCCEKQEPERNECFLSHKDDSPDLPKLKPDPNTLCDEFKADEKKFWGKYLYEIARRHPYFYAPELLYYANKYNGVFQECCQAEDKGACLLPKIET"
    )
    rng = random.Random(seed)
    peptides = []
    for _ in range(n):
        length = rng.randint(7, 15)
        start = rng.randrange(0, len(protein) - length)
        peptides.append(protein[start : start + length])

    return peptides


def test_dbg_cli_stops_with_dedicated_exit_code_on_noisy_input(tmp_path: Path) -> None:
    input_csv = tmp_path / "noisy.csv"
    pd.DataFrame({"cleaned_preds": noisy_peptides()}).to_csv(input_csv, index=False)

    # Without the limit this input does not finish within minutes
    proc = subprocess.run(
        [
            sys.executable,
            "-m",
            "instanexus.assembly",
            "--input-csv-path",
            str(input_csv),
            "--output-scaffolds-path",
            str(tmp_path / "scaffolds.fasta"),
            "--assembly-mode",
            "dbg",
            "--kmer-size",
            "7",
            "--min-overlap",
            "3",
            "--size-threshold",
            "10",
        ],
        capture_output=True,
        text=True,
        timeout=120,
    )

    assert proc.returncode == EXIT_SCAFFOLDING_LIMIT, proc.stderr[-2000:]
    assert "--max-scaffold-overlaps" in proc.stderr
    assert not (tmp_path / "scaffolds.fasta").exists()


def test_limit_raises_with_counts() -> None:
    assembler = asm.Assembler(mode="dbg", kmer_size=7, min_overlap=3, size_threshold=10, max_scaffold_overlaps=1_000)

    with pytest.raises(asm.ScaffoldingLimitExceeded) as excinfo:
        assembler.run(noisy_peptides())

    assert excinfo.value.limit == 1_000
    assert excinfo.value.n_overlaps > 1_000
    assert excinfo.value.min_overlap == 3


def test_limit_does_not_change_result_on_clean_input() -> None:
    peptides = clean_peptides()
    params = dict(mode="dbg", kmer_size=7, min_overlap=3, size_threshold=10)

    limited = asm.Assembler(**params).run(peptides)
    unlimited = asm.Assembler(**params, max_scaffold_overlaps=None).run(peptides)

    assert limited
    assert limited == unlimited


def _contained_reference(seqs: list) -> list:
    """Previous all-pairs implementation of merge_sequences_dbg."""
    merged = set(seqs)
    for c in seqs:
        for c2 in seqs:
            if c != c2 and c2 in c:
                merged.discard(c2)

    return asm.sort_by_length(merged)


def test_indexed_search_matches_all_pairs_search() -> None:
    rng = random.Random(0)
    for _ in range(500):
        alphabet = rng.choice(["AB", "ACDE", "ACDEFGHIKLMNPQRSTVWY"])
        seqs = ["".join(rng.choice(alphabet) for _ in range(rng.randint(0, 12))) for _ in range(rng.randint(0, 20))]
        seqs += seqs[: rng.randint(0, 3)]
        min_overlap = rng.randint(0, 5)

        assert asm.find_sliding_overlaps_indexed(seqs, min_overlap) == asm.find_sliding_overlaps(seqs, min_overlap)
        assert asm.merge_sequences_dbg(seqs, disable_tqdm=True) == _contained_reference(seqs)


def test_near_miss_overlaps_are_logged_at_debug_level_only(caplog: pytest.LogCaptureFixture) -> None:
    # suffix "ABCDE" of the first sequence differs from prefix "ABCDF" of the second in one residue
    seqs = ["XXABCDE", "ABCDFYY"]

    with caplog.at_level(logging.INFO, logger=asm.logger.name):
        asm.find_sliding_overlaps(seqs, 3)
    assert "POTENTIAL OVERLAP MISSED" not in caplog.text

    with caplog.at_level(logging.DEBUG, logger=asm.logger.name):
        asm.find_sliding_overlaps(seqs, 3)
    assert "POTENTIAL OVERLAP MISSED: ABCDE vs ABCDF" in caplog.text
