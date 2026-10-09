"""Isoleucine is reported in scaffolds, having been normalized away for assembly."""

import argparse
import logging
import subprocess
import sys
from pathlib import Path

import pandas as pd
from Bio import SeqIO

from instanexus import assembly
from instanexus.assembly import count_isoleucine_leucine, restore_isoleucine
from instanexus.helpers import compute_isoleucine_statistics
from instanexus.preprocessing import normalize_sequence, remove_modifications, strip_modifications


class TestStripModifications:
    """`strip_modifications` is `remove_modifications` without the I/L normalization."""

    def test_keeps_isoleucine_where_remove_modifications_collapses_it(self):
        assert remove_modifications("ALIVTQTMK") == "ALLVTQTMK"
        assert strip_modifications("ALIVTQTMK") == "ALIVTQTMK"

    def test_both_strip_the_same_modifications(self):
        for sequence in ("ALKALPM[UNIMOD:35]HIR", "AL(ox)KALPMHIR", "-ALKALPMHIR-"):
            assert normalize_sequence(strip_modifications(sequence)) == remove_modifications(sequence)

    def test_none_survives(self):
        assert strip_modifications(None) is None


class TestRestoreIsoleucine:
    """The reads vote on each I/L position they cover."""

    def test_a_single_read_restores_its_own_isoleucine(self):
        restored, decisions = restore_isoleucine(["AAALVTQTMKAA"], ["ALVTQTMK"], ["AIVTQTMK"])
        assert restored == ["AAAIVTQTMKAA"]
        assert decisions[0]["position_1based"] == 4
        assert decisions[0]["votes_isoleucine"] == 1
        assert decisions[0]["call"] == "I"

    def test_the_majority_decides(self):
        # Three reads of one position, two reading isoleucine.
        restored, _ = restore_isoleucine(["ALVTQTMK"], ["ALVTQTMK"] * 3, ["AIVTQTMK", "AIVTQTMK", "ALVTQTMK"])
        assert restored == ["AIVTQTMK"]

        # The same with the majority the other way.
        restored, _ = restore_isoleucine(["ALVTQTMK"], ["ALVTQTMK"] * 3, ["AIVTQTMK", "ALVTQTMK", "ALVTQTMK"])
        assert restored == ["ALVTQTMK"]

    def test_a_tie_stays_leucine(self):
        # No majority is not a reason to change the residue: leucine is both the
        # status quo and the commoner of the two.
        restored, decisions = restore_isoleucine(["ALVTQTMK"], ["ALVTQTMK"] * 2, ["AIVTQTMK", "ALVTQTMK"])
        assert restored == ["ALVTQTMK"]
        assert [d["call"] for d in decisions] == ["tie"]

    def test_a_position_no_read_covers_stays_leucine(self):
        restored, _ = restore_isoleucine(["ALVTQTMKLLLL"], ["ALVTQTMK"], ["AIVTQTMK"])
        assert restored == ["AIVTQTMKLLLL"]

    def test_positions_are_voted_on_independently(self):
        # Two I/L positions, each decided by the reads covering it alone.
        restored, _ = restore_isoleucine(["ALVTQTMKLR"], ["ALVTQTMKLR"] * 2, ["AIVTQTMKLR", "AIVTQTMKIR"])
        assert restored == ["AIVTQTMKLR"]

    def test_a_read_occurring_twice_votes_at_both_places(self):
        restored, _ = restore_isoleucine(["ALKGGGALK"], ["ALK"], ["AIK"])
        assert restored == ["AIKGGGAIK"]

    def test_a_read_whose_residues_do_not_line_up_cannot_vote(self):
        # A wrong offset is worse than no vote, so a mismatched pair is skipped.
        restored, decisions = restore_isoleucine(["ALVTQTMK"], ["ALVTQTMK"], ["AIVTQTM"])
        assert restored == ["ALVTQTMK"]
        assert [d["call"] for d in decisions] == ["no_reads"]

    def test_scaffolds_with_no_matching_read_are_returned_unchanged(self):
        restored, decisions = restore_isoleucine(["ALVTQTMK"], ["GGGGG"], ["GGGGG"])
        assert restored == ["ALVTQTMK"]
        assert [d["call"] for d in decisions] == ["no_reads"]

    def test_every_il_position_is_recorded_with_its_call(self):
        # Four I/L positions: isoleucine wins, leucine wins, a tie, and one no read covers.
        restored, decisions = restore_isoleucine(
            ["GLGLGLGGL"],
            ["GLGLGLG"] * 4,
            ["GIGLGIG", "GIGLGLG", "GLGLGIG", "GIGIGLG"],
        )
        assert restored == ["GIGLGLGGL"]
        assert decisions == [
            {"scaffold": "scaffold_1", "position_1based": 2, "votes_isoleucine": 3, "votes_leucine": 1, "call": "I"},
            {"scaffold": "scaffold_1", "position_1based": 4, "votes_isoleucine": 1, "votes_leucine": 3, "call": "L"},
            {"scaffold": "scaffold_1", "position_1based": 6, "votes_isoleucine": 2, "votes_leucine": 2, "call": "tie"},
            {
                "scaffold": "scaffold_1",
                "position_1based": 9,
                "votes_isoleucine": 0,
                "votes_leucine": 0,
                "call": "no_reads",
            },
        ]

    def test_nothing_other_than_i_and_l_is_touched(self):
        scaffold = "ALVTQTMKGGWYF"
        restored, _ = restore_isoleucine([scaffold], [scaffold], ["AIVTQTMKGGWYF"])
        assert restored[0].replace("I", "L") == scaffold


class TestIsoleucineStatistics:
    """Placement is normalized; the residues are then compared as they were read."""

    @staticmethod
    def _mapped(normalized, start):
        # What `map_to_protein` returns: (sequence, (start, end, mismatches, identity))
        # with `start` a 0-based index into the reference and `end` exclusive.
        return [(normalized, (start, start + len(normalized), [], 1.0))]

    def test_counts_only_reference_il_positions(self):
        reference = "GGAILGG"  # I at index 3, L at index 4
        stats = compute_isoleucine_statistics(self._mapped("GGALLGG", 0), {"GGALLGG": "GGAILGG"}, reference)
        assert stats["il_positions_covered"] == 2
        assert stats["il_correct"] == 2
        assert stats["il_accuracy"] == 1.0

    def test_a_wrong_residue_is_counted_as_wrong(self):
        reference = "GGAILGG"
        # Reported LL where the reference reads IL: one right, one wrong.
        stats = compute_isoleucine_statistics(self._mapped("GGALLGG", 0), {"GGALLGG": "GGALLGG"}, reference)
        assert stats["il_positions_covered"] == 2
        assert stats["il_correct"] == 1
        assert stats["il_accuracy"] == 0.5

    def test_this_is_what_the_normalized_comparison_could_not_see(self):
        # Both of these score identically on a normalized comparison -- that is the
        # defect. They must not score identically here.
        reference = "GGAIIGG"
        right = compute_isoleucine_statistics(self._mapped("GGALLGG", 0), {"GGALLGG": "GGAIIGG"}, reference)
        wrong = compute_isoleucine_statistics(self._mapped("GGALLGG", 0), {"GGALLGG": "GGALLGG"}, reference)
        assert right["il_accuracy"] == 1.0
        assert wrong["il_accuracy"] == 0.0

    def test_a_scaffold_placed_partway_in_uses_the_offset(self):
        reference = "GGGGIL"
        stats = compute_isoleucine_statistics(self._mapped("LL", 4), {"LL": "IL"}, reference)
        assert stats["il_positions_covered"] == 2
        assert stats["il_correct"] == 2

    def test_positions_past_the_end_of_the_reference_are_skipped(self):
        stats = compute_isoleucine_statistics(self._mapped("ILGG", 2), {"ILGG": "ILGG"}, "GGIL")
        assert stats["il_positions_covered"] == 2
        assert stats["il_correct"] == 2

    def test_no_il_under_a_scaffold_gives_zeros_rather_than_dividing(self):
        stats = compute_isoleucine_statistics(self._mapped("GGGG", 0), {"GGGG": "GGGG"}, "GGGG")
        assert stats == {
            "il_positions_covered": 0,
            "il_correct": 0,
            "il_accuracy": 0.0,
            "il_accuracy_all_leucine": 0.0,
        }

    def test_the_all_leucine_baseline_is_reported_beside_the_accuracy(self):
        # Three reference I/L positions, two of them L: reporting everything as L
        # scores 2/3 for free, which is what the accuracy has to be read against.
        reference = "GGILLGG"
        stats = compute_isoleucine_statistics(self._mapped("GGLLLGG", 0), {"GGLLLGG": "GGILLGG"}, reference)
        assert stats["il_positions_covered"] == 3
        assert stats["il_accuracy"] == 1.0
        assert round(stats["il_accuracy_all_leucine"], 4) == round(2 / 3, 4)


# A stretch of bovine serum albumin with both isoleucine and leucine
PROTEIN = "MKWVTFISLLLLFSSAYSRGVFRRDTHKSEIAHRFKDLGEEHFKGLVLIAFSQYLQQ"


def _run_command(tmp_path: Path, read_residues: list, **kwargs) -> tuple:
    """Assemble overlapping fragments of PROTEIN through ``assembly.main``, as the command line does.

    Args:
        tmp_path: Folder for input and outputs.
        read_residues: Predicted residues of each fragment, one per fragment.
        **kwargs: Extra arguments for ``assembly.main``.

    Returns:
        The reported scaffolds and the path of isoleucine_restoration.tsv.
    """
    input_csv = tmp_path / "cleaned.csv"
    pd.DataFrame(
        {"cleaned_preds": [normalize_sequence(r) for r in read_residues], "read_residues": read_residues}
    ).to_csv(input_csv, index=False)
    output_fasta = tmp_path / "out" / "scaffolds.fasta"
    assembly.main(
        input_csv_path=str(input_csv),
        output_scaffolds_path=str(output_fasta),
        metadata_json_path=None,
        assembly_mode="greedy",
        kmer_size=7,
        min_overlap=3,
        size_threshold=10,
        reference=False,
        chain="",
        min_identity=0.8,
        max_mismatches=10,
        **kwargs,
    )
    scaffolds = [str(record.seq) for record in SeqIO.parse(output_fasta, "fasta")]

    return scaffolds, output_fasta.parent / "isoleucine_restoration.tsv"


FRAGMENTS = [PROTEIN[start : start + 12] for start in range(0, len(PROTEIN) - 11, 4)]


class TestRestorationInTheCommand:
    """``assembly.main`` votes only when the reads can tell I from L, and can be told not to."""

    def test_reads_using_both_letters_are_voted_on(self, tmp_path):
        scaffolds, tsv = _run_command(tmp_path, FRAGMENTS)

        assert any("I" in s for s in scaffolds)
        assert tsv.exists()
        assert set(pd.read_csv(tsv, sep="\t")["call"]) <= {"I", "L", "tie", "no_reads"}

    def test_reads_with_only_isoleucine_are_not_voted_on(self, tmp_path, caplog):
        # Some models write I for every I/L: a vote would turn every position into I.
        only_isoleucine = [fragment.replace("L", "I") for fragment in FRAGMENTS]

        with caplog.at_level(logging.WARNING, logger=assembly.logger.name):
            scaffolds, tsv = _run_command(tmp_path, only_isoleucine)

        assert "do not distinguish isoleucine from leucine" in caplog.text
        assert scaffolds
        assert not any("I" in s for s in scaffolds)
        assert not tsv.exists()

    def test_count_isoleucine_leucine(self):
        assert count_isoleucine_leucine(["AIL", "IIG", None]) == (3, 1)
        assert count_isoleucine_leucine(["AIK", "IIG"]) == (3, 0)

    def test_restoration_can_be_disabled(self, tmp_path):
        scaffolds, tsv = _run_command(tmp_path, FRAGMENTS, isoleucine_restoration=False)

        assert scaffolds
        assert not any("I" in s for s in scaffolds)
        assert not tsv.exists()

    def test_option_is_on_by_default_and_turned_off_by_the_flag(self):
        parser = argparse.ArgumentParser()
        assembly.add_isoleucine_restoration_argument(parser)

        assert parser.parse_args([]).isoleucine_restoration is True
        assert parser.parse_args(["--no-isoleucine-restoration"]).isoleucine_restoration is False

    def test_flag_reaches_the_assembly_command(self, tmp_path):
        input_csv = tmp_path / "cleaned.csv"
        pd.DataFrame({"cleaned_preds": [normalize_sequence(f) for f in FRAGMENTS], "read_residues": FRAGMENTS}).to_csv(
            input_csv, index=False
        )
        args = [
            sys.executable,
            "-m",
            "instanexus.assembly",
            "--input-csv-path",
            str(input_csv),
            "--size-threshold",
            "10",
        ]

        for flag, expect_isoleucine in (([], True), (["--no-isoleucine-restoration"], False)):
            output_fasta = tmp_path / ("off" if flag else "on") / "scaffolds.fasta"
            proc = subprocess.run(
                args + ["--output-scaffolds-path", str(output_fasta)] + flag,
                capture_output=True,
                text=True,
                timeout=120,
            )
            assert proc.returncode == 0, proc.stderr[-2000:]
            reported = "".join(str(record.seq) for record in SeqIO.parse(output_fasta, "fasta"))
            assert ("I" in reported) is expect_isoleucine
