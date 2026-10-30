#!/usr/bin/env python3

"""Opt-in rank cutoff (--top-fraction) and read-support filter (--min-reads) in preprocessing (issue #67)."""

import json
import subprocess
import sys
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest

from instanexus import assembly, main, preprocessing

# PEPTIDEK is read three times: once with I spelled as L, once with a modification; ACDEFGHK twice
PSMS = pd.DataFrame(
    {
        "prediction": [
            "PEPTIDEK",
            "PEPTLDEK",
            "PEPT[UNIMOD:21]IDEK",
            "ACDEFGHK",
            "ACDEFGHK",
            "MNPQRSTK",
            "VWYAAGLR",
        ],
        "conf": [0.99, 0.70, 0.65, 0.90, 0.60, 0.98, 0.50],
    }
)


def run_preprocessing(tmp_path: Path, *args: str, psms: pd.DataFrame = PSMS) -> subprocess.CompletedProcess:
    """Run `python -m instanexus.preprocessing` on a PSM table.

    Args:
        tmp_path: Folder for the input and output files.
        *args: Extra command-line arguments.
        psms: PSMs to write as the input CSV.

    Returns:
        The completed process, with stdout and stderr captured as text.
    """
    input_csv = tmp_path / "sample.csv"
    psms.to_csv(input_csv, index=False)

    return subprocess.run(
        [
            sys.executable,
            "-m",
            "instanexus.preprocessing",
            "--input-csv",
            str(input_csv),
            "--output-csv-path",
            str(tmp_path / "cleaned.csv"),
            *args,
        ],
        capture_output=True,
        text=True,
        timeout=120,
    )


def kept_confidences(tmp_path: Path) -> list:
    """Return the confidences of the PSMs preprocessing wrote, highest first.

    Args:
        tmp_path: Folder holding cleaned.csv.

    Returns:
        The confidences.
    """
    return sorted(pd.read_csv(tmp_path / "cleaned.csv")["conf"], reverse=True)


def record_steps(tmp_path: Path) -> dict:
    """Return the steps of the filters record, by name.

    Args:
        tmp_path: Folder holding cleaned_filters.json.

    Returns:
        Each step as (applied, rows_in, rows_out).
    """
    record = json.loads((tmp_path / "cleaned_filters.json").read_text())

    return {s["name"]: (s["applied"], s["rows_in"], s["rows_out"]) for s in record["steps"]}


def test_top_fraction_mask_keeps_the_floor_of_the_fraction() -> None:
    df = pd.DataFrame({"conf": [0.1, 0.9, 0.5, 0.7, 0.3]})

    assert df[preprocessing.top_fraction_mask(df, 0.5)]["conf"].tolist() == [0.9, 0.7]
    assert df[preprocessing.top_fraction_mask(df, 0.6)]["conf"].tolist() == [0.9, 0.5, 0.7]
    assert preprocessing.top_fraction_mask(df, 1.0).all()


def test_top_fraction_mask_breaks_ties_by_row_order_and_drops_missing_confidence() -> None:
    df = pd.DataFrame({"conf": [0.8, np.nan, 0.8, 0.8, 0.2]})

    # four PSMs have a confidence, so half of them is two: the first two of the tied rows
    assert preprocessing.top_fraction_mask(df, 0.5).tolist() == [True, False, True, False, False]
    assert not preprocessing.top_fraction_mask(df, 1.0)[1]


@pytest.mark.parametrize("fraction", [0.0, 1.5])
def test_top_fraction_mask_rejects_fractions_outside_the_unit_interval(fraction: float) -> None:
    with pytest.raises(ValueError, match="must be in"):
        preprocessing.top_fraction_mask(pd.DataFrame({"conf": [0.5]}), fraction)


def test_off_by_default_and_recorded(tmp_path: Path) -> None:
    proc = run_preprocessing(tmp_path)

    assert proc.returncode == 0, proc.stderr[-2000:]
    assert len(kept_confidences(tmp_path)) == len(PSMS)
    record = json.loads((tmp_path / "cleaned_filters.json").read_text())
    assert (record["top_fraction"], record["min_reads"]) == (None, None)
    steps = record_steps(tmp_path)
    assert steps["min_reads"] == (False, len(PSMS), len(PSMS))
    assert steps["top_fraction"] == (False, len(PSMS), len(PSMS))


@pytest.mark.parametrize(("min_reads", "expected"), [(2, [0.99, 0.90, 0.70, 0.65, 0.60]), (3, [0.99, 0.70, 0.65])])
def test_min_reads_counts_peptides_without_modifications_and_with_i_as_l(
    tmp_path: Path, min_reads: int, expected: list
) -> None:
    proc = run_preprocessing(tmp_path, "--min-reads", str(min_reads))

    assert proc.returncode == 0, proc.stderr[-2000:]
    assert kept_confidences(tmp_path) == expected
    assert record_steps(tmp_path)["min_reads"] == (True, len(PSMS), len(expected))


def test_top_fraction_keeps_the_most_confident(tmp_path: Path) -> None:
    proc = run_preprocessing(tmp_path, "--top-fraction", "0.5")

    assert proc.returncode == 0, proc.stderr[-2000:]
    # half of 7 PSMs, rounded down
    assert kept_confidences(tmp_path) == [0.99, 0.98, 0.90]
    assert record_steps(tmp_path)["top_fraction"] == (True, len(PSMS), 3)


def test_both_filters_count_over_the_same_psms(tmp_path: Path) -> None:
    proc = run_preprocessing(tmp_path, "--min-reads", "2", "--top-fraction", "0.4")

    assert proc.returncode == 0, proc.stderr[-2000:]
    # 40% of 7 PSMs is 2: 0.99 and 0.98, and MNPQRSTK (0.98) is read once. Ranking only the 5 PSMs read twice
    # would keep 0.99 and 0.90; counting reads within the top 2 would keep nothing.
    assert kept_confidences(tmp_path) == [0.99]
    steps = record_steps(tmp_path)
    assert steps["min_reads"] == (True, 7, 5)
    assert steps["top_fraction"] == (True, 5, 1)


def test_reads_are_counted_after_the_conf_filter(tmp_path: Path) -> None:
    proc = run_preprocessing(tmp_path, "--conf", "0.62", "--min-reads", "2")

    assert proc.returncode == 0, proc.stderr[-2000:]
    # the conf filter removes the second ACDEFGHK, so only PEPTIDEK is still read twice
    assert kept_confidences(tmp_path) == [0.99, 0.70, 0.65]
    steps = record_steps(tmp_path)
    assert steps["conf"] == (True, 7, 5)
    assert steps["min_reads"] == (True, 5, 3)


def test_top_fraction_without_a_confidence_column_is_skipped(tmp_path: Path) -> None:
    proc = run_preprocessing(tmp_path, "--top-fraction", "0.5", psms=PSMS[["prediction"]])

    assert proc.returncode == 0, proc.stderr[-2000:]
    assert "'conf' column not found" in proc.stderr
    assert len(pd.read_csv(tmp_path / "cleaned.csv")) == len(PSMS)
    assert record_steps(tmp_path)["top_fraction"] == (False, len(PSMS), len(PSMS))


@pytest.mark.parametrize("args", [["--top-fraction", "0"], ["--top-fraction", "1.5"], ["--min-reads", "0"]])
def test_invalid_values_are_usage_errors(tmp_path: Path, args: list) -> None:
    proc = run_preprocessing(tmp_path, *args)

    assert proc.returncode == 2
    assert not (tmp_path / "cleaned.csv").exists()


def test_pipeline_cli_accepts_both_options() -> None:
    proc = subprocess.run(
        [sys.executable, "-m", "instanexus.main", "--help"], capture_output=True, text=True, timeout=120
    )

    assert "--top-fraction" in proc.stdout
    assert "--min-reads" in proc.stdout


def test_folder_name_and_summary_record_both_options(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    input_csv = tmp_path / "sample.csv"
    pd.DataFrame({"preds": ["PEPTIDEKAAGLR"] * 2 + ["AAGLRSTWYK"] * 2, "conf": [0.9] * 4}).to_csv(
        input_csv, index=False
    )
    folder_outputs = tmp_path / "out"

    # clustering, alignment and consensus need external tools; the summary is written without them
    monkeypatch.setattr(main.clustering, "main", lambda **kwargs: None)
    monkeypatch.setattr(main.alignment, "main", lambda **kwargs: None)
    monkeypatch.setattr(main.consensus, "main", lambda **kwargs: None)
    args = SimpleNamespace(
        input_csv=str(input_csv),
        output_dir=None,
        folder_outputs=str(folder_outputs),
        metadata_json_path=None,
        contaminants_fasta_path=None,
        chain=None,
        reference=False,
        conf=None,
        fdr=None,
        min_length=None,
        max_length=None,
        top_fraction=0.5,
        min_reads=2,
        assembly_mode="greedy",
        kmer_size=7,
        min_overlap=3,
        size_threshold=10,
        min_identity=0.8,
        max_mismatches=10,
        refine=False,
        min_seq_id=0.85,
        coverage=0.8,
        skip_plots=True,
        max_scaffold_overlaps=assembly.MAX_SCAFFOLD_OVERLAPS,
        isoleucine_restoration=True,
    )

    main.run_pipeline(args)

    experiment_folder = folder_outputs / "sample" / "greedy_reads2_top0.5"
    summary = pd.read_csv(experiment_folder / "summary.tsv", sep="\t")
    assert list(summary.columns)[-2:] == ["top_fraction", "min_reads"]
    assert (summary["top_fraction"][0], summary["min_reads"][0]) == (0.5, 2)
