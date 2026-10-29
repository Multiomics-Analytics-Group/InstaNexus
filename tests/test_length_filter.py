#!/usr/bin/env python3

"""Opt-in peptide length filter in preprocessing (issue #68)."""

import json
import subprocess
import sys
from pathlib import Path

import pandas as pd
import pytest

from instanexus import preprocessing

# 5 to 9 residues once modifications are removed; the second carries a modification that must not count
PEPTIDES = ["ACDEF", "C[UNIMOD:4]PEPTIDEK", "GHIKLMN", "PQRSTVWY", "ACDEFGHIK"]
LENGTHS = [5, 9, 7, 8, 9]


def run_preprocessing(tmp_path: Path, *args: str) -> subprocess.CompletedProcess:
    """Run `python -m instanexus.preprocessing` on PEPTIDES.

    Args:
        tmp_path: Folder for the input and output files.
        *args: Extra command-line arguments.

    Returns:
        The completed process, with stdout and stderr captured as text.
    """
    input_csv = tmp_path / "sample.csv"
    pd.DataFrame({"prediction": PEPTIDES, "conf": [0.9] * len(PEPTIDES)}).to_csv(input_csv, index=False)

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


def test_filter_by_length_bounds_are_inclusive() -> None:
    df = pd.DataFrame({"cleaned_preds": ["A" * n for n in range(5, 11)]})

    kept = preprocessing.filter_by_length(df, min_length=7, max_length=9)

    assert kept["cleaned_preds"].str.len().tolist() == [7, 8, 9]


def test_filter_by_length_one_sided() -> None:
    df = pd.DataFrame({"cleaned_preds": ["A" * n for n in range(5, 11)]})

    assert preprocessing.filter_by_length(df, 9, None)["cleaned_preds"].str.len().tolist() == [9, 10]
    assert preprocessing.filter_by_length(df, None, 6)["cleaned_preds"].str.len().tolist() == [5, 6]
    assert len(preprocessing.filter_by_length(df, None, None)) == len(df)


def test_filter_by_length_rejects_inverted_bounds() -> None:
    with pytest.raises(ValueError, match="greater than"):
        preprocessing.filter_by_length(pd.DataFrame({"cleaned_preds": ["AAAA"]}), 10, 5)


def test_off_by_default_and_recorded(tmp_path: Path) -> None:
    proc = run_preprocessing(tmp_path)

    assert proc.returncode == 0, proc.stderr[-2000:]
    cleaned = pd.read_csv(tmp_path / "cleaned.csv")
    assert sorted(cleaned["cleaned_preds"].str.len()) == sorted(LENGTHS)
    record = json.loads((tmp_path / "cleaned_filters.json").read_text())
    assert record["min_length"] is None
    assert record["max_length"] is None
    assert record["rows_removed_by_length"] == 0
    assert record["rows_output"] == len(PEPTIDES)


def test_bounds_apply_to_residues_and_are_recorded(tmp_path: Path) -> None:
    proc = run_preprocessing(tmp_path, "--min-length", "7", "--max-length", "8")

    assert proc.returncode == 0, proc.stderr[-2000:]
    cleaned = pd.read_csv(tmp_path / "cleaned.csv")
    # the modified peptide counts as 9 residues, not as its 19 characters
    assert sorted(cleaned["cleaned_preds"].str.len()) == [7, 8]
    record = json.loads((tmp_path / "cleaned_filters.json").read_text())
    assert (record["min_length"], record["max_length"]) == (7, 8)
    assert record["rows_after_cleaning"] == len(PEPTIDES)
    assert record["rows_removed_by_length"] == 3
    assert record["rows_output"] == 2


@pytest.mark.parametrize("args", [["--min-length", "10", "--max-length", "5"], ["--min-length", "0"]])
def test_invalid_bounds_are_usage_errors(tmp_path: Path, args: list) -> None:
    proc = run_preprocessing(tmp_path, *args)

    assert proc.returncode == 2
    assert not (tmp_path / "cleaned.csv").exists()


def test_pipeline_cli_accepts_the_bounds(tmp_path: Path) -> None:
    proc = subprocess.run(
        [sys.executable, "-m", "instanexus.main", "--help"], capture_output=True, text=True, timeout=120
    )

    assert "--min-length" in proc.stdout
    assert "--max-length" in proc.stdout
