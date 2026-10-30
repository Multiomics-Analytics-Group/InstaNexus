#!/usr/bin/env python3

"""Opt-in peptide length filter in preprocessing (issue #68)."""

import json
import subprocess
import sys
from pathlib import Path
from types import SimpleNamespace

import pandas as pd
import pytest

from instanexus import assembly, main, preprocessing

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
    length_step = next(step for step in record["steps"] if step["name"] == "length")
    assert length_step == {"name": "length", "applied": False, "rows_in": len(PEPTIDES), "rows_out": len(PEPTIDES)}
    assert record["rows_output"] == len(PEPTIDES)


def test_bounds_apply_to_residues_and_are_recorded(tmp_path: Path) -> None:
    proc = run_preprocessing(tmp_path, "--min-length", "7", "--max-length", "8")

    assert proc.returncode == 0, proc.stderr[-2000:]
    cleaned = pd.read_csv(tmp_path / "cleaned.csv")
    # the modified peptide counts as 9 residues, not as its 19 characters
    assert sorted(cleaned["cleaned_preds"].str.len()) == [7, 8]
    record = json.loads((tmp_path / "cleaned_filters.json").read_text())
    assert (record["min_length"], record["max_length"]) == (7, 8)
    length_step = next(step for step in record["steps"] if step["name"] == "length")
    assert length_step == {"name": "length", "applied": True, "rows_in": len(PEPTIDES), "rows_out": 2}
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


def test_record_lists_rows_after_every_filter_in_order(tmp_path: Path) -> None:
    input_csv = tmp_path / "sample.csv"
    # conf 0.5 on the last two peptides, so the conf filter removes them after the length filter
    pd.DataFrame({"prediction": PEPTIDES, "conf": [0.9, 0.9, 0.9, 0.5, 0.5]}).to_csv(input_csv, index=False)

    proc = subprocess.run(
        [
            sys.executable,
            "-m",
            "instanexus.preprocessing",
            "--input-csv",
            str(input_csv),
            "--output-csv-path",
            str(tmp_path / "cleaned.csv"),
            "--min-length",
            "6",
            "--conf",
            "0.8",
        ],
        capture_output=True,
        text=True,
        timeout=120,
    )

    assert proc.returncode == 0, proc.stderr[-2000:]
    record = json.loads((tmp_path / "cleaned_filters.json").read_text())
    assert record["input_csv"] == "sample.csv"
    assert [(s["name"], s["applied"], s["rows_in"], s["rows_out"]) for s in record["steps"]] == [
        ("cleaning", True, 5, 5),
        ("length", True, 5, 4),
        ("contaminants", False, 4, 4),
        ("conf", True, 4, 2),
        ("fdr", False, 2, 2),
        ("min_reads", False, 2, 2),
        ("top_fraction", False, 2, 2),
    ]
    assert record["rows_output"] == 2


def test_length_bounds_are_in_the_summary_columns(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    input_csv = tmp_path / "sample.csv"
    pd.DataFrame({"preds": ["PEPTIDEKAAGLR", "AAGLRSTWYK"], "conf": [0.9, 0.9]}).to_csv(input_csv, index=False)
    output_dir = tmp_path / "out"

    def write_consensus_stats(**kwargs: object) -> None:
        stats_path = output_dir / "scaffolds" / "consensus" / "consensus_stats.json"
        stats_path.parent.mkdir(parents=True, exist_ok=True)
        stats_path.write_text(json.dumps({"n_consensus_files": 0, "max_consensus_length": 0}))

    # clustering and alignment need external tools; the consensus step only has to leave its statistics
    monkeypatch.setattr(main.clustering, "main", lambda **kwargs: None)
    monkeypatch.setattr(main.alignment, "main", lambda **kwargs: None)
    monkeypatch.setattr(main.consensus, "main", write_consensus_stats)
    args = SimpleNamespace(
        input_csv=str(input_csv),
        output_dir=str(output_dir),
        folder_outputs=None,
        metadata_json_path=None,
        contaminants_fasta_path=None,
        chain=None,
        reference=False,
        conf=None,
        fdr=None,
        min_length=7,
        max_length=20,
        top_fraction=None,
        min_reads=None,
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

    summary = pd.read_csv(output_dir / "summary.tsv", sep="\t")
    assert list(summary.columns)[-4:-2] == ["min_length", "max_length"]
    assert "n_consensus_files" in summary.columns
    assert (summary["min_length"][0], summary["max_length"][0]) == (7, 20)
