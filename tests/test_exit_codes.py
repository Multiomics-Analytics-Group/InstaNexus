#!/usr/bin/env python3

"""The command-line interfaces must exit non-zero when a step fails (issue #60)."""

import subprocess
import sys
from pathlib import Path
from types import SimpleNamespace

import pandas as pd
import pytest

from instanexus import assembly, main
from tests.test_dbg_scaffolding_limit import noisy_peptides

# Exit codes documented in the README; hard-coded so these tests also fail without the fix
EXIT_STEP_FAILED = 1
EXIT_USAGE = 2
EXIT_SCAFFOLDING_LIMIT = 3


def run_cli(*args: str) -> subprocess.CompletedProcess:
    """Run `python -m instanexus.main` with the given arguments.

    Args:
        *args: Command-line arguments.

    Returns:
        The completed process, with stdout and stderr captured as text.
    """
    return subprocess.run(
        [sys.executable, "-m", "instanexus.main", *args],
        capture_output=True,
        text=True,
        timeout=120,
    )


def test_missing_input_exits_with_step_failed(tmp_path: Path) -> None:
    proc = run_cli("--input-csv", str(tmp_path / "does_not_exist.csv"), "--output-dir", str(tmp_path / "out"))

    assert proc.returncode == EXIT_STEP_FAILED, proc.stderr[-2000:]
    assert "Preprocessing failed" in proc.stderr


def test_input_without_sequence_column_exits_with_step_failed(tmp_path: Path) -> None:
    input_csv = tmp_path / "no_sequence_column.csv"
    input_csv.write_text("foo,bar\n1,2\n")

    proc = run_cli("--input-csv", str(input_csv), "--output-dir", str(tmp_path / "out"))

    assert proc.returncode == EXIT_STEP_FAILED, proc.stderr[-2000:]
    assert "No suitable column found for peptide sequences" in proc.stderr
    assert "finished successfully" not in proc.stderr


def test_invalid_arguments_exit_with_usage_error() -> None:
    proc = run_cli("--max-scaffold-overlaps", "-1")

    assert proc.returncode == EXIT_USAGE


def test_help_documents_exit_codes() -> None:
    proc = run_cli("--help")

    assert proc.returncode == 0
    assert "exit codes:" in proc.stdout
    assert f"{EXIT_SCAFFOLDING_LIMIT}  dbg scaffolding exceeded --max-scaffold-overlaps" in proc.stdout


def test_exit_code_constants_match_documentation() -> None:
    assert main.EXIT_STEP_FAILED == EXIT_STEP_FAILED
    assert assembly.EXIT_SCAFFOLDING_LIMIT == EXIT_SCAFFOLDING_LIMIT


@pytest.mark.parametrize("step", ["assembly", "clustering", "alignment", "consensus"])
def test_every_step_failure_raises(step: str, tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    input_csv = tmp_path / "sample.csv"
    pd.DataFrame({"preds": ["PEPTIDEK", "PEPTIDER"], "conf": [0.9, 0.9]}).to_csv(input_csv, index=False)
    args = SimpleNamespace(
        input_csv=str(input_csv),
        output_dir=str(tmp_path / "out"),
        folder_outputs=None,
        metadata_json_path=None,
        contaminants_fasta_path=None,
        chain=None,
        reference=False,
        conf=None,
        fdr=None,
        min_length=None,
        max_length=None,
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

    def fail(**kwargs: object) -> None:
        raise RuntimeError(f"{step} broke")

    # Steps before the failing one succeed without running external tools
    for name in ["assembly", "clustering", "alignment", "consensus"]:
        module = getattr(main, name)
        monkeypatch.setattr(module, "main", fail if name == step else lambda **kwargs: None)

    with pytest.raises(main.PipelineStepFailed) as excinfo:
        main.run_pipeline(args)

    assert excinfo.value.step == step.capitalize()
    assert not (tmp_path / "out" / "summary.tsv").exists()


def test_scaffolding_limit_keeps_its_own_exit_code(tmp_path: Path) -> None:
    peptides = noisy_peptides()
    input_csv = tmp_path / "noisy.csv"
    pd.DataFrame({"preds": peptides, "conf": [0.9] * len(peptides)}).to_csv(input_csv, index=False)

    proc = run_cli(
        "--input-csv",
        str(input_csv),
        "--output-dir",
        str(tmp_path / "out"),
        "--assembly-mode",
        "dbg",
        "--kmer-size",
        "7",
        "--min-overlap",
        "3",
        "--size-threshold",
        "10",
    )

    assert proc.returncode == EXIT_SCAFFOLDING_LIMIT, proc.stderr[-2000:]
    assert "--max-scaffold-overlaps" in proc.stderr
