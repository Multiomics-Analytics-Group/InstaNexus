#!/usr/bin/env python3

"""`python -m instanexus` must run the pipeline CLI, and --version must report the installed version (issue #74)."""

import subprocess
import sys
from importlib.metadata import version
from pathlib import Path

import pytest


def run(*args: str) -> subprocess.CompletedProcess:
    """Run a command with the test interpreter's environment.

    Args:
        *args: The command and its arguments.

    Returns:
        The completed process, with stdout and stderr captured as text.
    """
    return subprocess.run(list(args), capture_output=True, text=True, timeout=120)


def test_python_m_instanexus_runs_the_pipeline_cli() -> None:
    proc = run(sys.executable, "-m", "instanexus", "--help")

    assert proc.returncode == 0, proc.stderr[-2000:]
    assert "--input-csv" in proc.stdout
    assert "--max-scaffold-overlaps" in proc.stdout


@pytest.mark.parametrize("module", ["instanexus", "instanexus.main"])
def test_version_comes_from_package_metadata(module: str) -> None:
    proc = run(sys.executable, "-m", module, "--version")

    assert proc.returncode == 0, proc.stderr[-2000:]
    assert proc.stdout.strip() == f"InstaNexus {version('instanexus')}"
    assert proc.stdout.strip() != "InstaNexus 0.1.0"


def test_console_script_reports_version() -> None:
    script = Path(sys.executable).parent / "instanexus"
    if not script.exists():
        pytest.skip("instanexus console script not installed in this environment")

    proc = run(str(script), "--version")

    assert proc.returncode == 0, proc.stderr[-2000:]
    assert proc.stdout.strip() == f"InstaNexus {version('instanexus')}"
