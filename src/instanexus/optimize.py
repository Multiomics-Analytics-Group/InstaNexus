#!/usr/bin/env python
"""CLI entry point for the instanexus-optimize command."""

import sys
from pathlib import Path


def cli() -> None:
    """Launch the hyperparameter grid search optimizer.

    Delegates to scripts/optimization/grid_search.py, which contains
    the full implementation. The package entry point exists here so that
    `instanexus-optimize` is available after `uv sync`.
    """
    _scripts_opt = Path(__file__).resolve().parents[2] / "scripts" / "optimization"
    if str(_scripts_opt) not in sys.path:
        sys.path.insert(0, str(_scripts_opt))

    import grid_search  # type: ignore[import-not-found]

    grid_search.main()
