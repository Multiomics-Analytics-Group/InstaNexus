#!/usr/bin/env python3

"""compute_assembly_statistics must use the 0-based, end-exclusive coordinates of map_to_protein (issue #61)."""

from pathlib import Path

import pytest

from instanexus import helpers
from instanexus import visualization as viz

REFERENCE = "ABCDEFGHIJ"  # 10 residues, positions 0-9


def _statistics(sequences: list, tmp_path: Path) -> tuple:
    """Map exact matches to REFERENCE and compute the assembly statistics.

    Args:
        sequences: Sequences to map.
        tmp_path: Folder for the statistics JSON.

    Returns:
        The (start, end) mappings and the statistics dictionary.
    """
    mapped = viz.process_protein_contigs_scaffold(sequences, REFERENCE, max_mismatches=0, min_identity=1.0)
    df = viz.create_dataframe_from_mapped_sequences(mapped)
    statistics = helpers.compute_assembly_statistics(df, "test", str(tmp_path), REFERENCE)

    return list(zip(df["start"], df["end"], strict=True)), statistics


@pytest.mark.parametrize(
    "sequences, expected_mapping, expected_coverage, expected_lengths",
    [
        pytest.param(["DEF"], [(3, 6)], 3 / 10, [3], id="internal"),
        pytest.param(["ABC"], [(0, 3)], 3 / 10, [3], id="at-reference-start"),
        pytest.param(["ABC", "DEF"], [(0, 3), (3, 6)], 6 / 10, [3, 3], id="adjacent"),
        pytest.param(["ABC", "GHI"], [(0, 3), (6, 9)], 6 / 10, [3, 3], id="with-gap"),
        pytest.param([REFERENCE], [(0, 10)], 10 / 10, [10], id="whole-reference"),
        pytest.param(["ABCDE", "CDEFG"], [(0, 5), (2, 7)], 7 / 10, [5, 5], id="overlapping"),
    ],
)
def test_coverage_and_lengths_use_end_exclusive_coordinates(
    sequences: list,
    expected_mapping: list,
    expected_coverage: float,
    expected_lengths: list,
    tmp_path: Path,
) -> None:
    mapping, statistics = _statistics(sequences, tmp_path)

    # map_to_protein: 0-based start, exclusive end
    assert mapping == expected_mapping
    assert statistics["coverage"] == pytest.approx(expected_coverage)
    assert statistics["min_length"] == min(expected_lengths)
    assert statistics["max_length"] == max(expected_lengths)
    assert statistics["average_length"] == pytest.approx(sum(expected_lengths) / len(expected_lengths))
    assert statistics["N50"] == max(expected_lengths)
    assert statistics["reference_start"] == 0
    assert statistics["reference_end"] == len(REFERENCE)
