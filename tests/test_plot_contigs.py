#!/usr/bin/env python3

"""plot_contigs must put non-overlapping contigs on one track and overlapping ones on separate tracks (issue #62)."""

from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.colors as mcolors
import pytest

from instanexus import visualization as viz

REFERENCE = "ABCDEFGHIJ"
CONTIG_COLOR = mcolors.to_rgba("#007EA7")


def contig_tracks(contigs: list, tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> list:
    """Plot contigs mapped to REFERENCE and return the track of each contig, in input order.

    Args:
        contigs: Contig sequences; each must occur in REFERENCE.
        tmp_path: Folder for the plot file.
        monkeypatch: Keeps the figure open after plot_contigs saves it.

    Returns:
        One track index per contig, numbered from 0 in order of first use.
    """
    mapped = [(contig, viz.map_to_protein(contig, REFERENCE, max_mismatches=0, min_identity=1.0)) for contig in contigs]
    monkeypatch.setattr(viz.plt, "close", lambda *args, **kwargs: None)

    viz.plot_contigs(mapped, REFERENCE, "test", str(tmp_path / "contigs.png"))

    ax = viz.plt.gcf().axes[0]
    ys = [p.get_y() for p in ax.patches if p.get_facecolor() == CONTIG_COLOR]
    viz.plt.close("all")
    order = list(dict.fromkeys(ys))

    return [order.index(y) for y in ys]


@pytest.mark.parametrize(
    "contigs, expected",
    [
        (["ABC", "DEF"], [0, 0]),  # adjacent: (0, 3) and (3, 6) do not overlap
        (["ABC", "GHI"], [0, 0]),  # separated
        (["ABCD", "DEFG"], [0, 1]),  # share D
        (["CDE", "ABCDEFG"], [0, 1]),  # the second contains the first
        (["ABC", "DEF", "CD"], [0, 0, 1]),  # CD overlaps both contigs on track 0
    ],
)
def test_track_assignment(contigs: list, expected: list, tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    assert contig_tracks(contigs, tmp_path, monkeypatch) == expected
