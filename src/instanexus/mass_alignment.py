#!/usr/bin/env python

r"""Mass-aware alignment of a de novo sequence to a reference protein.

The scoring follows the peptide similarity alignment (PSA) of Cemil Can Saylan, Armin Suleymaniniya and
Mathias Wilhelm, "Beyond Best Guesses: A Multi-Level Confidence Assessment Framework for De Novo Peptide
Sequencing Results", HUPO poster PS02.146, Computational Mass Spectrometry, TUM School of Life Sciences,
Technical University of Munich (not yet published; cite the paper once it appears). No code was published with
the poster; this is a reimplementation from its Figure 1c, and work that uses it must cite the authors.

A Needleman-Wunsch alignment with:

- match +2, mismatch -1, gap -2 (linear);
- a *shuffle* step: blocks of 2 to 4 residues with the same residues in a different order, +1.0 per residue;
- an *isobaric* step: blocks of up to two residues on each side, with different residues but masses within
  0.01 Da (N and GG, Q and AG), +1.5 per residue of the longer block.

The alignment is semi-global: the whole sequence is aligned, to any region of the protein. That lets a scaffold
with one inserted or deleted residue, an isobaric block or a swap map where an ungapped comparison does not.

Three points the poster leaves open are settled here as follows. Blocks that are a permutation of each other
are shuffles, not isobaric blocks, even though their masses agree. "Per residue" of an isobaric block counts the
residues of its longer side. Residue masses are monoisotopic and unmodified, since modifications are removed in
preprocessing.

__authors__ = Jeroen Van Goey
__copyright__ = Copyright 2026
__research-group__ = DTU Biosustain (Multi-omics Network Analytics) and DTU Bioengineering
__date__ = 30 Oct 2026
__maintainer__ = Marco Reverenna
__email__ = marcor@dtu.dk
__status__ = Dev
"""

from collections import Counter
from dataclasses import dataclass, field
from typing import List, Optional, Tuple

import numpy as np

MATCH = 2.0
MISMATCH = -1.0
GAP = -2.0
SHUFFLE_PER_RESIDUE = 1.0
ISOBARIC_PER_RESIDUE = 1.5
ISOBARIC_TOLERANCE_DA = 0.01
SHUFFLE_BLOCK_SIZES = (2, 3, 4)
ISOBARIC_BLOCK_SIZES = ((1, 1), (1, 2), (2, 1), (2, 2))

# monoisotopic residue masses, unmodified
RESIDUE_MASSES = {
    "G": 57.021464,
    "A": 71.037114,
    "S": 87.032028,
    "P": 97.052764,
    "V": 99.068414,
    "T": 101.047679,
    "C": 103.009185,
    "L": 113.084064,
    "I": 113.084064,
    "N": 114.042927,
    "D": 115.026943,
    "Q": 128.058578,
    "K": 128.094963,
    "E": 129.042593,
    "M": 131.040485,
    "H": 137.058912,
    "F": 147.068414,
    "U": 150.953636,
    "R": 156.101111,
    "Y": 163.063329,
    "W": 186.079313,
    "O": 237.147727,
}

# the kinds of difference an alignment reports, as keys of `Alignment.differences`
DIFFERENCE_KINDS = ("substitutions", "insertions", "deletions", "isobaric_blocks", "shuffled_blocks")

# transitions of the dynamic programme: (kind, residues of the sequence, residues of the protein)
_TRANSITIONS: List[Tuple[str, int, int]] = (
    [("diagonal", 1, 1), ("insertion", 1, 0)]
    + [("isobaric", a, b) for a, b in ISOBARIC_BLOCK_SIZES]
    + [("shuffle", k, k) for k in SHUFFLE_BLOCK_SIZES]
)
_DELETION = len(_TRANSITIONS)


@dataclass
class Operation:
    """One step of an alignment.

    Attributes:
        kind: match, substitution, insertion (a residue of the sequence that is not in the protein), deletion
            (a residue of the protein that is not in the sequence), isobaric or shuffle.
        seq_start: First residue of the sequence in the step, 0-based.
        seq_end: End of the step in the sequence, exclusive.
        ref_start: First residue of the protein in the step, 0-based.
        ref_end: End of the step in the protein, exclusive.
    """

    kind: str
    seq_start: int
    seq_end: int
    ref_start: int
    ref_end: int

    @property
    def columns(self) -> int:
        """Alignment columns the step takes: the residues of its longer side."""
        return max(self.seq_end - self.seq_start, self.ref_end - self.ref_start)


@dataclass
class Alignment:
    """A sequence aligned to a region of a protein.

    Attributes:
        start: First residue of the protein region, 0-based.
        end: End of the protein region, exclusive. `end - start` differs from the sequence length when the
            alignment has gaps or isobaric blocks.
        score: Alignment score.
        operations: The steps, in order along the sequence.
    """

    start: int
    end: int
    score: float
    operations: List[Operation] = field(default_factory=list)

    @property
    def columns(self) -> int:
        """Alignment length, in columns."""
        return sum(op.columns for op in self.operations)

    @property
    def identity(self) -> float:
        """Share of the alignment columns that are exact matches."""
        matches = sum(op.columns for op in self.operations if op.kind == "match")

        return matches / self.columns if self.columns else 0.0

    @property
    def difference_offsets(self) -> List[int]:
        """One entry per alignment column that is not an exact match: its offset in the protein from `start`.

        A column without a protein residue, an inserted residue or the extra column of an isobaric block,
        takes the offset of the next protein residue, or of the last one at the end of the region.
        """
        offsets: List[int] = []
        for op in self.operations:
            if op.kind == "match":
                continue
            ref_residues = list(range(op.ref_start, op.ref_end))
            if not ref_residues:
                ref_residues = [min(op.ref_start, self.end - 1)]
            padding = [ref_residues[-1]] * (op.columns - len(ref_residues))
            offsets.extend(position - self.start for position in ref_residues + padding)

        return offsets

    @property
    def residue_positions(self) -> List[Optional[int]]:
        """For each residue of the sequence, the protein position it is aligned to.

        Only residues in a one-to-one column, a match or a substitution, have a position; residues of an
        insertion, an isobaric block or a shuffle have None.
        """
        positions: List[Optional[int]] = []
        for op in self.operations:
            if op.kind in ("match", "substitution"):
                positions.append(op.ref_start)
            else:
                positions.extend([None] * (op.seq_end - op.seq_start))

        return positions

    @property
    def differences(self) -> Counter:
        """Counts of each kind of difference: substituted, inserted and deleted residues, isobaric and shuffled
        blocks."""
        counts: Counter = Counter({kind: 0 for kind in DIFFERENCE_KINDS})
        for op in self.operations:
            if op.kind == "substitution":
                counts["substitutions"] += 1
            elif op.kind == "insertion":
                counts["insertions"] += 1
            elif op.kind == "deletion":
                counts["deletions"] += 1
            elif op.kind == "isobaric":
                counts["isobaric_blocks"] += 1
            elif op.kind == "shuffle":
                counts["shuffled_blocks"] += 1

        return counts


def _codes(sequence: str) -> np.ndarray:
    """Encode residues as integers 1-26, so that blocks can be compared as numbers."""
    return np.frombuffer(sequence.encode("ascii"), dtype=np.uint8).astype(np.int64) - 64


def _masses(sequence: str) -> np.ndarray:
    """Residue masses, NaN for residues without a known mass, which then match no isobaric block."""
    return np.array([RESIDUE_MASSES.get(residue, np.nan) for residue in sequence])


def _block_keys(codes: np.ndarray, size: int) -> Tuple[np.ndarray, np.ndarray]:
    """For every block of `size` residues: a number for the block, and one for its residues in sorted order.

    Args:
        codes: Residue codes, from `_codes`.
        size: Block length.

    Returns:
        The two keys of each block, by block start.
    """
    if len(codes) < size:
        return np.empty(0, dtype=np.int64), np.empty(0, dtype=np.int64)
    windows = np.lib.stride_tricks.sliding_window_view(codes, size)
    weights = 32 ** np.arange(size - 1, -1, -1)

    return windows @ weights, np.sort(windows, axis=1) @ weights


def _block_masses(masses: np.ndarray, size: int) -> np.ndarray:
    """Mass of every block of `size` residues, by block start."""
    if len(masses) < size:
        return np.empty(0)

    return np.asarray(np.lib.stride_tricks.sliding_window_view(masses, size).sum(axis=1))


def align(
    sequence: str,
    protein: str,
    isobaric_tolerance: float = ISOBARIC_TOLERANCE_DA,
) -> Optional[Alignment]:
    """Align the whole of `sequence` to the best-scoring region of `protein`.

    Ties between regions go to the one that ends first in the protein. A sequence that occurs in the protein
    is returned as its first occurrence without running the alignment, since no alignment scores higher.

    Args:
        sequence: Sequence to align, such as a scaffold.
        protein: Reference protein.
        isobaric_tolerance: Largest mass difference, in Da, between the blocks of an isobaric step.

    Returns:
        The alignment, or None when `sequence` is empty or aligns to no protein residue.
    """
    m, n = len(sequence), len(protein)
    if m == 0 or n == 0:
        return None
    first = protein.find(sequence)
    if first >= 0:
        return Alignment(
            first, first + m, MATCH * m, [Operation("match", i, i + 1, first + i, first + i + 1) for i in range(m)]
        )

    seq_codes, ref_codes = _codes(sequence), _codes(protein)
    seq_masses, ref_masses = _masses(sequence), _masses(protein)
    ref_blocks = {size: _block_keys(ref_codes, size) for size in (1, 2, 3, 4)}
    ref_block_masses = {size: _block_masses(ref_masses, size) for size in (1, 2)}

    score = np.full((m + 1, n + 1), -np.inf)
    # (m + 1, n + 1): index into _TRANSITIONS, or _DELETION
    pointer = np.zeros((m + 1, n + 1), dtype=np.int8)
    # the alignment may start anywhere in the protein
    score[0, :] = 0.0
    column = np.arange(n + 1)

    for i in range(1, m + 1):
        # (n + 1,) best score into each cell of row i from the rows above, and the transition that gives it
        best = np.full(n + 1, -np.inf)
        best_from = np.zeros(n + 1, dtype=np.int8)
        for t, (kind, a, b) in enumerate(_TRANSITIONS):
            if a > i or b > n:
                continue
            previous = score[i - a]
            # (n + 1,)
            candidate = np.full(n + 1, -np.inf)
            if kind == "insertion":
                candidate = previous + GAP
            else:
                seq_block = seq_codes[i - a : i]
                ref_key, ref_sorted = ref_blocks[b]
                if kind == "diagonal":
                    step = np.where(ref_codes == seq_codes[i - 1], MATCH, MISMATCH)
                elif kind == "shuffle":
                    seq_key, seq_sorted = _block_keys(seq_block, a)
                    step = np.where(
                        (ref_sorted == seq_sorted[0]) & (ref_key != seq_key[0]), SHUFFLE_PER_RESIDUE * a, -np.inf
                    )
                else:
                    _, seq_sorted = _block_keys(seq_block, a)
                    seq_mass = seq_masses[i - a : i].sum()
                    different = ref_sorted != seq_sorted[0] if a == b else np.ones(len(ref_key), dtype=bool)
                    isobaric = np.abs(ref_block_masses[b] - seq_mass) <= isobaric_tolerance
                    step = np.where(different & isobaric, ISOBARIC_PER_RESIDUE * max(a, b), -np.inf)
                # the block over protein residues j - b to j ends at column j
                candidate[b:] = previous[: n + 1 - b] + step
            better = candidate > best
            best = np.where(better, candidate, best)
            best_from = np.where(better, t, best_from)

        # a deleted protein residue moves along the row: running maximum of best[k] + GAP * (j - k)
        with_deletions = np.maximum.accumulate(best - GAP * column) + GAP * column
        score[i] = with_deletions
        pointer[i] = np.where(with_deletions > best, _DELETION, best_from)

    end = int(np.argmax(score[m]))
    operations = []
    i, j = m, end
    while i > 0:
        t = int(pointer[i, j])
        if t == _DELETION:
            operations.append(Operation("deletion", i, i, j - 1, j))
            j -= 1
            continue
        kind, a, b = _TRANSITIONS[t]
        if kind == "diagonal":
            kind = "match" if sequence[i - 1] == protein[j - 1] else "substitution"
        operations.append(Operation(kind, i - a, i, j - b, j))
        i, j = i - a, j - b
    operations.reverse()
    if j == end:
        return None

    return Alignment(j, end, float(score[m, end]), operations)
