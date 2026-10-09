# Changelog

Changes that are merged but not yet released are collected under **Unreleased** and copied into the
notes of the next [GitHub Release](https://github.com/Multiomics-Analytics-Group/InstaNexus/releases).
Notes for earlier versions are on the GitHub Releases page.

## Unreleased

### Changed: assembly statistics are not comparable with earlier versions

Reference-based statistics (`peptide_stats.json`, `contig_stats.json`, `scaffolds_stats.json`, written with
`--reference`) used the 0-based, end-exclusive mapping coordinates as if they were 1-based and inclusive
([#61](https://github.com/Multiomics-Analytics-Group/InstaNexus/issues/61)). They are now computed correctly,
so the following values change and **cannot be compared with results from 0.3.1 or earlier**:

| Metric | Direction | Size of the change |
|---|---|---|
| `coverage` | lower | Previously `(C + B) / (L + 1)` instead of `C / L` (`C` covered positions, `B` contiguous covered blocks, `L` reference length). On the nanobody demo data: 0.33-2.15 percentage points lower; the more fragmented the coverage, the larger the drop. Full coverage stays at 1.0. |
| `average_length`, `min_length`, `max_length` | lower | Exactly 1 residue lower. |
| `N50`, `N90` | lower or unchanged | Usually 1 residue lower. |
| `reference_end` | lower | Now the reference length (it was length + 1). |

`total_sequences`, `mean_identity`, `median_identity` and `perfect_matches` are not affected by this fix
(`total_mismatches` changes for a different reason, see below).
Parameter rankings from `scripts/optimization/grid_search.py`, which weights `coverage` and `N50`, may change.

### Changed: `total_mismatches` counts mismatches; new `mismatched_positions`

`total_mismatches` used to be the number of distinct mismatch offsets *within* the mapped sequences, so mismatches
of different sequences at the same offset were merged
([#63](https://github.com/Multiomics-Analytics-Group/InstaNexus/issues/63)). The same statistics files now report
two mismatch metrics:

| Metric | Meaning | Change |
|---|---|---|
| `total_mismatches` | Total number of mismatches over all mapped sequences. A reference position covered by several sequences with the same error counts once per sequence. | **New value, not comparable with earlier versions.** Higher or equal, never lower; on the nanobody demo data from 3-24 to 4-124, most for redundant sets such as peptides (12 → 90 at `conf > 0.5`). |
| `mismatched_positions` | Number of reference positions with a mismatch in at least one mapped sequence, using the same coordinates as `coverage`; at most the number of covered positions. | **New.** 3-22 on the nanobody demo data. |

`scripts/optimization/grid_search.py` also reports `mismatched_positions` in its results table; the composite score
does not use either metric and is unchanged.

### Changed: `dbg` stops with exit code 3 instead of running indefinitely on noisy input

On noisy, low-confidence input, `dbg` scaffolding merged every pair of overlapping contigs round after round, and
the number of overlaps grew combinatorially until the run never finished
([#58](https://github.com/Multiomics-Analytics-Group/InstaNexus/issues/58)).

- **New limit, `--max-scaffold-overlaps`** (on `instanexus` and `python -m instanexus.assembly`; also
  `Assembler(max_scaffold_overlaps=...)`): if one `dbg` scaffolding round finds more overlaps than this
  (default 100,000; `0` disables the limit), the assembly stops with an error naming the limit and the usual
  remedies (filter with `--conf`/`--fdr`, a larger `--min-overlap`, or `--assembly-mode dbg_weighted`).
- **Exit code 3** when the limit is reached, so pipelines can tell this case from other failures. This is a
  behaviour change for `dbg` users: inputs that exceed the limit used to run indefinitely.
- **`POTENTIAL OVERLAP MISSED` is logged at DEBUG** instead of INFO; it filled the log with tens of thousands of
  lines on noisy input.
- **`dbg` and `greedy` are faster**, with identical output. On the nanobody demo data, `dbg` at `conf > 0.3`
  went from 0.98 s to 0.02 s and `greedy` at `conf > 0.1` from 4.62 s to 0.72 s. At `conf > 0.1` (784 PSMs),
  `dbg` used to run past 600 s; it now stops at the limit after 0.34 s.
