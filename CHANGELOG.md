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

`total_sequences`, `mean_identity`, `median_identity`, `perfect_matches` and `total_mismatches` are unchanged.
Parameter rankings from `scripts/optimization/grid_search.py`, which weights `coverage` and `N50`, may change.
