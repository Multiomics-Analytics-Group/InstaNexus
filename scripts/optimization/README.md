# Optimization

Parallel hyperparameter optimization for the InstaNexus assembly module.

This folder sweeps assembly parameters across a grid, evaluates each combination against a
reference, and ranks them with a normalized **Composite Score** (Reverenna et al., bioRxiv
2025) combining Coverage, N50, scaffold count, and maximum contig length.

## Contents

| File | Description |
|---|---|
| `grid_search.py` | Core grid-search optimizer (also exposed as the `instanexus-optimize` CLI) |
| `analyze_optimization.py` | Aggregates results into summary tables and SVG heatmaps |
| `run_all_gridsearch.sh` | Batch-runs the grid search across all samples and assembly modes |
| `run_preprocessing.sh` | Batch-preprocesses raw input CSVs into `inputs/cleaned/` |
| `debug.py` | Minimal smoke test for imports and a single assembly run |

The parameter grid for each assembly mode (`greedy`, `dbg_weighted`, `multimodal_dbg`) is
defined in [`json/gridsearch_params.json`](../../json/gridsearch_params.json).

## Running a grid search

Use the installed CLI:

```bash
instanexus-optimize \
    --input-csv inputs/ma1_cleaned.csv \
    --metadata-json json/sample_metadata.json \
    --grid-json json/gridsearch_params.json \
    --mode dbg_weighted \
    --chain light \
    --workers 16
```

or run the script directly:

```bash
python scripts/optimization/grid_search.py \
    --input-csv inputs/ma1_cleaned.csv \
    --metadata-json json/sample_metadata.json \
    --grid-json json/gridsearch_params.json \
    --mode dbg_weighted \
    --chain light \
    --workers 16
```

### Arguments

| Flag | Default | Description |
|---|---|---|
| `--input-csv` | *(required)* | Raw or cleaned input CSV. Preprocessing runs automatically if no cleaned file exists in `--output-dir`. |
| `--metadata-json` | *(required)* | Path to `sample_metadata.json` (required for reference protein lookup). |
| `--grid-json` | *(required)* | Path to `gridsearch_params.json` defining the parameter grid. |
| `--mode` | *(required)* | Assembly mode (`greedy`, `dbg_weighted`, `multimodal_dbg`). |
| `--chain` | `""` | Chain type for antibodies (`light` / `heavy`); omit for single-chain samples. |
| `--workers` | `8` | Number of parallel worker processes. |
| `--output-dir` | `outputs/_grid_search` | Directory to save results. |
| `--eval-min-identity` | `0.8` | Minimum identity for evaluation mapping. |
| `--eval-max-mismatches` | `100` | Maximum mismatches for evaluation mapping. |

The preprocessing step (used when no cleaned CSV is present) also accepts optional
`--conf`, `--fdr`, and `--contaminants-fasta` flags.

## Typical workflow

1. **Preprocess** raw inputs into `inputs/cleaned/` (antibody samples get `--chain light`):

   ```bash
   bash scripts/optimization/run_preprocessing.sh
   ```

2. **Sweep** all samples across every assembly mode. Edit the `MODES`, `ANTIBODIES`,
   `OTHERS`, and `WORKERS` variables at the top of the script to match your data, then:

   ```bash
   bash scripts/optimization/run_all_gridsearch.sh
   ```

   Results land in `outputs/_grid_search/<mode>/` and per-run logs in `logs/`.

3. **Analyze** the results into publication-ready tables and heatmaps:

   ```bash
   python scripts/optimization/analyze_optimization.py
   ```

   Outputs are written to `outputs/_summary_tables/` and `outputs/_optimization_figures/`.

## Composite Score

Each parameter combination is scored by min-max normalizing the metrics across the grid
and combining them with fixed weights (Reverenna et al., bioRxiv 2025):

```
composite_score = 0.5 * coverage_norm
                + 0.3 * N50_norm
                + 0.1 * (1 - scaffolds_count_norm)   # inverted: fewer scaffolds = better
                + 0.1 * max_length_norm
```

Higher coverage, longer N50, fewer scaffolds, and longer maximum contigs all increase the
score. `mean_identity` is collected for reporting but is **not** part of this formula. The
top-ranked combination is the recommended parameter set for that sample and mode, and the
optimizer prints a ready-to-run `instanexus` command for it.
