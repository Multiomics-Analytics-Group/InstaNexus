#!/usr/bin/env python

r"""Grid Search Optimization for InstaNexus Assembly.

 ██████████   ███████████ █████  █████
░░███░░░░███ ░█░░░███░░░█░░███  ░░███
 ░███   ░░███░   ░███  ░  ░███   ░███
 ░███    ░███    ░███     ░███   ░███
 ░███    ░███    ░███     ░███   ░███
 ░███    ███     ░███     ░███   ░███
 ██████████      █████    ░░████████
░░░░░░░░░░      ░░░░░      ░░░░░░░░

__authors__ = Marco Reverenna
__copyright__ = Copyright 2025-2026
__research-group__ = DTU Biosustain (Multi-omics Network Analytics) and DTU Bioengineering
__date__ = 09 Dec 2025
__maintainer__ = Marco Reverenna
__email__ = marcor@dtu.dk
__status__ = Dev

This script performs hyperparameter optimization for the assembly module.
It runs in parallel, aggregates results, and calculates a normalized
Composite Score based on Coverage, N50, Identity, and Fragmentation.

Usage:
    python scripts/grid_search.py \
        --input-csv inputs/ma1_cleaned.csv \
        --metadata-json json/sample_metadata.json \
        --grid-json json/gridsearch_params.json \
        --mode dbg_weighted \
        --chain light \
        --workers 16
"""

import argparse
import itertools
import json
import logging
import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
from typing import Any, Dict, List

import pandas as pd
from sklearn.preprocessing import MinMaxScaler
from tqdm import tqdm

# Ensure project root is in path
SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = SCRIPT_DIR.parent
if str(PROJECT_ROOT) not in sys.path:
    sys.path.append(str(PROJECT_ROOT))

# Import InstaNexus modules
try:
    from instanexus import helpers, preprocessing, visualization
    from instanexus.assembly import Assembler
except ImportError as e:
    sys.exit(f"Error importing InstaNexus modules: {e}. Run from project root.")

# Setup Logging
logging.basicConfig(
    level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s", handlers=[logging.StreamHandler(sys.stdout)]
)
logger = logging.getLogger(__name__)


def load_grid_params(json_path: Path, mode: str) -> List[Dict[str, Any]]:
    """Loads and generates parameter combinations from JSON."""
    with open(json_path, "r") as f:
        all_grids = json.load(f)

    if mode not in all_grids:
        raise ValueError(f"Mode '{mode}' not found in {json_path}")

    grid = all_grids[mode]
    keys, values = zip(*grid.items(), strict=False)
    combinations = [dict(zip(keys, v, strict=False)) for v in itertools.product(*values)]
    return combinations


def compute_final_ranking(df_results: pd.DataFrame) -> pd.DataFrame:
    """
    Applies MinMax scaling to normalize metrics and computes the Composite Score.

    WEIGHTING STRATEGY: 'Aggressive Consolidation Split'
    ----------------------------------------------------
    1. Coverage (0.35): DOMINANT.
       Rationale: The primary goal is to recover the protein sequence. High N50
       is useless if we only recover 10% of the target.

    2. N50 (0.25) & Scaffold Count (0.25): STRUCTURAL (50% Combined).
       Rationale: We strongly penalize fragmentation. We want the algorithm to
       prioritize merging contigs into longer, fewer scaffolds over keeping
       them separate to maximize local identity.

    3. Mean Identity (0.15): QUALITY.
       Rationale: Lower weight because input data is usually pre-filtered
       (e.g., >80% identity during mapping). Differences between 95% and 99%
       are less critical than differences in coverage or fragmentation.
    """
    if df_results.empty:
        return df_results

    # Metrics to use for scoring
    metrics = ["coverage", "N50", "mean_identity", "scaffolds_count"]

    # Check if metrics exist
    available_metrics = [m for m in metrics if m in df_results.columns]
    if len(available_metrics) != len(metrics):
        logger.warning("Some metrics missing from results. Skipping scoring.")
        return df_results

    df_scoring = df_results[metrics].copy().fillna(0)

    # Normalize (0 to 1)
    scaler = MinMaxScaler()
    scaled_values = scaler.fit_transform(df_scoring)
    df_scaled = pd.DataFrame(scaled_values, columns=metrics, index=df_results.index)

    # Invert 'scaffolds_count' because fewer is better.
    # Formula: 1 - normalized_val (where 1 becomes best/fewest, 0 becomes worst/most)
    df_scaled["scaffolds_count"] = 1 - df_scaled["scaffolds_count"]

    # Define Weights
    weights = {"coverage": 0.35, "N50": 0.25, "scaffolds_count": 0.25, "mean_identity": 0.15}

    # Calculate Weighted Sum
    composite_scores = df_scaled[list(weights.keys())].dot(pd.Series(weights))

    # Merge back
    df_final = df_results.copy()
    df_final["composite_score"] = composite_scores

    # Sort by score descending
    return df_final.sort_values(by="composite_score", ascending=False)


def evaluate_combination(
    params: Dict[str, Any],
    df_original: pd.DataFrame,
    protein_norm: str,
    assembly_mode: str,
    mapping_cfg: Dict[str, Any],
) -> Dict[str, Any]:
    """
    Worker function:
    1. Filters input data based on dynamic FDR parameter only.
    2. Runs assembly on the filtered subset.
    3. Computes statistics.
    """
    try:
        df_subset = df_original.copy()

        if "fdr" in params:
            if "psm_q_value" in df_subset.columns:
                df_subset = df_subset[df_subset["psm_q_value"] <= params["fdr"]]
            else:
                return {**params, "error": "FDR requested but 'psm_q_value' column missing"}

        if "cleaned_preds" in df_subset.columns:
            sequences = df_subset["cleaned_preds"].dropna().tolist()
        else:
            return {**params, "error": "'cleaned_preds' column missing"}

        if not sequences:
            return {**params, "scaffolds_count": 0, "coverage": 0, "error": "No sequences after filtering"}

        # Exclude 'fdr' from params passed to Assembler class init
        assembler_params = {k: v for k, v in params.items() if k != "fdr"}

        assembler = Assembler(mode=assembly_mode, **assembler_params)

        start_time = time.time()
        # Note: We pass df_subset (filtered) because multimodal/weighted modes use it
        scaffolds = assembler.run(sequences=sequences, df_full=df_subset)
        duration = time.time() - start_time

        if not scaffolds:
            return {**params, "scaffolds_count": 0, "coverage": 0, "error": "No scaffolds produced"}

        mapped_scaffolds = visualization.process_protein_contigs_scaffold(
            assembled_contigs=scaffolds,
            target_protein=protein_norm,
            max_mismatches=mapping_cfg["max_mismatches"],
            min_identity=mapping_cfg["min_identity"],
        )

        if not mapped_scaffolds:
            return {**params, "scaffolds_count": len(scaffolds), "coverage": 0, "error": "Mapping failed"}

        df_mapped = visualization.create_dataframe_from_mapped_sequences(mapped_scaffolds)

        stats = helpers.compute_assembly_statistics(
            df=df_mapped,
            sequence_type="scaffolds",
            output_folder="",  # We don't save individual JSONs to save IO/Time
            reference=protein_norm,
        )

        return {
            **params,
            "scaffolds_count": len(scaffolds),
            "coverage": stats.get("coverage", 0),
            "N50": stats.get("N50", 0),
            "mean_identity": stats.get("mean_identity", 0),
            "total_mismatches": stats.get("total_mismatches", 0),
            "duration_sec": round(duration, 2),
            "input_sequences": len(sequences),
            "error": None,
        }

    except Exception as e:
        return {**params, "error": str(e)}


def main():
    parser = argparse.ArgumentParser(description="Hyperparameter Grid Search for InstaNexus.")
    parser.add_argument("--input-csv", required=True, help="Path to input cleaned CSV.")
    parser.add_argument("--metadata-json", required=True, help="Path to sample_metadata.json.")
    parser.add_argument("--grid-json", required=True, help="Path to gridsearch_params.json.")
    parser.add_argument("--mode", required=True, help="Assembly mode (greedy, dbg_weighted, etc.)")
    parser.add_argument("--output-dir", default="outputs/_grid_search", help="Directory to save results.")
    parser.add_argument("--chain", default="", help="Chain type (light/heavy).")
    parser.add_argument("--workers", type=int, default=8, help="Number of parallel workers.")

    parser.add_argument("--eval-min-identity", type=float, default=0.8, help="Min identity for evaluation mapping.")
    parser.add_argument("--eval-max-mismatches", type=int, default=100, help="Max mismatches for evaluation mapping.")

    args = parser.parse_args()

    input_path = Path(args.input_csv)
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    logger.info(f"Loading data from {input_path}...")
    run_name = input_path.stem.replace("_cleaned", "")

    df = pd.read_csv(input_path)

    try:
        meta = helpers.get_sample_metadata(run_name, chain=args.chain, json_path=args.metadata_json)
        protein_norm = preprocessing.normalize_sequence(meta["protein"])
        logger.info(f"Reference protein loaded for {run_name} (Chain: {args.chain}). Length: {len(protein_norm)}")
    except Exception as e:
        logger.error(f"Metadata error: {e}")
        return

    try:
        combinations = load_grid_params(Path(args.grid_json), args.mode)
        logger.info(f"Loaded {len(combinations)} combinations for mode '{args.mode}'.")
    except Exception as e:
        logger.error(f"Failed to load grid params: {e}")
        return

    results = []
    mapping_cfg = {"min_identity": args.eval_min_identity, "max_mismatches": args.eval_max_mismatches}

    logger.info(f"Starting execution with {args.workers} workers...")

    with ProcessPoolExecutor(max_workers=args.workers) as executor:
        futures = {
            executor.submit(evaluate_combination, params, df, protein_norm, args.mode, mapping_cfg): params
            for params in combinations
        }

        for future in tqdm(as_completed(futures), total=len(futures), desc="Evaluating"):
            res = future.result()
            results.append(res)

    df_results = pd.DataFrame(results)

    if not df_results.empty and "coverage" in df_results.columns:
        df_valid = df_results[df_results["error"].isnull()].copy()

        if not df_valid.empty:
            logger.info("Computing final ranking (Aggressive Consolidation Split)...")
            df_ranked = compute_final_ranking(df_valid)

            df_errors = df_results[df_results["error"].notnull()]
            df_final = pd.concat([df_ranked, df_errors])
            df_final["chain"] = args.chain if args.chain else "N/A"
            chain_suffix = f"_{args.chain}" if args.chain else ""
            timestamp = time.strftime("%Y%m%d-%H%M%S")

            csv_out = output_dir / f"grid_{args.mode}_{run_name}{chain_suffix}_{timestamp}.csv"
            df_final.to_csv(csv_out, index=False)

            logger.info(f"Results saved to: {csv_out}")

            best = df_ranked.iloc[0]
            logger.info(f" BEST CONFIG | Score: {best['composite_score']:.3f}")
            logger.info(
                f"   Cov: {best['coverage']*100:.1f}% | N50: {best['N50']} | Scaffolds: {best['scaffolds_count']}"
            )

            best_params = {k: best[k] for k in combinations[0].keys() if k in best}
            logger.info(f"   Params: {json.dumps(best_params, indent=2)}")
        else:
            logger.warning("All runs failed or produced no valid assemblies.")
            csv_out = output_dir / f"grid_{args.mode}_{run_name}_FAILED.csv"
            df_results.to_csv(csv_out, index=False)
    else:
        logger.error("No results produced.")


if __name__ == "__main__":
    main()
