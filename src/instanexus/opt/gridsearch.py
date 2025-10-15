#!/usr/bin/env python

r"""Hyperparameter optimization script for assembly analysis.
 _____  _______  _    _
|  __ \|__   __|| |  | |
| |  | |  | |   | |  | |
| |  | |  | |   | |  | |
| |__| |  | |   | |__| |
|_____/   |_|   |______|

__authors__ = Marco Reverenna and Pasquale D. Colaianni
__copyright__ = Copyright 2025-2026
__research-group__ = DTU Biosustain (Multi-omics Network Analytics) and DTU Bioengineering
__date__ = 06 Aug 2025
__maintainer__ = Marco Reverenna
__email__ = marcor@dtu.dk
__status__ = Dev
"""

import itertools
import json
import logging

# import libraries
import os
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

from tqdm import tqdm

from src.opt.opt_dbg import run_pipeline_dbg
from src.opt.opt_greedy import run_pipeline_greedy

BASE_DIR = Path(__file__).resolve().parents[2]


# Define the parameter grid and set values to test

# with open("../../json/gridsearch_params.json") as f:
#     all_grids = json.load(f)

with open(BASE_DIR / "json" / "gridsearch_params.json") as f:
    all_grids = json.load(f)

method = "greedy"  # Change to "greedy" for greedy method

selected_grid = all_grids[method]

keys, values = zip(*selected_grid.items())

combinations = [dict(zip(keys, v)) for v in itertools.product(*values)]
total_combinations = len(combinations)

os.makedirs("logs", exist_ok=True)

# Set up logging
handlers = [logging.FileHandler("logs/grid_search.log"), logging.StreamHandler()]
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(levelname)s - %(message)s",
    handlers=handlers,
)

logging.info(
    f"Starting hyperparameter optimization with {total_combinations} combinations."
)
print(f"Total combinations: {total_combinations}")


def run_analysis(params, iteration):
    """Wrapper function to run the main analysis with error handling."""
    try:
        logging.info(f"[ITER {iteration}] Starting with parameters: {params}")
        logging.info(f"Grid search started using method '{method}'.")
        if method == "dbg":
            run_pipeline_dbg(**params)
        elif method == "greedy":
            run_pipeline_greedy(**params)
        logging.info(f"[ITER {iteration}] Completed successfully.")
    except Exception as e:
        logging.error(f"[ITER {iteration}] Failed with parameters {params}: {str(e)}")


def grid_search_parallel():
    """Perform hyperparameter optimization in parallel."""
    with ProcessPoolExecutor(max_workers=64) as executor:
        futures = {
            executor.submit(run_analysis, params, idx + 1): idx + 1
            for idx, params in enumerate(combinations)
        }

        for _ in tqdm(as_completed(futures), total=len(futures), desc="Processing"):
            pass

    logging.info("Hyperparameter optimization completed.")


if __name__ == "__main__":
    grid_search_parallel()
