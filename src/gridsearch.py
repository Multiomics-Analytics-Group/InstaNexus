#!/usr/bin/env python

r""" Hyperparameter optimization script for assembly analysis.
 _____  _______  _    _
|  __ \|__   __|| |  | |
| |  | |  | |   | |  | |
| |  | |  | |   | |  | |
| |__| |  | |   | |__| |
|_____/   |_|   |______|

__authors__ = Marco Reverenna
__copyright__ = Copyright 2025-2026
__research-group__ = DTU Biosustain (Multi-omics Network Analytics) and DTU Bioengineering
__date__ = 22 Jun 2025
__maintainer__ = Marco Reverenna
__email__ = marcor@dtu.dk
__status__ = Dev
"""

# import libraries
import os
import logging
import itertools
from tqdm import tqdm
from complete_dbg import run_pipeline_dbg
#from complete_greedy import run_pipeline_greedy
from concurrent.futures import ProcessPoolExecutor, as_completed

# Define the parameter grid and set values to test

parameter_grid_dbg = {
       "kmer_size": [6, 7],
       "min_overlap": [3, 4],
       "size_threshold": [0, 5, 10],
       "max_mismatches": [8, 10, 12, 14],
       "min_identity": [0.6, 0.7, 0.8, 0.9],
       "conf": [0.86, 0.88, 0.90, 0.92]
       }


# parameter_grid_greedy = {
#       "min_overlap": [3, 4],
#       "size_threshold": [0, 5, 10],
#       "max_mismatches": [8, 10, 12, 14], 
#       "min_identity": [0.6, 0.7, 0.8, 0.9],
#       "conf": [0.86, 0.88, 0.90, 0.92]
#       }


keys, values = zip(*parameter_grid_dbg.items())
#keys, values = zip(*parameter_grid_dbg.items())

combinations = [dict(zip(keys, v)) for v in itertools.product(*values)]
total_combinations = len(combinations)

os.makedirs("logs", exist_ok=True)

# Set up logging
handlers = [
    logging.FileHandler(f"logs/grid_search.log"),
    logging.StreamHandler()
]
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(levelname)s - %(message)s",
    handlers=handlers
)

logging.info(f"Starting hyperparameter optimization with {total_combinations} combinations.")
print(f"Total combinations: {total_combinations}")

def run_analysis(params, iteration):
    """Wrapper function to run the main analysis with error handling."""
    try:
        logging.info(f"[ITER {iteration}] Starting with parameters: {params}")
        run_pipeline_dbg(**params)
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