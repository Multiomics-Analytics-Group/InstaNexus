#!/usr/bin/env python

r"""

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
__date__ = 01 Nov 2025
__maintainer__ = Marco Reverenna
__email__ = marcor@dtu.dk
__status__ = Dev
"""

# import libraries
import numpy as np
import pandas as pd
import os
import json
import plotly.express as px
import plotly.graph_objects as go
from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[2]
JSON_DIR = PROJECT_ROOT / "json"


def get_sample_metadata(run, chain="", json_path=JSON_DIR / "sample_metadata.json"):
    with open(json_path, "r") as f:
        all_meta = json.load(f)

    if run not in all_meta:
        raise ValueError(f"Run '{run}' not found in metadata.")

    entries = all_meta[run]

    if not chain:
        # If no chain is specified, return the first entry
        return entries[0]

    for entry in entries:
        if entry["chain"] == chain:
            return entry

    raise ValueError(f"No metadata found for run '{run}' with chain '{chain}'.")


# Define and create the necessary directories only if they don't exist
def create_directory(path):
    """Creates a directory if it does not already exist.
    Args:
        path (str): The path of the directory to create.
    """
    if not os.path.exists(path):
        os.makedirs(path)
        # print(f"Created: {path}")
    # else:
    # print(f"Already exists: {path}")


def create_subdirectories_outputs(folder):
    """Creates subdirectories within the specified folder.
    Args:
        folder (str): The path of the parent directory.
    """
    subdirectories = ["cleaned", "contigs", "scaffolds", "statistics"]
    for subdirectory in subdirectories:
        create_directory(f"{folder}/{subdirectory}")


def create_subdirectories_figures(folder):
    """Creates subdirectories within the specified folder.
    Args:
        folder (str): The path of the parent directory.
    """
    subdirectories = [
        "preprocessing",
        "contigs",
        "scaffolds",
        "consensus",
        "heatmap",
        "logo",
    ]
    for subdirectory in subdirectories:
        create_directory(f"{folder}/{subdirectory}")


def compute_assembly_statistics(df, sequence_type, output_folder, reference, **params):
    """Statistics for contigs and scaffolds

    Args:
        df: DataFrame with mapped values
        sequence_type: either 'contigs' or 'scaffold'
        output_folder: folder to save output
        reference: reference protein normalized
    """

    statistics = {}
    statistics.update(params)  # add the hyperparameters to the statistics

    df["sequence_length"] = df["end"] - df["start"] + 1

    # Reference coordinates
    statistics["reference_start"] = int(0)
    statistics["reference_end"] = int(len(reference) + 1)

    # Sequences statistics
    statistics["total_sequences"] = int(len(df))
    statistics["average_length"] = float(df["sequence_length"].mean())
    statistics["min_length"] = int(df["sequence_length"].min())
    statistics["max_length"] = int(df["sequence_length"].max())

    # Create a set of covered positions (adjusting for 0-based indexing)
    covered_positions = set()
    for start, end in zip(df["start"], df["end"]):
        covered_positions.update(
            range(start - 1, end)
        )  # Convert 1-based to 0-based indexing
    statistics["coverage"] = float(len(covered_positions) / statistics["reference_end"])

    # identity score statistics
    statistics["mean_identity"] = float(df["identity_score"].mean())
    statistics["median_identity"] = float(df["identity_score"].median())
    # statistics['std_identity'] = float(df['identity_score'].std())

    # mismatch statistics
    statistics["perfect_matches"] = int(
        sum(df["mismatches_pos"].apply(len) == 0)
    )  # sequences with no mismatches
    all_mismatches = [pos for mismatches in df["mismatches_pos"] for pos in mismatches]
    statistics["total_mismatches"] = int(len(set(all_mismatches)))

    # N50 and N90 calculations
    lengths = sorted(df["sequence_length"], reverse=True)
    total_length = sum(lengths)

    cumulative_length = 0
    n50 = None
    n90 = None
    for length in lengths:
        cumulative_length += length
        if n50 is None and cumulative_length >= total_length * 0.5:
            n50 = length
        if n90 is None and cumulative_length >= total_length * 0.9:
            n90 = length
        if n50 is not None and n90 is not None:
            break

    statistics["N50"] = int(n50)
    statistics["N90"] = int(n90)

    file_name = f"{sequence_type}_stats.json"
    output_path = os.path.join(output_folder, file_name)
    with open(output_path, "w") as file:
        json.dump(statistics, file, indent=4)

    return statistics