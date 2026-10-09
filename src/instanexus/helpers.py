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
import json
import os
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

    Reference positions are 0-based and end-exclusive, as returned by visualization.map_to_protein.
    Per-reference-position metrics count each position once, however many sequences cover it:

    - coverage: fraction of reference positions covered by at least one mapped sequence.
    - mismatched_positions: number of reference positions with a mismatch in at least one mapped
      sequence (start + offset of each mismatch).

    Per-sequence metrics add up over the mapped sequences:

    - perfect_matches: number of sequences without mismatches.
    - total_mismatches: total number of mismatches over all mapped sequences; a reference position
      covered by several sequences with the same error counts once per sequence.

    Args:
        df: DataFrame with mapped values
        sequence_type: either 'contigs' or 'scaffold'
        output_folder: folder to save output
        reference: reference protein normalized

    Returns:
        The statistics, also written to ``<output_folder>/<sequence_type>_stats.json``.
    """

    statistics = {}
    statistics.update(params)  # add the hyperparameters to the statistics

    # start is 0-based and end is exclusive, as returned by visualization.map_to_protein
    df["sequence_length"] = df["end"] - df["start"]

    # Reference coordinates (0-based, end exclusive)
    statistics["reference_start"] = int(0)
    statistics["reference_end"] = int(len(reference))

    # Sequences statistics
    statistics["total_sequences"] = int(len(df))
    statistics["average_length"] = float(df["sequence_length"].mean())
    statistics["min_length"] = int(df["sequence_length"].min())
    statistics["max_length"] = int(df["sequence_length"].max())

    # Set of covered reference positions; range(start, end) because start is 0-based and end exclusive
    covered_positions = set()
    for start, end in zip(df["start"], df["end"], strict=False):
        covered_positions.update(range(start, end))
    statistics["coverage"] = float(len(covered_positions) / statistics["reference_end"])

    # identity score statistics
    statistics["mean_identity"] = float(df["identity_score"].mean())
    statistics["median_identity"] = float(df["identity_score"].median())
    # statistics['std_identity'] = float(df['identity_score'].std())

    # mismatch statistics; mismatches_pos holds 0-based offsets within each sequence
    statistics["perfect_matches"] = int(sum(df["mismatches_pos"].apply(len) == 0))  # sequences with no mismatches
    statistics["total_mismatches"] = int(df["mismatches_pos"].apply(len).sum())
    # same reference coordinates as coverage: start (0-based) + offset
    mismatched_positions = set()
    for start, mismatches in zip(df["start"], df["mismatches_pos"], strict=False):
        mismatched_positions.update(start + offset for offset in mismatches)
    statistics["mismatched_positions"] = int(len(mismatched_positions))

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

    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    output_path = os.path.join(output_folder, file_name)

    with open(output_path, "w") as file:
        json.dump(statistics, file, indent=4)

    return statistics
