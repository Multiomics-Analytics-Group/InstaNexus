#!/usr/bin/env python

""" Calculates and saves statistics.
 _____  _______  _    _ 
|  __ \|__   __|| |  | |
| |  | |  | |   | |  | |
| |  | |  | |   | |  | |
| |__| |  | |   | |__| |
|_____/   |_|   |______|

__authors__ = Marco Reverenna
__copyright__ = Copyright 2024-2025
__reserach-group__ = DTU Biosustain (Multi-omics Network Analytics) and DTU Bioengineering
__date__ = 26 Jun 2024
__maintainer__ = Marco Reverenna
__email__ = marcor@dtu.dk
__status__ = Dev
"""

import os
import json
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt



def compute_assembly_statistics(df, sequence_type, output_folder, reference, **params):
    """ Statistics for contigs and scaffolds

    Args:
        df: DataFrame with mapped values
        sequence_type: either 'contigs' or 'scaffold'
        reference: reference protein normalized
    """

    statistics = {}
    statistics.update(params)  # add the hyperparameters to the statistics

    df['sequence_length'] = df['end'] - df['start'] + 1

    # Reference coordinates
    statistics['reference_start'] = 0
    statistics['reference_end'] = len(reference) + 1

    # Sequences statistics
    statistics['total_sequences'] = len(df)
    statistics['average_length'] = df['sequence_length'].mean()
    statistics['min_length'] = df['sequence_length'].min()
    statistics['max_length'] = df['sequence_length'].max()

    # Create a set of covered positions (adjusting for 0-based indexing)
    covered_positions = set()
    for start, end in zip(df["start"], df["end"]):
        covered_positions.update(range(start - 1, end))  # Convert 1-based to 0-based indexing
    statistics['coverage'] = len(covered_positions) / statistics['reference_end']

    # Identity score statistics
    statistics['mean_identity'] = df['identity_score'].mean()
    statistics['median_identity'] = df['identity_score'].median()
    #statistics['std_identity'] = df['identity_score'].std()

    # Mismatch statistics
    statistics['perfect_matches'] = sum(df['mismatches_pos'].apply(len) == 0)  # sequences with no mismatches
    all_mismatches = [pos for mismatches in df['mismatches_pos'] for pos in mismatches]
    statistics['total_mismatches'] = len(set(all_mismatches))

    # N50 and N90 calculations
    lengths = sorted(df['sequence_length'], reverse=True)
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
        # Break early if both values are computed
        if n50 is not None and n90 is not None:
            break

    statistics['N50'] = n50
    statistics['N90'] = n90

    # Save JSON file
    file_name = f"{sequence_type}_stats.json"
    output_path = os.path.join(output_folder, file_name)
    with open(output_path, "w") as file:
        json.dump(statistics, file, indent=4)

    return statistics