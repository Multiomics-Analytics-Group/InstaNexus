#!/usr/bin/env python

r"""
 _____  _______  _    _ 
|  __ \|__   __|| |  | |
| |  | |  | |   | |  | |
| |  | |  | |   | |  | |
| |__| |  | |   | |__| |
|_____/   |_|   |______|

__authors__ = Marco Reverenna & Konstantinos Kalogeropoulus
__copyright__ = Copyright 2024-2025
__research-group__ = DTU Biosustain (Multi-omics Network Analytics) and DTU Bioengineering
__date__ = 26 Nov 2024
__maintainer__ = Marco Reverenna
__email__ = marcor@dtu.dk
__status__ = Dev
"""

# import libraries
import os
import re
import numpy as np
import pandas as pd
from Bio import SeqIO
import plotly.express as px
import plotly.graph_objects as go


# Define and create the necessary directories only if they don't exist
def create_directory(path):
    """Creates a directory if it does not already exist.
    Args:
        path (str): The path of the directory to create.
    """
    if not os.path.exists(path):
        os.makedirs(path)
        #print(f"Created: {path}")
    #else:
        #print(f"Already exists: {path}")

def create_subdirectories_outputs(folder):
    """Creates subdirectories within the specified folder.
    Args:
        folder (str): The path of the parent directory.
    """
    subdirectories = ['preprocessing', 'contigs', 'scaffolds', 'consensus', 'logs', 'statistics']
    for subdirectory in subdirectories:
        create_directory(f"{folder}/{subdirectory}")

def create_subdirectories_figures(folder):
    """Creates subdirectories within the specified folder.
    Args:
        folder (str): The path of the parent directory.
    """
    subdirectories = ['preprocessing', 'contigs', 'scaffolds', 'consensus', 'heatmap', 'logo']
    for subdirectory in subdirectories:
        create_directory(f"{folder}/{subdirectory}")


def normalize_sequence(sequence):
    """
    Normalize the given amino acid sequence by replacing all occurrences of 'I' with 'L'.

    Parameters:
    sequence (str): The amino acid sequence to be normalized.

    Returns:
    str: The normalized amino acid sequence with 'I' replaced by 'L'.
    """
    return sequence.replace('I', 'L')



def remove_modifications(psm_column):
    """
    Remove any content within parentheses, including the parentheses, from a given string.

    Parameters:
    - psm_column (str): The string containing modifications in parentheses (e.g., "A(ox)BC(mod)D"). If the value is null, it returns None.

    Returns:
    - str: The string with all parenthetical modifications removed (e.g., "ABCD"), or None if the input was null.
    """

    if pd.notnull(psm_column):
        return re.sub(r'\(.*?\)', '', psm_column)  # Replace any content in parentheses with an empty string
    return None

def clean_dataframe(df):
    """
    Clean and preprocess a DataFrame for analysis by removing '(ox)' substrings from sequences in the 'seq' column.
    by replacing values of -1 with -10 in the 'log_probs' column, by dropping rows with missing values in the 'preds' column.
    by extracts a 'protease' value from the 'experiment_name' column based on a specific naming convention.
    by adding a 'conf' column, which is the exponentiated 'log_probs' to represent confidence and sorting 
    the DataFrame by the 'conf' column in descending order.

    Parameters:
    - df (DataFrame): The raw input DataFrame to clean.

    Returns:
    - DataFrame: The cleaned and processed DataFrame.
    """
    # update 1
    df = df.copy()

    df['log_probs'] = df['log_probs'].replace(-1, -10)
    # -10 is very low, replacing with -1 (so that we are sure is very low quality prediction)
    # check new update InstaNovo 

    df = df.dropna(subset=['preds'])

    #df['protease'] = df['experiment_name'].apply(lambda x: x.split('_')[-3] if isinstance(x, str) else None)
    #df.loc[:, 'protease'] = df['experiment_name'].apply(lambda x: x.split('_')[-3] if isinstance(x, str) else None)
    
    #df['conf'] = np.exp(df['log_probs'])
    
    df.loc[:, 'conf'] = np.exp(df['log_probs'])

    df = df.sort_values('conf', ascending=False)

    return df



def filter_contaminants(seqs, run, contaminants_fasta):
    """
    Filters out sequences from the input list `seqs` that are substrings of sequences
    in the contaminants file. If run == 'bsa', the Bovine serum albumin precursor is ignored.

    Parameters:
    - seqs (list of str): List of sequences to be filtered.
    - contaminants_fasta (str): Path to the FASTA file containing contaminant sequences.
    - run (str): Run identifier, used to control special filtering logic.
    """

    contam_records = []
    for record in SeqIO.parse(contaminants_fasta, "fasta"):
        if run == 'bsa' and "Bovine serum albumin precursor" in record.description:
            continue  # Skip BSA if run is 'bsa'
        contam_records.append(str(record.seq))

    filtered_seqs = []
    removed_count = 0

    for seq in seqs:
        if any(seq in contam_seq for contam_seq in contam_records):
            removed_count += 1
        else:
            filtered_seqs.append(seq)

    #print(f"Removed {removed_count} contaminant sequences, {len(filtered_seqs)} sequences remaining.")
    return filtered_seqs


def plot_confidence_distribution(df, folder_figures, min_conf=0, max_conf=1):
    """
    Plots the distribution of confidence scores from a DataFrame.

    Parameters:
    df (pandas.DataFrame): DataFrame containing a column 'conf' with confidence scores.
    min_conf (float, optional): Minimum value of confidence range (default is 0).
    max_conf (float, optional): Maximum value of confidence range (default is 1).
    """
    # Filter the data based on the specified range
    filtered_df = df[(df['conf'] >= min_conf) & (df['conf'] <= max_conf)]

    fig = go.Figure()

    fig.add_trace(go.Histogram(
        x=filtered_df['conf'],
        xbins=dict(start=min_conf, end=max_conf, size=(max_conf - min_conf) / 40),
        marker=dict(color='brown'),
        opacity=1  # Remove opacity
    ))

    # Configure layout
    fig.update_layout(title='Confidence score distribution between {} and {}'.format(min_conf, max_conf),
        xaxis_title='Values', yaxis_title='Frequency',
        bargap=0.1, height=700, width=1200,
        margin=dict(l=50, r=50, t=100, b=100)
    )

    # Configure axes
    fig.update_xaxes(showgrid=True, gridcolor='lightgray',
        ticklabelposition='outside bottom',  # Place labels outside the bottom of the axis
        dtick=0.02  # Set the distance between ticks
    )
    fig.update_yaxes(showgrid=True, gridcolor='lightgray')
    fig.write_image(f"{folder_figures}/confidence_distribution_range_{min_conf}_{max_conf}.png")
    #fig.show()


def extract_protease(experiment_name, proteases):
    """ Extracts the protease name from the given experiment name.
    Parameters:
        experiment_name (str): The name of the experiment.
        proteases (list or set): A list or set of known protease names.

    Returns:
        str or None: The matched protease name, or None if no match is found.
    """
    parts = experiment_name.split('_')
    for part in parts:
        if part in proteases:
            return part
    return None



def plot_protease_distribution(protease_counts, folder_figures):
    """Creates an interactive bar plot of protease distribution using Plotly.
    
    Parameters:
        protease_counts (pandas.Series): A Pandas Series with protease names as the index
                                         and their counts as the values.
    """
    # Convert the Series to a DataFrame for compatibility with Plotly
    protease_df = protease_counts.reset_index()
    protease_df.columns = ['Protease', 'Count']

    fig = px.bar(
        protease_df,
        x='Protease',
        y='Count',
        title='Proteases distribution',
        labels={'Protease': 'Proteases', 'Count': 'Counts'},
        text='Count'
    )

    # Reduce the width of the bars and set the text position outside the bars
    fig.update_traces(textposition='outside', width=0.4)

    # Define conversion factor from mm to pixels (approx. 3.78 px/mm at 96 DPI)
    mm_to_px = 3.78

    # Set the desired figure dimensions according to Nature's standards.
    # For a double-column figure:
    width_mm = 240   # width in mm (double column) 183 per Nature
    height_mm = 200  # full page depth in mm 247 per Nature

    fig.update_layout(
        width=int(width_mm * mm_to_px),
        height=int(height_mm * mm_to_px),
        xaxis_title='Proteases',
        yaxis_title='Counts',
        xaxis_tickangle=0,
        showlegend=False,
        title_font_size=12,
        font=dict(
            family="Arial, sans-serif",  # Change to your preferred font
            size=8,
            color="black"
        ),
        margin=dict(t=50, b=50, l=50, r=100),  # Adjust margins as needed
        plot_bgcolor='white',   # White background for the plot area
        paper_bgcolor='white'   # White background for the entire figure
    )

    # Export the figure in SVG vector format
    fig.write_image(f"{folder_figures}/proteases_distribution.svg")
    fig.show()
