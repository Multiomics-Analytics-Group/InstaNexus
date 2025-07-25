#!/usr/bin/env python

r""" 
 _____  _______  _    _ 
|  __ \|__   __|| |  | |
| |  | |  | |   | |  | |
| |  | |  | |   | |  | |
| |__| |  | |   | |__| |
|_____/   |_|   |______|

__authors__ = Marco Reverenna & Konstantinos Kalogeropoulus
__copyright__ = Copyright 2025-2026
__research-group__ = DTU Biosustain (Multi-omics Network Analytics) and DTU Bioengineering
__date__ = 21 Mar 2025
__maintainer__ = Marco Reverenna
__email__ = marcor@dtu.dk
__status__ = Dev
"""

import os
import json
import re
import statistics
import pandas as pd
import Bio.SeqIO
from collections import Counter
import matplotlib.pyplot as plt
import seaborn as sns
import logomaker
from tqdm import tqdm
from Bio import SeqIO
import plotly.express as px
import plotly.io as pio



def generate_pssm(aligned_records):
    """Generate a Position-Specific Scoring Matrix (PSSM) from aligned sequence records and calculates the frequency of each amino acid
    at each position across all sequences, resulting in a matrix where each position has the relative frequency of each amino acid.
    """
    pssm = {}
    
    for record in aligned_records:
        for i, aa in enumerate(record.seq):
            if i not in pssm:
                pssm[i] = Counter()
            if aa != '-':
                pssm[i][aa] += 1
    
    for i in pssm:
        total = sum(pssm[i].values())
        for aa in pssm[i]:
            pssm[i][aa] /= total
    
    pssm_df = pd.DataFrame(pssm).fillna(0).T
    pssm_df.index = pssm_df.index + 1
    pssm_df = pssm_df.sort_index(axis=1)
    
    return pssm_df



def generate_consensus(pssm_df, threshold=0.6):
    """ Generate a consensus sequence from a Position-Specific Scoring Matrix (PSSM) DataFrame. 
    It determine the most frequent amino acid at each position. If the highest frequency at a position exceeds the given threshold, that amino acid is added
    to the consensus sequence. Otherwise, a gap ('-') is added to indicate lack of consensus at that position.
    """
    
    consensus = ""
    
    for i in pssm_df.index:
        # Check if the highest frequency at this position is above the threshold
        if pssm_df.loc[i].max() > threshold:
            # If so, add the amino acid with the highest frequency to the consensus
            consensus += pssm_df.loc[i].idxmax()
        else:
            consensus += '-'
    
    return consensus



def plot_heatmap(pssm_df, output_file):
    """ Plots a heatmap of the given PSSM DataFrame and saves it to the specified output file.
    """
    # Set a maximum width for the figure to prevent excessively large images
    max_fig_width = 50  # Maximum width in inches
    fig_width = min(len(pssm_df) / 1.5, max_fig_width)

    plt.figure(figsize=(fig_width, 8))
    cmap = sns.cubehelix_palette(start=2, rot=0, dark=.15, light=.85, reverse=True, as_cmap=True)
    sns.heatmap(pssm_df.T, cmap=cmap, cbar=False, linewidths=0.1, linecolor='white')
    plt.yticks(rotation=0, fontsize=15)
    plt.xticks(fontsize=15)
    plt.tight_layout()
    plt.savefig(output_file)
    plt.close()
    #plt.show()



def plot_heatmap2(pssm_df, output_file):
    """
    Plots a heatmap of the given PSSM DataFrame using Plotly, styled in red tones,
    with Arial font, labeled axes and colorbar, and saves as high-resolution SVG.
    """
    df_t = pssm_df.T  # Transpose to have amino acids on Y and positions on X

    fig = px.imshow(
        df_t,
        color_continuous_scale='Reds',
        aspect='auto',
        labels={'x': 'Position', 'y': ' Amino acid', 'color': 'Frequency'}
    )

    fig.update_layout(
        font=dict(family="DejaVu Sans", size=15),
        xaxis_title="Position",
        yaxis_title=" Amino acid",
        margin=dict(l=40, r=40, t=40, b=40),
        coloraxis_colorbar=dict(title='Frequency')
    )

    # Show every 5th position on the x-axis
    x_tickvals = list(range(0, df_t.shape[1], 5))
    fig.update_xaxes(
        tickmode='array',
        tickvals=x_tickvals,
        tickangle=0
    )

    fig.update_yaxes(autorange="reversed")  # To match seaborn orientation
    pio.write_image(fig, output_file, format="svg", width=1200, height=400, scale=2)



def plot_logo(pssm_df, title, output_file):

    max_fig_width = 50  # Maximum width in inches
    fig_width = min(len(pssm_df) / 1.5, max_fig_width)

    _, ax = plt.subplots(1, 1, figsize=[fig_width, 3])
    logo = logomaker.Logo(
        pssm_df,
        ax=ax,
        font_name="DejaVu Sans",
        color_scheme="NajafabadiEtAl2017",
        stack_order="big_on_top",
        center_values=False,
        flip_below=False,
        fade_below=.5,
        shade_below=.5,
        fade_probabilities=False,
        vpad=0.05,
        vsep=0.0,
        width=0.85,
        baseline_width=0.5,
    )
    logo.style_xticks(anchor=0, rotation=0, spacing=1, fontsize=20, ha="center")
    plt.yticks([0, 0.5, 1], fontsize=20)
    logo.style_spines(spines=["left", "right", "bottom", "top"], visible=False)
    logo.ax.set_ylabel("Frequency", fontsize=20)
    logo.ax.set_xlabel("Position", fontsize=20)
    logo.ax.set_title(title, fontsize=20)
    plt.tight_layout()
    plt.savefig(output_file)
    #plt.show()
    plt.close()


def plot_logo2(pssm_df, output_file):
    """
    Plots a sequence logo from a PSSM DataFrame using Logomaker and saves it as a high-resolution SVG.
    """
    max_fig_width = 50  # Limit for very long sequences
    fig_width = min(len(pssm_df) / 1.5, max_fig_width)

    fig, ax = plt.subplots(figsize=[fig_width, 3])

    logo = logomaker.Logo(
        pssm_df,
        ax=ax,
        font_name="DejaVu Sans",
        color_scheme="NajafabadiEtAl2017",
        stack_order="big_on_top",
        center_values=False,
        flip_below=False,
        fade_below=0.5,
        shade_below=0.5,
        fade_probabilities=False,
        vpad=0.05,
        vsep=0.0,
        width=0.85,
        baseline_width=0.5,
    )

    # Style ticks and labels
    logo.style_xticks(anchor=0, rotation=0, spacing=1, fontsize=16, ha="center")
    ax.set_yticks([0, 0.5, 1])
    ax.set_yticklabels([0, 0.5, 1], fontsize=16)
    ax.set_ylabel("Frequency", fontsize=18)
    ax.set_xlabel("Position", fontsize=18)

    # Clean spines and remove title
    logo.style_spines(spines=["left", "right", "bottom", "top"], visible=False)

    # Tight layout and SVG save
    plt.tight_layout()
    plt.savefig(output_file, format="svg", dpi=300)
    plt.close()



def process_alignment_files(align_folder, consensus_folder):
    """ Process alignment files in the specified folder and generate consensus sequences, heatmaps, and logos"
    """
    for align_subfolder in os.listdir(align_folder):
        align_subfolder_path = os.path.join(align_folder, align_subfolder)
        print(f"Processing {align_subfolder_path}")
        if os.path.isdir(align_subfolder_path):
            # Create subfolders for consensus, heatmap, and logo
            consensus_subfolder = os.path.join(consensus_folder, align_subfolder)
            consensus_fasta_folder = os.path.join(consensus_subfolder, "consensus_fasta")
            heatmap_folder = os.path.join(consensus_subfolder, "heatmap")
            logo_folder = os.path.join(consensus_subfolder, "logo")
            os.makedirs(consensus_fasta_folder, exist_ok=True)
            os.makedirs(heatmap_folder, exist_ok=True)
            os.makedirs(logo_folder, exist_ok=True)

            for alignment_file in tqdm(os.listdir(align_subfolder_path)):
                if alignment_file.endswith('.afa'):
                    alignment_path = os.path.join(align_subfolder_path, alignment_file)
                    base_filename = os.path.splitext(alignment_file)[0]

                    aligned_records = list(Bio.SeqIO.parse(alignment_path, "fasta"))
                    pssm_df = generate_pssm(aligned_records)
                    consensus_sequence = generate_consensus(pssm_df)

                    # Write consensus sequence to fasta file
                    consensus_record = Bio.SeqRecord.SeqRecord(
                        Bio.Seq.Seq(consensus_sequence),
                        id=base_filename,
                        description="Consensus sequence"
                    )
                    consensus_fasta_path = os.path.join(consensus_fasta_folder, f"{base_filename}_consensus.fasta")
                    Bio.SeqIO.write([consensus_record], consensus_fasta_path, "fasta")

                    heatmap_path = os.path.join(heatmap_folder, f"{base_filename}_heatmap.svg")
                    plot_heatmap2(pssm_df, heatmap_path)

                    logo_path = os.path.join(logo_folder, f"{base_filename}_logo.svg")
                    #plot_logo(pssm_df, base_filename, logo_path)
                    plot_logo2(pssm_df, logo_path)


def generate_consensus_stats(consensus_base_folder):
    for cluster_folder in os.listdir(consensus_base_folder):
        cluster_path = os.path.join(consensus_base_folder, cluster_folder)

        if os.path.isdir(cluster_path):
            consensus_fasta_path = os.path.join(cluster_path, "consensus_fasta")

            if os.path.isdir(consensus_fasta_path):
                fasta_files = [f for f in os.listdir(consensus_fasta_path) if f.endswith('.fasta')]
                n_fasta_files = len(fasta_files)

                lengths = []
                gap_lengths_all = []
                sequences_without_gaps = 0

                for fasta_file in fasta_files:
                    fasta_path = os.path.join(consensus_fasta_path, fasta_file)
                    record = next(SeqIO.parse(fasta_path, "fasta"))
                    seq = str(record.seq)
                    seq_len = len(seq)
                    lengths.append(seq_len)

                    if "-" not in seq:
                        sequences_without_gaps += 1
                    else:
                        gap_lengths = [len(g.group()) for g in re.finditer(r'-+', seq)]
                        gap_lengths_all.extend(gap_lengths)

                longest_gap = max(gap_lengths_all) if gap_lengths_all else 0
                shortest_gap = min(gap_lengths_all) if gap_lengths_all else 0
                percent_no_gaps = (sequences_without_gaps / n_fasta_files * 100) if n_fasta_files > 0 else 0
                max_length = max(lengths) if lengths else 0
                min_length = min(lengths) if lengths else 0
                avg_length = statistics.mean(lengths) if lengths else 0

                stats = {
                    "n_fasta_files": n_fasta_files,
                    "longest_gap": longest_gap,
                    "shortest_gap": shortest_gap,
                    "percent_without_gaps": round(percent_no_gaps, 2),
                    "max_consensus_length": max_length,
                    "min_consensus_length": min_length,
                    "avg_consensus_length": round(avg_length, 2)
                }

                stats_folder = os.path.join(cluster_path, "consensus_stats")
                os.makedirs(stats_folder, exist_ok=True)

                stats_path = os.path.join(stats_folder, "stats.json")
                with open(stats_path, 'w') as f:
                    json.dump(stats, f, indent=4)

                print(f"Stats saved in: {stats_path}")


def load_all_consensus_sequences(consensus_base_folder):
    consensus_dict = {}

    for cluster_folder in os.listdir(consensus_base_folder):
        cluster_path = os.path.join(consensus_base_folder, cluster_folder)
        consensus_fasta_folder = os.path.join(cluster_path, "consensus_fasta")

        if os.path.isdir(consensus_fasta_folder):
            sequences = []

            for filename in os.listdir(consensus_fasta_folder):
                if filename.endswith((".fasta", ".fa")):
                    filepath = os.path.join(consensus_fasta_folder, filename)
                    
                    try:
                        for record in SeqIO.parse(filepath, "fasta"):
                            sequences.append(str(record.seq))
                    except Exception as e:
                        print(f"Error reading {filepath}: {e}")

            consensus_dict[cluster_folder] = sequences

    return consensus_dict