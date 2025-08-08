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
__date__ = 26 Jun 2025
__maintainer__ = Marco Reverenna
__email__ = marcor@dtu.dk
__status__ = Dev
"""

import os
import Bio.SeqIO
import matplotlib.pyplot as plt
import seaborn as sns
from tqdm import tqdm
import matplotlib.patches as patches
import plotly.graph_objects as go
import logging
import pandas as pd

def map_to_protein(seq, protein, max_mismatches, min_identity):
    """ Maps a sequence (`seq`) to a target protein sequence, allowing for mismatches, 
    and identifies the best match based on the maximum mismatches and minimum identity threshold.
    """

    best_match = None
    best_identity = 0

    # Slide `seq` across `protein` to check each possible alignment position
    for i in range(len(protein) - len(seq) + 1):
        mismatches_count = 0
        mismatch_positions = []
        
        # Compare `seq` to the substring of `protein` at the current position
        for j in range(len(seq)):
            if seq[j] != protein[i+j]:
                mismatches_count += 1
                mismatch_positions.append(j)
            
            # Stop checking this alignment if mismatches exceed the allowed threshold
            if mismatches_count > max_mismatches:
                break
        
        # If this alignment meets the mismatch requirement, calculate identity
        if mismatches_count <= max_mismatches:
            if len(seq) == 0:
                print("Zero length sequence found.")
                print(seq)
                continue
            
            identity = 1 - mismatches_count / len(seq)
            
            # Update the best match if this alignment has a higher identity and meets the minimum requirement
            if identity >= min_identity and identity > best_identity:
                best_match = (i, i + len(seq), mismatch_positions, identity)
                best_identity = identity

    return best_match


def process_protein_contigs_scaffold(assembled_contigs, target_protein, max_mismatches, min_identity):
    """ Maps each contig in `assembled_contigs` to a target protein sequence (`target_protein`)
    and identifies which contigs match based on specified mismatch and identity thresholds.
    """
    mapped_sequences = []

    # Map each contig to the target protein
    for contig in assembled_contigs:
        # Attempt to map the contig to the target protein
        target_mapping = map_to_protein(contig, target_protein, max_mismatches=max_mismatches, min_identity=min_identity)
        
        if target_mapping:
            mapped_sequences.append((contig, target_mapping))
    
    return mapped_sequences



def write_mapped_contigs(mapped_contigs, folder, filename_prefix):
    """
    Writes mapped contigs to a FASTA file with detailed annotations for each contig.

    This function takes a list of mapped contigs, creates a sequence record for each contig 
    with information about its length, alignment start and end positions, number of mismatches, 
    and identity. The records are then saved to a specified folder with a filename prefix.

    Parameters:
    - mapped_contigs (list of tuples): List of tuples where each tuple contains:
        - contig (str): The contig sequence.
        - mapping (tuple): A tuple with alignment details `(start, end, mismatches, identity)`,
          where:
            - start (int): Start position of the alignment in the target sequence.
            - end (int): End position of the alignment in the target sequence.
            - mismatches (list of int): List of mismatch positions in the contig.
            - identity (float): Identity of the alignment (fraction of matching positions).
    - folder (str): Directory path where the FASTA file will be saved.
    - filename_prefix (str): Prefix for the output FASTA filename (e.g., "light_chain_mapped").

    Returns:
    - None

    Output:
    - Writes a FASTA file containing each contig sequence with metadata about its alignment.

    Notes:
    - Each record in the FASTA file includes an identifier and a description. The description 
      includes the contig length, start and end positions of the alignment, number of mismatches, 
      and the alignment identity, formatted to two decimal places.
    """

    records = []
    for idx, (contig, mapping) in enumerate(mapped_contigs):
        start, end, mismatches, identity = mapping
        record = Bio.SeqRecord.SeqRecord(
            Bio.Seq.Seq(contig), 
            id=f"Contig {idx+1}", 
            description=f"length: {len(contig)}, start: {start}, end: {end}, mismatches: {len(mismatches)}, identity: {identity:.2f}"
        )
        records.append(record)
    Bio.SeqIO.write(records, os.path.join(folder, f"{filename_prefix}.fasta"), "fasta")



def plot_contigs(mapped_contigs, prot_seq, title, output_file):
    sns.set('paper', 'ticks', 'colorblind', font_scale=1.5)
    _, ax = plt.subplots(figsize=(12, 4))

    ax.add_patch(patches.Rectangle((0, 0), len(prot_seq), 0.2, facecolor='#e6f0ef', edgecolor='#e6f0ef'))

    tracks = {}
    ind = 0
    
    for _, (contig, mapping) in tqdm(enumerate(mapped_contigs), desc="Plotting contigs"):
        start_index, end_index, mismatches, _ = mapping
        
        ind += 1
        placed = False
        for track_num, track in tracks.items():
            if not any(s <= end_index <= e or s <= start_index <= e for s, e in track):
                track.append((start_index, end_index))
                ax.add_patch(patches.Rectangle((start_index, 0.3 + 0.1*track_num),
                                               len(contig),
                                               0.075,
                                               facecolor='#007EA7',
                                               edgecolor='#007EA7',
                                               label="Contig" if not placed else "")
                            )
                placed = True
                break

        if not placed:
            track_num = len(tracks) + 1
            tracks[track_num] = [(start_index, end_index)]
            ax.add_patch(patches.Rectangle((start_index, 0.3 + 0.1*track_num), len(contig), 0.075, facecolor='#007EA7', edgecolor='#007EA7', label="Contig" if not placed else ""))
        
        for mismatch in mismatches:
            ax.add_patch(patches.Rectangle((start_index + mismatch, 0.3 + 0.1*track_num), 1, 0.075, facecolor='#FCAB64', edgecolor='#FCAB64'))

    print(f"Plotted {ind} contigs.")
    
    ax.set_xlim(0, len(prot_seq))
    ax.set_ylim(0, 0.3 + 0.1*(len(tracks)+1))
    ax.get_yaxis().set_visible(False)

    ax.set_xlabel("Sequence range")
    ax.set_title(title)

    handles, labels = ax.get_legend_handles_labels()
    by_label = dict(zip(labels, handles))
    ax.legend(by_label.values(), by_label.keys(), loc="center right", frameon=False, bbox_to_anchor=(1.2, 0.8))

    for spine in ax.spines.values():
        spine.set_visible(False)

    plt.tight_layout(pad=1)
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.close()


def mapping_sequences(mapped_sequences, prot_seq, title, output_folder, output_file):
    """Plot sequences on a sequence using Plotly.
    Parameters:
    mapped_sequences (list): List of tuples containing contigs and their mappings.
    prot_seq (str): The protein sequence.
    title (str): Title of the plot.
    output_file (str): Path to save the output plot.
    """
    fig = go.Figure()

    # Add background rectangle for the protein sequence
    fig.add_shape(type="rect", x0=0, x1=len(prot_seq), y0=0, y1=0.2, fillcolor="#e6f0ef", line=dict(width=0))

    tracks = {}
    ind = 0

    for _, (_, mapping) in tqdm(enumerate(mapped_sequences), desc="Plotting contigs"):
        start_index, end_index, mismatches, _ = mapping

        ind += 1
        placed = False
        for track_num, track in tracks.items():
            if not any(s <= end_index <= e or s <= start_index <= e for s, e in track):
                track.append((start_index, end_index))
                fig.add_shape(type="rect", x0=start_index, x1=end_index, y0=0.3 + 0.1 * track_num,
                    y1=0.375 + 0.1 * track_num, fillcolor="#007EA7", line=dict(color="#007EA7"))
                placed = True
                break

        if not placed:
            track_num = len(tracks) + 1
            tracks[track_num] = [(start_index, end_index)]
            fig.add_shape(type="rect", x0=start_index, x1=end_index,
                y0=0.3 + 0.1 * track_num, y1=0.375 + 0.1 * track_num, fillcolor="#007EA7", line=dict(color="#007EA7"))

        for mismatch in mismatches:
            fig.add_shape(type="rect", x0=start_index + mismatch, x1=start_index + mismatch + 1,
                y0=0.3 + 0.1 * track_num, y1=0.375 + 0.1 * track_num, fillcolor="#FCAB64", line=dict(color="#FCAB64"))

    print(f"Plotted {ind} sequences.")

    fig.add_trace(go.Scatter(x=[None], y=[None], mode='markers', marker=dict(size=10, color="#007EA7"), name="Match"))
    fig.add_trace(go.Scatter(x=[None], y=[None], mode='markers', marker=dict(size=10, color="#FCAB64"), name="Mismatch"))

    fig.update_layout(
        title=title, xaxis=dict(title="Sequence range", range=[0, len(prot_seq)], showgrid=False),
        yaxis=dict(showgrid=False, zeroline=False, showticklabels=False, range=[0, 0.3 + 0.1 * (len(tracks) + 1)]),
        shapes=[], plot_bgcolor="white", width=1200, height=400)

    fig.write_image(f"{output_folder}/{output_file}", scale=2)
    #fig.show()



def create_dataframe_from_mapped_sequences(data):
    """Takes a list of tuples containing sequence data and returns a structured DataFrame.
    Parameters:
        data (list): A list of tuples, where each tuple contains:
                     - A sequence string
                     - A tuple of details (start, end, indices, score)
    Returns:
        pd.DataFrame: A DataFrame with columns for sequence, start, end, indices, and score.
    """
    # Create the initial DataFrame
    df = pd.DataFrame(data, columns=["Sequence", "Details"])

    # Expand the 'Details' column into separate columns
    df[["start", "end", "mismatches_pos", "identity_score"]] = pd.DataFrame(df["Details"].tolist(), index=df.index)

    # Drop the original 'Details' column
    df.drop(columns=["Details"], inplace=True)
    df.rename(columns={"Sequence": "sequence"}, inplace=True)

    return df


def mapping_substitutions(mapped_sequences, prot_seq, title,
                          bar_colors=None, output_file=None, output_folder=".", 
                          contig_colors="#6baed6", match_color="#6baed6"):

    default_colors = {
        "match": match_color,
        "mismatch": "#b30000",
        "D_to_N": "#000000",
        "E_to_Q": "#A8A29E"
    }
    colors = {**default_colors, **(bar_colors or {})}

    fig = go.Figure()

    fig.add_shape(
        type="rect",
        x0=0,
        x1=len(prot_seq),
        y0=0,
        y1=0.2,
        fillcolor="#e6f0ef",
        line=dict(width=0)
    )

    tracks = {}
    ind = 0

    for seq, mapping in tqdm(mapped_sequences, desc="Plotting contigs"):
        start_index, end_index, mismatches, _ = mapping
        ind += 1

        contig_color = contig_colors[ind % len(contig_colors)] if isinstance(contig_colors, list) else contig_colors

        placed = False
        for track_num, track in tracks.items():
            if not any(s <= end_index <= e or s <= start_index <= e for s, e in track):
                track.append((start_index, end_index))
                y0 = 0.3 + 0.1 * track_num
                y1 = y0 + 0.075
                fig.add_shape(type="rect", x0=start_index, x1=end_index, y0=y0, y1=y1,
                              fillcolor=contig_color, line=dict(color=contig_color))
                placed = True
                break

        if not placed:
            track_num = len(tracks)
            tracks[track_num] = [(start_index, end_index)]
            y0 = 0.3 + 0.1 * track_num
            y1 = y0 + 0.075
            fig.add_shape(type="rect", x0=start_index, x1=end_index, y0=y0, y1=y1,
                          fillcolor=contig_color, line=dict(color=contig_color))

        for mismatch in mismatches:
            abs_index = start_index + mismatch
            if abs_index >= len(prot_seq) or mismatch >= len(seq):
                continue

            ref_aa = prot_seq[abs_index]
            query_aa = seq[mismatch]

            if query_aa == 'D' and ref_aa == 'N':
                color = colors["D_to_N"]
            elif query_aa == 'E' and ref_aa == 'Q':
                color = colors["E_to_Q"]
            else:
                color = colors["mismatch"]

            fig.add_shape(
                type="rect",
                x0=abs_index,
                x1=abs_index + 1,
                y0=y0,
                y1=y1,
                fillcolor=color,
                line=dict(color=color)
            )

    print(f"Plotted {ind} sequences.")

    # Legend
    fig.add_trace(go.Scatter(x=[None], y=[None], mode='markers', marker=dict(size=10, color=colors["match"]), name="Match"))
    fig.add_trace(go.Scatter(x=[None], y=[None], mode='markers', marker=dict(size=10, color=colors["mismatch"]), name="Mismatch"))
    fig.add_trace(go.Scatter(x=[None], y=[None], mode='markers', marker=dict(size=10, color=colors["D_to_N"]), name="Seq:D → Ref:N"))
    fig.add_trace(go.Scatter(x=[None], y=[None], mode='markers', marker=dict(size=10, color=colors["E_to_Q"]), name="Seq:E → Ref:Q"))

    fig.update_layout(
        title=title,
        legend=dict(
            title=dict(text="Legend"),
            orientation="h",
            x=0.5, xanchor="center",
            y=1.05, yanchor="bottom"
        ),
        showlegend=True,
        xaxis=dict(title="Sequence range", range=[0, len(prot_seq)], showgrid=False),
        yaxis=dict(showgrid=False, zeroline=False, showticklabels=False,
                   range=[0, 0.3 + 0.1 * (len(tracks) + 1)]),
        plot_bgcolor="white",
        width=1200,
        height=400,
        font=dict(size=14, family="Arial, sans-serif", color='black')
    )

    if output_file:
        os.makedirs(output_folder, exist_ok=True)
        fig.write_image(os.path.join(output_folder, output_file), scale=2)

    # fig.show()


def mapping_psms_protease_associated(mapped_sequences, prot_seq, labels, palette, title, output_folder, output_file):

    fig = go.Figure()

    fig.add_shape(
        type="rect", 
        x0=0, 
        x1=len(prot_seq), 
        y0=0, 
        y1=0.2, 
        fillcolor="#e6f0ef", 
        line=dict(width=0)
    )
    tracks = {}
    ind = 0
    unique_labels = []
    for lab in labels:
        if lab not in unique_labels:
            unique_labels.append(lab)
    
    label_color = {lab: palette.get(lab, "#000000") for lab in unique_labels}
    
    for idx, (_, mapping) in tqdm(enumerate(mapped_sequences), desc="Plotting contigs"):
        start_index, end_index, mismatches, _ = mapping
        lab = labels[idx]
        ind += 1
        placed = False
        for track_num, track in tracks.items():
            if not any(s <= end_index <= e or s <= start_index <= e for s, e in track):
                track.append((start_index, end_index))
                y0 = 0.3 + 0.1 * track_num
                y1 = 0.375 + 0.1 * track_num
                fig.add_shape(
                    type="rect", 
                    x0=start_index, 
                    x1=end_index, 
                    y0=y0, 
                    y1=y1,
                    fillcolor=label_color[lab], 
                    line=dict(color=label_color[lab])
                )
                placed = True
                break
        if not placed:
            track_num = len(tracks) + 1
            tracks[track_num] = [(start_index, end_index)]
            y0 = 0.3 + 0.1 * track_num
            y1 = 0.375 + 0.1 * track_num
            fig.add_shape(
                type="rect", 
                x0=start_index, 
                x1=end_index, 
                y0=y0, 
                y1=y1,
                fillcolor=label_color[lab], 
                line=dict(color=label_color[lab])
            )
    
    print(f"Plotted {ind} sequences.")
    for lab, col in label_color.items():
        fig.add_trace(
            go.Scatter(
                x=[None], 
                y=[None], 
                mode='markers',
                marker=dict(size=10, color=col, symbol='square'),
                name=lab
            )
        )
    
    fig.update_layout(
        title=title,
        legend_title="Proteases",
        legend=dict(orientation='h', x=0.5, y=1.1, xanchor='center', yanchor='bottom', font=dict(size=16)),
        margin=dict(l=20, r=20, t=100, b=20),
        xaxis=dict(title="Reference", range=[0, len(prot_seq)], showgrid=False, dtick=50, tick0=0, tickfont=dict(size=16)),
        yaxis=dict(showgrid=False, zeroline=False, showticklabels=False, range=[0, 0.3 + 0.1 * (len(tracks) + 1)]),
        shapes=[], 
        plot_bgcolor="white", 
        width=1200, 
        height= 600, 
        showlegend=True,
        font=dict(size=14, family="Arial, sans-serif", color='black')
    )
    
    fig.write_image(f"{output_folder}/{output_file}", scale=2)
    # fig.show()