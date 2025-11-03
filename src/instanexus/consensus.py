#!/usr/bin/env python

r"""Consensus generation module for aligned scaffolds.

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
__date__ = 03 Nov 2025
__maintainer__ = Marco Reverenna
__email__ = marcor@dtu.dk
__status__ = Dev
"""

import argparse
import logging
import os
import json
import re
import statistics
from pathlib import Path
from collections import Counter
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import Bio.SeqRecord 
from tqdm import tqdm
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import plotly.express as px
import plotly.io as pio
import logomaker

# Setup logging
logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
logger = logging.getLogger(__name__)


def generate_pssm(aligned_records):
    """Generates a Position-Specific Scoring Matrix (PSSM) from aligned sequences."""
    pssm = {}
    for record in aligned_records:
        for i, aa in enumerate(record.seq):
            if i not in pssm:
                pssm[i] = Counter()
            if aa != "-":  # Ignore gaps
                pssm[i][aa] += 1
    
    # Normalize to frequencies
    for i in pssm:
        total = sum(pssm[i].values())
        if total > 0:
            for aa in pssm[i]:
                pssm[i][aa] /= total
                
    pssm_df = pd.DataFrame(pssm).fillna(0).T
    pssm_df.index = pssm_df.index + 1  # 1-based indexing for positions
    pssm_df = pssm_df.sort_index(axis=1) # Sort columns alphabetically (A, C, D...)
    return pssm_df


def generate_consensus(pssm_df, threshold=0.6):
    """Generates a consensus sequence from a PSSM DataFrame."""
    consensus = ""
    for i in pssm_df.index:
        if pssm_df.loc[i].max() > threshold:
            consensus += pssm_df.loc[i].idxmax()
        else:
            consensus += "-"  # Use gap if no amino acid exceeds threshold
    return consensus


def plot_heatmap2(pssm_df, output_file):
    """Generates and saves a Plotly heatmap from a PSSM."""
    df_t = pssm_df.T
    fig = px.imshow(
        df_t,
        color_continuous_scale="Reds",
        aspect="auto",
        labels={"x": "Position", "y": "Amino acid", "color": "Frequency"},
    )
    fig.update_layout(
        font=dict(family="DejaVu Sans", size=15),
        xaxis_title="Position",
        yaxis_title="Amino acid",
        margin=dict(l=40, r=40, t=40, b=40),
        coloraxis_colorbar=dict(title="Frequency"),
    )
    # Set x-ticks to be every 5 positions
    x_tickvals = list(range(0, df_t.shape[1], 5))
    fig.update_xaxes(tickmode="array", tickvals=x_tickvals, tickangle=0)
    fig.update_yaxes(autorange="reversed")
    pio.write_image(fig, output_file, format="svg", width=1200, height=400, scale=2)


def plot_logo2(pssm_df, output_file):
    """Generates and saves a Logomaker sequence logo from a PSSM."""
    max_fig_width = 50
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
    
    logo.style_xticks(anchor=0, rotation=0, spacing=1, fontsize=16, ha="center")
    ax.set_yticks([0, 0.5, 1])
    ax.set_yticklabels([0, 0.5, 1], fontsize=16)
    ax.set_ylabel("Frequency", fontsize=18)
    ax.set_xlabel("Position", fontsize=18)
    logo.style_spines(spines=["left", "right", "bottom", "top"], visible=False)
    plt.tight_layout()
    plt.savefig(output_file, format="svg", dpi=300)
    plt.close(fig)


def run_consensus_generation(align_folder: str, output_folder: str, run_id: str = ""):
    """
    Core logic: Process all .afa files from alignment folder.
    Generate consensus sequences, heatmaps, and logos.
    """
    align_path = Path(align_folder)
    output_path = Path(output_folder)
    output_path.mkdir(parents=True, exist_ok=True)

    if not align_path.exists():
        logger.error(f"Alignment folder not found: {align_path}")
        raise FileNotFoundError(f"Alignment folder not found: {align_path}")

    consensus_fasta_dir = output_path / "consensus_fasta"
    heatmap_dir = output_path / "heatmap"
    logo_dir = output_path / "logo"

    consensus_fasta_dir.mkdir(exist_ok=True)
    heatmap_dir.mkdir(exist_ok=True)
    logo_dir.mkdir(exist_ok=True)

    alignment_files = [f for f in sorted(os.listdir(align_path)) if f.endswith(".afa")]
    
    logger.info(f"Found {len(alignment_files)} aligned .afa files.")
    
    for alignment_file in tqdm(
        alignment_files,
        desc=f"[{run_id}] Generating Consensus" if run_id else "Generating Consensus",
    ):
        alignment_path = align_path / alignment_file
        base_filename = alignment_path.stem # e.g., 'scaffold_0001_aligned'

        aligned_records = list(SeqIO.parse(alignment_path, "fasta"))
        if not aligned_records:
            logger.warning(f"Skipping empty alignment file: {alignment_file}")
            continue

        pssm_df = generate_pssm(aligned_records)
        consensus_sequence = generate_consensus(pssm_df)

        # Write consensus FASTA
        consensus_record = Bio.SeqRecord.SeqRecord(
            Bio.Seq.Seq(consensus_sequence),
            id=base_filename.replace("_aligned", "_consensus"),
            description="Consensus sequence",
        )
        consensus_fasta_path = consensus_fasta_dir / f"{base_filename}_consensus.fasta"
        Bio.SeqIO.write([consensus_record], consensus_fasta_path, "fasta")

        # Plot heatmap
        heatmap_path = heatmap_dir / f"{base_filename}_heatmap.svg"
        plot_heatmap2(pssm_df, heatmap_path)

        # Plot logo
        logo_path = logo_dir / f"{base_filename}_logo.svg"
        plot_logo2(pssm_df, logo_path)

    logger.info("All consensus tasks completed.")
    return consensus_fasta_dir


def generate_consensus_stats(consensus_base_folder):
    """Calculates statistics on the generated consensus FASTA files."""
    consensus_base_folder = Path(consensus_base_folder)
    fasta_files = list(consensus_base_folder.glob("*_consensus.fasta"))
    
    if not fasta_files:
        logger.warning("No consensus FASTA files found, skipping stats.")
        return

    lengths = []
    gap_lengths_all = []
    sequences_without_gaps = 0
    
    for fasta_path in fasta_files:
        record = next(SeqIO.parse(fasta_path, "fasta"))
        seq = str(record.seq)
        seq_len = len(seq)
        lengths.append(seq_len)
        
        if "-" not in seq:
            sequences_without_gaps += 1
        else:
            gap_lengths = [len(g.group()) for g in re.finditer(r"-+", seq)]
            gap_lengths_all.extend(gap_lengths)
            
    longest_gap = max(gap_lengths_all) if gap_lengths_all else 0
    shortest_gap = min(gap_lengths_all) if gap_lengths_all else 0
    percent_no_gaps = (sequences_without_gaps / len(fasta_files) * 100) if fasta_files else 0
    max_length = max(lengths) if lengths else 0
    min_length = min(lengths) if lengths else 0
    avg_length = statistics.mean(lengths) if lengths else 0
    
    stats = {
        "n_consensus_files": len(fasta_files),
        "longest_gap": longest_gap,
        "shortest_gap": shortest_gap,
        "percent_without_gaps": round(percent_no_gaps, 2),
        "max_consensus_length": max_length,
        "min_consensus_length": min_length,
        "avg_consensus_length": round(avg_length, 2),
    }
    
    stats_path = consensus_base_folder.parent / "consensus_stats.json"
    with open(stats_path, "w") as f:
        json.dump(stats, f, indent=4)
    logger.info(f"Consensus statistics saved to: {stats_path}")


def main(
    combo_folder: str, 
    run_name: str, 
    assembly_mode: str, 
    conf: float
):
    """
    Main function to run the consensus generation script.
    """
    logger.info("--- Starting Step 5: Consensus Generation ---")
    
    # Reconstruct paths
    combo_folder_path = Path(combo_folder)
    scaffolds_folder_path = combo_folder_path / "scaffolds"
    align_folder_in = scaffolds_folder_path / "alignment"
    consensus_folder_out = scaffolds_folder_path / "consensus"

    logger.info(f"Alignment Folder (Input): {align_folder_in}")
    logger.info(f"Consensus Folder (Output): {consensus_folder_out}")

    # Recreate the run_id for the tqdm progress bar
    run_id = f"{assembly_mode}_{conf}_{run_name}"

    # --- Step 1: Generate consensus, heatmaps, and logos ---
    logger.info("Running consensus generation from alignment files...")
    consensus_fasta_dir = run_consensus_generation(
        align_folder=str(align_folder_in),
        output_folder=str(consensus_folder_out),
        run_id=run_id
    )
    
    # --- Step 2: Generate statistics on the consensus files ---
    if consensus_fasta_dir:
        logger.info("Running consensus statistics generation...")
        generate_consensus_stats(
            consensus_base_folder=consensus_fasta_dir
        )
    
    logger.info("--- Step 5: Consensus Generation Completed ---")


def cli():
    """
    Command-line interface (CLI) for the consensus generation script.
    """
    parser = argparse.ArgumentParser(
        description="Consensus generation script for aligned scaffolds."
    )
    
    parser.add_argument(
        "--combo-folder",
        type=str,
        required=True,
        help="Path to the 'comb_...' folder (output of the assembly step)."
    )
    parser.add_argument(
        "--run-name",
        type=str,
        required=True,
        help="The base name of the run (e.g., 'bsa')."
    )
    parser.add_argument(
        "--assembly-mode",
        type=str,
        default="greedy",
        help="Assembly mode used (e.g., 'greedy', 'dbg')."
    )
    parser.add_argument(
        "--conf",
        type=float,
        default=0.88,
        help="Confidence score used in the run."
    )
    
    args = parser.parse_args()
    
    main(
        combo_folder=args.combo_folder,
        run_name=args.run_name,
        assembly_mode=args.assembly_mode,
        conf=args.conf
    )

if __name__ == "__main__":
    cli()