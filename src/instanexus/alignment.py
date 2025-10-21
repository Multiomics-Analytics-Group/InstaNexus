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
__date__ = 21 Oct 2025
__maintainer__ = Marco Reverenna
__email__ = marcor@dtu.dk
__status__ = Dev
"""

import os
import shutil
import subprocess
from Bio import SeqIO
from pathlib import Path


def align_or_copy_fasta(fasta_file, output_file):
    sequences = list(SeqIO.parse(fasta_file, "fasta"))
    if len(sequences) == 1:
        shutil.copy(fasta_file, output_file)
    else:
        subprocess.run(
            ["clustalo", "-i", fasta_file, "-o", output_file, "--outfmt", "fa"],
            check=True,
        )


def process_alignment(input_folder: str):
    """
    Align all FASTA files in clustering/cluster_fasta and save results in alignment/.
    """
    clustering_dir = Path(input_folder) / "clustering"
    cluster_fasta_folder = clustering_dir / "cluster_fasta"
    alignment_folder = Path(input_folder) / "alignment"
    alignment_folder.mkdir(parents=True, exist_ok=True)

    if not cluster_fasta_folder.exists():
        raise FileNotFoundError(f"Cluster FASTA folder not found: {cluster_fasta_folder}")

    for fasta_file in sorted(os.listdir(cluster_fasta_folder)):
        if not fasta_file.endswith(".fasta"):
            continue
        fasta_path = cluster_fasta_folder / fasta_file
        output_path = alignment_folder / fasta_file.replace(".fasta", "_out.afa")
        align_or_copy_fasta(fasta_path, output_path)

    print("All alignment tasks completed.")

