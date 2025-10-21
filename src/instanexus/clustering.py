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
__date__ = 21 Oct 2025
__maintainer__ = Marco Reverenna
__email__ = marcor@dtu.dk
__status__ = Dev
"""

import os
import shutil
import subprocess
from tempfile import mkdtemp
from pathlib import Path
from Bio import SeqIO
import pandas as pd
from tqdm import tqdm


def cluster_fasta_files(input_folder: str):
    clustering_dir = Path(input_folder) / "clustering"
    cluster_fasta_dir = clustering_dir / "cluster_fasta"
    clustering_dir.mkdir(parents=True, exist_ok=True)
    cluster_fasta_dir.mkdir(exist_ok=True)

    temp_dir = mkdtemp(prefix="mmseqs-")

    for fasta_file in sorted(os.listdir(input_folder)):
        if not fasta_file.endswith(".fasta"):
            continue

        fasta_path = Path(input_folder) / fasta_file
        if not fasta_path.is_file():
            continue

        base_filename = fasta_path.stem
        prefix = clustering_dir / base_filename

        subprocess.run(
            [
                "mmseqs",
                "easy-cluster",
                str(fasta_path),
                str(prefix),
                temp_dir,
                "--min-seq-id",
                "0.85",
                "-c",
                "0.8",
                "--cov-mode",
                "1",
                "-v",
                "1",
            ],
            check=True,
        )

        tsv_src = Path(f"{prefix}_cluster.tsv")
        rep_src = Path(f"{prefix}_rep_seq.fasta")
        all_src = Path(f"{prefix}_all_seqs.fasta")


        if tsv_src.exists():
            shutil.move(str(tsv_src), clustering_dir / "scaffolds_cluster.tsv")
        if rep_src.exists():
            shutil.move(str(rep_src), clustering_dir / "scaffolds_rep_seq.fasta")
        if all_src.exists():
            shutil.move(str(all_src), clustering_dir / "scaffolds_all_seqs.fasta")

    shutil.rmtree(temp_dir)
    print("All clustering tasks completed.")


def process_fasta_and_clusters(fasta_file: str, input_folder: str):
    clustering_dir = Path(input_folder) / "clustering"
    cluster_tsv = clustering_dir / "scaffolds_cluster.tsv"
    cluster_fasta_dir = clustering_dir / "cluster_fasta"
    cluster_fasta_dir.mkdir(exist_ok=True)

    if not cluster_tsv.is_file():
        print(f"Cluster TSV file not found for {fasta_file}, skipping.")
        return

    cluster_df = pd.read_csv(cluster_tsv, sep="\t", header=None, names=["cluster", "contig"])
    records = list(SeqIO.parse(fasta_file, "fasta"))
    clusters = cluster_df["cluster"].unique()

    width = max(4, len(str(len(clusters))))

    for i, cluster in enumerate(tqdm(clusters, desc="Processing clusters")):
        contigs = cluster_df[cluster_df["cluster"] == cluster]["contig"].values
        contig_records = [r for r in records if r.id in contigs]
        cluster_name = f"scaffold_{str(i+1).zfill(width)}.fasta"
        SeqIO.write(contig_records, cluster_fasta_dir / cluster_name, "fasta")

    print("All cluster FASTA files created.")
