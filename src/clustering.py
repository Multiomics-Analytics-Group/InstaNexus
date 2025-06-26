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
__date__ = 21 Mar 2025
__maintainer__ = Marco Reverenna
__email__ = marcor@dtu.dk
__status__ = Dev
"""

import os
import subprocess
import shutil
from tempfile import mkdtemp
import pandas as pd
import Bio.SeqIO
from tqdm import tqdm
import logging

def cluster_fasta_files(input_folder):

    cluster_folder = os.path.join(input_folder, "cluster")
    os.makedirs(cluster_folder, exist_ok=True)

    temp_dir = mkdtemp(prefix="mmseqs-") # create a temporary directory for mmseqs

    # Iterate over all fasta files in the folder
    for fasta_file in os.listdir(input_folder):
        if fasta_file.endswith('.fasta'):
            fasta_path = os.path.join(input_folder, fasta_file)
            print(f"the current fasta path is: {fasta_path}")

            if os.path.isfile(fasta_path):

                base_filename = os.path.splitext(fasta_file)[0] # get the base filename with no extension

                prefix = os.path.join(cluster_folder, base_filename) # define the prefix for mmseqs easy-cluster
                
                print(f"Clustering {fasta_file}...") # run mmseqs easy-cluster
                subprocess.run(["mmseqs", "easy-cluster", fasta_path, prefix, temp_dir, "--min-seq-id", "0.85", "-c", "0.8", "--cov-mode", "1", "-v", "1"])
                print(f"Clustering completed for {fasta_file}, results stored with prefix {prefix}")

    shutil.rmtree(temp_dir) 
    print("All clustering tasks completed.")



def process_fasta_and_clusters(fasta_file, cluster_tsv_folder, output_base_folder):
        
    base_filename = os.path.basename(fasta_file).rsplit('.', 1)[0]
    
    cluster_tsv = os.path.join(cluster_tsv_folder, f"{base_filename}_cluster.tsv")

    if not os.path.isfile(cluster_tsv):
        print(f"Cluster TSV file not found for {fasta_file}, skipping.")
        return
    
    output_folder = os.path.join(output_base_folder, f"{base_filename}_cluster_fasta")
    os.makedirs(output_folder, exist_ok=True)
    
    cluster_df = pd.read_csv(cluster_tsv, sep='\t', header=None, names=['cluster', 'contig'])

    records = list(Bio.SeqIO.parse(fasta_file, 'fasta'))

    clusters = cluster_df['cluster'].unique()

    for cluster in tqdm(clusters, desc=f"Processing clusters for {base_filename}"):
        contigs = cluster_df[cluster_df['cluster'] == cluster]['contig'].values
        contig_records = [record for record in records if record.id in contigs]
        Bio.SeqIO.write(contig_records, os.path.join(output_folder, f'{cluster}.fasta'), 'fasta')