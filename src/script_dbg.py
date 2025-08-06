#!/usr/bin/env python

r""" Full assembly script for proteins.
 _____  _______  _    _ 
|  __ \|__   __|| |  | |
| |  | |  | |   | |  | |
| |  | |  | |   | |  | |
| |__| |  | |   | |__| |
|_____/   |_|   |______|

__authors__ = Marco Reverenna
__copyright__ = Copyright 2025-2026
__research-group__ = DTU Biosustain (Multi-omics Network Analytics) and DTU Bioengineering
__date__ = 06 Aug 2025
__maintainer__ = Marco Reverenna
__email__ = marcor@dtu.dk
__status__ = Dev
"""

# !pip install kaleido # to export plotly figures as png
# !pip install --upgrade nbformat # to avoid plotly error

# my modules
import greedy_method as greedy
import dbg
import mapping as map
import consensus as cons
import alignment as align
import clustering as clus
import preprocessing as prep
import compute_statistics as comp_stat

# import libraries
from tqdm import tqdm
from tempfile import mkdtemp
from itertools import combinations
from collections import defaultdict, Counter
from Bio import SeqIO

import sys
import os
import json
import re
import Bio
import shutil
import logging
import importlib
import statistics
import subprocess
import numpy as np
import pandas as pd
import seaborn as sns
import networkx as nx
import plotly.express as px
import matplotlib.pyplot as plt
import plotly.graph_objects as go

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s"
)
logger = logging.getLogger(__name__)


def main():
    """Main function to run the assembly script."""
    
    logger.info("Starting protein assembly pipeline.")

    # Set parameters for the assembly process
    protein = 'MKWVTFISLLLLFSSAYSRGVFRRDTHKSEIAHRFKDLGEEHFKGLVLIAFSQYLQQCPFDEHVKLVNELTEFAKTCVADESHAGCEKSLHTLFGDELCKVASLRETYGDMADCCEKQEPERNECFLSHKDDSPDLPKLKPDPNTLCDEFKADEKKFWGKYLYEIARRHPYFYAPELLYYANKYNGVFQECCQAEDKGACLLPKIETMREKVLASSARQRLRCASIQKFGERALKAWSVARLSQKFPKAEFVEVTKLVTDLTKVHKECCHGDLLECADDRADLAKYICDNQDTISSKLKECCDKPLLEKSHCIAEVEKDAIPENLPPLTADFAEDKDVCKNYQEAKDAFLGSFLYEYSRRHPEYAVSVLLRLAKEYEATLEECCAKDDPHACYSTVFDKLKHLVDEPQNLIKQNCDQFEKLGEYGFQNALIVRYTRKVPQVSTPTLVEVSRSLGKVGTRCCTKPESERMPCTEDYLSLILNRLCVLHEKTPVSEKVTKCCTESLVNRRPCFSALTPDETYVPKAFDEKLFTFHADICTLPDTEKQIKKQTALVELLKHKPKATEEQLKTVMENFVAFVDKCCAADDKEACFAVEGPKLVVSTQTALA'
    proteases = ['Chymotrypsin', 'Legumain', 'Krakatoa', 'Elastase', 'Trypsin', 'Papain', 'Thermo', 'ProtK', 'GluC', 'LysC']
    ass_method = 'dbg'
    run = "bsa"
    chain = ''
    conf = 0.8
    kmer_size = 7
    min_overlap = 3
    min_identity = 0.8
    max_mismatches = 20
    size_threshold = 20

    logger.info("Parameters loaded.")

    folder_outputs = f"../outputs/{run}{chain}"
    prep.create_directory(folder_outputs)
    combination_folder_out = os.path.join(folder_outputs, f"comb_{ass_method}_c{conf}_ks{kmer_size}_ts{size_threshold}_mo{min_overlap}_mi{min_identity}_mm{max_mismatches}")
    prep.create_subdirectories_outputs(combination_folder_out)
    
    logger.info(f"Output folders created at: {combination_folder_out}")

    # Data cleaning
    logger.info("Starting data cleaning...")
    
    protein_norm = prep.normalize_sequence(protein)
    df = pd.read_csv(f"../inputs/{run}.csv")
    df['protease'] = df['experiment_name'].apply(lambda name: prep.extract_protease(name, proteases))
    df = prep.clean_dataframe(df)
    df['cleaned_preds'] = df['preds'].apply(prep.remove_modifications)
    cleaned_psms = df['cleaned_preds'].tolist()
    filtered_psms = prep.filter_contaminants(cleaned_psms, run, "../fasta/contaminants.fasta")
    df = df[df['cleaned_preds'].isin(filtered_psms)]
    df["mapped"] = df["cleaned_preds"].apply(lambda x: "True" if x in protein_norm else "False")
    df = df[df['conf'] > conf]
    df.reset_index(drop=True, inplace=True)    
    final_psms = df['cleaned_preds'].tolist()
    logger.info("Data cleaning completed.")

    # Assembly
    kmers = dbg.get_kmers(final_psms, kmer_size=kmer_size)
    edges = dbg.get_debruijn_edges_from_kmers(kmers)
    assembled_contigs = dbg.assemble_contigs(edges)
    assembled_contigs = sorted(assembled_contigs, key=len, reverse=True)
    assembled_contigs = list(set(assembled_contigs))
    assembled_contigs = [seq for seq in assembled_contigs if len(seq) > size_threshold]
    assembled_contigs = sorted(assembled_contigs, key=len, reverse=True)
    records = [Bio.SeqRecord.SeqRecord(Bio.Seq.Seq(contig), id=f"contig_{idx+1}",
                                    description=f"length: {len(contig)}") for idx,
                                    contig in enumerate(assembled_contigs)]
    Bio.SeqIO.write(records, f"{combination_folder_out}/contigs/{ass_method}_contig_{conf}_{run}.fasta", "fasta")
    mapped_contigs = map.process_protein_contigs_scaffold(assembled_contigs, protein_norm, max_mismatches, min_identity)
    df_contigs = map.create_dataframe_from_mapped_sequences(data = mapped_contigs)
    comp_stat.compute_assembly_statistics(df = df_contigs, sequence_type='contigs', output_folder = f'{combination_folder_out}/statistics',
                                          reference = protein_norm)
    assembled_scaffolds = dbg.create_scaffolds(assembled_contigs, min_overlap)
    assembled_scaffolds = list(set(assembled_scaffolds))
    assembled_scaffolds = sorted(assembled_scaffolds, key=len, reverse=True)
    assembled_scaffolds = [scaffold for scaffold in assembled_scaffolds if len(scaffold) > size_threshold]
    assembled_scaffolds = dbg.merge_sequences(assembled_scaffolds)
    assembled_scaffolds = list(set(assembled_scaffolds))
    assembled_scaffolds = sorted(assembled_scaffolds, key=len, reverse=True)
    assembled_scaffolds = [scaffold for scaffold in assembled_scaffolds if len(scaffold) > size_threshold]

    records = []
    for i, seq in enumerate(assembled_scaffolds):
        record = Bio.SeqRecord.SeqRecord(Bio.Seq.Seq(seq), id=f"scaffold_{i+1}", description=f"length: {len(seq)}")
        records.append(record)

    Bio.SeqIO.write(records, f"{combination_folder_out}/scaffolds/{ass_method}_scaffold_{conf}_{run}.fasta", "fasta")
    mapped_scaffolds = map.process_protein_contigs_scaffold(assembled_contigs = assembled_scaffolds,
                                                            target_protein = protein_norm,
                                                            max_mismatches = max_mismatches,
                                                            min_identity = min_identity
                                                            )
    
    df_scaffolds_mapped = map.create_dataframe_from_mapped_sequences(data = mapped_scaffolds)
    comp_stat.compute_assembly_statistics(df = df_scaffolds_mapped, sequence_type='scaffolds',
                                          output_folder = f"{combination_folder_out}/statistics",
                                          reference = protein_norm)
    
    # Clustering
    scaffolds_folder_out = f"{combination_folder_out}/scaffolds"
    clus.cluster_fasta_files(input_folder = scaffolds_folder_out)

    cluster_tsv_folder = os.path.join(scaffolds_folder_out, "cluster")
    output_base_folder = os.path.join(scaffolds_folder_out, "cluster_fasta")

    for fasta_file in os.listdir(scaffolds_folder_out):
        if fasta_file.endswith('.fasta'):
            fasta_path = os.path.join(scaffolds_folder_out, fasta_file)
            clus.process_fasta_and_clusters(fasta_path, cluster_tsv_folder, output_base_folder)
    
    # Alignment
    cluster_fasta_folder = os.path.join(scaffolds_folder_out, "cluster_fasta") 
    align_folder = os.path.join(scaffolds_folder_out, "align")
    prep.create_directory(align_folder)

    for cluster_folder in os.listdir(cluster_fasta_folder): 
        cluster_folder_path = os.path.join(cluster_fasta_folder, cluster_folder) 
        if os.path.isdir(cluster_folder_path): 
            
            output_cluster_folder = os.path.join(align_folder, cluster_folder) 
            os.makedirs(output_cluster_folder, exist_ok=True) 
                
            for fasta_file in os.listdir(cluster_folder_path): 
                if fasta_file.endswith('.fasta'): 
                    fasta_file_path = os.path.join(cluster_folder_path, fasta_file)
                    base_filename = os.path.splitext(fasta_file)[0] 
                    output_file = os.path.join(output_cluster_folder, f"{base_filename}_out.afa")
                        
                    align.align_or_copy_fasta(fasta_file_path, output_file)

    logger.info("All alignment tasks completed.")

    # Consensus
    consensus_folder = os.path.join(scaffolds_folder_out, "consensus")
    cons.process_alignment_files(align_folder, consensus_folder)

if __name__ == "__main__":
    main()