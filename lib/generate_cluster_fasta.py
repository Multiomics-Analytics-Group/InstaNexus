import os
import argparse
import pandas as pd
import Bio.SeqIO
from tqdm import tqdm


def process_fasta_and_clusters(fasta_file, cluster_tsv_folder, output_base_folder):
        
    # Get the base filename without the extension
    base_filename = os.path.basename(fasta_file).rsplit('.', 1)[0]
    print(f"Processing {base_filename}")
    
    # Construct the corresponding cluster TSV file path
    cluster_tsv = os.path.join(cluster_tsv_folder, f"{base_filename}_cluster.tsv")
    print(f"Cluster TSV file: {cluster_tsv}")

    # Check if the cluster TSV file exists
    if not os.path.isfile(cluster_tsv):
        print(f"Cluster TSV file not found for {fasta_file}, skipping.")
        return
    
    # Create output directory for this fasta file
    output_folder = os.path.join(output_base_folder, f"{base_filename}_cluster_fasta")
    os.makedirs(output_folder, exist_ok=True)
    print(f"Output folder: {output_folder}")
    
    # Read the cluster TSV file
    cluster_df = pd.read_csv(cluster_tsv, sep='\t', header=None, names=['cluster', 'contig'])

    # Read the original fasta file
    records = list(Bio.SeqIO.parse(fasta_file, 'fasta'))

    # Get unique clusters
    clusters = cluster_df['cluster'].unique()
    print(f"Processing {base_filename}: Number of clusters: {len(clusters)}")

    # Generate fasta files for each cluster
    for cluster in tqdm(clusters, desc=f"Processing clusters for {base_filename}"):
        contigs = cluster_df[cluster_df['cluster'] == cluster]['contig'].values
        contig_records = [record for record in records if record.id in contigs]
        Bio.SeqIO.write(contig_records, os.path.join(output_folder, f'{cluster}.fasta'), 'fasta')

    print(f"Processing {base_filename} done, written to {output_folder} with {len(clusters)} clusters.")