import os
import shutil
import subprocess
from Bio import SeqIO


def align_or_copy_fasta(fasta_file, output_file):

    sequences = list(SeqIO.parse(fasta_file, "fasta"))
    
    if len(sequences) == 1:
        # If there's only one sequence, copy the file directly
        shutil.copy(fasta_file, output_file)
        #print(f"Copied {fasta_file} to {output_file}")
    else:
        # If there are multiple sequences, run clustalo
        subprocess.run(["clustalo", "-i", fasta_file, "-o", output_file, "--outfmt", "fa", "--force"])
        #print(f"Aligned {fasta_file} and saved to {output_file}")


def process_alignment(input_folder):
    """
    Process all fasta files in the cluster_fasta folder, align them if necessary,
    and save the results in the align folder.
    """
    cluster_fasta_folder = os.path.join(input_folder, "cluster_fasta")
    align_folder = os.path.join(input_folder, "align")

    # Create the align folder if it does not exist
    os.makedirs(align_folder, exist_ok=True)

    # Iterate over all folders in the cluster fasta folder
    for cluster_folder in os.listdir(cluster_fasta_folder):
        cluster_folder_path = os.path.join(cluster_fasta_folder, cluster_folder)
        if os.path.isdir(cluster_folder_path):
            # Create a corresponding folder in the align folder
            output_cluster_folder = os.path.join(align_folder, cluster_folder)
            os.makedirs(output_cluster_folder, exist_ok=True)
            
            # Iterate over all fasta files in the cluster folder
            for fasta_file in os.listdir(cluster_folder_path):
                if fasta_file.endswith('.fasta'):
                    fasta_file_path = os.path.join(cluster_folder_path, fasta_file)
                    base_filename = os.path.splitext(fasta_file)[0]
                    output_file = os.path.join(output_cluster_folder, f"{base_filename}_out.afa")
                    
                    align_or_copy_fasta(fasta_file_path, output_file)

    print("All alignment tasks completed.")