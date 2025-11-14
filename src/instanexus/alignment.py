#!/usr/bin/env python

r"""Alignment module for clustered scaffolds.

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
__date__ = 14 Nov 2025
__maintainer__ = Marco Reverenna
__email__ = marcor@dtu.dk
__status__ = Dev
"""

import argparse
import logging
import os
import shutil
import subprocess
from pathlib import Path
from Bio import SeqIO

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
logger = logging.getLogger(__name__)


def align_or_copy_fasta(fasta_file, output_file):
    """
    Aligns a FASTA file using clustalo if it has >1 sequence,
    otherwise just copies it.
    """
    try:
        sequences = list(SeqIO.parse(fasta_file, "fasta"))
    except Exception as e:
        logger.error(f"Could not parse FASTA file {fasta_file}: {e}")
        return

    if len(sequences) == 1:
        # Only one sequence, no alignment needed
        shutil.copy(fasta_file, output_file)
        logger.debug(f"Copied single-sequence file: {Path(fasta_file).name}")
    elif len(sequences) > 1:
        # Multiple sequences, run clustalo
        logger.debug(f"Aligning {len(sequences)} sequences from {Path(fasta_file).name}...")
        try:
            subprocess.run(
                ["clustalo", "-i", fasta_file, "-o", output_file, "--outfmt", "fa", "--force"],
                check=True,
                capture_output=True, # Suppress clustalo stdout
                text=True
            )
        except FileNotFoundError:
            logger.error("clustalo command not found. Please ensure it is in your system's PATH.")
            raise
        except subprocess.CalledProcessError as e:
            logger.error(f"Clustalo failed for {fasta_file}: {e.stderr}")
            raise
    else:
        # Zero sequences
        logger.warning(f"Skipping empty FASTA file: {Path(fasta_file).name}")


def process_alignment(input_dir: str, output_dir: str):
    """
    Align all FASTA files from input_dir and save results in output_dir.
    
    Args:
        input_dir (str): Path to the .../cluster_fasta/ directory.
        output_dir (str): Path to the .../alignment/ directory.
    """
    cluster_fasta_folder = Path(input_dir)
    alignment_folder = Path(output_dir)
    
    alignment_folder.mkdir(parents=True, exist_ok=True)

    if not cluster_fasta_folder.exists():
        logger.error(f"Cluster FASTA folder not found: {cluster_fasta_folder}")
        raise FileNotFoundError(f"Cluster FASTA folder not found: {cluster_fasta_folder}")

    fasta_files_to_align = [
        f for f in sorted(os.listdir(cluster_fasta_folder))
        if f.endswith(".fasta")
    ]

    logger.info(f"Found {len(fasta_files_to_align)} cluster FASTA files to align.")

    for fasta_file in fasta_files_to_align:
        fasta_path = cluster_fasta_folder / fasta_file
        # Save output as .afa (alignment FASTA)
        output_path = alignment_folder / fasta_file.replace(".fasta", "_aligned.afa") 
        
        align_or_copy_fasta(fasta_path, output_path)

    logger.info("All alignment tasks completed.")


def main(input_cluster_fasta_folder: str, 
    output_alignment_folder: str):
    """
    Main function to run the alignment script.
    """
    logger.info("--- Starting Step 4: Alignment ---")

    logger.info(f"Input Folder (cluster FASTA): {input_cluster_fasta_folder}")
    logger.info(f"Output Folder (Alignments): {output_alignment_folder}")

    process_alignment(
        input_dir=input_cluster_fasta_folder,
        output_dir=output_alignment_folder
    )
    
    logger.info("--- Step 4: Alignment Completed ---")


def cli():
    """
    Command-line interface (CLI) for the alignment script.
    """
    parser = argparse.ArgumentParser(
        description="Alignment script for clustered scaffolds."
    )
    
    parser.add_argument(
        "--input-folder",
        type=str,
        required=True,
        help="Path to the folder containing cluster FASTA files (e.g., .../cluster_fasta)."
    )
    parser.add_argument(
        "--output-folder",
        type=str,
        required=True,
        help="Path to the folder to save aligned .afa files (e.g., .../alignment)."
    )
    
    args = parser.parse_args()

    main(input_cluster_fasta_folder=args.input_folder,
         output_alignment_folder=args.output_folder
         )


if __name__ == "__main__":
    cli()

# python -m instanexus.alignment \
#     --input-folder outputs/bsa/scaffolds/clustering/cluster_fasta \
#     --output-folder outputs/bsa/scaffolds/alignment