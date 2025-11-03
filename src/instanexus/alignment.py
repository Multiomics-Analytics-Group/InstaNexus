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
__date__ = 03 Nov 2025
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

# Setup logging
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
        logger.info(f"Copied single-sequence file: {Path(fasta_file).name}")
    elif len(sequences) > 1:
        # Multiple sequences, run clustalo
        logger.info(f"Aligning {len(sequences)} sequences from {Path(fasta_file).name}...")
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


def process_alignment(scaffolds_folder: str):
    """
    Align all FASTA files in .../scaffolds/clustering/cluster_fasta
    and save results in .../scaffolds/alignment/.
    
    Args:
        scaffolds_folder (str): Path to the .../scaffolds/ directory.
    """
    scaffolds_folder_path = Path(scaffolds_folder)
    clustering_dir = scaffolds_folder_path / "clustering"
    cluster_fasta_folder = clustering_dir / "cluster_fasta"
    alignment_folder = scaffolds_folder_path / "alignment"
    
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


def main(combo_folder: str):
    """
    Main function to run the alignment script.
    """
    logger.info("--- Starting Step 4: Alignment ---")
    
    combo_folder_path = Path(combo_folder)
    scaffolds_folder_path = combo_folder_path / "scaffolds"

    if not scaffolds_folder_path.exists():
        logger.error(f"Scaffolds folder not found: {scaffolds_folder_path}")
        raise FileNotFoundError(f"Scaffolds folder not found: {scaffolds_folder_path}")

    logger.info(f"Input (Scaffolds Folder): {scaffolds_folder_path}")
    
    # Call the core logic function, passing the .../scaffolds/ path
    process_alignment(scaffolds_folder=str(scaffolds_folder_path))
    
    logger.info("--- Step 4: Alignment Completed ---")


def cli():
    """
    Command-line interface (CLI) for the alignment script.
    """
    parser = argparse.ArgumentParser(
        description="Alignment script for clustered scaffolds."
    )
    
    parser.add_argument(
        "--combo-folder",
        type=str,
        required=True,
        help="Path to the 'comb_...' folder (output of the assembly step)."
    )
    
    args = parser.parse_args()
    
    main(combo_folder=args.combo_folder)


if __name__ == "__main__":
    cli()