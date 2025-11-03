#!/usr/bin/env python

r"""Pipeline script for InstaNexus.

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

# !pip install kaleido # to export plotly figures as png
# !pip install --upgrade nbformat # to avoid plotly error

import argparse
import logging
from pathlib import Path
import preprocessing
import assembly 
import clustering
import alignment
import consensus

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
logger = logging.getLogger(__name__)

def cli():
    """Command-line interface for the master pipeline."""
    
    parser = argparse.ArgumentParser(
        description="Run the full InstaNexus preprocessing and assembly pipeline."
    )
    
    parser.add_argument(
        "--input-csv",
        type=str,
        required=True,
        help="Path to the RAW input CSV file containing PSM data.",
    )
    parser.add_argument(
        "--folder-outputs",
        type=str,
        default="outputs",
        help="Base folder to save all run outputs.",
    )
    parser.add_argument(
        "--chain",
        type=str,
        default="",
        help="Chain identifier for the sample (e.g., 'light', 'heavy').",
    )
    parser.add_argument(
        "--reference",
        action="store_true",
        help="Enable reference-based mode for statistics.",
    )
    parser.add_argument(
        "--conf",
        type=float,
        default=0.88,
        help="Confidence threshold for filtering (default: 0.88).",
    )
    parser.add_argument(
        "--assembly-mode",
        type=str,
        choices=["dbg", "greedy"],
        default="greedy",
        help="Assembly algorithm to use.",
    )
    parser.add_argument(
        "--kmer-size",
        type=int,
        default=6,
        help="K-mer size (only used if --assembly-mode dbg).",
    )
    parser.add_argument(
        "--min-overlap",
        type=int,
        default=4,
        help="Minimum overlap size between reads.",
    )
    parser.add_argument(
        "--size-threshold",
        type=int,
        default=10,
        help="Minimum contig size threshold.",
    )
    parser.add_argument(
        "--min-identity",
        type=float,
        default=0.8,
        help="Minimum identity threshold (only used if --reference).",
    )
    parser.add_argument(
        "--max-mismatches",
        type=int,
        default=14,
        help="Maximum allowed mismatches (only used if --reference).",
    )
    parser.add_argument(
        "--min-seq-id",
        type=float,
        default=0.85,
        help="Minimum sequence identity for mmseqs (default: 0.85)."
    )
    parser.add_argument(
        "--coverage",
        type=float,
        default=0.8,
        help="Coverage parameter (-c) for mmseqs (default: 0.8)."
    )
    
    args = parser.parse_args()
    
    run_pipeline(args)


def run_pipeline(args):
    """
    Orchestrates the preprocessing and assembly steps.
    """
    
    logger.info("--- InstaNexus Pipeline started ---")
    
    try:
        preprocessing.main(
            input_csv=args.input_csv,
            chain=args.chain,
            folder_outputs=args.folder_outputs,
            reference=args.reference,
            assembly_mode=args.assembly_mode,
            conf=args.conf,
            kmer_size=args.kmer_size,
            size_threshold=args.size_threshold,
            min_overlap=args.min_overlap,
            min_identity=args.min_identity,
            max_mismatches=args.max_mismatches
        )
    except Exception as e:
        logger.error(f"Preprocessing failed: {e}")
        return 

    run = Path(args.input_csv).stem
    base_output_folder = Path(args.folder_outputs) / run
    
    folder_name_parts = [f"comb_{args.assembly_mode}", f"c{args.conf}", f"ts{args.size_threshold}", f"mo{args.min_overlap}"]
    
    if args.assembly_mode == "dbg":
        if args.kmer_size is None:
             args.kmer_size = cli.get_default("kmer_size")
        folder_name_parts.insert(2, f"ks{args.kmer_size}")

    if args.reference:
        folder_name_parts.extend([f"mi{args.min_identity}", f"mm{args.max_mismatches}"])

    combination_folder_out = base_output_folder / "_".join(folder_name_parts)
    cleaned_csv_path = combination_folder_out / "cleaned" / "cleaned_data.csv"

    try:
        assembly.main(
            input_data=str(cleaned_csv_path),
            output_folder=str(combination_folder_out),
            assembly_mode=args.assembly_mode,
            kmer_size=args.kmer_size,
            min_overlap=args.min_overlap,
            size_threshold=args.size_threshold,
            reference=args.reference,
            chain=args.chain,
            min_identity=args.min_identity,
            max_mismatches=args.max_mismatches
        )
    except Exception as e:
        logger.error(f"Assembly failed: {e}")
        return

    try:
        clustering.main(
            combo_folder=str(combination_folder_out),
            run_name=run,
            assembly_mode=args.assembly_mode,
            conf=args.conf,
            min_seq_id=args.min_seq_id,
            coverage=args.coverage
        )
    except Exception as e:
        logger.error(f"Clustering failed: {e}")
        return
    
    try:
        alignment.main(
            combo_folder=str(combination_folder_out)
        )
    except Exception as e:
        logger.error(f"Alignment failed: {e}")
        return
    
    try:
        consensus.main(
            combo_folder=str(combination_folder_out),
            run_name=run,
            assembly_mode=args.assembly_mode,
            conf=args.conf
        )
    except Exception as e:
        logger.error(f"Consensus failed: {e}")
        return

    logger.info("--- InstaNexus Pipeline finished successfully! ---")


if __name__ == "__main__":
    cli()

# example command to run the pipeline:
# python main.py --input-csv ../../inputs/bsa.csv --folder-outputs ../../outputs --assembly-mode dbg --conf 0.92 --kmer-size 7 --size-threshold 12 --min-overlap 3 --min-seq-id 0.85