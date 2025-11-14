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
__date__ = 14 Nov 2025
__maintainer__ = Marco Reverenna
__email__ = marcor@dtu.dk
__status__ = Dev
"""

import argparse
import logging
from pathlib import Path
from . import preprocessing
from . import assembly 
from . import clustering
from . import alignment
from . import consensus

# Setup logging
logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
logging.getLogger("instanexus.preprocessing").setLevel(logging.WARNING)
logging.getLogger("instanexus.assembly").setLevel(logging.WARNING)
logging.getLogger("instanexus.clustering").setLevel(logging.WARNING)
logging.getLogger("instanexus.alignment").setLevel(logging.WARNING)
logging.getLogger("instanexus.consensus").setLevel(logging.WARNING)
logger = logging.getLogger(__name__)

def cli():
    """Command-line interface for the master pipeline."""
    
    parser = argparse.ArgumentParser(
        description="Run the full InstaNexus preprocessing and assembly pipeline."
    )
    
    # --- INPUT ARGUMENTS ---
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
        "--metadata-json-path",
        type=str,
        required=True,
        help="Path to the sample_metadata.json file (required by preprocessing and assembly)."
    )
    parser.add_argument(
        "--contaminants-fasta-path",
        type=str,
        required=True,
        help="Path to the contaminants.fasta file (required by preprocessing)."
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
        default=None, # Default to None for preprocessing logic
        help="Confidence threshold for filtering (optional).",
    )
    parser.add_argument(
        "--assembly-mode",
        type=str,
        choices=["greedy", "dbg", "dbg_weighted", "dbgX", "fusion"],
        default="greedy",
        help="Assembly algorithm to use.",
    )
    parser.add_argument(
        "--kmer-size",
        type=int,
        default=6,
        help="K-mer size (only used if --assembly-mode dbg*).",
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
    
    # Run the pipeline with the validated arguments
    run_pipeline(args)


def run_pipeline(args):
    """
    Orchestrates the full pipeline, calling each refactored module.
    """
    
    logger.info("--- InstaNexus Pipeline started ---")
    
    run_name = Path(args.input_csv).stem # e.g., 'bsa'
    base_output_folder = Path(args.folder_outputs) / run_name # e.g., 'outputs/bsa'
    
    # Build the experiment folder name based on parameters
    folder_name_parts = [f"{args.assembly_mode}"]
    
    if args.conf is not None:
        folder_name_parts.append(f"c{args.conf}")
    
    if "dbg" in args.assembly_mode:
        folder_name_parts.append(f"ks{args.kmer_size}")
    
    folder_name_parts.append(f"mo{args.min_overlap}")
    folder_name_parts.append(f"ts{args.size_threshold}")

    if args.reference:
        folder_name_parts.extend([f"mi{args.min_identity}", f"mm{args.max_mismatches}"])

    run_folder_name = "_".join(folder_name_parts)
    experiment_folder = base_output_folder / run_folder_name # e.g., 'outputs/bsa/greedy_c0.9_mo4_ts10'
    
    cleaned_csv_path = experiment_folder / "cleaned.csv"
    
    scaffolds_folder = experiment_folder / "scaffolds"
    scaffolds_fasta_path = scaffolds_folder / "scaffolds.fasta"
    
    clustering_folder = scaffolds_folder / "clustering" # Clustering output
    cluster_fasta_folder = clustering_folder / "cluster_fasta" # Input for alignment
    
    alignment_folder = scaffolds_folder / "alignment"
    
    consensus_folder = scaffolds_folder / "consensus"

    # ID for logs (optional)
    run_id_str = f"[{run_name} @ {run_folder_name}]"
    
    logger.info(f"Starting pipeline for run: {run_id_str}")
    logger.info(f"All results will be saved to: {experiment_folder}")


    try:
        logger.info(f"--- [Step 1/5] Running Preprocessing ---")
        preprocessing.main(
            input_csv=args.input_csv,
            metadata_json=args.metadata_json_path,
            contaminants_fasta=args.contaminants_fasta_path,
            chain=args.chain,
            reference=args.reference,
            conf=args.conf,
            output_csv_path=str(cleaned_csv_path) # Pass explicit path
        )
    except Exception as e:
        logger.error(f"Preprocessing failed: {e}")
        return 

    try:
        logger.info(f"--- [Step 2/5] Running Assembly ---")
        assembly.main(
            input_csv_path=str(cleaned_csv_path),
            output_scaffolds_path=str(scaffolds_fasta_path),
            metadata_json_path=args.metadata_json_path,
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
        logger.info(f"--- [Step 3/5] Running Clustering ---")
        clustering.main(
            input_scaffolds_folder=str(scaffolds_folder), # Pass the folder
            min_seq_id=args.min_seq_id,
            coverage=args.coverage
        )
    except Exception as e:
        logger.error(f"Clustering failed: {e}")
        return
    
    try:
        logger.info(f"--- [Step 4/5] Running Alignment ---")
        alignment.main(
            input_cluster_fasta_folder=str(cluster_fasta_folder),
            output_alignment_folder=str(alignment_folder)
        )
    except Exception as e:
        logger.error(f"Alignment failed: {e}")
        return
    
    try:
        logger.info(f"--- [Step 5/5] Running Consensus ---")
        consensus.main(
            input_alignment_folder=str(alignment_folder),
            output_consensus_folder=str(consensus_folder),
            run_id=run_id_str # Pass ID for logs
        )
    except Exception as e:
        logger.error(f"Consensus failed: {e}")
        return

    logger.info(f"--- InstaNexus Pipeline finished successfully! ---")
    logger.info(f"Final results in: {experiment_folder}")


if __name__ == "__main__":
    cli()

# Example command to run the full pipeline:
# python -m instanexus.main \
#    --input-csv inputs/bsa.csv \
#    --folder-outputs outputs \
#    --metadata-json-path json/sample_metadata.json \
#    --contaminants-fasta-path fasta/contaminants.fasta \
#    --assembly-mode dbg \
#    --conf 0.9 \
#    --kmer-size 7 \
#    --size-threshold 12 \
#    --min-overlap 3 \
#    --min-seq-id 0.85