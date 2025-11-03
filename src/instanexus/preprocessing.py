#!/usr/bin/env python

r""" Preprocessing module for InstaNexus.

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
__date__ = 29 Oct 2025
__maintainer__ = Marco Reverenna
__email__ = marcor@dtu.dk
__status__ = Dev
"""

# import libraries
import os
import re
import json
import logging
import numpy as np
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go

from pathlib import Path
from Bio import SeqIO


PROJECT_ROOT = Path(__file__).resolve().parents[2]
JSON_DIR = PROJECT_ROOT / "json"

def get_sample_metadata(run, chain="", json_path=JSON_DIR / "sample_metadata.json"):
    with open(json_path, "r") as f:
        all_meta = json.load(f)

    if run not in all_meta:
        raise ValueError(f"Run '{run}' not found in metadata.")

    entries = all_meta[run]

    if not chain:
        # If no chain is specified, return the first entry
        return entries[0]

    for entry in entries:
        if entry["chain"] == chain:
            return entry

    raise ValueError(f"No metadata found for run '{run}' with chain '{chain}'.")


# Define and create the necessary directories only if they don't exist
def create_directory(path):
    """Creates a directory if it does not already exist.
    Args:
        path (str): The path of the directory to create.
    """
    if not os.path.exists(path):
        os.makedirs(path)
        # print(f"Created: {path}")
    # else:
    # print(f"Already exists: {path}")


def create_subdirectories_outputs(folder):
    """Creates subdirectories within the specified folder.
    Args:
        folder (str): The path of the parent directory.
    """
    subdirectories = ["cleaned", "contigs", "scaffolds", "statistics"]
    for subdirectory in subdirectories:
        create_directory(f"{folder}/{subdirectory}")


def normalize_sequence(sequence):
    """Normalize the given amino acid sequence by replacing all occurrences of 'I' with
    'L'.

    Parameters:
        sequence (str): The amino acid sequence to be normalized.

    Returns:
        str: The normalized amino acid sequence with 'I' replaced by 'L'.
    """
    return sequence.replace("I", "L")


def extract_protease(experiment_name, proteases):
    """Extracts the protease name from the given experiment name.

    Parameters:
        experiment_name (str): The name of the experiment.
        proteases (list or set): A list or set of known protease names.

    Returns:
        str or None: The matched protease name, or None if no match is found.
    """
    parts = experiment_name.split("_")
    for part in parts:
        if part in proteases:
            return part
    return None


def remove_modifications(psm_column):
    """Remove any content within parentheses, including the parentheses, from a given
    string. Remove UNIMOD modifications and normalize I to L.

    Parameters:
        - psm_column (str): The string containing modifications in parentheses (e.g.,
          "A(ox)BC(mod)D"). If the value is null, it returns None.

    Returns:
        - str: The string with all parenthetical modifications removed (e.g., "ABCD"),
          or None if the input was null.
    """
    if pd.notnull(psm_column):
        ret = re.sub(r"\(.*?\)", "", psm_column)
        ret = re.sub(r"\[.*?\]", "", ret)  # replace UNIMOD modifications
        ret = normalize_sequence(ret)
        return ret
    return None


# ! needs to move once it is a package
def test_remove_modifications():
    assert remove_modifications("A(ox)BC(mod)D") == "ABCD"
    assert remove_modifications("A[UNIMOD:21]BC[UNIMOD:35]D") == "ABCD"
    assert (
        remove_modifications("A(ox)[UNIMOD:21]BC(mod)[UNIMOD:35]D") == "ABCD"
    )
    assert remove_modifications(None) is None
    assert remove_modifications("ACD") == "ACD"
    assert remove_modifications("A(I)BCD") == "ABCD"
    assert remove_modifications("A(ox)B(I)C(mod)D") == "ABCD"
    assert (
        remove_modifications("A(ox)[UNIMOD:21]B(I)C(mod)[UNIMOD:35]D") == "ABCD"
    )
    assert remove_modifications("AI BCD") == "AL BCD"
    assert remove_modifications("A(ox)I B(mod)CD") == "AL BCD"


def clean_dataframe(df):
    """Clean and preprocess a DataFrame for analysis by removing '(ox)' substrings
    from sequences in the 'seq' column. By replacing values of -1 with -10 in the
    'log_probs' column, by dropping rows with missing values in the 'preds' column.
    By extracts a 'protease' value from the 'experiment_name' column based on a
    specific naming convention. By adding a 'conf' column, which is the exponentiated
    'log_probs' to represent confidence and sorting the DataFrame by the 'conf' column
    in descending order.

    Parameters:
        - df (DataFrame): The raw input DataFrame to clean.

    Returns:
        - DataFrame: The cleaned and processed DataFrame.
    """
    df = df.copy()

    df["log_probs"] = df["log_probs"].replace(-1, -10)
    # -10 is very low, replacing with -1 (so that we are sure is very low quality)

    df = df.dropna(subset=["preds"])

    df.loc[:, "conf"] = np.exp(df["log_probs"])
    df = df.sort_values("conf", ascending=False)

    return df


def filter_contaminants(seqs, run, contaminants_fasta):
    """Filters out sequences from the input list `seqs` that are substrings of
    sequences in the contaminants file. If run == 'bsa', the Bovine serum albumin
    precursor is ignored.

    Parameters:
        - seqs (list of str): List of sequences to be filtered.
        - contaminants_fasta (str): Path to the FASTA file containing contaminant
          sequences.
        - run (str): Run identifier, used to control special filtering logic.
    """
    contam_records = []
    for record in SeqIO.parse(contaminants_fasta, "fasta"):
        if run == "bsa" and "Bovine serum albumin precursor" in record.description:
            continue  # skip BSA if run is 'bsa'
        contam_records.append(str(record.seq))

    filtered_seqs = []
    removed_count = 0

    for seq in seqs:
        if any(seq in contam_seq for contam_seq in contam_records):
            removed_count += 1
        else:
            filtered_seqs.append(seq)

    return filtered_seqs


logging.basicConfig(
    level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s"
)
logger = logging.getLogger(__name__)



def main(
    input_csv: str,
    chain: str = "",
    folder_outputs: str = "outputs",
    reference: bool = False,
    assembly_mode = "dbg",
    conf: float = 0.88,
    kmer_size: int = 6,
    size_threshold: int = 10,
    min_overlap: int = 4,
    min_identity: float = 0.8,
    max_mismatches: int = 14,
):
    """Main function to run the preprocessing script."""
    input_csv = Path(input_csv)

    print("Starting preprocessing pipeline.")

    run = input_csv.stem
    repo_folder = Path(__file__).resolve().parents[2]

    # load metadata
    if chain:
        meta = get_sample_metadata(run, chain=chain)
    else:
        meta = get_sample_metadata(run)

    proteases = meta["proteases"]

    if reference:
        protein = meta["protein"]

    if assembly_mode != "dbg":
        print("Ignoring kmer_size (only relevant for dbg mode)")
        kmer_size = None

    if not reference:
        print("Ignoring min_identity and max_mismatches (only relevant when reference=True)")
        min_identity = None
        max_mismatches = None

    print("Parameters loaded.")

    folder_outputs = Path(folder_outputs) / run
    folder_outputs.mkdir(parents=True, exist_ok=True)

    folder_name_parts = [f"comb_{assembly_mode}", f"c{conf}", f"ts{size_threshold}", f"mo{min_overlap}"]

    if assembly_mode == "dbg":
        folder_name_parts.insert(2, f"ks{kmer_size}")

    if reference:
        folder_name_parts.extend([f"mi{min_identity}", f"mm{max_mismatches}"])

    combination_folder_out = folder_outputs / "_".join(folder_name_parts)
    create_subdirectories_outputs(combination_folder_out)

    print(f"Output folders created at: {combination_folder_out}")

    # data cleaning

    logger.info("Starting data cleaning...")
    
    if reference:
        protein_norm = normalize_sequence(protein)
    df = pd.read_csv(input_csv)

    df["protease"] = df["experiment_name"].apply(
        lambda name: extract_protease(name, proteases)
    )
    
    df = clean_dataframe(df)
    
    df["cleaned_preds"] = df["preds"].apply(remove_modifications)
    
    cleaned_psms = df["cleaned_preds"].tolist()

    filtered_psms = filter_contaminants(
        cleaned_psms, run, repo_folder / "fasta/contaminants.fasta"
    )
    df = df[df["cleaned_preds"].isin(filtered_psms)]

    if reference:
        df["mapped"] = df["cleaned_preds"].apply(
            lambda x: "True" if x in protein_norm else "False"
        )

    # probably confidence trhreshold won't be necessary anymore
    df = df[df["conf"] > conf]
    
    df.reset_index(drop=True, inplace=True)
    
    logger.info("Data cleaning completed.")

    cleaned_csv_path = combination_folder_out / "cleaned" / "cleaned_data.csv"
    
    df.to_csv(cleaned_csv_path, index=False)

    logger.info("Cleaned data saved to: {}.".format(cleaned_csv_path))


def cli():
    """Command-line interface for the preprocessing script."""
    import argparse

    parser = argparse.ArgumentParser(
        description="Preprocess peptide-spectrum match (PSM) data."
    )
    parser.add_argument(
        "--input-csv",
        type=str,
        required=True,
        help="Path to the input CSV file containing PSM data.",
    )
    parser.add_argument(
        "--chain",
        type=str,
        default="",
        help="Chain identifier for the sample (optional).",
    )
    parser.add_argument(
        "--folder-outputs",
        type=str,
        default="outputs",
        help="Folder to save output files.",
    )
    parser.add_argument(
        "--reference",
        action="store_true",
        help="Whether to use reference protein sequence for mapping.",
    )

    parser.add_argument(
        "--assembly-mode",
        type=str,
        choices=["dbg", "greedy"],
        help="Assembly algorithm to use.",
    )
    parser.add_argument(
        "--kmer-size",
        type=int,
        default=6,
        help="K-mer size (only used if --assembly-mode dbg).",
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
        "--conf",
        type=float,
        default=0.88,
        help="Confidence threshold for filtering (default: 0.88).",
    )
    parser.add_argument(
        "--size-threshold",
        type=int,
        default=10,
        help="Minimum contig size threshold (default: 10).",
    )
    parser.add_argument(
        "--min-overlap",
        type=int,
        default=4,
        help="Minimum overlap size between reads (default: 4).",
    )

    args = parser.parse_args()

    main(
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


if __name__ == "__main__":
    cli()

# python -m instanexus.preprocessing --input-csv ../../inputs/bsa.csv --folder-outputs ../../outputs --assembly-mode dbg --conf 0.9 --kmer-size 7 --size-threshold 12 --min-overlap 5