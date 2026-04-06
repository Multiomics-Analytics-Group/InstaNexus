#!/usr/bin/env python

r"""Preprocessing module for InstaNexus.

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
__date__ = 20 Dec 2025
__maintainer__ = Marco Reverenna
__email__ = marcor@dtu.dk
__status__ = Dev
"""

# import libraries
import json
import logging
import re
from pathlib import Path
from typing import Optional

import numpy as np
import pandas as pd
from Bio import SeqIO

# PROJECT_ROOT = Path(__file__).resolve().parents[2]
# JSON_DIR = PROJECT_ROOT / "json"


def get_sample_metadata(run, chain="", json_path=None):
    """Retrieve sample metadata from a JSON file based on the run and optional chain."""
    if json_path is None:
        raise ValueError("json_path must be provided.")

    with open(json_path, "r") as f:
        all_meta = json.load(f)

    if run not in all_meta:
        raise ValueError(f"Run '{run}' not found in metadata.")

    entries = all_meta[run]

    if len(entries) == 1:
        return entries[0]

    if not chain:
        available_chains = [e.get("chain", "unknown") for e in entries]
        raise ValueError(
            f"Multiple entries found for run '{run}'. Available chains: {available_chains}. Please specify a chain."
        )

    for entry in entries:
        if entry["chain"] == chain:
            return entry

    raise ValueError(f"No metadata found for run '{run}' with chain '{chain}'.")


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
    """
    if pd.notnull(psm_column):
        # remove (ox) style modifications
        ret = re.sub(r"\(.*?\)-?", "", str(psm_column))
        # replace UNIMOD modifications
        ret = re.sub(r"\[.*?\]-?", "", ret)

        ret = ret.strip("-")
        ret = normalize_sequence(ret)
        return ret
    return None


def clean_instanovo_raw(df):
    """V1 Strategy: 'preds' -> 'cleaned_preds', 'log_probs' -> 'conf'."""
    logger.info("Detected V1 format (Standard).")

    if "log_probs" in df.columns:
        df["log_probs"] = df["log_probs"].replace(-1, -10)
        df.loc[:, "conf"] = np.exp(df["log_probs"])

    if "preds" in df.columns:
        df = df.dropna(subset=["preds"])
        df["cleaned_preds"] = df["preds"].apply(remove_modifications)
    else:
        logger.error("V1 Error: 'preds' column missing!")
        return pd.DataFrame()  # Empty DF on error

    if "conf" in df.columns:
        df = df.sort_values("conf", ascending=False)
    return df


def clean_winnow_rescored(df):
    """
    V2 Strategy (Robust):
    - Detects sequence column flexibly ('prediction_untokenised', 'prediction', etc.) -> 'cleaned_preds'
    - Maps 'calibrated_confidence' -> 'conf'
    - Handles missing columns gracefully (assigns None or skips).
    """
    logger.info("Detected V2 format (Winnow/Updated). Cleaning aggressively...")

    if "calibrated_confidence" in df.columns:
        df = df.rename(columns={"calibrated_confidence": "conf"})

    seq_candidates = ["prediction_untokenised", "prediction", "Peptide", "sequence", "cleaned_preds"]
    found_seq_col = None

    for col in seq_candidates:
        if col in df.columns:
            found_seq_col = col
            break

    if found_seq_col:
        logger.info(f"Using '{found_seq_col}' as sequence column.")
        df = df.dropna(subset=[found_seq_col])
        df["cleaned_preds"] = df[found_seq_col].apply(remove_modifications)
    else:
        logger.error("V2 Critical Error: No sequence column found (checked: prediction_untokenised, prediction, etc.)!")
        return pd.DataFrame()

    optional_cols = ["instanovo_token_log_probabilities", "ion_match_intensity", "iRT", "iRT error"]
    for col in optional_cols:
        if col not in df.columns:
            logger.warning(f"Optional column '{col}' missing. Filling with None.")
            df[col] = None

    cols_of_interest = [
        "experiment_name",
        "scan_number",
        "cleaned_preds",
        "conf",
        "psm_q_value",
        "delta_mass_ppm",
        "Mass Error",
        "is_missing_prosit_features",
        "instanovo_token_log_probabilities",
        "ion_match_intensity",
        "ion_matches",
        "iRT",
        "iRT error",
        "is_missing_irt_error",
        "predicted iRT",
        "margin",
        "entropy",
        "z-score",
        "correct",
        "protease",
        "mapped",
    ]

    existing_cols = [c for c in cols_of_interest if c in df.columns]
    df = df[existing_cols]

    if "psm_q_value" in df.columns:
        df = df.sort_values("psm_q_value", ascending=False)

    return df


def clean_dataframe(df):
    """
    Wrapper function that detects the data format
    and applies the correct cleaning strategy.
    """
    df = df.copy()

    # Auto-detection logic
    if "calibrated_confidence" in df.columns or "prediction_untokenised" in df.columns:
        return clean_winnow_rescored(df)

    elif "log_probs" in df.columns:
        return clean_instanovo_raw(df)

    else:
        logger.warning("Unknown data format! Attempting generic cleaning via V2 logic.")
        try:
            return clean_winnow_rescored(df)
        except Exception as e:
            logger.error(f"Generic cleaning failed: {e}. Returning dataframe as is.")
            return df


# def add_quantification_data(df_main, run_name, fdr_threshold, inputs_folder="inputs"):
#     """
#     Filters df_main by FDR, then looks for a quantification file ({run_name}_quant_scores.csv).
#     Merges the abundance data into the filtered dataframe.
#     """
#     if fdr_threshold is not None:
#         if "psm_q_value" in df_main.columns:
#             initial_len = len(df_main)
#             df_main = df_main[df_main['psm_q_value'] <= fdr_threshold].copy()
#             logger.info(f"FDR Filter applied inside merge function: {initial_len} -> {len(df_main)} rows (<= {fdr_threshold})")
#         else:
#             logger.warning("FDR threshold provided but 'psm_q_value' column missing. Skipping filter.")

#     # 2. Gestione Path
#     quant_file_name = f"{run_name}_quant_scores.csv"
#     quant_file_path = Path(inputs_folder) / quant_file_name

#     if not quant_file_path.exists():
#         logger.warning(f"Quantification file NOT FOUND: {quant_file_path}")
#         logger.warning("Skipping abundance merging. 'peptide_abundance' will be missing.")
#         return df_main

#     logger.info(f"Found quantification file: {quant_file_path}")

#     try:
#         df_quant = pd.read_csv(quant_file_path)

#         if "cleaned_preds" not in df_quant.columns or "total_abundance_norm" not in df_quant.columns:
#             logger.warning(f"Quantification file format error. Missing columns in {quant_file_path}")
#             return df_main

#         df_quant_summed = df_quant.groupby('cleaned_preds', as_index=False)['total_abundance_norm'].sum()
#         df_quant_summed.rename(columns={'total_abundance_norm': 'peptide_abundance'}, inplace=True)
#         df_merged = pd.merge(df_main, df_quant_summed, on='cleaned_preds', how='left')
#         df_merged['peptide_abundance'] = df_merged['peptide_abundance'].fillna(0)

#         logger.info(f"Quantification data merged successfully. Output rows: {len(df_merged)}")
#         return df_merged

#     except Exception as e:
#         logger.error(f"Error merging quantification data: {e}")
#         return df_main


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


def add_quantification_data(df_main, run_name, inputs_folder="inputs"):
    quant_file_path = Path(inputs_folder) / f"{run_name}_quant_scores.csv"

    if not quant_file_path.exists():
        return df_main

    try:
        df_quant = pd.read_csv(quant_file_path)

        required_cols = ["cleaned_preds", "total_abundance_norm"]
        if not all(col in df_quant.columns for col in required_cols):
            return df_main

        valid_peptides = set(df_main["cleaned_preds"].unique())
        df_quant = df_quant[df_quant["cleaned_preds"].isin(valid_peptides)]

        df_quant_summed = df_quant.groupby("cleaned_preds", as_index=False)["total_abundance_norm"].sum()
        df_quant_summed.rename(columns={"total_abundance_norm": "peptide_abundance"}, inplace=True)

        df_merged = pd.merge(df_main, df_quant_summed, on="cleaned_preds", how="left")
        df_merged["peptide_abundance"] = df_merged["peptide_abundance"].fillna(0)

        logger.info(f"Quantification merged. Rows with abundance: {len(df_quant_summed)}")
        return df_merged

    except Exception as e:
        logger.error(f"Error merging quantification: {e}")
        return df_main


logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
logger = logging.getLogger(__name__)


def main(
    input_csv: str,
    metadata_json: Optional[str],
    contaminants_fasta: Optional[str],
    chain: str,
    reference: bool,
    conf: float,
    fdr: float,
    output_csv_path: str,
):
    """Main function to run the preprocessing script."""
    input_csv_path = Path(input_csv)

    run = input_csv_path.stem  # stem gives the filename without suffix

    protein_norm = None
    if metadata_json is not None:
        if chain:
            meta = get_sample_metadata(run, chain=chain, json_path=metadata_json)
        else:
            meta = get_sample_metadata(run, json_path=metadata_json)
        proteases = meta["proteases"]
        if reference:
            protein = meta["protein"]
            protein_norm = normalize_sequence(protein)

    df = pd.read_csv(input_csv_path)

    if metadata_json is not None and "experiment_name" in df.columns:
        df["protease"] = df["experiment_name"].apply(lambda name: extract_protease(name, proteases))

    if "preds" in df.columns:
        df["cleaned_preds"] = df["preds"].apply(remove_modifications)
    elif "prediction_untokenised" in df.columns:
        df["cleaned_preds"] = df["prediction_untokenised"].apply(remove_modifications)
    else:
        raise ValueError("No suitable column found for peptide sequences.")

    df = clean_dataframe(df)

    if contaminants_fasta is not None:
        cleaned_psms = df["cleaned_preds"].tolist()
        filtered_psms = filter_contaminants(cleaned_psms, run, contaminants_fasta)
        df = df[df["cleaned_preds"].isin(filtered_psms)]

    if reference and protein_norm is not None:
        df["mapped"] = df["cleaned_preds"].apply(lambda x: True if x in protein_norm else False)

    if conf is not None:
        logger.info(f"Applying confidence threshold: > {conf}")
        if "conf" in df.columns:
            df = df[df["conf"] > conf]

    if fdr is not None:
        logger.info(f"Applying FDR threshold: <= {fdr}")
        if "psm_q_value" in df.columns:
            df = df[df["psm_q_value"] <= fdr]
        else:
            logger.warning("FDR filter requested but 'psm_q_value' column not found (V1 data?).")

    if conf is None and fdr is None:
        logger.info("No filters (Conf/FDR) applied.")

    input_dir = input_csv_path.parent
    df = add_quantification_data(df, run, inputs_folder=input_dir)

    df.reset_index(drop=True, inplace=True)

    logger.info("Data cleaning completed.")
    cleaned_csv_path = Path(output_csv_path)
    cleaned_csv_path.parent.mkdir(parents=True, exist_ok=True)

    df.to_csv(cleaned_csv_path, index=False)
    logger.info("Cleaned data saved to: {}.".format(cleaned_csv_path))


def cli():
    """Command-line interface for the preprocessing script."""
    import argparse

    parser = argparse.ArgumentParser(description="Preprocess peptide-spectrum match (PSM) data.")
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
        "--reference",
        action="store_true",
        help="Whether to use reference protein sequence for mapping.",
    )
    parser.add_argument(
        "--conf",
        type=float,
        default=None,
        help="Confidence threshold for filtering (default: 0.88).",
    )
    parser.add_argument(
        "--fdr",
        type=float,
        default=None,
        help="False discovery rate threshold for filtering (default: 0.1).",
    )
    parser.add_argument(
        "--metadata-json",
        type=str,
        default=None,
        help="Path to the sample_metadata.json file (optional).",
    )
    parser.add_argument(
        "--contaminants-fasta",
        type=str,
        default=None,
        help="Path to the contaminants.fasta file (optional).",
    )
    parser.add_argument(
        "--output-csv-path",
        type=str,
        required=True,
        help="Path to the output CSV file.",
    )

    args = parser.parse_args()

    main(**vars(args))


if __name__ == "__main__":
    cli()

# python -m instanexus.preprocessing --input-csv inputs/bsa.csv --metadata-json json/sample_metadata.json --contaminants-fasta fasta/contaminants.fasta --output-csv-path outputs/bsa_cleaned.csv --conf 0.9 --reference

# winnow
# python -m instanexus.preprocessing --input-csv inputs/bsa.csv --metadata-json json/sample_metadata.json --contaminants-fasta fasta/contaminants.fasta --output-csv-path outputs/bsa_cleaned.csv --fdr 0.02
