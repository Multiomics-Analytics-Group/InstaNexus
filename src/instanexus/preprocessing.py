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
import argparse
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


def strip_modifications(psm_column):
    """The bare residue letters of a prediction, with I and L left as they were read.

    Assembly works on the normalized form, where I and L are one residue, because two
    reads of the same region that spell it differently do not overlap on a literal
    comparison. The spelling is still what the model predicted, though, and the
    scaffold is the output of this tool, so the residues are kept here for the vote
    that puts isoleucine back (see ``assembly.restore_isoleucine``).
    """
    if pd.notnull(psm_column):
        # remove (ox) style modifications
        ret = re.sub(r"\(.*?\)-?", "", str(psm_column))
        # replace UNIMOD modifications
        ret = re.sub(r"\[.*?\]-?", "", ret)

        return ret.strip("-")
    return None


def remove_modifications(psm_column):
    """The bare residue letters, normalized so that I and L are one residue.

    This is what assembly overlaps on. ``strip_modifications`` is the same thing
    without the normalization.
    """
    stripped = strip_modifications(psm_column)
    return normalize_sequence(stripped) if stripped is not None else None


def clean_instanovo_raw(df):
    """V1 Strategy: 'preds' -> 'cleaned_preds', 'log_probs' -> 'conf'."""
    logger.info("Detected V1 format (Standard).")

    if "log_probs" in df.columns:
        df["log_probs"] = df["log_probs"].replace(-1, -10)
        df.loc[:, "conf"] = np.exp(df["log_probs"])

    if "preds" in df.columns:
        df = df.dropna(subset=["preds"])
        df["cleaned_preds"] = df["preds"].apply(remove_modifications)
        df["read_residues"] = df["preds"].apply(strip_modifications)
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
        df["read_residues"] = df[found_seq_col].apply(strip_modifications)
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
        "read_residues",
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


def filter_by_length(df: pd.DataFrame, min_length: Optional[int], max_length: Optional[int]) -> pd.DataFrame:
    """Keep peptides whose length, in residues after removing modifications, is within the bounds.

    Args:
        df: PSMs with a `cleaned_preds` column.
        min_length: Shortest peptide kept, inclusive; None for no lower bound.
        max_length: Longest peptide kept, inclusive; None for no upper bound.

    Returns:
        The rows within the bounds.
    """
    if min_length is not None and max_length is not None and min_length > max_length:
        raise ValueError(f"min_length ({min_length}) is greater than max_length ({max_length}).")
    lengths = df["cleaned_preds"].str.len()
    keep = pd.Series(True, index=df.index)
    if min_length is not None:
        keep &= lengths >= min_length
    if max_length is not None:
        keep &= lengths <= max_length

    return df[keep]


def read_support(df: pd.DataFrame) -> pd.Series:
    """Count how often each PSM's peptide is read in the sample.

    Peptides are compared after removing modifications and with I and L treated as one residue, so the count
    is taken over `cleaned_preds`.

    Args:
        df: PSMs with a `cleaned_preds` column.

    Returns:
        For each row, the number of rows with the same peptide.
    """
    return df.groupby("cleaned_preds")["cleaned_preds"].transform("size")


def top_fraction_mask(df: pd.DataFrame, fraction: float) -> pd.Series:
    """Select the most confident fraction of the PSMs.

    PSMs are ranked by `conf`, ties in the order of the rows, and the first `floor(fraction * n)` are kept,
    where `n` counts the PSMs with a confidence. PSMs without a confidence are not kept.

    Args:
        df: PSMs with a `conf` column.
        fraction: Share of the PSMs to keep, in (0, 1].

    Returns:
        True for the rows to keep.
    """
    if not 0 < fraction <= 1:
        raise ValueError(f"fraction ({fraction}) must be in (0, 1].")
    rank = df["conf"].rank(ascending=False, method="first", pct=True)

    return rank <= fraction


def _filter_step(name: str, applied: bool, rows_in: int, rows_out: int) -> dict:
    """One entry of the `steps` list in `<output stem>_filters.json`.

    Args:
        name: Filter name: cleaning, length, contaminants, conf, fdr, min_reads or top_fraction.
        applied: Whether the filter was requested and could run.
        rows_in: Rows before the filter.
        rows_out: Rows after the filter.

    Returns:
        The entry.
    """
    return {"name": name, "applied": applied, "rows_in": int(rows_in), "rows_out": int(rows_out)}


def filters_record_path(output_csv_path: Path) -> Path:
    """Return the path of the JSON file that records the filters applied to a preprocessing output.

    Args:
        output_csv_path: The cleaned CSV written by preprocessing.

    Returns:
        `<output stem>_filters.json`, next to the CSV.
    """
    return output_csv_path.with_name(f"{output_csv_path.stem}_filters.json")


def _peptide_length(value: str) -> int:
    """Parse --min-length / --max-length: a positive integer.

    Args:
        value: Command-line value.

    Returns:
        The length.
    """
    length = int(value)
    if length < 1:
        raise argparse.ArgumentTypeError("must be a positive integer")

    return length


def add_length_filter_arguments(parser: argparse.ArgumentParser) -> None:
    """Add --min-length and --max-length to a command-line parser.

    Args:
        parser: Parser to extend.
    """
    parser.add_argument(
        "--min-length",
        type=_peptide_length,
        default=None,
        help="Keep only peptides of at least this many residues, after removing modifications (default: no limit).",
    )
    parser.add_argument(
        "--max-length",
        type=_peptide_length,
        default=None,
        help="Keep only peptides of at most this many residues, after removing modifications (default: no limit).",
    )


def _top_fraction(value: str) -> float:
    """Parse --top-fraction: a number in (0, 1].

    Args:
        value: Command-line value.

    Returns:
        The fraction.
    """
    fraction = float(value)
    if not 0 < fraction <= 1:
        raise argparse.ArgumentTypeError("must be greater than 0 and at most 1")

    return fraction


def _min_reads(value: str) -> int:
    """Parse --min-reads: a positive integer.

    Args:
        value: Command-line value.

    Returns:
        The number of reads.
    """
    reads = int(value)
    if reads < 1:
        raise argparse.ArgumentTypeError("must be a positive integer")

    return reads


def add_rank_and_support_arguments(parser: argparse.ArgumentParser) -> None:
    """Add --top-fraction and --min-reads to a command-line parser.

    Args:
        parser: Parser to extend.
    """
    parser.add_argument(
        "--top-fraction",
        type=_top_fraction,
        default=None,
        help=(
            "Keep only this fraction of the sample's PSMs, the most confident, ranked after the conf and FDR "
            "filters (default: keep all)."
        ),
    )
    parser.add_argument(
        "--min-reads",
        type=_min_reads,
        default=None,
        help=(
            "Keep only PSMs whose peptide, with modifications removed and I = L, is read at least this many "
            "times in the sample, counted after the conf and FDR filters (default: no minimum)."
        ),
    )


def check_length_bounds(parser: argparse.ArgumentParser, args: argparse.Namespace) -> None:
    """Fail with a usage error when --min-length is greater than --max-length.

    Args:
        parser: The parser that produced args.
        args: Parsed arguments.
    """
    if args.min_length is not None and args.max_length is not None and args.min_length > args.max_length:
        parser.error(f"--min-length ({args.min_length}) is greater than --max-length ({args.max_length}).")


def main(
    input_csv: str,
    metadata_json: Optional[str],
    contaminants_fasta: Optional[str],
    chain: str,
    reference: bool,
    conf: float,
    fdr: float,
    output_csv_path: str,
    min_length: Optional[int] = None,
    max_length: Optional[int] = None,
    top_fraction: Optional[float] = None,
    min_reads: Optional[int] = None,
):
    """Main function to run the preprocessing script.

    Besides the cleaned CSV, writes `<output stem>_filters.json` next to it, recording the filters applied
    and, in `steps`, the rows before and after each of them, in the order they run: cleaning, length,
    contaminants, conf, FDR, min_reads, top_fraction. The read counts and the confidence rank are both taken
    over the PSMs left after FDR, so with both set a PSM is kept when it passes each of them.
    """
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

    seq_candidates = ["preds", "prediction_untokenised", "prediction", "Peptide", "sequence"]
    seq_col = next((c for c in seq_candidates if c in df.columns), None)
    if seq_col is not None:
        df["cleaned_preds"] = df[seq_col].apply(remove_modifications)
    else:
        raise ValueError("No suitable column found for peptide sequences.")

    # rows before and after each filter, in the order they are applied
    steps = []

    rows_in = len(df)
    df = clean_dataframe(df)
    steps.append(_filter_step("cleaning", True, rows_in, len(df)))

    rows_in = len(df)
    length_applied = min_length is not None or max_length is not None
    if length_applied:
        logger.info(f"Applying length filter: {min_length or 1} to {max_length or 'any'} residues")
        df = filter_by_length(df, min_length, max_length)
    steps.append(_filter_step("length", length_applied, rows_in, len(df)))

    rows_in = len(df)
    if contaminants_fasta is not None:
        cleaned_psms = df["cleaned_preds"].tolist()
        filtered_psms = filter_contaminants(cleaned_psms, run, contaminants_fasta)
        df = df[df["cleaned_preds"].isin(filtered_psms)]
    steps.append(_filter_step("contaminants", contaminants_fasta is not None, rows_in, len(df)))

    if reference and protein_norm is not None:
        df["mapped"] = df["cleaned_preds"].apply(lambda x: True if x in protein_norm else False)

    rows_in = len(df)
    conf_applied = conf is not None and "conf" in df.columns
    if conf is not None:
        logger.info(f"Applying confidence threshold: > {conf}")
        if "conf" in df.columns:
            df = df[df["conf"] > conf]
    steps.append(_filter_step("conf", conf_applied, rows_in, len(df)))

    rows_in = len(df)
    fdr_applied = fdr is not None and "psm_q_value" in df.columns
    if fdr is not None:
        logger.info(f"Applying FDR threshold: <= {fdr}")
        if "psm_q_value" in df.columns:
            df = df[df["psm_q_value"] <= fdr]
        else:
            logger.warning("FDR filter requested but 'psm_q_value' column not found (V1 data?).")
    steps.append(_filter_step("fdr", fdr_applied, rows_in, len(df)))

    if conf is None and fdr is None:
        logger.info("No filters (Conf/FDR) applied.")

    # counted over the same PSMs, so the order of these two filters does not change the result
    after_fdr = df
    rows_in = len(df)
    if min_reads is not None:
        logger.info(f"Applying read support: peptides read at least {min_reads} times")
        df = df[read_support(after_fdr) >= min_reads]
    steps.append(_filter_step("min_reads", min_reads is not None, rows_in, len(df)))

    rows_in = len(df)
    top_fraction_applied = top_fraction is not None and "conf" in df.columns
    if top_fraction is not None:
        logger.info(f"Applying rank cutoff: the most confident {top_fraction} of PSMs")
        if "conf" in df.columns:
            df = df[top_fraction_mask(after_fdr, top_fraction).loc[df.index]]
        else:
            logger.warning("Top-fraction filter requested but 'conf' column not found.")
    steps.append(_filter_step("top_fraction", top_fraction_applied, rows_in, len(df)))

    input_dir = input_csv_path.parent
    df = add_quantification_data(df, run, inputs_folder=input_dir)

    df.reset_index(drop=True, inplace=True)

    logger.info("Data cleaning completed.")
    cleaned_csv_path = Path(output_csv_path)
    cleaned_csv_path.parent.mkdir(parents=True, exist_ok=True)

    df.to_csv(cleaned_csv_path, index=False)
    logger.info("Cleaned data saved to: {}.".format(cleaned_csv_path))

    record = {
        "input_csv": input_csv_path.name,
        "min_length": min_length,
        "max_length": max_length,
        "conf": conf,
        "fdr": fdr,
        "top_fraction": top_fraction,
        "min_reads": min_reads,
        "contaminants_fasta": contaminants_fasta,
        "steps": steps,
        "rows_output": len(df),
    }
    with open(filters_record_path(cleaned_csv_path), "w") as f:
        json.dump(record, f, indent=4)


def cli():
    """Command-line interface for the preprocessing script."""
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
    add_length_filter_arguments(parser)
    add_rank_and_support_arguments(parser)

    args = parser.parse_args()
    check_length_bounds(parser, args)

    main(**vars(args))


if __name__ == "__main__":
    cli()

# python -m instanexus.preprocessing --input-csv inputs/bsa.csv --metadata-json json/sample_metadata.json --contaminants-fasta fasta/contaminants.fasta --output-csv-path outputs/bsa_cleaned.csv --conf 0.9 --reference

# winnow
# python -m instanexus.preprocessing --input-csv inputs/bsa.csv --metadata-json json/sample_metadata.json --contaminants-fasta fasta/contaminants.fasta --output-csv-path outputs/bsa_cleaned.csv --fdr 0.02
