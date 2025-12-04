#!/usr/bin/env python

r"""FDR analysis script: Aggregated Categories.

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
__date__ = 01 Dec 2025
__maintainer__ = Marco Reverenna
__email__ = marcor@dtu.dk
__status__ = Dev
"""

import json
import logging
import os
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

from instanexus import helpers, preprocessing, visualization
from instanexus.assembly import Assembler

SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = SCRIPT_DIR.parent
if str(PROJECT_ROOT) not in sys.path:
    sys.path.append(str(PROJECT_ROOT))

INPUTS_FOLDER = PROJECT_ROOT / "inputs"
METADATA_JSON = PROJECT_ROOT / "json/sample_metadata.json"
COLORS_JSON = PROJECT_ROOT / "json" / "colors.json"
CONTAMINANTS_FASTA = PROJECT_ROOT / "fasta/contaminants.fasta"
OUTPUT_FOLDER = PROJECT_ROOT / "outputs" / "_fdr_analysis"

# Define sample groups
SAMPLE_GROUPS = {
    "BSA": ["bsa.csv"],
    "Nanobodies": [f"nb{i}.csv" for i in range(1, 11)],  # Generates nb1..nb10
    "Antibodies": ["ma1.csv", "ma2.csv", "ma3.csv"],
    "Binders": ["bind1.csv", "bind2.csv", "bind3.csv"],
}

CATEGORY_COLOR_MAP = {
    "BSA": "bsa",
    "Nanobodies": "nanobodies",
    "Antibodies": "antibodies",
    "Binders": "binders",
}

FDR_THRESHOLDS = [0.01, 0.05, 0.10, 0.20, 0.40]
ASSEMBLY_MODE = "dbg_weighted"
KMER_SIZE = 7
MIN_OVERLAP = 3
SIZE_THRESHOLD = 10
MIN_IDENTITY = 0.8
MAX_MISMATCHES = 100

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(message)s")
logger = logging.getLogger(__name__)


def add_quantification_data(df_main, run_name, inputs_folder):
    """Merges peptide_abundance from quant files if available."""
    quant_file_name = f"{run_name}_quant_scores.csv"
    quant_file_path = Path(inputs_folder) / quant_file_name

    if not quant_file_path.exists():
        return df_main

    try:
        df_quant = pd.read_csv(quant_file_path)
        valid_peps = set(df_main["cleaned_preds"].unique())
        df_quant = df_quant[df_quant["cleaned_preds"].isin(valid_peps)]

        df_sum = df_quant.groupby("cleaned_preds", as_index=False)["total_abundance_norm"].sum()
        df_sum.rename(columns={"total_abundance_norm": "peptide_abundance"}, inplace=True)

        df_merged = pd.merge(df_main, df_sum, on="cleaned_preds", how="left")
        df_merged["peptide_abundance"] = df_merged["peptide_abundance"].fillna(0)
        return df_merged
    except Exception:
        return df_main


def load_custom_palette():
    default_palette = "viridis"
    if not os.path.exists(COLORS_JSON):
        return default_palette
    try:
        with open(COLORS_JSON, "r") as f:
            colors_data = json.load(f)
        custom_palette = {}
        for category_label, json_key in CATEGORY_COLOR_MAP.items():
            color = colors_data.get(json_key, {}).get("scaffold", "#333333")
            custom_palette[category_label] = color
        return custom_palette
    except Exception:
        return default_palette


def main():
    os.makedirs(OUTPUT_FOLDER, exist_ok=True)

    # 1. Load the full metadata JSON to handle it manually
    with open(METADATA_JSON, "r") as f:
        FULL_METADATA = json.load(f)

    assembler = Assembler(
        mode=ASSEMBLY_MODE,
        kmer_size=KMER_SIZE,
        min_overlap=MIN_OVERLAP,
        size_threshold=SIZE_THRESHOLD,
        min_identity=MIN_IDENTITY,
        max_mismatches=MAX_MISMATCHES,
    )

    all_results = []

    for category, file_list in SAMPLE_GROUPS.items():
        logger.info(f"=== Processing Category: {category} ({len(file_list)} samples) ===")

        for filename in file_list:
            csv_path = INPUTS_FOLDER / filename
            run_name = csv_path.stem
            clean_run_name = run_name.replace("_cleaned", "")

            logger.info(f"   -> Sample File: {run_name}")

            if not csv_path.exists():
                logger.warning(f"File not found: {csv_path}. Skipping.")
                continue

            # 2. Retrieve the LIST of entries (to handle antibodies with 2 chains)
            # We access the JSON directly instead of using the strict get_sample_metadata
            meta_entries = FULL_METADATA.get(clean_run_name, [])
            if isinstance(meta_entries, dict):
                meta_entries = [meta_entries]  # Normalize to list if it is a single dict

            if not meta_entries:
                logger.warning(f"Metadata not found for {clean_run_name}. Skipping.")
                continue

            # Read CSV once per Run
            df_original = pd.read_csv(csv_path)

            # 3. Iterate over EACH target (Chain) available for this Run
            for meta in meta_entries:
                target_protein = meta.get("protein", "")
                chain_type = meta.get("chain", "")
                proteases = meta.get("proteases", [])

                # Create a unique label for the report (e.g., "ma1 (heavy)")
                sample_label = f"{clean_run_name} ({chain_type})" if chain_type else clean_run_name

                logger.info(f"      Target: {sample_label}")

                protein_norm = preprocessing.normalize_sequence(target_protein)

                # Work on a copy of the dataframe
                df = df_original.copy()

                if "experiment_name" in df.columns:
                    # Use default argument binding (p=proteases) to fix the variable in lambda
                    df["protease"] = df["experiment_name"].apply(
                        lambda x, p=proteases: preprocessing.extract_protease(x, p)
                    )

                df = preprocessing.clean_dataframe(df)

                if "cleaned_preds" in df.columns:
                    df["cleaned_preds"] = df["cleaned_preds"].apply(preprocessing.remove_modifications)
                    df = df.dropna(subset=["cleaned_preds"])
                else:
                    logger.warning("No cleaned_preds column. Skipping.")
                    continue

                df = add_quantification_data(df, clean_run_name, inputs_folder=INPUTS_FOLDER)

                clean_list = df["cleaned_preds"].tolist()
                filtered = preprocessing.filter_contaminants(clean_list, clean_run_name, CONTAMINANTS_FASTA)
                df = df[df["cleaned_preds"].isin(filtered)]

                # Loop FDR
                for fdr in FDR_THRESHOLDS:
                    if "psm_q_value" in df.columns:
                        subset = df[df["psm_q_value"] <= fdr].copy()
                    else:
                        subset = df.copy()

                    input_seqs = subset["cleaned_preds"].tolist()

                    # Helper function to add empty/zero results
                    def add_result(
                        cov=0,
                        scaf_count=0,
                        # Capture loop variables here:
                        cat=category,
                        samp=sample_label,
                        run=clean_run_name,
                        ch=chain_type,
                        f=fdr,
                    ):
                        all_results.append(
                            {
                                "Category": cat,
                                "Sample": samp,
                                "Run": run,
                                "Chain": ch,
                                "FDR": f,
                                "Coverage": cov,
                                "Scaffolds": scaf_count,
                            }
                        )

                    if not input_seqs:
                        add_result()
                        continue

                    try:
                        # Assembly is always the same (Reference-free)
                        scaffolds = assembler.run(sequences=input_seqs, df_full=subset)
                    except Exception:
                        scaffolds = []

                    if not scaffolds:
                        add_result()
                        continue

                    # Mapping: Here we match scaffolds against the SPECIFIC current chain
                    mapped = visualization.process_protein_contigs_scaffold(
                        scaffolds, protein_norm, MAX_MISMATCHES, MIN_IDENTITY
                    )

                    cov = 0
                    if mapped:
                        df_map = visualization.create_dataframe_from_mapped_sequences(mapped)
                        stats = helpers.compute_assembly_statistics(
                            df=df_map,
                            sequence_type="temp",
                            output_folder=str(OUTPUT_FOLDER),
                            reference=protein_norm,
                        )
                        cov = stats["coverage"] * 100

                    add_result(cov=cov, scaf_count=len(scaffolds))

    if not all_results:
        logger.error("No results generated.")
        return

    results_df = pd.DataFrame(all_results)
    results_df.to_csv(OUTPUT_FOLDER / "aggregated_results.csv", index=False)

    # --- PLOTTING ---
    custom_palette = load_custom_palette()
    mode_output = OUTPUT_FOLDER / ASSEMBLY_MODE
    os.makedirs(mode_output, exist_ok=True)

    sns.set_style("whitegrid")
    sns.set_context("paper", font_scale=1.2)

    # Plotting: use 'style' to differentiate chains if present
    g = sns.relplot(
        data=results_df,
        x="FDR",
        y="Coverage",
        col="Category",
        col_wrap=2,
        kind="line",
        hue="Category",
        style="Chain",  # Differentiate Heavy/Light with different line styles
        markers=True,
        dashes=True,
        linewidth=2.5,
        palette=custom_palette,
        markersize=8,
        height=4,
        aspect=1.5,
        legend="full",  # Enable automatic legend to see chains
    )

    g.fig.suptitle(f"{ASSEMBLY_MODE} aggregated assembly performance (Mean ± 95% CI)", fontsize=16, y=1.02)
    g.fig.subplots_adjust(top=0.85, wspace=0.3, hspace=0.4)

    g.set_axis_labels("FDR threshold", "Sequence coverage (%)")
    g.set_titles("{col_name}")
    g.set(ylim=(0, 105))
    g.set(xticks=FDR_THRESHOLDS)

    for ax in g.axes.flat:
        ax.set_xticklabels([f"{int(x * 100)}%" for x in FDR_THRESHOLDS])
        ax.grid(True, which="major", color="#dddddd", linewidth=0.8)

    sns.move_legend(
        g,
        "lower center",
        bbox_to_anchor=(0.5, 1.02),
        ncol=4,
        title=None,
        frameon=False,
    )

    plt.savefig(mode_output / "aggregated_coverage_faceted.svg", bbox_inches="tight")
    plt.savefig(mode_output / "aggregated_coverage_faceted.png", dpi=300, bbox_inches="tight")

    logger.info(f"Aggregated plots saved to: {mode_output}")


if __name__ == "__main__":
    main()
