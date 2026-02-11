#!/usr/bin/env python

r"""FDR analysis script: Aggregated Categories (CLI enabled)."""

import argparse
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
PROJECT_ROOT = SCRIPT_DIR.parent.parent

if str(PROJECT_ROOT) not in sys.path:
    sys.path.append(str(PROJECT_ROOT))

INPUTS_FOLDER = PROJECT_ROOT / "inputs"
METADATA_JSON = PROJECT_ROOT / "json/sample_metadata.json"
COLORS_JSON = PROJECT_ROOT / "json" / "colors.json"
CONTAMINANTS_FASTA = PROJECT_ROOT / "fasta/contaminants.fasta"
BASE_OUTPUT_FOLDER = PROJECT_ROOT / "outputs" / "_fdr_analysis"

SAMPLE_GROUPS = {
    "BSA": ["bsa.csv"],
    "Nanobodies": [f"nb{i}.csv" for i in range(1, 11)],
    "Antibodies": ["ma1.csv", "ma2.csv", "ma3.csv"],
    "Binders": ["bind1.csv", "bind2.csv", "bind3.csv"],
}

CATEGORY_COLOR_MAP = {
    "BSA": "bsa",
    "Nanobodies": "nanobodies",
    "Antibodies": "antibodies",
    "Binders": "binders",
}

FDR_THRESHOLDS = [0.01, 0.05, 0.10, 0.20, 0.40, 1]
# FDR_THRESHOLDS = [0.01, 0.05]

DEFAULT_KMER = 7
MIN_OVERLAP = 3
SIZE_THRESHOLD = 10
MIN_IDENTITY = 0.8
MAX_MISMATCHES = 100
MAX_REFINE_ROUNDS = 10

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(message)s")
logger = logging.getLogger(__name__)


def add_quantification_data(df_main, run_name, inputs_folder):
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
    parser = argparse.ArgumentParser(description="Run FDR Analysis with specific Assembler.")
    parser.add_argument(
        "--mode", type=str, default="dbg_weighted", help="Assembly mode (e.g., greedy, dbg_weighted, multimodal_dbg)"
    )
    parser.add_argument("--refine", action="store_true", help="Enable iterative refinement (Overlap Graph)")
    parser.add_argument("--kmer", type=int, default=DEFAULT_KMER, help="K-mer size")

    args = parser.parse_args()

    assembly_mode = args.mode
    kmer_size = args.kmer
    refine_rounds = MAX_REFINE_ROUNDS if args.refine else 0

    folder_suffix = assembly_mode
    if args.refine:
        folder_suffix += "_refined"

    current_output_folder = BASE_OUTPUT_FOLDER / folder_suffix
    os.makedirs(current_output_folder, exist_ok=True)

    logger.info("--- Configuration ---")
    logger.info(f"Mode:   {assembly_mode}")
    logger.info(f"Refine: {'Enabled' if args.refine else 'Disabled'}")
    logger.info(f"Output: {current_output_folder}")
    logger.info("---------------------")

    with open(METADATA_JSON, "r") as f:
        FULL_METADATA = json.load(f)

    assembler = Assembler(
        mode=assembly_mode,
        kmer_size=kmer_size,
        min_overlap=MIN_OVERLAP,
        size_threshold=SIZE_THRESHOLD,
        min_identity=MIN_IDENTITY,
        max_mismatches=MAX_MISMATCHES,
        refine_rounds=refine_rounds,
    )

    all_results = []

    for category, file_list in SAMPLE_GROUPS.items():
        logger.info(f"=== Processing Category: {category} ===")

        for filename in file_list:
            csv_path = INPUTS_FOLDER / filename
            run_name = csv_path.stem
            clean_run_name = run_name.replace("_cleaned", "")

            if not csv_path.exists():
                continue

            meta_entries = FULL_METADATA.get(clean_run_name, [])
            if isinstance(meta_entries, dict):
                meta_entries = [meta_entries]
            if not meta_entries:
                continue

            df_original = pd.read_csv(csv_path)

            for meta in meta_entries:
                target_protein = meta.get("protein", "")
                chain_type = meta.get("chain", "")
                proteases = meta.get("proteases", [])
                sample_label = f"{clean_run_name} ({chain_type})" if chain_type else clean_run_name
                protein_norm = preprocessing.normalize_sequence(target_protein)

                df = df_original.copy()

                if "experiment_name" in df.columns:
                    df["protease"] = df["experiment_name"].apply(
                        lambda x, p=proteases: preprocessing.extract_protease(x, p)
                    )

                try:
                    df = preprocessing.clean_dataframe(df)
                except Exception:
                    continue

                if "cleaned_preds" in df.columns:
                    df["cleaned_preds"] = df["cleaned_preds"].apply(preprocessing.remove_modifications)
                    df = df.dropna(subset=["cleaned_preds"])
                else:
                    continue

                df = add_quantification_data(df, clean_run_name, inputs_folder=INPUTS_FOLDER)

                clean_list = df["cleaned_preds"].tolist()
                filtered = preprocessing.filter_contaminants(clean_list, clean_run_name, CONTAMINANTS_FASTA)
                df = df[df["cleaned_preds"].isin(filtered)]

                for fdr in FDR_THRESHOLDS:
                    if "psm_q_value" in df.columns:
                        subset = df[df["psm_q_value"] <= fdr].copy()
                    else:
                        subset = df.copy()

                    input_seqs = subset["cleaned_preds"].tolist()

                    def add_result(
                        cov=0,
                        scaf_count=0,
                        _category=category,
                        _sample_label=sample_label,
                        _clean_run_name=clean_run_name,
                        _chain_type=chain_type,
                        _fdr=fdr,
                    ):
                        all_results.append(
                            {
                                "Category": _category,
                                "Sample": _sample_label,
                                "Run": _clean_run_name,
                                "Chain": _chain_type,
                                "FDR": _fdr,
                                "Coverage": cov,
                                "Scaffolds": scaf_count,
                            }
                        )

                    if not input_seqs:
                        add_result()
                        continue

                    try:
                        scaffolds = assembler.run(sequences=input_seqs, df_full=subset)
                    except Exception:
                        scaffolds = []

                    if not scaffolds:
                        add_result()
                        continue

                    mapped = visualization.process_protein_contigs_scaffold(
                        scaffolds, protein_norm, MAX_MISMATCHES, MIN_IDENTITY
                    )

                    cov = 0
                    if mapped:
                        df_map = visualization.create_dataframe_from_mapped_sequences(mapped)
                        stats = helpers.compute_assembly_statistics(
                            df=df_map,
                            sequence_type="temp",
                            output_folder=str(current_output_folder),
                            reference=protein_norm,
                        )
                        cov = stats["coverage"] * 100

                    add_result(cov=cov, scaf_count=len(scaffolds))

    if not all_results:
        logger.error("No results generated.")
        return

    results_df = pd.DataFrame(all_results)
    results_df.to_csv(current_output_folder / "aggregated_results.csv", index=False)

    results_df["Chain"] = results_df["Chain"].replace("", "N/A").fillna("N/A")

    results_df = results_df.sort_values("FDR")
    results_df["FDR_Label"] = results_df["FDR"].apply(lambda x: f"{int(x * 100)}%")

    marker_map = {"heavy": "s", "light": "D", "N/A": "o"}

    visualization.set_publication_style()

    fig_w, fig_h = visualization.get_figsize(width_ratio=3)
    panel_height = 3.5
    calc_aspect = (fig_w / 2) / panel_height

    custom_palette = load_custom_palette()

    g = sns.relplot(
        data=results_df,
        x="FDR_Label",
        y="Coverage",
        col="Category",
        col_wrap=2,
        kind="line",
        hue="Category",
        style="Chain",
        markers=marker_map,
        dashes=True,
        linewidth=2.5,
        palette=custom_palette,
        markersize=9,
        height=panel_height,
        aspect=calc_aspect,
        legend="full",
    )

    g.fig.subplots_adjust(wspace=0.3, hspace=0.4)

    g.set_axis_labels("FDR threshold", "Sequence coverage (%)")
    g.set_titles("{col_name}")
    g.set(ylim=(0, 105))

    for ax in g.axes.flat:
        ax.grid(True, which="major", color="#dddddd", linewidth=0.8, alpha=0.5)

    sns.move_legend(
        g,
        "lower center",
        bbox_to_anchor=(0.5, 1.02),
        ncol=5,
        title=None,
        frameon=False,
    )

    plt.savefig(current_output_folder / "aggregated_coverage_faceted.svg", bbox_inches="tight")
    plt.savefig(current_output_folder / "aggregated_coverage_faceted.pdf", bbox_inches="tight")

    logger.info(f"Aggregated plots saved to: {current_output_folder}")


if __name__ == "__main__":
    main()
