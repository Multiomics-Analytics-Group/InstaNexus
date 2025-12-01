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
__date__ = 29 Nov 2025
__maintainer__ = Marco Reverenna
__email__ = marcor@dtu.dk
__status__ = Dev
"""

import matplotlib.pyplot as plt
import matplotlib.lines as mlines
import pandas as pd
import seaborn as sns
import os
import sys
import logging
import json
from pathlib import Path

from instanexus.assembly import Assembler
from instanexus import visualization, helpers, preprocessing

SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = SCRIPT_DIR.parent
if str(PROJECT_ROOT) not in sys.path:
    sys.path.append(str(PROJECT_ROOT))

INPUTS_FOLDER = PROJECT_ROOT / "inputs"
METADATA_JSON = PROJECT_ROOT / "json/sample_metadata.json"
COLORS_JSON = PROJECT_ROOT / "json" / "colors.json"
CONTAMINANTS_FASTA = PROJECT_ROOT / "fasta/contaminants.fasta"
OUTPUT_FOLDER = PROJECT_ROOT / "outputs" / "_fdr_analysis"

SAMPLE_GROUPS = {
    "BSA": ["bsa.csv"],
    "Nanobodies": [
        "nb1.csv",
        "nb2.csv",
        "nb3.csv",
        "nb4.csv",
        "nb5.csv",
        "nb6.csv",
        "nb7.csv",
        "nb8.csv",
        "nb9.csv",
        "nb10.csv",
    ],
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
ASSEMBLY_MODE = "multimodal_dbg"
KMER_SIZE = 7
MIN_OVERLAP = 4
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

        df_sum = df_quant.groupby("cleaned_preds", as_index=False)[
            "total_abundance_norm"
        ].sum()
        df_sum.rename(
            columns={"total_abundance_norm": "peptide_abundance"}, inplace=True
        )

        df_merged = pd.merge(df_main, df_sum, on="cleaned_preds", how="left")
        df_merged["peptide_abundance"] = df_merged["peptide_abundance"].fillna(0)
        return df_merged
    except:
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
    except:
        return default_palette


def main():
    os.makedirs(OUTPUT_FOLDER, exist_ok=True)

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
        logger.info(
            f"=== Processing Category: {category} ({len(file_list)} samples) ==="
        )

        for filename in file_list:
            csv_path = INPUTS_FOLDER / filename
            run_name = csv_path.stem

            logger.info(f"   -> Sample: {run_name}")

            if not csv_path.exists():
                logger.warning(f"File not found: {csv_path}. Skipping.")
                continue

            try:
                clean_run_name = run_name.replace("_cleaned", "")
                meta = preprocessing.get_sample_metadata(
                    clean_run_name, json_path=METADATA_JSON
                )
                protein_norm = preprocessing.normalize_sequence(meta.get("protein", ""))
                proteases = meta.get("proteases", [])
            except Exception as e:
                logger.warning(f"Metadata error for {run_name}: {e}. Skipping.")
                continue

            df = pd.read_csv(csv_path)

            if "experiment_name" in df.columns:
                df["protease"] = df["experiment_name"].apply(
                    lambda x: preprocessing.extract_protease(x, proteases)
                )

            df = preprocessing.clean_dataframe(df)

            if "cleaned_preds" in df.columns:
                df["cleaned_preds"] = df["cleaned_preds"].apply(
                    preprocessing.remove_modifications
                )
                df = df.dropna(subset=["cleaned_preds"])
            else:
                continue

            df = add_quantification_data(
                df, clean_run_name, inputs_folder=INPUTS_FOLDER
            )

            clean_list = df["cleaned_preds"].tolist()
            filtered = preprocessing.filter_contaminants(
                clean_list, clean_run_name, CONTAMINANTS_FASTA
            )
            df = df[df["cleaned_preds"].isin(filtered)]

            for fdr in FDR_THRESHOLDS:
                if "psm_q_value" in df.columns:
                    subset = df[df["psm_q_value"] <= fdr].copy()
                else:
                    subset = df.copy()

                input_seqs = subset["cleaned_preds"].tolist()

                if not input_seqs:
                    all_results.append(
                        {
                            "Category": category,
                            "Sample": run_name,
                            "FDR": fdr,
                            "Coverage": 0,
                            "Scaffolds": 0,
                        }
                    )
                    continue

                try:
                    scaffolds = assembler.run(sequences=input_seqs, df_full=subset)
                except Exception:
                    scaffolds = []

                if not scaffolds:
                    all_results.append(
                        {
                            "Category": category,
                            "Sample": run_name,
                            "FDR": fdr,
                            "Coverage": 0,
                            "Scaffolds": 0,
                        }
                    )
                    continue

                mapped = visualization.process_protein_contigs_scaffold(
                    scaffolds, protein_norm, MAX_MISMATCHES, MIN_IDENTITY
                )

                cov = 0
                if mapped:
                    df_map = visualization.create_dataframe_from_mapped_sequences(
                        mapped
                    )
                    stats = helpers.compute_assembly_statistics(
                        df=df_map,
                        sequence_type="temp",
                        output_folder=str(OUTPUT_FOLDER),
                        reference=protein_norm,
                    )
                    cov = stats["coverage"] * 100

                all_results.append(
                    {
                        "Category": category,
                        "Sample": run_name,
                        "FDR": fdr,
                        "Coverage": cov,
                        "Scaffolds": len(scaffolds),
                    }
                )

    if not all_results:
        logger.error("No results generated.")
        return

    results_df = pd.DataFrame(all_results)
    results_df.to_csv(OUTPUT_FOLDER / "aggregated_results.csv", index=False)

    custom_palette = load_custom_palette()

    mode_output = OUTPUT_FOLDER / ASSEMBLY_MODE
    os.makedirs(mode_output, exist_ok=True)

    sns.set_style("whitegrid")

    g = sns.relplot(
        data=results_df,
        x="FDR",
        y="Coverage",
        col="Category",
        col_wrap=2,
        kind="line",
        hue="Category",
        style="Category",
        markers=True,
        dashes=False,
        linewidth=3,
        palette=custom_palette,
        markersize=9,
        height=4,
        aspect=1.5,
        legend=False,
    )

    g.fig.subplots_adjust(top=0.82, wspace=0.3, hspace=0.4)
    g.fig.suptitle(
        "Aggregated assembly performance (Mean ± 95% CI)", fontsize=16, y=0.98
    )

    legend_handles = []
    for cat in SAMPLE_GROUPS.keys():
        if cat in custom_palette:
            h = mlines.Line2D(
                [],
                [],
                color=custom_palette[cat],
                marker="o",
                linewidth=3,
                markersize=9,
                label=cat,
            )
            legend_handles.append(h)

    g.fig.legend(
        handles=legend_handles,
        loc="upper center",
        bbox_to_anchor=(0.5, 0.92),
        ncol=4,
        frameon=False,
        fontsize=12,
    )
    g.set_axis_labels("FDR threshold", "Sequence coverage (%)")
    g.set_titles("{col_name}")
    g.set(ylim=(0, 105))
    g.set(xticks=FDR_THRESHOLDS)

    for ax in g.axes.flat:
        ax.set_xticklabels([f"{int(x * 100)}%" for x in FDR_THRESHOLDS])

    g.fig.subplots_adjust(top=0.82, wspace=0.3, hspace=0.4)

    plt.savefig(mode_output / "aggregated_coverage_faceted.svg", bbox_inches="tight")
    plt.savefig(
        mode_output / "aggregated_coverage_faceted.png", dpi=300, bbox_inches="tight"
    )

    logger.info(f"Aggregated plots saved to: {mode_output}")


if __name__ == "__main__":
    main()
