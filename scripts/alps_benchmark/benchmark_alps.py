#!/usr/bin/env python
import argparse
import json
import logging
import math
import os
import random
import re
import subprocess
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import pandas as pd
import seaborn as sns
from Bio import SeqIO

from instanexus import helpers, preprocessing, visualization
from instanexus.assembly import Assembler

# --- CONSTANTS ---
MAX_REFINE_ROUNDS = 10
COVERAGE_THRESHOLD = 0.9
ACCURACY_THRESHOLD = 0.7
PRECISION_THRESHOLD = 0.4
MAX_LEN_IDENTITY_THRESHOLD = 0.95

SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = SCRIPT_DIR.parent
if str(PROJECT_ROOT) not in sys.path:
    sys.path.append(str(PROJECT_ROOT))

ALPS_JAR_PATH = SCRIPT_DIR / "ALPS.jar"
if not ALPS_JAR_PATH.exists():
    ALPS_JAR_PATH = PROJECT_ROOT / "benchmark" / "ALPS.jar"

METADATA_PATH = PROJECT_ROOT / "json" / "sample_metadata.json"
CONTAMINANTS_FASTA = PROJECT_ROOT / "fasta/contaminants.fasta"
INPUTS_FOLDER = PROJECT_ROOT / "inputs"
BASE_OUTPUT_FOLDER = PROJECT_ROOT / "outputs/alps_benchmark_final"

SAMPLE_GROUPS = {
    "BSA": ["bsa.csv"],
    "Nanobodies": [f"nb{i}.csv" for i in range(1, 11)],
    "Antibodies": ["ma1.csv", "ma2.csv", "ma3.csv"],
    "Binders": ["bind1.csv", "bind2.csv", "bind3.csv"],
}

FDR_THRESHOLDS = [0.01, 0.05, 0.10]
ALPS_MAX_SCAFFOLDS = "5000"

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
logger = logging.getLogger(__name__)


def calculate_nx(seqs, x=50):
    if not seqs:
        return 0
    lengths = sorted([len(s) for s in seqs], reverse=True)
    target = sum(lengths) * (x / 100.0)
    cum_sum = 0
    for length in lengths:
        cum_sum += length
        if cum_sum >= target:
            return length
    return 0


# --- ROBUST PREPROCESSING (Local Override) ---
def robust_clean_dataframe(df):
    """Handles missing columns (V2 format) gracefully."""
    df = df.copy()
    seq_candidates = ["cleaned_preds", "prediction", "Peptide", "sequence", "prediction_untokenised"]
    found_col = None
    for col in seq_candidates:
        if col in df.columns:
            found_col = col
            break
    if found_col:
        df.rename(columns={found_col: "cleaned_preds"}, inplace=True)
    else:
        return pd.DataFrame()

    if "instanovo_token_log_probabilities" not in df.columns:
        df["instanovo_token_log_probabilities"] = None

    if "cleaned_preds" in df.columns:

        def clean_seq(seq):
            if not isinstance(seq, str):
                return ""
            seq = re.sub(r"\(.*?\)", "", seq)
            seq = re.sub(r"\[.*?\]", "", seq)
            return seq.strip()

        df["cleaned_preds"] = df["cleaned_preds"].apply(clean_seq)
        df = df[df["cleaned_preds"].str.len() > 0]
    return df


def get_advanced_stats(fasta_path, tool_name, protein_norm, fdr, run_label, category, input_psm_count):
    base_res = {
        "Category": category,
        "Tool": tool_name,
        "Sample": run_label,
        "FDR": fdr,
        "Coverage": 0,
        "N50": 0,
        "N90": 0,
        "Mean_Identity": 0,
        "Spectral_Depth": 0,
        "Num_Scaffolds": 0,
        "Max_Length": 0,
        "Precision": 0,
        "Assembly_Score": 0,
    }

    if not fasta_path or not Path(fasta_path).exists():
        return base_res
    seqs = [str(r.seq) for r in SeqIO.parse(fasta_path, "fasta")]
    if not seqs:
        return base_res

    total_scaffolds = len(seqs)
    n50 = calculate_nx(seqs, 50)
    n90 = calculate_nx(seqs, 90)

    mapped = visualization.process_protein_contigs_scaffold(
        assembled_contigs=seqs, target_protein=protein_norm, max_mismatches=1000, min_identity=PRECISION_THRESHOLD
    )

    if not mapped:
        return {**base_res, "Num_Scaffolds": total_scaffolds, "N50": n50, "N90": n90}

    df_map = visualization.create_dataframe_from_mapped_sequences(mapped)
    mapped_count = len(df_map)
    precision = (mapped_count / total_scaffolds) * 100

    df_hq_len = df_map[df_map["identity_score"] >= MAX_LEN_IDENTITY_THRESHOLD]
    max_hq_len = df_hq_len["end"].sub(df_hq_len["start"]).max() if not df_hq_len.empty else 0

    df_acc = df_map[df_map["identity_score"] >= ACCURACY_THRESHOLD]
    mean_identity = df_acc["identity_score"].mean() * 100 if not df_acc.empty else 0

    df_cov = df_map[df_map["identity_score"] >= COVERAGE_THRESHOLD].copy()
    coverage, spectral_depth = 0, 0
    if not df_cov.empty:
        parent_dir = str(Path(fasta_path).parent)
        s = helpers.compute_assembly_statistics(
            df=df_cov, sequence_type="temp_stats", output_folder=parent_dir, reference=protein_norm
        )
        coverage = s["coverage"] * 100
        covered_span = df_cov["end"].max() - df_cov["start"].min()
        if covered_span > 0:
            spectral_depth = input_psm_count / covered_span

    norm_coverage = coverage / 100.0
    norm_identity = mean_identity / 100.0
    norm_precision = precision / 100.0
    contiguity_factor = 1.0 / math.log10(total_scaffolds + 9) if total_scaffolds > 0 else 1.0
    assembly_score = (norm_coverage * norm_identity * norm_precision * contiguity_factor) * 100

    return {
        "Category": category,
        "Tool": tool_name,
        "Sample": run_label,
        "FDR": fdr,
        "Coverage": coverage,
        "N50": n50,
        "N90": n90,
        "Mean_Identity": mean_identity,
        "Spectral_Depth": spectral_depth,
        "Num_Scaffolds": total_scaffolds,
        "Max_Length": max_hq_len,
        "Precision": precision,
        "Assembly_Score": assembly_score,
    }


def prepare_alps_input(df_subset, output_csv_path):
    df = df_subset.copy()
    col_map = {"cleaned_preds": "Peptide", "prediction": "Peptide"}
    df.rename(columns=col_map, inplace=True)
    if "Peptide" not in df.columns:
        return False
    df["Peptide"] = df["Peptide"].apply(lambda x: re.sub(r"\(.*?\)", "", str(x)))
    df["Peptide"] = df["Peptide"].apply(lambda x: re.sub(r"\[.*?\]", "", str(x)))
    df["local confidence (%)"] = [
        " ".join([str(random.randint(85, 99)) for _ in range(len(str(p)))]) for p in df["Peptide"]
    ]
    df["Scan"] = [f"F1:{10000 + i}" for i in range(len(df))]
    df["Area"] = 100000
    df[["Scan", "Peptide", "local confidence (%)", "Area"]].to_csv(output_csv_path, index=False)
    return True


def plot_bar_metric(df, metric, y_label, title, output_path, subtitle=""):
    colors = {"InstaNexus": "#6BAED6", "ALPS": "#E74C3C"}
    sns.set_style("whitegrid")
    sns.set_context("paper", font_scale=1.4)

    g = sns.catplot(
        data=df,
        x="FDR",
        y=metric,
        hue="Tool",
        col="Category",
        kind="bar",
        palette=colors,
        height=5,
        aspect=0.9,
        legend=False,
        errorbar=("ci", 95),
        capsize=0.1,
        err_kws={"linewidth": 1.5, "color": "black"},
    )

    full_title = f"{title}\n" + (f"({subtitle})" if subtitle else "")
    g.fig.suptitle(full_title, fontsize=16, y=1.02, fontweight="bold")
    g.set_axis_labels("FDR Threshold", y_label)
    g.set_titles("{col_name}")
    for ax in g.axes.flat:
        ax.set_xticklabels([f"{int(x * 100)}%" for x in FDR_THRESHOLDS])

    handles = [
        mpatches.Patch(color=colors["InstaNexus"], label="InstaNexus"),
        mpatches.Patch(color=colors["ALPS"], label="ALPS"),
    ]
    g.fig.legend(
        handles=handles, loc="upper center", bbox_to_anchor=(0.5, 0.92), ncol=2, frameon=False, fontsize=14, title=""
    )
    g.fig.subplots_adjust(top=0.82, wspace=0.3, bottom=0.10)
    plt.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--mode", type=str, default="dbg_weighted")
    parser.add_argument("--kmer", type=int, default=7)
    parser.add_argument("--refine", action="store_true")
    args = parser.parse_args()

    assembly_mode = args.mode
    kmer_val = args.kmer
    refine_rounds = MAX_REFINE_ROUNDS if args.refine else 0

    print(f"--- Configuration: {assembly_mode}, k={kmer_val}, Refine={args.refine} ---")

    if not ALPS_JAR_PATH.exists():
        logger.error(f"ALPS.jar MISSING at {ALPS_JAR_PATH}")
        return

    with open(METADATA_PATH, "r") as f:
        FULL_METADATA = json.load(f)

    plot_title_suffix = f"{assembly_mode}, k={kmer_val}"
    if args.refine:
        plot_title_suffix += " [Refined]"

    out_name = f"{assembly_mode}_k{kmer_val}" + ("_refined" if args.refine else "")
    main_output_dir = BASE_OUTPUT_FOLDER / out_name
    dir_csv, dir_fasta, dir_svg = main_output_dir / "csv", main_output_dir / "fasta", main_output_dir / "svg"
    for d in [dir_csv, dir_fasta, dir_svg]:
        os.makedirs(d, exist_ok=True)

    assembler = Assembler(
        mode=assembly_mode, kmer_size=kmer_val, min_overlap=3, size_threshold=0, refine_rounds=refine_rounds
    )
    benchmark_results = []

    for category, file_list in SAMPLE_GROUPS.items():
        logger.info(f"=== Category: {category} ===")
        for filename in file_list:
            raw_path = INPUTS_FOLDER / filename
            file_stem = raw_path.stem
            if not raw_path.exists():
                continue

            meta_entries = FULL_METADATA.get(file_stem, [])
            if isinstance(meta_entries, dict):
                meta_entries = [meta_entries]
            if not meta_entries:
                continue

            df = pd.read_csv(raw_path)
            if "experiment_name" in df.columns:
                proteases = meta_entries[0].get("proteases", [])
                df["protease"] = df["experiment_name"].apply(
                    lambda x, p=proteases: preprocessing.extract_protease(x, p)
                )

            try:
                df = robust_clean_dataframe(df)
            except Exception:
                continue
            if "cleaned_preds" not in df.columns or df.empty:
                continue

            clean_list = df["cleaned_preds"].tolist()
            filtered = preprocessing.filter_contaminants(clean_list, file_stem, CONTAMINANTS_FASTA)
            df = df[df["cleaned_preds"].isin(filtered)]

            for fdr in FDR_THRESHOLDS:
                subset = df[df["psm_q_value"] <= fdr].copy() if "psm_q_value" in df.columns else df.copy()
                input_seqs = subset["cleaned_preds"].tolist()
                if not input_seqs:
                    continue

                # InstaNexus
                insta_fasta = dir_fasta / f"{file_stem}_FDR{fdr}_insta.fasta"
                try:
                    scaffolds = assembler.run(sequences=input_seqs, df_full=subset)
                    with open(insta_fasta, "w") as f:
                        for i, seq in enumerate(scaffolds):
                            f.write(f">scaffold_{i}\n{seq}\n")
                except Exception:
                    open(insta_fasta, "w").close()

                # ALPS
                alps_input = dir_csv / f"{file_stem}_FDR{fdr}_alps_input.csv"
                prepare_alps_input(subset, alps_input)
                # Cleanup old files
                for f in list(dir_csv.glob(f"{file_stem}_FDR{fdr}*alps*.fasta")):
                    os.remove(f)

                cmd = ["java", "-jar", str(ALPS_JAR_PATH), str(alps_input), str(kmer_val), ALPS_MAX_SCAFFOLDS]
                alps_raw_fasta = None
                try:
                    subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
                    found = list(dir_csv.glob(f"*{alps_input.stem}*.fasta"))
                    if found:
                        alps_raw_fasta = max(found, key=os.path.getctime)
                except Exception:
                    pass

                for meta in meta_entries:
                    chain, p_seq = meta.get("chain", ""), meta.get("protein", "")
                    p_norm = preprocessing.normalize_sequence(p_seq)
                    lbl = f"{file_stem} ({chain})" if chain else file_stem
                    benchmark_results.append(
                        get_advanced_stats(insta_fasta, "InstaNexus", p_norm, fdr, lbl, category, len(subset))
                    )
                    benchmark_results.append(
                        get_advanced_stats(alps_raw_fasta, "ALPS", p_norm, fdr, lbl, category, len(subset))
                    )

    if benchmark_results:
        df_res = pd.DataFrame(benchmark_results)
        df_res.to_csv(dir_csv / "benchmark_summary_raw.csv", index=False)
        metrics = {
            "Assembly_Score": ("Score", "Assembly Quality Score (AQS)"),
            "Coverage": ("Coverage (%)", "High-Confidence Coverage"),
            "N50": ("N50 (AA)", "Contiguity (N50)"),
            "N90": ("N90 (AA)", "Robustness (N90)"),
            "Precision": ("Mapped Scaffolds (%)", "Precision (Efficiency)"),
            "Mean_Identity": ("Mean Identity (%)", f"Accuracy (Mapped >{int(ACCURACY_THRESHOLD * 100)}%)"),
            "Num_Scaffolds": ("Count (Raw)", "Fragmentation"),
            "Max_Length": ("Length (AA)", f"Longest Sequence (>{int(MAX_LEN_IDENTITY_THRESHOLD * 100)}% Id)"),
        }
        for m, (ylab, tit) in metrics.items():
            plot_bar_metric(df_res, m, ylab, tit, dir_svg / f"{m.lower()}_barplot.svg", subtitle=plot_title_suffix)
        logger.info(f"Benchmark finished. Data saved in {dir_csv}")


if __name__ == "__main__":
    main()
