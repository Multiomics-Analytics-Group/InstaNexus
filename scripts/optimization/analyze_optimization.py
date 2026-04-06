#!/usr/bin/env python

"""
Optimization Analysis Script for InstaNexus.
Focus: Grid Search Analysis using Seaborn/Matplotlib.
Output: Publication-ready tables and SVG Heatmaps.
"""

import glob
import json
import os
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from matplotlib.colors import LinearSegmentedColormap

# --- CONFIGURATION ---

INPUT_DIR = Path("outputs/_grid_search")
OUTPUT_DIR = Path("outputs/_optimization_figures")
SUMMARY_DIR = Path("outputs/_summary_tables")

ABBREV = {
    "kmer_size": "k",
    "min_overlap": "mo",
    "size_threshold": "st",
    "min_weight": "mw",
    "fdr": "fdr",
    "refine_rounds": "rr",
}

# Definiamo i colori. Per Matplotlib, convertiremo queste liste in Colormap oggetti.
THEME_MAP_DATA = {
    "BSA": ["#fee0d2", "#fc9272", "#de2d26", "#67000d"],
    "Antibodies": ["#e0f3db", "#a8ddb5", "#41ab5d", "#00441b"],
    "Nanobodies": ["#eff3ff", "#bdd7e7", "#6baed6", "#2171b5"],
    "Binders": ["#F4F5F0", "#D8dbca", "#A3A694", "#46473E"],
    "Other": ["#f7f7f7", "#cccccc", "#969696", "#525252"],
}

HEATMAP_CONFIG = {
    "greedy": (["min_overlap"], ["size_threshold", "fdr"]),
    "dbg_weighted": (["kmer_size", "min_overlap"], ["size_threshold", "fdr"]),
    "multimodal_dbg": (["kmer_size", "min_weight"], ["size_threshold", "fdr"]),
    "default": (["min_overlap"], ["fdr"]),
}


def set_publication_style():
    """Sets the visual style for publication-quality figures."""
    sns.set_theme(style="ticks")

    plt.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.sans-serif": ["Arial"],
            "font.size": 14,
            "axes.titlesize": 16,
            "axes.labelsize": 15,
            "xtick.labelsize": 12,
            "ytick.labelsize": 12,
            "axes.linewidth": 1.5,
            "figure.dpi": 300,
            "legend.fontsize": 13,
            "legend.frameon": False,
            "legend.columnspacing": 1.5,
            "xtick.major.width": 1,
            "ytick.major.width": 1,
        }
    )


def get_figsize(width_ratio=1, total_width_inch=14.0):
    """Calculates figsize based on standard widths."""
    ratios = {
        1: (1, 1),
        2: (2, 1),
        3: (3, 1),
    }
    w_mult, h_mult = ratios.get(width_ratio, (1, 1))
    actual_width = (total_width_inch / 3) * w_mult
    actual_height = (actual_width / w_mult) * 0.85
    return (actual_width, actual_height)


def create_colormap(color_list, name="custom_cmap"):
    """Creates a LinearSegmentedColormap from a list of hex codes."""
    return LinearSegmentedColormap.from_list(name, color_list)


def plot_grid_search_clustermap(df, index_cols, column_cols, theme_colors, value_col, title="", output_file=None):
    """
    Generates heatmaps using Seaborn with White Grid.
    Replicates the exact look of the previous Plotly version.
    """
    try:
        pivot = df.pivot_table(values=value_col, index=index_cols, columns=column_cols, aggfunc="max")
    except KeyError:
        return

    pivot = pivot.sort_index(level=index_cols).sort_index(axis=1, level=column_cols).fillna(0)

    row_labels = [
        ", ".join(
            f"{ABBREV.get(c, c)}={v}"
            for c, v in zip(index_cols, idx if isinstance(idx, tuple) else [idx], strict=False)
        )
        for idx in pivot.index
    ]
    col_labels = [
        ", ".join(
            f"{ABBREV.get(c, c)}={v}"
            for c, v in zip(column_cols, col if isinstance(col, tuple) else [col], strict=False)
        )
        for col in pivot.columns
    ]

    fig_size = get_figsize(width_ratio=1)
    fig, ax = plt.subplots(figsize=fig_size)

    cmap = create_colormap(theme_colors)

    sns.heatmap(
        pivot,
        ax=ax,
        cmap=cmap,
        annot=False,
        linewidths=2,
        linecolor="white",
        vmin=0.0,
        vmax=1.0,
        cbar_kws={"label": value_col, "shrink": 0.8},
    )

    ax.set_xticklabels(col_labels, rotation=45, ha="right")
    ax.set_yticklabels(row_labels, rotation=0)
    ax.set_xlabel("")
    ax.set_ylabel("")
    ax.set_title(title, pad=20)

    if output_file:
        output_file.parent.mkdir(parents=True, exist_ok=True)
        plt.savefig(str(output_file), bbox_inches="tight", format="svg", transparent=False)
        print(f"Saved: {output_file}")

    plt.close(fig)


def get_category(run_name):
    name_lower = str(run_name).lower()
    if "bsa" in name_lower:
        return "BSA"
    if "nb" in name_lower:
        return "Nanobodies"
    if "bind" in name_lower:
        return "Binders"
    if "ma" in name_lower or "pama" in name_lower:
        return "Antibodies"
    return "Other"


def format_assembly_method(row):
    method = row.get("mode", "unknown")
    seq_type = row.get("seq_type", "Contigs")
    return f"{method} ({seq_type})"


def normalize_composite_score(df):
    if "composite_score" not in df.columns:
        return df
    if df["composite_score"].max() > 1.5:
        df["composite_score"] = df["composite_score"] / 1000.0
    return df


def save_detailed_rankings(df_best, output_dir, mode_name):
    for cat, sub in df_best.groupby("category"):
        sub = sub.copy()
        if "scaffolds_count" in sub.columns:
            sub = sub.rename(columns={"scaffolds_count": "total_sequences"})

        sub["assembly_method"] = sub.apply(format_assembly_method, axis=1)
        sub = sub.rename(columns={"display_name": "sample"})

        desired_order = [
            "category",
            "sample",
            "assembly_method",
            "total_sequences",
            "N50",
            "coverage",
            "mean_identity",
            "composite_score",
        ]
        final_cols = [c for c in desired_order if c in sub.columns]

        final_df = sub[final_cols].sort_values(by=["sample", "assembly_method"], ascending=[True, True])
        numeric_cols = final_df.select_dtypes(include=["float64", "float32"]).columns
        final_df[numeric_cols] = final_df[numeric_cols].round(3)

        out_path = output_dir / mode_name / f"best_results_{cat}_{mode_name}.csv"
        out_path.parent.mkdir(parents=True, exist_ok=True)
        final_df.to_csv(out_path, index=False)


def combine_json_to_csv(
    run: str,
    type_method: str,
    type_sequence: str,
    base_path: Path = Path("outputs"),
) -> None:
    """Walks output directories, reads JSON stats files, and saves a combined CSV.

    Args:
        run: Run identifier (e.g. 'bsa', 'ma1').
        type_method: Assembly method prefix used in the JSON filename (e.g. 'scaffolds').
        type_sequence: Sequence type suffix used in the JSON filename (e.g. 'contigs').
        base_path: Root outputs folder.
    """
    run_path = Path(base_path) / run
    dataframes = []
    files_added = 0

    for root, dirs, _ in os.walk(run_path):
        for dir_name in dirs:
            json_path = Path(root) / dir_name / "statistics" / f"{type_method}_{type_sequence}_stats.json"
            if json_path.exists():
                try:
                    with open(json_path) as f:
                        data = json.load(f)
                    df = pd.json_normalize(data)
                    df["source"] = dir_name
                    dataframes.append(df)
                    files_added += 1
                except Exception as e:
                    print(f"Error loading {json_path}: {e}")

    if not dataframes:
        print(f"No stats files found under {run_path}.")
        return

    combined_df = pd.concat(dataframes, ignore_index=True)

    if "ass_method" in combined_df.columns:
        combined_df["ass_method"] = combined_df["ass_method"].fillna("greedy")

    combined_df["sequence_type"] = type_sequence
    combined_df["method_type"] = type_method
    combined_df["run"] = run

    output_file = run_path / f"{type_sequence}_combined_stats.csv"
    combined_df.to_csv(output_file, index=False)
    print(f"Combined stats saved to: {output_file}  ({files_added} files merged)")


def plot_coverages_from_runs(
    runs: list,
    base_path: Path = Path("outputs"),
    combination_folder: str = "",
    contigs_json: str = "contigs_stats.json",
    scaffolds_json: str = "scaffolds_stats.json",
    save: bool = False,
    output_dir: Path = Path("."),
) -> None:
    """Plots coverage barplots for contigs and scaffolds across multiple runs.

    Args:
        runs: List of run identifiers to include (e.g. ['bsa', 'nb1']).
        base_path: Root outputs folder.
        combination_folder: Sub-folder name for the specific parameter combination.
        contigs_json: Filename of the contigs stats JSON (default: contigs_stats.json).
        scaffolds_json: Filename of the scaffolds stats JSON (default: scaffolds_stats.json).
        save: If True, saves plots as PNG files.
        output_dir: Directory where PNG files are saved when save=True.
    """
    base_path = Path(base_path)
    contig_coverages: list = []
    scaffold_coverages: list = []
    labels: list = []

    for run in runs:
        stats_path = base_path / run / combination_folder / "statistics"
        if not stats_path.exists():
            print(f"[{run}] Missing statistics folder: {stats_path}")
            continue

        for coverage_list, fname in [(contig_coverages, contigs_json), (scaffold_coverages, scaffolds_json)]:
            json_path = stats_path / fname
            if json_path.exists():
                try:
                    with open(json_path) as f:
                        coverage_list.append(json.load(f).get("coverage", 0))
                except Exception as e:
                    print(f"[{run}] Error reading {fname}: {e}")
                    coverage_list.append(0)
            else:
                print(f"[{run}] {fname} not found.")
                coverage_list.append(0)

        labels.append(run)

    for coverages, color, title, suffix in [
        (contig_coverages, "mediumslateblue", "Contigs Coverage per Run", "contigs"),
        (scaffold_coverages, "seagreen", "Scaffolds Coverage per Run", "scaffolds"),
    ]:
        plt.figure(figsize=(10, 4))
        plt.bar(labels, coverages, color=color)
        plt.ylabel("Coverage")
        plt.title(title)
        plt.xticks(rotation=45, ha="right")
        plt.tight_layout()

        if save:
            output_dir = Path(output_dir)
            output_dir.mkdir(parents=True, exist_ok=True)
            plt.savefig(output_dir / f"{suffix}_coverage.png", dpi=300)

        plt.show()


def main():
    print("--- Starting Optimization Analysis (Seaborn Edition) ---")

    set_publication_style()

    all_files = glob.glob(str(INPUT_DIR / "**/*.csv"), recursive=True)
    if not all_files:
        print(f"Error: No files found in {INPUT_DIR}")
        return

    dfs = []
    for f in all_files:
        try:
            temp = pd.read_csv(f)
            temp["source"] = os.path.basename(f)
            parts = temp["source"].iloc[0].split("_")

            if "weighted" in temp["source"].iloc[0]:
                mode, run_id = "dbg_weighted", parts[3]
            elif "multimodal" in temp["source"].iloc[0]:
                mode, run_id = "multimodal_dbg", parts[3]
            else:
                mode, run_id = "greedy", parts[2] if len(parts) > 2 else "unknown"

            temp["mode"] = mode
            temp["run_id"] = run_id
            if "chain" not in temp.columns:
                temp["chain"] = "N/A"
            dfs.append(temp)
        except Exception as e:
            print(f"Skipped {f}: {e}")

    if not dfs:
        return
    full_df = pd.concat(dfs, ignore_index=True)

    full_df = normalize_composite_score(full_df)
    full_df["category"] = full_df["run_id"].apply(get_category)
    full_df["seq_type"] = full_df["refine_rounds"].apply(lambda x: "Scaffolds" if x > 0 else "Contigs")

    mask = (full_df["chain"] != "N/A") & (full_df["chain"].notna())
    full_df["display_name"] = full_df["run_id"]
    full_df.loc[mask, "display_name"] = full_df.loc[mask, "run_id"] + " (" + full_df.loc[mask, "chain"] + ")"

    print("Generating Heatmaps...")
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    unique_combinations = full_df[["category", "mode", "display_name", "seq_type"]].drop_duplicates()

    for _, row in unique_combinations.iterrows():
        subset = full_df[
            (full_df["category"] == row["category"])
            & (full_df["mode"] == row["mode"])
            & (full_df["display_name"] == row["display_name"])
            & (full_df["seq_type"] == row["seq_type"])
        ]
        if subset.empty:
            continue

        target_y, target_x = HEATMAP_CONFIG.get(row["mode"], HEATMAP_CONFIG["default"])
        valid_y = [c for c in target_y if c in subset.columns]
        valid_x = [c for c in target_x if c in subset.columns]
        if not valid_y or not valid_x:
            continue

        safe_name = str(row["display_name"]).replace(" ", "_").replace("(", "").replace(")", "")
        out_subdir = OUTPUT_DIR / row["mode"]

        for metric in ["composite_score", "coverage"]:
            if metric not in subset.columns:
                continue
            fname = f"{safe_name}_{row['seq_type']}_heatmap_{metric.replace('_', '')}.svg"
            title = ""

            theme_colors = THEME_MAP_DATA.get(row["category"], THEME_MAP_DATA["Other"])

            plot_grid_search_clustermap(
                subset, valid_y, valid_x, theme_colors, metric, title, output_file=out_subdir / fname
            )

    print("Selecting best candidates...")
    group_cols = ["display_name", "mode", "seq_type"]
    full_df_sorted = full_df.sort_values("composite_score", ascending=False)
    df_best = full_df_sorted.drop_duplicates(subset=group_cols, keep="first").copy()

    print("Generating Summary Tables...")
    SUMMARY_DIR.mkdir(parents=True, exist_ok=True)
    for mode in df_best["mode"].unique():
        df_mode = df_best[df_best["mode"] == mode].copy()
        if df_mode.empty:
            continue
        save_detailed_rankings(df_mode, SUMMARY_DIR, mode)

    print("--- Analysis Complete. ---")


if __name__ == "__main__":
    main()
