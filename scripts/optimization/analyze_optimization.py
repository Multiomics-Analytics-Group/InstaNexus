#!/usr/bin/env python

"""
Optimization Analysis Script for InstaNexus.

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
__date__ = 11 Dec 2025
__maintainer__ = Marco Reverenna
__email__ = marcor@dtu.dk
__status__ = Dev

This script processes grid search results to generate manuscript-ready figures and tables:
1. Multi-Parameter Heatmaps (Clustermaps): Visualizes Composite Score and Coverage across
   parameter space (e.g., k-mer size vs overlap), separated by assembly mode.
2. Optimal Parameter Tables: Calculates the mode (most frequent value) of optimal parameters
   per category to identify robust settings (e.g., Supplementary Table 4).
3. Detailed Ranking Tables: Generates comprehensive CSVs with raw and scaled metrics
   (Coverage, N50, Total Sequences) for the best performing samples.
"""

import glob
import os
from pathlib import Path

import pandas as pd
import plotly.graph_objects as go
from sklearn.preprocessing import MinMaxScaler

INPUT_DIR = Path("outputs/_grid_search")
OUTPUT_DIR = Path("outputs/_optimization_figures")
SUMMARY_DIR = Path("outputs/_summary_tables")

ABBREV = {
    "kmer_size": "k",
    "min_overlap": "mo",
    "size_threshold": "st",
    "min_weight": "mw",
    "fdr": "fdr",
    "alpha_cov": "ac",
    "refine_rounds": "rr",
}

HEATMAP_CONFIG = {
    "greedy": (["min_overlap"], ["size_threshold", "fdr"]),
    "dbg_weighted": (["kmer_size", "min_overlap"], ["size_threshold", "fdr"]),
    "multimodal_dbg": (["kmer_size", "min_weight"], ["size_threshold", "fdr"]),
    "default": (["min_overlap"], ["fdr"]),
}

THEME_MAP = {
    "BSA": [[0.0, "#fdbb84"], [0.6, "#fdbb84"], [1.0, "#e34a33"]],
    "Antibodies": [[0.0, "#8dd3c7"], [0.6, "#8dd3c7"], [1.0, "#1b9e77"]],
    "Nanobodies": [[0.0, "#a6cee3"], [0.6, "#a6cee3"], [1.0, "#1f78b4"]],
    "Binders": [[0.0, "#F8ECA9"], [0.6, "#F8ECA9"], [1.0, "#FFC72C"]],
    "Other": [[0.0, "lightgrey"], [0.6, "lightgrey"], [1.0, "black"]],
}

PARAM_CONFIG = {
    "greedy": ["fdr", "min_overlap", "size_threshold", "refine_rounds"],
    "dbg_weighted": ["fdr", "kmer_size", "min_weight", "min_overlap", "size_threshold", "refine_rounds"],
    "multimodal_dbg": ["fdr", "kmer_size", "min_weight", "size_threshold", "refine_rounds"],
}


def get_category(run_name):
    """Maps run ID to biological category."""
    name_lower = run_name.lower()
    if "bsa" in name_lower:
        return "BSA"
    if "nb" in name_lower:
        return "Nanobodies"
    if "bind" in name_lower:
        return "Binders"
    if "ma" in name_lower or "pama" in name_lower:
        return "Antibodies"
    return "Other"


def classify_sequence_type(row):
    """Refine rounds == 0 -> Contigs; > 0 -> Scaffolds."""
    if row.get("refine_rounds", 0) == 0:
        return "Contigs"
    else:
        return "Scaffolds"


def plot_grid_search_clustermap(
    df, index_cols, column_cols, theme, value_col, title="", aggfunc="max", output_file=None
):
    try:
        pivot = df.pivot_table(values=value_col, index=index_cols, columns=column_cols, aggfunc=aggfunc)
    except KeyError:
        return

    pivot = pivot.sort_index(level=index_cols).sort_index(axis=1, level=column_cols)
    row_labels = [
        ", ".join(
            f"{ABBREV.get(c,c)}={v}" for c, v in zip(index_cols, idx if isinstance(idx, tuple) else [idx], strict=False)
        )
        for idx in pivot.index
    ]
    col_labels = [
        ", ".join(
            f"{ABBREV.get(c,c)}={v}"
            for c, v in zip(column_cols, col if isinstance(col, tuple) else [col], strict=False)
        )
        for col in pivot.columns
    ]
    pivot = pivot.fillna(0)

    fig = go.Figure(
        data=[
            go.Heatmap(
                z=pivot.values,
                x=col_labels,
                y=row_labels,
                colorscale=theme,
                zmin=0,
                zmax=1,
                showscale=True,
                colorbar=dict(title=value_col, len=0.75, thickness=20),
            )
        ]
    )

    n_rows, n_cols = pivot.shape
    shapes = []
    for i in range(n_rows + 1):
        shapes.append(
            dict(type="line", x0=-0.5, x1=n_cols - 0.5, y0=i - 0.5, y1=i - 0.5, line=dict(color="white", width=2))
        )
    for j in range(n_cols + 1):
        shapes.append(
            dict(type="line", x0=j - 0.5, x1=j - 0.5, y0=-0.5, y1=n_rows - 0.5, line=dict(color="white", width=2))
        )

    fig.update_layout(
        width=950,
        height=850,
        title=dict(text=title, x=0.5),
        xaxis=dict(
            tickangle=-45, showgrid=False, zeroline=False, title=", ".join([ABBREV.get(c, c) for c in column_cols])
        ),
        yaxis=dict(
            showgrid=False,
            zeroline=False,
            autorange="reversed",
            title=", ".join([ABBREV.get(c, c) for c in index_cols]),
        ),
        plot_bgcolor="white",
        paper_bgcolor="white",
        shapes=shapes,
        font=dict(family="Arial", size=14),
    )
    if output_file:
        output_file.parent.mkdir(parents=True, exist_ok=True)
        fig.write_image(str(output_file), format="svg", scale=2)


def compute_and_add_scaled_metrics(df):
    """Computes scaled metrics (0-1) for the final report to match Composite Score logic."""
    metrics_map = {
        "total_sequences": "total_sequences_scaled",
        "N50": "N50_scaled",
        "coverage": "coverage_scaled",
        "mean_identity": "mean_identity_scaled",
    }
    existing_metrics = [m for m in metrics_map.keys() if m in df.columns]

    if not existing_metrics:
        return df

    scaler = MinMaxScaler()
    for m in existing_metrics:
        values = df[[m]].fillna(0)
        scaled = scaler.fit_transform(values)

        if m == "total_sequences":
            df[metrics_map[m]] = 1 - scaled
        else:
            df[metrics_map[m]] = scaled
    return df


def generate_optimal_params_table(df_best, output_dir, mode_name):
    """Calculates the MODE (most frequent value) of optimal parameters per category."""
    target_params = PARAM_CONFIG.get(mode_name, [])
    available_params = [p for p in target_params if p in df_best.columns]

    summary_rows = []
    for cat, sub in df_best.groupby("category"):
        row = {"Category": cat}
        for p in available_params:
            try:
                row[p] = sub[p].mode().iloc[0]
            except Exception:
                row[p] = "-"
        summary_rows.append(row)

    summary_df = pd.DataFrame(summary_rows)
    out_path = output_dir / mode_name / f"optimal_params_{mode_name}.csv"
    out_path.parent.mkdir(parents=True, exist_ok=True)
    summary_df.to_csv(out_path, index=False)


def save_detailed_rankings(df_best, output_dir, mode_name):
    """Saves detailed ranking tables with raw and scaled metrics."""
    for cat, sub in df_best.groupby("category"):
        if "scaffolds_count" in sub.columns:
            sub = sub.rename(columns={"scaffolds_count": "total_sequences"})

        sub = compute_and_add_scaled_metrics(sub.copy())

        sub = sub.rename(columns={"display_name": "sample", "mode": "assembly_method"})

        desired_order = [
            "category",
            "sample",
            "assembly_method",
            "source",
            "total_sequences",
            "total_sequences_scaled",
            "N50",
            "N50_scaled",
            "coverage",
            "coverage_scaled",
            "composite_score",
        ]
        if "mean_identity" in sub.columns:
            desired_order.insert(-1, "mean_identity")
            desired_order.insert(-1, "mean_identity_scaled")

        final_cols = [c for c in desired_order if c in sub.columns]

        # 5. Sort and Save
        sub_sorted = sub[final_cols].sort_values("composite_score", ascending=False)
        numeric_cols = sub_sorted.select_dtypes(include=["float64", "float32"]).columns
        sub_sorted[numeric_cols] = sub_sorted[numeric_cols].round(3)

        out_path = output_dir / mode_name / f"best_results_{cat}_{mode_name}.csv"
        out_path.parent.mkdir(parents=True, exist_ok=True)
        sub_sorted.to_csv(out_path, index=False)


def main():
    print("--- Starting Optimization Analysis ---")

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

            if len(parts) >= 4 and parts[1] == "dbg" and parts[2] == "weighted":
                mode, run_id = "dbg_weighted", parts[3]
            elif len(parts) >= 4 and parts[1] == "multimodal":
                mode, run_id = "multimodal_dbg", parts[3]
            elif len(parts) >= 3:
                mode, run_id = parts[1], parts[2]
            else:
                mode, run_id = "unknown", "unknown"

            temp["mode"], temp["run_id"] = mode, run_id
            if "chain" not in temp.columns:
                temp["chain"] = "N/A"
            dfs.append(temp.fillna({"chain": "N/A"}))
        except Exception as e:
            print(f"Skipped {f}: {e}")

    if not dfs:
        return
    full_df = pd.concat(dfs, ignore_index=True)
    full_df["category"] = full_df["run_id"].apply(get_category)
    full_df["seq_type"] = full_df.apply(classify_sequence_type, axis=1)

    mask = (full_df["chain"] != "N/A") & (full_df["chain"] != "")
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

        safe_name = row["display_name"].replace(" ", "_").replace("(", "").replace(")", "")
        out_subdir = OUTPUT_DIR / row["mode"]

        for metric in ["composite_score", "coverage"]:
            if metric not in subset.columns:
                continue
            fname = f"{safe_name}_{row['seq_type']}_heatmap_{metric.replace('_', '')}.svg"
            title = f"{row['display_name']} - {row['seq_type']} ({row['mode']})<br>{metric.replace('_', ' ').title()}"
            plot_grid_search_clustermap(
                subset, valid_y, valid_x, THEME_MAP.get(row["category"]), metric, title, output_file=out_subdir / fname
            )

    print("Generating Summary Tables...")
    SUMMARY_DIR.mkdir(parents=True, exist_ok=True)

    df_best = full_df.loc[full_df.groupby(["display_name", "mode", "seq_type"])["composite_score"].idxmax()].copy()

    for mode in df_best["mode"].unique():
        df_mode = df_best[df_best["mode"] == mode].copy()
        if df_mode.empty:
            continue

        generate_optimal_params_table(df_mode, SUMMARY_DIR, mode)
        save_detailed_rankings(df_mode, SUMMARY_DIR, mode)

    print("--- Analysis Complete. Results saved to outputs/ ---")


if __name__ == "__main__":
    main()
