#!/usr/bin/env python

"""
Optimization Analysis Script for InstaNexus.
Focus: Grid Search Analysis, Parameter Selection, and Ranking.
Output: Publication-ready tables (Raw Metrics only) and Heatmaps.

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
"""

import glob
import os
from pathlib import Path
import numpy as np
import pandas as pd
import plotly.graph_objects as go

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

THEME_MAP = {
    "BSA": [[0.0, "#fdbb84"], [0.7, "#fdbb84"], [1.0, "#e34a33"]],
    "Antibodies": [[0.0, "#8dd3c7"], [0.7, "#8dd3c7"], [1.0, "#1b9e77"]],
    "Nanobodies": [[0.0, "#a6cee3"], [0.7, "#a6cee3"], [1.0, "#1f78b4"]],
    "Binders": [[0.0, "#F8ECA9"], [0.7, "#F8ECA9"], [1.0, "#FFC72C"]],
    "Other": [[0.0, "lightgrey"], [0.7, "lightgrey"], [1.0, "black"]],
}

HEATMAP_CONFIG = {
    "greedy": (["min_overlap"], ["size_threshold", "fdr"]),
    "dbg_weighted": (["kmer_size", "min_overlap"], ["size_threshold", "fdr"]),
    "multimodal_dbg": (["kmer_size", "min_weight"], ["size_threshold", "fdr"]),
    "default": (["min_overlap"], ["fdr"]),
}

def get_category(run_name):
    """Maps run ID to biological category."""
    name_lower = str(run_name).lower()
    if "bsa" in name_lower: return "BSA"
    if "nb" in name_lower: return "Nanobodies"
    if "bind" in name_lower: return "Binders"
    if "ma" in name_lower or "pama" in name_lower: return "Antibodies"
    return "Other"

def format_assembly_method(row):
    """
    Distinguishes between raw assembly (Contigs) and refined (Scaffolds).
    """
    method = row.get("mode", "unknown")
    seq_type = row.get("seq_type", "Contigs")
    return f"{method} ({seq_type})"

def normalize_composite_score(df):
    """
    Normalizes composite_score to 0-1 range if it appears to be in 0-1000 scale.
    """
    if "composite_score" not in df.columns:
        return df
    
    # Heuristic: if score > 1.5, assume it's out of 1000.
    if df["composite_score"].max() > 1.5:
        print(">> Normalizing composite_score (dividing by 1000) for consistency.")
        df["composite_score"] = df["composite_score"] / 1000.0
    
    return df

def save_detailed_rankings(df_best, output_dir, mode_name):
    """
    Saves detailed ranking tables with RAW metrics only.
    Source column removed.
    """
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
            "composite_score"
        ]
        
        final_cols = [c for c in desired_order if c in sub.columns]
        
        final_df = sub[final_cols].sort_values(by=["sample", "assembly_method"], ascending=[True, True])
        
        numeric_cols = final_df.select_dtypes(include=["float64", "float32"]).columns
        final_df[numeric_cols] = final_df[numeric_cols].round(3)

        out_path = output_dir / mode_name / f"best_results_{cat}_{mode_name}.csv"
        out_path.parent.mkdir(parents=True, exist_ok=True)
        final_df.to_csv(out_path, index=False)
        print(f"Saved cleaned table: {out_path}")

def plot_grid_search_clustermap(df, index_cols, column_cols, theme, value_col, title="", output_file=None):
    """Generates heatmaps for parameter exploration with White Grid."""
    try:
        pivot = df.pivot_table(values=value_col, index=index_cols, columns=column_cols, aggfunc="max")
    except KeyError:
        return

    pivot = pivot.sort_index(level=index_cols).sort_index(axis=1, level=column_cols).fillna(0)
    
    row_labels = [", ".join(f"{ABBREV.get(c,c)}={v}" for c, v in zip(index_cols, idx if isinstance(idx, tuple) else [idx])) for idx in pivot.index]
    col_labels = [", ".join(f"{ABBREV.get(c,c)}={v}" for c, v in zip(column_cols, col if isinstance(col, tuple) else [col])) for col in pivot.columns]

    fig = go.Figure(data=[go.Heatmap(
        z=pivot.values, x=col_labels, y=row_labels,
        colorscale=theme, zmin=0, zmax=1, showscale=True,
        colorbar=dict(title=value_col, len=0.75, thickness=20)
    )])
    
    n_rows, n_cols = pivot.shape
    shapes = []
    for i in range(n_rows + 1):
        shapes.append(dict(type="line", x0=-0.5, x1=n_cols-0.5, y0=i-0.5, y1=i-0.5, 
                           line=dict(color="white", width=2)))
    for j in range(n_cols + 1):
        shapes.append(dict(type="line", x0=j-0.5, x1=j-0.5, y0=-0.5, y1=n_rows-0.5, 
                           line=dict(color="white", width=2)))
    
    fig.update_layout(
        title=dict(text=title, x=0.5),
        width=950, height=850,
        plot_bgcolor="white",
        paper_bgcolor="white",
        shapes=shapes, # Apply the grid
        font=dict(family="Arial", size=14),
        xaxis=dict(tickangle=-45, showgrid=False, zeroline=False),
        yaxis=dict(showgrid=False, zeroline=False, autorange="reversed")
    )
    if output_file:
        output_file.parent.mkdir(parents=True, exist_ok=True)
        fig.write_image(str(output_file), format="svg")

def main():
    print("--- Starting Optimization Analysis (Final Supervisor Edition) ---")

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
            if "chain" not in temp.columns: temp["chain"] = "N/A"
            dfs.append(temp)
        except Exception as e:
            print(f"Skipped {f}: {e}")

    if not dfs: return
    full_df = pd.concat(dfs, ignore_index=True)

    full_df = normalize_composite_score(full_df)
    full_df["category"] = full_df["run_id"].apply(get_category)
    
    # Refine Rounds > 0 implies Scaffolding/Refinement occurred.
    full_df["seq_type"] = full_df["refine_rounds"].apply(lambda x: "Scaffolds" if x > 0 else "Contigs")

    mask = (full_df["chain"] != "N/A") & (full_df["chain"].notna())
    full_df["display_name"] = full_df["run_id"]
    full_df.loc[mask, "display_name"] = full_df.loc[mask, "run_id"] + " (" + full_df.loc[mask, "chain"] + ")"

    print("Generating Heatmaps...")
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    unique_combinations = full_df[["category", "mode", "display_name", "seq_type"]].drop_duplicates()
    
    for _, row in unique_combinations.iterrows():
        subset = full_df[
            (full_df["category"] == row["category"]) & 
            (full_df["mode"] == row["mode"]) & 
            (full_df["display_name"] == row["display_name"]) &
            (full_df["seq_type"] == row["seq_type"])
        ]
        if subset.empty: continue
        
        target_y, target_x = HEATMAP_CONFIG.get(row["mode"], HEATMAP_CONFIG["default"])
        valid_y = [c for c in target_y if c in subset.columns]
        valid_x = [c for c in target_x if c in subset.columns]
        if not valid_y or not valid_x: continue

        safe_name = str(row["display_name"]).replace(" ", "_").replace("(", "").replace(")", "")
        out_subdir = OUTPUT_DIR / row["mode"]
        
        for metric in ["composite_score", "coverage"]:
            if metric not in subset.columns: continue
            fname = f"{safe_name}_{row['seq_type']}_heatmap_{metric.replace('_', '')}.svg"
            title = f"{row['display_name']} - {row['seq_type']} ({row['mode']})<br>{metric.replace('_', ' ').title()}"
            plot_grid_search_clustermap(subset, valid_y, valid_x, THEME_MAP.get(row["category"]), metric, title, output_file=out_subdir/fname)

    print("Selecting best candidates...")
    group_cols = ["display_name", "mode", "seq_type"]
    full_df_sorted = full_df.sort_values("composite_score", ascending=False)
    df_best = full_df_sorted.drop_duplicates(subset=group_cols, keep="first").copy()

    print("Generating Summary Tables...")
    SUMMARY_DIR.mkdir(parents=True, exist_ok=True)

    for mode in df_best["mode"].unique():
        df_mode = df_best[df_best["mode"] == mode].copy()
        if df_mode.empty: continue
        save_detailed_rankings(df_mode, SUMMARY_DIR, mode)

    print("--- Analysis Complete. ---")

if __name__ == "__main__":
    main()