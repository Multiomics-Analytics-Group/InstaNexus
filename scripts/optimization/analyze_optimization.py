#!/usr/bin/env python

r""" Optimization Analysis Script for InstaNexus.
Aggregates Grid Search results, computes global composite scores, 
and generates manuscript-ready figures comparing Contigs vs Scaffolds.

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
__date__ = 10 Dec 2025
__maintainer__ = Marco Reverenna
__email__ = marcor@dtu.dk
__status__ = Dev
"""

import os
import glob
import pandas as pd
import numpy as np
import plotly.graph_objects as go
from sklearn.preprocessing import MinMaxScaler
from pathlib import Path

INPUT_DIR = Path("outputs/_grid_search")
OUTPUT_DIR = Path("outputs/_optimisation_figures")
SUMMARY_DIR = Path("outputs/_summary_tables")

# Reviewer-approved palette (Light=Contig, Dark=Scaffold)
COLORS = {
    "Antibodies": ("#8dd3c7", "#1b9e77"),
    "Nanobodies": ("#a6cee3", "#1f78b4"),
    "Binders":    ("#FFE066", "#E6AB02"),
    "BSA":        ("#fdbb84", "#e34a33"),
    "Other":      ("lightgrey", "darkgrey")
}

def get_category(run_name):
    """Maps run names to biological categories."""
    l = run_name.lower()
    if "bsa" in l: return "BSA"
    if "nb" in l: return "Nanobodies"
    if "bind" in l: return "Binders"
    if "ma" in l or "pama" in l: return "Antibodies"
    return "Other"

def classify_sequence_type(row):
    """
    Distinguishes Contigs from Scaffolds based on refinement rounds.
    refine_rounds == 0 -> Contigs (Raw Assembly)
    refine_rounds > 0  -> Scaffolds (Refined/Merged)
    """
    if row.get("refine_rounds", 0) == 0:
        return "Contigs"
    else:
        return "Scaffolds"

def compute_composite_score_global(df):
    """
    Computes a globally normalized composite score.
    Strategy: 'Aggressive Consolidation'
    - Coverage (35%): Dominant factor.
    - N50 (25%) + Scaffolds (25%): Structural quality.
    - Identity (15%): Sequence quality.
    """
    metrics = ["coverage", "N50", "scaffolds_count", "mean_identity"]
    weights = {"coverage": 0.35, "N50": 0.25, "scaffolds_count": 0.25, "mean_identity": 0.15}
    
    # Validation
    if df.empty or not all(m in df.columns for m in metrics):
        print("Warning: Missing metrics for scoring. Skipping.")
        return df

    # Normalize (0-1)
    scaler = MinMaxScaler()
    # FillNA with 0 to prevent errors on failed runs
    df_metrics = df[metrics].fillna(0)
    scaled = scaler.fit_transform(df_metrics)
    df_scaled = pd.DataFrame(scaled, columns=metrics, index=df.index)

    # Invert scaffolds count (Lower is better -> Higher Score)
    df_scaled["scaffolds_count"] = 1 - df_scaled["scaffolds_count"]

    # Weighted Sum
    df["composite_score"] = df_scaled[list(weights.keys())].dot(pd.Series(weights))
    
    return df

def plot_contig_vs_scaffold(df_best, category, mode, output_path):
    """
    Generates a grouped barplot: Best Contig vs Best Scaffold per run.
    """
    # Filter by category and mode (e.g., dbg_weighted)
    df_cat = df_best[
        (df_best["category"] == category) & 
        (df_best["mode"] == mode)
    ].copy()
    
    if df_cat.empty:
        return

    # Sort runs for consistent plotting
    df_cat = df_cat.sort_values("run")
    runs = df_cat["run"].unique()
    
    # Get Colors
    c_contig, c_scaffold = COLORS.get(category, COLORS["Other"])

    fig = go.Figure()

    # Trace 1: Contigs
    df_contigs = df_cat[df_cat["seq_type"] == "Contigs"]
    # Reindex ensures we have a value (or 0) for every run, even if data is missing
    y_contigs = df_contigs.set_index("run")["composite_score"].reindex(runs, fill_value=0)
    
    fig.add_trace(go.Bar(
        x=runs, y=y_contigs,
        name="Contigs", marker_color=c_contig
    ))

    # Trace 2: Scaffolds
    df_scaffolds = df_cat[df_cat["seq_type"] == "Scaffolds"]
    y_scaffolds = df_scaffolds.set_index("run")["composite_score"].reindex(runs, fill_value=0)
    
    fig.add_trace(go.Bar(
        x=runs, y=y_scaffolds,
        name="Scaffolds", marker_color=c_scaffold
    ))

    # Layout updates
    fig.update_layout(
        title=f"<b>{category}</b> Performance ({mode})",
        xaxis_title="",
        yaxis_title="Composite Score",
        yaxis=dict(range=[0, 1.05], showgrid=True, gridcolor='lightgrey'),
        barmode='group',
        template="plotly_white",
        width=800, height=500,
        legend=dict(
            orientation="h", 
            yanchor="bottom", y=1.02, 
            xanchor="right", x=1
        ),
        font=dict(family="Arial", size=14)
    )
    
    # Save
    filename = f"{category}_{mode}_benchmark.svg"
    fig.write_image(str(output_path / filename))
    print(f"   Saved: {filename}")

def summarize_best_params(df):
    """Extracts the optimal parameters for the best runs."""
    # Group by Category + Mode + Type to find the "Winning Configuration"
    # We take the mode (most frequent) parameter set across all runs in that category
    
    summary_data = []
    param_cols = ["kmer_size", "min_overlap", "size_threshold", "min_weight", "fdr", "alpha_cov"]
    
    # Filter relevant columns
    existing_params = [p for p in param_cols if p in df.columns]

    for (cat, mode, seq_type), sub_df in df.groupby(["category", "mode", "seq_type"]):
        row = {"Category": cat, "Mode": mode, "Type": seq_type}
        
        # Calculate mode for each parameter
        for p in existing_params:
            try:
                row[p] = sub_df[p].mode()[0]
            except:
                row[p] = "-"
        summary_data.append(row)
        
    return pd.DataFrame(summary_data)

def main():
    print("--- Starting Optimization Analysis ---")
    
    # 1. Load Data
    all_files = glob.glob(str(INPUT_DIR / "*.csv"))
    if not all_files:
        print(f"Error: No files found in {INPUT_DIR}")
        return

    dfs = []
    for f in all_files:
        try:
            temp = pd.read_csv(f)
            fname = os.path.basename(f)
            # Expected format: grid_MODE_RUNNAME_TIMESTAMP.csv
            parts = fname.split("_")
            if len(parts) >= 3:
                temp["mode"] = parts[1]
                temp["run"] = parts[2]
            dfs.append(temp)
        except Exception as e:
            print(f"Skipped {f}: {e}")

    full_df = pd.concat(dfs, ignore_index=True)    
    full_df["category"] = full_df["run"].apply(get_category)
    full_df["seq_type"] = full_df.apply(classify_sequence_type, axis=1)

    mask = (full_df["chain"].notna()) & (full_df["chain"] != "N/A") & (full_df["chain"] != "") & (full_df["chain"].str.lower() != "nan")
    full_df.loc[mask, "run"] = full_df.loc[mask, "run"] + " (" + full_df.loc[mask, "chain"] + ")"
    
    print("Computing Composite Scores...")
    full_df = compute_composite_score_global(full_df)

    best_runs = full_df.loc[
        full_df.groupby(["run", "mode", "seq_type"])["composite_score"].idxmax()
    ]
    
    print("Generating plots...")
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    
    modes_to_plot = best_runs["mode"].unique()
    categories = best_runs["category"].unique()

    for mode in modes_to_plot:
        for cat in categories:
            plot_contig_vs_scaffold(best_runs, cat, mode, OUTPUT_DIR)

    # 6. Save Tables
    print("Saving Summary Tables...")
    SUMMARY_DIR.mkdir(parents=True, exist_ok=True)
    
    # A. Best Parameters Table
    param_table = summarize_best_params(best_runs)
    param_table.to_csv(SUMMARY_DIR / "optimal_parameters.csv", index=False)
    
    # B. Full Ranked Results (for supplementary)
    full_df.sort_values("composite_score", ascending=False).to_csv(
        SUMMARY_DIR / "full_benchmark_ranked.csv", index=False
    )

    print("--- Analysis completed ---")

if __name__ == "__main__":
    main()