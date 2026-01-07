#!/usr/bin/env python

r"""Visualization module for InstaNexus.

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
__date__ = 25 Nov 2025
__maintainer__ = Marco Reverenna
__email__ = marcor@dtu.dk
__status__ = Dev
"""

import os
import json
import Bio.SeqIO
import matplotlib
matplotlib.use('Agg')
import matplotlib.patches as patches
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
import seaborn as sns

from tqdm import tqdm
from pathlib import Path
from instanexus import helpers


def set_publication_style():
    """
    Docstring for set_publication_style
    """
    sns.set_theme(style="ticks")
    
    plt.rcParams.update({
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
        "axes.linewidth": 1,
    })


def get_figsize(width_ratio=1, total_width_inch=14.0):
    """
    Calculates figsize based on a standard A4 width.
    Ratios: 
    3: 3x1 (Full width, short height)
    2: 2x1 (Medium)
    1: 1x1 (Square)  

    :param width_ratio: Description
    :param total_width_inch: Description
    
    """

    ratios = {
        1: (1, 1),
        2: (2, 1),
        3: (3, 1),
    }
    
    w_mult, h_mult = ratios.get(width_ratio, (1, 1))
    
    actual_width = (total_width_inch / 3) * w_mult
    actual_height = (actual_width / w_mult)
    
    return (actual_width, actual_height)


def plot_map_unmap_distribution(df, run, folder, conf_lim, ratio=1, unique_peps=True, title=False):
    """
    Docstring for plot_map_unmap_distribution
    
    :param df: Description
    :param run: Description
    :param folder: Description
    :param conf_lim: Description
    :param ratio: Description
    :param unique_peps: Description
    :param title: Description
    """
    set_publication_style()
    
    df_filtered = df[df["conf"] >= conf_lim].copy()
    if unique_peps:
        df_filtered = df_filtered.sort_values("conf", ascending=False).drop_duplicates(subset=["cleaned_preds"])

    _, ax = plt.subplots(figsize=get_figsize(width_ratio=ratio))
    
    c_maps, c_not = "#4A90E2", "#FF9F1C"
    bins = np.arange(0, 1.01, 0.01)

    sns.histplot(
        data=df_filtered, x="conf", hue="mapped",
        palette={True: c_maps, False: c_not},
        bins=bins, element="bars", fill=True, alpha=0.8,
        multiple="layer", ax=ax, legend=False
    )

    overlap_color = "#8c7c6d"
    handles = [
        patches.Patch(color=c_maps, alpha=0.8, label='Maps in protein'),
        patches.Patch(color=c_not, alpha=0.8, label='Not in protein'),
        patches.Patch(color=overlap_color, alpha=0.8, label='Overlap')
    ]
    ax.legend(
        handles=handles, loc='lower center', 
        bbox_to_anchor=(0.5, 1.02), ncol=3, borderaxespad=0.,
        frameon=False
    )

    ax.set_xlim(0, 1)
    ax.set_yscale("log")
    ax.set_xlabel("Confidence (calibrated)")
    ax.set_ylabel("Unique peptides" if unique_peps else "PSMs counts")
    
    if title:
        ax.set_title("Sequence distribution by confidence", pad=30)

    sns.despine()

    Path(folder).mkdir(parents=True, exist_ok=True)
    suffix = "unique" if unique_peps else "all"
    output_path = f"{folder}/fig2a_{run}_distribution_psm_{suffix}.svg"
    
    plt.savefig(output_path, format="svg", bbox_inches="tight")
    plt.close()

    return output_path



def plot_ratios_mapped_unmapped_threshold_set(run, df, reference, folder, ratio=1):
    """
    Plots confidence and FDR distributions using a square ratio (1:1).
    
    :param run: Description
    :param df: Description
    :param reference: Description
    :param folder: Description
    :param ratio: Description
    """
    
    
    set_publication_style()
    
    df_processed = df.copy()
    df_processed["mapped_status"] = df_processed["cleaned_preds"].apply(
        lambda x: "Maps to protein" if isinstance(x, str) and x in reference else "Not in protein"
    )
    
    palette = {
        "Maps to protein": "#4A90E2",
        "Not in protein": "#FF9F1C"
    }
    hue_order = ["Not in protein", "Maps to protein"]

    df_conf = df_processed[df_processed["conf"] >= 0.90].copy()
    
    _, ax1 = plt.subplots(figsize=get_figsize(width_ratio=ratio))
    
    sns.histplot(
        data=df_conf,
        x="conf",
        hue="mapped_status",
        multiple="stack",
        palette=palette,
        hue_order=hue_order,
        bins=np.arange(0.90, 1.01, 0.02),
        shrink=0.8,
        edgecolor="none",
        alpha=0.8,
        ax=ax1,
        legend=False
    )

    ax1.set_xlabel("Confidence (calibrated)")
    ax1.set_ylabel("PSMs counts")
    ax1.set_xlim(0.9, 1.0)
    
    ax1.legend(
        handles=[patches.Patch(color=palette[label], label=label) for label in hue_order],
        loc='lower center', 
        bbox_to_anchor=(0.5, 1.02), 
        ncol=2,
        frameon=False
    )
    
    sns.despine()
    
    if folder:
        Path(folder).mkdir(parents=True, exist_ok=True)
        plt.savefig(f"{folder}/fig2b_{run}_confidence_barplot.svg", format="svg", bbox_inches="tight")
    plt.show()

    df_fdr = df_processed[df_processed["psm_q_value"] <= 0.20].copy()
    
    fdr_bins = [0, 0.01, 0.05, 0.10, 0.20]
    fdr_labels = ["0-1%", "1-5%", "5-10%", "10-20%"]
    df_fdr["fdr_interval"] = pd.cut(df_fdr["psm_q_value"], bins=fdr_bins, labels=fdr_labels, include_lowest=True)

    _, ax2 = plt.subplots(figsize=get_figsize(width_ratio=ratio))
    
    sns.histplot(
        data=df_fdr,
        x="fdr_interval",
        hue="mapped_status",
        multiple="stack",
        palette=palette,
        hue_order=hue_order,
        shrink=0.7,
        edgecolor="none",
        alpha=0.8,
        discrete=True,
        ax=ax2,
        legend=False
    )

    ax2.set_xlabel("FDR interval (q-value)")
    ax2.set_ylabel("PSMs counts")
    
    ax2.legend(
        handles=[patches.Patch(color=palette[label], label=label) for label in hue_order],
        loc='lower center', 
        bbox_to_anchor=(0.5, 1.02), 
        ncol=2,
        frameon=False
    )

    sns.despine()

    if folder:
        plt.savefig(f"{folder}/fig2e_{run}_fdr_barplots.svg", format="svg", bbox_inches="tight")
    plt.show()



def plot_fdr_mapped_unmapped_ratio_across_all_range(run, df, reference, folder, ratio=1, unique_peps=False):
    """
    Plots the Mapped/Unmapped ratio across the full FDR range (0 to 1).
    
    :param run: Description
    :param df: Description
    :param reference: Description
    :param folder: Description
    :param ratio: Description
    :param unique_peps: Description
    """

    set_publication_style()
    
    df_copy = df.copy()
    
    df_copy["mapped_status"] = df_copy["cleaned_preds"].apply(
        lambda x: "mapped" if isinstance(x, str) and x in reference else "unmapped"
    )
    
    if unique_peps:
        df_copy = df_copy.sort_values("conf", ascending=False).drop_duplicates(subset=["cleaned_preds"])
    
    thresholds = np.linspace(0.01, 1.0, 50)
    data_points = []

    for t in thresholds:
        subset = df_copy[df_copy["psm_q_value"] <= t]
        if len(subset) == 0: continue
            
        counts = subset["mapped_status"].value_counts()
        n_mapped = counts.get("mapped", 0)
        n_unmapped = counts.get("unmapped", 0)
        
        ratio_val = n_mapped / n_unmapped if n_unmapped > 0 else np.nan
        data_points.append({"fdr": t, "ratio": ratio_val})

    plot_df = pd.DataFrame(data_points).dropna()

    _, ax = plt.subplots(figsize=get_figsize(width_ratio=ratio))

    VIBRANT_BLUE = "#4A90E2"
    VIBRANT_ORANGE = "#FF9F1C"
    
    sns.lineplot(
        data=plot_df, x='fdr', y='ratio', 
        color=VIBRANT_BLUE, linewidth=1.5, ax=ax, zorder=2
    )

    target_thresholds = [0.01, 0.05, 0.1, 0.2, 0.5, 1.0]
    highlights = []
    for target in target_thresholds:
        idx = (np.abs(plot_df['fdr'] - target)).argmin()
        highlights.append(plot_df.iloc[idx])
    
    highlights_df = pd.DataFrame(highlights)
    
    ax.scatter(
        highlights_df['fdr'], highlights_df['ratio'], 
        s=70, facecolors='white', edgecolors=VIBRANT_ORANGE, 
        linewidth=1.5, zorder=3
    )

    ax.set_xlabel("FDR Threshold (q-value)")
    ax.set_ylabel("Ratio (Mapped / Unmapped)")
    
    ax.set_xlim(0, 1.0)
    ax.set_xticks(np.arange(0, 1.1, 0.1))
    
    ax.grid(True, which='major', axis='both', color='#f0f0f0', linestyle='-', linewidth=0.5)
    
    sns.despine()

    if folder:
        Path(folder).mkdir(parents=True, exist_ok=True)
        out_name = f"fig2c_{run}_ratio_full_range.svg"
        plt.savefig(os.path.join(folder, out_name), format="svg", bbox_inches='tight')
    
    plt.show()



def plot_coverage_vs_fdr_curve(run, df, reference_protein, folder, min_identity=0.8, chain="Light", ratio=1):
    """
    Plot the PSM coverage across FDR thresholds
    
    :param run: sample run name
    :param df: filtered data
    :param reference_protein: reference protein normalised
    :param folder: figures folder
    :param min_identity: minimum identity required for mapping to the protein
    :param chain: specifically for antibodies
    :param ratio: set pubblication style
    """
    set_publication_style()
    
    df_copy = df.copy()
    fdr_thresholds = np.linspace(0, 1.0, 50)
    coverage_data = []

    temp_stats_folder = os.path.join(folder, "temp_stats_calc")
    os.makedirs(temp_stats_folder, exist_ok=True)

    for fdr in tqdm(fdr_thresholds):
        subset = df_copy[df_copy["psm_q_value"] <= fdr]
        if subset.empty:
            coverage_data.append({"fdr": fdr, "coverage": 0.0})
            continue

        sequences = subset["cleaned_preds"].tolist()
        try:
            mapped_psms = process_protein_contigs_scaffold(
                sequences, reference_protein, max_mismatches=10, min_identity=min_identity
            )
            df_mapped = create_dataframe_from_mapped_sequences(data=mapped_psms)
            stats = helpers.compute_assembly_statistics(
                df_mapped, f"temp_{chain}_{fdr:.2f}", temp_stats_folder, reference_protein
            )
            
            cov_value = 0
            for k in ["Coverage", "Protein Coverage", "coverage", "protein_coverage"]:
                if k in stats:
                    val = stats[k]
                    cov_value = float(val.replace("%", "")) if isinstance(val, str) else float(val)
                    if cov_value <= 1.0 and cov_value > 0:
                        cov_value *= 100
                    break
            coverage_data.append({"fdr": fdr, "coverage": cov_value})
        except:
            coverage_data.append({"fdr": fdr, "coverage": 0.0})

    plot_df = pd.DataFrame(coverage_data).sort_values("fdr")

    _, ax = plt.subplots(figsize=get_figsize(width_ratio=ratio))
    VIBRANT_BLUE = "#4A90E2"

    ax.plot(
        plot_df["fdr"], 
        plot_df["coverage"], 
        color=VIBRANT_BLUE, 
        linestyle='-', 
        linewidth=2, 
        zorder=2
    )

    critical_fdr = [0.01, 0.05, 0.1, 0.2, 0.5, 1.0]
    for val in critical_fdr:
        idx = (plot_df['fdr'] - val).abs().idxmin()
        ax.scatter(
            plot_df.loc[idx, 'fdr'], 
            plot_df.loc[idx, 'coverage'], 
            color='white', 
            edgecolor=VIBRANT_BLUE, 
            s=40, 
            zorder=3, 
            linewidth=1.2
        )

    ax.set_xlabel("FDR Threshold (q-value)", fontweight='normal')
    ax.set_ylabel("Protein Coverage (%)", fontweight='normal')
    ax.set_xlim(0, 1.0)
    ax.set_ylim(0, 105)
    ax.set_xticks(np.arange(0, 1.1, 0.2))
    
    ax.grid(True, which='major', axis='both', color='#f0f0f0', linestyle='-', linewidth=0.5)
    sns.despine()

    if folder:
        Path(folder).mkdir(parents=True, exist_ok=True)
        plt.savefig(os.path.join(folder, f"fig2d_{run}_{chain}_coverage_solid.svg"), format="svg", bbox_inches='tight')
    
    plt.show()



def plot_protease_confidence_ridges(df, colors_path='json/protease_colors.json', folder=None):
    """
    Docstring for plot_protease_confidence_ridges
    
    :param df: Description
    :param colors_path: Description
    :param folder: Description
    """
    set_publication_style()
    
    with open(colors_path, 'r') as f:
        protease_colors = json.load(f)

    proteases = sorted(df['protease'].unique())
    n_prot = len(proteases)

    total_w, total_h = get_figsize(width_ratio=1)
    row_height = total_h / n_prot

    sns.set(style="white", rc={"axes.facecolor": (0, 0, 0, 0)})
    
    g = sns.FacetGrid(
        df, row="protease", hue="protease", 
        aspect=1, 
        height=row_height, 
        palette=protease_colors,
        row_order=proteases
    )

    g.map(sns.kdeplot, "conf", log_scale=True, bw_adjust=.7, fill=True, alpha=0.8, linewidth=0)
    
    g.map(plt.axhline, y=0, lw=1, clip_on=False, color='#333333', alpha=0.8)

    def label(x, color, label):
        ax = plt.gca()
        ax.text(1.02, .1, label, fontweight="normal", color="black",
                ha="left", va="center", transform=ax.transAxes, fontsize=12)

    g.map(label, "conf")

    g.set_titles("")
    g.set(yticks=[], ylabel="")
    g.despine(bottom=True, left=True)
    
    g.figure.subplots_adjust(hspace=0.4) 
    g.figure.set_size_inches(total_w, total_h)

    tick_vals = [1e-10, 1e-8, 1e-6, 1e-4, 1e-2, 1] 
    tick_labels = ["0.0", "0.2", "0.4", "0.6", "0.8", "1.0"]

    for ax in g.axes.flat:
        ax.set_xlim(1e-10, 1)
        ax.set_xticks(tick_vals)
        ax.set_xticklabels(tick_labels)
        for label in ax.get_xticklabels():
            label.set_fontweight('normal')

    plt.xlabel("Confidence", fontsize=15, fontweight='normal', labelpad=20)
    
    if folder:
        Path(folder).mkdir(parents=True, exist_ok=True)
        g.savefig(f"{folder}/supp_fig1a_ridges_proteases.svg", format="svg", bbox_inches='tight')
    
    plt.show()



def plot_sunburst(df, output_folder, output_file, json_colors_path):
    """
    Docstring for plot_sunburst
    
    :param df: Description
    :param output_folder: Description
    :param output_file: Description
    :param json_colors_path: Description
    """
    set_publication_style()
    
    protease_colors_map = {}
    if os.path.exists(json_colors_path):
        with open(json_colors_path, 'r') as f:
            protease_colors_map = json.load(f)
    
    fallback_palette = ['#8dd3c7','#ffffb3','#bebada','#fb8072','#80b1d3','#fdb462','#b3de69','#fccde5']

    df = df.copy()
    df['is_mapped'] = df['mapped'].astype(str).str.lower().isin(['true', 'mapped', 'yes', '1'])

    counts = df.groupby(['protease', 'is_mapped']).size().unstack(fill_value=0).reset_index()
    counts.columns = [str(c) for c in counts.columns]
    
    if 'True' not in counts.columns: counts['True'] = 0
    if 'False' not in counts.columns: counts['False'] = 0
    counts.rename(columns={'True': 'mapped_count', 'False': 'unmapped_count'}, inplace=True)
    
    counts['total'] = counts['mapped_count'] + counts['unmapped_count']
    total_dataset_psms = counts['total'].sum()
    counts = counts.sort_values('total', ascending=False).reset_index(drop=True)

    colors_list = []
    for idx, row in counts.iterrows():
        colors_list.append(protease_colors_map.get(row['protease'], fallback_palette[idx % len(fallback_palette)]))
    counts['color'] = colors_list

    inner_sizes = counts['total'].values
    inner_colors = counts['color'].values
    inner_labels = [f"{row['total']/total_dataset_psms*100:.1f}%" for _, row in counts.iterrows()]
    
    outer_sizes, outer_colors, outer_labels = [], [], []
    for _, row in counts.iterrows():
        mapped, unmapped, total, c = row['mapped_count'], row['unmapped_count'], row['total'], row['color']
        outer_sizes.extend([mapped, unmapped])
        outer_colors.extend([c, (0,0,0,0)])
        eff_pct = (mapped / total * 100) if total > 0 else 0
        outer_labels.extend([f"{eff_pct:.0f}%" if (mapped/total_dataset_psms > 0.01) else "", ""])

    _, ax = plt.subplots(figsize=get_figsize(width_ratio=2))
    
    wedges_inner, texts_inner = ax.pie(
        inner_sizes, radius=0.7, colors=inner_colors, labels=inner_labels,
        labeldistance=0.5, startangle=90, counterclock=False, 
        wedgeprops=dict(width=0.4, edgecolor='white', linewidth=1.5) 
    )
    
    for t in texts_inner:
        t.set_size(9)
        t.set_fontweight("normal")

    wedges_outer, texts_outer = ax.pie(
        outer_sizes, radius=1.0, colors=outer_colors, labels=outer_labels,
        labeldistance=0.85, startangle=90, counterclock=False,
        wedgeprops=dict(width=0.25, edgecolor='white', linewidth=1)
    )

    for i, wedge in enumerate(wedges_outer):
        if i % 2 != 0: wedge.set_edgecolor('none')

    for t in texts_outer:
        t.set_size(9)
        t.set_fontweight("normal")

    ax.legend(
        wedges_inner, counts['protease'],
        title="Proteases",
        loc="center left",
        bbox_to_anchor=(1, 0, 0.5, 1),
        frameon=False
    )

    ax.add_artist(plt.Circle((0,0), 0.3, fc='white'))
    
    if output_folder and output_file:
        Path(output_folder).mkdir(parents=True, exist_ok=True)
        plt.savefig(os.path.join(output_folder, output_file), format='svg', bbox_inches='tight')
    
    plt.show()


def plot_barplot_composite_coverage_scores(
    csv_path, 
    category, 
    config_json_path="json/colors.json", 
    output_file_coverage=None, 
    output_file_composite=None,
    width_ratio=3
):

    set_publication_style()
    df = pd.read_csv(csv_path)
    
    try:
        with open(config_json_path, 'r') as f:
            color_data = json.load(f)
        palette = [
            color_data.get(category.lower(), {}).get("contig", "#a6cee3"),
            color_data.get(category.lower(), {}).get("scaffold", "#1f78b4")
        ]
    except Exception:
        palette = ["#a6cee3", "#1f78b4"]

    df['Type'] = df['assembly_method'].apply(lambda x: "Contigs" if "Contigs" in x else "Scaffolds")
    
    try:
        df['sample_num'] = df['sample'].str.extract(r'(\d+)').astype(int)
        df = df.sort_values('sample_num')
    except Exception:
        pass 

    plots_to_generate = [
        ("coverage", "Coverage", output_file_coverage),
        ("composite_score", "Composite score", output_file_composite)
    ]

    for col_name, y_label, out_file in plots_to_generate:
        fig_width, fig_height = get_figsize(width_ratio=width_ratio)
        plt.figure(figsize=(fig_width, fig_height))

        ax = sns.barplot(
            data=df,
            x='sample',
            y=col_name,
            hue='Type',
            palette=palette,
            edgecolor='black',
            linewidth=0.8
        )

        ax.set_ylim(0, 1.05)
        ax.set_xlabel("Samples")
        ax.set_ylabel(y_label)
        
        plt.legend(title="", loc='upper right', frameon=False)
        
        sns.despine()
        plt.tight_layout()

        if out_file:
            os.makedirs(os.path.dirname(out_file), exist_ok=True)
            plt.savefig(out_file, bbox_inches='tight')
            print(f"Saved {y_label} plot to {out_file}")
        
        plt.show()
        plt.close()



def fdr_ratio_mapped_unmapped(run, df, folder):
    bin_centers = []
    ratios = []

    for start in np.arange(0, 1, 0.05):
        end = start + 0.05
        subset_map = df[(df["mapped"] == "True") & (df["conf"] >= start) & (df["conf"] < end)]
        subset_unmap = df[(df["mapped"] == "False") & (df["conf"] >= start) & (df["conf"] < end)]

        count_map = len(subset_map)
        count_unmap = len(subset_unmap)

        ratio = count_map / count_unmap if count_unmap > 0 else np.nan

        bin_center = start + 0.025
        bin_centers.append(bin_center)
        ratios.append(ratio)

    bin_centers = np.array(bin_centers)
    ratios = np.array(ratios)

    y_horizontal = 1.3

    intersect_x = None
    for i in range(len(bin_centers) - 1):
        if (ratios[i] <= y_horizontal and ratios[i + 1] >= y_horizontal) or (
            ratios[i] >= y_horizontal and ratios[i + 1] <= y_horizontal
        ):
            x1, x2 = bin_centers[i], bin_centers[i + 1]
            y1, y2 = ratios[i], ratios[i + 1]
            intersect_x = x1 + (y_horizontal - y1) * (x2 - x1) / (y2 - y1)
            break

    fig = go.Figure()

    fig.add_trace(
        go.Scatter(
            x=bin_centers,
            y=ratios,
            mode="lines+markers",
            line=dict(color="#1f77b4", width=3.5),
            name="mapped/unmapped ratio",
        )
    )

    if intersect_x is not None:
        fig.add_shape(
            type="line",
            x0=-0.025,
            x1=intersect_x,
            xref="x",
            y0=y_horizontal,
            y1=y_horizontal,
            yref="y",
            line=dict(color="black", width=1.5, dash="dash"),
        )

    if intersect_x is not None:
        fig.add_shape(
            type="line",
            x0=intersect_x,
            x1=intersect_x,
            y0=0.001,
            y1=y_horizontal,
            yref="y",
            line=dict(color="black", width=1.5, dash="dash"),
        )

        fig.add_annotation(
            x=intersect_x - 0,
            y=y_horizontal - 0.8,
            text=f"Confidence: {intersect_x:.2f}",
            showarrow=False,
            font=dict(size=16, color="black"),
        )

    fig.update_layout(
        xaxis_title="Confidence",
        yaxis_title="Ratio mapped/unmapped",
        template="plotly_white",
        height=600,
        width=800,
        font=dict(size=22, family="Arial, sans-serif", color="black"),
        xaxis=dict(
            title=dict(font=dict(size=20)),
            tickmode="array",
            tickvals=np.arange(0, 1.1, 0.1),
            ticktext=[f"{x:.1f}" for x in np.arange(0, 1.1, 0.1)],
            zeroline=False,
            linewidth=1,
            linecolor="black",
            showline=True,
            showgrid=False,
        ),
        yaxis=dict(
            title=dict(font=dict(size=20)),
            type="log",
            range=[-3, 1.5],
            tickvals=[10**i for i in range(-3, 2)],
            ticktext=["10⁻³", "10⁻²", "10⁻¹", "10⁰", "10¹"],
            showline=True,
            linewidth=1,
            linecolor="black",
            zeroline=False,
            showgrid=False,
            side="left",
        ),
    )

    fig.write_image(f"{folder}/{run}_fdr_ratio_mapped_unmapped.svg")


def plot_relative_map_distribution(run, df, reference, folder, title=False):
    df = df[df["conf"] >= 0].copy()
    df["mapped"] = df["cleaned_preds"].apply(lambda x: "mapped" if x in reference else "unmapped")

    bins = np.arange(0, 1, 0.05)
    # bins = np.arange(0, 1.02, 0.02)
    df["bin"] = pd.cut(df["conf"], bins=bins, labels=bins[:-1])

    bin_counts = df.groupby("bin")["mapped"].count()
    mapped_counts = df[df["mapped"] == "mapped"].groupby("bin")["mapped"].count()
    unmapped_counts = df[df["mapped"] == "unmapped"].groupby("bin")["mapped"].count()

    mapped_percentages = (mapped_counts / bin_counts) * 100
    unmapped_percentages = (unmapped_counts / bin_counts) * 100

    hist_df = pd.DataFrame(
        {
            "confidence": bins[:-1],
            "Mapped": mapped_percentages.fillna(0).values,
            "Unmapped": unmapped_percentages.fillna(0).values,
        }
    )
    fig = px.line(
        hist_df,
        x="confidence",
        y=["Mapped", "Unmapped"],
        markers=False,
        line_shape="linear",
        color_discrete_map={"Mapped": "#A5C8E1", "Unmapped": "#FFC48C"},
    )
    # orange - unmapped: #ff7f0e, blue: #1f77b4, brown #AF6E7E

    title_text = "Relative distribution of mapped and unmapped peptides by confidence" if title else ""

    fig.update_layout(
        title=title_text,
        xaxis_title="Confidence",
        yaxis_title="Percentage (%)",
        template="plotly_white",
        height=600,
        width=800,
        font=dict(family="Arial, sans-serif", color="black"),
        showlegend=True,
        legend_title="",
        legend=dict(font=dict(size=16)),
    )

    fig.update_yaxes(
        range=[0, 100],
        title=dict(font=dict(size=20)),
        showline=True,
        linecolor="black",
        linewidth=1,
        showgrid=False,
    )

    fig.update_xaxes(
        range=[0, 1],
        title=dict(font=dict(size=20)),
        showline=True,
        linecolor="black",
        linewidth=1,
        showgrid=False,
        tickmode="linear",
        dtick=0.1,
        tickangle=0,
    )

    fig.add_vline(
        x=0.88,
        line_width=2,
        line_dash="dash",
        line_color="black",
        annotation_text="Cutoff",
        annotation_position="top",
        annotation_font_size=16,
    )

    for line_name in ["Mapped", "Unmapped"]:
        base_color = (165, 200, 225) if line_name == "Mapped" else (255, 196, 140)

        x_vals = hist_df["confidence"].astype(float).values
        y_vals = hist_df[line_name].astype(float).values
        fine_x = np.linspace(0.88, 1.0, 50)  # Start from 0.88 instead of intersection_x
        fine_y = np.interp(fine_x, x_vals, y_vals)

        for i in range(len(fine_x) - 1):
            x0 = fine_x[i]
            x1 = fine_x[i + 1]
            y_segment = (fine_y[i] + fine_y[i + 1]) / 2.0
            x_mid = (x0 + x1) / 2.0
            alpha = 1 - (x_mid - 0.88) / (1 - 0.88)
            fillcolor = f"rgba({base_color[0]}, {base_color[1]}, {base_color[2]}, {alpha:.2f})"

            fig.add_shape(
                type="rect",
                xref="x",
                yref="y",
                x0=x0,
                x1=x1,
                y0=0,
                y1=y_segment,
                fillcolor=fillcolor,
                line=dict(width=0),
                layer="below",
            )

    fig.update_traces(line=dict(width=3.5))
    fig.for_each_trace(lambda t: t.update(name="mapped") if t.name == "Mapped" else None)

    fig.for_each_trace(lambda t: t.update(name="unmapped") if t.name == "Unmapped" else None)

    fig.write_image(f"{folder}/{run}_relative_mapped_unmapped_distribution.svg")


def plot_map_distribution(run, df, reference, folder, threshold, title=False):
    df = df[df["conf"] >= threshold].copy()

    df["mapped"] = df["cleaned_preds"].apply(lambda x: "mapped" if x in reference else "unmapped")

    bins = np.arange(threshold, 1.002, 0.02)

    counts_mapped, edges = np.histogram(df[df["mapped"] == "mapped"]["conf"], bins=bins)
    counts_unmapped, _ = np.histogram(df[df["mapped"] == "unmapped"]["conf"], bins=bins)

    bin_centers = edges[:-1] + (0.5 * (edges[1] - edges[0]))

    hist_df = pd.DataFrame(
        {
            "confidence": np.tile(bin_centers, 2),
            "count": np.concatenate([counts_mapped, counts_unmapped]),
            "category": ["mapped"] * len(counts_mapped) + ["unmapped"] * len(counts_unmapped),
        }
    )

    fig = px.bar(
        hist_df,
        x="confidence",
        y="count",
        color="category",
        color_discrete_map={"mapped": "#A5C8E1", "unmapped": "#FFC48C"},
        barmode="stack",
    )

    title_text = "Distribution of mapped and unmapped sequences by confidence" if title else ""

    fig.update_layout(
        title=title_text,
        xaxis_title="Confidence",
        yaxis_title="PSMs counts",
        legend_title="",
        legend_font=dict(size=16),
        template="plotly_white",
        showlegend=True,
        height=600,
        width=800,
        font=dict(size=22, family="Arial, sans-serif", color="black"),
    )

    fig.update_yaxes(
        title=dict(font=dict(size=20)),
        showline=True,
        linecolor="black",
        linewidth=1,
        showgrid=False,
    )

    fig.update_xaxes(
        title=dict(font=dict(size=20)),
        showline=True,
        linecolor="black",
        linewidth=1,
        showgrid=False,
        tickmode="linear",
        dtick=0.02,
    )

    fig.for_each_trace(lambda t: t.update(name="mapped") if t.name == "Mapped" else None)

    fig.for_each_trace(lambda t: t.update(name="unmapped") if t.name == "Unmapped" else None)

    fig.write_image(f"{folder}/{run}_psms_mapped_unmapped_distribution.svg")


def plot_confidence_distribution(df, folder_figures, min_conf=0, max_conf=1):
    """Plots the distribution of confidence scores from a DataFrame."""
    # Filter the data based on the specified range
    filtered_df = df[(df["conf"] >= min_conf) & (df["conf"] <= max_conf)]

    fig = go.Figure()

    fig.add_trace(
        go.Histogram(
            x=filtered_df["conf"],
            xbins=dict(start=min_conf, end=max_conf, size=(max_conf - min_conf) / 40),
            marker=dict(color="brown"),
            opacity=1,  # remove opacity
        )
    )

    fig.update_layout(
        title="Confidence score distribution between {} and {}".format(min_conf, max_conf),
        xaxis_title="Values",
        yaxis_title="Frequency",
        bargap=0.1,
        height=700,
        width=1200,
        margin=dict(l=50, r=50, t=100, b=100),
    )

    fig.update_xaxes(
        showgrid=True,
        gridcolor="lightgray",
        ticklabelposition="outside bottom",
        dtick=0.02,
    )
    fig.update_yaxes(showgrid=True, gridcolor="lightgray")
    fig.write_image(f"{folder_figures}/confidence_distribution_range_{min_conf}_{max_conf}.png")


def plot_protease_distribution(protease_counts, folder_figures):
    """Creates an interactive bar plot of protease distribution using Plotly.

    Parameters:
        protease_counts (pandas.Series): A Pandas Series with protease names as the index
                                         and their counts as the values.
    """
    # Convert the Series to a DataFrame for compatibility with Plotly
    protease_df = protease_counts.reset_index()
    protease_df.columns = ["Protease", "Count"]

    fig = px.bar(
        protease_df,
        x="Protease",
        y="Count",
        title="Proteases distribution",
        labels={"Protease": "Proteases", "Count": "Counts"},
        text="Count",
    )

    fig.update_traces(textposition="outside", width=0.4)

    mm_to_px = 3.78
    width_mm = 240
    height_mm = 200

    fig.update_layout(
        width=int(width_mm * mm_to_px),
        height=int(height_mm * mm_to_px),
        xaxis_title="Proteases",
        yaxis_title="Counts",
        xaxis_tickangle=0,
        showlegend=False,
        title_font_size=12,
        font=dict(
            family="Arial, sans-serif",
            color="black",
        ),
        margin=dict(t=50, b=50, l=50, r=100),
        plot_bgcolor="white",
        paper_bgcolor="white",
    )

    fig.write_image(f"{folder_figures}/proteases_distribution.svg")


def map_to_protein(seq, protein, max_mismatches, min_identity):
    """Maps a sequence (`seq`) to a target protein sequence, allowing for mismatches,
    and identifies the best match based on the maximum mismatches and minimum identity threshold.
    """

    best_match = None
    best_identity = 0

    # Slide `seq` across `protein` to check each possible alignment position
    for i in range(len(protein) - len(seq) + 1):
        mismatches_count = 0
        mismatch_positions = []

        # Compare `seq` to the substring of `protein` at the current position
        for j in range(len(seq)):
            if seq[j] != protein[i + j]:
                mismatches_count += 1
                mismatch_positions.append(j)

            # Stop checking this alignment if mismatches exceed the allowed threshold
            if mismatches_count > max_mismatches:
                break

        # If this alignment meets the mismatch requirement, calculate identity
        if mismatches_count <= max_mismatches:
            if len(seq) == 0:
                print("Zero length sequence found.")
                print(seq)
                continue

            identity = 1 - mismatches_count / len(seq)

            # Update the best match if this alignment has a higher identity and meets the minimum requirement
            if identity >= min_identity and identity > best_identity:
                best_match = (i, i + len(seq), mismatch_positions, identity)
                best_identity = identity

    return best_match


def process_protein_contigs_scaffold(assembled_contigs, target_protein, max_mismatches, min_identity):
    """Maps each contig in `assembled_contigs` to a target protein sequence (`target_protein`)
    and identifies which contigs match based on specified mismatch and identity thresholds.
    """
    mapped_sequences = []

    # Map each contig to the target protein
    for contig in assembled_contigs:
        # Attempt to map the contig to the target protein
        target_mapping = map_to_protein(
            contig,
            target_protein,
            max_mismatches=max_mismatches,
            min_identity=min_identity,
        )

        if target_mapping:
            mapped_sequences.append((contig, target_mapping))

    return mapped_sequences


def write_mapped_contigs(mapped_contigs, folder, filename_prefix):
    """Writes mapped contigs to a FASTA file with detailed annotations for each contig."""

    records = []
    for idx, (contig, mapping) in enumerate(mapped_contigs):
        start, end, mismatches, identity = mapping
        record = Bio.SeqRecord.SeqRecord(
            Bio.Seq.Seq(contig),
            id=f"Contig {idx + 1}",
            description=f"length: {len(contig)}, start: {start}, end: {end}, mismatches: {len(mismatches)}, identity: {identity:.2f}",
        )
        records.append(record)
    Bio.SeqIO.write(records, os.path.join(folder, f"{filename_prefix}.fasta"), "fasta")


def plot_contigs(mapped_contigs, prot_seq, title, output_file):
    sns.set("paper", "ticks", "colorblind", font_scale=1.5)
    _, ax = plt.subplots(figsize=(12, 4))

    ax.add_patch(patches.Rectangle((0, 0), len(prot_seq), 0.2, facecolor="#e6f0ef", edgecolor="#e6f0ef"))

    tracks = {}
    ind = 0

    for _, (contig, mapping) in tqdm(enumerate(mapped_contigs), desc="Plotting contigs"):
        start_index, end_index, mismatches, _ = mapping

        ind += 1
        placed = False
        for track_num, track in tracks.items():
            if not any(s <= end_index <= e or s <= start_index <= e for s, e in track):
                track.append((start_index, end_index))
                ax.add_patch(
                    patches.Rectangle(
                        (start_index, 0.3 + 0.1 * track_num),
                        len(contig),
                        0.075,
                        facecolor="#007EA7",
                        edgecolor="#007EA7",
                        label="Contig" if not placed else "",
                    )
                )
                placed = True
                break

        if not placed:
            track_num = len(tracks) + 1
            tracks[track_num] = [(start_index, end_index)]
            ax.add_patch(
                patches.Rectangle(
                    (start_index, 0.3 + 0.1 * track_num),
                    len(contig),
                    0.075,
                    facecolor="#007EA7",
                    edgecolor="#007EA7",
                    label="Contig" if not placed else "",
                )
            )

        for mismatch in mismatches:
            ax.add_patch(
                patches.Rectangle(
                    (start_index + mismatch, 0.3 + 0.1 * track_num),
                    1,
                    0.075,
                    facecolor="#FCAB64",
                    edgecolor="#FCAB64",
                )
            )

    print(f"Plotted {ind} contigs.")

    ax.set_xlim(0, len(prot_seq))
    ax.set_ylim(0, 0.3 + 0.1 * (len(tracks) + 1))
    ax.get_yaxis().set_visible(False)

    ax.set_xlabel("Sequence range")
    ax.set_title(title)

    handles, labels = ax.get_legend_handles_labels()
    by_label = dict(zip(labels, handles, strict=False))
    ax.legend(
        by_label.values(),
        by_label.keys(),
        loc="center right",
        frameon=False,
        bbox_to_anchor=(1.2, 0.8),
    )

    for spine in ax.spines.values():
        spine.set_visible(False)

    plt.tight_layout(pad=1)
    plt.savefig(output_file, dpi=300, bbox_inches="tight")
    plt.close()
      


def create_dataframe_from_mapped_sequences(data):
    """Takes a list of tuples containing sequence data and returns a structured DataFrame."""
    # Create the initial DataFrame
    df = pd.DataFrame(data, columns=["Sequence", "Details"])

    # Expand the 'Details' column into separate columns
    df[["start", "end", "mismatches_pos", "identity_score"]] = pd.DataFrame(df["Details"].tolist(), index=df.index)

    df.drop(columns=["Details"], inplace=True)
    df.rename(columns={"Sequence": "sequence"}, inplace=True)

    return df


def mapping_sequences(
    mapped_sequences,
    prot_seq,
    category,
    config_json_path="json/colors.json",
    output_folder=".",
    output_file=None,
    show_figure=False,
):
    set_publication_style()
    
    try:
        with open(config_json_path, 'r') as f:
            color_data = json.load(f)
        main_color = color_data.get(category, {}).get("scaffold", "#1f78b4")
    except Exception as e:
        main_color = "#1f78b4"

    fig_width, fig_height = get_figsize(width_ratio=3)
    
    common_height = 0.3
    track_spacing = 0.45
    base_y_offset = 0.6

    _, ax = plt.subplots(figsize=(fig_width, fig_height))

    ax.add_patch(patches.Rectangle(
        (0, 0), len(prot_seq), common_height,
        linewidth=0, facecolor='#e6f0ef', zorder=0
    ))

    tracks = {}
    colors = {
        "match": main_color,
        "mismatch": "#b30000",
        "D_to_N": "#000000",
        "E_to_Q": "#A8A29E",
    }

    for seq, mapping in tqdm(mapped_sequences, desc=f"Mapping {category}"):
        start_index, end_index, mismatches, _ = mapping
        
        placed = False
        for track_num in sorted(tracks.keys()):
            if not any(max(s, start_index) < min(e, end_index) for s, e in tracks[track_num]):
                tracks[track_num].append((start_index, end_index))
                current_track_num = track_num
                placed = True
                break
        
        if not placed:
            current_track_num = len(tracks)
            tracks[current_track_num] = [(start_index, end_index)]

        current_y = base_y_offset + (current_track_num * track_spacing)
        
        ax.add_patch(patches.Rectangle(
            (start_index, current_y), end_index - start_index, common_height,
            linewidth=0.8, edgecolor='white', facecolor=colors["match"], alpha=0.9
        ))
        # showing mismatches
        for mismatch in mismatches:
            abs_index = start_index + mismatch
            if abs_index >= len(prot_seq) or mismatch >= len(seq): continue

            ref_aa = prot_seq[abs_index]
            query_aa = seq[mismatch]
            # showing modifications
            if query_aa == "D" and ref_aa == "N":
                mut_color = colors["D_to_N"]
            elif query_aa == "E" and ref_aa == "Q":
                mut_color = colors["E_to_Q"]
            else:
                mut_color = colors["mismatch"]

            ax.add_patch(patches.Rectangle(
                (abs_index, current_y), 1, common_height,
                linewidth=0, facecolor=mut_color, zorder=10 
            ))

    ax.set_xlim(0, len(prot_seq))
    max_y = base_y_offset + (len(tracks) * track_spacing) + 0.2
    ax.set_ylim(-0.1, max_y)
    
    ax.set_xlabel("Residue position")
    ax.set_yticks([])
    sns.despine(left=True)

    legend_labels = [
        patches.Patch(color=colors["match"], label=f"Match"),
        patches.Patch(color=colors["mismatch"], label="Mismatch"),
        patches.Patch(color=colors["D_to_N"], label="D \u2192 N"),
        patches.Patch(color=colors["E_to_Q"], label="E \u2192 Q"),
    ]
    ax.legend(handles=legend_labels, loc='upper center', 
              bbox_to_anchor=(0.5, 1.25), ncol=4)

    plt.tight_layout()

    if output_file:
        try:
            os.makedirs(output_folder, exist_ok=True)
            save_path = os.path.join(output_folder, output_file)
            plt.savefig(save_path, bbox_inches='tight')
        except PermissionError:
            plt.savefig(output_file, bbox_inches='tight')

    if show_figure:
        plt.show()
    plt.close()


def mapping_psms_protease_associated(
    mapped_sequences,
    prot_seq,
    labels,
    output_folder,
    output_file,
    json_colors_path=None,
    show_figure=False,
):
    """
    Docstring for mapping_psms_protease_associated_seaborn
    
    :param mapped_sequences: Description
    :param prot_seq: Description
    :param labels: Description
    :param output_folder: Description
    :param output_file: Description
    :param json_colors_path: Description
    :param show_figure: Description
    """
    
    set_publication_style()
    
    protease_colors_map = {}
    if json_colors_path and os.path.exists(json_colors_path):
        with open(json_colors_path, 'r') as f:
            protease_colors_map = json.load(f)
    
    unique_labels = sorted(list(set(labels)))
    fallback_palette = sns.color_palette("husl", len(unique_labels)).as_hex()
    label_color = {
        lab: protease_colors_map.get(lab, fallback_palette[i % len(fallback_palette)]) 
        for i, lab in enumerate(unique_labels)
    }
    
    fig_w, fig_h = get_figsize(width_ratio=3)
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))

    ax.add_patch(patches.Rectangle(
        (0, 0), len(prot_seq), 0.3, 
        linewidth=0, 
        facecolor='#e6f0ef', 
        zorder=0
    ))

    tracks = {}
    bar_height = 0.6
    track_spacing = 1.0
    base_y_offset = 1.2

    for idx, (_, mapping) in tqdm(enumerate(mapped_sequences), desc="Plotting peptides", total=len(mapped_sequences)):
        start_index, end_index, _, _ = mapping
        lab = labels[idx]
        placed = False
        
        for track_num in sorted(tracks.keys()):
            if not any(max(s, start_index) < min(e, end_index) for s, e in tracks[track_num]):
                tracks[track_num].append((start_index, end_index))
                current_y = base_y_offset + (track_num * track_spacing)
                placed = True
                break
        
        if not placed:
            new_track_num = len(tracks)
            tracks[new_track_num] = [(start_index, end_index)]
            current_y = base_y_offset + (new_track_num * track_spacing)

        ax.add_patch(patches.Rectangle(
            (start_index, current_y), end_index - start_index, bar_height,
            linewidth=0, facecolor=label_color[lab], alpha=0.85
        ))

    max_track = len(tracks) if tracks else 0
    ax.set_xlim(0, len(prot_seq))
    ax.set_ylim(-0.2, base_y_offset + (max_track * track_spacing) + 0.5)
    
    ax.set_xlabel("Residue position")
    ax.set_yticks([])
    
    sns.despine(ax=ax, left=True, top=True, right=True, bottom=False)

    legend_patches = [patches.Patch(color=label_color[lab], label=lab) for lab in unique_labels]
    
    ax.legend(
        handles=legend_patches, 
        title="Proteases",
        loc='upper center', 
        bbox_to_anchor=(0.5, 1.25),
        ncol=min(len(unique_labels), 5),
        frameon=False
    )

    plt.tight_layout()

    if output_folder and output_file:
        os.makedirs(output_folder, exist_ok=True)
        save_path = os.path.join(output_folder, output_file)
        plt.savefig(save_path, format='svg', bbox_inches='tight', transparent=True)

    if show_figure:
        plt.show()
    else:
        plt.close()