#!/usr/bin/env python

r""" Visualization module for InstaNexus.

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
__date__ = 14 Nov 2025
__maintainer__ = Marco Reverenna
__email__ = marcor@dtu.dk
__status__ = Dev
"""

import os
import numpy as np
import pandas as pd
import Bio.SeqIO
import matplotlib.patches as patches
import matplotlib.pyplot as plt
import plotly.express as px
import plotly.graph_objects as go
import seaborn as sns
from tqdm import tqdm


def missing_values_barplot(run, dataframe, folder):

    dataframe["missing_preds"] = dataframe["preds"].isna()

    missing_counts_df = dataframe["missing_preds"].value_counts().reset_index()
    missing_counts_df.columns = ["PSMs", "Count"]
    missing_counts_df["PSMs"] = missing_counts_df["PSMs"].map(
        {True: "Missing", False: "Valid"}
    )

    missing_counts_df = missing_counts_df.sort_values(
        "PSMs", key=lambda x: x.map({"Valid": 0, "Missing": 1})
    )

    valid_count = (
        missing_counts_df.loc[missing_counts_df["PSMs"] == "Valid", "Count"].values[0]
        if "Valid" in missing_counts_df["PSMs"].values
        else 0
    )
    missing_count = (
        missing_counts_df.loc[missing_counts_df["PSMs"] == "Missing", "Count"].values[0]
        if "Missing" in missing_counts_df["PSMs"].values
        else 0
    )

    total = valid_count + missing_count
    valid_pct = (valid_count / total * 100) if total > 0 else 0
    missing_pct = (missing_count / total * 100) if total > 0 else 0

    print(f"Valid PSMs: {valid_count} ({valid_pct:.2f}%)")
    print(f"Missing PSMs: {missing_count} ({missing_pct:.2f}%)")
    print(f"Total PSMs: {total}")

    color_map = {"Valid": "#337AB7", "Missing": "#FF5733"}

    fig = go.Figure()

    fig.add_trace(
        go.Bar(
            x=["Valid"],
            y=[valid_count],
            marker_color=color_map["Valid"],
            width=[0.4],
            name="Valid",
        )
    )

    fig.add_trace(
        go.Bar(
            x=["Missing"],
            y=[missing_count],
            marker_color=color_map["Missing"],
            width=[0.4],
            name="Missing",
        )
    )

    fig.update_layout(
        title="",
        xaxis_title="PSMs",
        yaxis_title="Count",
        barmode="group",
        width=800,
        height=600,
        margin=dict(t=50, l=10, r=10, b=10),
        legend_title_text="",
        paper_bgcolor="white",
        plot_bgcolor="white",
        font=dict(size=22, family="Arial, sans-serif", color="black"),
    )

    fig.update_xaxes(
        showline=True,
        linewidth=1,
        linecolor="black",
        title_font=dict(size=20),
        tickfont=dict(size=18),
    )

    fig.update_yaxes(
        showline=True,
        linewidth=1,
        linecolor="black",
        title_font=dict(size=20),
        tickfont=dict(size=18),
    )

    fig.write_image(f"{folder}/{run}_missing_value_bar.svg")


def plot_map_unmap_distribution(df, reference, run, folder, conf_lim, title=False):

    df = df[df["conf"] >= conf_lim]

    df["mapped"] = df["cleaned_preds"].apply(
        lambda x: "mapped" if x in reference else "unmapped"
    )

    fig = px.histogram(
        df,
        x="conf",
        color="mapped",
        nbins=int((1 - conf_lim) * 50),
        barmode="overlay",
        color_discrete_map={"mapped": "#1f77b4", "unmapped": "#ff7f0e"},
        # orange: #ff7f0e, blue: #1f77b4, brown #AF6E7E
    )

    fig.update_traces(xbins=dict(start=0, end=1, size=0.01))

    title_text = (
        "Distribution of mapped and unmapped sequences by confidence" if title else ""
    )

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
        font=dict(family="Arial,sans-serif", size=22, color="black"),
    )

    fig.update_xaxes(
        title=dict(font=dict(size=20)),
        dtick=0.1,
        range=[0, 1],
        showline=True,
        linecolor="black",
        linewidth=1,
    )

    fig.update_yaxes(
        title=dict(font=dict(size=20)),
        type="log",
        tickvals=[10, 100, 1000, 10000, 100000],
        ticktext=["10", "10²", "10³", "10⁴", "10⁵"],
        showgrid=False,
        gridwidth=1,
        gridcolor="white",
        showline=True,
        linecolor="black",
        linewidth=1,
        mirror=False,
    )

    fig.add_trace(
        go.Scatter(
            x=[None],
            y=[None],
            mode="markers",
            marker=dict(size=12, color="#AF6E7E", symbol="square"),
            name="overlap",
        )
    )

    fig.write_image(f"{folder}/{run}_confidence_distribution_range_mapped_unmapped.svg")


def fdr_ratio_mapped_unmapped(run, df, folder):

    bin_centers = []
    ratios = []

    for start in np.arange(0, 1, 0.05):
        end = start + 0.05
        subset_map = df[
            (df["mapped"] == "True") & (df["conf"] >= start) & (df["conf"] < end)
        ]
        subset_unmap = df[
            (df["mapped"] == "False") & (df["conf"] >= start) & (df["conf"] < end)
        ]

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
            # ticktext=[f"10^{i}" for i in range(-3, 2)],
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
    df["mapped"] = df["cleaned_preds"].apply(
        lambda x: "mapped" if x in reference else "unmapped"
    )

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

    # intersection_x = None
    # for i in range(1, len(hist_df)):
    #     mapped_prev, unmapped_prev = hist_df.iloc[i - 1][["Mapped", "Unmapped"]]
    #     mapped_curr, unmapped_curr = hist_df.iloc[i][["Mapped", "Unmapped"]]

    #     if mapped_prev < unmapped_prev and mapped_curr >= unmapped_curr:
    #         x0 = hist_df.iloc[i - 1]["confidence"]
    #         x1 = hist_df.iloc[i]["confidence"]

    #         y_diff_prev = mapped_prev - unmapped_prev
    #         y_diff_curr = mapped_curr - unmapped_curr

    #         intersection_x = x0 + (x1 - x0) * (-y_diff_prev) / (
    #             y_diff_curr - y_diff_prev
    #         )
    #         break

    fig = px.line(
        hist_df,
        x="confidence",
        y=["Mapped", "Unmapped"],
        markers=False,
        line_shape="linear",
        color_discrete_map={"Mapped": "#A5C8E1", "Unmapped": "#FFC48C"},
    )
    # orange - unmapped: #ff7f0e, blue: #1f77b4, brown #AF6E7E

    title_text = (
        "Relative distribution of mapped and unmapped peptides by confidence"
        if title
        else ""
    )

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
            alpha = 1 - (x_mid - 0.88) / (
                1 - 0.88
            )
            fillcolor = (
                f"rgba({base_color[0]}, {base_color[1]}, {base_color[2]}, {alpha:.2f})"
            )

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
    fig.for_each_trace(
        lambda t: t.update(name="mapped") if t.name == "Mapped" else None
    )

    fig.for_each_trace(
        lambda t: t.update(name="unmapped") if t.name == "Unmapped" else None
    )

    fig.write_image(f"{folder}/{run}_relative_mapped_unmapped_distribution.svg")


def plot_map_distribution(run, df, reference, folder, threshold, title=False):

    df = df[df["conf"] >= threshold].copy()

    df["mapped"] = df["cleaned_preds"].apply(
        lambda x: "mapped" if x in reference else "unmapped"
    )

    bins = np.arange(threshold, 1.002, 0.02)

    counts_mapped, edges = np.histogram(df[df["mapped"] == "mapped"]["conf"], bins=bins)
    counts_unmapped, _ = np.histogram(df[df["mapped"] == "unmapped"]["conf"], bins=bins)

    bin_centers = edges[:-1] + (0.5 * (edges[1] - edges[0]))

    hist_df = pd.DataFrame(
        {
            "confidence": np.tile(bin_centers, 2),
            "count": np.concatenate([counts_mapped, counts_unmapped]),
            "category": ["mapped"] * len(counts_mapped)
            + ["unmapped"] * len(counts_unmapped),
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

    title_text = (
        "Distribution of mapped and unmapped sequences by confidence" if title else ""
    )

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

    fig.for_each_trace(
        lambda t: t.update(name="mapped") if t.name == "Mapped" else None
    )

    fig.for_each_trace(
        lambda t: t.update(name="unmapped") if t.name == "Unmapped" else None
    )

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
        title="Confidence score distribution between {} and {}".format(
            min_conf, max_conf
        ),
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
    fig.write_image(
        f"{folder_figures}/confidence_distribution_range_{min_conf}_{max_conf}.png"
    )


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


def process_protein_contigs_scaffold(
    assembled_contigs, target_protein, max_mismatches, min_identity
):
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
            id=f"Contig {idx+1}",
            description=f"length: {len(contig)}, start: {start}, end: {end}, mismatches: {len(mismatches)}, identity: {identity:.2f}",
        )
        records.append(record)
    Bio.SeqIO.write(records, os.path.join(folder, f"{filename_prefix}.fasta"), "fasta")


def plot_contigs(mapped_contigs, prot_seq, title, output_file):
    sns.set("paper", "ticks", "colorblind", font_scale=1.5)
    _, ax = plt.subplots(figsize=(12, 4))

    ax.add_patch(
        patches.Rectangle(
            (0, 0), len(prot_seq), 0.2, facecolor="#e6f0ef", edgecolor="#e6f0ef"
        )
    )

    tracks = {}
    ind = 0

    for _, (contig, mapping) in tqdm(
        enumerate(mapped_contigs), desc="Plotting contigs"
    ):
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
    by_label = dict(zip(labels, handles))
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


def mapping_sequences(mapped_sequences, prot_seq, title, output_folder, output_file):
    """Plot sequences on a sequence using Plotly."""
    fig = go.Figure()

    fig.add_shape(
        type="rect",
        x0=0,
        x1=len(prot_seq),
        y0=0,
        y1=0.2,
        fillcolor="#e6f0ef",
        line=dict(width=0),
    )

    tracks = {}
    ind = 0

    for _, (_, mapping) in tqdm(enumerate(mapped_sequences), desc="Plotting contigs"):
        start_index, end_index, mismatches, _ = mapping

        ind += 1
        placed = False
        for track_num, track in tracks.items():
            if not any(s <= end_index <= e or s <= start_index <= e for s, e in track):
                track.append((start_index, end_index))
                fig.add_shape(
                    type="rect",
                    x0=start_index,
                    x1=end_index,
                    y0=0.3 + 0.1 * track_num,
                    y1=0.375 + 0.1 * track_num,
                    fillcolor="#007EA7",
                    line=dict(color="#007EA7"),
                )
                placed = True
                break

        if not placed:
            track_num = len(tracks) + 1
            tracks[track_num] = [(start_index, end_index)]
            fig.add_shape(
                type="rect",
                x0=start_index,
                x1=end_index,
                y0=0.3 + 0.1 * track_num,
                y1=0.375 + 0.1 * track_num,
                fillcolor="#007EA7",
                line=dict(color="#007EA7"),
            )

        for mismatch in mismatches:
            fig.add_shape(
                type="rect",
                x0=start_index + mismatch,
                x1=start_index + mismatch + 1,
                y0=0.3 + 0.1 * track_num,
                y1=0.375 + 0.1 * track_num,
                fillcolor="#FCAB64",
                line=dict(color="#FCAB64"),
            )

    print(f"Plotted {ind} sequences.")

    fig.add_trace(
        go.Scatter(
            x=[None],
            y=[None],
            mode="markers",
            marker=dict(size=10, color="#007EA7"),
            name="Match",
        )
    )
    fig.add_trace(
        go.Scatter(
            x=[None],
            y=[None],
            mode="markers",
            marker=dict(size=10, color="#FCAB64"),
            name="Mismatch",
        )
    )

    fig.update_layout(
        title=title,
        xaxis=dict(title="Sequence range", range=[0, len(prot_seq)], showgrid=False),
        yaxis=dict(
            showgrid=False,
            zeroline=False,
            showticklabels=False,
            range=[0, 0.3 + 0.1 * (len(tracks) + 1)],
        ),
        shapes=[],
        plot_bgcolor="white",
        width=1200,
        height=400,
    )

    fig.write_image(f"{output_folder}/{output_file}", scale=2)


def create_dataframe_from_mapped_sequences(data):
    """Takes a list of tuples containing sequence data and returns a structured DataFrame.
    """
    # Create the initial DataFrame
    df = pd.DataFrame(data, columns=["Sequence", "Details"])

    # Expand the 'Details' column into separate columns
    df[["start", "end", "mismatches_pos", "identity_score"]] = pd.DataFrame(
        df["Details"].tolist(), index=df.index
    )

    df.drop(columns=["Details"], inplace=True)
    df.rename(columns={"Sequence": "sequence"}, inplace=True)

    return df


def mapping_substitutions(
    mapped_sequences,
    prot_seq,
    title,
    bar_colors=None,
    output_file=None,
    output_folder=".",
    contig_colors="#6baed6",
    match_color="#6baed6",
):

    default_colors = {
        "match": match_color,
        "mismatch": "#b30000",
        "D_to_N": "#000000",
        "E_to_Q": "#A8A29E",
    }
    colors = {**default_colors, **(bar_colors or {})}

    fig = go.Figure()

    fig.add_shape(
        type="rect",
        x0=0,
        x1=len(prot_seq),
        y0=0,
        y1=0.2,
        fillcolor="#e6f0ef",
        line=dict(width=0),
    )

    tracks = {}
    ind = 0

    for seq, mapping in tqdm(mapped_sequences, desc="Plotting contigs"):
        start_index, end_index, mismatches, _ = mapping
        ind += 1

        contig_color = (
            contig_colors[ind % len(contig_colors)]
            if isinstance(contig_colors, list)
            else contig_colors
        )

        placed = False
        for track_num, track in tracks.items():
            if not any(s <= end_index <= e or s <= start_index <= e for s, e in track):
                track.append((start_index, end_index))
                y0 = 0.3 + 0.1 * track_num
                y1 = y0 + 0.075
                fig.add_shape(
                    type="rect",
                    x0=start_index,
                    x1=end_index,
                    y0=y0,
                    y1=y1,
                    fillcolor=contig_color,
                    line=dict(color=contig_color),
                )
                placed = True
                break

        if not placed:
            track_num = len(tracks)
            tracks[track_num] = [(start_index, end_index)]
            y0 = 0.3 + 0.1 * track_num
            y1 = y0 + 0.075
            fig.add_shape(
                type="rect",
                x0=start_index,
                x1=end_index,
                y0=y0,
                y1=y1,
                fillcolor=contig_color,
                line=dict(color=contig_color),
            )

        for mismatch in mismatches:
            abs_index = start_index + mismatch
            if abs_index >= len(prot_seq) or mismatch >= len(seq):
                continue

            ref_aa = prot_seq[abs_index]
            query_aa = seq[mismatch]

            if query_aa == "D" and ref_aa == "N":
                color = colors["D_to_N"]
            elif query_aa == "E" and ref_aa == "Q":
                color = colors["E_to_Q"]
            else:
                color = colors["mismatch"]

            fig.add_shape(
                type="rect",
                x0=abs_index,
                x1=abs_index + 1,
                y0=y0,
                y1=y1,
                fillcolor=color,
                line=dict(color=color),
            )

    print(f"Plotted {ind} sequences.")

    fig.add_trace(
        go.Scatter(
            x=[None],
            y=[None],
            mode="markers",
            marker=dict(size=10, color=colors["match"]),
            name="Match",
        )
    )
    fig.add_trace(
        go.Scatter(
            x=[None],
            y=[None],
            mode="markers",
            marker=dict(size=10, color=colors["mismatch"]),
            name="Mismatch",
        )
    )
    fig.add_trace(
        go.Scatter(
            x=[None],
            y=[None],
            mode="markers",
            marker=dict(size=10, color=colors["D_to_N"]),
            name="Seq:D → Ref:N",
        )
    )
    fig.add_trace(
        go.Scatter(
            x=[None],
            y=[None],
            mode="markers",
            marker=dict(size=10, color=colors["E_to_Q"]),
            name="Seq:E → Ref:Q",
        )
    )

    fig.update_layout(
        title=title,
        legend=dict(
            title=dict(text="Legend"),
            orientation="h",
            x=0.5,
            xanchor="center",
            y=1.05,
            yanchor="bottom",
        ),
        showlegend=True,
        xaxis=dict(title="Sequence range", range=[0, len(prot_seq)], showgrid=False),
        yaxis=dict(
            showgrid=False,
            zeroline=False,
            showticklabels=False,
            range=[0, 0.3 + 0.1 * (len(tracks) + 1)],
        ),
        plot_bgcolor="white",
        width=1200,
        height=400,
        font=dict(size=14, family="Arial, sans-serif", color="black"),
    )

    if output_file:
        os.makedirs(output_folder, exist_ok=True)
        fig.write_image(os.path.join(output_folder, output_file), scale=2)


def mapping_psms_protease_associated(
    mapped_sequences, prot_seq, labels, palette, title, output_folder, output_file
):

    fig = go.Figure()

    fig.add_shape(
        type="rect",
        x0=0,
        x1=len(prot_seq),
        y0=0,
        y1=0.2,
        fillcolor="#e6f0ef",
        line=dict(width=0),
    )
    tracks = {}
    ind = 0
    unique_labels = []
    for lab in labels:
        if lab not in unique_labels:
            unique_labels.append(lab)

    label_color = {lab: palette.get(lab, "#000000") for lab in unique_labels}

    for idx, (_, mapping) in tqdm(enumerate(mapped_sequences), desc="Plotting contigs"):
        start_index, end_index, mismatches, _ = mapping
        lab = labels[idx]
        ind += 1
        placed = False
        for track_num, track in tracks.items():
            if not any(s <= end_index <= e or s <= start_index <= e for s, e in track):
                track.append((start_index, end_index))
                y0 = 0.3 + 0.1 * track_num
                y1 = 0.375 + 0.1 * track_num
                fig.add_shape(
                    type="rect",
                    x0=start_index,
                    x1=end_index,
                    y0=y0,
                    y1=y1,
                    fillcolor=label_color[lab],
                    line=dict(color=label_color[lab]),
                )
                placed = True
                break
        if not placed:
            track_num = len(tracks) + 1
            tracks[track_num] = [(start_index, end_index)]
            y0 = 0.3 + 0.1 * track_num
            y1 = 0.375 + 0.1 * track_num
            fig.add_shape(
                type="rect",
                x0=start_index,
                x1=end_index,
                y0=y0,
                y1=y1,
                fillcolor=label_color[lab],
                line=dict(color=label_color[lab]),
            )

    print(f"Plotted {ind} sequences.")
    for lab, col in label_color.items():
        fig.add_trace(
            go.Scatter(
                x=[None],
                y=[None],
                mode="markers",
                marker=dict(size=10, color=col, symbol="square"),
                name=lab,
            )
        )

    fig.update_layout(
        title=title,
        legend_title="Proteases",
        legend=dict(
            orientation="h",
            x=0.5,
            y=1.1,
            xanchor="center",
            yanchor="bottom",
            font=dict(size=16),
        ),
        margin=dict(l=20, r=20, t=100, b=20),
        xaxis=dict(
            title="Reference",
            range=[0, len(prot_seq)],
            showgrid=False,
            dtick=50,
            tick0=0,
            tickfont=dict(size=16),
        ),
        yaxis=dict(
            showgrid=False,
            zeroline=False,
            showticklabels=False,
            range=[0, 0.3 + 0.1 * (len(tracks) + 1)],
        ),
        shapes=[],
        plot_bgcolor="white",
        width=1200,
        height=600,
        showlegend=True,
        font=dict(size=14, family="Arial, sans-serif", color="black"),
    )

    fig.write_image(f"{output_folder}/{output_file}", scale=2)