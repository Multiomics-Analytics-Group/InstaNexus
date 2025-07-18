#!/usr/bin/env python

r"""
 _____  _______  _    _ 
|  __ \|__   __|| |  | |
| |  | |  | |   | |  | |
| |  | |  | |   | |  | |
| |__| |  | |   | |__| |
|_____/   |_|   |______|

__authors__ = Marco Reverenna & Konstantinos Kalogeropoulus
__copyright__ = Copyright 2024-2025
__research-group__ = DTU Biosustain (Multi-omics Network Analytics) and DTU Bioengineering
__date__ = 26 Nov 2024
__maintainer__ = Marco Reverenna
__email__ = marcor@dtu.dk
__status__ = Dev
"""

# import libraries
import os
import re
import numpy as np
import pandas as pd
from Bio import SeqIO
import plotly.express as px
import plotly.graph_objects as go


# Define and create the necessary directories only if they don't exist
def create_directory(path):
    """Creates a directory if it does not already exist.
    Args:
        path (str): The path of the directory to create.
    """
    if not os.path.exists(path):
        os.makedirs(path)
        #print(f"Created: {path}")
    #else:
        #print(f"Already exists: {path}")

def create_subdirectories_outputs(folder):
    """Creates subdirectories within the specified folder.
    Args:
        folder (str): The path of the parent directory.
    """
    subdirectories = ['preprocessing', 'contigs', 'scaffolds', 'consensus', 'logs', 'statistics']
    for subdirectory in subdirectories:
        create_directory(f"{folder}/{subdirectory}")

def create_subdirectories_figures(folder):
    """Creates subdirectories within the specified folder.
    Args:
        folder (str): The path of the parent directory.
    """
    subdirectories = ['preprocessing', 'contigs', 'scaffolds', 'consensus', 'heatmap', 'logo']
    for subdirectory in subdirectories:
        create_directory(f"{folder}/{subdirectory}")


def normalize_sequence(sequence):
    """
    Normalize the given amino acid sequence by replacing all occurrences of 'I' with 'L'.

    Parameters:
    sequence (str): The amino acid sequence to be normalized.

    Returns:
    str: The normalized amino acid sequence with 'I' replaced by 'L'.
    """
    return sequence.replace('I', 'L')



def remove_modifications(psm_column):
    """
    Remove any content within parentheses, including the parentheses, from a given string.

    Parameters:
    - psm_column (str): The string containing modifications in parentheses (e.g., "A(ox)BC(mod)D"). If the value is null, it returns None.

    Returns:
    - str: The string with all parenthetical modifications removed (e.g., "ABCD"), or None if the input was null.
    """

    if pd.notnull(psm_column):
        return re.sub(r'\(.*?\)', '', psm_column)  # Replace any content in parentheses with an empty string
    return None

def clean_dataframe(df):
    """
    Clean and preprocess a DataFrame for analysis by removing '(ox)' substrings from sequences in the 'seq' column.
    by replacing values of -1 with -10 in the 'log_probs' column, by dropping rows with missing values in the 'preds' column.
    by extracts a 'protease' value from the 'experiment_name' column based on a specific naming convention.
    by adding a 'conf' column, which is the exponentiated 'log_probs' to represent confidence and sorting 
    the DataFrame by the 'conf' column in descending order.

    Parameters:
    - df (DataFrame): The raw input DataFrame to clean.

    Returns:
    - DataFrame: The cleaned and processed DataFrame.
    """
    # update 1
    df = df.copy()

    df['log_probs'] = df['log_probs'].replace(-1, -10)
    # -10 is very low, replacing with -1 (so that we are sure is very low quality prediction)
    # check new update InstaNovo 

    df = df.dropna(subset=['preds'])

    #df['protease'] = df['experiment_name'].apply(lambda x: x.split('_')[-3] if isinstance(x, str) else None)
    #df.loc[:, 'protease'] = df['experiment_name'].apply(lambda x: x.split('_')[-3] if isinstance(x, str) else None)
    
    #df['conf'] = np.exp(df['log_probs'])
    
    df.loc[:, 'conf'] = np.exp(df['log_probs'])

    df = df.sort_values('conf', ascending=False)

    return df



def filter_contaminants(seqs, run, contaminants_fasta):
    """
    Filters out sequences from the input list `seqs` that are substrings of sequences
    in the contaminants file. If run == 'bsa', the Bovine serum albumin precursor is ignored.

    Parameters:
    - seqs (list of str): List of sequences to be filtered.
    - contaminants_fasta (str): Path to the FASTA file containing contaminant sequences.
    - run (str): Run identifier, used to control special filtering logic.
    """

    contam_records = []
    for record in SeqIO.parse(contaminants_fasta, "fasta"):
        if run == 'bsa' and "Bovine serum albumin precursor" in record.description:
            continue  # Skip BSA if run is 'bsa'
        contam_records.append(str(record.seq))

    filtered_seqs = []
    removed_count = 0

    for seq in seqs:
        if any(seq in contam_seq for contam_seq in contam_records):
            removed_count += 1
        else:
            filtered_seqs.append(seq)

    #print(f"Removed {removed_count} contaminant sequences, {len(filtered_seqs)} sequences remaining.")
    return filtered_seqs


def plot_confidence_distribution(df, folder_figures, min_conf=0, max_conf=1):
    """
    Plots the distribution of confidence scores from a DataFrame.

    Parameters:
    df (pandas.DataFrame): DataFrame containing a column 'conf' with confidence scores.
    min_conf (float, optional): Minimum value of confidence range (default is 0).
    max_conf (float, optional): Maximum value of confidence range (default is 1).
    """
    # Filter the data based on the specified range
    filtered_df = df[(df['conf'] >= min_conf) & (df['conf'] <= max_conf)]

    fig = go.Figure()

    fig.add_trace(go.Histogram(
        x=filtered_df['conf'],
        xbins=dict(start=min_conf, end=max_conf, size=(max_conf - min_conf) / 40),
        marker=dict(color='brown'),
        opacity=1  # Remove opacity
    ))

    # Configure layout
    fig.update_layout(title='Confidence score distribution between {} and {}'.format(min_conf, max_conf),
        xaxis_title='Values', yaxis_title='Frequency',
        bargap=0.1, height=700, width=1200,
        margin=dict(l=50, r=50, t=100, b=100)
    )

    # Configure axes
    fig.update_xaxes(showgrid=True, gridcolor='lightgray',
        ticklabelposition='outside bottom',  # Place labels outside the bottom of the axis
        dtick=0.02  # Set the distance between ticks
    )
    fig.update_yaxes(showgrid=True, gridcolor='lightgray')
    fig.write_image(f"{folder_figures}/confidence_distribution_range_{min_conf}_{max_conf}.png")
    #fig.show()


def extract_protease(experiment_name, proteases):
    """ Extracts the protease name from the given experiment name.
    Parameters:
        experiment_name (str): The name of the experiment.
        proteases (list or set): A list or set of known protease names.

    Returns:
        str or None: The matched protease name, or None if no match is found.
    """
    parts = experiment_name.split('_')
    for part in parts:
        if part in proteases:
            return part
    return None



def plot_protease_distribution(protease_counts, folder_figures):
    """Creates an interactive bar plot of protease distribution using Plotly.
    
    Parameters:
        protease_counts (pandas.Series): A Pandas Series with protease names as the index
                                         and their counts as the values.
    """
    # Convert the Series to a DataFrame for compatibility with Plotly
    protease_df = protease_counts.reset_index()
    protease_df.columns = ['Protease', 'Count']

    fig = px.bar(
        protease_df,
        x='Protease',
        y='Count',
        title='Proteases distribution',
        labels={'Protease': 'Proteases', 'Count': 'Counts'},
        text='Count'
    )

    # Reduce the width of the bars and set the text position outside the bars
    fig.update_traces(textposition='outside', width=0.4)

    # Define conversion factor from mm to pixels (approx. 3.78 px/mm at 96 DPI)
    mm_to_px = 3.78

    # Set the desired figure dimensions according to Nature's standards.
    # For a double-column figure:
    width_mm = 240   # width in mm (double column) 183 per Nature
    height_mm = 200  # full page depth in mm 247 per Nature

    fig.update_layout(
        width=int(width_mm * mm_to_px),
        height=int(height_mm * mm_to_px),
        xaxis_title='Proteases',
        yaxis_title='Counts',
        xaxis_tickangle=0,
        showlegend=False,
        title_font_size=12,
        font=dict(
            family="Arial, sans-serif",  # Change to your preferred font
            size=8,
            color="black"
        ),
        margin=dict(t=50, b=50, l=50, r=100),  # Adjust margins as needed
        plot_bgcolor='white',   # White background for the plot area
        paper_bgcolor='white'   # White background for the entire figure
    )

    # Export the figure in SVG vector format
    fig.write_image(f"{folder_figures}/proteases_distribution.svg")
    fig.show()


def missing_values_barplot(run, dataframe, folder):

    dataframe["missing_preds"] = dataframe["preds"].isna()

    missing_counts_df = dataframe["missing_preds"].value_counts().reset_index()
    missing_counts_df.columns = ["PSMs", "Count"]
    missing_counts_df["PSMs"] = missing_counts_df["PSMs"].map({True: "Missing", False: "Valid"})
    
    missing_counts_df = missing_counts_df.sort_values("PSMs", key=lambda x: x.map({"Valid": 0, "Missing": 1}))

    valid_count = missing_counts_df.loc[missing_counts_df["PSMs"]=="Valid", "Count"].values[0] \
                  if "Valid" in missing_counts_df["PSMs"].values else 0
    missing_count = missing_counts_df.loc[missing_counts_df["PSMs"]=="Missing", "Count"].values[0] \
                    if "Missing" in missing_counts_df["PSMs"].values else 0

    total = valid_count + missing_count
    valid_pct = (valid_count / total * 100) if total > 0 else 0
    missing_pct = (missing_count / total * 100) if total > 0 else 0

    print(f"Valid PSMs: {valid_count} ({valid_pct:.2f}%)")
    print(f"Missing PSMs: {missing_count} ({missing_pct:.2f}%)")
    print(f"Total PSMs: {total}")

    color_map = {"Valid": "#337AB7", "Missing": "#FF5733"}  

    fig = go.Figure()
    
    fig.add_trace(go.Bar(
        x=["Valid"],
        y=[valid_count],
        marker_color=color_map["Valid"],
        width=[0.4],  
        name="Valid"
    ))
    
    fig.add_trace(go.Bar(
        x=["Missing"],
        y=[missing_count],
        marker_color=color_map["Missing"],
        width=[0.4],  
        name="Missing"
    ))
    
    fig.update_layout(
        title="",
        xaxis_title="PSMs",
        yaxis_title="Count",
        barmode="group",
        width=800,
        height=600,
        margin=dict(t=50, l=10, r=10, b=10),
        legend_title_text="",
        paper_bgcolor='white',
        plot_bgcolor='white',
        font=dict(size=22, family="Arial, sans-serif", color='black')
    )
    
    fig.update_xaxes(
        showline=True,
        linewidth=1,
        linecolor='black',
        title_font=dict(size=20),  
        tickfont=dict(size=18)  
    )
    
    fig.update_yaxes(
        showline=True,
        linewidth=1,
        linecolor='black',
        title_font=dict(size=20),  
        tickfont=dict(size=18)
    )
    
    fig.write_image(f"{folder}/{run}_missing_value_bar.svg")
    fig.show()


def plot_map_unmap_distribution(df, reference, run, folder, conf_lim, title=False):

    df = df[df['conf'] >= conf_lim]    

    df["mapped"] = df["cleaned_preds"].apply(lambda x: "mapped" if x in reference else "unmapped")

    fig = px.histogram(
        df, 
        x="conf", 
        color="mapped", 
        nbins=int((1 - conf_lim)*50), 
        barmode="overlay", 
        color_discrete_map={"mapped": "#1f77b4", "unmapped": "#ff7f0e"}
        # orange: #ff7f0e, blue: #1f77b4, brown #AF6E7E
    )
    
    fig.update_traces(xbins=dict(start=0, end=1, size=0.01))

    title_text = f"Distribution of mapped and unmapped sequences by confidence" if title else ''

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
        font=dict(family="Arial,sans-serif",size = 22, color='black'),
    )

    fig.update_xaxes(
        title=dict(font = dict(size=20)),
        dtick=0.1,  
        range=[0, 1],
        showline=True,
        linecolor='black',
        linewidth=1,  
    )

    fig.update_yaxes(
        title=dict(font = dict(size=20)),
        type="log",
        tickvals=[10, 100, 1000, 10000, 100000], 
        ticktext=["10", "10²", "10³", "10⁴", "10⁵"], 
        showgrid=False,
        gridwidth=1,
        gridcolor="white",
        showline=True,
        linecolor='black',
        linewidth=1,
        mirror = False 
    )

    fig.add_trace(go.Scatter(
        x=[None],
        y=[None],
        mode='markers',
        marker=dict(size=12, color="#AF6E7E", symbol = "square"
                    ),
        name="overlap"
    ))
    
    fig.write_image(f"{folder}/{run}_confidence_distribution_range_mapped_unmapped.svg")
    fig.show()


def fdr_ratio_mapped_unmapped(run, df, folder):

    bin_centers = []
    ratios = []

    for start in np.arange(0, 1, 0.05):
        end = start + 0.05
        subset_map = df[(df['mapped'] == 'True') & (df['conf'] >= start) & (df['conf'] < end)]
        subset_unmap = df[(df['mapped'] == 'False') & (df['conf'] >= start) & (df['conf'] < end)]
        
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
        if (ratios[i] <= y_horizontal and ratios[i + 1] >= y_horizontal) or \
        (ratios[i] >= y_horizontal and ratios[i + 1] <= y_horizontal):
            x1, x2 = bin_centers[i], bin_centers[i + 1]
            y1, y2 = ratios[i], ratios[i + 1]
            intersect_x = x1 + (y_horizontal - y1) * (x2 - x1) / (y2 - y1)
            break  

    fig = go.Figure()

    fig.add_trace(go.Scatter(
        x=bin_centers,
        y=ratios,
        mode='lines+markers',
        line=dict(color='#1f77b4', width=3.5),
        name='mapped/unmapped ratio'
    ))

    if intersect_x is not None:
        fig.add_shape(
            type="line",
            x0=-0.025,
            x1=intersect_x,  
            xref="x",
            y0=y_horizontal,
            y1=y_horizontal,
            yref="y",      
            line=dict(color="black", width=1.5, dash="dash")
        )

    if intersect_x is not None:
        fig.add_shape(
            type="line",
            x0=intersect_x,
            x1=intersect_x,
            y0=0.001, 
            y1=y_horizontal,
            yref="y",
            line=dict(color="black", width=1.5, dash="dash")
        )
        
        fig.add_annotation(
            x=intersect_x - 0,
            y=y_horizontal - 0.8, 
            text=f"Confidence: {intersect_x:.2f}",
            showarrow=False,
            font=dict(size=16, color="black")
        )

    fig.update_layout(
        xaxis_title="Confidence",
        yaxis_title="Ratio mapped/unmapped",
        template="plotly_white",
        height=600,
        width=800,
        font=dict(size=22, family="Arial, sans-serif", color="black"),
        xaxis=dict(
            title=dict(font = dict(size=20)),
            tickmode="array",
            tickvals=np.arange(0, 1.1, 0.1),  
            ticktext=[f"{x:.1f}" for x in np.arange(0, 1.1, 0.1)],
            zeroline=False,
            linewidth=1,
            linecolor='black',
            showline=True,
            showgrid=False,
        ),
        yaxis=dict(
            title=dict(font = dict(size=20)),
            type="log", 
            range=[-3, 1.5],  
            tickvals=[10**i for i in range(-3, 2)],
            ticktext=["10⁻³", "10⁻²", "10⁻¹", "10⁰", "10¹"],
            #ticktext=[f"10^{i}" for i in range(-3, 2)],  
            showline=True,  
            linewidth=1,
            linecolor='black',
            zeroline=False,
            showgrid=False,
            side="left" 
        )
    )

    fig.write_image(f"{folder}/{run}_fdr_ratio_mapped_unmapped.svg")
    fig.show()



def plot_relative_map_distribution(run, df, reference, folder, title=False):

    df = df[df['conf'] >= 0].copy()
    df["mapped"] = df["cleaned_preds"].apply(lambda x: "mapped" if x in reference else "unmapped")

    bins = np.arange(0, 1, 0.05)
    #bins = np.arange(0, 1.02, 0.02)
    df['bin'] = pd.cut(df['conf'], bins=bins, labels=bins[:-1])

    bin_counts = df.groupby('bin')['mapped'].count()
    mapped_counts = df[df['mapped'] == 'mapped'].groupby('bin')['mapped'].count()
    unmapped_counts = df[df['mapped'] == 'unmapped'].groupby('bin')['mapped'].count()

    mapped_percentages = (mapped_counts / bin_counts) * 100
    unmapped_percentages = (unmapped_counts / bin_counts) * 100

    hist_df = pd.DataFrame({
        "confidence": bins[:-1],
        "Mapped": mapped_percentages.fillna(0).values,
        "Unmapped": unmapped_percentages.fillna(0).values
    })

    intersection_x = None
    for i in range(1, len(hist_df)):
        mapped_prev, unmapped_prev = hist_df.iloc[i - 1][["Mapped", "Unmapped"]]
        mapped_curr, unmapped_curr = hist_df.iloc[i][["Mapped", "Unmapped"]]

        if mapped_prev < unmapped_prev and mapped_curr >= unmapped_curr:
            x0 = hist_df.iloc[i - 1]["confidence"]
            x1 = hist_df.iloc[i]["confidence"]

            y_diff_prev = mapped_prev - unmapped_prev
            y_diff_curr = mapped_curr - unmapped_curr

            intersection_x = x0 + (x1 - x0) * (-y_diff_prev) / (y_diff_curr - y_diff_prev)
            break

    fig = px.line(hist_df,
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
        legend=dict(font=dict(size=16))
    )


    fig.update_yaxes(range=[0, 100],
                     title=dict(font = dict(size=20)),
                     showline=True,
                     linecolor="black",
                     linewidth=1,
                     showgrid=False)

    fig.update_xaxes(range=[0, 1],
                     title=dict(font = dict(size=20)),
                     showline=True,
                     linecolor="black",
                     linewidth=1,
                     showgrid=False,
                     tickmode="linear",
                     dtick=0.1,
                     tickangle=0)

    fig.add_vline(
        x=0.88, 
        line_width=2,
        line_dash="dash",
        line_color="black",
        annotation_text="Cutoff",
        annotation_position="top",
        annotation_font_size=16
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
            alpha = 1 - (x_mid - 0.88) / (1 - 0.88)  # Use 0.88 as reference for transparency
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
                layer="below"
            )

    fig.update_traces(line=dict(width=3.5))
    fig.for_each_trace(lambda t: t.update(name="mapped") if t.name == "Mapped" else None)

    fig.for_each_trace(lambda t: t.update(name="unmapped") if t.name == "Unmapped" else None)


    fig.write_image(f"{folder}/{run}_relative_mapped_unmapped_distribution.svg")
    fig.show()


def plot_map_distribution(run, df, reference, folder, threshold, title=False):

    df = df[df['conf'] >= threshold].copy()
    
    df["mapped"] = df["cleaned_preds"].apply(lambda x: "mapped" if x in reference else "unmapped")
    
    bins = np.arange(threshold, 1.002, 0.02)  
    
    counts_mapped, edges = np.histogram(df[df["mapped"] == "mapped"]["conf"], bins=bins)
    counts_unmapped, _ = np.histogram(df[df["mapped"] == "unmapped"]["conf"], bins=bins)
    
    bin_centers = edges[:-1] + (0.5 * (edges[1] - edges[0]))

    hist_df = pd.DataFrame({
        "confidence": np.tile(bin_centers, 2),  # Repeat values for both categories
        "count": np.concatenate([counts_mapped, counts_unmapped]),  # Stack counts
        "category": ["mapped"] * len(counts_mapped) + ["unmapped"] * len(counts_unmapped)  # Category labels
    })
    
    fig = px.bar(hist_df,
                 x="confidence",
                 y="count",
                 color="category",
                 color_discrete_map={"mapped": "#A5C8E1", "unmapped": "#FFC48C"},
                 barmode="stack")  
    
    title_text = "Distribution of mapped and unmapped sequences by confidence" if title else ''
    
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
        font=dict(size=22, family="Arial, sans-serif", color='black')
    )
    
    fig.update_yaxes(
        title=dict(font=dict(size=20)),
        showline=True,
        linecolor='black',
        linewidth=1,
        showgrid=False
    )

    fig.update_xaxes(
        title=dict(font = dict(size=20)),
        showline=True,
        linecolor='black',
        linewidth=1,
        showgrid=False,
        tickmode="linear",
        dtick=0.02
    )

    fig.for_each_trace(
    lambda t: t.update(name="mapped") if t.name == "Mapped" else None
    )

    fig.for_each_trace(
    lambda t: t.update(name="unmapped") if t.name == "Unmapped" else None
    )

    fig.write_image(f"{folder}/{run}_psms_mapped_unmapped_distribution.svg")
    fig.show()