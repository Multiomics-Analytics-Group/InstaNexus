#!/usr/bin/env python

r"""Full assembly script for proteins.
 _____  _______  _    _
|  __ \|__   __|| |  | |
| |  | |  | |   | |  | |
| |  | |  | |   | |  | |
| |__| |  | |   | |__| |
|_____/   |_|   |______|

__authors__ = Marco Reverenna
__copyright__ = Copyright 2025-2026
__research-group__ = DTU Biosustain (Multi-omics Network Analytics) and DTU Bioengineering
__date__ = 26 Jun 2025
__maintainer__ = Marco Reverenna
__email__ = marcor@dtu.dk
__status__ = Dev
"""

# !pip install kaleido # to export plotly figures as png
# !pip install --upgrade nbformat # to avoid plotly error

import json
import os
from pathlib import Path

import Bio
import pandas as pd

# my modules
from src import compute_statistics as comp_stat
from src import dbg
from src import mapping as map
from src import preprocessing as prep

BASE_DIR = Path(__file__).resolve().parents[2]
JSON_DIR = BASE_DIR / "json"
INPUT_DIR = BASE_DIR / "inputs"
FASTA_DIR = BASE_DIR / "fasta"
OUTPUTS_DIR = BASE_DIR / "outputs"


def get_sample_metadata(run, chain="", json_path=JSON_DIR / "sample_metadata.json"):
    with open(json_path, "r") as f:
        all_meta = json.load(f)

    if run not in all_meta:
        raise ValueError(f"Run '{run}' not found in metadata.")

    entries = all_meta[run]

    for entry in entries:
        if entry["chain"] == chain:
            return entry

    raise ValueError(f"No metadata found for run '{run}' with chain '{chain}'.")


def run_pipeline_dbg(
    conf, kmer_size, min_overlap, max_mismatches, min_identity, size_threshold
):

    ass_method = "dbg"
    run = "ma1"

    meta = get_sample_metadata(run, chain="light")
    protein = meta["protein"]
    chain = meta["chain"]
    proteases = meta["proteases"]

    # run = "pa"
    # chain = "heavy"

    params = {
        "ass_method": "dbg",
        "conf": conf,
        "kmer_size": kmer_size,
        "min_overlap": min_overlap,
        "min_identity": min_identity,
        "max_mismatches": max_mismatches,
        "size_threshold": size_threshold,
    }

    # Directories
    folder_outputs = OUTPUTS_DIR / f"{run}{chain}"
    Path(folder_outputs).mkdir(parents=True, exist_ok=True)

    prep.create_directory(folder_outputs)
    combination_folder_out = os.path.join(
        folder_outputs,
        f"comb_{ass_method}_c{conf}_ks{kmer_size}_ts{size_threshold}_mo{min_overlap}_mi{min_identity}_mm{max_mismatches}",
    )

    prep.create_subdirectories_outputs(combination_folder_out)

    # Data cleaning
    protein_norm = prep.normalize_sequence(protein)

    df = pd.read_csv(INPUT_DIR / f"{run}.csv")

    df["protease"] = df["experiment_name"].apply(
        lambda name: prep.extract_protease(name, proteases)
    )

    df = prep.clean_dataframe(df)

    df["cleaned_preds"] = df["preds"].apply(prep.remove_modifications)

    cleaned_psms = df["cleaned_preds"].tolist()

    filtered_psms = prep.filter_contaminants(
        cleaned_psms, run, FASTA_DIR / "contaminants.fasta"
    )

    df = df[df["cleaned_preds"].isin(filtered_psms)]

    df["mapped"] = df["cleaned_preds"].apply(
        lambda x: "True" if x in protein_norm else "False"
    )

    df = df[df["conf"] > conf]

    df.reset_index(drop=True, inplace=True)

    final_psms = df["cleaned_preds"].tolist()

    # Assembly
    kmers = dbg.get_kmers(final_psms, kmer_size=kmer_size)

    edges = dbg.get_debruijn_edges_from_kmers(kmers)

    assembled_contigs = dbg.assemble_contigs(edges)

    assembled_contigs = sorted(assembled_contigs, key=len, reverse=True)

    assembled_contigs = list(set(assembled_contigs))

    assembled_contigs = [seq for seq in assembled_contigs if len(seq) > size_threshold]

    assembled_contigs = sorted(assembled_contigs, key=len, reverse=True)

    records = [
        Bio.SeqRecord.SeqRecord(
            Bio.Seq.Seq(contig),
            id=f"contig_{idx+1}",
            description=f"length: {len(contig)}",
        )
        for idx, contig in enumerate(assembled_contigs)
    ]

    Bio.SeqIO.write(
        records,
        f"{combination_folder_out}/contigs/{ass_method}_contig_{conf}_{run}.fasta",
        "fasta",
    )

    mapped_contigs = map.process_protein_contigs_scaffold(
        assembled_contigs, protein_norm, max_mismatches, min_identity
    )

    df_contigs = map.create_dataframe_from_mapped_sequences(data=mapped_contigs)

    comp_stat.compute_assembly_statistics(
        df=df_contigs,
        sequence_type="contigs",
        output_folder=f"{combination_folder_out}/statistics",
        reference=protein_norm,
        **params,
    )

    assembled_scaffolds = dbg.create_scaffolds(assembled_contigs, min_overlap)

    assembled_scaffolds = list(set(assembled_scaffolds))

    assembled_scaffolds = sorted(assembled_scaffolds, key=len, reverse=True)

    assembled_scaffolds = [
        scaffold for scaffold in assembled_scaffolds if len(scaffold) > size_threshold
    ]

    assembled_scaffolds = dbg.merge_sequences(assembled_scaffolds)

    assembled_scaffolds = list(set(assembled_scaffolds))

    assembled_scaffolds = sorted(assembled_scaffolds, key=len, reverse=True)

    assembled_scaffolds = [
        scaffold for scaffold in assembled_scaffolds if len(scaffold) > size_threshold
    ]

    records = []
    for i, seq in enumerate(assembled_scaffolds):
        record = Bio.SeqRecord.SeqRecord(
            Bio.Seq.Seq(seq), id=f"scaffold_{i+1}", description=f"length: {len(seq)}"
        )
        records.append(record)

    Bio.SeqIO.write(
        records,
        f"{combination_folder_out}/scaffolds/{ass_method}_scaffold_{conf}_{run}.fasta",
        "fasta",
    )

    mapped_scaffolds = map.process_protein_contigs_scaffold(
        assembled_contigs=assembled_scaffolds,
        target_protein=protein_norm,
        max_mismatches=max_mismatches,
        min_identity=min_identity,
    )

    df_scaffolds_mapped = map.create_dataframe_from_mapped_sequences(
        data=mapped_scaffolds
    )

    comp_stat.compute_assembly_statistics(
        df=df_scaffolds_mapped,
        sequence_type="scaffolds",
        output_folder=f"{combination_folder_out}/statistics",
        reference=protein_norm,
        **params,
    )
