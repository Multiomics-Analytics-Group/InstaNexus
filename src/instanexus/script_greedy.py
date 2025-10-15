#!/usr/bin/env python

r"""Full assembly script for proteins.
 _____  _______  _    _
|  __ \|__   __|| |  | |
| |  | |  | |   | |  | |
| |  | |  | |   | |  | |
| |__| |  | |   | |__| |
|_____/   |_|   |______|

__authors__ = Marco Reverenna & Konstantinos Kalogeropoulus
__copyright__ = Copyright 2024-2025
__research-group__ = DTU Biosustain (Multi-omics Network Analytics) and DTU Bioengineering
__date__ = 15 Oct 2025
__maintainer__ = Marco Reverenna
__email__ = marcor@dtu.dk
__status__ = Dev
"""

# !pip install kaleido # to export plotly figures as png
# !pip install --upgrade nbformat # to avoid plotly error

import argparse
import json
import logging
import os
from pathlib import Path

import Bio
import pandas as pd

# import libraries
from . import alignment as align
from . import clustering as clus
from . import compute_statistics as comp_stat
from . import consensus as cons
from . import greedy_method as greedy
from . import mapping as map
from . import preprocessing as prep

repo_folder = Path(__file__).resolve().parents[2]

parser = argparse.ArgumentParser(description="Protein Assembly Script")
parser.add_argument("--input_csv", type=str, help="Input file")
parser.add_argument(
    "--chain",
    type=str,
    choices=["light", "heavy"],
    default="",
    help="Specify chain type if applicable (light or heavy). Leave empty if not applicable.",
)
parser.add_argument(
    "--folder_outputs", default="outputs", type=str, help="Outputs folder"
)
parser.add_argument(
    "--reference",
    action="store_true",
    help="Enable reference-based mode (use protein reference and compute statistics)",
    )

logging.basicConfig(
    level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s"
)
logger = logging.getLogger(__name__)

BASE_DIR = Path(__file__).resolve().parents[1]
JSON_DIR = BASE_DIR / "json"
# INPUT_DIR = BASE_DIR / "inputs"
# FASTA_DIR = BASE_DIR / "fasta"
# OUTPUTS_DIR = BASE_DIR / "outputs"


def get_sample_metadata(run, chain="", json_path=JSON_DIR / "sample_metadata.json"):
    with open(json_path, "r") as f:
        all_meta = json.load(f)

    if run not in all_meta:
        raise ValueError(f"Run '{run}' not found in metadata.")

    entries = all_meta[run]

    if not chain:
        # If no chain is specified, return the first entry
        return entries[0]

    for entry in entries:
        if entry["chain"] == chain:
            return entry

    raise ValueError(f"No metadata found for run '{run}' with chain '{chain}'.")


def main(
    input_csv: str,
    chain: str = "",
    folder_outputs: str = "outputs",
    reference: bool = False,
):
    """Main function to run the assembly script."""

    input_csv = Path(input_csv)

    logger.info("Starting protein assembly pipeline.")

    run = input_csv.stem

    if chain:
        meta = get_sample_metadata(run, chain=chain)
    else:
        meta = get_sample_metadata(run)

    if reference:
        protein = meta["protein"]
        proteases = meta["proteases"]

    ass_method = "greedy"

    conf = 0.92
    min_overlap = 3
    min_identity = 0.6
    max_mismatches = 14
    size_threshold = 10

    logger.info("Parameters loaded.")

    folder_outputs = Path(folder_outputs) / run
    folder_outputs.mkdir(parents=True, exist_ok=True)

    combination_folder_out = (
        folder_outputs
        / f"comb_{ass_method}_c{conf}_ts{size_threshold}_mo{min_overlap}_mi{min_identity}_mm{max_mismatches}"
    )
    prep.create_subdirectories_outputs(combination_folder_out)

    logger.info(f"Output folders created at: {combination_folder_out}")

    # Data cleaning
    logger.info("Starting data cleaning...")
    if reference:
        protein_norm = prep.normalize_sequence(protein)
    df = pd.read_csv(input_csv)

    df["protease"] = df["experiment_name"].apply(
            lambda name: prep.extract_protease(name, proteases)
        )
    df = prep.clean_dataframe(df)
    df["cleaned_preds"] = df["preds"].apply(prep.remove_modifications)
    cleaned_psms = df["cleaned_preds"].tolist()
    filtered_psms = prep.filter_contaminants(
        cleaned_psms, run, repo_folder / "fasta/contaminants.fasta"
    )
    df = df[df["cleaned_preds"].isin(filtered_psms)]
    
    if reference:
        df["mapped"] = df["cleaned_preds"].apply(
            lambda x: "True" if x in protein_norm else "False"
        )
    df = df[df["conf"] > conf]
    df.reset_index(drop=True, inplace=True)
    final_psms = df["cleaned_preds"].tolist()
    logger.info("Data cleaning completed.")

    # Assembly
    assembled_contigs = greedy.assemble_contigs(final_psms, min_overlap)
    assembled_contigs = list(set(assembled_contigs))
    assembled_contigs = [
        contig for contig in assembled_contigs if len(contig) > size_threshold
    ]
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
    if reference:
        mapped_contigs = map.process_protein_contigs_scaffold(
            assembled_contigs, protein_norm, max_mismatches, min_identity
        )
        df_contigs = map.create_dataframe_from_mapped_sequences(data=mapped_contigs)
        comp_stat.compute_assembly_statistics(
            df=df_contigs,
            sequence_type="contigs",
            output_folder=f"{combination_folder_out}/statistics",
            reference=protein_norm,
        )

    assembled_scaffolds = greedy.combine_seqs_into_scaffolds(
        assembled_contigs, min_overlap
    )
    assembled_scaffolds = list(set(assembled_scaffolds))
    assembled_scaffolds = sorted(assembled_scaffolds, key=len, reverse=True)
    assembled_scaffolds = [
        scaffold for scaffold in assembled_scaffolds if len(scaffold) > size_threshold
    ]
    assembled_scaffolds = greedy.combine_seqs_into_scaffolds(
        assembled_scaffolds, min_overlap
    )
    assembled_scaffolds = list(set(assembled_scaffolds))
    assembled_scaffolds = sorted(assembled_scaffolds, key=len, reverse=True)
    assembled_scaffolds = [
        scaffold for scaffold in assembled_scaffolds if len(scaffold) > size_threshold
    ]
    assembled_scaffolds = greedy.merge_contigs(assembled_scaffolds)
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
    if reference:
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
        )

    # Clustering
    scaffolds_folder_out = f"{combination_folder_out}/scaffolds"
    clus.cluster_fasta_files(input_folder=scaffolds_folder_out)

    cluster_tsv_folder = os.path.join(scaffolds_folder_out, "cluster")
    output_base_folder = os.path.join(scaffolds_folder_out, "cluster_fasta")

    for fasta_file in os.listdir(scaffolds_folder_out):
        if fasta_file.endswith(".fasta"):
            fasta_path = os.path.join(scaffolds_folder_out, fasta_file)
            clus.process_fasta_and_clusters(
                fasta_path, cluster_tsv_folder, output_base_folder
            )

    # Alignment
    cluster_fasta_folder = os.path.join(scaffolds_folder_out, "cluster_fasta")
    align_folder = os.path.join(scaffolds_folder_out, "align")
    prep.create_directory(align_folder)

    for cluster_folder in os.listdir(cluster_fasta_folder):
        cluster_folder_path = os.path.join(cluster_fasta_folder, cluster_folder)
        if os.path.isdir(cluster_folder_path):

            output_cluster_folder = os.path.join(align_folder, cluster_folder)
            os.makedirs(output_cluster_folder, exist_ok=True)

            for fasta_file in os.listdir(cluster_folder_path):
                if fasta_file.endswith(".fasta"):
                    fasta_file_path = os.path.join(cluster_folder_path, fasta_file)
                    base_filename = os.path.splitext(fasta_file)[0]
                    output_file = os.path.join(
                        output_cluster_folder, f"{base_filename}_out.afa"
                    )

                    align.align_or_copy_fasta(fasta_file_path, output_file)

    logger.info("All alignment tasks completed.")

    # Consensus
    consensus_folder = os.path.join(scaffolds_folder_out, "consensus")
    cons.process_alignment_files(align_folder, consensus_folder)


def cli():
    """Command-line interface entry point for dbg."""
    args = parser.parse_args()
    main(
        input_csv=args.input_csv,
        chain=args.chain,
        folder_outputs=args.folder_outputs,
        reference=args.reference,
    )

if __name__ == "__main__":
    cli()