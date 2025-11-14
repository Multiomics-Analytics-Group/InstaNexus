#!/usr/bin/env python

r"""Assembly module for InstaNexus.

 ██████████   ███████████ █████  █████
░░███░░░░███ ░█░░░███░░░█░░███  ░░███ 
 ░███   ░░███░   ░███  ░  ░███   ░███ 
 ░███    ░███    ░███     ░███   ░███ 
 ░███    ░███    ░███     ░███   ░███ 
 ░███    ███     ░███     ░███   ░███ 
 ██████████      █████    ░░████████  
░░░░░░░░░░      ░░░░░      ░░░░░░░░   
                          
__authors__ = Marco Reverenna & Konstantinos Kalogeropoulus
__copyright__ = Copyright 2024-2025
__research-group__ = DTU Biosustain (Multi-omics Network Analytics) and DTU Bioengineering
__date__ = 14 Nov 2025
__maintainer__ = Marco Reverenna
__email__ = marcor@dtu.dk
__status__ = Dev
"""

# import libraries
import logging
import networkx as nx
import pandas as pd
import argparse
import Bio

from . import helpers
from . import visualization as viz

from tqdm import tqdm
from pathlib import Path
from Bio import SeqIO
from collections import defaultdict
from collections import Counter
from itertools import combinations
from dataclasses import dataclass
from typing import List, Tuple, Dict, Any, Iterable, Optional


logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)



def find_peptide_overlaps(peptides, min_overlap):
    """Finds overlaps between peptide sequences using a greedy approach."""
    overlaps = defaultdict(list) 

    for index_a, peptide_a in tqdm(enumerate(peptides), desc="Finding overlaps"):
        for index_b, peptide_b in enumerate(peptides):
            if index_a != index_b:  # Skip comparing the same peptide

                max_possible_overlap = min(len(peptide_a), len(peptide_b))

                for overlap_length in range(min_overlap, max_possible_overlap):
                    if peptide_a[-overlap_length:] == peptide_b[:overlap_length]:
                        overlaps[index_a].append(
                            (index_b, overlap_length)
                        )  # Add the overlap to the dictionary
                    if peptide_b[-overlap_length:] == peptide_a[:overlap_length]:
                        overlaps[index_b].append(
                            (index_a, overlap_length)
                        )  # Add the overlap to the dictionary
    return overlaps


def assemble_contigs_greedy(peptides, min_overlap):
    """Assembles contigs from peptide sequences using a greedy approach."""
    assembled_contigs = peptides[:]
    iteration = 0

    while True:
        iteration += 1
        overlaps = find_peptide_overlaps(assembled_contigs, min_overlap)
        if not overlaps:
            break
        new_contigs = []
        used_indices = set()

        # Process overlaps deterministically
        for i in sorted(overlaps.keys()):  # ensure deterministic order
            if i in used_indices:
                continue
            # Sort overlaps_list deterministically: prioritize longer overlap, then lower index
            overlaps_list = sorted(overlaps[i], key=lambda x: (-x[1], x[0]))
            best_match = overlaps_list[0] if overlaps_list else None

            if best_match:
                j, overlap_len = best_match
                if j not in used_indices:
                    new_contig = (
                        assembled_contigs[i] + assembled_contigs[j][overlap_len:]
                    )
                    new_contigs.append(new_contig)
                    used_indices.update([i, j])
        # Add unused peptides
        remaining_contigs = [
            contig
            for idx, contig in enumerate(assembled_contigs)
            if idx not in used_indices
        ]
        assembled_contigs = new_contigs + remaining_contigs
        if len(new_contigs) == 0:
            break

    return assembled_contigs


def merge_contigs_greedy(contigs):
    """Merges overlapping contigs into a set of unique contigs."""
    contigs = sorted(contigs, key=len, reverse=True)
    merged = set(contigs)
    for c in tqdm(contigs, desc="Merging contigs"):
        # print(c)
        for c2 in contigs:
            if c != c2 and c2 in c:  # if c is a substring of c2
                merged.discard(c2)

    return list(merged)


def combine_seqs_into_scaffolds(contigs, min_overlap):
    """Combine contigs based on a minimum overlap length."""
    overlaps = find_overlaps(contigs, min_overlap=min_overlap)
    combined_contigs = []

    for a, b, overlap in overlaps:
        combined = a + b[overlap:]
        combined_contigs.append(combined)

    return combined_contigs + contigs


def scaffold_iterative_greedy(contigs, min_overlap, size_threshold, disable_tqdm=False):
    """Iterative scaffolding using Greedy approach."""
    def clean(seqs):
        """Remove duplicates, filter by length, and sort by descending size."""
        seqs = list(set(seqs))
        seqs = [s for s in seqs if len(s) > size_threshold]
        return sorted(seqs, key=len, reverse=True)
    current = clean(contigs)
    prev = None
    
    while prev != current:
        prev = current
        next_round = combine_seqs_into_scaffolds(current, min_overlap)
        next_round = merge_contigs_greedy(next_round)
        next_round = clean(next_round)
        if next_round == current:
            break
        current = next_round
    
    return current


def get_kmers(seqs, kmer_size):
    """Generate k-mers of specified length from input sequences."""
    kmers = []
    for seq in seqs:
        kmers.extend(seq[i:i+kmer_size] for i in range(len(seq)-kmer_size+1))
    return kmers


def get_kmer_counts(kmers):
    """Count occurrences of each k-mer in a list of k-mers; it takes a list of k-mers and returns a dictionary where the keys
    are unique k-mers and the values are the counts of each k-mer's occurrence.
    """
    kmer_counts = {}
    for kmer in kmers:
        if kmer in kmer_counts:
            kmer_counts[kmer] += 1
        else:
            kmer_counts[kmer] = 1
    return kmer_counts


def get_debruijn_edges_from_kmers(kmers):
    """Generate edges of a De Bruijn graph from a list of k-mers.
    """
    edges = set()
    k_1mers = defaultdict(set)
    for kmer in kmers:
        k_1mers[kmer[:-1]].add(kmer[1:])
    for prefix in k_1mers:
        for suffix in k_1mers[prefix]:
            edges.add((prefix, suffix))
    return edges


def assemble_contigs_dbg(edges):
    """Assemble contigs from De Bruijn graph edges by traversing the graph; it takes a set of directed edges representing
    a De Bruijn graph and assembles contigs by performing a depth-first traversal.
    """
    graph = defaultdict(list)
    for start, end in edges:
        graph[start].append(end)
    # find starting nodes (nodes with no incoming edges)
    all_ends = set(e for _, e in edges)
    start_nodes = set(graph.keys()) - all_ends

    def traverse_iterative(start_node):
        """Traverse a graph iteratively to find paths (contigs) starting from a given node."""
        stack = [(start_node, start_node)]
        visited = set()
        while stack:
            node, path = stack.pop()
            if node not in visited:
                visited.add(node)
                if node not in graph or not graph[node]:  # end of a path
                    contigs.append(path)
                else:
                    for next_node in graph[node]:
                        stack.append((next_node, path + next_node[-1]))

    contigs = []
    for start_node in tqdm(start_nodes, desc="Traversing nodes"):
        traverse_iterative(start_node)
    contigs = sorted(contigs, key=len, reverse=True)
    contigs = list(set(contigs))
    return contigs


def find_overlaps(contigs, min_overlap, disable_tqdm=False):
    """Find overlaps between pairs of contigs based on specified minimum overlap."""
    overlaps = []
    total_pairs = sum(
        1 for _ in combinations(contigs, 2)
    )
    with tqdm(total=total_pairs, desc="Finding overlaps", disable=disable_tqdm) as pbar:
        for a, b in combinations(
            contigs, 2
        ):  # combinations() generates all pairs of contigs
            for i in range(
                min_overlap, min(len(a), len(b)) + 1
            ):  # Check overlaps of different lengths
                if a[-i:] == b[:i]:
                    overlaps.append((a, b, i))
                if b[-i:] == a[:i]:
                    overlaps.append((b, a, i))
            pbar.update(1)

    return overlaps


def create_scaffolds(contigs, min_overlap, disable_tqdm=False):
    """Create scaffolds from a list of contigs by merging overlapping sequences."""
    overlaps = find_overlaps(
        contigs, min_overlap=min_overlap, disable_tqdm=disable_tqdm
    )
    combined_contigs = []
    for a, b, overlap in tqdm(
        overlaps, desc="Merging overlaps", total=len(overlaps), disable=disable_tqdm
    ):
        combined = a + b[overlap:]
        combined_contigs.append(combined)

    return combined_contigs + contigs


def merge_sequences_dbg(contigs, disable_tqdm=False):
    """Merges overlapping sequences."""
    contigs = sorted(contigs, key=len, reverse=True)
    merged = set(contigs)
    for c in tqdm(contigs, desc="Merging contigs", disable=disable_tqdm):
        for c2 in contigs:
            if c != c2 and c2 in c:  # if c2 is a substring of c
                merged.discard(c2)
    return list(merged)


def scaffold_iterative_dbg(contigs, min_overlap, size_threshold, disable_tqdm=False):
    """Iterative scaffolding using DBG approach."""
    prev = None
    current = contigs
    while prev != current:
        prev = current
        current = create_scaffolds(current, min_overlap, disable_tqdm)
        current = merge_sequences_dbg(current, disable_tqdm)
        current = list(set(current))
        current = [s for s in current if len(s) > size_threshold]
    return sorted(current, key=len, reverse=True)



def get_kmers(sequences: Iterable[str], kmer_size: int) -> List[str]:
    """Generate all k-mers from a list of sequences."""
    kmers = []
    for seq in sequences:
        if not seq:
            continue
        L = len(seq)
        if L < kmer_size:
            continue
        kmers.extend(seq[i:i+kmer_size] for i in range(L - kmer_size + 1))
    return kmers


def get_kmer_counts(kmers: Iterable[str]) -> Counter:
    """Return a Counter of k-mer frequencies."""
    return Counter(kmers)


def build_dbg_from_kmers(kmers: Iterable[str]) -> nx.DiGraph:
    """
    Build a weighted De Bruijn graph:
    - Nodes: (k-1)-mers
    - Edges: k-mers, with attributes:
        weight: int (count / coverage)
        kmers:  list of the k-mers mapping to this edge (optional, kept for debug)
    """
    G = nx.DiGraph()
    kmer_counts = Counter(kmers)
    for kmer, count in tqdm(kmer_counts.items(), desc="Building DBG"):
        prefix, suffix = kmer[:-1], kmer[1:]
        if G.has_edge(prefix, suffix):
            G[prefix][suffix]["weight"] += count
            G[prefix][suffix]["kmers"].append(kmer)
        else:
            G.add_edge(prefix, suffix, weight=count, kmers=[kmer])
    return G


def filter_low_weight_edges(G: nx.DiGraph, min_weight: int = 2) -> nx.DiGraph:
    """Remove edges with weight < min_weight (light error correction)."""
    to_remove = [(u, v) for u, v, d in G.edges(data=True) if d.get("weight", 0) < min_weight]
    G.remove_edges_from(to_remove)
    # drop isolated nodes
    iso = [n for n in G.nodes if G.in_degree(n) == 0 and G.out_degree(n) == 0]
    G.remove_nodes_from(iso)
    return G


@dataclass
class ContigPath:
    nodes: List[str]         # list of (k-1)-mer node labels in path order
    seq: str                 # assembled sequence
    weights: List[int]       # edge weights along the path

def _extend_linear_path(G: nx.DiGraph, start: str, succ: str) -> ContigPath:
    """Extend from start→succ while in/out-degree == 1 (unbranched)."""
    path_nodes = [start, succ]
    weights = [G[start][succ]["weight"]]
    cur = succ
    while G.in_degree(cur) == 1 and G.out_degree(cur) == 1:
        nxt = next(iter(G.successors(cur)), None)
        if nxt is None:
            break
        weights.append(G[cur][nxt]["weight"])
        path_nodes.append(nxt)
        cur = nxt

    # build sequence from node labels
    seq = path_nodes[0]
    for n in path_nodes[1:]:
        seq += n[-1]
    return ContigPath(nodes=path_nodes, seq=seq, weights=weights)


def assemble_contigs_dbgx(G: nx.DiGraph, min_length: int = 0) -> List[ContigPath]:
    """
    Collapse unbranched paths into contigs.
    Returns ContigPath items (with sequence + per-edge weights).
    """
    contigs: List[ContigPath] = []

    # start from "branch" nodes (sources/sinks/branch points)
    for node in tqdm(list(G.nodes), desc="Assembling contigs"):
        if G.out_degree(node) == 0:
            continue
        if G.in_degree(node) != 1 or G.out_degree(node) != 1:
            for succ in G.successors(node):
                cp = _extend_linear_path(G, node, succ)
                if len(cp.seq) >= min_length:
                    contigs.append(cp)

    # edge case: pure cycles (every node deg=1/1). Traverse any cycle once.
    if not contigs and len(G) > 0:
        # pick arbitrary node and walk until it closes
        start = next(iter(G.nodes))
        succs = list(G.successors(start))
        if succs:
            cp = _extend_linear_path(G, start, succs[0])
            if len(cp.seq) >= min_length:
                contigs.append(cp)

    # deduplicate by sequence
    unique: Dict[str, ContigPath] = {}
    for cp in contigs:
        if cp.seq not in unique or len(cp.seq) > len(unique[cp.seq].seq):
            unique[cp.seq] = cp

    contigs = sorted(unique.values(), key=lambda c: len(c.seq), reverse=True)
    return contigs


@dataclass
class ContigScore:
    seq: str
    length: int
    mean_weight: float
    min_weight: int
    max_weight: int
    score: float

def score_contig(cp: ContigPath, alpha_len: float = 1.0, alpha_cov: float = 1.0, alpha_min: float = 0.2) -> ContigScore:
    """
    Simple reference-free score combining length and coverage:
      score = alpha_len * log(length) + alpha_cov * mean_weight + alpha_min * min_weight
    Adjust alphas to your data. You can also plug-in intensity-based terms later.
    """
    import math
    if cp.weights:
        mean_w = sum(cp.weights) / len(cp.weights)
        min_w = min(cp.weights)
        max_w = max(cp.weights)
    else:
        mean_w = min_w = max_w = 0.0
    L = len(cp.seq)
    composite = alpha_len * math.log(max(L, 2)) + alpha_cov * mean_w + alpha_min * min_w
    return ContigScore(seq=cp.seq, length=L, mean_weight=mean_w, min_weight=min_w, max_weight=max_w, score=composite)


def rank_contigs_by_score(contigs: List[ContigPath],
                          alpha_len: float = 1.0,
                          alpha_cov: float = 1.0,
                          alpha_min: float = 0.2) -> List[ContigScore]:
    scored = [score_contig(c, alpha_len, alpha_cov, alpha_min) for c in contigs]
    return sorted(scored, key=lambda s: (s.score, s.length), reverse=True)


def scaffold_iterative_dbgx(seqs: List[str],
                            kmer_size: int,
                            size_threshold: int = 10,
                            min_weight: int = 2,
                            max_rounds: int = 5,
                            patience: int = 2,
                            alpha_len: float = 1.0,
                            alpha_cov: float = 1.0,
                            alpha_min: float = 0.2) -> List[str]:
    """
    Optional refinement:
      rebuild DBG from current contigs → collapse → filter by size → repeat
    Stops when no improvement for `patience` rounds or `max_rounds` reached.
    """
    best: List[str] = list(seqs)
    no_improve = 0

    for rnd in range(1, max_rounds + 1):
        kmers = get_kmers(seqs, kmer_size)
        if not kmers:
            break
        G = build_dbg_from_kmers(kmers)
        G = filter_low_weight_edges(G, min_weight=min_weight)
        contigs = assemble_contigs_dbgx(G, min_length=size_threshold)

        if not contigs:
            break

        ranked = rank_contigs_by_score(contigs, alpha_len, alpha_cov, alpha_min)
        seqs_new = [r.seq for r in ranked]

        # improvement heuristic: fewer contigs or longer top contig
        improved = (len(seqs_new) < len(best)) or (seqs_new and best and len(seqs_new[0]) > len(best[0]))
        if improved:
            best = seqs_new
            no_improve = 0
        else:
            no_improve += 1

        seqs = seqs_new
        if no_improve >= patience:
            break

    # final unique & size filter
    uniq = []
    seen = set()
    for s in best:
        if len(s) >= size_threshold and s not in seen:
            seen.add(s)
            uniq.append(s)
    return uniq


def extend_path_dbg(G, contig, k, min_weight=1):
    """
    Extend a contig in both directions along the DBG G, following dominant edges
    """
    seq = contig
    extended = True

    # Extension forwards
    while extended:
        extended = False
        suffix = seq[-(k-1):]
        if suffix not in G or G.out_degree(suffix) == 0:
            break

        successors = list(G.successors(suffix))
        if len(successors) > 1:
            # choose only if there is a clear dominant edge
            best_succ, best_w = None, 0
            for s in successors:
                w = G[suffix][s].get("weight", 0)
                if w > best_w:
                    best_succ, best_w = s, w
            if best_succ and best_w >= min_weight:
                seq += best_succ[-1]
                extended = True
            else:
                break
        else:
            nxt = successors[0]
            w = G[suffix][nxt].get("weight", 0)
            if w >= min_weight:
                seq += nxt[-1]
                extended = True

    # Extension backwards
    extended = True
    while extended:
        extended = False
        prefix = seq[:k-1]
        if prefix not in G or G.in_degree(prefix) == 0:
            break

        predecessors = list(G.predecessors(prefix))
        if len(predecessors) > 1:
            best_pred, best_w = None, 0
            for p in predecessors:
                w = G[p][prefix].get("weight", 0)
                if w > best_w:
                    best_pred, best_w = p, w
            if best_pred and best_w >= min_weight:
                seq = best_pred[0] + seq
                extended = True
            else:
                break
        else:
            p = predecessors[0]
            w = G[p][prefix].get("weight", 0)
            if w >= min_weight:
                seq = p[0] + seq
                extended = True

    return seq


class Assembler:
    """Unified assembler supporting 'greedy', 'dbg' and 'dbgx' (NetworkX) modes."""

    def __init__(self,
                 mode: str = "greedy",
                 min_overlap: int = 4,
                 size_threshold: int = 10,
                 kmer_size: int = 6,
                 min_identity: float = 0.8,
                 max_mismatches: int = 10,
                 # dbgx-specific:
                 min_weight: int = 2,
                 refine_rounds: int = 0,   # 0 = no refine, >0 enables iterative refine
                 refine_patience: int = 2,
                 alpha_len: float = 1.0,
                 alpha_cov: float = 1.0,
                 alpha_min: float = 0.2):
        if mode not in ["greedy", "dbg", "dbg_weighted", "dbgX", "fusion"]:
            raise ValueError("mode must be 'greedy', 'dbg', 'dbg_weighted', 'dbgX' or 'fusion'")

        self.mode = mode
        self.min_overlap = min_overlap
        self.size_threshold = size_threshold
        self.kmer_size = kmer_size
        self.min_identity = min_identity
        self.max_mismatches = max_mismatches

        # dbgx params
        self.min_weight = min_weight
        self.refine_rounds = refine_rounds
        self.refine_patience = refine_patience
        self.alpha_len = alpha_len
        self.alpha_cov = alpha_cov
        self.alpha_min = alpha_min

    def assemble_greedy(self, sequences):

        logger.info(f"[Assembler] Running Greedy assembly (min_overlap={self.min_overlap})")
        contigs = assemble_contigs_greedy(sequences, self.min_overlap)
        contigs = list(set(contigs))
        contigs = sorted(contigs, key=len, reverse=True)

        scaffolds =  scaffold_iterative_greedy(contigs, self.min_overlap, self.size_threshold)

        return scaffolds

    def assemble_dbg(self, sequences):
        logger.info(f"[Assembler] Running DBG assembly (kmer_size={self.kmer_size})")

        kmers = get_kmers(sequences, self.kmer_size)
        edges = get_debruijn_edges_from_kmers(kmers)
        contigs = assemble_contigs_dbg(edges)
        contigs = list(set(contigs))
        contigs = sorted(contigs, key=len, reverse=True)
        contigs = [seq for seq in contigs if len(seq) > self.size_threshold]

        scaffolds = scaffold_iterative_dbg(contigs, self.min_overlap, self.size_threshold)

        return scaffolds
    

    def assemble_dbg_weighted(self, sequences: List[str]) -> List[str]:
        logger.info(f"[Assembler] Running DBG weighted (k={self.kmer_size}, min_weight={self.min_weight}, refine_rounds={self.refine_rounds})")

        kmers = get_kmers(sequences, self.kmer_size)
        if not kmers:
            logger.warning("No kmers generated; returning empty result.")
            return []

        G = build_dbg_from_kmers(kmers)
        G = filter_low_weight_edges(G, min_weight=self.min_weight)

        contigs_cp = assemble_contigs_dbgx(G, min_length=self.size_threshold)
        if not contigs_cp:
            logger.warning("No contigs assembled from DBGX; returning empty result.")
            return []

        ranked = rank_contigs_by_score(contigs_cp, self.alpha_len, self.alpha_cov, self.alpha_min)
        contigs = [r.seq for r in ranked]

        if self.refine_rounds and self.refine_rounds > 0:
            contigs = scaffold_iterative_dbgx(
                contigs,
                kmer_size=self.kmer_size,
                size_threshold=self.size_threshold,
                min_weight=self.min_weight,
                max_rounds=self.refine_rounds,
                patience=self.refine_patience,
                alpha_len=self.alpha_len,
                alpha_cov=self.alpha_cov,
                alpha_min=self.alpha_min,
            )

        scaffolds = list(contigs)

        return scaffolds
    
    def assemble_dbgX(self, sequences):
        logger.info(f"[Assembler] Running DBG-Extension (k={self.kmer_size})")

        kmers = get_kmers(sequences, self.kmer_size)
        G = build_dbg_from_kmers(kmers)
        G = filter_low_weight_edges(G, min_weight=self.min_weight)

        contigs_cp = assemble_contigs_dbgx(G, min_length=self.size_threshold)
        contigs = [c.seq for c in contigs_cp]

        logger.info("Extending contigs using DBG paths (coverage-aware)...")
        extended_contigs = [extend_path_dbg(G, c, self.kmer_size, self.min_weight) for c in contigs]
        extended_contigs = sorted(set(extended_contigs), key=len, reverse=True)

        return extended_contigs
    
    def assemble_fusion(self, sequences):
        logger.info(f"[Assembler] Running FUSION (DBG weighted + greedy merge)")

        contigs_dbg_weighted = self.assemble_dbg_weighted(sequences)

        logger.info("Running greedy merge on DBG weighted contigs...")
        contigs_greedy = assemble_contigs_greedy(sequences, self.min_overlap)
        contigs_greedy = merge_contigs_greedy(contigs_greedy)

        combined = list(set(contigs_dbg_weighted + contigs_greedy))
        combined = [s for s in combined if len(s) > self.size_threshold]
        logger.info(f"Combined {len(combined)} contigs from DBG weighted + Greedy")

        fused = assemble_contigs_greedy(combined, self.min_overlap)
        fused = merge_contigs_greedy(fused)
        fused = [s for s in fused if len(s) > self.size_threshold]
        fused = sorted(set(fused), key=len, reverse=True)

        return fused


    def run(self, sequences: List[str]):
        if not sequences:
            logger.error("No valid sequences provided for assembly.")
            raise ValueError("Input sequences list is empty.")
        
        if self.mode == "greedy":
            return self.assemble_greedy(sequences)
        elif self.mode == "dbg":
            return self.assemble_dbg(sequences)
        elif self.mode == "dbg_weighted":
            return self.assemble_dbg_weighted(sequences)
        elif self.mode == "dbgX":
            return self.assemble_dbgX(sequences)
        elif self.mode == "fusion":
            return self.assemble_fusion(sequences)


def main(
    input_csv_path: str,
    output_scaffolds_path: str,
    metadata_json_path: str,
    assembly_mode: str,
    kmer_size: int,
    min_overlap: int,
    size_threshold: int,
    reference: bool,
    chain: str,    
    min_identity: float,
    max_mismatches: int,
):
    """Main function for standalone assembly."""

    protein_norm = None # None means no reference mode
    if reference:
        logger.info(f"Reference mode enabled. Loading reference protein...")
        if not metadata_json_path:
            raise ValueError("metadata_json_path is required when reference mode is enabled.")
        
        try:
            run_name = Path(input_csv_path).stem # extract run name from input file
            meta = helpers.get_sample_metadata(run=run_name, chain=chain, json_path=metadata_json_path)
            protein= meta["protein"]
            protein_norm = helpers.get_normalized_protein(protein)
            logger.info("Reference protein loaded and normalized successfully.")
        
        except Exception as e:
            logger.error(f"Failed to get reference protein: {e}")
            logger.warning("Disabling reference mode.")
            reference = False # if fails, disable reference mode

    input_data = Path(input_csv_path)

    print(f"Starting assembly pipeline with {assembly_mode}...")

    # Load input sequences
    df = pd.read_csv(input_data)

    if "cleaned_preds" in df.columns:
        sequences = df["cleaned_preds"].dropna().tolist()
    else:
        raise ValueError("CSV must contain a 'cleaned_preds' column.")

    assembler = Assembler(
        mode=assembly_mode,
        min_overlap=min_overlap,
        size_threshold=size_threshold,
        kmer_size=kmer_size,
        min_identity=min_identity,
        max_mismatches=max_mismatches
    )

    scaffolds = assembler.run(sequences = sequences)

    output_path = Path(output_scaffolds_path)
    output_folder = output_path.parent.mkdir(parents=True, exist_ok=True)

    records = [
        Bio.SeqRecord.SeqRecord(
            Bio.Seq.Seq(seq), id=f"scaffold_{i+1}", description=f"length: {len(seq)}"
        )
        for i, seq in enumerate(scaffolds)
    ]

    Bio.SeqIO.write(
        records,
        output_path,
        "fasta",
    )

    logger.info(f"Assembly completed — {len(scaffolds)} scaffolds saved to {output_path}")

    if protein_norm:
        logger.info("Reference mode: calculating statistics...")
        stats_path = output_path.parent / "statistics"
        stats_path.mkdir(parents=True, exist_ok=True)

        mapped_scaffolds = viz.process_protein_contigs_scaffold(
            assembled_contigs=scaffolds,
            target_protein=protein_norm,
            max_mismatches=max_mismatches,
            min_identity=min_identity, 
        )
        df_scaffolds_mapped = viz.create_dataframe_from_mapped_sequences(
            data=mapped_scaffolds
        )
        helpers.compute_assembly_statistics(
            df=df_scaffolds_mapped,
            sequence_type="scaffolds",
            output_folder=str(stats_path),
            reference=protein_norm,
        )
        logger.info(f"Reference mode: Statistics saved to {stats_path}")


def cli():
    """Command-line interface for the assembly module."""

    parser = argparse.ArgumentParser(
        description="Run Greedy or DBG assembly on peptide sequences."
    )
    parser.add_argument(
        "--input-csv-path",
        type=str,
        required=True,
        help="Path to input CSV or FASTA file containing sequences.",
    )
    parser.add_argument(
        "--output-scaffolds-path",
        type=str,
        required=True,
        help="Path to save the output scaffolds FASTA file.",
    )
    parser.add_argument(
        "--metadata-json-path",
        type=str,
        default=None, # Optional, but required by --reference
        help="Path to sample_metadata.json (required for --reference).",
    )
    parser.add_argument(
        "--assembly-mode",
        type=str,
        choices=["greedy", "dbg"],
        default="greedy",
        help="Assembly mode to use: greedy or dbg.",
    )
    parser.add_argument(
        "--kmer-size",
        type=int,
        default=6,
        help="K-mer size (used only for DBG mode).",
    )
    parser.add_argument(
        "--min-overlap",
        type=int,
        default=3,
        help="Minimum overlap for merging sequences.",
    )
    parser.add_argument(
        "--size-threshold",
        type=int,
        default=10,
        help="Minimum contig length to retain after assembly.",
    )
    parser.add_argument(
        "--reference",
        action="store_true",
        help="Enables reference-based statistics.",
    )

    parser.add_argument(
        "--chain",
        type=str,
        default="",
        help="Specify chain type (light/heavy) required for reference lookup."
    )
    parser.add_argument(
        "--min-identity",
        type=float,
        default=0.8,
        help="Minimum identity for reference mapping."
    )
    parser.add_argument(
        "--max-mismatches",
        type=int,
        default=10,
        help="Maximum mismatches for reference mapping."
    )

    args = parser.parse_args()

    # in case of non-DBG mode, ignore kmer_size
    if args.assembly_mode != "dbg":
        args.kmer_size = 0
        logger.info("Ignoring kmer_size (used only for DBG mode).")
    
    if args.reference and not args.metadata_json_path:
        parser.error("--metadata-json-path is required when --reference is enabled.")

    main(**vars(args))


if __name__ == "__main__":
    cli()

# python -m instanexus.assembly --input-csv-path outputs/bsa/bsa_cleaned.csv --output-scaffolds-path outputs/bsa/bsa_scaffolds.fasta --metadata-json-path json/sample_metadata.json --assembly-mode dbg_weighted --kmer-size 7 --min-overlap 3 --size-threshold 10 --min-identity 0.8 --max-mismatches 10
