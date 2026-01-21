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
import argparse
import ast
import logging
import math
from collections import Counter, defaultdict
from dataclasses import dataclass
from itertools import combinations
from pathlib import Path
from typing import Dict, Iterable, List, Optional

import Bio
import networkx as nx
import pandas as pd
from tqdm import tqdm

from . import helpers, preprocessing
from . import visualization as viz

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

MAX_REFINE_ROUNDS = 10


# def find_peptide_overlaps(peptides, min_overlap):
#     """Finds overlaps between peptide sequences using a greedy approach."""
#     overlaps = defaultdict(list)

#     for index_a, peptide_a in tqdm(enumerate(peptides), desc="Finding overlaps"):
#         for index_b, peptide_b in enumerate(peptides):
#             if index_a != index_b:  # Skip comparing the same peptide
#                 max_possible_overlap = min(len(peptide_a), len(peptide_b))

#                 for overlap_length in range(min_overlap, max_possible_overlap):
#                     if peptide_a[-overlap_length:] == peptide_b[:overlap_length]:
#                         overlaps[index_a].append((index_b, overlap_length))  # Add the overlap to the dictionary
#                     if peptide_b[-overlap_length:] == peptide_a[:overlap_length]:
#                         overlaps[index_b].append((index_a, overlap_length))  # Add the overlap to the dictionary
#     return overlaps


def find_sliding_overlaps(sequences: list, min_overlap: int):
    overlaps = []
    for i, seq_a in enumerate(sequences):
        for j, seq_b in enumerate(sequences):
            if i == j: continue
                        
            max_search = min(len(seq_a), len(seq_b))
            for length in range(max_search, min_overlap - 1, -1):
                s1 = seq_a[-length:]
                s2 = seq_b[:length]
                
                if s1 == s2:
                    overlaps.append((i, j, length))
                    break
                
                diff = sum(1 for a, b in zip(s1, s2) if a != b)
                if diff == 1 and length >= 5: 
                    logger.info(f"POTENTIAL OVERLAP MISSED: {s1} vs {s2}")
    return overlaps


def merge_with_overhang(seq_a, seq_b, overlap_len):
    """
    Finds the exact alignment point and merges.
    This handles: A='...KGR', B='SVKGR...' -> Result: 'SVKGR...' 
    (if B contains A) or the correct concatenation.
    """
    if seq_a in seq_b:
        return seq_b
    if seq_b in seq_a:
        return seq_a
        
    for i in range(len(seq_a) - overlap_len + 1):
        suffix = seq_a[i:]
        if seq_b.startswith(suffix):
            return seq_a[:i] + seq_b
            
    return seq_a + seq_b[overlap_len:]


def assemble_contigs_greedy(peptides, min_overlap):
    assembled_contigs = peptides[:]
    iteration = 0
    MAX_ITERATIONS = 50

    while iteration < MAX_ITERATIONS:
        iteration += 1
        overlaps_list = find_sliding_overlaps(assembled_contigs, min_overlap)
        
        if not overlaps_list:
            break
            
        overlaps_list.sort(key=lambda x: x[2], reverse=True)
        
        new_contigs = []
        used_indices = set()

        for i, j, overlap_len in overlaps_list:
            if i in used_indices or j in used_indices:
                continue
            
            new_contig = merge_with_overhang(assembled_contigs[i], assembled_contigs[j], overlap_len)
            new_contigs.append(new_contig)
            used_indices.update([i, j])

        if not new_contigs:
            break

        remaining = [c for idx, c in enumerate(assembled_contigs) if idx not in used_indices]
        assembled_contigs = new_contigs + remaining

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


# def combine_seqs_into_scaffolds(contigs, min_overlap):
#     """Combine contigs based on a minimum overlap length."""
#     overlaps = find_overlaps(contigs, min_overlap=min_overlap, disable_tqdm=True)
#     combined_contigs = []

#     for a, b, overlap in overlaps:
#         combined = a + b[overlap:]
#         combined_contigs.append(combined)

#     return combined_contigs + contigs

def combine_seqs_into_scaffolds(contigs, min_overlap):
    """Combine contigs using the new sliding logic and overhang-aware merging."""
    # CAMBIO 1: Usa la logica sliding del tuo branch
    overlaps = find_sliding_overlaps(contigs, min_overlap=min_overlap)
    combined_contigs = []

    for i, j, overlap_len in overlaps:
        # CAMBIO 2: Usa il merger intelligente che abbiamo scritto
        a = contigs[i]
        b = contigs[j]
        combined = merge_with_overhang(a, b, overlap_len)
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
    MAX_ROUNDS = 10

    logger.info(f"Starting iterative scaffolding (Max rounds: {MAX_ROUNDS})...")

    for i in range(MAX_ROUNDS):
        next_round = combine_seqs_into_scaffolds(current, min_overlap)
        next_round = merge_contigs_greedy(next_round)
        next_round = clean(next_round)

        if len(next_round) == len(current):
            break

        current = next_round
        logger.info(f"  Round {i + 1}: {len(current)} contigs")

    return current


def get_weighted_kmers_from_df(
    df: pd.DataFrame,
    kmer_size: int,
    use_abundance: bool = True,
    use_quality: bool = True,
) -> Counter:
    """
    Generates k-mer weights integrating multiple data modalities.

    SCORING LOGIC:
    1. Sequence Confidence (Deep Learning): Geometric Mean of 'instanovo_token_log_probabilities'.
    2. Abundance (MS1): Logarithmic scaling of 'peptide_abundance'.
    3. Quality (MS2): Linear boost from 'ion_match_intensity'.
    4. Physics (iRT): Exponential penalty for 'iRT error'.
    """
    kmer_weights = Counter()

    col_tokens = "instanovo_token_log_probabilities"
    col_ms1_abundance = "peptide_abundance" 
    col_ms2_intensity = "ion_match_intensity"
    col_irt = "iRT error"

    for _, row in df.iterrows():
        sequence = row.get("cleaned_preds")

        if not isinstance(sequence, str) or len(sequence) < kmer_size:
            continue

        log_probs = []
        raw_probs_str = row.get(col_tokens)

        if isinstance(raw_probs_str, str) and raw_probs_str.startswith("["):
            try:
                parsed_log_probs = ast.literal_eval(raw_probs_str)

                if len(parsed_log_probs) == len(sequence) + 2:
                    log_probs = parsed_log_probs[1:-1]
                elif len(parsed_log_probs) == len(sequence):
                    log_probs = parsed_log_probs
                else:
                    log_probs = [0.0] * len(sequence)
            except Exception:
                log_probs = [0.0] * len(sequence)
        else:
            log_probs = [0.0] * len(sequence)

        global_multiplier = 1.0

        if use_abundance:
            ms1_val = row.get(col_ms1_abundance, 0)
            if pd.notnull(ms1_val) and ms1_val > 0:
                base_weight = math.log10(float(ms1_val) + 10)
            else:
                base_weight = 1.0

            ms2_val = row.get(col_ms2_intensity, 0)
            ms2_boost = 1.0
            if pd.notnull(ms2_val):
                # Linear boost: 0.1 -> 1.3x, 0.9 -> 3.7x
                ms2_boost = 1.0 + (float(ms2_val) * 3.0)

            global_multiplier = base_weight * ms2_boost

        # C. Quality Penalty: iRT Error
        if use_quality:
            irt_err = row.get(col_irt, None)
            is_missing = row.get("is_missing_irt_error", False)
            pred_irt = row.get("predicted iRT", 0)

            # Apply penalty ONLY if data is valid and error exists
            # We ignore missing flags and dummy values (-30)
            if pd.notnull(irt_err) and not is_missing and pred_irt > -20:
                try:
                    # Exponential decay (Sigma = 20.0 to handle long tails)
                    global_multiplier *= math.exp(-abs(float(irt_err)) / 20.0)
                except ValueError:
                    pass

        # --- 3. K-MER WEIGHTING (Geometric Mean) ---
        for i in range(len(sequence) - kmer_size + 1):
            kmer = sequence[i : i + kmer_size]

            # Extract log-probs for this specific k-mer
            if i + kmer_size <= len(log_probs):
                sub_logs = log_probs[i : i + kmer_size]

                # Geometric Mean Calculation
                # Math: GeoMean(p1...pn) = exp( average(ln(p1)...ln(pn)) )
                if sub_logs:
                    avg_log_prob = sum(sub_logs) / len(sub_logs)
                    local_confidence = math.exp(avg_log_prob)
                else:
                    local_confidence = 1.0
            else:
                local_confidence = 1.0

            # Final Weight = Global (MS1/MS2/iRT) * Local (Token Prob)
            weight = global_multiplier * local_confidence

            kmer_weights[kmer] += weight

    return kmer_weights


def get_debruijn_edges_from_kmers(kmers):
    """Generate edges of a De Bruijn graph from a list of k-mers."""
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
    total_pairs = sum(1 for _ in combinations(contigs, 2))
    with tqdm(total=total_pairs, desc="Finding overlaps", disable=disable_tqdm) as pbar:
        for a, b in combinations(contigs, 2):  # combinations() generates all pairs of contigs
            for i in range(min_overlap, min(len(a), len(b)) + 1):  # Check overlaps of different lengths
                if a[-i:] == b[:i]:
                    overlaps.append((a, b, i))
                if b[-i:] == a[:i]:
                    overlaps.append((b, a, i))
            pbar.update(1)

    return overlaps


def create_scaffolds(contigs, min_overlap, disable_tqdm=False):
    """
    Improved version: uses sliding overlaps and overhang-aware merging.
    """
    # Usa find_sliding_overlaps che restituisce (i, j, overlap_len)
    overlaps = find_sliding_overlaps(contigs, min_overlap=min_overlap)
    combined_contigs = []
    
    for i, j, overlap_len in tqdm(overlaps, desc="Merging overlaps", disable=disable_tqdm):
        a = contigs[i]
        b = contigs[j]
        # Usa il merger intelligente che gestisce gli overhangs
        combined = merge_with_overhang(a, b, overlap_len)
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
    """Generate all k-mers from a list of sequences (preserves duplicates)."""
    kmers = []
    for seq in sequences:
        if not seq:
            continue
        L = len(seq)
        if L < kmer_size:
            continue
        kmers.extend(seq[i : i + kmer_size] for i in range(L - kmer_size + 1))
    return kmers


def get_kmer_counts(kmers: Iterable[str]) -> Counter:
    """Return a Counter of k-mer frequencies."""
    return Counter(kmers)


def build_dbg_from_kmers(kmers: Iterable[str], weights: Counter = None) -> nx.DiGraph:
    """
    Build a De Bruijn graph.
    If 'weights' (Counter) is provided, uses those values for edges.
    Otherwise, counts occurrences from the list.
    """
    G = nx.DiGraph()

    if weights:
        iterator = weights.items()  # (kmer, calculated_weight)
    else:
        iterator = Counter(kmers).items()  # (kmer, count)

    for kmer, weight in iterator:
        prefix, suffix = kmer[:-1], kmer[1:]
        if G.has_edge(prefix, suffix):
            G[prefix][suffix]["weight"] += weight
        else:
            G.add_edge(prefix, suffix, weight=weight)
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
    nodes: List[str]  # list of (k-1)-mer node labels in path order
    seq: str  # assembled sequence
    weights: List[int]  # edge weights along the path


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


def score_contig(
    cp: ContigPath,
    alpha_len: float = 1.0,
    alpha_cov: float = 1.0,
    alpha_min: float = 0.2,
) -> ContigScore:
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
    return ContigScore(
        seq=cp.seq,
        length=L,
        mean_weight=mean_w,
        min_weight=min_w,
        max_weight=max_w,
        score=composite,
    )


def rank_contigs_by_score(
    contigs: List[ContigPath],
    alpha_len: float = 1.0,
    alpha_cov: float = 1.0,
    alpha_min: float = 0.2,
) -> List[ContigScore]:
    scored = [score_contig(c, alpha_len, alpha_cov, alpha_min) for c in contigs]
    return sorted(scored, key=lambda s: (s.score, s.length), reverse=True)


def build_overlap_graph(contigs: List[str], min_overlap: int) -> nx.DiGraph:
    G = nx.DiGraph()
    for i, seq in enumerate(contigs):
        G.add_node(i, seq=seq, length=len(seq))

    n_contigs = len(contigs)
    for i in range(n_contigs):
        for j in range(n_contigs):
            if i == j: continue

            seq_a, seq_b = contigs[i], contigs[j]
            # Cerchiamo l'overlap usando la logica sliding che abbiamo validato
            # Questo permette a A='...KGR' di connettersi a B='SVKGR...'
            best_overlap = 0
            max_ov = min(len(seq_a), len(seq_b))
            
            for k in range(max_ov, min_overlap - 1, -1):
                # Caso 1: La fine di A è contenuta nell'inizio di B (Sliding)
                if seq_b.startswith(seq_a[-k:]) or seq_a.endswith(seq_b[:k]):
                    best_overlap = k
                    break

            if best_overlap > 0:
                G.add_edge(i, j, weight=best_overlap)
    return G


def merge_paths_from_overlap_graph(G: nx.DiGraph) -> List[str]:
    """
    Traverses the overlap graph to find and merge the optimal paths
    (heaviest overlaps), resolving branches by prioritizing longer overlaps.
    """
    merged_contigs = []
    G_work = G.copy()

    while G_work.number_of_nodes() > 0:
        start_nodes = [n for n in G_work.nodes if G_work.in_degree(n) == 0]

        if not start_nodes:
            start_node = max(G_work.nodes, key=lambda n: len(G_work.nodes[n]["seq"]))
        else:
            start_node = max(start_nodes, key=lambda n: len(G_work.nodes[n]["seq"]))

        path = [start_node]
        current = start_node

        while True:
            if G_work.out_degree(current) == 0:
                break

            neighbors = list(G_work.successors(current))
            best_next = max(neighbors, key=lambda n: G_work[current][n]["weight"])

            if best_next in path:
                break

            path.append(best_next)
            current = best_next

        if len(path) == 1:
            merged_contigs.append(G_work.nodes[start_node]["seq"])
        else:
            first_idx = path[0]
            super_seq = G_work.nodes[first_idx]["seq"]

            for i in range(len(path) - 1):
                u, v = path[i], path[i + 1]
                overlap_len = G_work[u][v]["weight"]
                seq_v = G_work.nodes[v]["seq"]
                super_seq = merge_with_overhang(super_seq, seq_v, overlap_len)

            merged_contigs.append(super_seq)

        G_work.remove_nodes_from(path)

    return merged_contigs


def refine_using_overlap_graph(contigs: List[str], min_overlap: int) -> List[str]:
    """
    Wrapper function: Builds graph -> Merges paths -> Cleans up substrings.
    """
    if not contigs:
        return []

    G = build_overlap_graph(contigs, min_overlap)
    logger.info(f"DEBUG: Overlap Graph has {G.number_of_nodes()} nodes and {G.number_of_edges()} edges.")
    
    if G.number_of_edges() == 0:
        logger.warning("No overlaps found between scaffolds! Check sliding logic.")    
    
    refined = merge_paths_from_overlap_graph(G)

    refined = sorted(list(set(refined)), key=len, reverse=True)
    final_set = []
    for seq in refined:
        if not any(seq in other and seq != other for other in refined):
            final_set.append(seq)

    return final_set


def scaffold_iterative_dbgx(
    seqs: List[str],
    kmer_size: int,
    size_threshold: int = 10,
    min_weight: int = 2,
    max_rounds: int = 5,
    patience: int = 2,
    alpha_len: float = 1.0,
    alpha_cov: float = 1.0,
    alpha_min: float = 0.2,
) -> List[str]:
    """
    Optional refinement:
      rebuild DBG from current contigs → collapse → filter by size → repeat
    Stops when no improvement for `patience` rounds or `max_rounds` reached.
    """
    best: List[str] = list(seqs)
    no_improve = 0

    for _rnd in range(1, max_rounds + 1):
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

    while extended:
        extended = False
        suffix = seq[-(k - 1) :]
        if suffix not in G or G.out_degree(suffix) == 0:
            break

        successors = list(G.successors(suffix))
        if len(successors) > 1:
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

    extended = True
    while extended:
        extended = False
        prefix = seq[: k - 1]
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


def get_hybrid_kmer_weights(
    df: pd.DataFrame,
    kmer_size: int,
) -> Counter:
    """
    Calculates k-mer weights using ONLY:
    1. Frequency (Implicit via accumulation)
    2. MS1 Abundance (Peptide Area)
    3. AI Confidence (Token Probabilities)
    """
    kmer_weights = Counter()

    col_tokens = "instanovo_token_log_probabilities"
    col_abundance = "peptide_abundance"

    for _, row in df.iterrows():
        sequence = row.get("cleaned_preds")
        if not isinstance(sequence, str) or len(sequence) < kmer_size:
            continue

        
        raw_probs_str = row.get(col_tokens)
        log_probs = [0.0] * len(sequence)  # Default neutral

        if isinstance(raw_probs_str, str) and raw_probs_str.startswith("["):
            try:
                parsed = ast.literal_eval(raw_probs_str)
                # Adjust for SOS/EOS tokens if present
                if len(parsed) == len(sequence) + 2:
                    log_probs = parsed[1:-1]
                elif len(parsed) == len(sequence):
                    log_probs = parsed
            except Exception:
                pass  # Keep default

        # --- 2. MS1 Abundance Score ---
        # Log scale: log10(Area + 10).
        abundance_val = row.get(col_abundance, 0)
        if pd.notnull(abundance_val) and abundance_val > 0:
            abundance_score = math.log10(float(abundance_val) + 10)
        else:
            abundance_score = 1.0

        # --- 3. Accumulate Weights ---
        for i in range(len(sequence) - kmer_size + 1):
            kmer = sequence[i : i + kmer_size]

            # Calculate local AI confidence for this specific k-mer
            sub_logs = log_probs[i : i + kmer_size] if i + kmer_size <= len(log_probs) else []
            if sub_logs:
                # Geometric mean of probabilities in the k-mer
                avg_log_prob = sum(sub_logs) / len(sub_logs)
                ai_score = math.exp(avg_log_prob)
            else:
                ai_score = 1.0

            # Final Weight: Abundance * AI Confidence
            # Frequency is handled implicitly because we += this value every time we see the k-mer
            weight = abundance_score * ai_score

            kmer_weights[kmer] += weight

    return kmer_weights


class Assembler:
    """
    Unified assembler supporting:
    - 'greedy': Overlap-Layout-Consensus style.
    - 'dbg': Standard De Bruijn Graph.
    - 'dbg_weighted': DBG with node/edge filtering and scoring.
    - 'dbgX': DBG with extension heuristics.
    - 'fusion': Hybrid DBG + Greedy.
    - 'multimodal': Heuristic DBG using MS1/MS2/AI/iRT features.
    """

    def __init__(
        self,
        mode: str = "greedy",
        min_overlap: int = 4,
        size_threshold: int = 10,
        kmer_size: int = 6,
        min_identity: float = 0.8,
        max_mismatches: int = 10,
        min_weight: int = 2,
        refine_rounds: int = 0,
        refine_patience: int = 2,
        alpha_len: float = 1.0,
        alpha_cov: float = 1.0,
        alpha_min: float = 0.2,
        reference_protein: str = None,
        stats_output_folder: str = None,
    ):
        if mode not in ["greedy", "dbg", "dbg_weighted", "dbgX", "fusion", "multimodal_dbg", "hybrid_dbg"]:
            raise ValueError(
                "mode must be 'greedy', 'dbg', 'dbg_weighted', 'dbgX', 'fusion', 'multimodal_dbg' or 'hybrid_dbg'"
            )

        self.mode = mode
        self.min_overlap = min_overlap
        self.size_threshold = size_threshold
        self.kmer_size = kmer_size
        self.min_identity = min_identity
        self.max_mismatches = max_mismatches
        self.min_weight = min_weight
        self.refine_rounds = refine_rounds
        self.refine_patience = refine_patience
        self.alpha_len = alpha_len
        self.alpha_cov = alpha_cov
        self.alpha_min = alpha_min
        self.reference_protein = reference_protein
        self.stats_output_folder = stats_output_folder

    def _compute_intermediate_stats(self, contigs, label):
        """Internal wrapper for statistics."""
        if self.reference_protein and self.stats_output_folder:
            logger.info(f"Computing intermediate statistics for {label}...")
            try:
                mapped = viz.process_protein_contigs_scaffold(
                    assembled_contigs=contigs,
                    target_protein=self.reference_protein,
                    max_mismatches=self.max_mismatches,
                    min_identity=self.min_identity
                )
                df_mapped = viz.create_dataframe_from_mapped_sequences(data=mapped)
                if not df_mapped.empty:
                    helpers.compute_assembly_statistics(
                        df=df_mapped,
                        sequence_type=label,
                        output_folder=self.stats_output_folder,
                        reference=self.reference_protein,
                    )
            except Exception as e:
                logger.warning(f"Could not compute intermediate stats: {e}")

    def assemble_greedy(self, sequences):
        logger.info(f"[Assembler] Running Greedy assembly (min_overlap={self.min_overlap})")
        contigs = assemble_contigs_greedy(sequences, self.min_overlap)
        contigs = list(set(contigs))
        contigs = sorted(contigs, key=len, reverse=True)
        self._compute_intermediate_stats(contigs, label="contig")

        scaffolds = scaffold_iterative_greedy(contigs, self.min_overlap, self.size_threshold)

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
        logger.info(f"[Assembler] Running DBG weighted (k={self.kmer_size}, min_weight={self.min_weight})")

        # 1. Standard Weighted DBG Assembly
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

        logger.info(f"DBG produced {len(contigs)} initial contigs.")

        self._compute_intermediate_stats(contigs, label="contig")

        # 2. OVERLAP GRAPH REFINEMENT (The new logic)
        # Only runs if refine_rounds is > 0
        if self.refine_rounds > 0:
            logger.info("Refining contigs using Overlap Graph (Bird's Eye View)...")

            # Use a slightly safer/larger overlap for this final merge to avoid false positives
            # e.g., max(min_overlap, 2) or just self.min_overlap
            #safe_overlap = max(self.min_overlap, 2)
            safe_overlap = self.min_overlap
            
            iteration = 0

            while iteration < self.refine_rounds:
                iteration += 1
                prev_count = len(contigs)

                # Core logic
                contigs = refine_using_overlap_graph(contigs, min_overlap=safe_overlap)

                new_count = len(contigs)

                # CONVERGENZA: Se non abbiamo fuso nulla, ci fermiamo.
                if new_count >= prev_count:
                    logger.info(f"Refinement converged at round {iteration}.")
                    break

                logger.info(f"Round {iteration}: reduced to {new_count} scaffolds.")

            if iteration >= self.refine_rounds:
                logger.warning(f"Refinement stopped hit max rounds limit ({self.refine_rounds}).")

        return contigs

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
        logger.info("[Assembler] Running FUSION (DBG weighted + greedy merge)")

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

    def assemble_multimodal_dbg(self, sequences: List[str], df_full: Optional[pd.DataFrame] = None) -> List[str]:
        """
        Multimodal Heuristic Assembly Strategy.

        Logic:
        1. Landscape Construction: Weights nodes by MS1 Abundance, MS2 Intensity, iRT, and AI Confidence.
        2. Seed Selection: Picks the 'Heaviest Seed' (highest confidence node).
        3. Smart Navigation: Extends forward/backward using Lookahead Score (Edge * Node).
        4. Path Burning: Removes assembled nodes to uncover lower-abundance variants.
        """
        logger.info(f"[Assembler] Running Multimodal DBG (Heuristic) k={self.kmer_size}")

        # --- 1. WEIGHT CALCULATION (NODES & EDGES) ---
        if df_full is not None:
            logger.info("Using Multimodal Features (Token Probs, MS1/MS2, iRT) for weighting.")

            # Nodes are (k-1)-mers
            node_weights = get_weighted_kmers_from_df(df_full, self.kmer_size - 1, use_abundance=True, use_quality=True)
            # Edges are k-mers
            edge_weights = get_weighted_kmers_from_df(df_full, self.kmer_size, use_abundance=True, use_quality=True)
            # Build graph using calculated edge weights
            G = build_dbg_from_kmers([], weights=edge_weights)
        else:
            # Fallback: Simple counts (if no dataframe provided)
            logger.warning("No DataFrame provided for Multimodal DBG. Falling back to raw counts.")
            node_kmers = get_kmers(sequences, self.kmer_size - 1)
            node_weights = Counter(node_kmers)

            edge_kmers = get_kmers(sequences, self.kmer_size)
            G = build_dbg_from_kmers(edge_kmers)

        # Apply weights to nodes
        nx.set_node_attributes(G, node_weights, name="weight")

        # --- 2. FILTERING ---
        if self.min_weight > 1:
            # Note: With multimodal weights (floats), min_weight might need tuning (e.g. 5.0)
            # but usually >1 filters out single-observation errors effectively.
            nodes_to_remove = [n for n, w in node_weights.items() if w < self.min_weight]
            G.remove_nodes_from(nodes_to_remove)

        contigs = []

        # --- 3. ASSEMBLY LOOP (Heaviest Seed + Smart Greedy) ---
        while G.number_of_nodes() > 0:
            # Find Seed: Node with highest weight
            try:
                seed_node = max(G.nodes, key=lambda n: G.nodes[n].get("weight", 0))
            except ValueError:
                break

            # Stop if seed is too weak
            if G.nodes[seed_node].get("weight", 0) < self.min_weight:
                break

            # Helper for Smart Greedy Decision
            def get_best_neighbor(current_node, direction="successors"):
                if direction == "successors":
                    neighbors = list(G.successors(current_node))
                else:
                    neighbors = list(G.predecessors(current_node))

                if not neighbors:
                    return None

                # Heuristic Score = EdgeWeight * DestNodeWeight
                def score(n):
                    if direction == "successors":
                        edge_w = G.get_edge_data(current_node, n).get("weight", 1)
                    else:
                        edge_w = G.get_edge_data(n, current_node).get("weight", 1)
                    node_w = G.nodes[n].get("weight", 1)
                    return edge_w * node_w

                return max(neighbors, key=score)

            # Extend Forward
            path_fwd = [seed_node]
            curr = seed_node
            while True:
                best_next = get_best_neighbor(curr, direction="successors")
                if not best_next or best_next in path_fwd:
                    break
                path_fwd.append(best_next)
                curr = best_next

            # Extend Backward
            path_bwd = []
            curr = seed_node
            while True:
                best_prev = get_best_neighbor(curr, direction="predecessors")
                if not best_prev or best_prev in path_bwd or best_prev in path_fwd:
                    break
                path_bwd.append(best_prev)
                curr = best_prev

            # Reconstruct
            full_path_nodes = path_bwd[::-1] + path_fwd

            if not full_path_nodes:
                G.remove_node(seed_node)
                continue

            sequence = full_path_nodes[0]
            for kmer in full_path_nodes[1:]:
                sequence += kmer[-1]

            if len(sequence) >= self.size_threshold:
                contigs.append(sequence)

            # Burn Path
            G.remove_nodes_from(full_path_nodes)

        return contigs

    def assemble_hybrid_dbg(self, sequences: List[str], df_full: pd.DataFrame) -> List[str]:
        """
        Streamlined Weighted DBG using only Abundance + AI Confidence + Frequency.
        Uses heuristic traversal (smart greedy) to resolve branches based on score.
        """
        logger.info(f"[Assembler] Running Hybrid DBG (Freq + Abundance + AI) k={self.kmer_size}")

        if df_full is None or df_full.empty:
            logger.warning("No DataFrame provided for Hybrid DBG. Falling back to standard DBG.")
            return self.assemble_dbg(sequences)

        # 1. Calculate Weights (Node & Edge)
        node_weights = get_hybrid_kmer_weights(df_full, self.kmer_size - 1)
        edge_weights = get_hybrid_kmer_weights(df_full, self.kmer_size)

        # 2. Build Graph
        G = build_dbg_from_kmers([], weights=edge_weights)
        nx.set_node_attributes(G, node_weights, name="weight")

        # 3. Filter Noise
        if self.min_weight > 1:
            to_remove = [n for n, w in node_weights.items() if w < self.min_weight]
            G.remove_nodes_from(to_remove)

        contigs = []

        # 4. Assembly Loop (Heuristic Traversal - Heaviest Path)
        while G.number_of_nodes() > 0:
            try:
                seed_node = max(G.nodes, key=lambda n: G.nodes[n].get("weight", 0))
            except ValueError:
                break

            if G.nodes[seed_node].get("weight", 0) < self.min_weight:
                break

            # Helper to pick best neighbor based on Edge * Node weight
            def get_best_neighbor(curr, direction="successors"):
                neighbors = list(G.successors(curr)) if direction == "successors" else list(G.predecessors(curr))
                if not neighbors:
                    return None

                def score(n):
                    edge_w = (
                        G.get_edge_data(curr, n)["weight"]
                        if direction == "successors"
                        else G.get_edge_data(n, curr)["weight"]
                    )
                    return edge_w * G.nodes[n].get("weight", 0)

                return max(neighbors, key=score)

            # Extend Forward
            path_fwd = [seed_node]
            curr = seed_node
            while True:
                nxt = get_best_neighbor(curr, "successors")
                if not nxt or nxt in path_fwd:
                    break
                path_fwd.append(nxt)
                curr = nxt

            # Extend Backward
            path_bwd = []
            curr = seed_node
            while True:
                prev = get_best_neighbor(curr, "predecessors")
                if not prev or prev in path_bwd or prev in path_fwd:
                    break
                path_bwd.append(prev)
                curr = prev

            # Reconstruct Sequence
            full_path = path_bwd[::-1] + path_fwd
            if not full_path:
                G.remove_node(seed_node)
                continue

            seq = full_path[0]
            for node in full_path[1:]:
                seq += node[-1]

            if len(seq) >= self.size_threshold:
                contigs.append(seq)

            # Burn Path
            G.remove_nodes_from(full_path)

        return contigs

    def run(self, sequences: List[str], df_full: Optional[pd.DataFrame] = None):
        """
        Main entry point.
        Args:
            sequences: List of peptide strings (required for basic modes).
            df_full: Optional DataFrame containing advanced features (required for optimized dbg_greedy).
        """
        if not sequences and (df_full is None or df_full.empty):
            logger.error("No valid input provided for assembly.")
            raise ValueError("Input sequences list or DataFrame is empty.")

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
        elif self.mode == "multimodal_dbg":
            return self.assemble_multimodal_dbg(sequences, df_full=df_full)
        elif self.mode == "hybrid_dbg":
            return self.assemble_hybrid_dbg(sequences, df_full=df_full)


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
    refine_rounds: int = 0,
):
    """Main function for standalone assembly."""
    output_path = Path(output_scaffolds_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    stats_folder = output_path.parent / "statistics"

    protein_norm = None  # None means no reference mode

    if reference:
        logger.info("Reference mode enabled. Loading reference protein...")
        if not metadata_json_path:
            raise ValueError("metadata_json_path is required.")
        try:
            path_obj = Path(input_csv_path)
            
            if path_obj.name == "cleaned.csv":
                run_name = path_obj.parent.parent.name
            else:
                run_name = path_obj.stem.replace("_cleaned", "")
    
            meta = helpers.get_sample_metadata(run=run_name, chain=chain, json_path=metadata_json_path)
            protein = meta["protein"]
            protein_norm = preprocessing.normalize_sequence(protein)
            logger.info("Reference protein loaded and normalized successfully.")
            stats_folder.mkdir(parents=True, exist_ok=True)

        except Exception as e:
            logger.error(f"Failed to get reference protein: {e}")
            logger.warning("Disabling reference mode.")
            reference = False  # if fails, disable reference mode

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
        max_mismatches=max_mismatches,
        refine_rounds=refine_rounds,
        reference_protein=protein_norm,
        stats_output_folder=str(stats_folder) if protein_norm else None
    )

    scaffolds = assembler.run(sequences=sequences, df_full=df)

    output_path = Path(output_scaffolds_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    # output_folder = output_path.parent.mkdir(parents=True, exist_ok=True)

    records = [
        Bio.SeqRecord.SeqRecord(Bio.Seq.Seq(seq), id=f"scaffold_{i + 1}", description=f"length: {len(seq)}")
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
        df_scaffolds_mapped = viz.create_dataframe_from_mapped_sequences(data=mapped_scaffolds)
        helpers.compute_assembly_statistics(
            df=df_scaffolds_mapped,
            sequence_type="scaffolds",
            output_folder=str(stats_path),
            reference=protein_norm,
        )
        logger.info(f"Reference mode: Statistics saved to {stats_path}")


def cli():
    """Command-line interface for the assembly module."""

    parser = argparse.ArgumentParser(description="Run Greedy or DBG assembly on peptide sequences.")
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
        default=None,
        help="Path to sample_metadata.json (required for --reference).",
    )
    parser.add_argument(
        "--assembly-mode",
        type=str,
        choices=["greedy", "dbg", "dbg_weighted", "dbgX", "fusion", "multimodal_dbg"],
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
        "--refine",
        action="store_true",
        help="Enables iterative refinement (Overlap Graph) until convergence.",
    )
    parser.add_argument(
        "--chain",
        type=str,
        default="",
        help="Specify chain type (light/heavy) required for reference lookup.",
    )
    parser.add_argument(
        "--min-identity",
        type=float,
        default=0.8,
        help="Minimum identity for reference mapping.",
    )
    parser.add_argument(
        "--max-mismatches",
        type=int,
        default=10,
        help="Maximum mismatches for reference mapping.",
    )

    args = parser.parse_args()

    refine_rounds_val = MAX_REFINE_ROUNDS if args.refine else 0

    # in case of non-DBG mode, ignore kmer_size
    if args.assembly_mode == "greedy":
        args.kmer_size = 0
        logger.info("Ignoring kmer_size (used only for DBG mode).")

    if args.reference and not args.metadata_json_path:
        parser.error("--metadata-json-path is required when --reference is enabled.")

    args_dict = vars(args)

    if "refine" in args_dict:
        del args_dict["refine"]

    main(refine_rounds=refine_rounds_val, **args_dict)


if __name__ == "__main__":
    cli()

# python -m instanexus.assembly --input-csv-path outputs/bsa/bsa_cleaned.csv --output-scaffolds-path outputs/bsa/bsa_scaffolds.fasta --metadata-json-path json/sample_metadata.json --assembly-mode dbg_weighted --kmer-size 7 --min-overlap 3 --size-threshold 10 --min-identity 0.8 --max-mismatches 10
