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
from collections import defaultdict
from itertools import combinations
from tqdm import tqdm
import pandas as pd
import networkx as nx


def get_kmers(seqs, kmer_size):
    """ Generate k-mers of specified length from a list of sequences; a k-mer is a substring of length `kmer_size` extracted from each input sequence.
    """
    
    kmers = []

    for seq in seqs:
        kmers.extend(seq[i:i+kmer_size] for i in range(len(seq) - kmer_size + 1))
    
    return kmers


def get_kmer_counts(kmers):
    """ Count occurrences of each k-mer in a list of k-mers; it takes a list of k-mers and returns a dictionary where the keys 
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
    """ Generate edges of a De Bruijn graph from a list of k-mers.
    
    A De Bruijn graph is a directed graph used in sequence assembly where each k-mer represents an edge, 
    and the nodes are (k-1)-mers. This function takes a list of k-mers and generates unique edges between 
    (k-1)-mers, avoiding duplicate edges using a set.
    """

    edges = set()
    k_1mers = defaultdict(set)
    
    for kmer in kmers:
        k_1mers[kmer[:-1]].add(kmer[1:])
    
    for prefix in k_1mers:
        for suffix in k_1mers[prefix]:
            edges.add((prefix, suffix))
    
    return edges


def assemble_contigs(edges):
    """ Assemble contigs from De Bruijn graph edges by traversing the graph; it takes a set of directed edges representing
    a De Bruijn graph and assembles contigs by performing a depth-first traversal. Each contig is a path in the graph where 
    each node is connected by an edge. The function uses an iterative approach to avoid recursion depth limits.
    """

    graph = defaultdict(list)
    for start, end in edges:
        graph[start].append(end)
    
    # Find starting nodes (nodes with no incoming edges)
    all_ends = set(e for _, e in edges)
    start_nodes = set(graph.keys()) - all_ends

    def traverse_iterative(start_node):
        """ Traverse a graph iteratively to find paths (contigs) starting from a given node.
        """
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


def get_kmers_from_df(df, kmer_size):
    """ Generate k-mers of specified length from a DataFrame, preserving metadata.
    """
    kmers_list = []
    
    for _, row in df.iterrows():
        sequence = row['preds']  
        for i in range(len(sequence) - kmer_size + 1):
            kmer = sequence[i:i+kmer_size]
            kmers_list.append({**row.to_dict(), 'kmer': kmer})
    
    return pd.DataFrame(kmers_list)


def find_overlaps(contigs, min_overlap, disable_tqdm=False):
    """Find overlaps between pairs of contigs based on specified minimum overlap.
    """
    overlaps = []
    total_pairs = sum(1 for _ in combinations(contigs, 2))  # Calculate total number of pairs

    with tqdm(total=total_pairs, desc="Finding overlaps", disable=disable_tqdm) as pbar:
        for a, b in combinations(contigs, 2): # combinations() generates all pairs of contigs
            for i in range(min_overlap, min(len(a), len(b)) + 1): # Check overlaps of different lengths
                if a[-i:] == b[:i]:
                    overlaps.append((a, b, i))
                if b[-i:] == a[:i]:
                    overlaps.append((b, a, i))
            pbar.update(1) 
    
    return overlaps


def create_scaffolds(contigs, min_overlap, disable_tqdm=False):
    """ Create scaffolds from a list of contigs by merging overlapping sequences.
    """

    overlaps = find_overlaps(contigs, min_overlap=min_overlap, disable_tqdm=disable_tqdm)
    combined_contigs = []
    for a, b, overlap in tqdm(overlaps, desc="Merging overlaps", total=len(overlaps), disable=disable_tqdm):
        combined = a + b[overlap:]
        combined_contigs.append(combined)

    return combined_contigs + contigs



def merge_sequences(contigs, disable_tqdm=False):
   """ Merges overlapping sequences. """
   contigs = sorted(contigs, key=len, reverse=True)
   merged = set(contigs)
   for c in tqdm(contigs, desc="Merging contigs", disable=disable_tqdm):
       for c2 in contigs:
           if c != c2 and c2 in c:  # if c2 is a substring of c
               merged.discard(c2)
   return list(merged)



def find_best_overlap(seq1, seq2, min_overlap):
    """ Identifies the longest overlap between two sequences, seq1 and seq2, with a minimum required overlap.
    """

    best_overlap_len = 0  # Length of the longest overlap found
    best_overlap_index = -1  # Starting index of the best overlap in seq1
    best_overlap_index2 = -1  # Starting index of the best overlap in seq2

    # Iterate through possible starting positions in seq2 for an overlap
    for i in range(len(seq2) - min_overlap + 1):
        
        # Find the starting index in seq1 where a potential overlap with seq2[i:i+min_overlap] begins
        ind = seq1.find(seq2[i:i+min_overlap])
        
        # Check if this substring exists in seq1
        if ind != -1:
            overlap_len = min_overlap  # Start with the minimum overlap length

            # Extend the overlap as long as characters continue to match
            for k in range(min_overlap, len(seq2) - i):
                # Stop if the overlap extends beyond the end of seq1 or characters no longer match
                if ind + k >= len(seq1) or seq1[ind + k] != seq2[i + k]:
                    break
                overlap_len += 1  # Increase overlap length for each matching character

            # Update best overlap if a longer one has been found
            if overlap_len > best_overlap_len:
                best_overlap_len = overlap_len
                best_overlap_index = ind  # Update starting index in seq1
                best_overlap_index2 = i  # Update starting index in seq2

    return best_overlap_len, best_overlap_index, best_overlap_index2



def recombine_sequences(sequences, min_overlap):
    """ Recombine a list of contigs (sequences) by identifying overlapping regions and generating new recombined sequences
    from these overlaps.
    """
    
    recombined_sequences = []

    # iterate over all unique pairs of contigs to check for possible overlaps
    for seq1, seq2 in tqdm(combinations(sequences, 2), desc="Recombining contigs"):
        
        # Find the best overlap between seq1 and seq2
        overlap_len, overlap_index1, overlap_index2 = find_best_overlap(seq1, seq2, min_overlap)
        
        # Add the original sequences to the list of recombined sequences
        recombined_sequences.append(seq1)
        recombined_sequences.append(seq2)
        
        # If an overlap was found, proceed with recombination
        if overlap_len != -1:
            # Extract the overlapping segment from seq1 based on the indices
            overlap = seq1[overlap_index1:overlap_index1 + overlap_len]
            
            # Split the sequences into overlapping and non-overlapping regions
            nterm_seq1 = seq1[:overlap_index1]        # N-terminal (start) of seq1 before overlap
            nterm_seq2 = seq2[:overlap_index2]        # N-terminal (start) of seq2 before overlap
            cterm_seq1 = seq1[overlap_index1 + overlap_len:]  # C-terminal (end) of seq1 after overlap
            cterm_seq2 = seq2[overlap_index2 + overlap_len:]  # C-terminal (end) of seq2 after overlap

            # Create new recombined sequences using the non-overlapping ends and the stable overlap
            recombined_sequences.append(nterm_seq1 + overlap + cterm_seq2)  # Seq1 start + overlap + Seq2 end
            recombined_sequences.append(nterm_seq2 + overlap + cterm_seq1)  # Seq2 start + overlap + Seq1 end

    return recombined_sequences


def calculate_overlap(s1, s2, min_overlap):
    max_olap = 0
    max_len = min(len(s1), len(s2))
    for k in range(min_overlap, max_len + 1):
        if s1[-k:] == s2[:k]:
            max_olap = k
    return max_olap

def remove_contained_sequences(sequences):
    unique_seqs = list(set(sequences))
    to_remove = set()
    for i, s1 in enumerate(unique_seqs):
        for j, s2 in enumerate(unique_seqs):
            if i != j and s1 in s2:
                to_remove.add(s1)
    return [s for s in unique_seqs if s not in to_remove]

def build_overlap_graph(sequences, min_overlap):
    G = nx.DiGraph()
    for seq in sequences:
        G.add_node(seq)
    for s1 in sequences:
        for s2 in sequences:
            if s1 != s2:
                olap = calculate_overlap(s1, s2, min_overlap)
                if olap >= min_overlap:
                    G.add_edge(s1, s2, weight=olap)
    return G

def merge_path_sequences(path, min_overlap):
    scaffold = path[0]
    for i in range(1, len(path)):
        olap = calculate_overlap(path[i-1], path[i], min_overlap)
        scaffold += path[i][olap:]
    return scaffold

def find_all_paths_networkx_include_isolated(graph):
    all_paths = []
    nodes = list(graph.nodes())
    included_nodes = set()
    for source in nodes:
        for target in nodes:
            if source != target:
                try:
                    paths = list(nx.all_simple_paths(graph, source=source, target=target))
                    all_paths.extend(paths)
                    for p in paths:
                        included_nodes.update(p)
                except nx.NetworkXNoPath:
                    continue
    isolated_nodes = [n for n in nodes if n not in included_nodes]
    for n in isolated_nodes:
        all_paths.append([n])
    return all_paths


def dfs_custom(graph, start, path=None, visited=None):
    if path is None:
        path = [start]
    if visited is None:
        visited = set([start])
    paths = [path]
    for neighbor in graph.successors(start):
        if neighbor not in visited:
            new_paths = dfs_custom(graph, neighbor, path + [neighbor], visited | {neighbor})
            paths.extend(new_paths)
    return paths

def find_all_paths_dfs(graph):
    all_paths = []
    for node in graph.nodes():
        all_paths.extend(dfs_custom(graph, node))
    return all_paths


def assemble_scaffolds_networkx(sequences, min_overlap):
    """ Assemble scaffolds using NetworkX to find all paths in the overlap graph."""

    filtered_seqs = remove_contained_sequences(sequences)

    graph = build_overlap_graph(filtered_seqs, min_overlap)
    
    paths_nx = find_all_paths_networkx_include_isolated(graph)

    scaffolds_nx = [merge_path_sequences(path, min_overlap) for path in paths_nx]
    
    return scaffolds_nx


def assemble_scaffolds_dfs(sequences, min_overlap):
    """ Assemble scaffolds using a custom DFS approach to find all paths in the overlap graph."""
    
    filtered_seqs = remove_contained_sequences(sequences)
    
    graph = build_overlap_graph(filtered_seqs, min_overlap)
    
    paths_dfs = find_all_paths_dfs(graph)
    
    scaffolds_dfs = [merge_path_sequences(path, min_overlap) for path in paths_dfs]
    
    return scaffolds_dfs


def scaffold_iterative(contigs, min_overlap, size_threshold, disable_tqdm=False):
    prev = None
    current = contigs
    while prev != current:
        prev = current
        current = create_scaffolds(current, min_overlap, disable_tqdm)
        current = merge_sequences(current, disable_tqdm)

        current = list(set(current))
        current = [s for s in current if len(s) > size_threshold]
    return sorted(current, key=len, reverse=True)