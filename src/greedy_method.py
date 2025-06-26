#!/usr/bin/env python

r""" Greedy assembly script.
 _____  _______  _    _ 
|  __ \|__   __|| |  | |
| |  | |  | |   | |  | |
| |  | |  | |   | |  | |
| |__| |  | |   | |__| |
|_____/   |_|   |______|

__authors__ = Marco Reverenna & Konstantinos Kalogeropoulus
__copyright__ = Copyright 2024-2025
__research-group__ = DTU Biosustain (Multi-omics Network Analytics) and DTU Bioengineering
__date__ = 25 Feb 2025
__maintainer__ = Marco Reverenna
__email__ = marcor@dtu.dk
__status__ = Dev
"""

import Bio
from itertools import combinations
from collections import defaultdict
import numpy as np
from tqdm import tqdm
import Bio.Seq
import Bio.SeqRecord



# def find_overlaps_greedy(peptides, min_overlap):
#     overlaps = defaultdict(list) # dict to store overlaps
#     for i, pep1 in tqdm(enumerate(peptides), desc="Finding overlaps"): # iterate over peptides
#         for j, pep2 in enumerate(peptides): # iterate over peptides
#             if i != j: # if peptides are not the same
#                 # avoid the compare the same peptides twice
#                 for k in range(min_overlap, min(len(pep1), len(pep2))): # 
#                     if pep1[-k:] == pep2[:k]:
#                         overlaps[i].append((j, k))
#                     if pep2[-k:] == pep1[:k]:
#                         overlaps[j].append((i, k))
#     return overlaps


def find_peptide_overlaps(peptides, min_overlap):
    """Finds overlaps between peptide sequences using a greedy approach.
    
    Args:
        peptides (list of str): A list of peptide sequences.
        min_overlap (int): The minimum length required for an overlap, but 
                           the overlap can also be done with 3 amino acids.
        
    Returns:
        dict: A dictionary where the keys are the indices of the peptides and 
              the values are lists of tuples. Each tuple contains the index of 
              the overlapping peptide and the length of the overlap.

    Example:
        >>> peptides = ['AABB', 'AABB']
        >>> find_peptide_overlaps(peptides, 2)
        {0: [(1, 2)], 1: [(0, 2)]}
        # peptide at index 0 (first one) has an overlap 

    """
    overlaps = defaultdict(list)  # Dictionary to store overlaps

    for index_a, peptide_a in tqdm(enumerate(peptides), desc="Finding overlaps"):
        # index_a is the index of the first peptide
        # peptide_a is the sequence of the first peptide
        for index_b, peptide_b in enumerate(peptides):
            # index_b is the index of the second peptide
            # peptide_b is the sequence of the second peptide

            if index_a != index_b:  # Skip comparing the same peptide

                max_possible_overlap = min(len(peptide_a), len(peptide_b))
                
                # Check all possible overlap lengths starting from min_overlap
                for overlap_length in range(min_overlap, max_possible_overlap):
                    # If the suffix of peptide_a matches the prefix of peptide_b
                    if peptide_a[-overlap_length:] == peptide_b[:overlap_length]:
                        # peptide_a[-overlap_length:] is the suffix of peptide_a (ABCD -> CD)
                        # peptide_b[:overlap_length] is the prefix of peptide_b (CDEF -> CD)
                    
                        overlaps[index_a].append((index_b, overlap_length)) # Add the overlap to the dictionary

                    # If the suffix of peptide_b matches the prefix of peptide_a
                    if peptide_b[-overlap_length:] == peptide_a[:overlap_length]:
                        # peptide_b[-overlap_length:] is the suffix of peptide_b (CDEF -> EF)
                        # peptide_a[:overlap_length] is the prefix of peptide_a (ABCD -> AB)

                        overlaps[index_b].append((index_a, overlap_length)) # Add the overlap to the dictionary

    return overlaps



def assemble_contigs(peptides, min_overlap):
    """Assemble contigs from a list of peptides using a greedy approach."""

    assembled_contigs = peptides[:] # copy the list of peptides
    iteration = 0
    
    while True: # loop until no more overlaps are found
        iteration += 1
        #print(f"Iteration {iteration}:")

        overlaps = find_peptide_overlaps(assembled_contigs, min_overlap) # find overlaps between contigs
        
        if not overlaps: # if no overlaps are found,
            break
        new_contigs = set() # store new contigs
        used_indices = set() # store indices of used contigs
        
        for i, overlaps_list in tqdm(overlaps.items(), desc="Processing overlaps"):
            if i in used_indices:
                continue
            best_match = max(overlaps_list, key=lambda x: x[1], default=None)
            if best_match:
                j, overlap_len = best_match
                if j not in used_indices:
                    new_contig = assembled_contigs[i] + assembled_contigs[j][overlap_len:]
                    new_contigs.add(new_contig)
                    used_indices.update([i, j])
        
        # add contigs that were not merged
        remaining_contigs = [contig for idx, contig in enumerate(assembled_contigs) if idx not in used_indices] 
        # remaining_contigs are peptides which are not merged with any other peptide
        
        # update assembled contigs
        assembled_contigs = list(new_contigs) + remaining_contigs
        
        if len(new_contigs) == 0:
            #print(f"No more overlaps found after {iteration} iterations. Assembly complete.")
            break

    return assembled_contigs



def find_contig_overlap(seq1, seq2, min_overlap):
    """ 
    Finds the maximum overlap between two sequences with a minimum overlap length.
    
    Parameters:
        seq1 (str): The first sequence.
        seq2 (str): The second sequence.
        min_overlap (int): The minimum required overlap length.
        
    Returns:
        tuple: (smallest_seq, largest_seq, max_overlap_len, max_overlap_pos_small, max_overlap_pos_large)
               where smallest_seq and largest_seq are the shorter and longer sequences respectively,
               max_overlap_len is the length of the maximum overlap found,
               max_overlap_pos_small is the starting position of the overlap in the smaller sequence,
               max_overlap_pos_large is the starting position of the overlap in the larger sequence.
        None: If no overlap of at least min_overlap is found.
    """
    max_overlap_len = 0  # Maximum overlap length found so far
    max_overlap_pos_small = None  # Starting position of the overlap in the smaller sequence
    max_overlap_pos_large = None  # Starting position of the overlap in the larger sequence
    
    # Determine which sequence is smaller and which is larger.
    # In case of equal lengths, seq2 is chosen as the smaller sequence.
    smallest_seq = seq1 if len(seq1) < len(seq2) else seq2
    largest_seq = seq1 if len(seq1) >= len(seq2) else seq2
    
    # Iterate over each possible starting position in the smaller sequence.
    for i in range(len(smallest_seq)):
        # Iterate over each possible starting position in the larger sequence.
        for k in range(len(largest_seq)):
            # Check if the substrings of length min_overlap starting at positions i and k are equal.
            if smallest_seq[i:i + min_overlap] == largest_seq[k:k + min_overlap]:
                overlap_len = min_overlap
                # Extend the overlap while the next characters match in both sequences.
                while (i + overlap_len < len(smallest_seq) and 
                       k + overlap_len < len(largest_seq) and 
                       smallest_seq[i + overlap_len] == largest_seq[k + overlap_len]):
                    overlap_len += 1
                
                # Update the maximum overlap if the current one is longer.
                if overlap_len > max_overlap_len:
                    max_overlap_len = overlap_len
                    max_overlap_pos_small = i
                    max_overlap_pos_large = k

    # If a valid overlap is found that meets the minimum requirement, return the results.
    if max_overlap_len >= min_overlap:
        return smallest_seq, largest_seq, max_overlap_len, max_overlap_pos_small, max_overlap_pos_large
    else:
        # Otherwise, return None indicating no sufficient overlap was found.
        return None



def combine_sequences(seq1, seq2, min_overlap):
    """ Combine two sequences based on a minimum overlap.
    This function attempts to find an overlap between two sequences and combines them in various ways
    based on the overlap. If an overlap is found, it returns a list of possible combinations of the sequences.
    If no overlap is found, it returns None.
    """
    overlap = find_contig_overlap(seq1, seq2, min_overlap)

    if overlap:
        smallest_seq, largest_seq, overlap_len, overlap_pos_small, overlap_pos_large = overlap
        nterm_smallest = smallest_seq[:overlap_pos_small]
        cterm_smallest = smallest_seq[overlap_pos_small + overlap_len:]
        nterm_largest = largest_seq[:overlap_pos_large]
        cterm_largest = largest_seq[overlap_pos_large + overlap_len:]
        seq_overlap = smallest_seq[overlap_pos_small:overlap_pos_small + overlap_len]

        combinations = [
            nterm_smallest + seq_overlap + cterm_largest,
            nterm_largest + seq_overlap + cterm_smallest,
            smallest_seq,
            largest_seq
        ]
        return combinations
    else:
        return None



def recombine_contigs(contigs, min_overlap):
    """ Recombines a list of contigs by finding overlaps and merging them.
    """
    combined = []

    for i, contig1 in tqdm(enumerate(contigs), desc="Combining contigs"):
        for j, contig2 in enumerate(contigs):
            if i != j:
                combinations = combine_sequences(contig1, contig2, min_overlap)
                if combinations:
                    combined.extend(combinations)

    return combined


def merge_contigs(contigs):
   """ Merges overlapping contigs into a set of unique contigs.
   Takes a list of contigs and merges them by checking if one contig is a 
   substring of another. If a contig is found to be a substring of another,
   it is discarded, and the larger contig is kept.
   """
   contigs = sorted(contigs, key=len, reverse=True)
   merged = set(contigs)
   for c in tqdm(contigs, desc="Merging contigs"):
       #print(c)
       for c2 in contigs:
           if c != c2 and c2 in c: # if c is a substring of c2
               merged.discard(c2)

   return list(merged)



def find_overlaps(contigs, min_overlap):
    """ Find overlaps between pairs of contigs.
    This function takes a list of contigs and a minimum overlap length, and 
    returns a list of tuples representing the overlaps between pairs of contigs. 
    Each tuple contains two contigs and the length of their overlap.
    """
    overlaps = []

    for a, b in combinations(contigs, 2):
        for i in range(min_overlap, min(len(a), len(b)) + 1):
            if a[-i:] == b[:i]:
                overlaps.append((a, b, i))
            if b[-i:] == a[:i]:
                overlaps.append((b, a, i))

    return overlaps



def combine_seqs_into_scaffolds(contigs, min_overlap):
    """ Combine contigs based on a minimum overlap length.
    This function takes a list of contigs and a minimum overlap length, finds
    the overlaps between the contigs, and combines them. The combined
    contigs are then returned along with the original contigs.
    """
    overlaps = find_overlaps(contigs, min_overlap=min_overlap)
    combined_contigs = []

    for a, b, overlap in overlaps:
        combined = a + b[overlap:]
        combined_contigs.append(combined)

    return combined_contigs + contigs