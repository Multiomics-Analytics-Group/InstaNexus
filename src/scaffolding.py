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
from tqdm import tqdm
from itertools import combinations


def find_overlaps(contigs, min_overlap):

    """
    Find overlaps between pairs of contigs based on specified minimum overlap.

    This function takes a list of contigs and identifies overlapping
    regions between pairs of contigs where the overlap is at least `min_overlap` 
    nucleotides. For each pair of contigs that overlap, it records the contigs and 
    the length of the overlap.

    Parameters:
    - contigs (list of str): A list of contig sequences (strings) to check for overlaps.
    - min_overlap (int): The minimum number of nucleotides that must overlap between two contigs.

    Returns:
    - overlaps (list of tuples): A list of tuples, where each tuple contains:
        - contig1 (str): The first contig in the overlapping pair.
        - contig2 (str): The second contig in the overlapping pair.
        - overlap_length (int): The length of the overlap between the two contigs.

    Example:
    ```
    contigs = ["ATCG", "CGTA", "GTAC", "TACG"]
    min_overlap = 2
    overlaps = find_overlaps(contigs, min_overlap)
    # Output: [('ATCG', 'CGTA', 2), ('CGTA', 'TACG', 3), ('GTAC', 'TACG', 3)]
    ```
    """
    overlaps = []
    total_pairs = sum(1 for _ in combinations(contigs, 2))  # Calculate total number of pairs

    with tqdm(total=total_pairs, desc="Finding overlaps") as pbar:
        for a, b in combinations(contigs, 2): # combinations() generates all pairs of contigs
            for i in range(min_overlap, min(len(a), len(b)) + 1): # Check overlaps of different lengths
                if a[-i:] == b[:i]:
                    overlaps.append((a, b, i))
                if b[-i:] == a[:i]:
                    overlaps.append((b, a, i))
            pbar.update(1)  # Update progress bar after each pair is processed
    
    return overlaps


# def filter_contained_sequences(sequences):
#     """
#     Filters out sequences that are fully contained within other sequences in a list, retaining only the most complete or longest versions.

#     This function iterates through a list of sequences (e.g., contigs or scaffolds) and identifies any sequence that is entirely contained within another sequence.
#     It removes the contained sequence from the result to leave only unique, non-redundant sequences.

#     Parameters:
#     - sequences (list of str): A list of sequence strings (e.g., contigs or scaffolds) to be filtered.

#     Returns:
#     - list of str: A list of sequences with fully contained, redundant sequences removed.

#     Notes:
#     - Useful for both contig merging and scaffolding, as it eliminates sequences that are subsequences of larger ones.

#     """

#     merged = set(sequences)  # Use a set to handle unique sequences
    
#     # Compare each sequence with every other to remove redundant contained sequences
#     for seq in tqdm(sequences, desc="Filtering contained sequences"):
#         for other_seq in sequences:
#             # Check if `seq` is fully contained within `other_seq`
#             if seq != other_seq and seq in other_seq:
#                 merged.discard(seq)  # Remove `seq` if it's contained in a larger sequence
#                 merged.add(other_seq)  # Ensure the larger sequence remains in the set
#                 break

#     return list(merged)


def create_scaffolds(contigs, min_overlap):

    """
    Create scaffolds from a list of contigs by merging overlapping sequences.

    This function finds overlaps between contigs based on a specified minimum overlap 
    and merges them to form longer sequences (scaffolds). Overlapping regions are 
    combined, and the resulting merged contigs are added to the final list of scaffolds.

    Parameters:
    - contigs (list of str): A list of contig sequences to be merged based on overlaps.
    - min_overlap (int): The minimum length of overlap required to merge two contigs.

    Returns:
    - combined_contigs (list of str): A list of merged contig sequences (scaffolds).
    """

    overlaps = find_overlaps(contigs, min_overlap=min_overlap)
    combined_contigs = []
    for a, b, overlap in tqdm(overlaps, desc="Merging overlaps", total=len(overlaps)):
        combined = a + b[overlap:]
        combined_contigs.append(combined)
    return combined_contigs + contigs


def filter_and_sort_scaffolds(scaffolds, size_threshold):
    """
    Filters and sorts a list of scaffolds.
    
    Parameters:
    - scaffolds (list of str): List of scaffolds (sequences).
    - size_threshold (int): Minimum length of scaffolds to keep.
    
    Returns:
    - list of str: Filtered and sorted scaffolds.
    """
    # Remove duplicates
    scaffolds = list(set(scaffolds))
    
    # Filter scaffolds based on size threshold
    filtered_scaffolds = [contig for contig in scaffolds if len(contig) > size_threshold]
    
    # Sort scaffolds by length in descending order
    return sorted(filtered_scaffolds, key=len, reverse=True)



def find_best_overlap(seq1, seq2, min_overlap):
    """
    Identifies the longest overlap between two sequences, seq1 and seq2, with a minimum required overlap.

    This function goes through seq2 to find any substring of length `min_overlap` or more 
    that overlaps with a portion of seq1. It returns the length of the best overlap found, 
    along with the starting indices of the overlap in both sequences.

    Parameters:
    - seq1 (str): The first sequence to compare.
    - seq2 (str): The second sequence to compare.
    - min_overlap (int): The minimum required length for an overlap to be considered valid.

    Returns:
    - tuple: A tuple containing:
        - best_overlap_len (int): Length of the longest overlap found. If no overlap, returns 0.
        - best_overlap_index (int): Starting index of the overlap in seq1. If no overlap, returns -1.
        - best_overlap_index2 (int): Starting index of the overlap in seq2. If no overlap, returns -1.
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

    # Return the details of the best overlap found
    return best_overlap_len, best_overlap_index, best_overlap_index2


# dbg makes mistakes in ends
# min_overlap in this case is the minimum length required between two sequences to be considered overlapping

def recombine_sequences(sequences, min_overlap):
    """ Recombine a list of contigs (sequences) by identifying overlapping regions and generating 
    new recombined sequences from these overlaps.

    This function goes through all pairs of contigs and attempts to find overlaps of at least 
    `min_overlap` length. If an overlap is found, the overlapping region is retained, while 
    the non-overlapping ends are recombined to create new sequence variations.
    """
    
    recombined_sequences = []  # List to store original and recombined sequences

    # Iterate over all unique pairs of contigs to check for possible overlaps
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

    # Return the list of recombined sequences
    return recombined_sequences