#!/usr/bin/env python

r""" Full assembly script for proteins.
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

# my modules
import greedy_method as greedy
import dbg
import mapping as map
import consensus as cons
import alignment as align
import clustering as clus
import preprocessing as prep
import compute_statistics as comp_stat

# import libraries
from tqdm import tqdm
from tempfile import mkdtemp
from itertools import combinations
from collections import defaultdict, Counter
from Bio import SeqIO

import sys
import os
import json
import re
import Bio
import shutil
import logging
import importlib
import statistics
import subprocess
import numpy as np
import pandas as pd
import seaborn as sns
import networkx as nx
import plotly.express as px
import matplotlib.pyplot as plt
import plotly.graph_objects as go

def run_pipeline_dbg(conf, kmer_size, min_overlap, max_mismatches, min_identity, size_threshold):

    ass_method = 'dbg'

    # run = "bsa"
    # protein = 'MKWVTFISLLLLFSSAYSRGVFRRDTHKSEIAHRFKDLGEEHFKGLVLIAFSQYLQQCPFDEHVKLVNELTEFAKTCVADESHAGCEKSLHTLFGDELCKVASLRETYGDMADCCEKQEPERNECFLSHKDDSPDLPKLKPDPNTLCDEFKADEKKFWGKYLYEIARRHPYFYAPELLYYANKYNGVFQECCQAEDKGACLLPKIETMREKVLASSARQRLRCASIQKFGERALKAWSVARLSQKFPKAEFVEVTKLVTDLTKVHKECCHGDLLECADDRADLAKYICDNQDTISSKLKECCDKPLLEKSHCIAEVEKDAIPENLPPLTADFAEDKDVCKNYQEAKDAFLGSFLYEYSRRHPEYAVSVLLRLAKEYEATLEECCAKDDPHACYSTVFDKLKHLVDEPQNLIKQNCDQFEKLGEYGFQNALIVRYTRKVPQVSTPTLVEVSRSLGKVGTRCCTKPESERMPCTEDYLSLILNRLCVLHEKTPVSEKVTKCCTESLVNRRPCFSALTPDETYVPKAFDEKLFTFHADICTLPDTEKQIKKQTALVELLKHKPKATEEQLKTVMENFVAFVDKCCAADDKEACFAVEGPKLVVSTQTALA'
    # chain = ''
    # proteases = ['Chymotrypsin', 'Legumain', 'Krakatoa', 'Elastase', 'Trypsin', 'Papain', 'Thermo', 'ProtK', 'GluC', 'LysC']

    ##################
    ### ANTOBODIES ###
    ##################

    #run = "ma1"
    #protein = "EVQLVQSGAEVKKPGASVKVSCKASGYTFTSYGLSWVRQAPGQGLEWMGWLSAYNGNTNYAQKLQGRVTMTTDTSTSTAYMELRSLRSDDTAVYYCAREEYCSSTSCYVGNYYGMDVWGQGTLVTVSSASTKGPSVFPLAPSSKSTSGGTAALGCLVKDYFPEPVTVSWNSGALTSGVHTFPAVLQSSGLYSLSSVVTVPSSSLGTQTYLCNVNHKPSNTKVDKKVEPKSCDKTHTCPPCPAPEAAGGPSVFLFPPKPKDTLYLTREPEVTCVVVDVSHEDPEVKFNWYVDGVEVHNAKTKPREEQYNSTYRVVSVLTVLHQDWLNGKEYKCKVSNKALPAPLEKTLSKAKGQPREPQVYTLPPSRDELTKNQVSLTCLVKGFYPSDLAVEWESNGQPENNYKTTPPVLDSDGSFFLYSKLTVDKSRWQQGNVFSCSVMHEALHNHYTQKSLSLSPGK"
    #chain = "heavy"
    #proteases = ["Thermo", "Papain", "Chemo", "Trypsin", "Elastase", "ProtK", "GluC"]

    #run = "ma2"
    #protein = "EVQLVQSGAEVKKPGSSVKVSCKASGGTFSSYALSWVRQAPGQGLEWMGGLLPLFGTANYAQKFQGRVTLTADESTSTAYMELRSLRSDDTAVYYCARDNLGYCSGGSCYSDYYYYYMDVWGQGTLVTVSSASTKGPSVFPLAPSSKSTSGGTAALGCLVKDYFPEPVTVSWNSGALTSGVHTFPAVLQSSGLYSLSSVVTVPSSSLGTQTYLCNVNHKPSNTKVDKKVEPKSCDKTHTCPPCPAPELLGGPSVFLFPPKPKDTLMLSRTPEVTCVVVDVSHEDPEVKFNWYVDGVEVHNAKTKPREEQYNSTYRVVSVLTVLHQDWLNGKEYKCKVSNKALPAPLEKTLSKAKGQPREPQVYTLPPSRDELTKNQVSLTCLVKGFYPSDLAVEWESNGQPENNYKTTPPVLDSDGSFFLYSKLTVDKSRWQQGNVFSCSVMHEALHNHYTQKSLSLSPGK"
    #chain = "heavy"
    #proteases = ["Thermo", "Papain", "Chemo", "Trypsin", "Elastase", "ProtK", "GluC"]

    #run = "ma3"
    #protein = "EVQLVQSGAEVKKPGSSVKVSCKASGGTFSSYALSWVRQAPGQGLEWMGGLLPLFGTANYAQKFQGRVTLTADESTSTAYMELRSLRSDDTAVYYCARDNLGYCSGGSCYSDYYYYYMDVWGQGTLVTVSSASTKGPSVFPLAPSSKSTSGGTAALGCLVKDYFPEPVTVSWNSGALTSGVHTFPAVLQSSGLYSLSSVVTVPSSSLGTQTYLCNVNHKPSNTKVDKKVEPKSCDKTHTCPPCPAPEAAGGPSVFLFPPKPKDTLMLSRTPEVTCVVVDVSHEDPEVKFNWYVDGVEVHNAKTKPREEQYNSTYRVVSVLTVLHQDWLNGKEYKCKVSNKALPAPLEKTLSKAKGQPREPQVYTLPPSRDELTKNQVSLTCLVKGFYPSDLAVEWESNGQPENNYKTTPPVLDSDGSFFLYSKLTVDKSRWQQGNVFSCSVMHEALHNHYTQKSLSLSPGK"
    #chain = "heavy"
    #proteases = ["Thermo", "Papain", "Chemo", "Trypsin", "Elastase", "ProtK", "GluC"]

    #run = "ma1"
    #protein = "DLVMTQSPSSLSASVGDRVTLTCQASQDLSNYLNWYQQKPGKAPKLLLYAASSLESGVPSRFSGSGSGTEFTLTLTNLQVDDFATYYCQRYDSNFAFGQGTKVELKRTVAAPSVFLFPPSDEQLKSGTASVVCLLNNFYPREAKVQWKVDNALQSGNSQESVTEQDSKDSTYSLSSTLTLSKADYEKHKVYACEVTHQGLSSPVTKSFNRGEC"
    #chain = "light"
    #proteases = ["Thermo", "Papain", "Chemo", "Trypsin", "Elastase", "ProtK", "GluC"]

    #run = "ma2"
    #protein = "NFMLTQPRSVSESPGKTVTLSCTRSSGSLGSDYVHWYQQRPGSSPTTVLYEDNQRPSGVPDRFSGSLDSSSNSASLTLSGLKTEDEADYYCQSYDRSNHEVVFGGGTKLTVLGQPKAAPSVTLFPPSSEELQANKATLVCLLSDFYPGAVTVAWKADSSPVKAGVETTTPSKQSNNKYAASSYLSLTPEQWKSHRSYSCQVTHEGSTVEKTVAPTECS"
    #chain = "light"
    #proteases = ["Thermo", "Papain", "Chemo", "Trypsin", "Elastase", "ProtK", "GluC"]

    #run = "ma3"
    #protein = "NFMLTQPPSVSVAPGRTATLTCEGDNLGQQLVHWYQQKPGQAPVAVLSSDSDRPSGLPERFSGSNSGNTATLTLSRVEAGDEADYYCQVWDSGSDHVVFGGGTKLTVLGQPKAAPSVTLFPPSSEELQANKATLVCLLSDFYPGAVTVAWKADSSPVKAGVETTTPSKQSNNKYAASSYLSLTPEQWKSHRSYSCQVTHEGSTVEKTVAPTECS"
    #chain = "light"
    #proteases = ["Thermo", "Papain", "Chemo", "Trypsin", "Elastase", "ProtK", "GluC"]

    #run = "NB1"
    #protein = 'QVQLQESGGGLVQPGGSLRLSCAASGSIVNINYMRWYRQAPGKQRELVAVITSAGNTNYAESVKGRFTISIDNAKKMVFLQMNSLKPDDTAVYYCHADLRVRDGVRGDYWGQGTQVTVSSAAADYKDHDGDYKDHDIDYKDDDDKGAAHHHHHH'
    #proteases = ["Vesuvius", "Krakatoa", "Elastase", "Trypsin", "GluC", "Chymotrypsin", "Papain", "ProteinaseK", "Thermolysin"]
    #chain = ""

    #run = "NB2"
    #protein = 'QVQLQESGGGLVQPGGSLRLSCAASGFTVSSVTLSWLRQAPGKGLEWVSDITSNGQTYYADSVKGRFTISRDNAKNTIYLQMNSLKADDSAVYFCAEDRWRSSNHPRGQGTQVTVSSAAADYKDHDGDYKDHDIDYKDDDDKGAAHHHHHH'
    #proteases = ["Vesuvius", "Krakatoa", "Elastase", "Trypsin", "GluC", "Chymotrypsin", "Papain", "ProteinaseK", "Thermolysin"]
    #chain = ""

    #run = "NB3"
    #protein = 'QVQLQESGGGLVQAGGSLRLSCLASGRTFSDYRIGWFRQAPGKEREFVSTIRNDDANTYYADSVKGRFTISRDNAKNTVYLQMNSLKPEDTAVYYCAAGARHTAQTMAAGKGIDYWGQGTQVTVSSAAADYKDHDGDYKDHDIDYKDDDDKGAAHHHHHH'
    #proteases = ["Vesuvius", "Krakatoa", "Elastase", "Trypsin", "GluC", "Chymotrypsin", "Papain", "ProteinaseK", "Thermolysin"]
    #chain = ""

    #run = "NB4"
    #protein = 'QVQLQESGGGLVEPGGSLRLSCAASGFTFINYRMSWVRQAPGKGLEWVSGINPDGGTSYSDSVKGRFTISRDNAKNTLYLQMNSLKVEDTAVYYCIQSGTSRRGQGTQVTVSSAAADYKDHDGDYKDHDIDYKDDDDKGAAHHHHHH'
    #proteases = ["Vesuvius", "Krakatoa", "Elastase", "Trypsin", "GluC", "Chymotrypsin", "Papain", "ProteinaseK", "Thermolysin"]
    #chain = ""

    #run = "NB5"
    #protein = 'QVQLQESGGGLVQPGGSLRLSCAASGNIFSINYMKWYRQAPGKQRELVAVITDGGRTNYADSVKGRFAISRDNAKNTTYLQMSDLQPEDTAVYYCYADLRVVDGRHLPRGDYWGQGTQVTVSSAAADYKDHDGDYKDHDIDYKDDDDKGAAHHHHHH'
    #proteases = ["Vesuvius", "Krakatoa", "Elastase", "Trypsin", "GluC", "Chymotrypsin", "Papain", "ProteinaseK", "Thermolysin"]
    #chain = ""

    #run = "NB6"
    #protein = 'QVQLQESGGGLVQPGGSLRLSCTASLNIFSINAMGWYRQAPGKQRELVAAITSGGSTNYADSVKGRFTISRDNAKSTVYLQMNSLKPEDTAVYYCHAEGPFNIATKEQYDYWGQGTQVTVSSAAADYKDHDGDYKDHDIDYKDDDDKGAAHHHHHH'
    #proteases = ["Vesuvius", "Krakatoa", "Elastase", "Trypsin", "GluC", "Chymotrypsin", "Papain", "ProteinaseK", "Thermolysin"]
    #chain = ""

    #run = "NB8"
    #protein = 'QVQLQESGGGLVEPGGSLRLSCAVSGGSLNHYAMAWFRQAPGQEREGVACINRSGISTTYADSVKGRFTISRDNTKNTVWLQMNSLKPEDTGVYYCSAGKYYTVGHCDQDDYRGQGTQVTVSSAAADYKDHDGDYKDHDIDYKDDDDKGAAHHHHHH'
    #proteases = ["Vesuvius", "Krakatoa", "Elastase", "Trypsin", "GluC", "Chymotrypsin", "Papain", "ProteinaseK", "Thermolysin"]
    #chain = ""

    #run = "NB10"
    #protein = 'QVQLQESGGGLVQAGGSLRLSCAASGRTFDDYSMGWFRQAPGKEYEFVASINWSGSYTYYTDSVKGRFTISRDNAKNTVYLQMNSLKPDDTAVYYCAARDSIGVAVRRIDYDYWGQGTQVTVSSAAADYKDHDGDYKDHDIDYKDDDDKGAAHHHHHH'
    #proteases = ["Vesuvius", "Krakatoa", "Elastase", "Trypsin", "GluC", "Chymotrypsin", "Papain", "ProteinaseK", "Thermolysin"]
    #chain = ""

    #run = "NB12"
    #protein = 'QVQLQESGGGLVQAGGSLRLSCAASGRTFSSYAMAWFRQAPGKEREFVASISWSGDSTYYADSVKGRFTISRDNAKNTWYLQMNSLKPEDTAVYYCNTEEESTGTYYEWGQGTQVTVSSAAADYKDHDGDYKDHDIDYKDDDDKGAAHHHHHH'
    #proteases = ["Vesuvius", "Krakatoa", "Elastase", "Trypsin", "GluC", "Chymotrypsin", "Papain", "ProteinaseK", "Thermolysin"]
    #chain = ""

    #run = "NB13"
    #protein = 'QVQLQESGGGLVQPGGSLRLSCAASGSASSMYTLAWYRQAPGKQRELVALITSGHMTHYEDSVKGRFTISRDNAKEVLYLQMNSLKPEDTAVYFCNLHRLTSSDDDGRTWGQGTQVTVSSAAADYKDHDGDYKDHDIDYKDDDDKGAAHHHHHH'
    #proteases = ["Vesuvius", "Krakatoa", "Elastase", "Trypsin", "GluC", "Chymotrypsin", "Papain", "ProteinaseK", "Thermolysin"]
    #chain = ""

    #run = "BIND15"
    #protein = 'TEEELQKEIDEFKEKLLEKTRELVARAIEAQRAGEELWEQALNLAAGIYMQSISLLDMLELGKLAEAEEILKSIEELRKLLKEQAEKLKAKDPSKAPEMEELLKEVEETLKELEELVEKVKKFGGGGSEAAAKGGGGSGLNDIFEAQKIEWHEHHHHHH'
    #proteases = ["Vesuvius","Chymo","GluC","Papain","Krakatoa","ProtK","Thermo","Trypsin","Elastase"]
    #chain = ""

    #run = "BIND16"
    #protein = "SLSEEEKAELEELEKKALEKAAELVRRAIEAERAGRDLLAEALNLAAGILMLLVSAAEQAKRGDLAALKEQLEAAEKLLEILKEVAEQIKAGTPEERKAAEEALKLAEEAVKEIKRVLKLAEKAGGGGSEAAAKGGGGSGLNDIFEAQKIEWHEHHHHHHHH"
    #proteases = ["Vesuvius","Chymo","GluC","Papain","Krakatoa","ProtK","Thermo","Trypsin","Elastase"]
    #chain = ""

    #run = "BIND17"
    #protein = "KEELRAAAAELLAAAEALAEELRRLGLEEAAAHVLAAARHVAAALELIAATPASELNPELKREVAAHLREAAAHFEAAAEIVAAEDPLAGAMLREAALAARSMAAYVLHSSPEEALQQAAVFATGLAGAMLTMTGLVRERLAARGLNDIFEAQKIEWHEHHHHHH"
    #proteases = ["Vesuvius","Chymo","GluC","Papain","Krakatoa","ProtK","Thermo","Trypsin","Elastase"]
    #chain = ""

    run = "pa"
    chain = "heavy"

    params = {"ass_method": 'dbg',
              "conf": conf,
              "kmer_size": kmer_size,
              "min_overlap": min_overlap,
              "min_identity": min_identity,
              "max_mismatches": max_mismatches,
              "size_threshold": size_threshold
              }

    # Directories
    folder_outputs = f"../outputs/{run}{chain}"

    prep.create_directory(folder_outputs)
    
    combination_folder_out = os.path.join(folder_outputs, f"comb_{ass_method}_c{conf}_ks{kmer_size}_ts{size_threshold}_mo{min_overlap}_mi{min_identity}_mm{max_mismatches}")
    
    prep.create_subdirectories_outputs(combination_folder_out)

    # Data cleaning
    protein_norm = prep.normalize_sequence(protein)
    
    df = pd.read_csv(f"../input/{run}.csv")
    
    df['protease'] = df['experiment_name'].apply(lambda name: prep.extract_protease(name, proteases))
    
    df = prep.clean_dataframe(df)
    
    df['cleaned_preds'] = df['preds'].apply(prep.remove_modifications)
    
    cleaned_psms = df['cleaned_preds'].tolist()
    
    filtered_psms = prep.filter_contaminants(cleaned_psms, run, "../fasta/contaminants.fasta")
    
    df = df[df['cleaned_preds'].isin(filtered_psms)]
    
    df["mapped"] = df["cleaned_preds"].apply(lambda x: "True" if x in protein_norm else "False")
    
    df = df[df['conf'] > conf]
    
    df.reset_index(drop=True, inplace=True)
    
    final_psms = df['cleaned_preds'].tolist()

    # Assembly
    kmers = dbg.get_kmers(final_psms, kmer_size=kmer_size)
    
    edges = dbg.get_debruijn_edges_from_kmers(kmers)
    
    assembled_contigs = dbg.assemble_contigs(edges)
    
    assembled_contigs = sorted(assembled_contigs, key=len, reverse=True)
    
    assembled_contigs = list(set(assembled_contigs))
    
    assembled_contigs = [seq for seq in assembled_contigs if len(seq) > size_threshold]
    
    assembled_contigs = sorted(assembled_contigs, key=len, reverse=True)

    records = [Bio.SeqRecord.SeqRecord(Bio.Seq.Seq(contig), id=f"contig_{idx+1}",
                                    description=f"length: {len(contig)}") for idx,
                                    contig in enumerate(assembled_contigs)]
    
    Bio.SeqIO.write(records, f"{combination_folder_out}/contigs/{ass_method}_contig_{conf}_{run}.fasta", "fasta")

    mapped_contigs = map.process_protein_contigs_scaffold(assembled_contigs, protein_norm, max_mismatches, min_identity)
    
    df_contigs = map.create_dataframe_from_mapped_sequences(data = mapped_contigs)
    
    comp_stat.compute_assembly_statistics(df = df_contigs,
                                          sequence_type='contigs',
                                          output_folder = f'{combination_folder_out}/statistics',
                                          reference = protein_norm,
                                          **params
                                          )

    assembled_scaffolds = dbg.create_scaffolds(assembled_contigs, min_overlap)

    assembled_scaffolds = list(set(assembled_scaffolds))
    
    assembled_scaffolds = sorted(assembled_scaffolds, key=len, reverse=True)
    
    assembled_scaffolds = [scaffold for scaffold in assembled_scaffolds if len(scaffold) > size_threshold]

    assembled_scaffolds = dbg.merge_sequences(assembled_scaffolds)
    
    assembled_scaffolds = list(set(assembled_scaffolds))
    
    assembled_scaffolds = sorted(assembled_scaffolds, key=len, reverse=True)
    
    assembled_scaffolds = [scaffold for scaffold in assembled_scaffolds if len(scaffold) > size_threshold]

    records = []
    for i, seq in enumerate(assembled_scaffolds):
        record = Bio.SeqRecord.SeqRecord(Bio.Seq.Seq(seq), id=f"scaffold_{i+1}", description=f"length: {len(seq)}")
        records.append(record)

    Bio.SeqIO.write(records, f"{combination_folder_out}/scaffolds/{ass_method}_scaffold_{conf}_{run}.fasta", "fasta")

    mapped_scaffolds = map.process_protein_contigs_scaffold(assembled_contigs = assembled_scaffolds,
                                                            target_protein = protein_norm,
                                                            max_mismatches = max_mismatches,
                                                            min_identity = min_identity
                                                            )
    
    df_scaffolds_mapped = map.create_dataframe_from_mapped_sequences(data = mapped_scaffolds)
    
    comp_stat.compute_assembly_statistics(df = df_scaffolds_mapped, sequence_type='scaffolds',
                                          output_folder = f"{combination_folder_out}/statistics",
                                          reference = protein_norm, **params)    