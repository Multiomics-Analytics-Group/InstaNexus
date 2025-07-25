#!/usr/bin/env python

r"""
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

import os
import pandas as pd


base_directory = "/home/marcor/works/assembly/outputs/ma1/light"  # Change this to your actual path

merged_data = []

log_file_path = os.path.join(base_directory, "missing_statistics_log.txt")
missing_folders = []

for folder in os.listdir(base_directory):
    folder_path = os.path.join(base_directory, folder)
    statistics_path = os.path.join(folder_path, "statistics")
    
    if os.path.isdir(statistics_path):
        contigs_file = os.path.join(statistics_path, "contigs_stats_default.tsv")
        scaffolds_file = os.path.join(statistics_path, "scaffolds_stats_default.tsv")
        
        if os.path.exists(contigs_file) and os.path.exists(scaffolds_file):
            # Read both files
            df_contigs = pd.read_csv(contigs_file, sep="\t")
            df_scaffolds = pd.read_csv(scaffolds_file, sep="\t")
            
            # Add a column indicating the source folder
            df_contigs["Source_Folder"] = folder
            df_scaffolds["Source_Folder"] = folder
            
            # Append to the merged list
            merged_data.append(df_contigs)
            merged_data.append(df_scaffolds)
        else:
            # Log missing files
            missing_folders.append(folder)

if merged_data:
    final_df = pd.concat(merged_data, ignore_index=True)
    final_df.to_csv(os.path.join(base_directory, "merged_statistics.tsv"), sep="\t", index=False)
    print("Merged DataFrame saved as 'merged_statistics.tsv'.")
else:
    print("No valid statistics files found.")

# Write the log file
if missing_folders:
    with open(log_file_path, "w") as log_file:
        log_file.write("Missing statistics files in the following folders:\n")
        log_file.write("\n".join(missing_folders))
    print(f"Log file saved as '{log_file_path}'.")
else:
    print("All statistics files were found.")
