# InstaNexus demo

This folder contains a self-contained, teaching-oriented walkthrough of the
InstaNexus assembly pipeline applied to a nanobody dataset.

## Contents

- [`nanobody_demo.ipynb`](nanobody_demo.ipynb) — end-to-end notebook that takes
  raw de novo PSM predictions (`inputs/nb1.csv`), runs the InstaNexus CLI, and
  inspects the resulting scaffolds, consensus sequences, coverage, and sequence
  logos.

## Running the notebook in a Codespace

1. On GitHub, go to **Code → Codespaces → New codespace** on the
   `Multiomics-Analytics-Group/InstaNexus` repository.
2. Wait for the container to build — the `.devcontainer/` setup installs `uv`,
   `MMseqs2`, `Clustal Omega`, and syncs the Python environment automatically.
3. Open `demo/nanobody_demo.ipynb` and select the `.venv` kernel
   (`${workspaceFolder}/.venv/bin/python`).
4. Run the cells from top to bottom.

## Running the notebook on Binder

Click the badge below (or open the link directly) to launch the notebook in a
temporary, browser-based environment — no installation required:

[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/Multiomics-Analytics-Group/InstaNexus/main?filepath=demo/nanobody_demo.ipynb)

<https://mybinder.org/v2/gh/Multiomics-Analytics-Group/InstaNexus/main?filepath=demo/nanobody_demo.ipynb>

## Parameters you can experiment with

The "Parameters" section near the top of the notebook collects everything a
student would typically want to change:

| Parameter | Meaning |
|---|---|
| `INPUT_CSV` | Raw PSM predictions to assemble (default: `inputs/nb1.csv`) |
| `METADATA_JSON` | Sample metadata (protein sequence, chain, proteases) |
| `CONTAM_FASTA` | Contaminant sequences to filter out before assembly |
| `OUTPUT_DIR` | Where the pipeline writes its results |
| `assembly_mode` | Assembly algorithm (`"dbg"`, `"greedy"`, `"dbg_weighted"`, ...) |
| `conf` | Minimum confidence score kept during preprocessing |
| `kmer_size` | k-mer length used by the De Bruijn graph assembler |
| `min_overlap` | Minimum overlap required to merge two reads/contigs |
| `size_threshold` | Minimum contig length kept after assembly |
| `min_seq_id` | Minimum sequence identity used by MMseqs2 clustering |
| `coverage` | Minimum alignment coverage used by MMseqs2 clustering |

Changing these values and re-running the pipeline cell is the fastest way to
see how each parameter affects the assembled scaffolds and consensus
sequences.
