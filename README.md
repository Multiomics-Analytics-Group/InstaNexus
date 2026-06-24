<p align="center">
  <img src="https://raw.githubusercontent.com/Multiomics-Analytics-Group/InstaNexus/main/docs/source/assets/instanexus_logo%202.svg" width="600" alt="InstaNexus logo">
</p>

<p align="center"><em>A de novo protein sequencing workflow</em></p>

<p align="center">
  <a href="https://github.com/pre-commit/pre-commit"><img src="https://img.shields.io/badge/pre--commit-enabled-brightgreen?logo=pre-commit" alt="pre-commit"></a>
  <a href="https://github.com/astral-sh/ruff"><img src="https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/astral-sh/ruff/main/assets/badge/v2.json" alt="Ruff"></a>
  <img src="https://img.shields.io/badge/license-MIT-green" alt="License">
  <img src="https://img.shields.io/badge/python-3.10+-blue" alt="Python">
</p>

---

## Table of Contents
- [Introduction](#introduction)
- [Features](#features)
- [Workflow Diagram](#workflow-diagram)
- [Repository Structure](#repository-structure)
- [Installation](#installation)
- [Command-Line Usage](#command-line-usage)
- [Hyperparameter Optimization](#hyperparameter-optimization)
- [License](#license)
- [Acknowledgments](#acknowledgments)
- [References](#references)
- [Citation](#citation)

---

## Introduction

InstaNexus is a generalizable, end-to-end workflow for direct protein sequencing, tailored to reconstruct full-length protein therapeutics such as antibodies and nanobodies. It integrates AI-driven de novo peptide sequencing with optimized assembly and scoring strategies to maximize accuracy, coverage, and functional relevance.

This pipeline enables robust reconstruction of critical protein regions, advancing applications in therapeutic discovery, immune profiling, and protein engineering.

---

## Features

- 🧬 Supports De Bruijn Graph and Greedy-based assembly
- ⚗️ Handles multiple protease digestions (Trypsin, LysC, GluC, etc.)
- 🧹 Integrated contaminant removal and confidence filtering
- 🧩 Clustering, alignment, and consensus sequence reconstruction
- 🔗 Integrates with external tools:
  - [MMseqs2](https://github.com/soedinglab/MMseqs2) for fast clustering
  - [Clustal Omega](https://www.ebi.ac.uk/Tools/msa/clustalo/) for high-quality alignment
- 📊 Output-ready for downstream analysis and visualization

---

## Workflow Diagram

<p align="center">
  <img src="https://raw.githubusercontent.com/Multiomics-Analytics-Group/InstaNexus/main/docs/source/assets/instanexus_panel.png" width="900" alt="InstaNexus Workflow">
</p>

---

## Repository Structure


| Folder / File | Description |
|----------------|-------------|
| `docs/` | Sphinx documentation, tutorials, and images |
| `fasta/` | FASTA reference and contaminant sequences |
| `inputs/` | Example input CSV files |
| `json/` | Metadata and parameter configuration files |
| `outputs/` | Generated results (created during execution) |
| `src/instanexus/` | Core InstaNexus package |
| `src/instanexus/main.py` | Runs the full pipeline |
| `src/instanexus/preprocessing.py` | Module for data cleaning |
| `src/instanexus/assembly.py` | Module for sequence assembly |
| `src/instanexus/clustering.py` | Module for clustering (mmseqs2) |
| `src/instanexus/alignment.py` | Module for alignment (clustalo) |
| `src/instanexus/consensus.py` | Module for consensus generation |
| `src/instanexus/opt/` | Grid search and optimization workflows |
| `tests/` | Pytest unit and integration tests |
| `pyproject.toml` | Package metadata, dependencies, and entry point |
| `.pre-commit-config.yaml` | Pre-commit hook configuration |

---

## Installation

InstaNexus requires Python 3.10+, [uv](https://docs.astral.sh/uv/), **MMseqs2**, and **Clustal Omega**.

- [uv](https://docs.astral.sh/uv/) — fast Python package manager
- [MMseqs2](https://github.com/soedinglab/MMseqs2)
- [Clustal Omega](https://www.ebi.ac.uk/Tools/msa/clustalo/)

---

## Getting Started

### Option 1: Install from PyPI

```bash
pip install instanexus
```

### Option 2: Install from Source (for Developers)

#### Clone the repository:

```bash
git clone git@github.com:Multiomics-Analytics-Group/InstaNexus.git
cd InstaNexus
```

#### Install uv (if not already installed):
```bash
curl -LsSf https://astral.sh/uv/install.sh | sh
```

#### Sync the environment:
```bash
uv sync --all-extras
```

#### Set up pre-commit hooks:
```bash
uv run pre-commit install --hook-type pre-commit --hook-type commit-msg
```

#### Verify the installation:
```bash
uv run instanexus --help
```

---

## Command-line usage

After installation (and adding the `[project.scripts]` entry point), you can run the entire InstaNexus pipeline using the `instanexus` command.

All parameters for preprocessing, assembly, clustering, and consensus are provided in a single call. The pipeline will automatically create a unique, timestamped output folder for that specific combination of parameters.

```bash
instanexus --help
```

Example: Run the full pipeline
This command runs the complete workflow:

Preprocesses the input CSV.

Assembles using dbg (De Bruijn graph).

Clusters the resulting scaffolds.

Aligns the clusters.

Generates consensus sequences.

```bash
instanexus \
    --input-csv inputs/bsa.csv \
    --folder-outputs outputs \
    --metadata-json-path json/sample_metadata.json \
    --contaminants-fasta-path fasta/contaminants.fasta \
    --assembly-mode dbg \
    --conf 0.9 \
    --kmer-size 7 \
    --size-threshold 12 \
    --min-overlap 3 \
    --min-seq-id 0.85 \
    --coverage 0.8
```

The results for this specific run will be saved in a unique directory, such as:```outputs/bsa/dbg_c0.9_ks7_mo3_ts12/```



---

## License

This project is licensed under the [MIT License](LICENSE).

---

## Acknowledgments

InstaNexus was developed at **DTU Biosustain** and **DTU Bioengineering**.

We are grateful to the **DTU Bioengineering Proteomics Core Facility** for maintenance and operation of mass spectrometry instrumentation.

We also thank the **Informatics Platform at DTU Biosustain** for their support during the development and optimization of InstaNexus.

Special thanks to the users and developers of:
- [MMseqs2](https://github.com/soedinglab/MMseqs2)
- [Clustal Omega](https://www.ebi.ac.uk/Tools/msa/clustalo/)

---

## References

1. Hauser, M., et al. **MMseqs2: ultra fast and sensitive sequence searching**. *Nature Biotechnology* 35, 1026–1028 (2016). https://doi.org/10.1038/nbt.3988  
2. Sievers, F., et al. **Fast, scalable generation of high-quality protein multiple sequence alignments using Clustal Omega**. *Molecular Systems Biology* 7, 539 (2011). https://doi.org/10.1038/msb.2011.75
3. Eloff, K., Kalogeropoulos, K., Mabona, A., Morell, O., Catzel, R., Rivera-de-Torre, E., ... & Jenkins, T. P. (2025). **InstaNovo enables diffusion-powered de novo peptide sequencing in large-scale proteomics experiments.** Nature Machine Intelligence, 1-15.

---

## Citation

If you find this project useful in your research or work, please cite our publication in Molecular & Cellular Proteomics: [Generalizable direct protein sequencing with InstaNexus](https://doi.org/10.1016/j.mcpro.2026.101547).

```bibtex
@article{reverenna2026generalizable,
  title={Generalizable direct protein sequencing with InstaNexus},
  author={Reverenna, Marco and Nielsen, Maike Wennekers and Wolff, Darian Stephan and Daniel, Jemma and Lytra, Elpida and 
          Thumtecho, Suthimon and Colaianni, Pasquale D and Ljungars, Anne and Laustsen, Andreas H and Schoof, Erwin M and
          Van Goey, Jeroen and Jenkins, Timothy P and Lukassen, Marie V and Santos, Alberto and Kalogeropoulos, Konstantinos},
  journal={Molecular \& Cellular Proteomics},
  volume={25},
  number={4},
  pages={101547},
  year={2026},
  doi={10.1016/j.mcpro.2026.101547},
  publisher={Elsevier}
}
```
