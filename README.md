<p align="center">
  <img src="images/instanexus_logo 2.svg" width="600" alt="InstaNexus logo">
</p>

<p align="center"><em>A de novo protein sequencing workflow</em></p>

<p align="center">
  <img src="https://img.shields.io/badge/environment-conda-blue" alt="Conda">
  <img src="https://img.shields.io/badge/license-MIT-green" alt="License">
  <img src="https://img.shields.io/badge/python-3.9+-blue" alt="Python">
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
  <img src="images/instanexus_panel.png" width="900" alt="InstaNexus Workflow">
</p>

---

## Repository Structure


| Folder / File | Description |
|----------------|-------------|
| `environment.linux.yml` | Conda environment for Linux systems |
| `environment.osx-arm64.yaml` | Conda environment for macOS (Apple Silicon) |
| `src/instanexus/` | Core InstaNexus package (modules + CLI) |
| `src/instanexus/__main__.py` | Entry point for CLI (`instanexus` command) |
| `src/instanexus/script_dbg.py` | De Bruijn Graph-based assembly |
| `src/instanexus/script_greedy.py` | Greedy-based peptide assembly |
| `src/opt/` | Grid search and optimization workflows |
| `fasta/` | FASTA reference and contaminant sequences |
| `inputs/` | Example input CSV files |
| `json/` | Metadata and parameter configuration files |
| `notebooks/` | Jupyter notebooks for analysis and visualization |
| `images/` | Logos and workflow figures |
| `outputs/` | Generated results (created during execution) |

---

## Installation

- [Conda](https://docs.conda.io/en/latest/)
- [MMseqs2](https://github.com/soedinglab/MMseqs2)
- [Clustal Omega](https://www.ebi.ac.uk/Tools/msa/clustalo/)

> [!IMPORTANT]
> MMseqs2 and Clustal Omega are available through Conda, but compatibility depends on your system architecture.
> - 🔍 [Clustal Omega on Anaconda.org](https://anaconda.org/search?q=clustalo)   

---

## Getting Started

Follow these steps to clone the repository and set up the environment using Conda:

### 1. Clone the repository

To clone and set up the environment:

```bash
git clone git@github.com:Multiomics-Analytics-Group/InstaNexus.git
cd instanexus
```

### 2. Create and activate the Conda environment

Create instanexus conda environment for linux.

```bash
conda env create -f environment.linux.yml
```

Create instanexus conda environment for OS.

```bash
conda env create -f environment.osx-arm64.yaml
```

Activate:

```bash
conda activate instanexus
```

---

### 3. Install InstaNexus as a local package

```
pip install -e .
```

Then verify the CLI installation:

```
instanexus --version
```

---

## Command-line usage

After activating the environment, you can run InstaNexus directly from the terminal:
```bash
instanexus --help
```

### Run De Bruijn graph assembly

```
instanexus dbg --input_csv inputs/sample.csv --chain light --folder_outputs outputs --reference
```

### Run greedy assembly

```
instanexus greedy --input_csv inputs/sample.csv --folder_outputs outputs
```




---

## Hyperparameter Optimization

To launch the hyperparameter grid search, run the following command from the project root (the folder containing ```src/``` and ```json/```):

```bash
python -m src.opt.gridsearch
```
**Adjusting Parameters**

Grid search parameters for both the De Bruijn graph (dbg) and Greedy (greedy) assembly methods are defined in:

```bash
json/gridsearch_params.json
```

To test more (or fewer) combinations, edit the arrays for each parameter in this file.

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

If you find this project useful in your research or work, please cite it as:

Reverenna M., Nielsen M. W., Wolff D. S., Lytra E., Colaianni P. D., Ljungars A., Laustsen A. H., Schoof E. M., Van Goey J., Jenkins T. P., Lukassen M. V., Santos A., Kalogeropoulos K. (2025). *Generalizable direct protein sequencing with InstaNexus* [Preprint]. bioRxiv. https://doi.org/10.1101/2025.07.25.666861
