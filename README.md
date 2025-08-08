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
- [Prerequisites and Installation](#prerequisites-and-installation)
- [Getting Started](#getting-started)
- [License](#license)
- [Acknowledgments](#acknowledgments)
- [References](#references)

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

| File / Folder       | Description                                                                  |
|---------------------|------------------------------------------------------------------------------|
| `environment.linux.yml`        | Conda environment definition with required dependencies for linux |
| `environment.osx-arm64.yaml`   | Conda environment definition with required dependencies for OS    |
| `README.md`         | Project documentation                                                        |
| `examples/`         |                                                                              |
| `fasta/`            | Known contaminants and example FASTA sequences                               |
| `images/`           | Logos and workflow diagrams (PNG, SVG, PDF)                                  |
| `inputs/`           | Example datasets (e.g., BSA, antibody, nanobody)                             |
| `json/`             | JSON metadata for peptide color coding and analysis                          |
| `notebooks/`        | Jupyter notebooks for visualization and exploration                          |
| `src/`              | Core scripts to run the InstaNexus pipeline                                  |

---

## Prerequisites and Installation

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
git clone https://github.com/your-username/instanexus.git
cd instanexus
```

### 2. Create the conda environment

Create instanexus conda environment for linux

```bash
conda env create -f environment.linux.yml
```

Create instanexus conda environment for OS

```bash
conda env create -f environment.osx-arm64.yaml
```

### 3. Activate the environment

```bash
conda activate instanexus
```

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
