# LD-scaling genome scans for selection

This repository contains code to reproduce the analyses presented in:

> **Linkage disequilibrium scaling improves robustness and power to detect genomic regions under selection** (<https://www.biorxiv.org/cgi/content/short/2026.01.19.700334v1>)

The study introduces an LD-scaled association statistic (**F′**), a region-based outlier framework, and a consistency-based approach for integrating uncertainty across parameter choices and inference methods.

------------------------------------------------------------------------

## Overview

Genome-scan and genotype–environment association methods often lose power or generate false positives when linkage disequilibrium (LD) and population structure are strong. This project provides:

-   **LD-scaled association statistics (F′)** for LFMM and EMMAX
-   **Outlier-region (OR)–based inference**, rather than SNP-level tests
-   **Consistency scores (C)** integrating uncertainty across parameter settings
-   **Simulation-based benchmarking** using forward-in-time Nemo simulations
-   **Empirical re-analyses** of three- and nine-spined stickleback datasets

All analyses are implemented in **R**. This repository contains to code to analyse data, parsed simulated data (500 individuals with MAF\>0.05 and 125 individuals with MAF\>0.1) and intermediate data files are available from Zenondo (10.5281/zenodo.18466057). Raw data files from simulations are \>30gb so are only available upon request.

------------------------------------------------------------------------

## Repository structure

```         
├── README.md
├── Nemo/                               # Input files related to simulations (exl. raw data)
│   ├── generate_cmd_for_simulations.R/ # Main script to generate files for simulations       
├── R_sim/                              # Scripts to parse and analyse simulated data
│   ├── parse_sim_data.R/               # Parsing of simulate data  
│   ├── outlier_analyses_sim.R/         # Main analyses of simulated data
├── empirical_data/                     # 
│   ├── 3sp/                            # raw data - dowload from Zenondo
│   ├── 9sp/                            # raw data - dowload from Zenondo
│   ├── R/                              
│     ├── 3sp_sticklebacks.R/           # R-code to analyses three-spined sticklebacks
│     ├── 9sp_sticklebacks.R/           # R-code to analyses nine-spined sticklebacks
├── R/                                  # R-functions that are common for many analyses
├── figures/                            # Figures for manuscript
│   └── R/                              # Scripts to reproduce manuscript figures
└── data/                               # Intermediate data files (Dowload from Zenondo)
└── parsed_data/                        # Parsed data form simulations (Dowload from Zenondo)
└── sim_results/                        # Outlier results from simulations (Zenondo)
```

------------------------------------------------------------------------

## Requirements

Analyses were run in **R ≥ 4.2**.

Main R package dependencies include:

-   **LEA** (LFMM)
-   **SNPRelate**
-   **igraph**
-   **gstat**
-   **effectsize**
-   **tidyverse**

Clone the repository:

``` bash
git clone https://github.com/<username>/LD-scaling-genome-scans.git
```

## Reproducing the analyses

### 1. Simulations

Forward-in-time simulations were generated using **Nemo** and post-processed in R.

-   Simulation design: continuous landscape, spatially autocorrelated optima, variable gene flow and selection strength
-   Recombination maps derived from stickleback Marey maps

Script:

```         
Nemo/generate_cmd_for_simulations.R
```

Produces all necessary input files for simulations. Creates file containing all command-line arguments:

```         
Nemo/cmds.sh
```

Each line starts a separate run on a cluster (sbatch).

> **Note:** Raw Nemo outputs are large and only one is included in the Nemo/folder. Some of the input files used for simulations are provided as well as processed summary files used in the manuscript.

------------------------------------------------------------------------

### 2. Parsing simulated data

-   Takes output from simulations (one example data set is provided: ./Nemo/chr1_V0.5_c1_rep4.tar.gz) extracts neutral and causal loci for 500 fixed individuals (./data/keep_500.rds) with maf\>0.05 and produces files genotype files (*GTs*) and map files (*map*) and the environment (*env*) as well as the Nemo output containing population genetic info (*sim_data*).

### 3. Simulated data: LD-scaled outlier detection and benchmarking

```         
/R_sim/outlier_analyses_sim.R
```

This script performs the full downstream analysis of simulated genomic data used to benchmark LD-scaled genome-scan methods. Parsed simulation outputs containing genotypes for 500 individuals are subsampled to a fixed set of 125 individuals to match the empirical sampling design and to reduce computational cost while preserving LD structure.

For each simulation replicate, the pipeline estimates chromosome-specific LD-decay, performs genotype–environment association analyses using **EMMAX** and **LFMM**, applies LD-scaling to association statistics, and defines **outlier regions (ORs)** across thousands of random parameter combinations. Evidence is integrated across parameter draws using **consistency scores (C)**, and performance is evaluated at the level of ORs rather than individual SNPs.

True and false positives are defined using known causal QTNs, allowing calculation of OR-level precision–recall metrics. These are aggregated into **AUC-PR\*** and **AUC-PRC** scores to quantify robustness and power across demographic and selective scenarios. All intermediate and final results are saved to (./data/) and reused by downstream plotting and summary scripts.

### 4. Core R-functions

Core R-functions are found in:

```         
/R/
```

### 4. Analyses of empirical stickleback data

Empirical stickleback datasets originate from Fang et al. (2021):

-   Zenodo accession: **4722879**

The analyses of three- and nine-spined sticklebacks follow the same pipeline as the simulated data, except the AUC-analyses that are only possible for data where the ground truth is known.

```         
/empirical_data/R/3sp_sticklebacks.R
/empirical_data/R/9sp_sticklebacks.R
```

Code for Manhattan plots (for Fig. 5) are given in these files as well. Intermediate data files (LD-decay stats, and draws for different parameter combinations etc) are saved to:

```         
/empirical_data/3sp/
/empirical_data/3sp/
```

Manhattanplots for the empirical data are joined with simulated data in the code for Fig5.

```         
/figures/R/Fig5.R
```

### 5. Figures

## Relation to LDscanR

Core methods introduced here will be implemented in the R package **LDscanR**:

> <https://github.com/genomeinsights/LDscanR>

This repository provides the **full analysis and benchmarking pipeline** used in the manuscript, rather than a minimal package interface.

------------------------------------------------------------------------

## Citation

If you use this code, please cite:

> Kemppainen, P. & Guillaume, F. *Linkage disequilibrium scaling improves robustness and power to detect genomic regions under selection*. <https://www.biorxiv.org/cgi/content/short/2026.01.19.700334v1>

------------------------------------------------------------------------

## Contact

**Petri Kemppainen** 📧 [petri\@genomeinsights.fi](mailto:petri.kemppainen@helsinki.fi)
