---
editor_options: 
  markdown: 
    wrap: 72
---

# LD-scaling genome scans for selection

This repository contains code to reproduce the analyses presented in:

> **Linkage disequilibrium scaling improves robustness and power to
> detect genomic regions under selection** Petri Kemppainen & Frédéric
> Guillaume *Proceedings of the National Academy of Sciences (PNAS)*, in
> review

The study introduces an LD-scaled association statistic (**F′**), a
region-based outlier framework, and a consistency-based approach for
integrating uncertainty across parameter choices and inference methods.

------------------------------------------------------------------------

## Overview

Genome-scan and genotype–environment association methods often lose
power or generate false positives when linkage disequilibrium (LD) and
population structure are strong. This project provides:

-   **LD-scaled association statistics (F′)** for LFMM and EMMAX
-   **Outlier-region (OR)–based inference**, rather than SNP-level tests
-   **Consistency scores (C)** integrating uncertainty across parameter
    settings
-   **Simulation-based benchmarking** using forward-in-time Nemo
    simulations
-   **Empirical re-analyses** of three- and nine-spined stickleback
    datasets

All analyses are implemented in **R**. This repository contains to code
to analyse data, all raw data and intermediate data files are available
from Zenondo (10.5281/zenodo.18466057). Place holders for folders
containing this data exist but are empty.

------------------------------------------------------------------------

## Repository structure

```         
├── README.md
├── Nemo/                               # Input files related to simulations (partially empty)
│   ├── generate_cmd_for_simulations.R/ # Main script to generate files for simulations       
├── R_sim/                              # Scripts to parse and analyse simulated data
│   ├── parse_sim_data.R/                 # Parsing of simulate data  
├── R_emp/                              # Scripts analyse empirical data
├── empirical_data/                     # Empirical stickleback data (empty)
├── R/                                  # R-functions that are common for all analyses
├── figures/                            # Figures for manuscript
│   └── R/                              # Scripts to reproduce manuscript figures
└── data/                               # Intermediate data produced by analyses (empty)
│   ├── parsed_data/                    # Parsed data from simulations (empty)
│   ├── sim_results/                    # Results for analyses of simulated data (empty)
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

Exact package versions are recorded in `renv.lock` (recommended).

------------------------------------------------------------------------

## Installation

Clone the repository:

``` bash
git clone https://github.com/<username>/LD-scaling-genome-scans.git
cd LD-scaling-genome-scans
```

Restore the R environment:

``` r
renv::restore()
```

------------------------------------------------------------------------

## Reproducing the analyses

### 1. Simulations

Forward-in-time simulations were generated using **Nemo** and
post-processed in R.

-   Simulation design: continuous landscape, spatially autocorrelated
    optima, variable gene flow and selection strength
-   Recombination maps derived from stickleback Marey maps

Script:

``` r
Nemo/generate_cmd_for_simulations.R
```

Produces all necessary input files for simulations. Creates file
containing all command-line arguments:

``` r
Nemo/cmds.sh
```

Each line starts a separate run on a cluster (sbatch).

> **Note:** Raw Nemo outputs are large and not included here. Some of
> the input files used for simulations are provided as well as processed
> summary files used in the manuscript.

------------------------------------------------------------------------

### 2. Parsing and analysing simulated data

### 2. Association analyses

Genotype–environment associations were computed using:

-   **LFMM**
-   **EMMAX**

Scripts:

``` r
scripts/association/run_lfmm.R
scripts/association/run_emmax.R
```

------------------------------------------------------------------------

### 3. LD scaling and F′ statistic

Local LD was summarized in sliding windows and used to compute the
LD-scaled statistic **F′** via a permutation-based quantile
transformation.

Scripts:

``` r
scripts/ld_scaling/estimate_ld_decay.R
scripts/ld_scaling/calc_Fprime.R
```

------------------------------------------------------------------------

### 4. Outlier regions and consistency scores

Outlier SNPs were clustered into **outlier regions (ORs)** based on LD
and physical distance thresholds. Consistency scores (**C**) integrate
results across parameter draws.

Scripts:

``` r
scripts/outlier_regions/define_ORs.R
scripts/benchmarking/consistency_scores.R
```

------------------------------------------------------------------------

### 5. Benchmarking and figures

Performance was evaluated using OR-level precision–recall metrics
(**AUC-PR**\* and **AUC-PRC**).

Scripts:

``` r
scripts/benchmarking/auc_pr.R
scripts/figures/make_figures.R
```

------------------------------------------------------------------------

## Empirical data availability

Empirical stickleback datasets originate from Fang et al. (2021):

-   Zenodo accession: **4722879**

Due to size, genotype data are **not redistributed** here. The
repository assumes the user downloads these datasets separately and
updates the paths in `scripts/config.R`.

------------------------------------------------------------------------

## Relation to LDscanR

Core methods introduced here are implemented in the R package
**LDscanR**:

> <https://github.com/genomeinsights/LDscanR>

This repository provides the **full analysis and benchmarking pipeline**
used in the manuscript, rather than a minimal package interface.

------------------------------------------------------------------------

## Citation

If you use this code, please cite:

> Kemppainen, P. & Guillaume, F. *Linkage disequilibrium scaling
> improves robustness and power to detect genomic regions under
> selection*. PNAS.

(Updated with DOI upon acceptance.)

------------------------------------------------------------------------

## Contact

**Petri Kemppainen** University of Helsinki 📧
[petri.kemppainen\@helsinki.fi](mailto:petri.kemppainen@helsinki.fi)
