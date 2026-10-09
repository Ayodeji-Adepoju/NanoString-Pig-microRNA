# ============================================================================
# NanoString renal miRNA analysis pipeline: human NanoString panel -> pig miRNAs
# ============================================================================
#
# This single script combines:
#   1. miRBase human/pig sequence retrieval and exact mature-sequence mapping
#   2. NanoString RCC import and normalization
#   3. human-to-Sus scrofa miRNA annotation
#   4. differential-expression analysis and visualization
#   5. target retrieval, enrichment, ranking, and Cytoscape-ready exports
#
# GitHub / reproducibility notes
# ------------------------------
# No user-specific absolute file paths are stored in this script.
#
# Recommended repository structure:
#
# project/
# |-- NanoString_miRNA_pipeline_GitHub.R
# |-- data/
# |   |-- rcc/                 # NanoString .RCC files
# |   `-- metadata.csv         # optional but recommended
# |-- reference/               # downloaded miRBase files and mapping tables
# `-- results/                 # all analysis outputs
#
# Run from the repository root:
#   Rscript NanoString_miRNA_pipeline_GitHub.R
#
# Alternatively set a project location without editing this script:
#   NANOSTRING_PROJECT_DIR=/path/to/project Rscript NanoString_miRNA_pipeline_GitHub.R
#
# Optional sample prefix:
#   NANOSTRING_SAMPLE_PREFIX=<prefix> Rscript NanoString_miRNA_pipeline_GitHub.R
#
# ============================================================================


# ============================================================================
# PART 0 - Packages
# ============================================================================

cran_pkgs <- c(
  "remotes", "nanostringr", "dplyr", "tidyr", "tibble", "readr",
  "ggplot2", "ggrepel", "pheatmap", "matrixStats", "httr2", "stringr"
)

to_install_cran <- cran_pkgs[
  !cran_pkgs %in% rownames(installed.packages())
]

if (length(to_install_cran) > 0) {
  install.packages(to_install_cran)
}

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

bioc_pkgs <- c(
  "Biostrings",
  "limma",
  "ComplexHeatmap",
  "circlize",
  "clusterProfiler",
  "org.Hs.eg.db",
  "ReactomePA",
  "multiMiR"
)

to_install_bioc <- bioc_pkgs[
  !bioc_pkgs %in% rownames(installed.packages())
]

if (length(to_install_bioc) > 0) {
  BiocManager::install(
    to_install_bioc,
    ask = FALSE,
    update = FALSE
  )
}

library(nanostringr)
library(dplyr)
library(tidyr)
library(tibble)
library(readr)
library(ggplot2)
library(ggrepel)
library(pheatmap)
library(matrixStats)
library(httr2)
library(stringr)
library(Biostrings)
library(limma)
library(ComplexHeatmap)
library(circlize)
