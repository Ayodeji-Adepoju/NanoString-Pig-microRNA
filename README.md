# NanoString miRNA Analysis Pipeline for Neonatal Pig Kidney

## Overview

This repository contains an R-based analysis pipeline for profiling renal microRNA (miRNA) expression in neonatal pig kidney samples assayed with a **human NanoString nCounter miRNA panel**.

The workflow combines:

1. Retrieval of human and pig miRNA reference sequences from **miRBase**.
2. Exact mature-sequence matching between **human (hsa)** and **pig (ssc)** miRNAs.
3. Import and normalization of NanoString RCC files.
4. Exclusion of human miRNAs without a supported porcine counterpart.
5. Differential expression analysis of mapped porcine miRNAs.
6. Heatmap and volcano plot generation.
7. miRNA target-gene retrieval using **multiMiR**.
8. GO, KEGG, and Reactome enrichment analyses.
9. Generation of Cytoscape-compatible edge and node tables.
10. Export of intermediate and final analysis tables for reproducibility.

The pipeline was developed for a neonatal porcine model of sepsis-associated acute kidney injury (SA-AKI), but the structure can be adapted for other studies using human NanoString miRNA panels on porcine tissue.

---

## Main Script

The complete workflow is contained in:

```text
NanoString_miRNA_pipeline_GitHub.R
```

The script combines the original miRNA mapping workflow and the downstream NanoString analysis into a single reproducible pipeline.

---

## Recommended Repository Structure

```text
project/
├── NanoString_miRNA_pipeline_GitHub.R
├── README.md
│
├── data/
│   ├── rcc/
│   │   ├── sample_01.RCC
│   │   ├── sample_02.RCC
│   │   └── ...
│   │
│   └── metadata.csv
│
├── reference/
│
└── results/
```

### `data/rcc/`

Place all NanoString `.RCC` files in this directory.

### `data/metadata.csv`

A metadata file is recommended and should contain at least:

```text
SampleID,Group
sample_01,Control
sample_02,Control
sample_03,Control
sample_07,Sepsis
sample_08,Sepsis
sample_09,Sepsis
```

`SampleID` values must match the sample names detected from the NanoString RCC files.

If `metadata.csv` is not present, the pipeline attempts to infer group assignments from sample names containing:

- `Control` or `Sham`
- `Sepsis` or `CLP`

If groups cannot be inferred and exactly 12 samples are present, the script falls back to the original study layout of six Control and six Sepsis samples.

For public or independent reuse, providing `metadata.csv` is strongly recommended.

---

## Requirements

The pipeline requires R and installs missing packages automatically.

### CRAN packages

- `remotes`
- `nanostringr`
- `dplyr`
- `tidyr`
- `tibble`
- `readr`
- `ggplot2`
- `ggrepel`
- `pheatmap`
- `matrixStats`
- `httr2`
- `stringr`

### Bioconductor packages

- `Biostrings`
- `limma`
- `ComplexHeatmap`
- `circlize`
- `clusterProfiler`
- `org.Hs.eg.db`
- `ReactomePA`
- `multiMiR`

The script installs `BiocManager` if it is not already available.

---

## Running the Pipeline

From the repository root:

```bash
Rscript NanoString_miRNA_pipeline_GitHub.R
```

By default, the script assumes that the repository root is the working project directory.

A custom project directory can also be specified using an environment variable:

```bash
NANOSTRING_PROJECT_DIR=/path/to/project \
Rscript NanoString_miRNA_pipeline_GitHub.R
```

An optional NanoString sample-column prefix can also be provided:

```bash
NANOSTRING_SAMPLE_PREFIX=sample_prefix \
Rscript NanoString_miRNA_pipeline_GitHub.R
```

If no prefix is provided, all detected sample columns are used.

---

## Analysis Workflow

### 1. miRBase sequence retrieval

The pipeline downloads the current miRBase:

```text
mature.fa
hairpin.fa
```

files and extracts:

- human mature miRNAs (`hsa-`)
- pig mature miRNAs (`ssc-`)
- pig precursor miRNAs

Downloaded reference files are stored in:

```text
reference/
```

---

### 2. Human-to-pig miRNA mapping

Because the NanoString panel is human-based, the original assay features are initially annotated with human miRNA identifiers.

Human and pig mature miRNAs are matched using **exact mature-sequence identity**.

When multiple pig miRNAs share the same mature sequence, the pipeline resolves the preferred pig assignment using:

1. exact human-to-pig name correspondence where available;
2. mature-arm consistency (`-3p` or `-5p`);
3. the first supported sequence match if no more specific assignment is available.

The final mapping table is written to:

```text
reference/hsa_to_ssc_best_exact_sequence_match.csv
```

Only human miRNAs with a supported porcine sequence match are retained for downstream analysis.

---

### 3. NanoString RCC import

NanoString RCC files are imported using:

```r
nanostringr::read_rcc()
```

The pipeline identifies:

- endogenous miRNA features;
- negative controls;
- positive controls;
- housekeeping controls.

---

### 4. NanoString normalization

Normalization includes:

#### Background correction

Background is estimated for each sample from the negative controls as:

```text
mean negative-control count + 2 × standard deviation
```

Background-corrected counts are constrained to a minimum value of 1.

#### Positive-control normalization

Positive-control geometric means are calculated per sample and used to generate sample-specific scaling factors.

#### Housekeeping normalization

Where housekeeping probes are available, their geometric mean is used for an additional normalization step.

---

### 5. Porcine miRNA expression matrix

NanoString features are converted from human miRNA labels to mapped `ssc-miRNA` identifiers.

Features without a supported pig equivalent are excluded.

If multiple human features map to the same porcine miRNA, normalized expression values are averaged.

Outputs include:

```text
endogenous_normalized_counts_matrix_ssc_labels_mapped_only.csv
endogenous_log2_matrix_ssc_labels_mapped_only.csv
final_feature_id_mapping_mapped_only.csv
mapping_summary.csv
```

---

## Differential Expression Analysis

Differential expression between Sepsis and Control groups is performed using **limma**.

The model contrast is:

```text
Sepsis - Control
```

The resulting table includes:

- `logFC`
- `AveExpr`
- `t`
- `P.Value`
- `adj.P.Val`
- `B`
- mapped `ssc_miRNA`
- corresponding `hsa_query_miRNA`

The pipeline retains adjusted P values in the result table, but the study-specific candidate prioritization can use predefined raw P-value and fold-change thresholds as appropriate.

Primary result file:

```text
limma_results_ssc_labels_mapped_only.csv
```

---

## Heatmap

Significant miRNAs are visualized with `ComplexHeatmap`.

The heatmap:

- uses row-wise z-scored expression;
- groups samples by experimental condition;
- displays individual biological replicates;
- clusters miRNAs;
- preserves experimental group annotation.

Outputs include:

```text
ComplexHeatmap_grouped_ordered_samples_ssc_mapped_only.pdf
ComplexHeatmap_grouped_ordered_samples_ssc_mapped_only.png
```

---

## Volcano Plot

The volcano plot displays:

```text
x-axis: log2 fold change
y-axis: -log10(P value)
```

Default candidate thresholds in the script are:

```text
P < 0.05
|log2FC| >= 0.7
```

Significantly upregulated and downregulated miRNAs are highlighted separately, and selected miRNAs are labeled.

---

## miRNA Target-Gene Analysis

Significant porcine miRNAs are mapped back to their matched human miRNA identifiers for target retrieval because the target databases queried through `multiMiR` are primarily human annotated.

Validated interactions are retrieved using:

```r
multiMiR::get_multimir(
  table = "validated",
  org = "hsa"
)
```

Target information may include:

- miRNA identifier;
- target gene symbol;
- Entrez ID;
- source database;
- support type;
- experimental evidence.

The porcine miRNA identifier is retained alongside the human query miRNA.

---

## Functional Enrichment

Target genes are mapped to human Entrez identifiers for downstream functional annotation.

The pipeline performs:

### Gene Ontology Biological Process

```r
clusterProfiler::enrichGO()
```

### KEGG pathway enrichment

```r
clusterProfiler::enrichKEGG()
```

### Reactome pathway enrichment

```r
ReactomePA::enrichPathway()
```

These analyses provide pathway-level context for miRNAs altered during experimental sepsis.

---

## Cytoscape Export

The pipeline generates files suitable for import into Cytoscape.

### Edge table

Contains miRNA-to-target-gene relationships:

```text
miRNA_target_edges_for_Cytoscape_ssc_mapped_only.csv
```

### Node table

Contains node-level annotations including:

- miRNA/gene node type;
- miRNA log2 fold change;
- P value;
- adjusted P value;
- matched human miRNA;
- target counts;
- regulator counts.

Output:

```text
cytoscape_nodes_ssc_mapped_only.csv
```

---

## Output Directories

The pipeline automatically creates:

```text
reference/
results/
```

`reference/` contains downloaded miRBase files and human-to-pig annotation tables.

`results/` contains normalized matrices, differential-expression results, target-gene tables, pathway-enrichment outputs, figures, Cytoscape files, and session information.

---

## Reproducibility

At the end of the analysis, the pipeline saves:

```text
sessionInfo_mapped_only.txt
```

This records the R version and package versions used during the analysis.

For full reproducibility, users should also preserve:

- raw RCC files;
- `metadata.csv`;
- the generated miRBase mapping tables;
- the R session information;
- the exact version of this script used for analysis.

---

## Important Methodological Note

The NanoString assay used in this workflow is a **human miRNA panel applied to porcine tissue**. Therefore, the pipeline does not attempt to interpret all human assay probes as porcine miRNAs.

Only features with a sequence-supported `Sus scrofa` counterpart are retained.

Consequently, the final porcine miRNA dataset represents the subset of miRNAs that can be confidently recovered from the human NanoString panel and should not be interpreted as the complete porcine miRNome. Pig-specific miRNAs lacking homologous probes on the human assay may not be detected.

---

## Citation

If this pipeline is used in a publication, please cite the associated manuscript when available.

Suggested repository citation format:

```text
Adepoju AE et al. NanoString miRNA analysis pipeline for cross-species
human-to-pig renal miRNA profiling in neonatal sepsis-associated acute
kidney injury.
```

The formal citation can be updated after publication.

---

## Data Availability

Raw NanoString RCC files should only be included in the public repository if their release is permitted by the study's data-sharing requirements.

For repositories where raw RCC files cannot be distributed, the repository may instead contain:

```text
data/
└── README.md
```

describing how authorized users can obtain the raw data.

Processed, non-identifiable output tables can be provided separately where appropriate.

---

## License

Add the license appropriate for your repository before public release.

Common options include:

- MIT License
- GNU General Public License v3.0
- BSD 3-Clause License

For an academic analysis repository intended primarily for reproducibility, select the license that best matches your institution's data and software-sharing requirements.

---

## Contact

For questions regarding the analysis workflow, please use the Issues section of the GitHub repository or provide the corresponding author's preferred contact information here.
