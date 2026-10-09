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
library(clusterProfiler)
library(org.Hs.eg.db)
library(ReactomePA)
library(multiMiR)


# ============================================================================
# PART 1 - Portable project paths
# ============================================================================

project_dir <- Sys.getenv(
  "NANOSTRING_PROJECT_DIR",
  unset = "."
)

project_dir <- normalizePath(
  project_dir,
  winslash = "/",
  mustWork = FALSE
)

data_dir <- file.path(project_dir, "data")
rcc_dir <- file.path(data_dir, "rcc")
reference_dir <- file.path(project_dir, "reference")
out_dir <- file.path(project_dir, "results")

metadata_file <- file.path(data_dir, "metadata.csv")

# Optional NanoString sample-column prefix. Blank means no prefix filtering.
sample_prefix <- Sys.getenv(
  "NANOSTRING_SAMPLE_PREFIX",
  unset = ""
)

dir.create(reference_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

if (!dir.exists(rcc_dir)) {
  stop(
    "RCC directory not found: ", rcc_dir, "\n",
    "Place the NanoString .RCC files in data/rcc/ or set NANOSTRING_PROJECT_DIR."
  )
}


# ============================================================================
# PART 2 - Build human-to-pig miRNA sequence mapping from miRBase
# ============================================================================

mature_url  <- "https://www.mirbase.org/download/CURRENT/mature.fa"
hairpin_url <- "https://www.mirbase.org/download/CURRENT/hairpin.fa"

mature_file  <- file.path(reference_dir, "mirbase_mature.fa")
hairpin_file <- file.path(reference_dir, "mirbase_hairpin.fa")

mapping_file <- file.path(
  reference_dir,
  "hsa_to_ssc_best_exact_sequence_match.csv"
)


# ----------------------------------------------------------------------------
# Robust miRBase FASTA downloader
# ----------------------------------------------------------------------------

download_mirbase_fasta <- function(url, destfile) {

  resp <- httr2::request(url) |>
    httr2::req_user_agent("R httr2 miRBase FASTA downloader") |>
    httr2::req_retry(max_tries = 3) |>
    httr2::req_perform()

  txt <- httr2::resp_body_string(resp)

  # miRBase may occasionally return FASTA wrapped in HTML.
  txt <- gsub("&gt;", ">", txt, fixed = TRUE)
  txt <- gsub("&lt;", "<", txt, fixed = TRUE)
  txt <- gsub("<br\\s*/?>", "\n", txt, perl = TRUE)
  txt <- gsub("</p>", "\n", txt, ignore.case = TRUE)
  txt <- gsub("<p>", "", txt, ignore.case = TRUE)
  txt <- gsub("<[^>]+>", "", txt)

  txt <- trimws(txt)

  lines <- unlist(strsplit(txt, "\n", fixed = TRUE))
  lines <- trimws(lines)
  lines <- lines[lines != ""]

  if (length(lines) == 0 || !startsWith(lines[1], ">")) {
    stop(
      "Downloaded content from ", url,
      " could not be converted to FASTA."
    )
  }

  writeLines(lines, destfile, useBytes = TRUE)

  invisible(destfile)
}


# Download only when local cached copies do not already exist.
if (!file.exists(mature_file)) {
  download_mirbase_fasta(mature_url, mature_file)
}

if (!file.exists(hairpin_file)) {
  download_mirbase_fasta(hairpin_url, hairpin_file)
}


# ----------------------------------------------------------------------------
# Parse mature and precursor miRNAs
# ----------------------------------------------------------------------------

parse_mirbase_fasta <- function(
  fasta_path,
  species_prefix,
  seq_type = c("mature", "precursor")
) {

  seq_type <- match.arg(seq_type)

  fa <- Biostrings::readRNAStringSet(fasta_path)
  hdr <- names(fa)

  keep <- grepl(paste0("^", species_prefix), hdr)
  fa <- fa[keep]
  hdr <- names(fa)

  first_space <- regexpr(" ", hdr)
  has_space <- first_space > 0

  mirna_name <- ifelse(
    has_space,
    substr(hdr, 1, first_space - 1),
    hdr
  )

  rest <- ifelse(
    has_space,
    substr(hdr, first_space + 1, nchar(hdr)),
    ""
  )

  second_space <- regexpr(" ", rest)
  has_second <- second_space > 0

  accession <- ifelse(
    has_second,
    substr(rest, 1, second_space - 1),
    rest
  )

  description <- ifelse(
    has_second,
    substr(rest, second_space + 1, nchar(rest)),
    ""
  )

  data.frame(
    miRNA = mirna_name,
    accession = accession,
    description = description,
    sequence = as.character(fa),
    type = seq_type,
    stringsAsFactors = FALSE
  )
}


human_mature <- parse_mirbase_fasta(
  mature_file,
  species_prefix = "hsa-",
  seq_type = "mature"
)

pig_mature <- parse_mirbase_fasta(
  mature_file,
  species_prefix = "ssc-",
  seq_type = "mature"
)

pig_precursor <- parse_mirbase_fasta(
  hairpin_file,
  species_prefix = "ssc-",
  seq_type = "precursor"
)

pig_all <- dplyr::bind_rows(
  pig_mature,
  pig_precursor
)

write.csv(
  human_mature,
  file.path(reference_dir, "miRBase_hsa_mature.csv"),
  row.names = FALSE
)

write.csv(
  pig_mature,
  file.path(reference_dir, "miRBase_ssc_mature.csv"),
  row.names = FALSE
)

write.csv(
  pig_precursor,
  file.path(reference_dir, "Pig_miRNA_precursor.csv"),
  row.names = FALSE
)

write.csv(
  pig_all,
  file.path(reference_dir, "Pig_miRNA_all.csv"),
  row.names = FALSE
)


# ----------------------------------------------------------------------------
# Exact human-to-pig mature-sequence matching
# ----------------------------------------------------------------------------

hsa_ssc_seq_match <- human_mature %>%
  dplyr::select(-type) %>%
  dplyr::inner_join(
    pig_mature %>% dplyr::select(-type),
    by = "sequence",
    suffix = c("_hsa", "_ssc"),
    relationship = "many-to-many"
  ) %>%
  dplyr::arrange(miRNA_hsa, miRNA_ssc)

write.csv(
  hsa_ssc_seq_match,
  file.path(reference_dir, "hsa_ssc_mature_exact_sequence_matches.csv"),
  row.names = FALSE
)


# ----------------------------------------------------------------------------
# Resolve multiple pig matches sharing the same mature sequence
# ----------------------------------------------------------------------------

choose_best_ssc_from_seq <- function(hsa_name, ssc_names) {

  ssc_names <- unique(ssc_names)

  exact_match <- sub("^hsa-", "ssc-", hsa_name)

  if (exact_match %in% ssc_names) {
    return(exact_match)
  }

  if (grepl("-3p$", hsa_name, ignore.case = TRUE)) {

    arm_hits <- ssc_names[
      grepl("-3p$", ssc_names, ignore.case = TRUE)
    ]

    if (length(arm_hits) > 0) {
      return(arm_hits[1])
    }
  }

  if (grepl("-5p$", hsa_name, ignore.case = TRUE)) {

    arm_hits <- ssc_names[
      grepl("-5p$", ssc_names, ignore.case = TRUE)
    ]

    if (length(arm_hits) > 0) {
      return(arm_hits[1])
    }
  }

  ssc_names[1]
}


hsa_to_ssc_best <- hsa_ssc_seq_match %>%
  dplyr::group_by(
    miRNA_hsa,
    accession_hsa,
    sequence
  ) %>%
  dplyr::summarise(
    miRNA_ssc_best = choose_best_ssc_from_seq(
      miRNA_hsa[1],
      miRNA_ssc
    ),
    accession_ssc_best = accession_ssc[
      match(miRNA_ssc_best, miRNA_ssc)
    ],
    all_ssc_matches = paste(
      unique(miRNA_ssc),
      collapse = "; "
    ),
    n_ssc_matches = dplyr::n_distinct(miRNA_ssc),
    .groups = "drop"
  )


write.csv(
  hsa_to_ssc_best,
  mapping_file,
  row.names = FALSE
)


# Create the mapping object expected by the NanoString analysis.
mir_map <- hsa_to_ssc_best %>%
  dplyr::transmute(
    hsa_miRNA = as.character(miRNA_hsa),
    ssc_miRNA = as.character(miRNA_ssc_best),
    hsa_accession = as.character(accession_hsa),
    ssc_accession = as.character(accession_ssc_best),
    sequence = as.character(sequence)
  ) %>%
  dplyr::filter(
    !is.na(hsa_miRNA),
    hsa_miRNA != "",
    !is.na(ssc_miRNA),
    ssc_miRNA != ""
  ) %>%
  dplyr::distinct(
    hsa_miRNA,
    .keep_all = TRUE
  )

write.csv(
  mir_map,
  file.path(
    reference_dir,
    "hsa_to_ssc_selected_mapping_from_sequence.csv"
  ),
  row.names = FALSE
)

valid_hsa <- unique(mir_map$hsa_miRNA)

message(
  "miRBase mapping complete: ",
  nrow(hsa_to_ssc_best),
  " human miRNA entries have an exact-sequence pig match."
)


# ============================================================================
# PART 3 - NanoString analysis
# ============================================================================

# ============================================================
# STEP 3 - Helpers for cleaning NanoString miRNA labels
# ============================================================

# Remove "|" and everything after it
strip_pipe_suffix <- function(x) {
  sub("\\|.*$", "", x)
}

# Split on "+" and keep the first candidate that exists in valid_hsa
choose_best_hsa_from_feature <- function(x, valid_hsa_vec) {
  x_clean <- strip_pipe_suffix(x)
  candidates <- trimws(unlist(strsplit(x_clean, "\\+")))
  candidates <- candidates[candidates != ""]
  
  in_map <- candidates[candidates %in% valid_hsa_vec]
  if (length(in_map) >= 1) return(in_map[1])
  
  return(NA_character_)
}


# ============================================================
# STEP 4 - Read NanoString RCC files
# ============================================================

rcc_obj <- nanostringr::read_rcc(path = rcc_dir)

counts_df <- as_tibble(rcc_obj[["raw"]])

if (ncol(counts_df) > 0 &&
    is.numeric(counts_df[[1]]) &&
    all(!is.na(counts_df[[1]])) &&
    all(counts_df[[1]] == seq_len(nrow(counts_df)))) {
  counts_df <- counts_df %>% dplyr::select(-1)
}

stopifnot(all(c("Code.Class", "Name", "Accession") %in% names(counts_df)))

id_cols <- c("Code.Class", "Name", "Accession")
sample_cols <- setdiff(names(counts_df), id_cols)

# Optional sample-name prefix. Leave blank to use all sample columns.
if (nzchar(sample_prefix)) {
  sample_cols <- sample_cols[grepl(paste0("^", sample_prefix), sample_cols)]
}

stopifnot(length(sample_cols) > 0)


# ============================================================
# STEP 5 - Convert wide counts to long format
# ============================================================

long <- counts_df %>%
  tidyr::pivot_longer(
    cols = dplyr::all_of(sample_cols),
    names_to = "SampleID",
    values_to = "Count"
  ) %>%
  dplyr::transmute(
    SampleID  = SampleID,
    CodeClass = as.character(.data[["Code.Class"]]),
    Target    = as.character(.data[["Name"]]),
    Accession = as.character(.data[["Accession"]]),
    Count     = as.numeric(Count)
  ) %>%
  dplyr::filter(!is.na(Count))


# ============================================================
# STEP 6 - Identify control probes and endogenous miRNAs
# ============================================================

neg <- long %>% dplyr::filter(grepl("NEG", Target, ignore.case = TRUE))
pos <- long %>% dplyr::filter(grepl("POS", Target, ignore.case = TRUE))
hk  <- long %>% dplyr::filter(tolower(CodeClass) %in% c("housekeeping", "hk"))

endo <- long %>%
  dplyr::filter(!(Target %in% unique(neg$Target))) %>%
  dplyr::filter(!(Target %in% unique(pos$Target))) %>%
  dplyr::filter(!(Target %in% unique(hk$Target)))


# ============================================================
# STEP 7 - Normalize the NanoString data
# ============================================================

# 7A. Background correction
bg_tbl <- neg %>%
  dplyr::group_by(SampleID) %>%
  dplyr::summarise(
    bg = mean(Count, na.rm = TRUE) + 2 * sd(Count, na.rm = TRUE),
    .groups = "drop"
  )

endo2 <- endo %>%
  dplyr::left_join(bg_tbl, by = "SampleID") %>%
  dplyr::mutate(
    bg = ifelse(is.na(bg), 0, bg),
    Count_bg = pmax(Count - bg, 1)
  )

# 7B. Positive-control normalization
pos_gm <- pos %>%
  dplyr::group_by(SampleID) %>%
  dplyr::summarise(
    pos_geo = exp(mean(log(pmax(Count, 1)), na.rm = TRUE)),
    .groups = "drop"
  )

pos_med <- median(pos_gm$pos_geo, na.rm = TRUE)

endo3 <- endo2 %>%
  dplyr::left_join(pos_gm, by = "SampleID") %>%
  dplyr::mutate(
    pos_sf = ifelse(is.na(pos_geo) | pos_geo == 0, 1, pos_med / pos_geo),
    Count_posnorm = Count_bg * pos_sf
  )

# 7C. Housekeeping normalization
if (nrow(hk) > 0) {
  hk_tbl <- hk %>%
    dplyr::group_by(SampleID) %>%
    dplyr::summarise(
      hk_geo = exp(mean(log(pmax(Count, 1)), na.rm = TRUE)),
      .groups = "drop"
    )
  
  hk_med <- median(hk_tbl$hk_geo, na.rm = TRUE)
  
  endo4 <- endo3 %>%
    dplyr::left_join(hk_tbl, by = "SampleID") %>%
    dplyr::mutate(
      hk_sf = ifelse(is.na(hk_geo) | hk_geo == 0, 1, hk_med / hk_geo),
      Count_norm = Count_posnorm * hk_sf
    )
} else {
  message("No housekeeping probes detected. Skipping HK normalization.")
  endo4 <- endo3 %>%
    dplyr::mutate(Count_norm = Count_posnorm)
}


# ============================================================
# STEP 8 - Build feature annotation
# original NanoString label -> cleaned hsa -> mapped ssc
# IMPORTANT: exclude unmapped hsa miRNAs
# ============================================================

feature_annot <- endo4 %>%
  dplyr::distinct(Target) %>%
  dplyr::mutate(
    hsa_clean = strip_pipe_suffix(Target),
    hsa_query_miRNA = vapply(Target, choose_best_hsa_from_feature, character(1), valid_hsa_vec = valid_hsa),
    ssc_miRNA = mir_map$ssc_miRNA[match(hsa_query_miRNA, mir_map$hsa_miRNA)]
  ) %>%
  dplyr::mutate(
    mapping_status = ifelse(!is.na(ssc_miRNA) & ssc_miRNA != "", "mapped_to_ssc", "unmapped_excluded")
  ) %>%
  dplyr::filter(!is.na(hsa_query_miRNA), hsa_query_miRNA != "") %>%
  dplyr::filter(!is.na(ssc_miRNA), ssc_miRNA != "") %>%
  dplyr::mutate(feature_id = ssc_miRNA)

write.csv(feature_annot, file.path(out_dir, "feature_annotation_hsa_to_ssc_sequence_based.csv"), row.names = FALSE)

mapping_summary <- data.frame(
  mapped_features = nrow(feature_annot),
  total_unique_targets = n_distinct(endo4$Target),
  excluded_unmapped = n_distinct(endo4$Target) - nrow(feature_annot)
)
write.csv(mapping_summary, file.path(out_dir, "mapping_summary.csv"), row.names = FALSE)


# ============================================================
# STEP 9 - Build normalized expression matrix using ONLY mapped ssc labels
# If multiple hsa probes collapse to same ssc label, average them
# ============================================================

expr_counts <- endo4 %>%
  dplyr::select(Target, SampleID, Count_norm) %>%
  dplyr::inner_join(feature_annot, by = "Target") %>%
  dplyr::select(feature_id, hsa_query_miRNA, ssc_miRNA, SampleID, Count_norm)

expr_counts_wide <- expr_counts %>%
  dplyr::select(feature_id, SampleID, Count_norm) %>%
  tidyr::pivot_wider(names_from = SampleID, values_from = Count_norm)

expr_counts_collapsed <- expr_counts_wide %>%
  dplyr::group_by(feature_id) %>%
  dplyr::summarise(
    dplyr::across(dplyr::all_of(sample_cols), ~ mean(.x, na.rm = TRUE)),
    .groups = "drop"
  )

expr_mat <- expr_counts_collapsed %>%
  tibble::column_to_rownames("feature_id") %>%
  as.matrix()

log2_mat <- log2(expr_mat)

write.csv(expr_mat, file.path(out_dir, "endogenous_normalized_counts_matrix_ssc_labels_mapped_only.csv"))
write.csv(log2_mat, file.path(out_dir, "endogenous_log2_matrix_ssc_labels_mapped_only.csv"))

feature_map_final <- expr_counts %>%
  dplyr::distinct(feature_id, hsa_query_miRNA, ssc_miRNA)

write.csv(feature_map_final, file.path(out_dir, "final_feature_id_mapping_mapped_only.csv"), row.names = FALSE)


# ============================================================
# STEP 10 - Metadata
# ============================================================

# Preferred GitHub/reproducible format:
#   data/metadata.csv
# with columns:
#   SampleID,Group
#
# Example:
#   SampleID,Group
#   sample_01,Control
#   sample_02,Control
#   sample_07,Sepsis
#
# If metadata.csv is absent, the script first attempts to infer groups from
# sample names containing Control/Sham or Sepsis/CLP. As a final fallback,
# for exactly 12 samples it assigns the first 6 to Control and the next 6
# to Sepsis to reproduce the original study layout.

if (file.exists(metadata_file)) {

  metadata <- read.csv(
    metadata_file,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  stopifnot(all(c("SampleID", "Group") %in% colnames(metadata)))

  metadata <- metadata[
    match(colnames(log2_mat), metadata$SampleID),
    ,
    drop = FALSE
  ]

  if (anyNA(metadata$SampleID)) {
    stop("metadata.csv does not contain all sample IDs present in the expression matrix.")
  }

} else {

  sample_ids <- colnames(log2_mat)

  inferred_group <- dplyr::case_when(
    grepl("control|sham", sample_ids, ignore.case = TRUE) ~ "Control",
    grepl("sepsis|clp", sample_ids, ignore.case = TRUE) ~ "Sepsis",
    TRUE ~ NA_character_
  )

  if (all(!is.na(inferred_group))) {

    metadata <- data.frame(
      SampleID = sample_ids,
      Group = inferred_group,
      stringsAsFactors = FALSE
    )

  } else if (length(sample_ids) == 12) {

    warning(
      "No metadata.csv found and groups could not be inferred from sample names. ",
      "Using the original study layout: first 6 samples = Control, next 6 = Sepsis. ",
      "For public reuse, provide data/metadata.csv."
    )

    metadata <- data.frame(
      SampleID = sample_ids,
      Group = c(rep("Control", 6), rep("Sepsis", 6)),
      stringsAsFactors = FALSE
    )

  } else {

    stop(
      "Group labels could not be inferred. Create data/metadata.csv ",
      "with columns SampleID and Group."
    )
  }
}

metadata$Group <- factor(metadata$Group, levels = c("Control", "Sepsis"))

if (anyNA(metadata$Group)) {
  stop("Group values must resolve to 'Control' or 'Sepsis'. Check data/metadata.csv.")
}

stopifnot(all(metadata$SampleID == colnames(log2_mat)))


# ============================================================
# STEP 11 - Differential expression with limma
# Row names are SSC labels only
# ============================================================

design <- model.matrix(~ 0 + Group, data = metadata)
colnames(design) <- levels(metadata$Group)

fit <- lmFit(log2_mat, design)

contrast.matrix <- makeContrasts(
  Sepsis_vs_Control = Sepsis - Control,
  levels = design
)

fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)

results <- topTable(
  fit2,
  coef = "Sepsis_vs_Control",
  number = Inf,
  adjust.method = "BH",
  sort.by = "P"
)

results <- as.data.frame(results, stringsAsFactors = FALSE)
results$ssc_miRNA <- rownames(results)
results$hsa_query_miRNA <- feature_map_final$hsa_query_miRNA[
  match(results$ssc_miRNA, feature_map_final$feature_id)
]

# remove any residual NA-mapped rows just in case
results <- results[!is.na(results$hsa_query_miRNA) & results$hsa_query_miRNA != "", , drop = FALSE]

write.csv(results, file.path(out_dir, "limma_results_ssc_labels_mapped_only.csv"), row.names = TRUE)


# ============================================================
# STEP 12 - QC summary
# ============================================================

qc_summary <- tibble::tibble(
  n_samples = ncol(log2_mat),
  n_targets_after_mapping = nrow(log2_mat),
  n_neg_targets = dplyr::n_distinct(neg$Target),
  n_pos_targets = dplyr::n_distinct(pos$Target),
  n_hk_targets  = dplyr::n_distinct(hk$Target)
)

readr::write_csv(qc_summary, file.path(out_dir, "qc_summary_ssc_labels_mapped_only.csv"))


# ============================================================
# STEP 13 - Heatmap
# ============================================================

sig <- results[!is.na(results$P.Value) & results$P.Value < 0.05, , drop = FALSE]

if (nrow(sig) > 1) {
  feat <- rownames(sig)[1:min(50, nrow(sig))]
  mat_hm <- log2_mat[feat, , drop = FALSE]
  
  mat_z <- t(scale(t(mat_hm)))
  mat_z[is.na(mat_z)] <- 0
  
  metadata$SampleNum <- as.numeric(sub(".*_(\\d+)$", "\\1", metadata$SampleID))
  metadata_ord <- metadata[order(metadata$Group, metadata$SampleNum), , drop = FALSE]
  mat_z <- mat_z[, metadata_ord$SampleID, drop = FALSE]
  colnames(mat_z) <- sprintf("%02d", metadata_ord$SampleNum)
  
  # make sure group order is consistent
  metadata_ord$Group <- factor(metadata_ord$Group, levels = c("Control", "Sepsis"))
  
  ha <- HeatmapAnnotation(
    Group = metadata_ord$Group,
    col = list(
      Group = c(
        "Control" = "blue",
        "Sepsis" = "red"
      )
    ),
    annotation_name_side = "left",
    annotation_legend_param = list(
      Group = list(title = "Group")
    )
  )
  
  col_fun <- colorRamp2(c(-2, 0, 2), c("blue", "white", "red"))
  
  ht <- Heatmap(
    mat_z,
    name = "z-score",
    top_annotation = ha,
    col = col_fun,
    show_row_names = TRUE,
    show_column_names = TRUE,
    cluster_rows = TRUE,
    cluster_columns = TRUE,
    column_split = metadata_ord$Group,
    row_names_gp = gpar(fontsize = 12),
    column_names_gp = gpar(fontsize = 12),
    column_names_rot = 90,
    column_title = NULL
  )
  
  # Save PDF
  pdf(
    file.path(out_dir, "ComplexHeatmap_grouped_ordered_samples_ssc_mapped_only.pdf"),
    width = 10,
    height = 8
  )
  draw(
    ht,
    heatmap_legend_side = "right",
    annotation_legend_side = "right"
  )
  dev.off()
  
  # Save PNG
  png(
    filename = file.path(out_dir, "ComplexHeatmap_grouped_ordered_samples_ssc_mapped_only.png"),
    width = 2000,
    height = 1600,
    res = 300
  )
  draw(
    ht,
    heatmap_legend_side = "right",
    annotation_legend_side = "right"
  )
  dev.off()
}


# ============================================================
# STEP 14 - Volcano plot
# ============================================================

p_cut <- 0.05
fc_cut <- 0.7

vol <- as.data.frame(results, stringsAsFactors = FALSE)
vol$miRNA <- vol$ssc_miRNA
vol$P.Value <- as.numeric(vol$P.Value)
vol$logFC <- as.numeric(vol$logFC)
vol$negLog10P <- -log10(pmax(vol$P.Value, 1e-300))

# Define groups
vol$ColorClass <- ifelse(
  vol$P.Value < p_cut & vol$logFC >= fc_cut,
  "Up (Sepsis)",
  ifelse(
    vol$P.Value < p_cut & vol$logFC <= -fc_cut,
    "Down (Sepsis)",
    "Not sig"
  )
)

# Count significant up/down miRNAs
n_up <- sum(vol$P.Value < p_cut & vol$logFC >= fc_cut, na.rm = TRUE)
n_down <- sum(vol$P.Value < p_cut & vol$logFC <= -fc_cut, na.rm = TRUE)

# Select top 10 upregulated
lab_up <- vol[
  vol$P.Value < p_cut & vol$logFC >= fc_cut,
  ,
  drop = FALSE
]
lab_up <- lab_up[order(lab_up$P.Value, -lab_up$logFC, na.last = NA), , drop = FALSE]
lab_up <- lab_up[1:min(10, nrow(lab_up)), , drop = FALSE]

# Select top 10 downregulated
lab_down <- vol[
  vol$P.Value < p_cut & vol$logFC <= -fc_cut,
  ,
  drop = FALSE
]
lab_down <- lab_down[order(lab_down$P.Value, lab_down$logFC, na.last = NA), , drop = FALSE]
lab_down <- lab_down[1:min(10, nrow(lab_down)), , drop = FALSE]

# Combine labels
lab <- rbind(lab_up, lab_down)

p <- ggplot(vol, aes(x = logFC, y = negLog10P)) +
  geom_point(
    aes(color = ColorClass),
    alpha = 0.9,
    size = 3.8
  ) +
  geom_vline(xintercept = c(-fc_cut, fc_cut), linetype = "dashed", linewidth = 0.7) +
  geom_hline(yintercept = -log10(p_cut), linetype = "dashed", linewidth = 0.7) +
  geom_text_repel(
    data = lab,
    aes(label = miRNA),
    size = 3.5,
    box.padding = 0.4,
    point.padding = 0.3,
    max.overlaps = Inf
  ) +
  scale_color_manual(
    values = c(
      "Up (Sepsis)" = "red",
      "Down (Sepsis)" = "blue",
      "Not sig" = "grey70"
    ),
    name = "Direction"
  ) +
  theme_minimal(base_size = 14) +
  labs(
    title = paste0(
      "Volcano plot (Sepsis vs Control)\n",
      "Upregulated: ", n_up, "   |   Downregulated: ", n_down
    ),
    x = "log2 Fold Change",
    y = "-log10(P.Value)"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.title = element_text(face = "bold"),
    legend.title = element_text(face = "bold")
  )

ggsave(
  filename = file.path(out_dir, "volcano_limma_top10_up_top10_down.png"),
  plot = p,
  width = 9,
  height = 7,
  dpi = 300
)

# ============================================================
# STEP 15 - Significant miRNAs for target analysis
# Use matched HUMAN miRNAs only
# ============================================================

sig_miRNAs <- results %>%
  dplyr::filter(P.Value < 0.05) %>%
  dplyr::filter(!is.na(hsa_query_miRNA), hsa_query_miRNA != "") %>%
  dplyr::arrange(P.Value)

write.csv(sig_miRNAs, file.path(out_dir, "significant_ssc_miRNAs_rawP_lt_0.05_mapped_only.csv"), row.names = TRUE)

miRNA_list_hsa <- unique(sig_miRNAs$hsa_query_miRNA)
miRNA_list_hsa <- miRNA_list_hsa[!is.na(miRNA_list_hsa) & miRNA_list_hsa != ""]


# ============================================================
# STEP 16 - Retrieve target genes with multiMiR
# ============================================================

if (length(miRNA_list_hsa) > 0) {
  
  mm <- multiMiR::get_multimir(
    mirna = miRNA_list_hsa,
    table = "validated",
    summary = FALSE,
    legacy.out = FALSE,
    org = "hsa"
  )
  
  target_tbl <- mm@data
  
  keep_cols <- intersect(
    c("mature_mirna_id", "target_symbol", "target_entrez", "database", "support_type", "experiment"),
    colnames(target_tbl)
  )
  
  target_tbl <- target_tbl[, keep_cols, drop = FALSE]
  
  if ("target_symbol" %in% colnames(target_tbl)) {
    target_tbl <- target_tbl[!is.na(target_tbl$target_symbol) & target_tbl$target_symbol != "", , drop = FALSE]
  }
  
  target_tbl$ssc_miRNA <- sig_miRNAs$ssc_miRNA[match(target_tbl$mature_mirna_id, sig_miRNAs$hsa_query_miRNA)]
  target_tbl <- target_tbl[!is.na(target_tbl$ssc_miRNA) & target_tbl$ssc_miRNA != "", , drop = FALSE]
  
  write.csv(target_tbl, file.path(out_dir, "validated_miRNA_targets_human_query_with_ssc_labels_mapped_only.csv"), row.names = FALSE)
  
  target_genes_human <- unique(target_tbl$target_symbol)
  target_genes_human <- target_genes_human[!is.na(target_genes_human) & target_genes_human != ""]
  
  write.csv(
    data.frame(target_gene_human = target_genes_human),
    file.path(out_dir, "unique_target_genes_human_mapped_only.csv"),
    row.names = FALSE
  )
  
} else {
  target_tbl <- data.frame()
  target_genes_human <- character(0)
}


# ============================================================
# STEP 17 - Convert HUMAN target genes to Entrez IDs
# ============================================================

if (length(target_genes_human) > 0) {
  
  gene_map <- clusterProfiler::bitr(
    target_genes_human,
    fromType = "SYMBOL",
    toType = c("ENTREZID"),
    OrgDb = org.Hs.eg.db
  )
  
  write.csv(gene_map, file.path(out_dir, "target_gene_ID_mapping_human_mapped_only.csv"), row.names = FALSE)
  
  entrez_ids <- unique(gene_map$ENTREZID)
  
} else {
  gene_map <- data.frame()
  entrez_ids <- character(0)
}


# ============================================================
# STEP 18 - GO Biological Process enrichment
# ============================================================

if (length(entrez_ids) > 0) {
  
  ego_bp <- enrichGO(
    gene = entrez_ids,
    OrgDb = org.Hs.eg.db,
    keyType = "ENTREZID",
    ont = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.01,
    readable = TRUE
  )
  
  ego_bp_df <- as.data.frame(ego_bp)
  write.csv(ego_bp_df, file.path(out_dir, "GO_BP_enrichment_human_targets_mapped_only.csv"), row.names = FALSE)
  
  if (nrow(ego_bp_df) > 0) {
    pdf(file.path(out_dir, "GO_BP_dotplot_human_targets_mapped_only.pdf"), width = 10, height = 7)
    print(dotplot(ego_bp, showCategory = 20, title = "GO Biological Process"))
    dev.off()
  }
  
} else {
  ego_bp <- NULL
}


# ============================================================
# STEP 19 - KEGG pathway enrichment
# ============================================================

if (length(entrez_ids) > 0) {
  
  ekegg <- enrichKEGG(
    gene = entrez_ids,
    organism = "hsa",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.01,
    pAdjustMethod = "BH"
  )
  
  ekegg_df <- as.data.frame(ekegg)
  write.csv(ekegg_df, file.path(out_dir, "KEGG_enrichment_human_targets_mapped_only.csv"), row.names = FALSE)
}

  library(dplyr)
  library(ggplot2)
  
  # --------------------------------------------------
  # 1. Start from your enrichment results table
  # Example: ego_bp_df, ekegg_df, ereact_df, or similar
  # --------------------------------------------------
  enrich_df <- ekegg_df   # change this to your enrichment dataframe
  
  # --------------------------------------------------
  # 2. Detect the subcategory column automatically
  # --------------------------------------------------
  subcat_col <- c("subcategory", "Subcategory", "Category")
  subcat_col <- subcat_col[subcat_col %in% colnames(enrich_df)][1]
  
  if (is.na(subcat_col) || length(subcat_col) == 0) {
    stop("No subcategory column found. Add a column such as 'subcategory' or 'Category' first.")
  }
  
  # --------------------------------------------------
  # 3. Make sure FoldEnrichment exists
  # If not, calculate it from GeneRatio and BgRatio
  # --------------------------------------------------
  if (!"FoldEnrichment" %in% colnames(enrich_df)) {
    
    ratio_to_num <- function(x) {
      sapply(strsplit(as.character(x), "/"), function(z) as.numeric(z[1]) / as.numeric(z[2]))
    }
    
    enrich_df$GeneRatio_num <- ratio_to_num(enrich_df$GeneRatio)
    enrich_df$BgRatio_num   <- ratio_to_num(enrich_df$BgRatio)
    enrich_df$FoldEnrichment <- enrich_df$GeneRatio_num / enrich_df$BgRatio_num
  }
  
  # --------------------------------------------------
  # 4. Filter for FoldEnrichment >= 1.3
  # --------------------------------------------------
  enrich_df_filt <- enrich_df %>%
    filter(!is.na(FoldEnrichment), FoldEnrichment >= 1.3)
  
  # --------------------------------------------------
  # 5. Aggregate by subcategory
  # - mean FoldEnrichment
  # - total gene counts across terms
  # - number of descriptions/terms in subcategory
  # --------------------------------------------------
  subcat_summary <- enrich_df_filt %>%
    group_by(.data[[subcat_col]]) %>%
    summarise(
      MeanFoldEnrichment = mean(FoldEnrichment, na.rm = TRUE),
      TotalGeneCount = sum(Count, na.rm = TRUE),
      N_Descriptions = n(),
      .groups = "drop"
    ) %>%
    arrange(MeanFoldEnrichment)
  
  # --------------------------------------------------
  # 6. Plot
  # --------------------------------------------------
  p <- ggplot(subcat_summary,
              aes(x = MeanFoldEnrichment,
                  y = reorder(.data[[subcat_col]], MeanFoldEnrichment))) +
    geom_point(aes(size = N_Descriptions, color = TotalGeneCount), alpha = 0.9) +
    
    geom_vline(xintercept = 1.3, linetype = "dashed", linewidth = 0.7) +
    
    scale_x_continuous(
      limits = c(1.25, NA),
      expand = expansion(mult = c(0, 0.05))
    ) +
    
    scale_size_continuous(name = "No. of Pathways") +
    scale_color_gradient(low = "lightblue", high = "red", name = "Total Gene Count") +
    
    theme_minimal(base_size = 14) +
    labs(
      title = "Enriched Subcategories",
      x = "Fold Enrichment",
      y = "Subcategory"
    ) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      axis.title = element_text(face = "bold")
    )
  
  print(p)
  
  # --------------------------------------------------
  # 7. Save outputs
  # --------------------------------------------------
  write.csv(subcat_summary,
            file.path(out_dir, "enrichment_subcategory_summary.csv"),
            row.names = FALSE)
  
  ggsave(
    file.path(out_dir, "enrichment_subcategory_bubbleplot.png"),
    plot = p,
    width = 10,
    height = 7,
    dpi = 300
  )

# ============================================================
# STEP 20 - Reactome pathway enrichment
# ============================================================

if (length(entrez_ids) > 0) {
  
  ereact <- enrichPathway(
    gene = entrez_ids,
    organism = "human",
    pvalueCutoff = 0.05,
    pAdjustMethod = "BH",
    readable = TRUE
  )
  
  ereact_df <- as.data.frame(ereact)
  write.csv(ereact_df, file.path(out_dir, "Reactome_enrichment_human_targets_mapped_only.csv"), row.names = FALSE)
  
  if (nrow(ereact_df) > 0) {
    pdf(file.path(out_dir, "Reactome_dotplot_human_targets_mapped_only.pdf"), width = 10, height = 7)
    print(dotplot(ereact, showCategory = 20, title = "Reactome pathways"))
    dev.off()
  }
  
} else {
  ereact <- NULL
}


# ============================================================
# STEP 21 - Cytoscape edge table
# ============================================================

if (exists("target_tbl") && nrow(target_tbl) > 0) {
  
  stopifnot("mature_mirna_id" %in% colnames(target_tbl))
  stopifnot("target_symbol" %in% colnames(target_tbl))
  
  edge_tbl <- unique(
    data.frame(
      miRNA = as.character(target_tbl[["ssc_miRNA"]]),
      TargetGene = as.character(target_tbl[["target_symbol"]]),
      hsa_query_miRNA = as.character(target_tbl[["mature_mirna_id"]]),
      stringsAsFactors = FALSE
    )
  )
  
  edge_tbl <- edge_tbl[
    !is.na(edge_tbl$miRNA) & edge_tbl$miRNA != "" &
      !is.na(edge_tbl$TargetGene) & edge_tbl$TargetGene != "",
    ,
    drop = FALSE
  ]
  
  write.csv(
    edge_tbl,
    file.path(out_dir, "miRNA_target_edges_for_Cytoscape_ssc_mapped_only.csv"),
    row.names = FALSE
  )
  
  mirna_gene_counts <- aggregate(
    TargetGene ~ miRNA,
    data = edge_tbl,
    FUN = function(x) length(unique(x))
  )
  colnames(mirna_gene_counts)[colnames(mirna_gene_counts) == "TargetGene"] <- "n_targets"
  mirna_gene_counts <- mirna_gene_counts[order(-mirna_gene_counts$n_targets), , drop = FALSE]
  
  write.csv(
    mirna_gene_counts,
    file.path(out_dir, "miRNA_target_gene_counts_ssc_mapped_only.csv"),
    row.names = FALSE
  )
}


# ============================================================
# STEP 22 - Cytoscape node table
# ============================================================

if (exists("edge_tbl") && nrow(edge_tbl) > 0) {
  
  mirna_nodes <- results %>%
    dplyr::transmute(
      node = ssc_miRNA,
      type = "miRNA",
      logFC = logFC,
      P.Value = P.Value,
      adj.P.Val = adj.P.Val,
      hsa_query_miRNA = hsa_query_miRNA
    )
  
  target_counts <- aggregate(
    TargetGene ~ miRNA,
    data = edge_tbl,
    FUN = function(x) length(unique(x))
  )
  colnames(target_counts) <- c("node", "target_count")
  
  mirna_nodes <- merge(mirna_nodes, target_counts, by = "node", all.x = TRUE)
  mirna_nodes$target_count[is.na(mirna_nodes$target_count)] <- 0
  
  gene_counts <- aggregate(
    miRNA ~ TargetGene,
    data = edge_tbl,
    FUN = function(x) length(unique(x))
  )
  colnames(gene_counts) <- c("node", "regulator_count")
  
  gene_nodes <- data.frame(
    node = gene_counts$node,
    type = "Gene",
    logFC = NA,
    P.Value = NA,
    adj.P.Val = NA,
    hsa_query_miRNA = NA,
    target_count = NA,
    regulator_count = gene_counts$regulator_count,
    stringsAsFactors = FALSE
  )
  
  mirna_nodes$regulator_count <- NA
  
  mirna_nodes <- mirna_nodes[, c(
    "node", "type", "logFC", "P.Value", "adj.P.Val",
    "hsa_query_miRNA", "target_count", "regulator_count"
  )]
  
  node_tbl <- rbind(mirna_nodes, gene_nodes)
  
  write.csv(
    node_tbl,
    file.path(out_dir, "cytoscape_nodes_ssc_mapped_only.csv"),
    row.names = FALSE
  )
}


# ============================================================
# STEP 23 - Save session summary
# ============================================================

capture.output(sessionInfo(), file = file.path(out_dir, "sessionInfo_mapped_only.txt"))

message("Pipeline complete. All outputs saved to: ", out_dir)



  
#CHAPTER 2
  # ============================================================
  # Integrated miRNA target ranking + Cytoscape export
  # ============================================================
  
  library(dplyr)
  
  # ------------------------------------------------------------
  # STEP 1 - Check required input
  # ------------------------------------------------------------
  stopifnot(exists("results"))
  stopifnot(exists("target_tbl"))
  stopifnot(exists("out_dir"))
  
  stopifnot(all(c("ssc_miRNA", "hsa_query_miRNA", "logFC", "P.Value", "adj.P.Val") %in% colnames(results)))
  stopifnot(all(c("mature_mirna_id", "target_symbol", "ssc_miRNA") %in% colnames(target_tbl)))
  
  # ------------------------------------------------------------
  # STEP 2 - Clean target table
  # ------------------------------------------------------------
  target_tbl2 <- target_tbl %>%
    dplyr::filter(!is.na(ssc_miRNA), ssc_miRNA != "",
                  !is.na(mature_mirna_id), mature_mirna_id != "",
                  !is.na(target_symbol), target_symbol != "") %>%
    dplyr::distinct()
  
  # optional: inspect available columns
  print(colnames(target_tbl2))
  
  # ------------------------------------------------------------
  # STEP 3 - Summarize support across databases
  # ------------------------------------------------------------
  target_support <- target_tbl2 %>%
    dplyr::group_by(ssc_miRNA, mature_mirna_id, target_symbol) %>%
    dplyr::summarise(
      n_databases = dplyr::n_distinct(database),
      databases = paste(sort(unique(database)), collapse = "; "),
      n_records = dplyr::n(),
      .groups = "drop"
    )
  
  # ------------------------------------------------------------
  # STEP 4 - Add validation evidence
  # Here, because target_tbl came from multiMiR validated query in your pipeline,
  # all entries can be considered validated. If later you combine predicted+validated,
  # this still works.
  # ------------------------------------------------------------
  validated_db_keywords <- c("mirtarbase", "tarbase", "validated")
  
  target_support <- target_support %>%
    dplyr::mutate(
      validated = ifelse(
        grepl(paste(validated_db_keywords, collapse = "|"),
              tolower(databases)),
        1, 1   # set to 1 because your current target_tbl comes from validated table
      )
    )
  
  # ------------------------------------------------------------
  # STEP 5 - Add miRNA differential expression metrics
  # ------------------------------------------------------------
  mirna_info <- results %>%
    dplyr::select(ssc_miRNA, hsa_query_miRNA, logFC, P.Value, adj.P.Val) %>%
    dplyr::distinct()
  
  target_rank <- target_support %>%
    dplyr::left_join(mirna_info, by = "ssc_miRNA")
  
  # ------------------------------------------------------------
  # STEP 6 - Optional: add gene-level expression evidence
  # If you have gene_results with columns:
  #   gene_symbol, logFC_gene, P.Value_gene, adj.P.Val_gene
  # then this section will add opposite-direction scoring.
  # If not available, gene evidence is set to NA / 0.
  # ------------------------------------------------------------
  if (exists("gene_results")) {
    
    # expected columns: gene_symbol, logFC_gene or logFC
    gene_results2 <- gene_results
    
    if (!"gene_symbol" %in% colnames(gene_results2) && "target_symbol" %in% colnames(gene_results2)) {
      gene_results2$gene_symbol <- gene_results2$target_symbol
    }
    
    if (!"gene_symbol" %in% colnames(gene_results2)) {
      stop("gene_results exists, but no gene_symbol column was found.")
    }
    
    if (!"logFC_gene" %in% colnames(gene_results2) && "logFC" %in% colnames(gene_results2)) {
      gene_results2$logFC_gene <- gene_results2$logFC
    }
    
    target_rank <- target_rank %>%
      dplyr::left_join(
        gene_results2 %>%
          dplyr::select(gene_symbol, logFC_gene),
        by = c("target_symbol" = "gene_symbol")
      ) %>%
      dplyr::mutate(
        opposite_direction = ifelse(
          !is.na(logFC) & !is.na(logFC_gene) &
            ((logFC > 0 & logFC_gene < 0) | (logFC < 0 & logFC_gene > 0)),
          1, 0
        )
      )
    
  } else {
    target_rank$logFC_gene <- NA_real_
    target_rank$opposite_direction <- 0
  }
  
  # ------------------------------------------------------------
  # STEP 7 - Optional: add correlation evidence
  # If you have:
  #   miRNA_mat  = matrix of miRNA expression (rows = ssc_miRNA, cols = samples)
  #   gene_mat   = matrix of gene expression (rows = gene symbols, cols = samples)
  # then this section computes pairwise correlations.
  # Otherwise correlation is set to NA / 0.
  # ------------------------------------------------------------
  if (exists("miRNA_mat") && exists("gene_mat")) {
    
    get_pair_cor <- function(mir, gene, miRNA_mat, gene_mat) {
      if (!(mir %in% rownames(miRNA_mat))) return(NA_real_)
      if (!(gene %in% rownames(gene_mat))) return(NA_real_)
      suppressWarnings(
        cor(
          as.numeric(miRNA_mat[mir, ]),
          as.numeric(gene_mat[gene, ]),
          method = "pearson",
          use = "pairwise.complete.obs"
        )
      )
    }
    
    target_rank$correlation <- mapply(
      get_pair_cor,
      target_rank$ssc_miRNA,
      target_rank$target_symbol,
      MoreArgs = list(miRNA_mat = miRNA_mat, gene_mat = gene_mat)
    )
    
  } else {
    target_rank$correlation <- NA_real_
  }
  
  # ------------------------------------------------------------
  # STEP 8 - Optional: add pathway relevance
  # If enrichment target genes are available, mark them.
  # This uses any of:
  #   target_genes_human, kegg_df, ego_bp_df, ereact_df
  # ------------------------------------------------------------
  pathway_gene_set <- character(0)
  
  if (exists("target_genes_human")) {
    pathway_gene_set <- unique(c(pathway_gene_set, target_genes_human))
  }
  
  if (exists("kegg_df") && "geneID" %in% colnames(kegg_df)) {
    kg <- unique(unlist(strsplit(kegg_df$geneID, "/")))
    pathway_gene_set <- unique(c(pathway_gene_set, kg))
  }
  
  target_rank$pathway_relevant <- ifelse(
    target_rank$target_symbol %in% pathway_gene_set,
    1, 0
  )
  
  # ------------------------------------------------------------
  # STEP 9 - Composite score
  # You can tune these weights.
  # ------------------------------------------------------------
  target_rank <- target_rank %>%
    dplyr::mutate(
      score_validated = validated * 5,
      score_databases = pmin(n_databases, 5),
      score_correlation = ifelse(
        is.na(correlation), 0,
        ifelse(correlation <= -0.6, 3,
               ifelse(correlation <= -0.4, 1, 0))
      ),
      score_direction = opposite_direction * 2,
      score_pathway = pathway_relevant * 2,
      target_score = score_validated + score_databases + score_correlation + score_direction + score_pathway
    ) %>%
    dplyr::arrange(dplyr::desc(target_score), P.Value)
  
  # ------------------------------------------------------------
  # STEP 10 - Save ranked target table
  # ------------------------------------------------------------
  write.csv(
    target_rank,
    file.path(out_dir, "ranked_miRNA_targets_integrated.csv"),
    row.names = FALSE
  )
  
  # top 10 targets per miRNA
  top_targets_per_miRNA <- target_rank %>%
    dplyr::group_by(ssc_miRNA) %>%
    dplyr::slice_max(order_by = target_score, n = 10, with_ties = FALSE) %>%
    dplyr::ungroup()
  
  write.csv(
    top_targets_per_miRNA,
    file.path(out_dir, "top10_targets_per_miRNA_ranked.csv"),
    row.names = FALSE
  )
  
  # ------------------------------------------------------------
  # STEP 11 - Build Cytoscape edge table
  # Score is attached to each edge
  # ------------------------------------------------------------
  cyto_edges <- target_rank %>%
    dplyr::transmute(
      source = ssc_miRNA,
      target = target_symbol,
      interaction = "miRNA-target",
      hsa_query_miRNA = mature_mirna_id,
      n_databases = n_databases,
      databases = databases,
      validated = validated,
      correlation = correlation,
      opposite_direction = opposite_direction,
      pathway_relevant = pathway_relevant,
      target_score = target_score
    ) %>%
    dplyr::distinct()
  
  write.csv(
    cyto_edges,
    file.path(out_dir, "cytoscape_edges_ranked.csv"),
    row.names = FALSE
  )
  
  # ------------------------------------------------------------
  # STEP 12 - Build Cytoscape node table
  # miRNA nodes carry logFC
  # gene nodes carry regulator counts
  # ------------------------------------------------------------
  
  # miRNA node counts
  miRNA_target_counts <- cyto_edges %>%
    dplyr::group_by(source) %>%
    dplyr::summarise(
      target_count = dplyr::n_distinct(target),
      max_target_score = max(target_score, na.rm = TRUE),
      .groups = "drop"
    )
  
  miRNA_nodes <- results %>%
    dplyr::transmute(
      id = ssc_miRNA,
      type = "miRNA",
      logFC = logFC,
      P.Value = P.Value,
      adj.P.Val = adj.P.Val,
      hsa_query_miRNA = hsa_query_miRNA
    ) %>%
    dplyr::left_join(miRNA_target_counts, by = c("id" = "source")) %>%
    dplyr::mutate(
      target_count = ifelse(is.na(target_count), 0, target_count),
      max_target_score = ifelse(is.na(max_target_score), 0, max_target_score),
      regulator_count = NA_real_
    )
  
  # gene node counts
  gene_nodes <- cyto_edges %>%
    dplyr::group_by(target) %>%
    dplyr::summarise(
      regulator_count = dplyr::n_distinct(source),
      max_target_score = max(target_score, na.rm = TRUE),
      mean_target_score = mean(target_score, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::transmute(
      id = target,
      type = "Gene",
      logFC = if ("logFC_gene" %in% colnames(target_rank)) NA_real_ else NA_real_,
      P.Value = NA_real_,
      adj.P.Val = NA_real_,
      hsa_query_miRNA = NA_character_,
      target_count = NA_real_,
      regulator_count = regulator_count,
      max_target_score = max_target_score
    )
  
  cyto_nodes <- dplyr::bind_rows(
    miRNA_nodes %>%
      dplyr::select(id, type, logFC, P.Value, adj.P.Val, hsa_query_miRNA,
                    target_count, regulator_count, max_target_score),
    gene_nodes %>%
      dplyr::select(id, type, logFC, P.Value, adj.P.Val, hsa_query_miRNA,
                    target_count, regulator_count, max_target_score)
  )
  
  write.csv(
    cyto_nodes,
    file.path(out_dir, "cytoscape_nodes_ranked.csv"),
    row.names = FALSE
  )
  
  # ------------------------------------------------------------
  # STEP 13 - Cytoscape style guide export
  # Useful notes for mapping inside Cytoscape
  # ------------------------------------------------------------
  style_guide <- data.frame(
    table = c("Edges", "Edges", "Nodes", "Nodes", "Nodes"),
    attribute = c("target_score", "validated", "logFC", "type", "max_target_score"),
    recommended_visual = c(
      "Edge color or width",
      "Edge line style",
      "Node fill color",
      "Node shape",
      "Node border width or color"
    ),
    suggested_mapping = c(
      "Low=light grey, High=red",
      "0=dashed, 1=solid",
      "Negative=blue, Positive=red",
      "miRNA=diamond, Gene=ellipse",
      "Low=thin border, High=thick border"
    ),
    stringsAsFactors = FALSE
  )
  
  write.csv(
    style_guide,
    file.path(out_dir, "cytoscape_style_guide_ranked.csv"),
    row.names = FALSE
  )
  
message("Integrated miRNA target ranking and Cytoscape export complete.")
  
  
  # ============================================================
  # Volcano + network integration
  # Requires:
  #   results
  #   target_rank
  #   out_dir
  # ============================================================
  
  # ----------------------------
  # 0. Install/load packages
  # ----------------------------
  cran_pkgs <- c("dplyr", "ggplot2", "ggrepel", "igraph", "ggraph", "tidyr")
  to_install <- cran_pkgs[!cran_pkgs %in% rownames(installed.packages())]
  if (length(to_install) > 0) install.packages(to_install)
  
  library(dplyr)
  library(ggplot2)
  library(ggrepel)
  library(igraph)
  library(ggraph)
  library(tidyr)
  
  # ----------------------------
  # 1. Prepare volcano table
  # ----------------------------
  p_cut <- 0.05
  fc_cut <- 0.7
  
  vol <- as.data.frame(results, stringsAsFactors = FALSE)
  vol$miRNA <- vol$ssc_miRNA
  vol$P.Value <- as.numeric(vol$P.Value)
  vol$logFC <- as.numeric(vol$logFC)
  vol$adj.P.Val <- as.numeric(vol$adj.P.Val)
  vol$negLog10P <- -log10(pmax(vol$P.Value, 1e-300))
  
  vol$ColorClass <- ifelse(
    vol$P.Value < p_cut & vol$logFC >= fc_cut,
    "Up (Sepsis)",
    ifelse(
      vol$P.Value < p_cut & vol$logFC <= -fc_cut,
      "Down (Sepsis)",
      "Not sig"
    )
  )
  
  # ----------------------------
  # 2. Select miRNAs to carry into the network
  # Strategy:
  #   top 10 upregulated + top 10 downregulated by volcano criteria
  # ----------------------------
  vol_up <- vol %>%
    dplyr::filter(P.Value < p_cut, logFC >= fc_cut) %>%
    dplyr::arrange(P.Value, dplyr::desc(logFC)) %>%
    dplyr::slice(1:min(10, n()))
  
  vol_down <- vol %>%
    dplyr::filter(P.Value < p_cut, logFC <= -fc_cut) %>%
    dplyr::arrange(P.Value, logFC) %>%
    dplyr::slice(1:min(10, n()))
  
  vol_net_miRNAs <- dplyr::bind_rows(vol_up, vol_down) %>%
    dplyr::distinct(miRNA, .keep_all = TRUE)
  
  network_miRNAs <- unique(vol_net_miRNAs$miRNA)
  
  # ----------------------------
  # 3. Label these miRNAs on the volcano
  # ----------------------------
  vol$NetworkHighlight <- ifelse(vol$miRNA %in% network_miRNAs, "In network", "Background")
  
  n_up <- sum(vol$P.Value < p_cut & vol$logFC >= fc_cut, na.rm = TRUE)
  n_down <- sum(vol$P.Value < p_cut & vol$logFC <= -fc_cut, na.rm = TRUE)
  
  p_volcano <- ggplot(vol, aes(x = logFC, y = negLog10P)) +
    geom_point(
      aes(color = ColorClass),
      alpha = 0.85,
      size = 3.2
    ) +
    geom_point(
      data = vol %>% dplyr::filter(miRNA %in% network_miRNAs),
      shape = 21,
      stroke = 0.9,
      size = 4.8,
      fill = NA,
      color = "black"
    ) +
    geom_vline(xintercept = c(-fc_cut, fc_cut), linetype = "dashed", linewidth = 0.7) +
    geom_hline(yintercept = -log10(p_cut), linetype = "dashed", linewidth = 0.7) +
    geom_text_repel(
      data = vol %>% dplyr::filter(miRNA %in% network_miRNAs),
      aes(label = miRNA),
      size = 3.4,
      max.overlaps = Inf,
      box.padding = 0.4,
      point.padding = 0.3
    ) +
    scale_color_manual(
      values = c(
        "Up (Sepsis)" = "red",
        "Down (Sepsis)" = "blue",
        "Not sig" = "grey70"
      ),
      name = "Direction"
    ) +
    theme_minimal(base_size = 14) +
    labs(
      title = paste0(
        "Volcano plot with network-selected miRNAs\n",
        "Upregulated: ", n_up, "   |   Downregulated: ", n_down
      ),
      x = "log2 Fold Change",
      y = "-log10(P.Value)"
    ) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      axis.title = element_text(face = "bold"),
      legend.title = element_text(face = "bold")
    )
  
  # ----------------------------
  # 4. Build a focused network
  # Keep top 5 targets per selected miRNA by target_score
  # ----------------------------
  stopifnot(exists("target_rank"))
  
  network_edges <- target_rank %>%
    dplyr::filter(ssc_miRNA %in% network_miRNAs) %>%
    dplyr::group_by(ssc_miRNA) %>%
    dplyr::arrange(dplyr::desc(target_score), P.Value, .by_group = TRUE) %>%
    dplyr::slice(1:min(5, dplyr::n())) %>%
    dplyr::ungroup() %>%
    dplyr::transmute(
      from = ssc_miRNA,
      to = target_symbol,
      target_score = target_score,
      validated = validated,
      correlation = correlation,
      opposite_direction = opposite_direction,
      pathway_relevant = pathway_relevant
    ) %>%
    dplyr::distinct()
  
  # ----------------------------
  # 5. Build node table for plotting
  # miRNA nodes keep volcano metrics
  # gene nodes keep regulator count
  # ----------------------------
  miRNA_nodes <- vol %>%
    dplyr::filter(miRNA %in% network_miRNAs) %>%
    dplyr::transmute(
      name = miRNA,
      type = "miRNA",
      logFC = logFC,
      P.Value = P.Value,
      adj.P.Val = adj.P.Val
    )
  
  gene_nodes <- network_edges %>%
    dplyr::count(to, name = "regulator_count") %>%
    dplyr::transmute(
      name = to,
      type = "Gene",
      logFC = NA_real_,
      P.Value = NA_real_,
      adj.P.Val = NA_real_,
      regulator_count = regulator_count
    )
  
  plot_nodes <- dplyr::bind_rows(
    miRNA_nodes %>%
      dplyr::mutate(regulator_count = NA_real_),
    gene_nodes
  ) %>%
    dplyr::distinct(name, .keep_all = TRUE)
  
  # save network tables
  write.csv(network_edges, file.path(out_dir, "network_edges_for_plot.csv"), row.names = FALSE)
  write.csv(plot_nodes, file.path(out_dir, "network_nodes_for_plot.csv"), row.names = FALSE)
  
  # ----------------------------
  # 6. Build igraph object
  # ----------------------------
  g <- igraph::graph_from_data_frame(
    d = network_edges,
    vertices = plot_nodes,
    directed = TRUE
  )
  
  # node aesthetics
  V(g)$node_type <- plot_nodes$type[match(V(g)$name, plot_nodes$name)]
  V(g)$logFC <- plot_nodes$logFC[match(V(g)$name, plot_nodes$name)]
  V(g)$regulator_count <- plot_nodes$regulator_count[match(V(g)$name, plot_nodes$name)]
  
  V(g)$node_size <- ifelse(
    V(g)$node_type == "miRNA",
    8 + 2 * pmin(abs(V(g)$logFC), 4),
    4 + 2 * pmin(ifelse(is.na(V(g)$regulator_count), 1, V(g)$regulator_count), 6)
  )
  
  V(g)$node_color <- ifelse(
    V(g)$node_type == "miRNA" & !is.na(V(g)$logFC) & V(g)$logFC > 0, "red",
    ifelse(
      V(g)$node_type == "miRNA" & !is.na(V(g)$logFC) & V(g)$logFC < 0, "blue",
      "grey80"
    )
  )
  
  E(g)$edge_width <- 0.5 + 0.35 * pmin(network_edges$target_score, 10)
  E(g)$edge_color <- ifelse(network_edges$pathway_relevant == 1, "firebrick", "grey70")
  
  # ----------------------------
  # 7. Network plot
  # ----------------------------
  p_network <- ggraph(g, layout = "fr") +
    geom_edge_link(
      aes(width = target_score, color = factor(pathway_relevant)),
      alpha = 0.7,
      show.legend = TRUE
    ) +
    geom_node_point(
      aes(size = node_size, fill = I(node_color), shape = node_type),
      color = "black",
      stroke = 0.3
    ) +
    geom_node_text(
      aes(label = name),
      repel = TRUE,
      size = 3
    ) +
    scale_shape_manual(values = c("miRNA" = 23, "Gene" = 21)) +
    scale_edge_width(range = c(0.4, 2.2), name = "Target score") +
    scale_edge_color_manual(
      values = c("0" = "grey75", "1" = "firebrick"),
      labels = c("0" = "Not pathway-marked", "1" = "Pathway-relevant"),
      name = "Gene relevance"
    ) +
    theme_void() +
    labs(title = "Integrated miRNA-target network") +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold")
    )
  
  # ----------------------------
  # 8. Save separate network plot
  # ----------------------------
  ggsave(
    filename = file.path(out_dir, "network_top_miRNA_targets.png"),
    plot = p_network,
    width = 11,
    height = 8,
    dpi = 300
  )
  
  # ----------------------------
  # 9. Create a combined volcano + network figure
  # Uses patchwork if available; otherwise saves separately
  # ----------------------------
  if (!requireNamespace("patchwork", quietly = TRUE)) {
    install.packages("patchwork")
  }
  library(patchwork)
  
  combined_plot <- p_volcano + p_network + patchwork::plot_layout(widths = c(1.1, 1))
  
  ggsave(
    filename = file.path(out_dir, "volcano_network_integrated.png"),
    plot = combined_plot,
    width = 17,
    height = 8,
    dpi = 300
  )
  
  message("Volcano + network integration complete.")  
  
  
  # ============================================================
  # Heatmap of miRNA–target network
  # Uses network_edges from your integration pipeline
  # ============================================================
  
  # ============================================================
  # Safe miRNA–target heatmap
  # ============================================================
  
  library(pheatmap)
  library(dplyr)
  library(tidyr)
  library(matrixStats)
  
  # 1. Build matrix
  heat_df <- network_edges %>%
    dplyr::select(from, to, target_score) %>%
    tidyr::pivot_wider(
      names_from = to,
      values_from = target_score,
      values_fill = 0
    )
  
  heat_mat <- as.matrix(heat_df[, -1])
  rownames(heat_mat) <- heat_df$from
  
  # 2. Order miRNAs
  miRNA_order <- vol_net_miRNAs$miRNA
  heat_mat <- heat_mat[intersect(miRNA_order, rownames(heat_mat)), , drop = FALSE]
  
  # 3. Keep top variable target genes
  gene_var <- apply(heat_mat, 2, var, na.rm = TRUE)
  gene_var <- gene_var[!is.na(gene_var)]
  top_genes <- names(sort(gene_var, decreasing = TRUE))[1:min(30, length(gene_var))]
  heat_mat <- heat_mat[, top_genes, drop = FALSE]
  
  # 4. Remove rows/columns with zero variance
  row_var <- apply(heat_mat, 1, var, na.rm = TRUE)
  col_var <- apply(heat_mat, 2, var, na.rm = TRUE)
  
  heat_mat <- heat_mat[row_var > 0 & !is.na(row_var), , drop = FALSE]
  heat_mat <- heat_mat[, col_var > 0 & !is.na(col_var), drop = FALSE]
  
  # 5. Row-wise z-score safely
  safe_scale_row <- function(x) {
    s <- sd(x, na.rm = TRUE)
    m <- mean(x, na.rm = TRUE)
    if (is.na(s) || s == 0) {
      return(rep(0, length(x)))
    } else {
      return((x - m) / s)
    }
  }
  
  heat_mat_scaled <- t(apply(heat_mat, 1, safe_scale_row))
  
  # keep dimnames
  rownames(heat_mat_scaled) <- rownames(heat_mat)
  colnames(heat_mat_scaled) <- colnames(heat_mat)
  
  # 6. Remove any remaining bad values
  heat_mat_scaled[!is.finite(heat_mat_scaled)] <- 0
  
  # 7. Annotation
  miRNA_annot <- vol %>%
    dplyr::filter(miRNA %in% rownames(heat_mat_scaled)) %>%
    dplyr::select(miRNA, logFC)
  
  rownames(miRNA_annot) <- miRNA_annot$miRNA
  miRNA_annot <- miRNA_annot[, "logFC", drop = FALSE]
  
  # 8. Plot
  pheatmap(
    heat_mat_scaled,
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    show_rownames = TRUE,
    show_colnames = TRUE,
    fontsize_row = 8,
    fontsize_col = 7,
    border_color = NA,
    main = "miRNA–target interaction heatmap (Z-score)",
    annotation_row = miRNA_annot,
    color = colorRampPalette(c("blue", "white", "red"))(100)
  )
  
  # 9. Save
  png(
    filename = file.path(out_dir, "miRNA_target_heatmap.png"),
    width = 1400,
    height = 1000,
    res = 150
  )
  
  pheatmap(
    heat_mat_scaled,
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    show_rownames = TRUE,
    show_colnames = TRUE,
    fontsize_row = 8,
    fontsize_col = 7,
    border_color = NA,
    main = "miRNA–target interaction heatmap (Z-score)",
    annotation_row = miRNA_annot,
    color = colorRampPalette(c("blue", "white", "red"))(100)
  )
  
  dev.off()
    
# ============================================================
# KEGG enrichment for target genes of each miRNA
# - per-miRNA enrichment
# - save top 20 pathways for each miRNA
# - build heatmap across miRNAs and pathways
# ============================================================

# Required packages
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")

cran_pkgs <- c("dplyr", "tidyr", "readr", "stringr")
to_install_cran <- cran_pkgs[!cran_pkgs %in% rownames(installed.packages())]
if (length(to_install_cran) > 0) install.packages(to_install_cran)

bioc_pkgs <- c("clusterProfiler", "org.Hs.eg.db", "ComplexHeatmap", "circlize")
to_install_bioc <- bioc_pkgs[!bioc_pkgs %in% rownames(installed.packages())]
if (length(to_install_bioc) > 0) BiocManager::install(to_install_bioc, ask = FALSE, update = FALSE)

library(dplyr)
library(tidyr)
library(readr)
library(stringr)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ComplexHeatmap)
library(circlize)

# ------------------------------------------------------------
# 1. Check required input
# ------------------------------------------------------------
stopifnot(exists("target_tbl"))
stopifnot("ssc_miRNA" %in% colnames(target_tbl))
stopifnot("target_symbol" %in% colnames(target_tbl))
stopifnot(exists("out_dir"))

# ------------------------------------------------------------
# 2. Prepare miRNA -> target gene table
# ------------------------------------------------------------
mirna_targets <- target_tbl %>%
  dplyr::select(ssc_miRNA, target_symbol) %>%
  dplyr::filter(!is.na(ssc_miRNA), ssc_miRNA != "",
                !is.na(target_symbol), target_symbol != "") %>%
  dplyr::distinct()

write.csv(
  mirna_targets,
  file.path(out_dir, "miRNA_target_gene_table_for_per_miRNA_KEGG.csv"),
  row.names = FALSE
)

# ------------------------------------------------------------
# 3. Run KEGG enrichment separately for each miRNA
# ------------------------------------------------------------
mirna_list <- sort(unique(mirna_targets$ssc_miRNA))

kegg_results_list <- list()
top20_list <- list()

for (mir in mirna_list) {
  
  genes_this <- mirna_targets %>%
    dplyr::filter(ssc_miRNA == mir) %>%
    dplyr::pull(target_symbol) %>%
    unique()
  
  # Skip if too few genes
  if (length(genes_this) < 5) next
  
  # Convert gene symbols to ENTREZ IDs
  gene_map_this <- tryCatch({
    clusterProfiler::bitr(
      genes_this,
      fromType = "SYMBOL",
      toType = "ENTREZID",
      OrgDb = org.Hs.eg.db
    )
  }, error = function(e) NULL)
  
  if (is.null(gene_map_this) || nrow(gene_map_this) == 0) next
  
  entrez_this <- unique(gene_map_this$ENTREZID)
  if (length(entrez_this) < 5) next
  
  # KEGG enrichment
  kegg_this <- tryCatch({
    clusterProfiler::enrichKEGG(
      gene = entrez_this,
      organism = "hsa",
      pvalueCutoff = 0.05,
      pAdjustMethod = "BH"
    )
  }, error = function(e) NULL)
  
  if (is.null(kegg_this)) next
  
  kegg_df <- as.data.frame(kegg_this)
  if (nrow(kegg_df) == 0) next
  
  kegg_df$ssc_miRNA <- mir
  kegg_df$n_target_genes <- length(genes_this)
  
  kegg_results_list[[mir]] <- kegg_df
  
  top20_df <- kegg_df %>%
    dplyr::arrange(p.adjust, pvalue) %>%
    dplyr::slice(1:min(20, n()))
  
  top20_list[[mir]] <- top20_df
}

# Combine all results
all_kegg_results <- dplyr::bind_rows(kegg_results_list)
all_top20_kegg <- dplyr::bind_rows(top20_list)

# Save full and top20 tables
write.csv(
  all_kegg_results,
  file.path(out_dir, "per_miRNA_KEGG_all_results.csv"),
  row.names = FALSE
)

write.csv(
  all_top20_kegg,
  file.path(out_dir, "per_miRNA_KEGG_top20_pathways.csv"),
  row.names = FALSE
)



kidney_keywords <- c(
  "renal", "kidney", "nephro", "glomer", "tubule", "tubular", "podocyte",
  "fibrosis", "extracellular matrix", "ecm", "hypoxia", "oxidative stress",
  "apoptosis", "inflammation", "tgf", "pi3k", "mapk", "hif", "nf-kappa",
  "tnf", "ferroptosis", "focal adhesion", "cytokine"
)

kidney_pattern <- paste(kidney_keywords, collapse = "|")

all_top20_kegg_kidney <- all_top20_kegg[
  grepl(kidney_pattern, all_top20_kegg$Description, ignore.case = TRUE),
  ,
  drop = FALSE
]

write.csv(
  all_top20_kegg_kidney,
  file.path(out_dir, "per_miRNA_KEGG_top20_kidney_related.csv"),
  row.names = FALSE
)

# ------------------------------------------------------------
# 4. Build a heatmap matrix
# Rows = pathways
# Columns = miRNAs
# Values = -log10(p.adjust)
# ------------------------------------------------------------

if (nrow(all_top20_kegg) > 0) {
  
  heat_df <- all_top20_kegg %>%
    dplyr::mutate(
      pathway_label = Description,
      score = -log10(p.adjust)
    ) %>%
    dplyr::select(ssc_miRNA, pathway_label, score) %>%
    dplyr::distinct()
  
  # Keep pathways that appear in at least one miRNA top20
  heat_wide <- heat_df %>%
    tidyr::pivot_wider(
      names_from = ssc_miRNA,
      values_from = score,
      values_fill = 0
    )
  
  heat_mat <- heat_wide %>%
    tibble::column_to_rownames("pathway_label") %>%
    as.matrix()
  
  # Optional: keep only pathways appearing in >= 2 miRNAs
  pathway_freq <- rowSums(heat_mat > 0)
  heat_mat_filt <- heat_mat[pathway_freq >= 2, , drop = FALSE]
  
  # If filtering removes too much, fall back to all
  if (nrow(heat_mat_filt) < 2) {
    heat_mat_filt <- heat_mat
  }
  
  # Order rows/columns by clustering-ready summaries
  row_ord <- order(rowSums(heat_mat_filt), decreasing = TRUE)
  col_ord <- order(colSums(heat_mat_filt), decreasing = TRUE)
  
  heat_mat_filt <- heat_mat_filt[row_ord, col_ord, drop = FALSE]
  
  # Save matrix
  write.csv(
    heat_mat_filt,
    file.path(out_dir, "per_miRNA_KEGG_heatmap_matrix.csv"),
    row.names = TRUE
  )
  
  # ----------------------------------------------------------
  # 5. Plot heatmap
  # ----------------------------------------------------------
  col_fun <- circlize::colorRamp2(
    c(0, 1.3, 3, 5),
    c("white", "lightyellow", "orange", "red")
  )
  
  # 9. Save
  png(
    filename = file.path(out_dir, "miRNA_target_heatmap.png"),
    width = 1400,
    height = 1200,
    res = 150
  )
  
  draw(
    Heatmap(
      heat_mat_filt,
      name = "-log10(adj.P)",
      col = col_fun,
      cluster_rows = TRUE,
      cluster_columns = TRUE,
      show_row_names = TRUE,
      show_column_names = TRUE,
      row_names_gp = grid::gpar(fontsize = 8),
      column_names_gp = grid::gpar(fontsize = 8),
      column_names_rot = 45,
      row_title = "KEGG pathways",
      column_title = "miRNAs"
    ),
    heatmap_legend_side = "right"
  )
  
  dev.off()
}



# ------------------------------------------------------------
# 6. Optional: save top 20 pathways for each miRNA separately
# ------------------------------------------------------------
if (length(top20_list) > 0) {
  per_mirna_dir <- file.path(out_dir, "per_miRNA_top20_KEGG_tables")
  dir.create(per_mirna_dir, showWarnings = FALSE, recursive = TRUE)
  
  for (mir in names(top20_list)) {
    safe_name <- gsub("[^A-Za-z0-9_\\-]", "_", mir)
    write.csv(
      top20_list[[mir]],
      file.path(per_mirna_dir, paste0(safe_name, "_top20_KEGG.csv")),
      row.names = FALSE
    )
  }
}

message("Per-miRNA KEGG enrichment complete. Files saved to: ", out_dir)


# ============================================================
# 2D PCA plot for NanoString miRNA data
# Requires:
# - log2_mat
# - metadata with SampleID and Group
# - out_dir
# ============================================================

# install.packages(c("ggplot2", "ggrepel", "matrixStats"))
library(ggplot2)
library(ggrepel)
library(matrixStats)

# ----------------------------
# 1. Prepare PCA input
# ----------------------------
stopifnot(is.matrix(log2_mat))

# remove zero-variance miRNAs
var_genes <- matrixStats::rowVars(log2_mat)
log2_mat_filt <- log2_mat[var_genes > 0, , drop = FALSE]

# transpose so samples are rows
pca_input <- t(log2_mat_filt)

# ----------------------------
# 2. Run PCA
# ----------------------------
pca_res <- prcomp(
  pca_input,
  center = TRUE,
  scale. = TRUE
)

# ----------------------------
# 3. Build PCA dataframe
# ----------------------------
pca_df <- as.data.frame(pca_res$x[, 1:2, drop = FALSE])
pca_df$SampleID <- rownames(pca_df)
pca_df$Group <- metadata$Group[match(pca_df$SampleID, metadata$SampleID)]
pca_df$SampleNum <- as.numeric(sub(".*_(\\d+)$", "\\1", pca_df$SampleID))

# make sure order is consistent
pca_df$Group <- factor(pca_df$Group, levels = c("Control", "Sepsis"))

# variance explained
var_explained <- (pca_res$sdev^2) / sum(pca_res$sdev^2)
pc1_var <- round(var_explained[1] * 100, 1)
pc2_var <- round(var_explained[2] * 100, 1)

# ----------------------------
# 4. Plot PCA
# ----------------------------
p <- ggplot(pca_df, aes(x = PC1, y = PC2, color = Group)) +
  geom_point(size = 5, alpha = 0.9) +
  geom_text_repel(
    aes(label = SampleNum),
    size = 5,
    max.overlaps = Inf
  ) +
  stat_ellipse(
    aes(fill = Group),
    geom = "polygon",
    alpha = 0.15,
    color = NA
  ) +
  scale_color_manual(
    values = c("Control" = "blue", "Sepsis" = "red")
  ) +
  scale_fill_manual(
    values = c("Control" = "blue", "Sepsis" = "red")
  ) +
  labs(
    title = "PCA of Renal miRNA Expression",
    x = paste0("PC1 (", pc1_var, "%)"),
    y = paste0("PC2 (", pc2_var, "%)")
  ) +
  theme_classic(base_size = 16) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5, size = 20),
    axis.title = element_text(face = "bold", size = 18),
    axis.text = element_text(size = 18),
    legend.title = element_text(size = 18, face = "bold"),
    legend.text = element_text(size = 16)
  )

print(p)

# ----------------------------
# 5. Save as PNG and PDF
# ----------------------------
ggsave(
  filename = file.path(out_dir, "PCA_2D_miRNA.png"),
  plot = p,
  width = 8,
  height = 6,
  dpi = 300
)

ggsave(
  filename = file.path(out_dir, "PCA_2D_miRNA.pdf"),
  plot = p,
  width = 8,
  height = 6
)



# ============================================================
# 3D PCA plot for NanoString miRNA data
# Requires:
# - log2_mat
# - metadata with SampleID and Group
# - out_dir
# ============================================================

# Install packages if needed
cran_pkgs <- c("plotly", "htmlwidgets", "matrixStats")
to_install <- cran_pkgs[!cran_pkgs %in% rownames(installed.packages())]
if (length(to_install) > 0) install.packages(to_install)

library(plotly)
library(htmlwidgets)
library(matrixStats)

# ----------------------------
# 1. Prepare PCA input
# ----------------------------
# ============================================================
# 3D PCA with group ellipsoids
# Requires:
# - log2_mat
# - metadata with SampleID and Group
# - out_dir
# ============================================================

# Install packages if needed
cran_pkgs <- c("plotly", "htmlwidgets", "matrixStats", "MASS")
to_install <- cran_pkgs[!cran_pkgs %in% rownames(installed.packages())]
if (length(to_install) > 0) install.packages(to_install)

library(plotly)
library(htmlwidgets)
library(matrixStats)
library(MASS)

# ----------------------------
# 1. Prepare PCA input
# ----------------------------
stopifnot(is.matrix(log2_mat))

var_genes <- matrixStats::rowVars(log2_mat)
log2_mat_filt <- log2_mat[var_genes > 0, , drop = FALSE]

pca_input <- t(log2_mat_filt)

# ----------------------------
# 2. Run PCA
# ----------------------------
pca_res <- prcomp(
  pca_input,
  center = TRUE,
  scale. = TRUE
)

# ----------------------------
# 3. Build PCA dataframe
# ----------------------------
pca_df <- as.data.frame(pca_res$x[, 1:3, drop = FALSE])
pca_df$SampleID <- rownames(pca_df)
pca_df$Group <- metadata$Group[match(pca_df$SampleID, metadata$SampleID)]
pca_df$SampleNum <- as.numeric(sub(".*_(\\d+)$", "\\1", pca_df$SampleID))

var_explained <- (pca_res$sdev^2) / sum(pca_res$sdev^2)
pc1_var <- round(var_explained[1] * 100, 1)
pc2_var <- round(var_explained[2] * 100, 1)
pc3_var <- round(var_explained[3] * 100, 1)

# ----------------------------
# 4. Function to generate ellipsoid mesh
# ----------------------------
make_ellipsoid <- function(center, covmat, scale_factor = 2, n = 40) {
  u <- seq(0, 2 * pi, length.out = n)
  v <- seq(0, pi, length.out = n)
  
  x <- outer(cos(u), sin(v))
  y <- outer(sin(u), sin(v))
  z <- outer(rep(1, length(u)), cos(v))
  
  sphere <- array(0, dim = c(length(u), length(v), 3))
  sphere[, , 1] <- x
  sphere[, , 2] <- y
  sphere[, , 3] <- z
  
  eig <- eigen(covmat)
  transform_mat <- eig$vectors %*% diag(sqrt(pmax(eig$values, 0))) * scale_factor
  
  ellipsoid <- array(0, dim = dim(sphere))
  for (i in seq_len(dim(sphere)[1])) {
    for (j in seq_len(dim(sphere)[2])) {
      pt <- c(sphere[i, j, 1], sphere[i, j, 2], sphere[i, j, 3])
      new_pt <- center + transform_mat %*% pt
      ellipsoid[i, j, ] <- new_pt
    }
  }
  
  list(
    x = ellipsoid[, , 1],
    y = ellipsoid[, , 2],
    z = ellipsoid[, , 3]
  )
}

# ----------------------------
# 5. Base 3D PCA scatter
# ----------------------------
group_colors <- c("Control" = "blue", "Sepsis" = "red")

p3d <- plot_ly()

# add sample points
for (grp in levels(factor(pca_df$Group))) {
  df_grp <- pca_df[pca_df$Group == grp, , drop = FALSE]
  
  p3d <- p3d %>%
    add_trace(
      data = df_grp,
      x = ~PC1,
      y = ~PC2,
      z = ~PC3,
      type = "scatter3d",
      mode = "markers+text",
      text = ~SampleNum,
      textposition = "top center",
      textfont = list(size = 16),
      hovertext = ~paste(
        "SampleID:", SampleID,
        "<br>Group:", Group,
        "<br>PC1:", round(PC1, 2),
        "<br>PC2:", round(PC2, 2),
        "<br>PC3:", round(PC3, 2)
      ),
      hoverinfo = "text",
      marker = list(
        size = 10,   # increased dot size
        color = group_colors[[grp]]
      ),
      name = grp
    )
}

# ----------------------------
# 6. Add ellipsoids
# ----------------------------
for (grp in levels(factor(pca_df$Group))) {
  df_grp <- pca_df[pca_df$Group == grp, c("PC1", "PC2", "PC3"), drop = FALSE]
  
  if (nrow(df_grp) >= 4) {
    center <- colMeans(df_grp)
    covmat <- stats::cov(df_grp)
    
    ell <- make_ellipsoid(center = center, covmat = covmat, scale_factor = 2, n = 35)
    
    p3d <- p3d %>%
      add_surface(
        x = ell$x,
        y = ell$y,
        z = ell$z,
        opacity = 0.18,
        showscale = FALSE,
        surfacecolor = matrix(1, nrow = nrow(ell$x), ncol = ncol(ell$x)),
        colorscale = list(c(0, group_colors[[grp]]), c(1, group_colors[[grp]])),
        hoverinfo = "skip",
        name = paste0(grp, " ellipsoid"),
        inherit = FALSE,
        showlegend = FALSE
      )
  }
}

# ----------------------------
# 7. Layout
# ----------------------------
p3d <- p3d %>%
  layout(
    title = list(
      text = "3D PCA of Renal miRNA Expression",
      font = list(size = 22)
    ),
    legend = list(
      title = list(text = "Group"),
      font = list(size = 24)
    ),
    scene = list(
      xaxis = list(
        title = list(
          text = paste0("PC1 (", pc1_var, "%)"),
          font = list(size = 18)
        ),
        tickfont = list(size = 18)
      ),
      yaxis = list(
        title = list(
          text = paste0("PC2 (", pc2_var, "%)"),
          font = list(size = 18)
        ),
        tickfont = list(size = 18)
      ),
      zaxis = list(
        title = list(
          text = paste0("PC3 (", pc3_var, "%)"),
          font = list(size = 18)
        ),
        tickfont = list(size = 18)
      )
    )
  )

# ----------------------------
# 8. Save outputs
# ----------------------------
htmlwidgets::saveWidget(
  p3d,
  file = file.path(out_dir, "PCA_3D_miRNA_with_ellipsoids.html"),
  selfcontained = TRUE
)

try({
  plotly::save_image(
    p3d,
    file = file.path(out_dir, "PCA_3D_miRNA_with_ellipsoids.png"),
    width = 1200,
    height = 900
  )
}, silent = TRUE)

p3d



install.packages("webshot2")
library(webshot2)

htmlwidgets::saveWidget(
  p3d,
  file = file.path(out_dir, "PCA_3D_miRNA_with_ellipsoids.html"),
  selfcontained = TRUE
)

webshot2::webshot(
  url = file.path(out_dir, "PCA_3D_miRNA_with_ellipsoids.html"),
  file = file.path(out_dir, "PCA_3D_miRNA_with_ellipsoids_webshot.png"),
  vwidth = 1400,
  vheight = 1000
)



# ============================================================
# Enrichment dot plot
# - y-axis: term id + term name
# - x-axis: -log10(p-value)
# - dot color: chart color
# - dot size: intersection size
# ============================================================

library(readr)
library(dplyr)
library(ggplot2)

# ----------------------------------------------------------------
# 1. Read your enrichment file
#    Replace with your actual file path if needed
# ----------------------------------------------------------------
enrich_df <- read_csv(file.path(out_dir, "miR_374b-5p_enrichment.csv"))

enrich_df<- as.data.frame(enrich_df, stringsAsFactors = FALSE)

# ----------------------------------------------------------------
# 2. Clean and prepare data
# ----------------------------------------------------------------
colnames(enrich_df) <- c(
  "background_size",
  "chart_color",
  "description",
  "evidence_codes",
  "intersecting_genes",
  "intersection_size",
  "network_suid",
  "nodes_suid",
  "p_value",
  "precision",
  "query_size",
  "recall",
  "source",
  "term_id",
  "term_name",
  "term_size"
)

# ----------------------------------------------------------------
# 3. Order terms so the most significant are at the top
# ----------------------------------------------------------------

plot_df <- enrich_df %>%
  mutate(
    p_value = as.numeric(p_value),
    intersection_size = as.numeric(intersection_size),
    term_label = paste(term_id, term_name, sep = ": "),
    neg_log10_p = -log10(p_value)
  ) %>%
  filter(!is.na(p_value), p_value > 0, !is.na(intersection_size)) %>%
  arrange(p_value)

# take top 50 using base indexing
plot_df <- plot_df[1:min(20, nrow(plot_df)), , drop = FALSE]

plot_df$plot_color <- "#FFCDCD"
plot_df$plot_color[1:min(8, nrow(plot_df))] <- plot_df$chart_color[1:min(8, nrow(plot_df))]


plot_df$term_label <- factor(
  plot_df$term_label,
  levels = plot_df$term_label[order(plot_df$neg_log10_p)]
)

# ----------------------------------------------------------------
# 4. Plot
# ----------------------------------------------------------------
p <- ggplot(plot_df, aes(x = neg_log10_p, y = term_label)) +
  geom_point(
    aes(size = intersection_size, fill = chart_color),
    shape = 21,
    color = "black",
    stroke = 0.4,
    alpha = 0.9
  ) +
  scale_fill_identity() +
  scale_size_continuous(name = "Intersection size") +
  scale_x_continuous(
    breaks = c(1.8, 2.2, 2.6, 3.0, 3.4),
    limits = c(1.8, 3.4)
  ) +
  labs(
    title = "GO Molecular Function Enrichment Dot Plot",
    x = expression(-log[10](p-value)),
    y = "Term ID: Term Name"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.y = element_text(size = 9),
    plot.title = element_text(face = "bold", hjust = 0.5)
  )

print(p)

ggsave(
  filename = file.path(out_dir, "GO_MF_enrichment_dotplot_374-5p.png"),
  plot = p,
  width = 15,
  height = 7,
  dpi = 300
)



# ----------------------------------------------------------------
# 1. Read your enrichment file
#    Replace with your actual file path if needed
# ----------------------------------------------------------------
enrich_df2 <- read_csv(file.path(out_dir, "miR-196b_enrichment.csv"))

# ----------------------------------------------------------------
# 2. Clean and prepare data
# ----------------------------------------------------------------
colnames(enrich_df2) <- c(
  "background_size",
  "chart_color",
  "description",
  "evidence_codes",
  "intersecting_genes",
  "intersection_size",
  "network_suid",
  "nodes_suid",
  "p_value",
  "precision",
  "query_size",
  "recall",
  "source",
  "term_id",
  "term_name",
  "term_size"
)

# ----------------------------------------------------------------
# 3. Order terms so the most significant are at the top
# ----------------------------------------------------------------
enrich_df2 <- as.data.frame(enrich_df, stringsAsFactors = FALSE)

plot_df2 <- enrich_df2 %>%
  mutate(
    p_value = as.numeric(p_value),
    intersection_size = as.numeric(intersection_size),
    term_label = paste(term_id, term_name, sep = ": "),
    neg_log10_p = -log10(p_value)
  ) %>%
  filter(!is.na(p_value), p_value > 0, !is.na(intersection_size)) %>%
  arrange(p_value)

# take top 50 using base indexing
plot_df2 <- plot_df[1:min(20, nrow(plot_df)), , drop = FALSE]

plot_df2$plot_color <- "#FFCDCD"
plot_df2$plot_color[1:min(8, nrow(plot_df))] <- plot_df$chart_color[1:min(8, nrow(plot_df))]


plot_df2$term_label <- factor(
  plot_df2$term_label,
  levels = plot_df2$term_label[order(plot_df2$neg_log10_p)]
)

# ----------------------------------------------------------------
# 4. Plot
# ----------------------------------------------------------------
p2 <- ggplot(plot_df2, aes(x = neg_log10_p, y = term_label)) +
  geom_point(
    aes(size = intersection_size, fill = chart_color),
    shape = 21,
    color = "black",
    stroke = 0.4,
    alpha = 0.9
  ) +
  scale_fill_identity() +
  scale_size_continuous(name = "Intersection size") +
  scale_x_continuous(
    breaks = c(2.2, 2.6, 3.0, 3.4),
    limits = c(2.2, 3.4)
  ) +
  labs(
    title = "GO Molecular Function Enrichment Dot Plot",
    x = expression(-log[10](p-value)),
    y = "Term ID: Term Name"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.y = element_text(size = 9),
    plot.title = element_text(face = "bold", hjust = 0.5)
  )

print(p2)

ggsave(
  filename = file.path(out_dir, "GO_MF_enrichment_dotplot_miR196b.png"),
  plot = p,
  width = 15,
  height = 7,
  dpi = 300
)