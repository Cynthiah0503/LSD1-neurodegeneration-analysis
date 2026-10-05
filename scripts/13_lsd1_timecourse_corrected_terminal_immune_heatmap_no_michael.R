# 13_lsd1_timecourse_corrected_terminal_immune_heatmap_no_michael.R
#
# Purpose: Corrected LSD1 time-course immune heatmap excluding Michael terminal samples.
# Inputs are expected under data/processed/ or data/external/ relative to this repository.
# Outputs are written under results/ or script-defined subfolders.

# Resolve paths relative to the repository root when the script is run with Rscript.
get_script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) > 0) {
    return(dirname(normalizePath(sub("^--file=", "", file_arg[[1]]), mustWork = FALSE)))
  }
  getwd()
}
repo_root <- normalizePath(file.path(get_script_dir(), ".."), mustWork = FALSE)
suppressPackageStartupMessages({
  library(ComplexHeatmap)
  library(dplyr)
  library(grid)
  library(readxl)
  library(circlize)
  library(tibble)
  library(sva)
  library(ggplot2)
})

katz_dir <- repo_root

lsd1_dir  <- file.path(katz_dir, "timepoint_lsd1")
trem2_dir <- file.path(katz_dir, "Trem2 terminal experiment")

wm_out_dir <- file.path(
  lsd1_dir,
  "combine terminal_CORRECTED_LSD1terminal_plus_Trem2WT_LSD1WT_and_LSD1KO_no_Michael"
)
if (!dir.exists(wm_out_dir)) dir.create(wm_out_dir, recursive = TRUE)

# LSD1 timepoint files
wk3_file    <- file.path(lsd1_dir, "3 weeks feature counts.xlsx")
wk4_file    <- file.path(lsd1_dir, "4 weeks feature counts.xlsx")
early_file  <- file.path(lsd1_dir, "early_onset_featurecounts.xlsx")
term_file   <- file.path(lsd1_dir, "Terminal feature counts.xlsx")
allc_file   <- file.path(lsd1_dir, "NIHMS472534-supplement-02.csv")
module_file <- file.path(lsd1_dir, "modules.csv")

# Trem2 raw heatmap featureCounts file.
# This is the same import logic as your working Trem2 heatmap code:
#   read Trem2heatmap feature counts.xlsx
#   use sample columns starting with LSD1
#   parse groups from LSD1WT/LSD1KO + Trem2WW/Trem2MM
#   select Trem2WW background only: LSD1WT_Trem2WW and LSD1KO_Trem2WW, displayed as Trem2WT_LSD1WT and Trem2WT_LSD1KO.
trem2_heatmap_file <- file.path(
  trem2_dir,
  "Trem2heatmap feature counts.xlsx"
)

# Hard safety block: never use Michael in this corrected workflow.
michael_file_that_must_not_be_used <- file.path(lsd1_dir, "featurecountsLSD1del vs. LSD1wt .xlsx")

required_files <- c(
  wk3_file,
  wk4_file,
  early_file,
  term_file,
  allc_file,
  module_file,
  trem2_heatmap_file
)

missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files) > 0) {
  stop("Missing required input files:\n", paste(missing_files, collapse = "\n"))
}

if (file.exists(michael_file_that_must_not_be_used)) {
  message("Michael file exists on disk but will NOT be read: ", michael_file_that_must_not_be_used)
}

cat("Output folder:\n", wm_out_dir, "\n\n")

wm_wk3        <- readxl::read_excel(wk3_file, sheet = 1)
wm_wk4        <- readxl::read_excel(wk4_file, sheet = 1)
wm_early      <- readxl::read_excel(early_file, sheet = 1)
wm_term       <- readxl::read_excel(term_file, sheet = 1)
wm_trem2_term <- readxl::read_excel(trem2_heatmap_file, sheet = 1)

wm_allc <- read.csv(allc_file, header = TRUE, check.names = FALSE)
wm_mod  <- read.csv(module_file, header = TRUE, check.names = FALSE)

# 2) Fix names / BOM

wm_fix_names <- function(df) {
  names(df) <- sub("^\\ufeff", "", names(df))
  bad <- is.na(names(df)) | names(df) == ""
  if (any(bad)) names(df)[bad] <- paste0("V", seq_len(sum(bad)))
  names(df) <- make.unique(names(df))
  df
}

wm_wk3        <- wm_fix_names(wm_wk3)
wm_wk4        <- wm_fix_names(wm_wk4)
wm_early      <- wm_fix_names(wm_early)
wm_term       <- wm_fix_names(wm_term)
wm_trem2_term <- wm_fix_names(wm_trem2_term)
wm_allc       <- wm_fix_names(wm_allc)
wm_mod        <- wm_fix_names(wm_mod)

wm_allc$Gene_Symbol <- toupper(trimws(as.character(wm_allc$Gene_Symbol)))

# 2.5) Display-name helpers

wm_make_display_labels <- function(samples) {
  samples <- as.character(samples)

  genotype <- dplyr::case_when(
    grepl("Trem2WT_LSD1KO", samples, ignore.case = TRUE) ~ "LSD1KO",
    grepl("Trem2WT_LSD1WT", samples, ignore.case = TRUE) ~ "LSD1WT",
    grepl("LSD1\\s*KO|LSD1[_ ]?KO|LSD1[_ ]?Del|LSD1\\s*Del|Cre\\+", samples, ignore.case = TRUE) ~ "LSD1KO",
    grepl("LSD1\\s*WT|LSD1[_ ]?WT|Cre\\-", samples, ignore.case = TRUE) ~ "LSD1WT",
    TRUE ~ "UNK"
  )

  dataset <- dplyr::case_when(
    grepl("^3w__", samples)           ~ "3W",
    grepl("^4w__", samples)           ~ "4W",
    grepl("^early__", samples)        ~ "Early",
    grepl("^term_lsd1__", samples)    ~ "Term",
    grepl("^term_trem2wt__", samples) ~ "Trem2WT",
    TRUE ~ "Sample"
  )

  labels <- character(length(samples))

  non_terminal <- !(grepl("^term_lsd1__", samples) | grepl("^term_trem2wt__", samples))
  if (any(non_terminal)) {
    idx_non <- ave(seq_along(samples)[non_terminal],
                   paste(dataset[non_terminal], genotype[non_terminal], sep = "__"),
                   FUN = seq_along)
    labels[non_terminal] <- paste0(dataset[non_terminal], "_", genotype[non_terminal], "_", sprintf("%02d", idx_non))
  }

  terminal <- !non_terminal
  if (any(terminal)) {
    for (g in c("LSD1WT", "LSD1KO", "UNK")) {
      idx <- which(terminal & genotype == g)
      if (length(idx) > 0) {
        seq_idx <- seq_along(idx)
        labels[idx] <- ifelse(
          grepl("^term_trem2wt__", samples[idx]) & g == "LSD1WT",
          paste0("Trem2WT_LSD1WT_", g, "_", sprintf("%02d", seq_idx)),
          ifelse(
            grepl("^term_trem2wt__", samples[idx]) & g == "LSD1KO",
            paste0("Trem2WT_LSD1KO_", g, "_", sprintf("%02d", seq_idx)),
            paste0("Term_", g, "_", sprintf("%02d", seq_idx))
          )
        )
      }
    }
  }

  labels
}
wm_dataset_display <- function(samples) {
  dplyr::case_when(
    grepl("^3w__", samples)                 ~ "3 weeks",
    grepl("^4w__", samples)                 ~ "4 weeks",
    grepl("^early__", samples)              ~ "Early onset",
    grepl("^term_lsd1__", samples)           ~ "Terminal",
    grepl("^term_trem2wt__", samples) ~ "Terminal",
    TRUE ~ NA_character_
  )
}

wm_terminal_source <- function(samples) {
  dplyr::case_when(
    grepl("^term_lsd1__", samples)            ~ "LSD1_terminal",
    grepl("^term_trem2wt__", samples) ~ "Trem2WT_terminal_addon",
    TRUE ~ NA_character_
  )
}

# 3) Immune functions -> modules

wm_target_term <- "Immune functions"

wm_target_modules <- wm_mod %>%
  dplyr::filter(trimws(CategoryTerm) == wm_target_term) %>%
  dplyr::pull(Module) %>%
  unique()

cat("Immune modules:", paste(wm_target_modules, collapse = ", "), "\n")

# 4) Gene -> module map

wm_gene2module <- wm_allc %>%
  dplyr::select(Gene_Symbol, Module) %>%
  dplyr::filter(!is.na(Gene_Symbol), !is.na(Module)) %>%
  dplyr::distinct()

wm_immune_genes <- wm_gene2module %>%
  dplyr::filter(Module %in% wm_target_modules) %>%
  dplyr::pull(Gene_Symbol) %>%
  unique()

# 5) Generic helpers for featureCounts-like tables

wm_pick_gene_col <- function(df, context = "input file") {
  priority <- c(
    "Gene_name", "Gene_Name", "GeneSymbol", "Gene_Symbol", "gene_symbol",
    "Symbol", "SYMBOL", "GeneID", "Gene_ID", "Geneid", "gene_id",
    "Gene", "gene", "external_gene_name"
  )
  hit <- priority[priority %in% names(df)]
  if (length(hit) > 0) return(hit[1])

  # In your older LSD1 files, the last column is the gene symbol.
  tail(names(df), 1)
}

wm_numeric_sample_cols <- function(df, exclude_cols = character()) {
  anno_patterns <- paste(
    c(
      "^Geneid$", "^Gene_ID$", "^Gene_name$", "^Gene_Name$", "^Gene_Symbol$",
      "^Symbol$", "^Chr$", "^Start$", "^End$", "^Strand$", "^Length$",
      "^gene$", "^gene_id$", "^external_gene_name$"
    ),
    collapse = "|"
  )

  candidates <- setdiff(names(df), exclude_cols)
  candidates <- candidates[!grepl(anno_patterns, candidates, ignore.case = TRUE)]

  is_num_like <- vapply(candidates, function(x) {
    vals <- suppressWarnings(as.numeric(df[[x]]))
    mean(!is.na(vals)) > 0.80
  }, logical(1))

  candidates[is_num_like]
}

wm_aggregate_counts <- function(df, gene_col, sample_cols, context = "input file") {
  if (length(sample_cols) == 0) {
    stop("No sample columns selected for ", context)
  }

  df$Gene_name <- toupper(trimws(as.character(df[[gene_col]])))

  out <- df[, c("Gene_name", sample_cols), drop = FALSE] %>%
    dplyr::filter(!is.na(Gene_name), Gene_name != "") %>%
    dplyr::group_by(Gene_name) %>%
    dplyr::summarise(
      dplyr::across(dplyr::all_of(sample_cols), ~ sum(as.numeric(.x), na.rm = TRUE)),
      .groups = "drop"
    )

  list(df_agg = out, sample_cols = sample_cols)
}

# 6) Prep standard LSD1 feature-count files

wm_prep_featurecounts_standard <- function(df, context = "standard LSD1 file") {
  gene_col <- wm_pick_gene_col(df, context)

  # Old LSD1 timepoint files use sample names starting with LSD1.
  sample_cols <- grep("^LSD1", names(df), value = TRUE)

  # Fallback in case a file was saved with numeric sample IDs instead.
  if (length(sample_cols) == 0) {
    sample_cols <- wm_numeric_sample_cols(df, exclude_cols = gene_col)
  }

  wm_aggregate_counts(df, gene_col, sample_cols, context)
}

# 7) Prep early onset merged file

wm_prep_featurecounts_early <- function(df) {
  if (!("Gene_ID" %in% names(df))) {
    stop("For early_onset_featurecounts.xlsx, column 'Gene_ID' was not found.")
  }

  df$Gene_name <- toupper(trimws(as.character(df[["Gene_ID"]])))

  sample_cols <- c(
    "0501_LSD1WT",
    "0502_LSD1WT",
    "0503_LSD1WT",
    "0504_LSD1KO",
    "0505_LSD1KO",
    "0506_LSD1KO"
  )

  missing_cols <- setdiff(sample_cols, names(df))
  if (length(missing_cols) > 0) {
    stop("Missing sample columns in early onset file: ", paste(missing_cols, collapse = ", "))
  }

  wm_aggregate_counts(df, "Gene_ID", sample_cols, "early_onset_featurecounts.xlsx")
}

# 8) Prep Trem2WT_LSD1KO terminal group

wm_make_group_id_4_from_trem2_samples <- function(samples) {
  s <- as.character(samples)

  is_lsd1wt <- grepl("LSD1\\s*WT|LSD1[_ ]?WT", s, ignore.case = TRUE)
  is_lsd1ko <- grepl("LSD1\\s*KO|LSD1[_ ]?KO", s, ignore.case = TRUE)

  lsd1_part <- dplyr::case_when(
    is_lsd1wt ~ "LSD1WT",
    is_lsd1ko ~ "LSD1KO",
    TRUE ~ NA_character_
  )

  is_trem2ww <- grepl("Trem2WW", s, ignore.case = TRUE)
  is_trem2mm <- grepl("Trem2MM", s, ignore.case = TRUE)

  trem2_part <- dplyr::case_when(
    is_trem2ww ~ "Trem2WW",
    is_trem2mm ~ "Trem2MM",
    TRUE ~ NA_character_
  )

  ifelse(
    !is.na(lsd1_part) & !is.na(trem2_part),
    paste(lsd1_part, trem2_part, sep = "_"),
    NA_character_
  )
}

wm_prep_trem2wt_terminal_addon <- function(df) {
  # EXACT Trem2 import logic from your working heatmap:
  # gene symbol = last column; sample columns = names starting with LSD1.
  gene_col <- tail(names(df), 1)
  sample_cols_all <- grep("^LSD1", names(df), value = TRUE)

  if (length(sample_cols_all) == 0) {
    stop("No Trem2 sample columns found starting with 'LSD1' in Trem2heatmap feature counts.xlsx.")
  }

  gid4 <- wm_make_group_id_4_from_trem2_samples(sample_cols_all)

  mapping_df <- data.frame(
    sample_col = sample_cols_all,
    group_internal = gid4,
    group_display = gsub("Trem2WW", "Trem2WT", gsub("Trem2MM", "Trem2KO", gid4, fixed = TRUE), fixed = TRUE),
    stringsAsFactors = FALSE
  )
  write.csv(
    mapping_df,
    file.path(wm_out_dir, "CHECK_Trem2heatmap_all_sample_group_mapping.csv"),
    row.names = FALSE
  )

  # Corrected target:
  #   keep the Trem2WT background, BOTH LSD1WT and LSD1KO.
  #   internal Trem2WW = display Trem2WT.
  target_groups <- c("LSD1WT_Trem2WW", "LSD1KO_Trem2WW")
  keep <- gid4 %in% target_groups
  sample_cols <- sample_cols_all[keep]
  selected_gid <- gid4[keep]

  if (length(sample_cols) == 0) {
    stop(
      "No samples mapped to LSD1WT_Trem2WW or LSD1KO_Trem2WW in Trem2heatmap feature counts.xlsx.\n",
      "Check: ", file.path(wm_out_dir, "CHECK_Trem2heatmap_all_sample_group_mapping.csv")
    )
  }

  write.csv(
    data.frame(
      selected_Trem2WT_background_columns = sample_cols,
      selected_internal_group = selected_gid,
      selected_display_group = gsub("Trem2WW", "Trem2WT", selected_gid, fixed = TRUE),
      source_file = trem2_heatmap_file,
      stringsAsFactors = FALSE
    ),
    file.path(wm_out_dir, "CHECK_selected_Trem2WT_background_terminal_columns.csv"),
    row.names = FALSE
  )

  out <- wm_aggregate_counts(df, gene_col, sample_cols, "Trem2heatmap feature counts.xlsx: Trem2WT background")

  clean_names <- gsub("[^A-Za-z0-9]+", "_", out$sample_cols)
  clean_names <- gsub("_+", "_", clean_names)
  clean_names <- gsub("^_|_$", "", clean_names)

  new_cols <- ifelse(
    selected_gid == "LSD1WT_Trem2WW",
    paste0("Trem2WT_LSD1WT_", clean_names),
    paste0("Trem2WT_LSD1KO_", clean_names)
  )

  names(out$df_agg)[match(out$sample_cols, names(out$df_agg))] <- new_cols
  out$sample_cols <- new_cols

  out
}

# 9) Sample annotation

wm_group_colors_cre <- c(
  "LSD1KO" = "firebrick",
  "LSD1WT" = "steelblue"
)

wm_make_sample_anno_cre <- function(samples) {
  is_ko <- grepl("Trem2WT_LSD1KO", samples, ignore.case = TRUE) |
    grepl("LSD1\\s*KO|LSD1[_ ]?KO|LSD1[_ ]?Del|LSD1\\s*Del|Cre\\+", samples, ignore.case = TRUE)
  is_wt <- grepl("Trem2WT_LSD1WT", samples, ignore.case = TRUE) |
    grepl("LSD1\\s*WT|LSD1[_ ]?WT|Cre\\-", samples, ignore.case = TRUE)

  grp <- dplyr::case_when(
    is_ko ~ "LSD1KO",
    is_wt ~ "LSD1WT",
    TRUE  ~ NA_character_
  )

  HeatmapAnnotation(
    Group  = grp,
    col    = list(Group = wm_group_colors_cre),
    na_col = "grey90",
    show_annotation_name = TRUE
  )
}
wm_make_groupmean_anno_cre <- function(group_cols) {
  grp <- ifelse(grepl("LSD1KO", group_cols, ignore.case = TRUE), "LSD1KO", "LSD1WT")
  HeatmapAnnotation(
    Group = grp,
    col   = list(Group = wm_group_colors_cre),
    na_col = "grey90",
    show_annotation_name = TRUE
  )
}

# 10) df_agg -> matrix

wm_dfagg_to_mat <- function(df_agg, sample_cols) {
  mat <- as.matrix(df_agg[, sample_cols, drop = FALSE])
  rownames(mat) <- df_agg$Gene_name
  storage.mode(mat) <- "numeric"
  mat
}

# 11) Build full union matrix

wm_make_full_union_mat <- function(mat, all_genes) {
  out <- matrix(
    0,
    nrow = length(all_genes),
    ncol = ncol(mat),
    dimnames = list(all_genes, colnames(mat))
  )
  out[rownames(mat), ] <- mat
  out
}

# 12) Row z-score

wm_row_zscore <- function(mat) {
  if (nrow(mat) == 0) return(mat)
  z <- t(scale(t(mat)))
  z[!is.finite(z)] <- 0
  z
}

# 13) Collapse group mean

wm_collapse_group_mean_from_z <- function(mat_z) {
  if (nrow(mat_z) == 0) return(mat_z)

  groups <- ifelse(
    grepl("Trem2WT_LSD1KO", colnames(mat_z), ignore.case = TRUE) |
      grepl("LSD1\\s*KO|LSD1[_ ]?KO|LSD1[_ ]?Del|LSD1\\s*Del|Cre\\+", colnames(mat_z), ignore.case = TRUE),
    "LSD1KO",
    ifelse(
      grepl("Trem2WT_LSD1WT", colnames(mat_z), ignore.case = TRUE) |
        grepl("LSD1\\s*WT|LSD1[_ ]?WT|Cre\\-", colnames(mat_z), ignore.case = TRUE),
      "LSD1WT",
      NA
    )
  )

  if (any(is.na(groups))) {
    bad <- colnames(mat_z)[is.na(groups)]
    stop("Could not assign some columns to LSD1WT/LSD1KO: ", paste(bad, collapse = ", "))
  }

  out <- sapply(c("LSD1WT", "LSD1KO"), function(g) {
    idx <- which(groups == g)
    if (length(idx) == 0) {
      rep(NA_real_, nrow(mat_z))
    } else if (length(idx) == 1) {
      mat_z[, idx]
    } else {
      rowMeans(mat_z[, idx, drop = FALSE], na.rm = TRUE)
    }
  })

  rownames(out) <- rownames(mat_z)
  out
}

# 14) Prepare all datasets

wm_wk3_p        <- wm_prep_featurecounts_standard(wm_wk3, "3 weeks feature counts.xlsx")
wm_wk4_p        <- wm_prep_featurecounts_standard(wm_wk4, "4 weeks feature counts.xlsx")
wm_early_p      <- wm_prep_featurecounts_early(wm_early)
wm_term_p       <- wm_prep_featurecounts_standard(wm_term, "Terminal feature counts.xlsx")
wm_trem2_term_p <- wm_prep_trem2wt_terminal_addon(wm_trem2_term)

wm_wk3_mat        <- wm_dfagg_to_mat(wm_wk3_p$df_agg,        wm_wk3_p$sample_cols)
wm_wk4_mat        <- wm_dfagg_to_mat(wm_wk4_p$df_agg,        wm_wk4_p$sample_cols)
wm_early_mat      <- wm_dfagg_to_mat(wm_early_p$df_agg,      wm_early_p$sample_cols)
wm_term_mat       <- wm_dfagg_to_mat(wm_term_p$df_agg,       wm_term_p$sample_cols)
wm_trem2_term_mat <- wm_dfagg_to_mat(wm_trem2_term_p$df_agg, wm_trem2_term_p$sample_cols)

colnames(wm_wk3_mat)        <- paste0("3w__", colnames(wm_wk3_mat))
colnames(wm_wk4_mat)        <- paste0("4w__", colnames(wm_wk4_mat))
colnames(wm_early_mat)      <- paste0("early__", colnames(wm_early_mat))
colnames(wm_term_mat)       <- paste0("term_lsd1__", colnames(wm_term_mat))
colnames(wm_trem2_term_mat) <- paste0("term_trem2wt__", colnames(wm_trem2_term_mat))

# Final safety check: no Michael samples.
all_raw_sample_names <- c(
  colnames(wm_wk3_mat),
  colnames(wm_wk4_mat),
  colnames(wm_early_mat),
  colnames(wm_term_mat),
  colnames(wm_trem2_term_mat)
)

if (any(grepl("mich|michael", all_raw_sample_names, ignore.case = TRUE))) {
  stop("Michael/Mich sample detected. Stopping because this corrected workflow must exclude Michael.")
}

# 15) Build GLOBAL matrix

wm_all_genes_union <- unique(c(
  rownames(wm_wk3_mat),
  rownames(wm_wk4_mat),
  rownames(wm_early_mat),
  rownames(wm_term_mat),
  rownames(wm_trem2_term_mat)
))

wm_wk3_full        <- wm_make_full_union_mat(wm_wk3_mat,        wm_all_genes_union)
wm_wk4_full        <- wm_make_full_union_mat(wm_wk4_mat,        wm_all_genes_union)
wm_early_full      <- wm_make_full_union_mat(wm_early_mat,      wm_all_genes_union)
wm_term_full       <- wm_make_full_union_mat(wm_term_mat,       wm_all_genes_union)
wm_trem2_term_full <- wm_make_full_union_mat(wm_trem2_term_mat, wm_all_genes_union)

wm_combined_counts <- cbind(
  wm_wk3_full,
  wm_wk4_full,
  wm_early_full,
  wm_term_full,
  wm_trem2_term_full
)

if (any(grepl("mich|michael", colnames(wm_combined_counts), ignore.case = TRUE))) {
  stop("Michael/Mich sample detected after combining. Stopping.")
}

# 15.5) Sample information for QC and plotting

wm_sample_info <- data.frame(
  Sample = colnames(wm_combined_counts),
  stringsAsFactors = FALSE
)

wm_sample_info$Dataset <- dplyr::case_when(
  grepl("^3w__", wm_sample_info$Sample)                  ~ "3weeks",
  grepl("^4w__", wm_sample_info$Sample)                  ~ "4weeks",
  grepl("^early__", wm_sample_info$Sample)               ~ "early_onset",
  grepl("^term_lsd1__", wm_sample_info$Sample)            ~ "terminal",
  grepl("^term_trem2wt__", wm_sample_info$Sample) ~ "terminal",
  TRUE ~ NA_character_
)

wm_sample_info$DatasetDisplay <- wm_dataset_display(wm_sample_info$Sample)
wm_sample_info$TerminalSource <- wm_terminal_source(wm_sample_info$Sample)
wm_sample_info$DisplayName    <- wm_make_display_labels(wm_sample_info$Sample)

wm_sample_info$Genotype <- dplyr::case_when(
  grepl("Trem2WT_LSD1KO", wm_sample_info$Sample, ignore.case = TRUE) ~ "LSD1KO",
  grepl("Trem2WT_LSD1WT", wm_sample_info$Sample, ignore.case = TRUE) ~ "LSD1WT",
  grepl("LSD1\\s*KO|LSD1[_ ]?KO|LSD1[_ ]?Del|LSD1\\s*Del|Cre\\+", wm_sample_info$Sample, ignore.case = TRUE) ~ "LSD1KO",
  grepl("LSD1\\s*WT|LSD1[_ ]?WT|Cre\\-", wm_sample_info$Sample, ignore.case = TRUE) ~ "LSD1WT",
  TRUE ~ NA_character_
)

if (any(is.na(wm_sample_info$Dataset))) {
  stop("Some samples could not be assigned to Dataset.")
}
if (any(is.na(wm_sample_info$Genotype))) {
  write.csv(wm_sample_info, file.path(wm_out_dir, "DEBUG_sample_info_unassigned_genotype.csv"), row.names = FALSE)
  stop("Some samples could not be assigned to LSD1WT/LSD1KO. See DEBUG_sample_info_unassigned_genotype.csv")
}

rownames(wm_sample_info) <- wm_sample_info$Sample

write.csv(wm_sample_info, file.path(wm_out_dir, "CHECK_corrected_sample_info_NO_MICHAEL.csv"), row.names = FALSE)

wm_dataset_colors <- c(
  "3 weeks" = "#1b9e77",
  "4 weeks" = "#d95f02",
  "Early onset" = "#7570b3",
  "Terminal" = "#e7298a"
)

# 16) GLOBAL normalize + ComBat

wm_combined_log <- log2(wm_combined_counts + 1)

wm_combat_info <- wm_sample_info
wm_combat_info$Dataset <- factor(
  wm_combat_info$Dataset,
  levels = c("3weeks", "4weeks", "early_onset", "terminal")
)
wm_combat_info$Genotype <- factor(
  wm_combat_info$Genotype,
  levels = c("LSD1WT", "LSD1KO")
)

wm_combat_mod <- model.matrix(~ Genotype, data = wm_combat_info)

wm_combined_combat <- sva::ComBat(
  dat = wm_combined_log,
  batch = wm_combat_info$Dataset,
  mod = wm_combat_mod,
  par.prior = TRUE,
  prior.plots = FALSE
)

wm_combined_z <- wm_row_zscore(wm_combined_combat)

wm_combined_immune_z <- wm_combined_z[
  intersect(rownames(wm_combined_z), wm_immune_genes),
  ,
  drop = FALSE
]

if (nrow(wm_combined_immune_z) == 0) {
  stop("No immune genes found in combined global matrix.")
}

# 16.25) AFTER-ComBat PCA / clustering QC

wm_qc_combat <- wm_combined_combat[
  intersect(rownames(wm_combined_combat), wm_immune_genes),
  ,
  drop = FALSE
]

if (nrow(wm_qc_combat) == 0) {
  stop("No immune genes found for after-ComBat QC PCA.")
}

wm_gene_var_combat <- apply(wm_qc_combat, 1, var, na.rm = TRUE)
wm_qc_combat <- wm_qc_combat[wm_gene_var_combat > 0, , drop = FALSE]

if (nrow(wm_qc_combat) == 0) {
  stop("No variable immune genes left for after-ComBat QC PCA.")
}

wm_pca_combat <- prcomp(t(wm_qc_combat), center = TRUE, scale. = TRUE)

wm_pca_combat_df <- as.data.frame(wm_pca_combat$x[, 1:min(4, ncol(wm_pca_combat$x)), drop = FALSE])
wm_pca_combat_df$Sample <- rownames(wm_pca_combat_df)
wm_pca_combat_df <- dplyr::left_join(wm_pca_combat_df, wm_sample_info, by = "Sample")

wm_percent_var_combat <- 100 * (wm_pca_combat$sdev^2 / sum(wm_pca_combat$sdev^2))

wm_p_pca_dataset_combat <- ggplot(
  wm_pca_combat_df,
  aes(x = PC1, y = PC2, color = DatasetDisplay, shape = Genotype)
) +
  geom_point(size = 4) +
  geom_text(aes(label = DisplayName), vjust = -0.8, size = 3, show.legend = FALSE) +
  labs(
    title = "AFTER ComBat: PCA of immune genes (NO Michael; Terminal = LSD1 terminal + Trem2WT_LSD1WT/LSD1KO)",
    x = paste0("PC1 (", round(wm_percent_var_combat[1], 1), "%)"),
    y = paste0("PC2 (", round(wm_percent_var_combat[2], 1), "%)")
  ) +
  theme_bw(base_size = 14)

wm_p_pca_genotype_combat <- ggplot(
  wm_pca_combat_df,
  aes(x = PC1, y = PC2, color = Genotype, shape = DatasetDisplay)
) +
  geom_point(size = 4) +
  geom_text(aes(label = DisplayName), vjust = -0.8, size = 3, show.legend = FALSE) +
  labs(
    title = "AFTER ComBat: PCA colored by genotype (NO Michael)",
    x = paste0("PC1 (", round(wm_percent_var_combat[1], 1), "%)"),
    y = paste0("PC2 (", round(wm_percent_var_combat[2], 1), "%)")
  ) +
  theme_bw(base_size = 14)

wm_sample_cor_combat <- cor(wm_qc_combat, method = "pearson", use = "pairwise.complete.obs")

wm_ha_qc_combat <- HeatmapAnnotation(
  Dataset = wm_sample_info[colnames(wm_sample_cor_combat), "DatasetDisplay"],
  Genotype = wm_sample_info[colnames(wm_sample_cor_combat), "Genotype"],
  TerminalSource = wm_sample_info[colnames(wm_sample_cor_combat), "TerminalSource"],
  col = list(
    Dataset = wm_dataset_colors,
    Genotype = c("LSD1WT" = "steelblue", "LSD1KO" = "firebrick")
  ),
  na_col = "grey90"
)

wm_ht_cor_combat <- Heatmap(
  wm_sample_cor_combat,
  name = "cor",
  top_annotation = wm_ha_qc_combat,
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  row_labels = wm_sample_info[rownames(wm_sample_cor_combat), "DisplayName"],
  column_labels = wm_sample_info[colnames(wm_sample_cor_combat), "DisplayName"],
  column_title = "AFTER ComBat: sample correlation heatmap (immune genes; NO Michael)"
)

wm_sample_dist_combat <- dist(t(wm_qc_combat))
wm_sample_dist_mat_combat <- as.matrix(wm_sample_dist_combat)

wm_ht_dist_combat <- Heatmap(
  wm_sample_dist_mat_combat,
  name = "distance",
  top_annotation = wm_ha_qc_combat,
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  row_labels = wm_sample_info[rownames(wm_sample_dist_mat_combat), "DisplayName"],
  column_labels = wm_sample_info[colnames(wm_sample_dist_mat_combat), "DisplayName"],
  column_title = "AFTER ComBat: sample distance heatmap (immune genes; NO Michael)"
)

pdf(
  file.path(wm_out_dir, "CORRECTED_COMBINE_terminal_QC_PCA_and_sample_clustering_immune_genes_AFTER_ComBat_NO_MICHAEL.pdf"),
  width = 15,
  height = 12
)
print(wm_p_pca_dataset_combat)
print(wm_p_pca_genotype_combat)
draw(wm_ht_cor_combat, newpage = TRUE)
draw(wm_ht_dist_combat, newpage = TRUE)
dev.off()

# 16.5) PCA BEFORE vs AFTER ComBat

wm_pca_genes <- intersect(rownames(wm_combined_log), wm_immune_genes)

if (length(wm_pca_genes) == 0) {
  stop("No immune genes found for before/after ComBat PCA.")
}

wm_mat_before <- wm_combined_log[wm_pca_genes, , drop = FALSE]
wm_mat_after  <- wm_combined_combat[wm_pca_genes, , drop = FALSE]

wm_var_before <- apply(wm_mat_before, 1, var, na.rm = TRUE)
wm_var_after  <- apply(wm_mat_after, 1, var, na.rm = TRUE)

wm_keep_genes <- names(wm_var_before)[wm_var_before > 0 & wm_var_after > 0]

wm_mat_before <- wm_mat_before[wm_keep_genes, , drop = FALSE]
wm_mat_after  <- wm_mat_after[wm_keep_genes, , drop = FALSE]

if (nrow(wm_mat_before) == 0 || nrow(wm_mat_after) == 0) {
  stop("No variable immune genes left for PCA after filtering.")
}

wm_pca_before <- prcomp(t(wm_mat_before), center = TRUE, scale. = TRUE)
wm_pca_before_df <- as.data.frame(wm_pca_before$x[, 1:min(4, ncol(wm_pca_before$x)), drop = FALSE])
wm_pca_before_df$Sample <- rownames(wm_pca_before_df)
wm_pca_before_df <- dplyr::left_join(wm_pca_before_df, wm_sample_info, by = "Sample")

wm_percent_var_before <- 100 * (wm_pca_before$sdev^2 / sum(wm_pca_before$sdev^2))

wm_p_before_dataset <- ggplot(wm_pca_before_df, aes(x = PC1, y = PC2, color = DatasetDisplay, shape = Genotype)) +
  geom_point(size = 4) +
  geom_text(aes(label = DisplayName), vjust = -0.8, size = 3, show.legend = FALSE) +
  labs(
    title = "BEFORE ComBat: PCA of immune genes (NO Michael)",
    x = paste0("PC1 (", round(wm_percent_var_before[1], 1), "%)"),
    y = paste0("PC2 (", round(wm_percent_var_before[2], 1), "%)")
  ) +
  theme_bw(base_size = 14)

wm_p_before_genotype <- ggplot(wm_pca_before_df, aes(x = PC1, y = PC2, color = Genotype, shape = DatasetDisplay)) +
  geom_point(size = 4) +
  geom_text(aes(label = DisplayName), vjust = -0.8, size = 3, show.legend = FALSE) +
  labs(
    title = "BEFORE ComBat: PCA colored by genotype (NO Michael)",
    x = paste0("PC1 (", round(wm_percent_var_before[1], 1), "%)"),
    y = paste0("PC2 (", round(wm_percent_var_before[2], 1), "%)")
  ) +
  theme_bw(base_size = 14)

wm_pca_after <- prcomp(t(wm_mat_after), center = TRUE, scale. = TRUE)
wm_pca_after_df <- as.data.frame(wm_pca_after$x[, 1:min(4, ncol(wm_pca_after$x)), drop = FALSE])
wm_pca_after_df$Sample <- rownames(wm_pca_after_df)
wm_pca_after_df <- dplyr::left_join(wm_pca_after_df, wm_sample_info, by = "Sample")

wm_percent_var_after <- 100 * (wm_pca_after$sdev^2 / sum(wm_pca_after$sdev^2))

wm_p_after_dataset <- ggplot(wm_pca_after_df, aes(x = PC1, y = PC2, color = DatasetDisplay, shape = Genotype)) +
  geom_point(size = 4) +
  geom_text(aes(label = DisplayName), vjust = -0.8, size = 3, show.legend = FALSE) +
  labs(
    title = "AFTER ComBat: PCA of immune genes (NO Michael)",
    x = paste0("PC1 (", round(wm_percent_var_after[1], 1), "%)"),
    y = paste0("PC2 (", round(wm_percent_var_after[2], 1), "%)")
  ) +
  theme_bw(base_size = 14)

wm_p_after_genotype <- ggplot(wm_pca_after_df, aes(x = PC1, y = PC2, color = Genotype, shape = DatasetDisplay)) +
  geom_point(size = 4) +
  geom_text(aes(label = DisplayName), vjust = -0.8, size = 3, show.legend = FALSE) +
  labs(
    title = "AFTER ComBat: PCA colored by genotype (NO Michael)",
    x = paste0("PC1 (", round(wm_percent_var_after[1], 1), "%)"),
    y = paste0("PC2 (", round(wm_percent_var_after[2], 1), "%)")
  ) +
  theme_bw(base_size = 14)

pdf(
  file.path(wm_out_dir, "CORRECTED_COMBINE_terminal_PCA_before_vs_after_ComBat_immune_genes_NO_MICHAEL.pdf"),
  width = 10,
  height = 8
)
print(wm_p_before_dataset)
print(wm_p_before_genotype)
print(wm_p_after_dataset)
print(wm_p_after_genotype)
dev.off()

# 17) GLOBAL immune row order

if (nrow(wm_combined_immune_z) == 1) {
  wm_global_immune_row_order <- rownames(wm_combined_immune_z)
} else {
  wm_hc_global <- hclust(dist(wm_combined_immune_z), method = "ward.D2")
  wm_global_immune_row_order <- rownames(wm_combined_immune_z)[wm_hc_global$order]
}

write.csv(
  data.frame(Gene = wm_global_immune_row_order, RowIndex = seq_along(wm_global_immune_row_order)),
  file.path(wm_out_dir, "CHECK_global_immune_gene_order_CORRECTED_NO_MICHAEL.csv"),
  row.names = FALSE
)

# 18) Split back to datasets

wm_wk3_z        <- wm_combined_z[, colnames(wm_wk3_full),        drop = FALSE]
wm_wk4_z        <- wm_combined_z[, colnames(wm_wk4_full),        drop = FALSE]
wm_early_z      <- wm_combined_z[, colnames(wm_early_full),      drop = FALSE]
wm_lsd1_term_z   <- wm_combined_z[, colnames(wm_term_full),       drop = FALSE]
wm_trem2_term_z <- wm_combined_z[, colnames(wm_trem2_term_full), drop = FALSE]

wm_term_z <- cbind(wm_lsd1_term_z, wm_trem2_term_z)

# Keep terminal columns grouped by genotype:
# LSD1WT first, then LSD1KO. Trem2WT_LSD1KO samples go inside LSD1KO.
wm_terminal_col_order <- c(
  colnames(wm_term_z)[
    grepl("Trem2WT_LSD1WT", colnames(wm_term_z), ignore.case = TRUE) |
      grepl("LSD1\\s*WT|LSD1[_ ]?WT|Cre\\-", colnames(wm_term_z), ignore.case = TRUE)
  ],
  colnames(wm_term_z)[
    grepl("Trem2WT_LSD1KO", colnames(wm_term_z), ignore.case = TRUE) |
      grepl("LSD1\\s*KO|LSD1[_ ]?KO|LSD1[_ ]?Del|LSD1\\s*Del|Cre\\+", colnames(wm_term_z), ignore.case = TRUE)
  ]
)
wm_terminal_col_order <- unique(wm_terminal_col_order)
wm_term_z <- wm_term_z[, wm_terminal_col_order, drop = FALSE]

# 19) One shared legend

wm_global_col_fun <- circlize::colorRamp2(
  c(-4, -2, 0, 2, 4),
  c("blue", "lightblue", "white", "salmon", "red")
)

# 20) Module display names

wm_module_display_name <- function(x) {
  dplyr::case_when(
    x == "Yellow"      ~ "Microglial and immune genes",
    x == "Light cyan"  ~ "Light cyan",
    x == "Forestgreen" ~ "Forestgreen",
    x == "Gold"        ~ "Gold",
    TRUE ~ x
  )
}

# 21) Output PDF

wm_out_pdf <- file.path(
  wm_out_dir,
  "CORRECTED_COMBINE_terminal_Aligned_Timepoints_Immune_LSD1terminal_plus_Trem2WT_LSD1WT_and_LSD1KO_GLOBAL_withGroupMean_NO_MICHAEL.pdf"
)

pdf(wm_out_pdf, width = 24, height = 8)

for (wm_m in wm_target_modules) {
  message("-> Plotting module: ", wm_m)

  wm_m_display <- wm_module_display_name(wm_m)

  wm_genes_m <- wm_gene2module %>%
    dplyr::filter(Module == wm_m) %>%
    dplyr::pull(Gene_Symbol) %>%
    unique()

  wm_shared_module_order <- wm_global_immune_row_order[wm_global_immune_row_order %in% wm_genes_m]
  if (length(wm_shared_module_order) == 0) next

  wm_z3 <- wm_wk3_z[intersect(wm_shared_module_order, rownames(wm_wk3_z)), , drop = FALSE]
  wm_z4 <- wm_wk4_z[intersect(wm_shared_module_order, rownames(wm_wk4_z)), , drop = FALSE]
  wm_ze <- wm_early_z[intersect(wm_shared_module_order, rownames(wm_early_z)), , drop = FALSE]
  wm_zt <- wm_term_z[intersect(wm_shared_module_order, rownames(wm_term_z)), , drop = FALSE]

  if (nrow(wm_z3) == 0 && nrow(wm_z4) == 0 && nrow(wm_ze) == 0 && nrow(wm_zt) == 0) next

  if (nrow(wm_z3) == 0) wm_z3 <- matrix(0, nrow = 1, ncol = ncol(wm_wk3_z),   dimnames = list("NO_GENE", colnames(wm_wk3_z)))
  if (nrow(wm_z4) == 0) wm_z4 <- matrix(0, nrow = 1, ncol = ncol(wm_wk4_z),   dimnames = list("NO_GENE", colnames(wm_wk4_z)))
  if (nrow(wm_ze) == 0) wm_ze <- matrix(0, nrow = 1, ncol = ncol(wm_early_z), dimnames = list("NO_GENE", colnames(wm_early_z)))
  if (nrow(wm_zt) == 0) wm_zt <- matrix(0, nrow = 1, ncol = ncol(wm_term_z),  dimnames = list("NO_GENE", colnames(wm_term_z)))

  wm_z3_g <- wm_collapse_group_mean_from_z(wm_z3)
  wm_z4_g <- wm_collapse_group_mean_from_z(wm_z4)
  wm_ze_g <- wm_collapse_group_mean_from_z(wm_ze)
  wm_zt_g <- wm_collapse_group_mean_from_z(wm_zt)

  wm_ha3 <- wm_make_sample_anno_cre(colnames(wm_z3))
  wm_ha4 <- wm_make_sample_anno_cre(colnames(wm_z4))
  wm_hae <- wm_make_sample_anno_cre(colnames(wm_ze))
  wm_hat <- wm_make_sample_anno_cre(colnames(wm_zt))

  wm_ha3_g <- wm_make_groupmean_anno_cre(colnames(wm_z3_g))
  wm_ha4_g <- wm_make_groupmean_anno_cre(colnames(wm_z4_g))
  wm_hae_g <- wm_make_groupmean_anno_cre(colnames(wm_ze_g))
  wm_hat_g <- wm_make_groupmean_anno_cre(colnames(wm_zt_g))

  wm_ht3 <- Heatmap(
    wm_z3,
    name = "z-score",
    col = wm_global_col_fun,
    cluster_rows = FALSE,
    row_order = seq_len(nrow(wm_z3)),
    cluster_columns = FALSE,
    column_order = colnames(wm_z3),
    column_labels = wm_make_display_labels(colnames(wm_z3)),
    show_column_dend = FALSE,
    top_annotation = wm_ha3,
    show_row_names = TRUE,
    row_names_gp = grid::gpar(fontsize = 6),
    column_title = paste0("3 weeks - ", wm_m_display)
  )

  wm_ht4 <- Heatmap(
    wm_z4,
    name = "z-score",
    col = wm_global_col_fun,
    cluster_rows = FALSE,
    row_order = seq_len(nrow(wm_z4)),
    cluster_columns = FALSE,
    column_order = colnames(wm_z4),
    column_labels = wm_make_display_labels(colnames(wm_z4)),
    show_column_dend = FALSE,
    top_annotation = wm_ha4,
    show_row_names = FALSE,
    column_title = paste0("4 weeks - ", wm_m_display)
  )

  wm_hte <- Heatmap(
    wm_ze,
    name = "z-score",
    col = wm_global_col_fun,
    cluster_rows = FALSE,
    row_order = seq_len(nrow(wm_ze)),
    cluster_columns = FALSE,
    column_order = colnames(wm_ze),
    column_labels = wm_make_display_labels(colnames(wm_ze)),
    show_column_dend = FALSE,
    top_annotation = wm_hae,
    show_row_names = FALSE,
    column_title = paste0("Early onset - ", wm_m_display)
  )

  wm_htt <- Heatmap(
    wm_zt,
    name = "z-score",
    col = wm_global_col_fun,
    cluster_rows = FALSE,
    row_order = seq_len(nrow(wm_zt)),
    cluster_columns = FALSE,
    column_order = colnames(wm_zt),
    column_labels = wm_make_display_labels(colnames(wm_zt)),
    show_column_dend = FALSE,
    top_annotation = wm_hat,
    show_row_names = FALSE,
    column_title = paste0("Terminal = LSD1 terminal + Trem2WT_LSD1WT/LSD1KO; NO Michael - ", wm_m_display)
  )

  draw(
    wm_ht3 + wm_ht4 + wm_hte + wm_htt,
    newpage = TRUE,
    merge_legends = TRUE,
    heatmap_legend_side = "right",
    annotation_legend_side = "right"
  )

  wm_ht3_g <- Heatmap(
    wm_z3_g,
    name = "z-score",
    col = wm_global_col_fun,
    cluster_rows = FALSE,
    row_order = seq_len(nrow(wm_z3_g)),
    cluster_columns = FALSE,
    column_order = colnames(wm_z3_g),
    show_column_dend = FALSE,
    top_annotation = wm_ha3_g,
    show_row_names = TRUE,
    row_names_gp = grid::gpar(fontsize = 6),
    column_title = paste0("3 weeks - ", wm_m_display)
  )

  wm_ht4_g <- Heatmap(
    wm_z4_g,
    name = "z-score",
    col = wm_global_col_fun,
    cluster_rows = FALSE,
    row_order = seq_len(nrow(wm_z4_g)),
    cluster_columns = FALSE,
    column_order = colnames(wm_z4_g),
    show_column_dend = FALSE,
    top_annotation = wm_ha4_g,
    show_row_names = FALSE,
    column_title = paste0("4 weeks - ", wm_m_display)
  )

  wm_hte_g <- Heatmap(
    wm_ze_g,
    name = "z-score",
    col = wm_global_col_fun,
    cluster_rows = FALSE,
    row_order = seq_len(nrow(wm_ze_g)),
    cluster_columns = FALSE,
    column_order = colnames(wm_ze_g),
    show_column_dend = FALSE,
    top_annotation = wm_hae_g,
    show_row_names = FALSE,
    column_title = paste0("Early onset - ", wm_m_display)
  )

  wm_htt_g <- Heatmap(
    wm_zt_g,
    name = "z-score",
    col = wm_global_col_fun,
    cluster_rows = FALSE,
    row_order = seq_len(nrow(wm_zt_g)),
    cluster_columns = FALSE,
    column_order = colnames(wm_zt_g),
    show_column_dend = FALSE,
    top_annotation = wm_hat_g,
    show_row_names = FALSE,
    column_title = paste0("Terminal = LSD1 terminal + Trem2WT_LSD1WT/LSD1KO; NO Michael - ", wm_m_display)
  )

  draw(
    wm_ht3_g + wm_ht4_g + wm_hte_g + wm_htt_g,
    newpage = TRUE,
    merge_legends = TRUE,
    heatmap_legend_side = "right",
    annotation_legend_side = "right"
  )
}

dev.off()

cat("Saved main heatmap PDF:\n", wm_out_pdf, "\n")
cat("Saved QC/sample-check files in:\n", wm_out_dir, "\n")
cat("DONE. Corrected terminal = LSD1 terminal + Trem2WT_LSD1WT and Trem2WT_LSD1KO terminal add-on samples. Michael was not imported.\n")
