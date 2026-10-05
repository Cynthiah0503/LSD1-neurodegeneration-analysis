# 14_lsd1_timecourse_corrected_terminal_reactive_astrocyte_heatmap_no_michael.R
#
# Purpose: Corrected LSD1 time-course reactive-astrocyte heatmap excluding Michael terminal samples.
# Inputs are expected under data/processed/ or data/external/ relative to this repository.
# Outputs are written under results/ or script-defined subfolders.

# CORRECTED COMBINE TERMINAL reactive-astrocyte marker heatmap
#
# FIXED LOGIC:
#   3w + 4w + early onset + terminal are still from timepoint_lsd1.
#   Terminal is now:
#     1) timepoint_lsd1 Terminal feature counts.xlsx
#     2) Trem2heatmap feature counts.xlsx, Trem2WT background only:
#          LSD1WT_Trem2WW + LSD1KO_Trem2WW
#        displayed as Trem2WT_LSD1WT + Trem2WT_LSD1KO.
#
#   Michael terminal is completely excluded.
#
#   This script intentionally DOES NOT read:
#     featurecountsLSD1del vs. LSD1wt .xlsx
#   and it stops if any "mich" / "michael" sample accidentally enters.
#
# Output:
#     combine terminal_REACTIVE_ASTRO_CORRECTED_LSD1terminal_plus_Trem2WT_LSD1WT_and_LSD1KO_no_Michael/


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
  "combine terminal_REACTIVE_ASTRO_CORRECTED_LSD1terminal_plus_Trem2WT_LSD1WT_and_LSD1KO_no_Michael"
)
if (!dir.exists(wm_out_dir)) dir.create(wm_out_dir, recursive = TRUE)

wk3_file   <- file.path(lsd1_dir, "3 weeks feature counts.xlsx")
wk4_file   <- file.path(lsd1_dir, "4 weeks feature counts.xlsx")
early_file <- file.path(lsd1_dir, "early_onset_featurecounts.xlsx")
term_file  <- file.path(lsd1_dir, "Terminal feature counts.xlsx")

trem2_heatmap_file <- file.path(trem2_dir, "Trem2heatmap feature counts.xlsx")

# Hard safety block: never use Michael in this corrected workflow.
michael_file_that_must_not_be_used <- file.path(lsd1_dir, "featurecountsLSD1del vs. LSD1wt .xlsx")

required_files <- c(wk3_file, wk4_file, early_file, term_file, trem2_heatmap_file)
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

# 3) Reactive astrocyte markers
wm_reactive_astro_markers <- c(
  "Gfap", "Nes", "Synm", "Vim",
  "Aldoc", "Fabp7", "Maob", "Tspo",
  "Cryab", "Hspb1",
  "C3", "Chi3l1", "Lcn2", "Serpina3n", "Mt1", "Mt2", "Thbs1",
  "Nfatc3", "Nfatc4", "Ntrk2", "Il17ra", "S100b", "Sox9", "Stat3",
  "Slc1a3", "Slc1a2", "Kcnj10"
)
wm_reactive_astro_markers_up <- toupper(wm_reactive_astro_markers)

cat("Reactive astrocyte markers:\n")
print(wm_reactive_astro_markers_up)

# 4) Display-name helpers
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
    idx_non <- ave(
      seq_along(samples)[non_terminal],
      paste(dataset[non_terminal], genotype[non_terminal], sep = "__"),
      FUN = seq_along
    )
    labels[non_terminal] <- paste0(
      dataset[non_terminal], "_", genotype[non_terminal], "_", sprintf("%02d", idx_non)
    )
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
    grepl("^3w__", samples)           ~ "3 weeks",
    grepl("^4w__", samples)           ~ "4 weeks",
    grepl("^early__", samples)        ~ "Early onset",
    grepl("^term_lsd1__", samples)    ~ "Terminal",
    grepl("^term_trem2wt__", samples) ~ "Terminal",
    TRUE ~ NA_character_
  )
}

wm_terminal_source <- function(samples) {
  dplyr::case_when(
    grepl("^term_lsd1__", samples)    ~ "LSD1_terminal",
    grepl("^term_trem2wt__", samples) ~ "Trem2WT_terminal_addon",
    TRUE ~ NA_character_
  )
}

# 5) Generic featureCounts helpers
wm_aggregate_counts <- function(df, gene_col, sample_cols, context = "input file") {
  if (length(sample_cols) == 0) stop("No sample columns selected for ", context)

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

wm_prep_featurecounts_standard <- function(df, context = "standard LSD1 file") {
  gene_col <- tail(names(df), 1)
  sample_cols <- grep("^LSD1", names(df), value = TRUE)
  if (length(sample_cols) == 0) stop("No sample columns found starting with 'LSD1' in ", context)
  wm_aggregate_counts(df, gene_col, sample_cols, context)
}

wm_prep_featurecounts_early <- function(df) {
  if (!("Gene_ID" %in% names(df))) {
    stop("For early_onset_featurecounts.xlsx, column 'Gene_ID' was not found.")
  }

  sample_cols <- c(
    "0501_LSD1WT", "0502_LSD1WT", "0503_LSD1WT",
    "0504_LSD1KO", "0505_LSD1KO", "0506_LSD1KO"
  )

  missing_cols <- setdiff(sample_cols, names(df))
  if (length(missing_cols) > 0) {
    stop("Missing sample columns in early onset file: ", paste(missing_cols, collapse = ", "))
  }

  wm_aggregate_counts(df, "Gene_ID", sample_cols, "early_onset_featurecounts.xlsx")
}

# 6) Trem2heatmap group parsing; keep Trem2WT background only
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

  ifelse(!is.na(lsd1_part) & !is.na(trem2_part), paste(lsd1_part, trem2_part, sep = "_"), NA_character_)
}

wm_prep_trem2wt_terminal_addon <- function(df) {
  # EXACT Trem2 import logic from your working Trem2 heatmap:
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
  write.csv(mapping_df, file.path(wm_out_dir, "CHECK_Trem2heatmap_all_sample_group_mapping.csv"), row.names = FALSE)

  # Correct target for this combine-terminal correction:
  # keep the Trem2WT background, BOTH LSD1WT and LSD1KO.
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

# 7) Matrix helpers
wm_dfagg_to_mat <- function(df_agg, sample_cols) {
  mat <- as.matrix(df_agg[, sample_cols, drop = FALSE])
  rownames(mat) <- df_agg$Gene_name
  storage.mode(mat) <- "numeric"
  mat
}

wm_make_full_union_mat <- function(mat, all_genes) {
  out <- matrix(0, nrow = length(all_genes), ncol = ncol(mat), dimnames = list(all_genes, colnames(mat)))
  out[rownames(mat), ] <- mat
  out
}

wm_row_zscore <- function(mat) {
  if (nrow(mat) == 0) return(mat)
  z <- t(scale(t(mat)))
  z[!is.finite(z)] <- 0
  z
}

# 8) Annotation and group mean helpers
wm_group_colors_cre <- c("LSD1KO" = "firebrick", "LSD1WT" = "steelblue")

wm_make_sample_anno_cre <- function(samples) {
  is_ko <- grepl("Trem2WT_LSD1KO", samples, ignore.case = TRUE) |
    grepl("LSD1\\s*KO|LSD1[_ ]?KO|LSD1[_ ]?Del|LSD1\\s*Del|Cre\\+", samples, ignore.case = TRUE)
  is_wt <- grepl("Trem2WT_LSD1WT", samples, ignore.case = TRUE) |
    grepl("LSD1\\s*WT|LSD1[_ ]?WT|Cre\\-", samples, ignore.case = TRUE)

  grp <- dplyr::case_when(is_ko ~ "LSD1KO", is_wt ~ "LSD1WT", TRUE ~ NA_character_)

  HeatmapAnnotation(
    Group  = grp,
    col    = list(Group = wm_group_colors_cre),
    na_col = "grey90",
    show_annotation_name = TRUE
  )
}

wm_get_lsd1_group <- function(samples) {
  dplyr::case_when(
    grepl("Trem2WT_LSD1KO", samples, ignore.case = TRUE) ~ "LSD1KO",
    grepl("Trem2WT_LSD1WT", samples, ignore.case = TRUE) ~ "LSD1WT",
    grepl("LSD1\\s*KO|LSD1[_ ]?KO|LSD1[_ ]?Del|LSD1\\s*Del|Cre\\+", samples, ignore.case = TRUE) ~ "LSD1KO",
    grepl("LSD1\\s*WT|LSD1[_ ]?WT|Cre\\-", samples, ignore.case = TRUE) ~ "LSD1WT",
    TRUE ~ NA_character_
  )
}

wm_make_groupmean_anno_cre <- function(group_cols) {
  grp <- wm_get_lsd1_group(group_cols)
  HeatmapAnnotation(
    Group = grp,
    col   = list(Group = wm_group_colors_cre),
    na_col = "grey90",
    show_annotation_name = TRUE
  )
}

wm_collapse_group_mean_from_z <- function(mat_z) {
  if (nrow(mat_z) == 0) return(mat_z)

  groups <- wm_get_lsd1_group(colnames(mat_z))
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

# 9) Prepare all datasets
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

all_raw_sample_names <- c(colnames(wm_wk3_mat), colnames(wm_wk4_mat), colnames(wm_early_mat), colnames(wm_term_mat), colnames(wm_trem2_term_mat))
if (any(grepl("mich|michael", all_raw_sample_names, ignore.case = TRUE))) {
  stop("Michael/Mich sample detected. Stopping because this corrected workflow must exclude Michael.")
}

# 10) Build GLOBAL matrix
wm_all_genes_union <- unique(c(
  rownames(wm_wk3_mat), rownames(wm_wk4_mat), rownames(wm_early_mat),
  rownames(wm_term_mat), rownames(wm_trem2_term_mat)
))

wm_wk3_full        <- wm_make_full_union_mat(wm_wk3_mat,        wm_all_genes_union)
wm_wk4_full        <- wm_make_full_union_mat(wm_wk4_mat,        wm_all_genes_union)
wm_early_full      <- wm_make_full_union_mat(wm_early_mat,      wm_all_genes_union)
wm_term_full       <- wm_make_full_union_mat(wm_term_mat,       wm_all_genes_union)
wm_trem2_term_full <- wm_make_full_union_mat(wm_trem2_term_mat, wm_all_genes_union)

wm_combined_counts <- cbind(wm_wk3_full, wm_wk4_full, wm_early_full, wm_term_full, wm_trem2_term_full)

if (any(grepl("mich|michael", colnames(wm_combined_counts), ignore.case = TRUE))) {
  stop("Michael/Mich sample detected after combining. Stopping.")
}

wm_reactive_genes_present <- intersect(rownames(wm_combined_counts), wm_reactive_astro_markers_up)
if (length(wm_reactive_genes_present) == 0) {
  stop("None of the reactive astrocyte marker genes were found in the combined matrix.")
}

cat("Reactive astrocyte markers found:\n")
print(wm_reactive_genes_present)
write.csv(
  data.frame(
    marker_requested = wm_reactive_astro_markers_up,
    found_in_combined_matrix = wm_reactive_astro_markers_up %in% wm_reactive_genes_present
  ),
  file.path(wm_out_dir, "CHECK_reactive_astrocyte_markers_found.csv"),
  row.names = FALSE
)

# 11) Sample information
wm_sample_info <- data.frame(Sample = colnames(wm_combined_counts), stringsAsFactors = FALSE)

wm_sample_info$Dataset <- dplyr::case_when(
  grepl("^3w__", wm_sample_info$Sample)           ~ "3weeks",
  grepl("^4w__", wm_sample_info$Sample)           ~ "4weeks",
  grepl("^early__", wm_sample_info$Sample)        ~ "early_onset",
  grepl("^term_lsd1__", wm_sample_info$Sample)    ~ "terminal",
  grepl("^term_trem2wt__", wm_sample_info$Sample) ~ "terminal",
  TRUE ~ NA_character_
)
wm_sample_info$DatasetDisplay <- wm_dataset_display(wm_sample_info$Sample)
wm_sample_info$TerminalSource <- wm_terminal_source(wm_sample_info$Sample)
wm_sample_info$DisplayName    <- wm_make_display_labels(wm_sample_info$Sample)
wm_sample_info$Genotype       <- wm_get_lsd1_group(wm_sample_info$Sample)

if (any(is.na(wm_sample_info$Dataset))) stop("Some samples could not be assigned to Dataset.")
if (any(is.na(wm_sample_info$Genotype))) {
  write.csv(wm_sample_info, file.path(wm_out_dir, "DEBUG_sample_info_unassigned_genotype.csv"), row.names = FALSE)
  stop("Some samples could not be assigned to LSD1WT/LSD1KO. See DEBUG_sample_info_unassigned_genotype.csv")
}

rownames(wm_sample_info) <- wm_sample_info$Sample
write.csv(wm_sample_info, file.path(wm_out_dir, "CHECK_corrected_reactive_astro_sample_info_NO_MICHAEL.csv"), row.names = FALSE)

wm_dataset_colors <- c(
  "3 weeks" = "#1b9e77",
  "4 weeks" = "#d95f02",
  "Early onset" = "#7570b3",
  "Terminal" = "#e7298a"
)

# 12) GLOBAL normalize + ComBat + row z-score
wm_combined_log <- log2(wm_combined_counts + 1)

wm_combat_info <- wm_sample_info
wm_combat_info$Dataset <- factor(wm_combat_info$Dataset, levels = c("3weeks", "4weeks", "early_onset", "terminal"))
wm_combat_info$Genotype <- factor(wm_combat_info$Genotype, levels = c("LSD1WT", "LSD1KO"))

wm_combat_mod <- model.matrix(~ Genotype, data = wm_combat_info)

wm_combined_combat <- sva::ComBat(
  dat = wm_combined_log,
  batch = wm_combat_info$Dataset,
  mod = wm_combat_mod,
  par.prior = TRUE,
  prior.plots = FALSE
)

wm_combined_z <- wm_row_zscore(wm_combined_combat)
wm_combined_reactive_z <- wm_combined_z[wm_reactive_genes_present, , drop = FALSE]

if (nrow(wm_combined_reactive_z) == 0) {
  stop("No reactive astrocyte markers found after processing.")
}

# 13) AFTER-ComBat QC PCA / clustering using reactive markers
wm_qc_combat <- wm_combined_combat[wm_reactive_genes_present, , drop = FALSE]
wm_gene_var_combat <- apply(wm_qc_combat, 1, var, na.rm = TRUE)
wm_qc_combat <- wm_qc_combat[wm_gene_var_combat > 0, , drop = FALSE]

if (nrow(wm_qc_combat) >= 2) {
  wm_pca_combat <- prcomp(t(wm_qc_combat), center = TRUE, scale. = TRUE)
  wm_pca_combat_df <- as.data.frame(wm_pca_combat$x[, 1:min(4, ncol(wm_pca_combat$x)), drop = FALSE])
  wm_pca_combat_df$Sample <- rownames(wm_pca_combat_df)
  wm_pca_combat_df <- dplyr::left_join(wm_pca_combat_df, wm_sample_info, by = "Sample")
  wm_percent_var_combat <- 100 * (wm_pca_combat$sdev^2 / sum(wm_pca_combat$sdev^2))

  wm_p_pca_dataset_combat <- ggplot(wm_pca_combat_df, aes(x = PC1, y = PC2, color = DatasetDisplay, shape = Genotype)) +
    geom_point(size = 4) +
    geom_text(aes(label = DisplayName), vjust = -0.8, size = 3, show.legend = FALSE) +
    labs(
      title = "AFTER ComBat: PCA of reactive astrocyte markers (NO Michael; Terminal = LSD1 terminal + Trem2WT_LSD1WT/LSD1KO)",
      x = paste0("PC1 (", round(wm_percent_var_combat[1], 1), "%)"),
      y = paste0("PC2 (", round(wm_percent_var_combat[2], 1), "%)")
    ) +
    theme_bw(base_size = 14)

  wm_p_pca_genotype_combat <- ggplot(wm_pca_combat_df, aes(x = PC1, y = PC2, color = Genotype, shape = DatasetDisplay)) +
    geom_point(size = 4) +
    geom_text(aes(label = DisplayName), vjust = -0.8, size = 3, show.legend = FALSE) +
    labs(
      title = "AFTER ComBat: PCA colored by genotype (reactive astrocyte markers; NO Michael)",
      x = paste0("PC1 (", round(wm_percent_var_combat[1], 1), "%)"),
      y = paste0("PC2 (", round(wm_percent_var_combat[2], 1), "%)")
    ) +
    theme_bw(base_size = 14)

  wm_sample_cor_combat <- cor(wm_qc_combat, method = "pearson", use = "pairwise.complete.obs")
  wm_ha_qc_combat <- HeatmapAnnotation(
    Dataset = wm_sample_info[colnames(wm_sample_cor_combat), "DatasetDisplay"],
    Genotype = wm_sample_info[colnames(wm_sample_cor_combat), "Genotype"],
    TerminalSource = wm_sample_info[colnames(wm_sample_cor_combat), "TerminalSource"],
    col = list(Dataset = wm_dataset_colors, Genotype = c("LSD1WT" = "steelblue", "LSD1KO" = "firebrick")),
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
    column_title = "AFTER ComBat: sample correlation heatmap (reactive astrocyte markers; NO Michael)"
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
    column_title = "AFTER ComBat: sample distance heatmap (reactive astrocyte markers; NO Michael)"
  )

  pdf(
    file.path(wm_out_dir, "CORRECTED_COMBINE_terminal_QC_PCA_and_sample_clustering_reactive_astro_markers_AFTER_ComBat_NO_MICHAEL.pdf"),
    width = 15,
    height = 12
  )
  print(wm_p_pca_dataset_combat)
  print(wm_p_pca_genotype_combat)
  draw(wm_ht_cor_combat, newpage = TRUE)
  draw(wm_ht_dist_combat, newpage = TRUE)
  dev.off()
}

# 14) Split back to datasets
wm_wk3_z        <- wm_combined_z[, colnames(wm_wk3_full),        drop = FALSE]
wm_wk4_z        <- wm_combined_z[, colnames(wm_wk4_full),        drop = FALSE]
wm_early_z      <- wm_combined_z[, colnames(wm_early_full),      drop = FALSE]
wm_lsd1_term_z  <- wm_combined_z[, colnames(wm_term_full),       drop = FALSE]
wm_trem2_term_z <- wm_combined_z[, colnames(wm_trem2_term_full), drop = FALSE]
wm_term_z <- cbind(wm_lsd1_term_z, wm_trem2_term_z)

# Keep terminal columns grouped by genotype: LSD1WT first, then LSD1KO.
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

# 15) Extract reactive marker matrices
wm_z3 <- wm_wk3_z[intersect(wm_reactive_astro_markers_up, rownames(wm_wk3_z)), , drop = FALSE]
wm_z4 <- wm_wk4_z[intersect(wm_reactive_astro_markers_up, rownames(wm_wk4_z)), , drop = FALSE]
wm_ze <- wm_early_z[intersect(wm_reactive_astro_markers_up, rownames(wm_early_z)), , drop = FALSE]
wm_zt <- wm_term_z[intersect(wm_reactive_astro_markers_up, rownames(wm_term_z)), , drop = FALSE]

wm_z3_g <- wm_collapse_group_mean_from_z(wm_z3)
wm_z4_g <- wm_collapse_group_mean_from_z(wm_z4)
wm_ze_g <- wm_collapse_group_mean_from_z(wm_ze)
wm_zt_g <- wm_collapse_group_mean_from_z(wm_zt)

wm_common_genes_all <- Reduce(
  intersect,
  list(rownames(wm_z3), rownames(wm_z4), rownames(wm_ze), rownames(wm_zt), rownames(wm_z3_g), rownames(wm_z4_g), rownames(wm_ze_g), rownames(wm_zt_g))
)

if (length(wm_common_genes_all) == 0) {
  stop("No shared reactive astrocyte marker genes across all datasets/pages.")
}

wm_z3   <- wm_z3[wm_common_genes_all, , drop = FALSE]
wm_z4   <- wm_z4[wm_common_genes_all, , drop = FALSE]
wm_ze   <- wm_ze[wm_common_genes_all, , drop = FALSE]
wm_zt   <- wm_zt[wm_common_genes_all, , drop = FALSE]
wm_z3_g <- wm_z3_g[wm_common_genes_all, , drop = FALSE]
wm_z4_g <- wm_z4_g[wm_common_genes_all, , drop = FALSE]
wm_ze_g <- wm_ze_g[wm_common_genes_all, , drop = FALSE]
wm_zt_g <- wm_zt_g[wm_common_genes_all, , drop = FALSE]

# 16) ONE shared gene order across individual + group mean pages
wm_all_mix <- cbind(wm_z3, wm_z4, wm_ze, wm_zt, wm_z3_g, wm_z4_g, wm_ze_g, wm_zt_g)

if (nrow(wm_all_mix) == 1) {
  wm_gene_names_shared <- rownames(wm_all_mix)
} else {
  wm_hc_all <- hclust(dist(wm_all_mix), method = "ward.D2")
  wm_gene_names_shared <- rownames(wm_all_mix)[wm_hc_all$order]
}

write.csv(
  data.frame(Gene = wm_gene_names_shared, RowIndex = seq_along(wm_gene_names_shared)),
  file.path(wm_out_dir, "CHECK_reactive_astrocyte_gene_order_CORRECTED_NO_MICHAEL.csv"),
  row.names = FALSE
)

wm_z3   <- wm_z3[wm_gene_names_shared, , drop = FALSE]
wm_z4   <- wm_z4[wm_gene_names_shared, , drop = FALSE]
wm_ze   <- wm_ze[wm_gene_names_shared, , drop = FALSE]
wm_zt   <- wm_zt[wm_gene_names_shared, , drop = FALSE]
wm_z3_g <- wm_z3_g[wm_gene_names_shared, , drop = FALSE]
wm_z4_g <- wm_z4_g[wm_gene_names_shared, , drop = FALSE]
wm_ze_g <- wm_ze_g[wm_gene_names_shared, , drop = FALSE]
wm_zt_g <- wm_zt_g[wm_gene_names_shared, , drop = FALSE]

# 17) Shared color scale + annotations
wm_global_col_fun <- circlize::colorRamp2(
  c(-4, -2, 0, 2, 4),
  c("blue", "lightblue", "white", "salmon", "red")
)

wm_ha3   <- wm_make_sample_anno_cre(colnames(wm_z3))
wm_ha4   <- wm_make_sample_anno_cre(colnames(wm_z4))
wm_hae   <- wm_make_sample_anno_cre(colnames(wm_ze))
wm_hat   <- wm_make_sample_anno_cre(colnames(wm_zt))

wm_ha3_g <- wm_make_groupmean_anno_cre(colnames(wm_z3_g))
wm_ha4_g <- wm_make_groupmean_anno_cre(colnames(wm_z4_g))
wm_hae_g <- wm_make_groupmean_anno_cre(colnames(wm_ze_g))
wm_hat_g <- wm_make_groupmean_anno_cre(colnames(wm_zt_g))

# 18) Heatmaps
wm_ht3 <- Heatmap(
  wm_z3,
  name = "z-score",
  col = wm_global_col_fun,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  row_order = seq_len(nrow(wm_z3)),
  column_order = colnames(wm_z3),
  column_labels = wm_make_display_labels(colnames(wm_z3)),
  show_column_dend = FALSE,
  top_annotation = wm_ha3,
  show_row_names = TRUE,
  row_names_gp = grid::gpar(fontsize = 8),
  column_title = "3 weeks - Reactive astrocyte markers"
)

wm_ht4 <- Heatmap(
  wm_z4,
  name = "z-score",
  col = wm_global_col_fun,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  row_order = seq_len(nrow(wm_z4)),
  column_order = colnames(wm_z4),
  column_labels = wm_make_display_labels(colnames(wm_z4)),
  show_column_dend = FALSE,
  top_annotation = wm_ha4,
  show_row_names = FALSE,
  column_title = "4 weeks - Reactive astrocyte markers"
)

wm_hte <- Heatmap(
  wm_ze,
  name = "z-score",
  col = wm_global_col_fun,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  row_order = seq_len(nrow(wm_ze)),
  column_order = colnames(wm_ze),
  column_labels = wm_make_display_labels(colnames(wm_ze)),
  show_column_dend = FALSE,
  top_annotation = wm_hae,
  show_row_names = FALSE,
  column_title = "Early onset - Reactive astrocyte markers"
)

wm_htt <- Heatmap(
  wm_zt,
  name = "z-score",
  col = wm_global_col_fun,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  row_order = seq_len(nrow(wm_zt)),
  column_order = colnames(wm_zt),
  column_labels = wm_make_display_labels(colnames(wm_zt)),
  show_column_dend = FALSE,
  top_annotation = wm_hat,
  show_row_names = TRUE,
  row_names_gp = grid::gpar(fontsize = 8),
  column_title = "Terminal - Reactive astrocyte markers"
)

wm_ht3_g <- Heatmap(
  wm_z3_g,
  name = "z-score",
  col = wm_global_col_fun,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  row_order = seq_len(nrow(wm_z3_g)),
  column_order = colnames(wm_z3_g),
  show_column_dend = FALSE,
  top_annotation = wm_ha3_g,
  show_row_names = TRUE,
  row_names_gp = grid::gpar(fontsize = 8),
  column_title = "3 weeks - Reactive astrocyte markers"
)

wm_ht4_g <- Heatmap(
  wm_z4_g,
  name = "z-score",
  col = wm_global_col_fun,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  row_order = seq_len(nrow(wm_z4_g)),
  column_order = colnames(wm_z4_g),
  show_column_dend = FALSE,
  top_annotation = wm_ha4_g,
  show_row_names = FALSE,
  column_title = "4 weeks - Reactive astrocyte markers"
)

wm_hte_g <- Heatmap(
  wm_ze_g,
  name = "z-score",
  col = wm_global_col_fun,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  row_order = seq_len(nrow(wm_ze_g)),
  column_order = colnames(wm_ze_g),
  show_column_dend = FALSE,
  top_annotation = wm_hae_g,
  show_row_names = FALSE,
  column_title = "Early onset - Reactive astrocyte markers"
)

wm_htt_g <- Heatmap(
  wm_zt_g,
  name = "z-score",
  col = wm_global_col_fun,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  row_order = seq_len(nrow(wm_zt_g)),
  column_order = colnames(wm_zt_g),
  show_column_dend = FALSE,
  top_annotation = wm_hat_g,
  show_row_names = FALSE,
  column_title = "Terminal - Reactive astrocyte markers"
)

wm_out_pdf <- file.path(
  wm_out_dir,
  "CORRECTED_COMBINE_terminal_Aligned_Timepoints_ReactiveAstrocytes_LSD1terminal_plus_Trem2WT_LSD1WT_and_LSD1KO_GLOBAL_withGroupMean_NO_MICHAEL.pdf"
)

pdf(wm_out_pdf, width = 24, height = 8)

draw(
  wm_ht3 + wm_ht4 + wm_hte + wm_htt,
  newpage = TRUE,
  merge_legends = TRUE,
  heatmap_legend_side = "right",
  annotation_legend_side = "right"
)

draw(
  wm_ht3_g + wm_ht4_g + wm_hte_g + wm_htt_g,
  newpage = TRUE,
  merge_legends = TRUE,
  heatmap_legend_side = "right",
  annotation_legend_side = "right"
)

dev.off()

cat("Saved:", wm_out_pdf, "\n")
cat("Reactive markers found:", length(wm_reactive_genes_present), "/", length(wm_reactive_astro_markers_up), "\n")
cat("Check files written to:", wm_out_dir, "\n")
