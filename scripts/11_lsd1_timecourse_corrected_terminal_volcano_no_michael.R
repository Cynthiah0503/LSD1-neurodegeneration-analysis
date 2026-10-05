# 11_lsd1_timecourse_corrected_terminal_volcano_no_michael.R
#
# Purpose: Corrected LSD1 time-course DESeq2 volcano workflow excluding Michael terminal samples.
# Inputs are expected under data/processed/ or data/external/ relative to this repository.
# Outputs are written under results/ or script-defined subfolders.

# CORRECTED GALAXY-STYLE FC COMBINE TERMINAL VOLCANO
#
#   Same volcano logic as the older combine-terminal script:
#     3 weeks KO vs WT
#     4 weeks KO vs WT
#     Early onset KO vs WT
#     Terminal KO vs WT
#
#   Corrected terminal = timepoint_lsd1 Terminal + Trem2WT-background terminal
#
# Corrected Terminal source:
#   1) timepoint_lsd1/Terminal feature counts.xlsx
#   2) Trem2 terminal experiment/Trem2heatmap feature counts.xlsx
#      KEEP ONLY:
#        LSD1WT_Trem2WW  = LSD1WT_Trem2WT control
#        LSD1KO_Trem2WW  = LSD1KO_Trem2WT experiment
#      EXCLUDE:
#        Trem2MM / Trem2KO background
#        Michael completely
#
# Direction:
#   log2FC = LSD1KO / LSD1WT
#   Positive log2FC = higher in LSD1KO
#   Negative log2FC = lower in LSD1KO
#
# DESeq2 settings:
#   DESeq(fitType = "parametric", betaPrior = TRUE)
#   results(contrast = c("Genotype", "LSD1KO", "LSD1WT"))
#
# Model:
#   3 weeks     -> ~ Genotype
#   4 weeks     -> ~ Genotype
#   Early onset -> ~ Genotype
#   Terminal    -> ~ Dataset + Genotype
#
# Output folder:
#     CORRECTED_volcano_LSD1terminal_plus_Trem2WTterminal_no_Michael_GALAXY_style_FC/


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
  library(dplyr)
  library(ggplot2)
  library(ggrepel)
  library(patchwork)
  library(readxl)
  library(readr)
  library(DESeq2)
})

select  <- dplyr::select
filter  <- dplyr::filter
mutate  <- dplyr::mutate
arrange <- dplyr::arrange

katz_dir <- repo_root

base_dir <- file.path(katz_dir, "timepoint_lsd1")
trem2_dir <- file.path(katz_dir, "Trem2 terminal experiment")

corrected_root <- file.path(base_dir, "正确combine terminal")

out_dir <- file.path(
  corrected_root,
  "CORRECTED_volcano_LSD1terminal_plus_Trem2WTterminal_no_Michael_GALAXY_style_FC"
)

if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}

# Raw featureCounts inputs
wk3_file   <- file.path(base_dir, "3 weeks feature counts.xlsx")
wk4_file   <- file.path(base_dir, "4 weeks feature counts.xlsx")
early_file <- file.path(base_dir, "early_onset_featurecounts.xlsx")
term_file  <- file.path(base_dir, "Terminal feature counts.xlsx")

# Correct Trem2 terminal raw matrix
trem2_heatmap_file <- file.path(trem2_dir, "Trem2heatmap feature counts.xlsx")

# Michael exists but must never be read
mich_file_that_must_not_be_used <- file.path(base_dir, "featurecountsLSD1del vs. LSD1wt .xlsx")

# Module annotation files
allc_file <- file.path(base_dir, "NIHMS472534-supplement-02.csv")
mod_file  <- file.path(base_dir, "modules.csv")

required_files <- c(
  wk3_file,
  wk4_file,
  early_file,
  term_file,
  trem2_heatmap_file,
  allc_file,
  mod_file
)

cat("Output folder:\n", out_dir, "\n\n")
cat("Required input files:\n")
print(required_files)
cat("\nFile existence check:\n")
print(file.exists(required_files))

missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files) > 0) {
  stop("Missing required input files:\n", paste(missing_files, collapse = "\n"))
}

if (file.exists(mich_file_that_must_not_be_used)) {
  message("Michael file exists on disk but will NOT be read: ", mich_file_that_must_not_be_used)
}

# 1) Helpers
fix_names <- function(df) {
  names(df) <- sub("^\\ufeff", "", names(df))
  bad <- is.na(names(df)) | names(df) == ""
  if (any(bad)) names(df)[bad] <- paste0("V", seq_len(sum(bad)))
  names(df) <- make.unique(names(df))
  as.data.frame(df)
}

read_xlsx_fixed <- function(path) {
  df <- readxl::read_excel(path, sheet = 1)
  fix_names(df)
}

clean_colname <- function(x) {
  x <- gsub("[^A-Za-z0-9]+", "_", x)
  x <- gsub("_+", "_", x)
  x <- gsub("^_|_$", "", x)
  x
}

pick_gene_symbol_col <- function(df, context = "input") {
  priority <- c(
    "Gene_name", "Gene_Name", "GeneSymbol", "Gene_Symbol", "gene_symbol",
    "Symbol", "SYMBOL", "external_gene_name", "Gene", "gene"
  )
  hit <- priority[priority %in% names(df)]
  if (length(hit) > 0) return(hit[1])

  # Older featureCounts-style files usually keep gene symbol as the last column.
  tail(names(df), 1)
}

aggregate_by_gene_symbol <- function(df, gene_col, sample_cols, context = "input") {
  if (length(sample_cols) == 0) {
    stop("No sample columns selected for ", context)
  }

  missing_cols <- setdiff(c(gene_col, sample_cols), names(df))
  if (length(missing_cols) > 0) {
    stop("Missing columns in ", context, ": ", paste(missing_cols, collapse = ", "))
  }

  df$Gene_name <- toupper(trimws(as.character(df[[gene_col]])))

  df_agg <- df[, c("Gene_name", sample_cols), drop = FALSE] %>%
    filter(!is.na(Gene_name), Gene_name != "") %>%
    group_by(Gene_name) %>%
    summarise(
      across(all_of(sample_cols), ~ sum(as.numeric(.x), na.rm = TRUE)),
      .groups = "drop"
    )

  mat <- as.matrix(df_agg[, sample_cols, drop = FALSE])
  rownames(mat) <- df_agg$Gene_name
  storage.mode(mat) <- "numeric"
  mat <- round(mat)
  mat[is.na(mat)] <- 0
  mat[mat < 0] <- 0

  mat
}

make_full_union_mat <- function(mat, all_genes) {
  out <- matrix(
    0,
    nrow = length(all_genes),
    ncol = ncol(mat),
    dimnames = list(all_genes, colnames(mat))
  )
  out[rownames(mat), ] <- mat
  out
}

make_group_id_4_from_samples <- function(samples) {
  s <- as.character(samples)

  is_lsd1wt <- grepl("LSD1\\s*WT|LSD1[_ ]?WT", s, ignore.case = TRUE)
  is_lsd1ko <- grepl(
    "LSD1\\s*KO|LSD1[_ ]?KO|LSD1[_ ]?Del|LSD1\\s*Del|LDS1[_ ]?Del|Cre\\+",
    s,
    ignore.case = TRUE
  )

  lsd1_part <- dplyr::case_when(
    is_lsd1wt ~ "LSD1WT",
    is_lsd1ko ~ "LSD1KO",
    TRUE ~ NA_character_
  )

  is_trem2ww <- grepl("Trem2WW|Trem2[_ ]?WT|Trem2WT", s, ignore.case = TRUE)
  is_trem2mm <- grepl("Trem2MM|Trem2[_ ]?KO|Trem2KO", s, ignore.case = TRUE)

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

wk3_raw   <- read_xlsx_fixed(wk3_file)
wk4_raw   <- read_xlsx_fixed(wk4_file)
early_raw <- read_xlsx_fixed(early_file)
term_raw  <- read_xlsx_fixed(term_file)
trem2_raw <- read_xlsx_fixed(trem2_heatmap_file)

# 3) Prep standard LSD1 featureCounts files
prep_featurecounts_standard <- function(df, dataset_label, timepoint_label) {
  gene_col <- tail(names(df), 1)
  sample_cols <- grep("^LSD1", names(df), value = TRUE)

  if (length(sample_cols) == 0) {
    stop("No sample columns starting with 'LSD1' were found for dataset: ", dataset_label)
  }

  mat <- aggregate_by_gene_symbol(
    df,
    gene_col = gene_col,
    sample_cols = sample_cols,
    context = dataset_label
  )

  colnames(mat) <- paste0(dataset_label, "__", colnames(mat))

  sample_info <- data.frame(
    Sample = colnames(mat),
    Dataset = dataset_label,
    Timepoint = timepoint_label,
    stringsAsFactors = FALSE
  )

  sample_info$Genotype <- dplyr::case_when(
    grepl("LSD1\\s*KO|LSD1[_ ]?KO", sample_info$Sample, ignore.case = TRUE) ~ "LSD1KO",
    grepl("LSD1\\s*WT|LSD1[_ ]?WT", sample_info$Sample, ignore.case = TRUE) ~ "LSD1WT",
    TRUE ~ NA_character_
  )

  if (any(is.na(sample_info$Genotype))) {
    stop(
      "Could not assign genotype for samples in ", dataset_label, ": ",
      paste(sample_info$Sample[is.na(sample_info$Genotype)], collapse = ", ")
    )
  }

  list(mat = mat, sample_info = sample_info)
}

# 4) Prep early onset merged featureCounts file
prep_featurecounts_early <- function(df, dataset_label = "early_onset", timepoint_label = "Early onset") {
  if (!("Gene_ID" %in% names(df))) {
    stop("For early_onset_featurecounts.xlsx, column 'Gene_ID' was not found.")
  }

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

  mat <- aggregate_by_gene_symbol(
    df,
    gene_col = "Gene_ID",
    sample_cols = sample_cols,
    context = "early_onset_featurecounts.xlsx"
  )

  colnames(mat) <- paste0(dataset_label, "__", colnames(mat))

  sample_info <- data.frame(
    Sample = colnames(mat),
    Dataset = dataset_label,
    Timepoint = timepoint_label,
    stringsAsFactors = FALSE
  )

  sample_info$Genotype <- dplyr::case_when(
    grepl("LSD1\\s*KO|LSD1[_ ]?KO", sample_info$Sample, ignore.case = TRUE) ~ "LSD1KO",
    grepl("LSD1\\s*WT|LSD1[_ ]?WT", sample_info$Sample, ignore.case = TRUE) ~ "LSD1WT",
    TRUE ~ NA_character_
  )

  if (any(is.na(sample_info$Genotype))) {
    stop(
      "Could not assign genotype for early onset samples: ",
      paste(sample_info$Sample[is.na(sample_info$Genotype)], collapse = ", ")
    )
  }

  list(mat = mat, sample_info = sample_info)
}

# 5) Prep corrected Trem2WT-background terminal samples
prep_trem2wt_terminal <- function(df, dataset_label = "Trem2WT_terminal", timepoint_label = "Terminal") {
  gene_col <- pick_gene_symbol_col(df, "Trem2heatmap featureCounts")
  sample_cols_all <- grep("^LSD1", names(df), value = TRUE)

  if (length(sample_cols_all) == 0) {
    stop("No Trem2heatmap sample columns found starting with 'LSD1'.")
  }

  group_id <- make_group_id_4_from_samples(sample_cols_all)

  mapping <- data.frame(
    original_column = sample_cols_all,
    group_id = group_id,
    stringsAsFactors = FALSE
  )

  write.csv(
    mapping,
    file.path(out_dir, "CHECK_Trem2heatmap_all_sample_group_mapping.csv"),
    row.names = FALSE
  )

  target_groups <- c("LSD1WT_Trem2WW", "LSD1KO_Trem2WW")

  selected_cols <- mapping$original_column[mapping$group_id %in% target_groups]

  if (length(selected_cols) == 0) {
    stop(
      "No Trem2WT-background samples were selected from Trem2heatmap.\n",
      "Open CHECK_Trem2heatmap_all_sample_group_mapping.csv and inspect group_id."
    )
  }

  selected_mapping <- mapping %>%
    filter(group_id %in% target_groups)

  write.csv(
    selected_mapping,
    file.path(out_dir, "CHECK_selected_Trem2WT_background_terminal_columns.csv"),
    row.names = FALSE
  )

  if (!all(target_groups %in% selected_mapping$group_id)) {
    warning(
      "Not all target Trem2WT-background groups were detected. Detected groups: ",
      paste(unique(selected_mapping$group_id), collapse = ", ")
    )
  }

  mat <- aggregate_by_gene_symbol(
    df,
    gene_col = gene_col,
    sample_cols = selected_cols,
    context = "Trem2heatmap Trem2WT-background featureCounts"
  )

  selected_group_id <- selected_mapping$group_id[
    match(selected_cols, selected_mapping$original_column)
  ]

  genotype <- dplyr::case_when(
    selected_group_id == "LSD1WT_Trem2WW" ~ "LSD1WT",
    selected_group_id == "LSD1KO_Trem2WW" ~ "LSD1KO",
    TRUE ~ NA_character_
  )

  if (any(is.na(genotype))) {
    stop("Could not assign Trem2WT-background selected columns to LSD1WT/LSD1KO.")
  }

  colnames(mat) <- paste0(
    dataset_label,
    "__",
    genotype,
    "__",
    clean_colname(selected_cols)
  )

  sample_info <- data.frame(
    Sample = colnames(mat),
    Dataset = dataset_label,
    Timepoint = timepoint_label,
    stringsAsFactors = FALSE
  )

  sample_info$Genotype <- genotype

  list(mat = mat, sample_info = sample_info)
}

# 6) Build matrices
wk3_p   <- prep_featurecounts_standard(wk3_raw,  dataset_label = "3weeks",        timepoint_label = "3 weeks")
wk4_p   <- prep_featurecounts_standard(wk4_raw,  dataset_label = "4weeks",        timepoint_label = "4 weeks")
early_p <- prep_featurecounts_early(early_raw,   dataset_label = "early_onset",   timepoint_label = "Early onset")
term_p  <- prep_featurecounts_standard(term_raw, dataset_label = "LSD1_terminal", timepoint_label = "Terminal")
trem2_p <- prep_trem2wt_terminal(trem2_raw,      dataset_label = "Trem2WT_terminal", timepoint_label = "Terminal")

all_genes_union <- unique(c(
  rownames(wk3_p$mat),
  rownames(wk4_p$mat),
  rownames(early_p$mat),
  rownames(term_p$mat),
  rownames(trem2_p$mat)
))

wk3_full   <- make_full_union_mat(wk3_p$mat,   all_genes_union)
wk4_full   <- make_full_union_mat(wk4_p$mat,   all_genes_union)
early_full <- make_full_union_mat(early_p$mat, all_genes_union)
term_full  <- make_full_union_mat(term_p$mat,  all_genes_union)
trem2_full <- make_full_union_mat(trem2_p$mat, all_genes_union)

combined_counts <- cbind(
  wk3_full,
  wk4_full,
  early_full,
  term_full,
  trem2_full
)

combined_sample_info <- bind_rows(
  wk3_p$sample_info,
  wk4_p$sample_info,
  early_p$sample_info,
  term_p$sample_info,
  trem2_p$sample_info
)

rownames(combined_sample_info) <- combined_sample_info$Sample
combined_sample_info <- combined_sample_info[colnames(combined_counts), , drop = FALSE]

# Hard safety block
if (any(grepl("mich|michael", combined_sample_info$Sample, ignore.case = TRUE)) ||
    any(grepl("mich|michael", combined_sample_info$Dataset, ignore.case = TRUE))) {
  stop("Michael/mich sample entered corrected volcano workflow. Stop.")
}

write.csv(
  combined_counts,
  file.path(out_dir, "CORRECTED_input_raw_featureCounts_gene_union_NO_MICHAEL.csv")
)
write.csv(
  combined_sample_info,
  file.path(out_dir, "CORRECTED_input_sample_info_NO_MICHAEL.csv"),
  row.names = FALSE
)

cat("\nCombined sample table:\n")
print(table(combined_sample_info$Timepoint, combined_sample_info$Dataset, combined_sample_info$Genotype))

# 7) DESeq2 runner
run_deseq_for_timepoint <- function(count_mat, sample_info, timepoint_label, output_prefix) {
  message("\n============================")
  message("Running Galaxy-style DESeq2 for: ", timepoint_label)
  message("============================")

  sample_info <- as.data.frame(sample_info)
  sample_info$Genotype <- factor(sample_info$Genotype, levels = c("LSD1KO", "LSD1WT"))
  sample_info$Dataset  <- factor(sample_info$Dataset)

  keep_samples <- sample_info$Timepoint == timepoint_label
  si <- droplevels(sample_info[keep_samples, , drop = FALSE])
  mat <- count_mat[, si$Sample, drop = FALSE]

  mat <- round(as.matrix(mat))
  mat[is.na(mat)] <- 0
  mat[mat < 0] <- 0
  storage.mode(mat) <- "integer"

  if (nrow(mat) == 0) {
    stop("No genes found for ", timepoint_label)
  }

  if (!all(c("LSD1WT", "LSD1KO") %in% unique(as.character(si$Genotype)))) {
    stop("Both LSD1WT and LSD1KO must be present for ", timepoint_label)
  }

  dataset_counts <- table(si$Dataset, si$Genotype)
  message("Sample table for ", timepoint_label, ":")
  print(dataset_counts)

  if (any(dataset_counts == 0) && nlevels(si$Dataset) > 1) {
    stop(
      "At least one Dataset x Genotype cell has zero samples for ", timepoint_label,
      ". Cannot safely run ~ Dataset + Genotype."
    )
  }

  if (nlevels(si$Dataset) > 1) {
    design_formula <- ~ Dataset + Genotype
    model_text <- "~ Dataset + Genotype"
    message("Model: ~ Dataset + Genotype")
  } else {
    design_formula <- ~ Genotype
    model_text <- "~ Genotype"
    message("Model: ~ Genotype  [single dataset; no identifiable batch term]")
  }

  dds <- DESeqDataSetFromMatrix(
    countData = mat,
    colData = si,
    design = design_formula
  )

  dds <- DESeq(
    dds,
    fitType   = "parametric",
    betaPrior = TRUE
  )

  res <- results(
    dds,
    contrast = c("Genotype", "LSD1KO", "LSD1WT"),
    alpha    = 0.1
  )

  res_df <- as.data.frame(res)
  res_df$Gene_name <- rownames(res_df)

  out <- res_df %>%
    transmute(
      Gene_name = Gene_name,
      baseMean = baseMean,
      log2.FC. = log2FoldChange,
      lfcSE = lfcSE,
      stat = stat,
      P.Value = pvalue,
      P.adj = padj
    ) %>%
    filter(!is.na(Gene_name), !is.na(log2.FC.), !is.na(P.adj)) %>%
    mutate(
      Significance = case_when(
        P.adj < 0.05 & log2.FC. > 1  ~ "Up regulated",
        P.adj < 0.05 & log2.FC. < -1 ~ "Down regulated",
        TRUE ~ "Not Significant"
      ),
      Timepoint = timepoint_label,
      Model = model_text,
      Contrast = "LSD1KO_vs_LSD1WT",
      DESeq2_style = "Galaxy_2.11.40.8_style_parametric_betaPrior_TRUE"
    ) %>%
    arrange(P.Value)

  write.csv(
    out,
    file.path(out_dir, paste0(output_prefix, "_GALAXY_STYLE_DESeq2_KO_vs_WT.csv")),
    row.names = FALSE
  )

  out
}

wk3_combine <- run_deseq_for_timepoint(
  combined_counts,
  combined_sample_info,
  timepoint_label = "3 weeks",
  output_prefix = "3weeks"
)

wk4_combine <- run_deseq_for_timepoint(
  combined_counts,
  combined_sample_info,
  timepoint_label = "4 weeks",
  output_prefix = "4weeks"
)

early_combine <- run_deseq_for_timepoint(
  combined_counts,
  combined_sample_info,
  timepoint_label = "Early onset",
  output_prefix = "Early_onset"
)

term_combine <- run_deseq_for_timepoint(
  combined_counts,
  combined_sample_info,
  timepoint_label = "Terminal",
  output_prefix = "Terminal_LSD1terminal_plus_Trem2WTterminal_NO_MICHAEL"
)

# 9) Label helper
get_top_labels_combine <- function(df, n_each = 8) {
  up <- df %>%
    filter(P.adj < 0.05, log2.FC. > 1) %>%
    arrange(P.adj) %>%
    slice_head(n = n_each)

  down <- df %>%
    filter(P.adj < 0.05, log2.FC. < -1) %>%
    arrange(P.adj) %>%
    slice_head(n = n_each)

  bind_rows(up, down) %>%
    distinct(Gene_name, .keep_all = TRUE)
}

# 10) Axis limits
get_ymax_combine <- function(df) {
  y <- -log10(pmax(df$P.adj, 1e-300))
  y <- y[is.finite(y)]
  if (length(y) == 0) return(10)
  max(10, ceiling(max(y) + 2))
}

global_xlim_combine <- c(-10, 10)

global_ymax_combine <- max(
  get_ymax_combine(early_combine),
  get_ymax_combine(wk3_combine),
  get_ymax_combine(wk4_combine),
  get_ymax_combine(term_combine)
)

global_ylim_combine <- c(-2, global_ymax_combine)

# 11) Volcano function
make_volcano_combine <- function(df, title_text, xlab_text,
                                 xlims = c(-10, 10),
                                 ylims = c(-2, 30),
                                 n_each = 8) {
  df2 <- df %>%
    mutate(P.adj.plot = pmax(P.adj, 1e-300))

  top_genes <- get_top_labels_combine(df2, n_each = n_each)

  ggplot(df2, aes(x = log2.FC., y = -log10(P.adj.plot), color = Significance)) +
    geom_point(alpha = 0.65, size = 1.6) +
    geom_text_repel(
      data = top_genes,
      aes(label = Gene_name),
      size = 3.2,
      max.overlaps = Inf,
      show.legend = FALSE
    ) +
    scale_color_manual(values = c(
      "Up regulated" = "firebrick",
      "Down regulated" = "steelblue",
      "Not Significant" = "grey70"
    )) +
    geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
    coord_cartesian(xlim = xlims, ylim = ylims) +
    labs(
      title = title_text,
      x = xlab_text,
      y = "-log10(P.adj)"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      legend.position = "right",
      plot.title = element_text(size = 11),
      axis.title = element_text(size = 9),
      axis.text = element_text(size = 8)
    )
}

# 12) Full plots
p_wk3_combine <- make_volcano_combine(
  wk3_combine,
  title_text = "3 weeks",
  xlab_text  = "log2FC",
  xlims = global_xlim_combine,
  ylims = global_ylim_combine
)

p_wk4_combine <- make_volcano_combine(
  wk4_combine,
  title_text = "4 weeks",
  xlab_text  = "log2FC",
  xlims = global_xlim_combine,
  ylims = global_ylim_combine
)

p_early_combine <- make_volcano_combine(
  early_combine,
  title_text = "Early onset",
  xlab_text  = "log2FC",
  xlims = global_xlim_combine,
  ylims = global_ylim_combine
)

p_term_combine <- make_volcano_combine(
  term_combine,
  title_text = "Terminal corrected: LSD1 terminal + Trem2WT terminal",
  xlab_text  = "log2FC",
  xlims = global_xlim_combine,
  ylims = global_ylim_combine
)

# 13) MODULE / Microglial and immune genes overlay
#     Keeps the same old logic: yellow module only.
allc_combine <- read.csv(allc_file, header = TRUE, check.names = FALSE)
mod_combine  <- read.csv(mod_file,  header = TRUE, check.names = FALSE)

names(mod_combine)  <- sub("^\\ufeff", "", names(mod_combine))
names(allc_combine) <- sub("^\\ufeff", "", names(allc_combine))

bad_combine <- is.na(names(mod_combine)) | names(mod_combine) == ""
if (any(bad_combine)) {
  names(mod_combine)[bad_combine] <- paste0("V", seq_len(sum(bad_combine)))
}

names(mod_combine) <- make.unique(names(mod_combine))
allc_combine$Gene_Symbol <- toupper(as.character(allc_combine$Gene_Symbol))

target_term_combine <- "Immune functions"

immune_modules_combine <- mod_combine %>%
  filter(trimws(CategoryTerm) == target_term_combine) %>%
  pull(Module) %>%
  unique()

if (length(immune_modules_combine) == 0) {
  stop("No modules found for CategoryTerm == 'Immune functions'.")
}

gene2module_combine <- allc_combine %>%
  select(Gene_Symbol, Module) %>%
  filter(!is.na(Gene_Symbol), !is.na(Module)) %>%
  distinct()

all_genes_universe_combine <- unique(toupper(c(
  early_combine$Gene_name,
  wk3_combine$Gene_name,
  wk4_combine$Gene_name,
  term_combine$Gene_name
)))

setA_yellow_combine <- gene2module_combine %>%
  filter(Module %in% immune_modules_combine) %>%
  filter(tolower(trimws(Module)) == "yellow") %>%
  pull(Gene_Symbol) %>%
  unique()

setA_yellow_combine <- intersect(setA_yellow_combine, all_genes_universe_combine)

if (length(setA_yellow_combine) == 0) {
  stop("Set A (Microglial and immune genes / yellow module) is empty.")
}

setB_no_yellow_combine <- setdiff(all_genes_universe_combine, setA_yellow_combine)

cat("\n[CORRECTED COMBINE TERMINAL] Microglial and immune genes:", length(setA_yellow_combine), "\n")
cat("[CORRECTED COMBINE TERMINAL] Other genes:", length(setB_no_yellow_combine), "\n")
cat("[CORRECTED COMBINE TERMINAL] Overlap check:", length(intersect(setA_yellow_combine, setB_no_yellow_combine)), "\n\n")

filter_to_gene_set_combine <- function(df, gene_set) {
  df %>%
    mutate(Gene_UP = toupper(as.character(Gene_name))) %>%
    filter(Gene_UP %in% gene_set) %>%
    select(-Gene_UP)
}

make_volcano_overlay_yellow_combine <- function(df, title_text, xlab_text,
                                                xlims = c(-10, 10),
                                                ylims = c(-2, 30),
                                                yellow_genes,
                                                n_each = 8,
                                                base_point_size = 1.6,
                                                yellow_point_size = 1.2) {
  df2 <- df %>%
    mutate(
      P.adj.plot = pmax(P.adj, 1e-300),
      Gene_UP = toupper(as.character(Gene_name)),
      Significance = case_when(
        P.adj < 0.05 & log2.FC. > 1  ~ "Up regulated",
        P.adj < 0.05 & log2.FC. < -1 ~ "Down regulated",
        TRUE ~ "Not Significant"
      ),
      is_yellow = Gene_UP %in% toupper(yellow_genes)
    )

  top_genes <- get_top_labels_combine(df2, n_each = n_each)

  ggplot(df2, aes(x = log2.FC., y = -log10(P.adj.plot))) +
    geom_point(aes(color = Significance), alpha = 0.65, size = base_point_size) +
    geom_point(
      data = df2 %>% filter(is_yellow),
      shape = 21,
      fill = "gold",
      color = "black",
      stroke = 0.35,
      size = yellow_point_size,
      alpha = 0.95
    ) +
    geom_text_repel(
      data = top_genes,
      aes(label = Gene_name),
      size = 3.2,
      max.overlaps = Inf,
      show.legend = FALSE
    ) +
    scale_color_manual(values = c(
      "Up regulated" = "firebrick",
      "Down regulated" = "steelblue",
      "Not Significant" = "grey70"
    )) +
    geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
    coord_cartesian(xlim = xlims, ylim = ylims) +
    labs(
      title = title_text,
      x = xlab_text,
      y = "-log10(P.adj)"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      legend.position = "right",
      plot.title = element_text(size = 11),
      axis.title = element_text(size = 9),
      axis.text = element_text(size = 8)
    )
}

early_A_combine <- filter_to_gene_set_combine(early_combine, setA_yellow_combine)
wk3_A_combine   <- filter_to_gene_set_combine(wk3_combine, setA_yellow_combine)
wk4_A_combine   <- filter_to_gene_set_combine(wk4_combine, setA_yellow_combine)
term_A_combine  <- filter_to_gene_set_combine(term_combine, setA_yellow_combine)

early_B_combine <- filter_to_gene_set_combine(early_combine, setB_no_yellow_combine)
wk3_B_combine   <- filter_to_gene_set_combine(wk3_combine, setB_no_yellow_combine)
wk4_B_combine   <- filter_to_gene_set_combine(wk4_combine, setB_no_yellow_combine)
term_B_combine  <- filter_to_gene_set_combine(term_combine, setB_no_yellow_combine)

p_wk3_A_combine <- make_volcano_combine(wk3_A_combine, "3 weeks (ONLY Microglial and immune genes)", "log2FC", global_xlim_combine, global_ylim_combine)
p_wk4_A_combine <- make_volcano_combine(wk4_A_combine, "4 weeks (ONLY Microglial and immune genes)", "log2FC", global_xlim_combine, global_ylim_combine)
p_early_A_combine <- make_volcano_combine(early_A_combine, "Early onset (ONLY Microglial and immune genes)", "log2FC", global_xlim_combine, global_ylim_combine)
p_term_A_combine <- make_volcano_combine(term_A_combine, "Terminal corrected (ONLY Microglial and immune genes)", "log2FC", global_xlim_combine, global_ylim_combine)

p_wk3_B_combine <- make_volcano_combine(wk3_B_combine, "3 weeks (WITHOUT Microglial and immune genes)", "log2FC", global_xlim_combine, global_ylim_combine)
p_wk4_B_combine <- make_volcano_combine(wk4_B_combine, "4 weeks (WITHOUT Microglial and immune genes)", "log2FC", global_xlim_combine, global_ylim_combine)
p_early_B_combine <- make_volcano_combine(early_B_combine, "Early onset (WITHOUT Microglial and immune genes)", "log2FC", global_xlim_combine, global_ylim_combine)
p_term_B_combine <- make_volcano_combine(term_B_combine, "Terminal corrected (WITHOUT Microglial and immune genes)", "log2FC", global_xlim_combine, global_ylim_combine)

p_wk3_overlay_combine <- make_volcano_overlay_yellow_combine(wk3_combine, "3 weeks", "log2FC", global_xlim_combine, global_ylim_combine, setA_yellow_combine)
p_wk4_overlay_combine <- make_volcano_overlay_yellow_combine(wk4_combine, "4 weeks", "log2FC", global_xlim_combine, global_ylim_combine, setA_yellow_combine)
p_early_overlay_combine <- make_volcano_overlay_yellow_combine(early_combine, "Early onset", "log2FC", global_xlim_combine, global_ylim_combine, setA_yellow_combine)
p_term_overlay_combine <- make_volcano_overlay_yellow_combine(term_combine, "Terminal corrected", "log2FC", global_xlim_combine, global_ylim_combine, setA_yellow_combine)

# 14) 4-panel layout
one_page_4panel_combine <- function(p_wk3, p_wk4, p_early, p_term) {
  (p_wk3 | p_wk4 | p_early | p_term) +
    plot_layout(guides = "collect") &
    theme(
      legend.position = "right",
      plot.title = element_text(size = 11),
      axis.title = element_text(size = 9),
      axis.text  = element_text(size = 8),
      legend.title = element_text(size = 9),
      legend.text  = element_text(size = 8)
    )
}

# 15) PDF output
final_pdf_combine <- file.path(
  out_dir,
  "VOLCANO_CORRECTED_combine_terminal_LSD1terminal_plus_Trem2WTterminal_NO_MICHAEL_GALAXY_STYLE_FC.pdf"
)

pdf(final_pdf_combine, width = 23, height = 6.5, onefile = TRUE)

print(one_page_4panel_combine(
  p_wk3_combine,
  p_wk4_combine,
  p_early_combine,
  p_term_combine
))

print(one_page_4panel_combine(
  p_wk3_A_combine,
  p_wk4_A_combine,
  p_early_A_combine,
  p_term_A_combine
))

print(one_page_4panel_combine(
  p_wk3_B_combine,
  p_wk4_B_combine,
  p_early_B_combine,
  p_term_B_combine
))

print(one_page_4panel_combine(
  p_wk3_overlay_combine,
  p_wk4_overlay_combine,
  p_early_overlay_combine,
  p_term_overlay_combine
))

dev.off()

message("Done. Wrote: ", final_pdf_combine)

# 16) Export UP / DOWN genes
get_up_table_combine <- function(df, timepoint_name) {
  df %>%
    filter(P.adj < 0.05, log2.FC. > 1) %>%
    mutate(
      Timepoint = timepoint_name,
      Gene_name = as.character(Gene_name)
    ) %>%
    select(Timepoint, Gene_name, log2.FC., P.adj, Significance, baseMean, P.Value, Model, Contrast) %>%
    arrange(P.adj)
}

get_down_table_combine <- function(df, timepoint_name) {
  df %>%
    filter(P.adj < 0.05, log2.FC. < -1) %>%
    mutate(
      Timepoint = timepoint_name,
      Gene_name = as.character(Gene_name)
    ) %>%
    select(Timepoint, Gene_name, log2.FC., P.adj, Significance, baseMean, P.Value, Model, Contrast) %>%
    arrange(P.adj)
}

up_wk3_tbl_combine    <- get_up_table_combine(wk3_combine, "3 weeks")
up_wk4_tbl_combine    <- get_up_table_combine(wk4_combine, "4 weeks")
up_early_tbl_combine  <- get_up_table_combine(early_combine, "Early onset")
up_term_tbl_combine   <- get_up_table_combine(term_combine, "Terminal corrected")

down_wk3_tbl_combine   <- get_down_table_combine(wk3_combine, "3 weeks")
down_wk4_tbl_combine   <- get_down_table_combine(wk4_combine, "4 weeks")
down_early_tbl_combine <- get_down_table_combine(early_combine, "Early onset")
down_term_tbl_combine  <- get_down_table_combine(term_combine, "Terminal corrected")

up_all_tbl_combine <- bind_rows(
  up_wk3_tbl_combine,
  up_wk4_tbl_combine,
  up_early_tbl_combine,
  up_term_tbl_combine
)

down_all_tbl_combine <- bind_rows(
  down_wk3_tbl_combine,
  down_wk4_tbl_combine,
  down_early_tbl_combine,
  down_term_tbl_combine
)

write.csv(up_wk3_tbl_combine,    file.path(out_dir, "UP_genes_3weeks_CORRECTED_COMBINE_terminal.csv"), row.names = FALSE)
write.csv(up_wk4_tbl_combine,    file.path(out_dir, "UP_genes_4weeks_CORRECTED_COMBINE_terminal.csv"), row.names = FALSE)
write.csv(up_early_tbl_combine,  file.path(out_dir, "UP_genes_Early_onset_CORRECTED_COMBINE_terminal.csv"), row.names = FALSE)
write.csv(up_term_tbl_combine,   file.path(out_dir, "UP_genes_Terminal_CORRECTED_LSD1terminal_plus_Trem2WTterminal_NO_MICHAEL.csv"), row.names = FALSE)
write.csv(up_all_tbl_combine,    file.path(out_dir, "UP_genes_ALL_4_timepoints_CORRECTED_COMBINE_terminal.csv"), row.names = FALSE)

write.csv(down_wk3_tbl_combine,   file.path(out_dir, "DOWN_genes_3weeks_CORRECTED_COMBINE_terminal.csv"), row.names = FALSE)
write.csv(down_wk4_tbl_combine,   file.path(out_dir, "DOWN_genes_4weeks_CORRECTED_COMBINE_terminal.csv"), row.names = FALSE)
write.csv(down_early_tbl_combine, file.path(out_dir, "DOWN_genes_Early_onset_CORRECTED_COMBINE_terminal.csv"), row.names = FALSE)
write.csv(down_term_tbl_combine,  file.path(out_dir, "DOWN_genes_Terminal_CORRECTED_LSD1terminal_plus_Trem2WTterminal_NO_MICHAEL.csv"), row.names = FALSE)
write.csv(down_all_tbl_combine,   file.path(out_dir, "DOWN_genes_ALL_4_timepoints_CORRECTED_COMBINE_terminal.csv"), row.names = FALSE)

# 17) Export cleaned final tables used for volcano
write.csv(wk3_combine,   file.path(out_dir, "CLEANED_3weeks_CORRECTED_COMBINE_terminal.csv"), row.names = FALSE)
write.csv(wk4_combine,   file.path(out_dir, "CLEANED_4weeks_CORRECTED_COMBINE_terminal.csv"), row.names = FALSE)
write.csv(early_combine, file.path(out_dir, "CLEANED_Early_onset_CORRECTED_COMBINE_terminal.csv"), row.names = FALSE)
write.csv(term_combine,  file.path(out_dir, "CLEANED_Terminal_CORRECTED_LSD1terminal_plus_Trem2WTterminal_NO_MICHAEL.csv"), row.names = FALSE)

summary_df <- data.frame(
  Timepoint = c("3 weeks", "4 weeks", "Early onset", "Terminal corrected"),
  Model = c(
    unique(wk3_combine$Model),
    unique(wk4_combine$Model),
    unique(early_combine$Model),
    unique(term_combine$Model)
  ),
  n_UP = c(
    nrow(up_wk3_tbl_combine),
    nrow(up_wk4_tbl_combine),
    nrow(up_early_tbl_combine),
    nrow(up_term_tbl_combine)
  ),
  n_DOWN = c(
    nrow(down_wk3_tbl_combine),
    nrow(down_wk4_tbl_combine),
    nrow(down_early_tbl_combine),
    nrow(down_term_tbl_combine)
  )
)

write.csv(
  summary_df,
  file.path(out_dir, "SUMMARY_CORRECTED_volcano_DESeq2_UP_DOWN_counts.csv"),
  row.names = FALSE
)

cat("\n[CORRECTED COMBINE TERMINAL] Exported UP-regulated genes:\n")
cat("3 weeks:", nrow(up_wk3_tbl_combine), "\n")
cat("4 weeks:", nrow(up_wk4_tbl_combine), "\n")
cat("Early onset:", nrow(up_early_tbl_combine), "\n")
cat("Terminal corrected:", nrow(up_term_tbl_combine), "\n")
cat("All combined:", nrow(up_all_tbl_combine), "\n")

cat("\n[CORRECTED COMBINE TERMINAL] Exported DOWN-regulated genes:\n")
cat("3 weeks:", nrow(down_wk3_tbl_combine), "\n")
cat("4 weeks:", nrow(down_wk4_tbl_combine), "\n")
cat("Early onset:", nrow(down_early_tbl_combine), "\n")
cat("Terminal corrected:", nrow(down_term_tbl_combine), "\n")
cat("All combined:", nrow(down_all_tbl_combine), "\n")

message("Done. Wrote corrected volcano PDF, DE tables, UP gene tables, DOWN gene tables, and check files.")
