# 08_plx_terminal_microglia_immune_heatmap_by_sex.R
#
# Purpose: PLX terminal sex-stratified microglia and immune-response heatmap using Trem2 gene order.
# Inputs are expected under data/processed/ or data/external/ relative to this repository.
# Outputs are written under results/ or script-defined subfolders.

# PLX terminal microglia / immune response heatmap - WITHIN-SEX NORMALIZED
# 32 single-sample featureCounts files: 1301-1332
# Genotype logic:
#   cre- = LSD1WT
#   cre+ = LSD1KO
#   ctrl/plx = treatment annotation
# No ComBat / no batch correction / no top-variable-gene filtering
#
# Gene-category logic follows the Trem2 immune heatmap workflow:
#   modules.csv: CategoryTerm == "Immune functions" -> target modules
#   NIHMS472534-supplement-02.csv: Gene_Symbol -> Module
#   Yellow module display title = "Microglial and immune genes"
#
# Output:
#   ONE PDF with pages ordered by module:
#     Yellow group mean: Female left + Male right
#     Yellow individual samples: Female left + Male right
#     Light cyan group mean: Female left + Male right
#     Light cyan individual samples: Female left + Male right
#     remaining Immune-functions modules in the same side-by-side format
#   Z-score normalization is calculated WITHIN SEX only.
#   Group mean is calculated FROM within-sex individual sample z-scores.
#
# NEW requested version:
#   Use the exact Trem2 heatmap gene order from:
#     ALL_HEATMAP_GENE_ORDER_with_LeftCleanRowIndex.csv
#   Add the same left-side clean 0-based row index axis.
#   PLX will NOT re-cluster rows for the immune module heatmaps.


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
  library(readr)
})

katz_root <- repo_root

project_dir <- file.path(katz_root, "PLX terminal experiment/results")
featurecounts_dir <- file.path(project_dir, "featurecounts")
genotype_file <- file.path(project_dir, "genotypes.xlsx")

# Use the same immune-module annotation files used in previous heatmap workflows.
# These paths were present in the file inventory.
allc_file <- file.path(katz_root, "PLX RNA sequencing/NIHMS472534-supplement-02.csv")
mod_file  <- file.path(katz_root, "PLX RNA sequencing/modules.csv")

# Exact row order exported from the Trem2 immune heatmap script.
# This is the file that makes PLX use the SAME gene order as Trem2.
trem2_gene_order_file <- file.path(
  katz_root,
  "Trem2 terminal experiment/ALL_HEATMAP_GENE_ORDER_with_LeftCleanRowIndex.csv"
)

out_dir <- file.path(project_dir, "heatmap_microglia_immune_response_by_sex_Trem2_gene_order_SIDE_BY_SIDE_WITHIN_SEX_NORMALIZED")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# 1) Basic file checks
if (!file.exists(genotype_file)) stop("Missing genotype file: ", genotype_file)
if (!dir.exists(featurecounts_dir)) stop("Missing featureCounts folder: ", featurecounts_dir)
if (!file.exists(allc_file)) stop("Missing allc/module gene file: ", allc_file)
if (!file.exists(mod_file)) stop("Missing modules file: ", mod_file)
if (!file.exists(trem2_gene_order_file)) {
  stop(
    "Missing Trem2 gene-order file: ", trem2_gene_order_file, "\n",
    "Run the Trem2 heatmap script first, or move ALL_HEATMAP_GENE_ORDER_with_LeftCleanRowIndex.csv to this path."
  )
}

sample_ids <- as.character(1301:1332)
featurecount_files <- file.path(featurecounts_dir, paste0(sample_ids, " featurecounts.xlsx"))
missing_fc <- featurecount_files[!file.exists(featurecount_files)]
if (length(missing_fc) > 0) {
  stop("Missing featureCounts files:\n", paste(missing_fc, collapse = "\n"))
}

fix_names <- function(df) {
  names(df) <- sub("^\\ufeff", "", names(df))
  bad <- is.na(names(df)) | names(df) == ""
  if (any(bad)) names(df)[bad] <- paste0("V", seq_len(sum(bad)))
  names(df) <- make.unique(names(df))
  df
}

geno_raw <- readxl::read_excel(genotype_file, sheet = 1, col_names = TRUE)
geno_raw <- fix_names(geno_raw)

# Expected first columns from current file:
#   well/animal label, full sample ID, short sample number, genotype, Sex
# Make this robust by taking the first 5 columns.
if (ncol(geno_raw) < 5) stop("genotypes.xlsx has fewer than 5 columns; check file format.")

geno <- geno_raw[, 1:5]
names(geno) <- c("Animal_or_well", "Full_ID", "Sample", "Raw_genotype", "Sex")

geno <- geno %>%
  dplyr::mutate(
    Sample = as.character(Sample),
    Raw_genotype = trimws(as.character(Raw_genotype)),
    Sex = trimws(as.character(Sex)),
    Treatment = dplyr::case_when(
      grepl("^ctrl", Raw_genotype, ignore.case = TRUE) ~ "Ctrl",
      grepl("^plx",  Raw_genotype, ignore.case = TRUE) ~ "PLX",
      TRUE ~ NA_character_
    ),
    LSD1_status = dplyr::case_when(
      grepl("cre\\-", Raw_genotype, ignore.case = TRUE) ~ "LSD1WT",
      grepl("cre\\+", Raw_genotype, ignore.case = TRUE) ~ "LSD1KO",
      TRUE ~ NA_character_
    ),
    Group4 = paste(Treatment, LSD1_status, sep = "_")
  ) %>%
  dplyr::filter(Sample %in% sample_ids)

if (nrow(geno) != length(sample_ids)) {
  stop("genotype file did not map all 1301-1332 samples. Mapped n = ", nrow(geno))
}
if (any(is.na(geno$Treatment)) || any(is.na(geno$LSD1_status))) {
  print(geno %>% dplyr::filter(is.na(Treatment) | is.na(LSD1_status)))
  stop("Some genotype rows could not be mapped to Treatment or LSD1_status.")
}

# fixed plot order: control first, then PLX; within each, cre-/WT first, then cre+/KO
group_order <- c("Ctrl_LSD1WT", "Ctrl_LSD1KO", "PLX_LSD1WT", "PLX_LSD1KO")
geno$Group4 <- factor(geno$Group4, levels = group_order)
geno <- geno %>% dplyr::arrange(Group4, Sex, Sample)

cat("\n=== Sample counts by 4 groups ===\n")
print(table(geno$Group4, useNA = "ifany"))
cat("\n=== Sample counts by 4 groups and sex ===\n")
print(table(geno$Group4, geno$Sex, useNA = "ifany"))

readr::write_csv(geno, file.path(out_dir, "PLX_terminal_heatmap_sample_info.csv"))

read_one_featurecounts <- function(file, sample_id) {
  df <- readxl::read_excel(file, sheet = 1, col_names = FALSE)
  df <- as.data.frame(df, stringsAsFactors = FALSE)

  # Keep first 8 columns; expected:
  # Geneid, Count, Chr, Start, End, Strand, Extra/Length, Gene_name
  if (ncol(df) < 8) {
    stop("FeatureCounts file has fewer than 8 columns: ", file)
  }
  df <- df[, 1:8, drop = FALSE]
  names(df) <- c("Geneid", sample_id, "Chr", "Start", "End", "Strand", "Extra", "Gene_name")

  # Remove fake/header-like first row if present.
  first_gene <- tolower(trimws(as.character(df$Geneid[1])))
  first_count <- suppressWarnings(as.numeric(df[[sample_id]][1]))
  if (first_gene %in% c("geneid", "gene_id", "gene id") || is.na(first_count)) {
    df <- df[-1, , drop = FALSE]
  }

  df <- df %>%
    dplyr::mutate(
      Geneid = trimws(as.character(Geneid)),
      Gene_name = toupper(trimws(as.character(Gene_name))),
      !!sample_id := suppressWarnings(as.numeric(.data[[sample_id]]))
    ) %>%
    dplyr::filter(!is.na(Geneid), Geneid != "")

  # Save cleaned single-sample version with correct column names.
  readr::write_csv(df, file.path(out_dir, paste0("CLEANED_featureCounts_", sample_id, ".csv")))

  df %>% dplyr::select(Geneid, Gene_name, dplyr::all_of(sample_id))
}

fc_list <- Map(read_one_featurecounts, featurecount_files, sample_ids)

# Merge all samples by Geneid + Gene_name.
merged_geneid <- Reduce(function(x, y) {
  dplyr::full_join(x, y, by = c("Geneid", "Gene_name"))
}, fc_list)

for (sid in sample_ids) {
  merged_geneid[[sid]][is.na(merged_geneid[[sid]])] <- 0
}

# If the same Geneid has multiple Gene_name labels, keep all rows at geneid level,
# then collapse by symbol below for heatmap.
merged_geneid <- merged_geneid %>% dplyr::arrange(Geneid)
readr::write_csv(merged_geneid, file.path(out_dir, "PLX_terminal_raw_featureCounts_geneid_1301_1332.csv"))

# Collapse to gene-symbol level because immune module annotation uses Gene_Symbol.
merged_symbol <- merged_geneid %>%
  dplyr::filter(!is.na(Gene_name), Gene_name != "", Gene_name != "NA") %>%
  dplyr::group_by(Gene_name) %>%
  dplyr::summarise(
    dplyr::across(dplyr::all_of(sample_ids), ~ sum(as.numeric(.x), na.rm = TRUE)),
    .groups = "drop"
  ) %>%
  dplyr::arrange(Gene_name)

readr::write_csv(merged_symbol, file.path(out_dir, "PLX_terminal_raw_featureCounts_symbol_1301_1332.csv"))

count_mat <- as.matrix(merged_symbol[, sample_ids, drop = FALSE])
rownames(count_mat) <- merged_symbol$Gene_name
storage.mode(count_mat) <- "numeric"

# Reorder columns by group order.
sample_order <- geno$Sample
count_mat <- count_mat[, sample_order, drop = FALSE]

allc <- read.csv(allc_file, header = TRUE, check.names = FALSE)
mod  <- read.csv(mod_file,  header = TRUE, check.names = FALSE)
allc <- fix_names(allc)
mod  <- fix_names(mod)

if (!("Gene_Symbol" %in% names(allc))) stop("Gene_Symbol column not found in: ", allc_file)
if (!("Module" %in% names(allc))) stop("Module column not found in: ", allc_file)
if (!("CategoryTerm" %in% names(mod))) stop("CategoryTerm column not found in: ", mod_file)
if (!("Module" %in% names(mod))) stop("Module column not found in: ", mod_file)

allc$Gene_Symbol <- toupper(trimws(as.character(allc$Gene_Symbol)))
mod$CategoryTerm <- trimws(as.character(mod$CategoryTerm))
mod$Module <- trimws(as.character(mod$Module))

# Find modules classified as immune-function modules.
target_term <- "Immune functions"
target_modules <- mod %>%
  dplyr::filter(CategoryTerm == target_term) %>%
  dplyr::pull(Module) %>%
  unique()

if (length(target_modules) == 0) {
  stop("No modules found for CategoryTerm == 'Immune functions'. Check modules.csv.")
}

cat("\n=== Immune-function modules ===\n")
print(target_modules)

gene2module <- allc %>%
  dplyr::select(Gene_Symbol, Module) %>%
  dplyr::filter(!is.na(Gene_Symbol), Gene_Symbol != "", !is.na(Module), Module != "") %>%
  dplyr::distinct()

# Save matched gene lists for checking.
immune_gene_summary <- gene2module %>%
  dplyr::filter(Module %in% target_modules) %>%
  dplyr::mutate(Detected_in_featureCounts = Gene_Symbol %in% rownames(count_mat))
readr::write_csv(immune_gene_summary, file.path(out_dir, "PLX_terminal_immune_module_gene_mapping_check.csv"))

# 6) Expression transforms
log_mat <- log2(count_mat + 1)
readr::write_csv(
  as.data.frame(log_mat) %>% tibble::rownames_to_column("Gene_name"),
  file.path(out_dir, "PLX_terminal_log2_count_plus1_symbol_1301_1332.csv")
)

row_zscore <- function(mat) {
  z <- t(scale(t(mat), center = TRUE, scale = TRUE))
  z[!is.finite(z)] <- 0
  z
}

sex_display_name <- function(sex_value) {
  dplyr::case_when(
    sex_value == "F" ~ "Female",
    sex_value == "M" ~ "Male",
    TRUE ~ as.character(sex_value)
  )
}

# IMPORTANT NORMALIZATION LOGIC:
#   z-score is calculated WITHIN SEX.
#
#   Female samples:
#     each gene is centered/scaled only across Female samples.
#
#   Male samples:
#     each gene is centered/scaled only across Male samples.
#
#   Female and Male are NOT normalized together.
#   This prevents sex baseline expression differences from affecting
#   within-sex heatmap patterns.
make_sex_specific_z_mat <- function(log_mat, sample_info, sex_col = "Sex", sample_col = "Sample") {
  z_list <- list()
  z_check_list <- list()

  for (sx in c("F", "M")) {
    sx_samples <- sample_info %>%
      dplyr::filter(.data[[sex_col]] == sx) %>%
      dplyr::pull(.data[[sample_col]]) %>%
      as.character()

    sx_samples <- intersect(colnames(log_mat), sx_samples)

    if (length(sx_samples) == 0) {
      stop("No samples found in log_mat for sex = ", sx)
    }

    log_mat_sx <- log_mat[, sx_samples, drop = FALSE]
    z_mat_sx <- row_zscore(log_mat_sx)

    z_list[[sx]] <- z_mat_sx

    sx_name <- sex_display_name(sx)

    readr::write_csv(
      as.data.frame(z_mat_sx) %>% tibble::rownames_to_column("Gene_name"),
      file.path(out_dir, paste0("PLX_terminal_row_zscore_WITHIN_SEX_", sx_name, ".csv"))
    )

    z_check_list[[sx]] <- data.frame(
      Sex = sx_name,
      N_samples_used_for_zscore = length(sx_samples),
      Samples_used_for_zscore = paste(sx_samples, collapse = ";")
    )
  }

  z_check_df <- dplyr::bind_rows(z_check_list)
  readr::write_csv(
    z_check_df,
    file.path(out_dir, "CHECK_WITHIN_SEX_zscore_sample_sets.csv")
  )

  # Combine Female and Male z-score matrices back together.
  # Columns are then reordered to match sample_order used by the heatmap.
  z_combined <- cbind(z_list[["F"]], z_list[["M"]])
  z_combined <- z_combined[, sample_order, drop = FALSE]

  return(z_combined)
}

z_mat <- make_sex_specific_z_mat(
  log_mat = log_mat,
  sample_info = geno,
  sex_col = "Sex",
  sample_col = "Sample"
)

readr::write_csv(
  as.data.frame(z_mat) %>% tibble::rownames_to_column("Gene_name"),
  file.path(out_dir, "PLX_terminal_row_zscore_WITHIN_SEX_symbol_1301_1332.csv")
)

cat("\n=== Z-score normalization logic ===\n")
cat("PASS: z-score was calculated separately within Female and Male samples.\n")
cat("Female and Male samples were NOT normalized together.\n")

group_colors <- c(
  "LSD1WT" = "steelblue",
  "LSD1KO" = "firebrick"
)

plx_colors <- c(
  "Ctrl" = "grey55",
  "PLX"  = "darkorange2"
)

sex_colors <- c(
  "F" = "orchid3",
  "M" = "seagreen4"
)

col_fun <- circlize::colorRamp2(
  c(-4, -2, 0, 2, 4),
  c("blue", "lightblue", "white", "salmon", "red")
)

make_group_mean_from_individual_z_by_sex <- function(ind_z, sample_info, group_order, sex_value) {
  stopifnot(all(colnames(ind_z) %in% sample_info$Sample))
  info <- sample_info[match(colnames(ind_z), sample_info$Sample), , drop = FALSE]

  keep <- info$Sex == sex_value
  ind_z_sex <- ind_z[, keep, drop = FALSE]
  info_sex <- info[keep, , drop = FALSE]

  if (ncol(ind_z_sex) == 0) {
    stop("No samples found for sex = ", sex_value)
  }

  gid <- as.character(info_sex$Group4)

  out <- sapply(group_order, function(g) {
    idx <- which(gid == g)
    if (length(idx) == 0) return(rep(NA_real_, nrow(ind_z_sex)))
    if (length(idx) == 1) return(ind_z_sex[, idx])
    rowMeans(ind_z_sex[, idx, drop = FALSE], na.rm = TRUE)
  })
  out <- as.matrix(out)
  rownames(out) <- rownames(ind_z_sex)
  colnames(out) <- group_order

  list(mean_z = out, ind_z = ind_z_sex, sample_info = info_sex)
}

make_group_labels <- function(group_order) {
  group_order
}

make_individual_labels_by_sex <- function(sample_ids_ordered, sample_info_sex) {
  info <- sample_info_sex[match(sample_ids_ordered, sample_info_sex$Sample), , drop = FALSE]
  labels <- character(nrow(info))
  for (g in group_order) {
    idx <- which(as.character(info$Group4) == g)
    if (length(idx) > 0) {
      labels[idx] <- paste0(g, "_", sprintf("%02d", seq_along(idx)))
    }
  }
  labels
}

make_group_anno <- function(cols_group) {
  s <- as.character(cols_group)
  treatment <- ifelse(grepl("^Ctrl_", s), "Ctrl", ifelse(grepl("^PLX_", s), "PLX", NA_character_))
  lsd1 <- ifelse(grepl("_LSD1WT$", s), "LSD1WT", ifelse(grepl("_LSD1KO$", s), "LSD1KO", NA_character_))

  HeatmapAnnotation(
    Treatment = factor(treatment, levels = c("Ctrl", "PLX")),
    Group = factor(lsd1, levels = c("LSD1WT", "LSD1KO")),
    col = list(
      Treatment = plx_colors,
      Group = group_colors
    ),
    na_col = "grey90",
    show_annotation_name = TRUE,
    annotation_name_gp = grid::gpar(fontsize = 10)
  )
}

make_individual_anno_by_sex <- function(sample_ids_ordered, sample_info_sex) {
  info <- sample_info_sex[match(sample_ids_ordered, sample_info_sex$Sample), , drop = FALSE]

  HeatmapAnnotation(
    Treatment = factor(info$Treatment, levels = c("Ctrl", "PLX")),
    Group = factor(info$LSD1_status, levels = c("LSD1WT", "LSD1KO")),
    col = list(
      Treatment = plx_colors,
      Group = group_colors
    ),
    na_col = "grey90",
    show_annotation_name = TRUE,
    annotation_name_gp = grid::gpar(fontsize = 10)
  )
}

make_module_title <- function(module_name, sex_value, page_type = c("group mean", "individual samples")) {
  page_type <- match.arg(page_type)
  main_title <- if (module_name == "Yellow") {
    "Microglial and immune genes"
  } else {
    paste0("PLX terminal - Immune module: ", module_name)
  }
  paste0(main_title, " - ", sex_display_name(sex_value), " (", page_type, ")")
}

# Page order like old Trem2 code: Yellow first, Light cyan second, then remaining immune modules.
priority_modules <- c("Yellow", "Light cyan")
modules_to_plot <- c(intersect(priority_modules, target_modules), setdiff(target_modules, priority_modules))

cat("\n=== Modules to plot ===\n")
print(modules_to_plot)

# 8) Trem2 gene-order + left clean row-index helpers
trem2_gene_order <- readr::read_csv(trem2_gene_order_file, show_col_types = FALSE) %>%
  fix_names()

required_order_cols <- c("Module", "Row_index_0_based", "Row_order_1_based", "Gene")
missing_order_cols <- setdiff(required_order_cols, names(trem2_gene_order))
if (length(missing_order_cols) > 0) {
  stop(
    "Trem2 gene-order file is missing required columns: ",
    paste(missing_order_cols, collapse = ", ")
  )
}

trem2_gene_order <- trem2_gene_order %>%
  dplyr::mutate(
    Module = trimws(as.character(Module)),
    Gene = toupper(trimws(as.character(Gene))),
    Row_index_0_based = as.integer(Row_index_0_based),
    Row_order_1_based = as.integer(Row_order_1_based)
  ) %>%
  dplyr::filter(!is.na(Module), Module != "", !is.na(Gene), Gene != "")

get_trem2_row_order_for_module <- function(module_name, available_genes, strict_missing = TRUE) {
  ord_df <- trem2_gene_order %>%
    dplyr::filter(Module == module_name) %>%
    dplyr::arrange(Row_index_0_based, Row_order_1_based)

  if (nrow(ord_df) == 0) {
    stop("No Trem2 gene order found for module: ", module_name)
  }

  row_ord_full <- unique(ord_df$Gene)
  available_genes <- toupper(as.character(available_genes))

  missing_from_plx <- setdiff(row_ord_full, available_genes)
  extra_in_plx <- setdiff(available_genes, row_ord_full)

  if (length(missing_from_plx) > 0 && strict_missing) {
    stop(
      "PLX is missing genes that exist in the Trem2 row-order file for module ", module_name, ".\n",
      "Because you requested the SAME gene order, I am stopping instead of silently dropping rows.\n",
      "Missing genes: ", paste(head(missing_from_plx, 50), collapse = ", "),
      ifelse(length(missing_from_plx) > 50, " ...", "")
    )
  }

  row_ord <- intersect(row_ord_full, available_genes)

  list(
    row_ord = row_ord,
    trem2_n = length(row_ord_full),
    plx_n_after_ordering = length(row_ord),
    missing_from_plx = missing_from_plx,
    extra_in_plx_not_used = extra_in_plx
  )
}

choose_clean_tick_step <- function(n) {
  if (n <= 60) return(10)
  if (n <= 500) return(50)
  return(100)
}

get_clean_axis_ticks <- function(n) {
  tick_step <- choose_clean_tick_step(n)
  max_clean_tick <- floor((n - 1) / tick_step) * tick_step
  ticks <- seq(0, max_clean_tick, by = tick_step)
  unique(ticks)
}

make_row_index_axis <- function(n) {
  row_index_0_based <- seq_len(n) - 1
  ticks <- get_clean_axis_ticks(n)

  axis_labels <- rep("", n)
  axis_labels[row_index_0_based %in% ticks] <- as.character(row_index_0_based[row_index_0_based %in% ticks])

  rowAnnotation(
    `Row index` = anno_text(
      axis_labels,
      just = "right",
      location = 0.5,
      gp = grid::gpar(fontsize = 6)
    ),
    annotation_name_side = "top",
    annotation_name_gp = grid::gpar(fontsize = 8),
    width = grid::unit(11, "mm")
  )
}

# 9) Main immune heatmap PDF
#    Side-by-side layout:
#      For each module:
#        Page 1 = Female group mean  |  Male group mean
#        Page 2 = Female individual  |  Male individual
#      Female and Male are still separated.
#      They are only arranged side-by-side on the same PDF page.
#      Both sides use the exact SAME Trem2 row order.
#      The clean row-index axis is shown once, on the far left.
out_pdf <- file.path(
  out_dir,
  "PLX_terminal_Microglia_Immune_Response_Heatmap_BY_SEX_SIDE_BY_SIDE_Trem2GeneOrder_LeftCleanRowIndex_WITHIN_SEX_NORMALIZED.pdf"
)

# Store the exact row order actually used in the PLX heatmaps.
plx_gene_orders_used <- list()
plx_gene_order_checks <- list()

sex_order <- c("F", "M")

# Helper: prepare one sex-specific matrix set for one module using a fixed Trem2 row order.
prepare_plx_module_by_sex <- function(ind_z_all, row_ord, sex_value) {
  sex_result <- make_group_mean_from_individual_z_by_sex(
    ind_z = ind_z_all,
    sample_info = geno,
    group_order = group_order,
    sex_value = sex_value
  )

  mean_z <- sex_result$mean_z
  ind_z  <- sex_result$ind_z
  info_sex <- sex_result$sample_info

  # Order same-sex individual columns by 4 groups, then sample number.
  info_sex <- info_sex %>%
    dplyr::mutate(
      Group4 = factor(Group4, levels = group_order),
      Sample_numeric = as.numeric(Sample)
    ) %>%
    dplyr::arrange(Group4, Sample_numeric)

  sex_sample_order <- info_sex$Sample
  ind_z <- ind_z[, sex_sample_order, drop = FALSE]

  # Force the SAME Trem2 row order.
  mean_z <- mean_z[row_ord, group_order, drop = FALSE]
  ind_z  <- ind_z[row_ord, sex_sample_order, drop = FALSE]

  list(
    mean_z = mean_z,
    ind_z = ind_z,
    info_sex = info_sex,
    group_labels = make_group_labels(colnames(mean_z)),
    ind_labels = make_individual_labels_by_sex(colnames(ind_z), info_sex)
  )
}

# Female annotation: carries the visible annotation legends.
make_group_anno_female <- function(cols_group) {
  make_group_anno(cols_group)
}

make_individual_anno_female <- function(sample_ids_ordered, sample_info_sex) {
  make_individual_anno_by_sex(sample_ids_ordered, sample_info_sex)
}

# Male annotation: unique internal names so ComplexHeatmap does not complain
# when Female and Male heatmaps are combined into one HeatmapList.
# Legends are hidden here because Female already provides the same legends.
make_group_anno_male <- function(cols_group) {
  s <- as.character(cols_group)
  treatment <- ifelse(grepl("^Ctrl_", s), "Ctrl", ifelse(grepl("^PLX_", s), "PLX", NA_character_))
  lsd1 <- ifelse(grepl("_LSD1WT$", s), "LSD1WT", ifelse(grepl("_LSD1KO$", s), "LSD1KO", NA_character_))

  HeatmapAnnotation(
    Treatment_Male = factor(treatment, levels = c("Ctrl", "PLX")),
    Group_Male = factor(lsd1, levels = c("LSD1WT", "LSD1KO")),
    col = list(
      Treatment_Male = plx_colors,
      Group_Male = group_colors
    ),
    na_col = "grey90",
    show_annotation_name = FALSE,
    show_legend = FALSE
  )
}

make_individual_anno_male <- function(sample_ids_ordered, sample_info_sex) {
  info <- sample_info_sex[match(sample_ids_ordered, sample_info_sex$Sample), , drop = FALSE]

  HeatmapAnnotation(
    Treatment_Male = factor(info$Treatment, levels = c("Ctrl", "PLX")),
    Group_Male = factor(info$LSD1_status, levels = c("LSD1WT", "LSD1KO")),
    col = list(
      Treatment_Male = plx_colors,
      Group_Male = group_colors
    ),
    na_col = "grey90",
    show_annotation_name = FALSE,
    show_legend = FALSE
  )
}

pdf(out_pdf, width = 16, height = 8)

for (m in modules_to_plot) {
  message("-> plotting module side-by-side: ", m)

  genes_m <- gene2module %>%
    dplyr::filter(Module == m) %>%
    dplyr::pull(Gene_Symbol) %>%
    unique()

  keep_genes <- intersect(rownames(z_mat), genes_m)
  if (length(keep_genes) == 0) {
    message("Skipped module with no detected genes: ", m)
    next
  }

  ind_z_all <- z_mat[keep_genes, sample_order, drop = FALSE]

  # Use the exact Trem2 heatmap row order for this module.
  # Do NOT re-cluster PLX rows here.
  trem2_order_result <- get_trem2_row_order_for_module(
    module_name = m,
    available_genes = rownames(ind_z_all),
    strict_missing = TRUE
  )
  row_ord <- trem2_order_result$row_ord

  # Clean row-axis ticks for this module.
  axis_ticks <- get_clean_axis_ticks(length(row_ord))
  axis_tick_step <- choose_clean_tick_step(length(row_ord))
  axis_last_labeled_index <- max(axis_ticks)

  # Prepare Female and Male matrices separately, but with identical row order.
  female <- prepare_plx_module_by_sex(ind_z_all, row_ord, "F")
  male   <- prepare_plx_module_by_sex(ind_z_all, row_ord, "M")

  # Store/check row order once per sex, so the CSV proves Female and Male are identical.
  for (sx in sex_order) {
    plx_gene_orders_used[[paste(m, sx, sep = "__")]] <- data.frame(
      Module = m,
      Sex = sex_display_name(sx),
      Total_genes_in_module = length(row_ord),
      Row_index_0_based = seq_along(row_ord) - 1,
      Row_order_1_based = seq_along(row_ord),
      Axis_tick_step = axis_tick_step,
      Axis_last_labeled_index = axis_last_labeled_index,
      Gene = row_ord
    )

    plx_gene_order_checks[[paste(m, sx, sep = "__")]] <- data.frame(
      Module = m,
      Sex = sex_display_name(sx),
      Trem2_order_gene_count = trem2_order_result$trem2_n,
      PLX_ordered_gene_count = trem2_order_result$plx_n_after_ordering,
      Missing_from_PLX_count = length(trem2_order_result$missing_from_plx),
      Extra_PLX_module_genes_not_used_count = length(trem2_order_result$extra_in_plx_not_used),
      Extra_PLX_module_genes_not_used = paste(trem2_order_result$extra_in_plx_not_used, collapse = ";")
    )
  }

  # Page 1: group mean, Female left + Male right
  row_axis_mean <- make_row_index_axis(length(row_ord))

  ht_mean_female <- Heatmap(
    female$mean_z,
    name = "z-score",
    col = col_fun,
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    show_row_dend = FALSE,
    show_column_dend = FALSE,
    column_order = colnames(female$mean_z),
    column_labels = female$group_labels,
    top_annotation = make_group_anno_female(colnames(female$mean_z)),
    left_annotation = row_axis_mean,
    show_row_names = FALSE,
    show_column_names = TRUE,
    column_names_gp = grid::gpar(fontsize = 10),
    column_title = paste0(sex_display_name("F"), " (group mean)"),
    column_title_gp = grid::gpar(fontsize = 12, fontface = "bold"),
    heatmap_legend_param = list(
      title = "z-score",
      title_gp = grid::gpar(fontsize = 10),
      labels_gp = grid::gpar(fontsize = 9)
    )
  )

  ht_mean_male <- Heatmap(
    male$mean_z,
    name = "z-score_male_group_mean",
    col = col_fun,
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    show_row_dend = FALSE,
    show_column_dend = FALSE,
    column_order = colnames(male$mean_z),
    column_labels = male$group_labels,
    top_annotation = make_group_anno_male(colnames(male$mean_z)),
    show_row_names = FALSE,
    show_column_names = TRUE,
    column_names_gp = grid::gpar(fontsize = 10),
    column_title = paste0(sex_display_name("M"), " (group mean)"),
    column_title_gp = grid::gpar(fontsize = 12, fontface = "bold"),
    show_heatmap_legend = FALSE
  )

  group_mean_page_title <- make_module_title(m, "F", "group mean")
  group_mean_page_title <- sub(" - Female \\(group mean\\)$", " - Female vs Male (group mean)", group_mean_page_title)

  draw(
    ht_mean_female + ht_mean_male,
    newpage = TRUE,
    heatmap_legend_side = "right",
    annotation_legend_side = "right",
    merge_legends = TRUE,
    column_title = group_mean_page_title,
    column_title_gp = grid::gpar(fontsize = 15, fontface = "bold")
  )

  # Page 2: individual samples, Female left + Male right
  row_axis_ind <- make_row_index_axis(length(row_ord))

  ht_ind_female <- Heatmap(
    female$ind_z,
    name = "z-score_individual",
    col = col_fun,
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    show_row_dend = FALSE,
    show_column_dend = FALSE,
    column_order = colnames(female$ind_z),
    column_labels = female$ind_labels,
    top_annotation = make_individual_anno_female(colnames(female$ind_z), female$info_sex),
    left_annotation = row_axis_ind,
    show_row_names = FALSE,
    show_column_names = TRUE,
    column_names_gp = grid::gpar(fontsize = 6),
    column_title = paste0(sex_display_name("F"), " (individual samples)"),
    column_title_gp = grid::gpar(fontsize = 12, fontface = "bold"),
    heatmap_legend_param = list(
      title = "z-score",
      title_gp = grid::gpar(fontsize = 10),
      labels_gp = grid::gpar(fontsize = 9)
    )
  )

  ht_ind_male <- Heatmap(
    male$ind_z,
    name = "z-score_male_individual",
    col = col_fun,
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    show_row_dend = FALSE,
    show_column_dend = FALSE,
    column_order = colnames(male$ind_z),
    column_labels = male$ind_labels,
    top_annotation = make_individual_anno_male(colnames(male$ind_z), male$info_sex),
    show_row_names = FALSE,
    show_column_names = TRUE,
    column_names_gp = grid::gpar(fontsize = 6),
    column_title = paste0(sex_display_name("M"), " (individual samples)"),
    column_title_gp = grid::gpar(fontsize = 12, fontface = "bold"),
    show_heatmap_legend = FALSE
  )

  individual_page_title <- make_module_title(m, "F", "individual samples")
  individual_page_title <- sub(" - Female \\(individual samples\\)$", " - Female vs Male (individual samples)", individual_page_title)

  draw(
    ht_ind_female + ht_ind_male,
    newpage = TRUE,
    heatmap_legend_side = "right",
    annotation_legend_side = "right",
    merge_legends = TRUE,
    column_title = individual_page_title,
    column_title_gp = grid::gpar(fontsize = 15, fontface = "bold")
  )
}

dev.off()

plx_gene_order_used_df <- dplyr::bind_rows(plx_gene_orders_used)
plx_gene_order_check_df <- dplyr::bind_rows(plx_gene_order_checks)

readr::write_csv(
  plx_gene_order_used_df,
  file.path(out_dir, "PLX_terminal_gene_order_USED_from_Trem2_SIDE_BY_SIDE_LeftCleanRowIndex_WITHIN_SEX_NORMALIZED.csv")
)

readr::write_csv(
  plx_gene_order_check_df,
  file.path(out_dir, "CHECK_PLX_vs_Trem2_gene_order_by_module_and_sex_SIDE_BY_SIDE_WITHIN_SEX_NORMALIZED.csv")
)

cat("\nSaved main side-by-side heatmap PDF:\n", out_pdf, "\n")
cat("Saved PLX gene order used file:\n", file.path(out_dir, "PLX_terminal_gene_order_USED_from_Trem2_SIDE_BY_SIDE_LeftCleanRowIndex_WITHIN_SEX_NORMALIZED.csv"), "\n")
cat("Saved PLX vs Trem2 gene order check file:\n", file.path(out_dir, "CHECK_PLX_vs_Trem2_gene_order_by_module_and_sex_SIDE_BY_SIDE_WITHIN_SEX_NORMALIZED.csv"), "\n")

# 9) Extra outputs for checking / reuse
# Sex-specific matrices for all detected immune genes.
# Group means are calculated from individual z-scores WITHIN sex only.
all_immune_genes <- gene2module %>%
  dplyr::filter(Module %in% target_modules) %>%
  dplyr::pull(Gene_Symbol) %>%
  unique()
all_immune_keep <- intersect(rownames(z_mat), all_immune_genes)
immune_ind_z_all <- z_mat[all_immune_keep, sample_order, drop = FALSE]

for (sx in sex_order) {
  sex_result <- make_group_mean_from_individual_z_by_sex(
    ind_z = immune_ind_z_all,
    sample_info = geno,
    group_order = group_order,
    sex_value = sx
  )

  info_sex <- sex_result$sample_info %>%
    dplyr::mutate(Group4 = factor(Group4, levels = group_order)) %>%
    dplyr::arrange(Group4, Sample)

  sex_sample_order <- info_sex$Sample
  immune_ind_z_sex <- sex_result$ind_z[, sex_sample_order, drop = FALSE]
  immune_mean_z_sex <- sex_result$mean_z[, group_order, drop = FALSE]

  sx_name <- sex_display_name(sx)

  readr::write_csv(
    as.data.frame(immune_ind_z_sex) %>% tibble::rownames_to_column("Gene_name"),
    file.path(out_dir, paste0("PLX_terminal_immune_individual_row_zscore_matrix_", sx_name, ".csv"))
  )

  readr::write_csv(
    as.data.frame(immune_mean_z_sex) %>% tibble::rownames_to_column("Gene_name"),
    file.path(out_dir, paste0("PLX_terminal_immune_group_mean_from_individual_zscore_matrix_", sx_name, ".csv"))
  )
}

cat("Saved output folder:\n", out_dir, "\n")
cat("Done.\n")

# 10) PRINT / SAVE sample label mapping used in heatmap
# This block tells you exactly which raw sample number corresponds
# to each displayed heatmap label.
#
# Logic must match the heatmap:
#   Female and Male are separated.
#   Within each sex:
#     Ctrl_LSD1WT -> Ctrl_LSD1KO -> PLX_LSD1WT -> PLX_LSD1KO
#   Within each group:
#     samples are ordered by original sample number.

make_heatmap_label_check_by_sex <- function(sample_info, sex_value, group_order) {
  info <- sample_info %>%
    dplyr::filter(Sex == sex_value) %>%
    dplyr::mutate(
      Group4 = factor(Group4, levels = group_order),
      Sample_numeric = as.numeric(Sample)
    ) %>%
    dplyr::arrange(Group4, Sample_numeric) %>%
    dplyr::group_by(Group4) %>%
    dplyr::mutate(
      Figure_label = paste0(as.character(Group4), "_", sprintf("%02d", dplyr::row_number()))
    ) %>%
    dplyr::ungroup() %>%
    dplyr::select(
      Figure_label,
      Sample,
      Animal_or_well,
      Full_ID,
      Raw_genotype,
      Sex,
      Treatment,
      LSD1_status,
      Group4
    )

  info
}

female_label_check <- make_heatmap_label_check_by_sex(
  sample_info = geno,
  sex_value = "F",
  group_order = group_order
)

male_label_check <- make_heatmap_label_check_by_sex(
  sample_info = geno,
  sex_value = "M",
  group_order = group_order
)

cat("\n\n============================================================\n")
cat("HEATMAP SAMPLE LABEL CHECK: FEMALE\n")
cat("These are the exact raw sample numbers behind Female heatmap labels.\n")
cat("============================================================\n")
print(female_label_check, n = Inf)

cat("\n\n============================================================\n")
cat("HEATMAP SAMPLE LABEL CHECK: MALE\n")
cat("These are the exact raw sample numbers behind Male heatmap labels.\n")
cat("============================================================\n")
print(male_label_check, n = Inf)

cat("\n\n============================================================\n")
cat("COUNT CHECK BY SEX AND GROUP\n")
cat("============================================================\n")
print(table(geno$Sex, geno$Group4, useNA = "ifany"))

cat("\n\n============================================================\n")
cat("CHECK: Female and Male are separated\n")
cat("============================================================\n")
cat("Female samples in Female label table:\n")
print(unique(female_label_check$Sex))

cat("Male samples in Male label table:\n")
print(unique(male_label_check$Sex))

if (!all(female_label_check$Sex == "F")) {
  stop("ERROR: Female label table contains non-Female samples.")
}

if (!all(male_label_check$Sex == "M")) {
  stop("ERROR: Male label table contains non-Male samples.")
}

cat("\nPASS: Female and Male samples are separated correctly.\n")

cat("\n\n============================================================\n")
cat("CHECK: cre- / cre+ mapping\n")
cat("============================================================\n")

mapping_check <- geno %>%
  dplyr::select(Sample, Raw_genotype, Sex, Treatment, LSD1_status, Group4) %>%
  dplyr::arrange(Sex, Group4, as.numeric(Sample))

print(mapping_check, n = Inf)

bad_cre_minus <- mapping_check %>%
  dplyr::filter(grepl("cre\\-", Raw_genotype, ignore.case = TRUE), LSD1_status != "LSD1WT")

bad_cre_plus <- mapping_check %>%
  dplyr::filter(grepl("cre\\+", Raw_genotype, ignore.case = TRUE), LSD1_status != "LSD1KO")

if (nrow(bad_cre_minus) > 0) {
  print(bad_cre_minus)
  stop("ERROR: Some cre- samples were not mapped to LSD1WT.")
}

if (nrow(bad_cre_plus) > 0) {
  print(bad_cre_plus)
  stop("ERROR: Some cre+ samples were not mapped to LSD1KO.")
}

cat("\nPASS: cre- = LSD1WT and cre+ = LSD1KO mapping is correct.\n")

# Save these checks as CSV files.
readr::write_csv(
  female_label_check,
  file.path(out_dir, "CHECK_heatmap_individual_sample_labels_Female.csv")
)

readr::write_csv(
  male_label_check,
  file.path(out_dir, "CHECK_heatmap_individual_sample_labels_Male.csv")
)

readr::write_csv(
  mapping_check,
  file.path(out_dir, "CHECK_cre_mapping_all_samples.csv")
)

cat("\n\nSaved sample-label check files:\n")
cat(file.path(out_dir, "CHECK_heatmap_individual_sample_labels_Female.csv"), "\n")
cat(file.path(out_dir, "CHECK_heatmap_individual_sample_labels_Male.csv"), "\n")
cat(file.path(out_dir, "CHECK_cre_mapping_all_samples.csv"), "\n")
cat("\nSample label checking complete.\n")

# 11) Mini heatmap: CSF1R / KDM1A / ESR1 only

mini_genes <- c("CSF1R", "KDM1A", "ESR1", "TCF24")

mini_genes_found <- intersect(mini_genes, rownames(z_mat))
mini_genes_missing <- setdiff(mini_genes, rownames(z_mat))

cat("\n============================================================\n")
cat("3-gene mini heatmap check\n")
cat("Requested genes:\n")
print(mini_genes)
cat("Found genes:\n")
print(mini_genes_found)
cat("Missing genes:\n")
print(mini_genes_missing)
cat("============================================================\n")

if (length(mini_genes_found) == 0) {
  stop("None of these genes were found in z_mat: CSF1R, KDM1A, ESR1, TCF24")
}

mini_z_all <- z_mat[mini_genes_found, sample_order, drop = FALSE]

mini_pdf <- file.path(out_dir, "PLX_terminal_3genes_CSF1R_KDM1A_ESR1_TCF24_heatmap_BY_SEX_WITHIN_SEX_NORMALIZED.pdf")
pdf(mini_pdf, width = 10, height = 6)

for (sx in sex_order) {
  cat("\nNow plotting mini heatmap for sex =", sx, "\n")

  sex_result <- make_group_mean_from_individual_z_by_sex(
    ind_z = mini_z_all,
    sample_info = geno,
    group_order = group_order,
    sex_value = sx
  )

  mini_mean_z <- sex_result$mean_z
  mini_ind_z  <- sex_result$ind_z
  info_sex    <- sex_result$sample_info

  info_sex <- info_sex %>%
    dplyr::mutate(
      Group4 = factor(Group4, levels = group_order),
      Sample_numeric = as.numeric(Sample)
    ) %>%
    dplyr::arrange(Group4, Sample_numeric)

  sex_sample_order <- info_sex$Sample
  mini_ind_z <- mini_ind_z[, sex_sample_order, drop = FALSE]

  # shared row order for same-sex group mean + individual
  mix_all <- cbind(mini_mean_z, mini_ind_z)
  if (nrow(mix_all) > 1) {
    hc_all <- hclust(dist(mix_all), method = "ward.D2")
    row_ord <- rownames(mix_all)[hc_all$order]
  } else {
    row_ord <- rownames(mix_all)
  }

  mini_mean_z <- mini_mean_z[row_ord, group_order, drop = FALSE]
  mini_ind_z  <- mini_ind_z[row_ord, sex_sample_order, drop = FALSE]

  group_labels <- make_group_labels(colnames(mini_mean_z))
  ind_labels <- make_individual_labels_by_sex(colnames(mini_ind_z), info_sex)

  ha_mean <- make_group_anno(colnames(mini_mean_z))
  ha_ind  <- make_individual_anno_by_sex(colnames(mini_ind_z), info_sex)

  # group mean page
  ht_mean <- Heatmap(
    mini_mean_z,
    name = "z-score",
    col = col_fun,
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    show_row_dend = FALSE,
    show_column_dend = FALSE,
    column_order = colnames(mini_mean_z),
    column_labels = group_labels,
    top_annotation = ha_mean,
    show_row_names = TRUE,
    row_names_gp = grid::gpar(fontsize = 12, fontface = "bold"),
    show_column_names = TRUE,
    column_names_gp = grid::gpar(fontsize = 10),
    column_title = paste0("CSF1R / KDM1A / ESR1 / TCF24 - ", sex_display_name(sx), " (group mean)"),
    column_title_gp = grid::gpar(fontsize = 14, fontface = "bold"),
    heatmap_legend_param = list(
      title = "z-score",
      title_gp = grid::gpar(fontsize = 10),
      labels_gp = grid::gpar(fontsize = 9)
    )
  )

  draw(
    ht_mean,
    newpage = TRUE,
    heatmap_legend_side = "right",
    annotation_legend_side = "right",
    merge_legends = TRUE
  )

  # individual page
  ht_ind <- Heatmap(
    mini_ind_z,
    name = "z-score",
    col = col_fun,
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    show_row_dend = FALSE,
    show_column_dend = FALSE,
    column_order = colnames(mini_ind_z),
    column_labels = ind_labels,
    top_annotation = ha_ind,
    show_row_names = TRUE,
    row_names_gp = grid::gpar(fontsize = 12, fontface = "bold"),
    show_column_names = TRUE,
    column_names_gp = grid::gpar(fontsize = 6),
    column_title = paste0("CSF1R / KDM1A / ESR1 / TCF24 - ", sex_display_name(sx), " (individual samples)"),
    column_title_gp = grid::gpar(fontsize = 14, fontface = "bold"),
    heatmap_legend_param = list(
      title = "z-score",
      title_gp = grid::gpar(fontsize = 10),
      labels_gp = grid::gpar(fontsize = 9)
    )
  )

  draw(
    ht_ind,
    newpage = TRUE,
    heatmap_legend_side = "right",
    annotation_legend_side = "right",
    merge_legends = TRUE
  )
}

dev.off()

cat("\nSaved 3-gene mini heatmap PDF:\n")
cat(mini_pdf, "\n")
