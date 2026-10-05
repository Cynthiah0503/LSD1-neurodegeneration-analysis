# 09_plx_terminal_reactive_astrocyte_heatmap_by_sex.R
#
# Purpose: PLX terminal sex-stratified reactive-astrocyte heatmap.
# Inputs are expected under data/processed/ or data/external/ relative to this repository.
# Outputs are written under results/ or script-defined subfolders.

# PLX terminal reactive astrocyte marker heatmap
# FULL FIXED VERSION
#
# Fixed from previous error:
#   The previous print-check block selected columns before arrange(),
#   which dropped SexDisplay and caused:
#     object 'SexDisplay' not found
#
# Rules:
# 1) Single-batch project: NO ComBat / NO batch correction
# 2) 1301-1332 featureCounts files have no reliable column names:
#    read with col_names = FALSE and manually assign names
# 3) cre- = LSD1WT; cre+ = LSD1KO
# 4) Gene order is fixed to match the Trem22 reference figure
# 5) Female and Male are shown on the SAME page:
#    Page 1 = group mean: Female left, Male right
#    Page 2 = individual samples: Female left, Male right
# 6) Group mean is calculated FROM individual sample z-scores within each sex
# 7) Gene names are shown ONLY on the far right side of the Male heatmap
# 8) Female and Male heatmaps are separated by a larger middle gap


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

base_dir <- file.path(repo_root, "data", "processed", "PLX terminal experiment/results")
featurecounts_dir <- file.path(base_dir, "featurecounts")
genotype_file <- file.path(base_dir, "genotypes.xlsx")

out_dir <- file.path(base_dir, "heatmap_reactive_astrocyte_shared_geneorder_side_by_side")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

if (!dir.exists(featurecounts_dir)) {
  stop("Missing featureCounts folder: ", featurecounts_dir)
}
if (!file.exists(genotype_file)) {
  stop("Missing genotype file: ", genotype_file)
}

# 1) Fixed reactive astrocyte gene order
#    EXACTLY follows the Trem22 reference figure.
reactive_gene_order_ref <- c(
  "CRYAB",
  "SOX9",
  "MT1",
  "MT2",
  "TSPO",
  "GFAP",
  "VIM",
  "LCN2",
  "THBS1",
  "IL17RA",
  "STAT3",
  "NES",
  "HSPB1",
  "SERPINA3N",
  "ALDOC",
  "S100B",
  "SLC1A2",
  "MAOB",
  "SLC1A3",
  "KCNJ10",
  "C3",
  "FABP7",
  "NFATC3",
  "NFATC4",
  "NTRK2",
  "SYNM"
)

# CHI3L1 is checked but NOT plotted because it is not in the reference figure order.
reactive_marker_check_full <- c(reactive_gene_order_ref, "CHI3L1")

geno_raw <- readxl::read_excel(genotype_file, sheet = 1, col_names = FALSE)

sample_info <- geno_raw %>%
  dplyr::select(
    Animal_or_Well = 1,
    FullSampleID = 2,
    Sample = 3,
    RawGenotype = 4,
    Sex = 5
  ) %>%
  dplyr::filter(!is.na(Sample), !is.na(RawGenotype), !is.na(Sex)) %>%
  dplyr::mutate(
    Sample = as.character(Sample),
    RawGenotype_original = trimws(as.character(RawGenotype)),
    RawGenotype = trimws(tolower(as.character(RawGenotype))),
    Sex = toupper(trimws(as.character(Sex))),
    Treatment = dplyr::case_when(
      grepl("^ctrl", RawGenotype) ~ "Ctrl",
      grepl("^plx", RawGenotype) ~ "PLX",
      TRUE ~ NA_character_
    ),
    LSD1 = dplyr::case_when(
      grepl("cre-", RawGenotype, fixed = TRUE) ~ "LSD1WT",
      grepl("cre+", RawGenotype, fixed = TRUE) ~ "LSD1KO",
      TRUE ~ NA_character_
    ),
    Group = paste(Treatment, LSD1, sep = "_"),
    SexDisplay = dplyr::case_when(
      Sex == "F" ~ "Female",
      Sex == "M" ~ "Male",
      TRUE ~ Sex
    )
  ) %>%
  dplyr::filter(Sample %in% as.character(1301:1332))

if (nrow(sample_info) != 32) {
  warning("Expected 32 samples in genotype file, but found ", nrow(sample_info), ". Check sample_info output.")
}

if (any(is.na(sample_info$Treatment)) || any(is.na(sample_info$LSD1))) {
  print(sample_info %>% dplyr::filter(is.na(Treatment) | is.na(LSD1)))
  stop("Some genotype rows could not be parsed into Treatment/LSD1.")
}

group_order <- c("Ctrl_LSD1WT", "Ctrl_LSD1KO", "PLX_LSD1WT", "PLX_LSD1KO")

sample_info <- sample_info %>%
  dplyr::mutate(
    Treatment = factor(Treatment, levels = c("Ctrl", "PLX")),
    LSD1 = factor(LSD1, levels = c("LSD1WT", "LSD1KO")),
    Group = factor(Group, levels = group_order),
    SexDisplay = factor(SexDisplay, levels = c("Female", "Male"))
  ) %>%
  dplyr::arrange(SexDisplay, Group, as.numeric(Sample))

readr::write_csv(sample_info, file.path(out_dir, "ReactiveAstrocyte_sample_info.csv"))

cat("\n=== Sample counts by Sex and Group ===\n")
print(table(sample_info$SexDisplay, sample_info$Group, useNA = "ifany"))

cat("\n=== cre mapping check ===\n")
mapping_check <- sample_info %>%
  dplyr::arrange(SexDisplay, Group, as.numeric(Sample)) %>%
  dplyr::select(Sample, RawGenotype_original, SexDisplay, Sex, Treatment, LSD1, Group)

print(mapping_check, n = Inf)

readr::write_csv(mapping_check, file.path(out_dir, "CHECK_cre_mapping_all_samples.csv"))

bad_cre_minus <- mapping_check %>%
  dplyr::filter(grepl("cre\\-", RawGenotype_original, ignore.case = TRUE), LSD1 != "LSD1WT")

bad_cre_plus <- mapping_check %>%
  dplyr::filter(grepl("cre\\+", RawGenotype_original, ignore.case = TRUE), LSD1 != "LSD1KO")

if (nrow(bad_cre_minus) > 0) {
  print(bad_cre_minus)
  stop("ERROR: Some cre- samples were not mapped to LSD1WT.")
}

if (nrow(bad_cre_plus) > 0) {
  print(bad_cre_plus)
  stop("ERROR: Some cre+ samples were not mapped to LSD1KO.")
}

cat("\nPASS: cre- = LSD1WT and cre+ = LSD1KO mapping is correct.\n")

read_one_featurecounts <- function(sample_id) {
  f <- file.path(featurecounts_dir, paste0(sample_id, " featurecounts.xlsx"))
  if (!file.exists(f)) {
    stop("Missing featureCounts file: ", f)
  }

  df <- readxl::read_excel(f, sheet = 1, col_names = FALSE)
  df <- as.data.frame(df, stringsAsFactors = FALSE)

  if (ncol(df) < 8) {
    stop("Expected at least 8 columns in ", f, ", but found ", ncol(df))
  }

  df <- df[, 1:8, drop = FALSE]
  names(df) <- c("Geneid", "Count", "Chr", "Start", "End", "Strand", "Extra", "Gene_name")

  first_gene <- tolower(trimws(as.character(df$Geneid[1])))
  first_count <- suppressWarnings(as.numeric(df$Count[1]))
  if (first_gene %in% c("geneid", "gene_id", "gene id") || is.na(first_count)) {
    df <- df[-1, , drop = FALSE]
  }

  df_clean <- df %>%
    dplyr::transmute(
      Geneid = trimws(as.character(Geneid)),
      Chr = as.character(Chr),
      Start = suppressWarnings(as.numeric(Start)),
      End = suppressWarnings(as.numeric(End)),
      Strand = as.character(Strand),
      Gene_name = toupper(trimws(as.character(Gene_name))),
      Count = suppressWarnings(as.numeric(Count))
    ) %>%
    dplyr::filter(!is.na(Geneid), Geneid != "", !is.na(Gene_name), Gene_name != "", Gene_name != "NA") %>%
    dplyr::mutate(Count = ifelse(is.na(Count), 0, Count))

  cleaned_file <- file.path(out_dir, paste0("CLEANED_featureCounts_", sample_id, ".csv"))
  readr::write_csv(df_clean, cleaned_file)

  df_clean %>%
    dplyr::select(Geneid, Gene_name, Count) %>%
    dplyr::rename(!!sample_id := Count)
}

sample_ids <- as.character(1301:1332)
featurecount_files <- file.path(featurecounts_dir, paste0(sample_ids, " featurecounts.xlsx"))
missing_fc <- featurecount_files[!file.exists(featurecount_files)]

if (length(missing_fc) > 0) {
  stop("Missing featureCounts files:\n", paste(missing_fc, collapse = "\n"))
}

fc_list <- lapply(sample_ids, read_one_featurecounts)

merged_geneid <- Reduce(function(x, y) {
  dplyr::full_join(x, y, by = c("Geneid", "Gene_name"))
}, fc_list)

count_cols <- sample_ids
merged_geneid[count_cols] <- lapply(merged_geneid[count_cols], function(x) {
  x <- suppressWarnings(as.numeric(x))
  x[is.na(x)] <- 0
  x
})

readr::write_csv(
  merged_geneid,
  file.path(out_dir, "ReactiveAstrocyte_raw_featureCounts_geneid_1301_1332.csv")
)

merged_symbol <- merged_geneid %>%
  dplyr::filter(!is.na(Gene_name), Gene_name != "", Gene_name != "NA") %>%
  dplyr::group_by(Gene_name) %>%
  dplyr::summarise(
    dplyr::across(dplyr::all_of(count_cols), ~ sum(as.numeric(.x), na.rm = TRUE)),
    .groups = "drop"
  ) %>%
  dplyr::arrange(Gene_name)

readr::write_csv(
  merged_symbol,
  file.path(out_dir, "ReactiveAstrocyte_raw_featureCounts_symbol_1301_1332.csv")
)

expr_mat <- as.matrix(merged_symbol[, count_cols, drop = FALSE])
rownames(expr_mat) <- merged_symbol$Gene_name
storage.mode(expr_mat) <- "numeric"

expr_mat <- expr_mat[, sample_info$Sample, drop = FALSE]

# 5) Keep reference reactive astrocyte genes only
reactive_genes_present <- reactive_gene_order_ref[reactive_gene_order_ref %in% rownames(expr_mat)]
reactive_genes_missing <- setdiff(reactive_gene_order_ref, reactive_genes_present)

marker_check <- data.frame(
  Gene_name = reactive_marker_check_full,
  In_reference_plot_order = reactive_marker_check_full %in% reactive_gene_order_ref,
  Detected_in_featureCounts = reactive_marker_check_full %in% rownames(expr_mat),
  stringsAsFactors = FALSE
)

readr::write_csv(marker_check, file.path(out_dir, "ReactiveAstrocyte_marker_detection_check.csv"))

cat("\n============================================================\n")
cat("Reactive astrocyte fixed gene order check\n")
cat("Reference gene order:\n")
print(reactive_gene_order_ref)
cat("\nGenes found in this PLX dataset and used for plotting:\n")
print(reactive_genes_present)
cat("\nGenes missing from this PLX dataset:\n")
print(reactive_genes_missing)
cat("\nNote: CHI3L1 is checked but not plotted because it is not in the reference figure order.\n")
cat("============================================================\n")

if (length(reactive_genes_present) == 0) {
  stop("None of the reference reactive astrocyte genes were found in the merged matrix.")
}

reactive_counts <- expr_mat[reactive_genes_present, , drop = FALSE]
reactive_log <- log2(reactive_counts + 1)

readr::write_csv(
  data.frame(Gene_name = rownames(reactive_counts), reactive_counts, check.names = FALSE),
  file.path(out_dir, "ReactiveAstrocyte_raw_counts_reference_order_present_markers.csv")
)

readr::write_csv(
  data.frame(Gene_name = rownames(reactive_log), reactive_log, check.names = FALSE),
  file.path(out_dir, "ReactiveAstrocyte_log2_count_plus1_reference_order_present_markers.csv")
)

# 6) Row z-score helper
row_zscore <- function(mat) {
  if (nrow(mat) == 0) return(mat)
  z <- t(scale(t(mat), center = TRUE, scale = TRUE))
  z[!is.finite(z)] <- 0
  z
}

# 7) Group-mean helper: mean from individual z-scores
collapse_group_mean_from_z <- function(mat_z, info_sex) {
  stopifnot(ncol(mat_z) == nrow(info_sex))

  out <- sapply(group_order, function(g) {
    idx <- which(as.character(info_sex$Group) == g)
    if (length(idx) == 0) {
      return(rep(NA_real_, nrow(mat_z)))
    }
    if (length(idx) == 1) {
      return(mat_z[, idx])
    }
    rowMeans(mat_z[, idx, drop = FALSE], na.rm = TRUE)
  })

  out <- as.matrix(out)
  rownames(out) <- rownames(mat_z)
  colnames(out) <- group_order
  out
}

# 8) Display labels
make_individual_display_labels <- function(info_sex) {
  labels <- character(nrow(info_sex))

  for (g in group_order) {
    idx <- which(as.character(info_sex$Group) == g)
    if (length(idx) > 0) {
      labels[idx] <- paste0(g, "_", sprintf("%02d", seq_along(idx)))
    }
  }

  labels
}

# 9) Colors and annotations
lsd1_colors <- c(
  "LSD1WT" = "steelblue",
  "LSD1KO" = "firebrick"
)

treatment_colors <- c(
  "Ctrl" = "grey55",
  "PLX" = "darkorange2"
)

col_fun <- circlize::colorRamp2(
  c(-4, -2, 0, 2, 4),
  c("blue", "lightblue", "white", "salmon", "red")
)

make_anno_from_info <- function(info_sex) {
  HeatmapAnnotation(
    Treatment = factor(info_sex$Treatment, levels = c("Ctrl", "PLX")),
    Group = factor(info_sex$LSD1, levels = c("LSD1WT", "LSD1KO")),
    col = list(
      Treatment = treatment_colors,
      Group = lsd1_colors
    ),
    na_col = "grey90",
    show_annotation_name = TRUE,
    annotation_name_gp = grid::gpar(fontsize = 9)
  )
}

make_anno_from_group_cols <- function(group_cols) {
  treatment <- ifelse(grepl("^Ctrl_", group_cols), "Ctrl",
                      ifelse(grepl("^PLX_", group_cols), "PLX", NA_character_))
  lsd1 <- ifelse(grepl("LSD1WT$", group_cols), "LSD1WT",
                 ifelse(grepl("LSD1KO$", group_cols), "LSD1KO", NA_character_))

  HeatmapAnnotation(
    Treatment = factor(treatment, levels = c("Ctrl", "PLX")),
    Group = factor(lsd1, levels = c("LSD1WT", "LSD1KO")),
    col = list(
      Treatment = treatment_colors,
      Group = lsd1_colors
    ),
    na_col = "grey90",
    show_annotation_name = TRUE,
    annotation_name_gp = grid::gpar(fontsize = 9)
  )
}

# 10) Build sex-specific matrices
build_sex_mats <- function(sex_label) {
  info_sex <- sample_info %>%
    dplyr::filter(as.character(SexDisplay) == sex_label) %>%
    dplyr::arrange(Group, as.numeric(Sample))

  if (nrow(info_sex) == 0) {
    stop("No samples found for sex: ", sex_label)
  }

  log_sex <- reactive_log[reactive_genes_present, info_sex$Sample, drop = FALSE]

  # z-score within each sex-specific matrix.
  z_sex <- row_zscore(log_sex)

  # keep exact reference gene order; do NOT cluster rows.
  z_sex <- z_sex[reactive_genes_present, , drop = FALSE]

  mean_z <- collapse_group_mean_from_z(z_sex, info_sex)
  mean_z <- mean_z[reactive_genes_present, group_order, drop = FALSE]

  colnames(z_sex) <- make_individual_display_labels(info_sex)

  label_check <- info_sex %>%
    dplyr::mutate(Figure_label = make_individual_display_labels(info_sex)) %>%
    dplyr::select(Figure_label, Sample, Animal_or_Well, FullSampleID, RawGenotype_original, SexDisplay, Sex, Treatment, LSD1, Group)

  readr::write_csv(
    label_check,
    file.path(out_dir, paste0("CHECK_heatmap_individual_sample_labels_", sex_label, ".csv"))
  )

  cat("\n=== Heatmap sample labels for ", sex_label, " ===\n", sep = "")
  print(label_check, n = Inf)

  readr::write_csv(
    data.frame(Gene_name = rownames(z_sex), z_sex, check.names = FALSE),
    file.path(out_dir, paste0("ReactiveAstrocyte_", sex_label, "_individual_row_zscore_REFERENCE_ORDER.csv"))
  )

  readr::write_csv(
    data.frame(Gene_name = rownames(mean_z), mean_z, check.names = FALSE),
    file.path(out_dir, paste0("ReactiveAstrocyte_", sex_label, "_group_mean_from_individual_zscore_REFERENCE_ORDER.csv"))
  )

  list(
    info = info_sex,
    individual_z = z_sex,
    group_mean_z = mean_z
  )
}

female <- build_sex_mats("Female")
male <- build_sex_mats("Male")

# 11) Heatmap object builders
make_group_mean_heatmap <- function(mat, sex_label, show_gene_names = FALSE) {
  Heatmap(
    mat,
    name = "z-score",
    col = col_fun,
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    show_row_dend = FALSE,
    show_column_dend = FALSE,
    row_order = seq_len(nrow(mat)),
    column_order = colnames(mat),
    top_annotation = make_anno_from_group_cols(colnames(mat)),
    show_row_names = show_gene_names,
    row_names_side = "right",
    row_names_gp = grid::gpar(fontsize = 8),
    show_column_names = TRUE,
    column_names_gp = grid::gpar(fontsize = 9),
    column_title = paste0("Reactive astrocyte markers - ", sex_label, " (group mean)"),
    column_title_gp = grid::gpar(fontsize = 13, fontface = "bold"),
    heatmap_legend_param = list(
      title = "z-score",
      title_gp = grid::gpar(fontsize = 10),
      labels_gp = grid::gpar(fontsize = 9)
    )
  )
}

make_individual_heatmap <- function(mat, info_sex, sex_label, show_gene_names = FALSE) {
  Heatmap(
    mat,
    name = "z-score",
    col = col_fun,
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    show_row_dend = FALSE,
    show_column_dend = FALSE,
    row_order = seq_len(nrow(mat)),
    column_order = colnames(mat),
    top_annotation = make_anno_from_info(info_sex),
    show_row_names = show_gene_names,
    row_names_side = "right",
    row_names_gp = grid::gpar(fontsize = 8),
    show_column_names = TRUE,
    column_names_gp = grid::gpar(fontsize = 6),
    column_title = paste0("Reactive astrocyte markers - ", sex_label, " (individual samples)"),
    column_title_gp = grid::gpar(fontsize = 13, fontface = "bold"),
    heatmap_legend_param = list(
      title = "z-score",
      title_gp = grid::gpar(fontsize = 10),
      labels_gp = grid::gpar(fontsize = 9)
    )
  )
}

out_pdf <- file.path(
  out_dir,
  "PLX_terminal_ReactiveAstrocyte_SHARED_GENEORDER_Female_Male_side_by_side.pdf"
)

pdf(out_pdf, width = 18, height = 8)

# Page 1: group mean
# Gene names are shown only on the FAR RIGHT side of the Male heatmap.
ht_group_female <- make_group_mean_heatmap(
  mat = female$group_mean_z,
  sex_label = "Female",
  show_gene_names = TRUE
)

ht_group_male <- make_group_mean_heatmap(
  mat = male$group_mean_z,
  sex_label = "Male",
  show_gene_names = TRUE
)

draw(
  ht_group_female + ht_group_male,
  newpage = TRUE,
  merge_legends = TRUE,
  heatmap_legend_side = "right",
  annotation_legend_side = "right",
  ht_gap = grid::unit(2.5, "cm")
)

# Page 2: individual samples
# Gene names are shown only on the FAR RIGHT side of the Male heatmap.
ht_ind_female <- make_individual_heatmap(
  mat = female$individual_z,
  info_sex = female$info,
  sex_label = "Female",
  show_gene_names = TRUE
)

ht_ind_male <- make_individual_heatmap(
  mat = male$individual_z,
  info_sex = male$info,
  sex_label = "Male",
  show_gene_names = TRUE
)

draw(
  ht_ind_female + ht_ind_male,
  newpage = TRUE,
  merge_legends = TRUE,
  heatmap_legend_side = "right",
  annotation_legend_side = "right",
  ht_gap = grid::unit(2.5, "cm")
)

dev.off()

cat("\nSaved PDF:\n", out_pdf, "\n")
cat("\nOutput folder:\n", out_dir, "\n")
cat("\nDone.\n")
