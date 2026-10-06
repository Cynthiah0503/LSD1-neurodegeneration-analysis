# 09_trem2_terminal_reactive_astrocyte_heatmap.R
#
# Purpose: Trem2/LSD1 terminal reactive-astrocyte marker heatmap.
# Inputs are expected under data/processed/ or data/external/ relative to this repository.
# Outputs are written under results/ or script-defined subfolders.

# Trem2 reactive astrocyte markers heatmap
# Single file: Trem2heatmap feature counts.xlsx
# + ONE output PDF only
# + page 1 = 4-group mean
# + page 2 = individual samples
# + SAME shared gene order across both pages
# + individual labels use genotype_01/02/03...
# + group mean is calculated FROM individual z-scores
#   (same logic as previous code)


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
data_processed_dir <- file.path(repo_root, "data", "processed")
data_external_dir <- file.path(repo_root, "data", "external")
results_dir <- file.path(repo_root, "results", "trem2_terminal")
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
library(ComplexHeatmap)
library(dplyr)
library(grid)
library(readxl)
library(circlize)

trem22 <- readxl::read_excel(file.path(data_processed_dir, "Trem2heatmap feature counts.xlsx"), sheet = 1)

# 2) Fix names / BOM
fix_names <- function(df) {
  names(df) <- sub("^\\ufeff", "", names(df))
  bad <- is.na(names(df)) | names(df) == ""
  if (any(bad)) names(df)[bad] <- paste0("V", seq_len(sum(bad)))
  names(df) <- make.unique(names(df))
  df
}

trem22 <- fix_names(trem22)

# 3) Prep featureCounts
prep_featurecounts <- function(df) {
  gene_col <- tail(names(df), 1)
  df$Gene_name <- toupper(trimws(as.character(df[[gene_col]])))

  sample_cols <- grep("^LSD1", names(df), value = TRUE)
  if (length(sample_cols) == 0) {
    stop("No sample columns found starting with 'LSD1'.")
  }

  df_agg <- df[, c("Gene_name", sample_cols), drop = FALSE] %>%
    dplyr::filter(!is.na(Gene_name), Gene_name != "") %>%
    dplyr::group_by(Gene_name) %>%
    dplyr::summarise(
      dplyr::across(dplyr::all_of(sample_cols), ~ sum(as.numeric(.x), na.rm = TRUE)),
      .groups = "drop"
    )

  list(df_agg = df_agg, sample_cols = sample_cols)
}

trem22_p <- prep_featurecounts(trem22)

trem22_mat <- as.matrix(trem22_p$df_agg[, trem22_p$sample_cols, drop = FALSE])
rownames(trem22_mat) <- trem22_p$df_agg$Gene_name
storage.mode(trem22_mat) <- "numeric"

# 4) Reactive astrocyte marker list
reast <- c(
  "Gfap", "Nes", "Synm", "Vim",
  "Aldoc", "Fabp7", "Maob", "Tspo",
  "Cryab", "Hspb1",
  "C3", "Chi3l1", "Lcn2", "Serpina3n", "Mt1", "Mt2", "Thbs1",
  "Nfatc3", "Nfatc4", "Ntrk2", "Il17ra", "S100b", "Sox9", "Stat3",
  "Slc1a3", "Slc1a2", "Kcnj10"
)
reast_up <- toupper(reast)

# 5) Parse sample names into 4 groups
#    internal uses WW/MM
#    display uses WT/KO
make_group_id_4_from_samples <- function(samples) {
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

  out <- ifelse(!is.na(lsd1_part) & !is.na(trem2_part),
                paste(lsd1_part, trem2_part, sep = "_"),
                NA_character_)
  out
}

convert_internal_to_display <- function(x) {
  x <- gsub("Trem2WW", "Trem2WT", x, fixed = TRUE)
  x <- gsub("Trem2MM", "Trem2KO", x, fixed = TRUE)
  x
}

fixed_group_order_internal <- c(
  "LSD1WT_Trem2WW",
  "LSD1WT_Trem2MM",
  "LSD1KO_Trem2WW",
  "LSD1KO_Trem2MM"
)

fixed_group_order_display <- convert_internal_to_display(fixed_group_order_internal)

# 6) Collapse to 4-group mean
collapse_mat_by_group_mean <- function(mat, group_id, group_order) {
  stopifnot(length(group_id) == ncol(mat))

  ok <- !is.na(group_id)
  mat2 <- mat[, ok, drop = FALSE]
  gid  <- group_id[ok]

  if (ncol(mat2) == 0) {
    stop("No samples mapped to 4 groups. Check sample names.")
  }

  out <- sapply(group_order, function(g) {
    idx <- which(gid == g)
    if (length(idx) == 0) return(rep(NA_real_, nrow(mat2)))
    if (length(idx) == 1) return(mat2[, idx])
    rowMeans(mat2[, idx, drop = FALSE], na.rm = TRUE)
  })

  out <- as.matrix(out)
  rownames(out) <- rownames(mat2)
  colnames(out) <- group_order
  out
}

gid4 <- make_group_id_4_from_samples(colnames(trem22_mat))

cat("\n=== 4-group mapping counts ===\n")
print(table(gid4, useNA = "ifany"))

if (any(is.na(gid4))) {
  cat("\nUnmapped sample names:\n")
  print(colnames(trem22_mat)[is.na(gid4)])
}

# 7) Build individual labels: genotype_01/02/03...
make_individual_order_and_labels <- function(samples) {
  gid_internal <- make_group_id_4_from_samples(samples)

  ord <- unlist(lapply(fixed_group_order_internal, function(g) which(gid_internal == g)))
  if (length(ord) == 0) {
    stop("No samples matched the fixed 4-group order.")
  }

  ordered_samples <- samples[ord]
  ordered_gid_internal <- gid_internal[ord]
  ordered_gid_display  <- convert_internal_to_display(ordered_gid_internal)

  labels <- character(length(ordered_samples))

  for (g in fixed_group_order_display) {
    idx <- which(ordered_gid_display == g)
    if (length(idx) > 0) {
      labels[idx] <- paste0(g, "_", sprintf("%02d", seq_along(idx)))
    }
  }

  list(
    ordered_samples = ordered_samples,
    ordered_gid_internal = ordered_gid_internal,
    ordered_gid_display = ordered_gid_display,
    display_labels = labels
  )
}

# 8) Keep reactive markers only
keep_ra <- intersect(rownames(trem22_mat), reast_up)
if (length(keep_ra) == 0) {
  stop("No Reactive Astrocyte markers found in trem22_mat.")
}

ra_ind_mat <- trem22_mat[keep_ra, , drop = FALSE]

# 9) log2 + row z-score
#    do this ONCE on individual samples
row_zscore <- function(mat) {
  z <- t(scale(t(mat), center = TRUE, scale = TRUE))
  z[!is.finite(z)] <- 0
  z
}

ra_ind_log <- log2(ra_ind_mat + 1)
ra_ind_z   <- row_zscore(ra_ind_log)

# 10) group mean FROM z-score matrix
#     do NOT z-score group mean again
ra_mean_z <- collapse_mat_by_group_mean(
  ra_ind_z,
  gid4,
  group_order = fixed_group_order_internal
)

common_ra_genes <- intersect(rownames(ra_ind_z), rownames(ra_mean_z))
if (length(common_ra_genes) == 0) {
  stop("No shared Reactive Astrocyte genes between individual and 4-group mean matrices.")
}

ra_ind_z2  <- ra_ind_z[common_ra_genes, , drop = FALSE]
ra_mean_z2 <- ra_mean_z[common_ra_genes, , drop = FALSE]

# 11) Shared row order
#     built from mean + individual together
mix_ra <- cbind(
  ra_mean_z2[, fixed_group_order_internal, drop = FALSE],
  ra_ind_z2
)

if (nrow(mix_ra) == 1) {
  ra_gene_order_shared <- rownames(mix_ra)
} else {
  hc_ra <- hclust(dist(mix_ra), method = "ward.D2")
  ra_gene_order_shared <- rownames(mix_ra)[hc_ra$order]
}

ra_mean_z_shared <- ra_mean_z2[ra_gene_order_shared, fixed_group_order_internal, drop = FALSE]
ra_ind_z_shared  <- ra_ind_z2[ra_gene_order_shared, , drop = FALSE]

# 12) Fixed column order
ind_info <- make_individual_order_and_labels(colnames(ra_ind_z_shared))
ra_ind_z_shared <- ra_ind_z_shared[, ind_info$ordered_samples, drop = FALSE]

# 13) Column labels for display
mean_col_labels <- fixed_group_order_display
ind_col_labels  <- ind_info$display_labels

# 14) Top annotations
group_colors_cre <- c(
  "LSD1WT" = "steelblue",
  "LSD1KO" = "firebrick"
)

trem2_colors <- c(
  "Trem2WT" = "#8E9AAF",
  "Trem2KO" = "#C77DFF"
)

make_sample_anno_4groups <- function(display_cols) {
  s <- as.character(display_cols)

  grp <- ifelse(grepl("^LSD1WT_", s), "LSD1WT",
                ifelse(grepl("^LSD1KO_", s), "LSD1KO", NA_character_))

  trem2 <- ifelse(grepl("_Trem2WT$", s), "Trem2WT",
                  ifelse(grepl("_Trem2KO$", s), "Trem2KO", NA_character_))

  HeatmapAnnotation(
    Group = factor(grp, levels = c("LSD1WT", "LSD1KO")),
    Trem2 = factor(trem2, levels = c("Trem2WT", "Trem2KO")),
    col = list(
      Group = group_colors_cre,
      Trem2 = trem2_colors
    ),
    na_col = "grey90",
    show_annotation_name = TRUE,
    annotation_name_gp = grid::gpar(fontsize = 10)
  )
}

make_sample_anno_individual <- function(display_cols) {
  s <- as.character(display_cols)

  grp <- ifelse(grepl("^LSD1WT_", s), "LSD1WT",
                ifelse(grepl("^LSD1KO_", s), "LSD1KO", NA_character_))

  trem2 <- ifelse(grepl("_Trem2WT_", s), "Trem2WT",
                  ifelse(grepl("_Trem2KO_", s), "Trem2KO", NA_character_))

  HeatmapAnnotation(
    Group = factor(grp, levels = c("LSD1WT", "LSD1KO")),
    Trem2 = factor(trem2, levels = c("Trem2WT", "Trem2KO")),
    col = list(
      Group = group_colors_cre,
      Trem2 = trem2_colors
    ),
    na_col = "grey90",
    show_annotation_name = TRUE,
    annotation_name_gp = grid::gpar(fontsize = 10)
  )
}

ha_ra_mean <- make_sample_anno_4groups(mean_col_labels)
ha_ra_ind  <- make_sample_anno_individual(ind_col_labels)

# 15) Shared color scale
col_fun <- circlize::colorRamp2(
  c(-4, -2, 0, 2, 4),
  c("blue", "lightblue", "white", "salmon", "red")
)

# 16) Heatmaps
ht_ra_mean <- Heatmap(
  ra_mean_z_shared,
  name = "z-score",
  col = col_fun,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_row_dend = FALSE,
  show_column_dend = FALSE,
  row_order = seq_len(nrow(ra_mean_z_shared)),
  column_order = fixed_group_order_internal,
  column_labels = mean_col_labels,
  top_annotation = ha_ra_mean,
  show_row_names = TRUE,
  row_names_side = "right",
  row_names_gp = grid::gpar(fontsize = 8),
  column_names_gp = grid::gpar(fontsize = 10),
  column_title = "Reactive Astrocyte markers (group mean)",
  column_title_gp = grid::gpar(fontsize = 14, fontface = "bold"),
  heatmap_legend_param = list(
    title = "z-score",
    title_gp = grid::gpar(fontsize = 10),
    labels_gp = grid::gpar(fontsize = 9)
  )
)

ht_ra_ind <- Heatmap(
  ra_ind_z_shared,
  name = "z-score",
  col = col_fun,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_row_dend = FALSE,
  show_column_dend = FALSE,
  row_order = seq_len(nrow(ra_ind_z_shared)),
  column_order = colnames(ra_ind_z_shared),
  column_labels = ind_col_labels,
  top_annotation = ha_ra_ind,
  show_row_names = TRUE,
  row_names_side = "right",
  row_names_gp = grid::gpar(fontsize = 8),
  column_names_gp = grid::gpar(fontsize = 8),
  column_title = "Reactive Astrocyte markers (individual samples)",
  column_title_gp = grid::gpar(fontsize = 14, fontface = "bold"),
  heatmap_legend_param = list(
    title = "z-score",
    title_gp = grid::gpar(fontsize = 10),
    labels_gp = grid::gpar(fontsize = 9)
  )
)

out_pdf <- "Trem22_ReactiveAstrocyte_SHARED_GENEORDER_groupmean_and_individual.pdf"

pdf(out_pdf, width = 12, height = 6)

draw(
  ht_ra_mean,
  newpage = TRUE,
  heatmap_legend_side = "right",
  annotation_legend_side = "right",
  merge_legends = TRUE
)

draw(
  ht_ra_ind,
  newpage = TRUE,
  heatmap_legend_side = "right",
  annotation_legend_side = "right",
  merge_legends = TRUE
)

dev.off()

cat("Saved:", out_pdf, "\n")
cat("Reactive markers found:", length(common_ra_genes), "/", length(reast_up), "\n")
cat("4-group mean display order:", paste(mean_col_labels, collapse = " | "), "\n")
cat("Individual sample display labels:\n")
print(ind_col_labels)
