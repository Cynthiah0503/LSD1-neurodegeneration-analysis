# 02_trem2_terminal_microglia_immune_heatmap.R
#
# Purpose: Trem2/LSD1 terminal microglia and immune-module heatmap from featureCounts matrices.
# Inputs are expected under data/processed/ or data/external/ relative to this repository.
# Outputs are written under results/ or script-defined subfolders.

# Trem2 immune heatmap
# Single file: Trem2heatmap feature counts.xlsx
#
# 1) Change title for Yellow to:
#    Microglial and immune genes
# 2) Change display labels:
#    LSD1WT_Trem2WW -> LSD1WT_Trem2WT
#    LSD1WT_Trem2MM -> LSD1WT_Trem2KO
#    LSD1KO_Trem2WW -> LSD1KO_Trem2WT
#    LSD1KO_Trem2MM -> LSD1KO_Trem2KO
# 3) Delete:
#    - gene names
#    - yellow/module side bar on the left
#    - module legend
# 4) Individual sample version:
#    use genotype_01/02/03... naming logic
# 5) Generate ONE PDF ONLY
# 6) Page order:
#    Yellow (group mean)
#    Yellow (individual sample)
#    Light cyan (group mean)
#    Light cyan (individual sample)
#    ...
# 7) IMPORTANT LOGIC:
#    group mean is calculated FROM individual sample z-scores
# 8) NEW:
#    export exact heatmap gene order to CSV
# 9) NEW:
#    add left-side clean row index axis, using 0-based row numbers
#    axis labels use clean rounded ticks only; the final row is not labeled unless it is a clean tick


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
library(ComplexHeatmap)
library(dplyr)
library(grid)
library(readxl)
library(circlize)

trem22 <- readxl::read_excel("Trem2heatmap feature counts.xlsx", sheet = 1)

allc <- read.csv("NIHMS472534-supplement-02.csv", header = TRUE, check.names = FALSE)
mod  <- read.csv("modules.csv", header = TRUE, check.names = FALSE)

# 修复 modules.csv 里 NA/空列名
names(mod) <- sub("^\\ufeff", "", names(mod))
bad <- is.na(names(mod)) | names(mod) == ""
if (any(bad)) names(mod)[bad] <- paste0("V", seq_len(sum(bad)))
names(mod) <- make.unique(names(mod))

names(allc) <- sub("^\\ufeff", "", names(allc))
allc$Gene_Symbol <- toupper(as.character(allc$Gene_Symbol))

# 3) "Immune functions" -> target modules
target_term <- "Immune functions"

target_modules <- mod %>%
  dplyr::filter(trimws(CategoryTerm) == target_term) %>%
  dplyr::pull(Module) %>%
  unique()

# 4) Gene -> Module map
gene2module <- allc %>%
  dplyr::select(Gene_Symbol, Module) %>%
  dplyr::filter(!is.na(Gene_Symbol), !is.na(Module)) %>%
  dplyr::distinct()

# 5) Prep featureCounts
prep_featurecounts <- function(df) {
  names(df) <- sub("^\\ufeff", "", names(df))
  bad <- is.na(names(df)) | names(df) == ""
  if (any(bad)) names(df)[bad] <- paste0("V", seq_len(sum(bad)))
  names(df) <- make.unique(names(df))

  gene_col <- tail(names(df), 1)
  df$Gene_name <- toupper(trimws(as.character(df[[gene_col]])))

  sample_cols <- grep("^LSD1", names(df), value = TRUE)

  df_agg <- df[, c("Gene_name", sample_cols), drop = FALSE] %>%
    dplyr::filter(!is.na(Gene_name), Gene_name != "") %>%
    dplyr::group_by(Gene_name) %>%
    dplyr::summarise(
      dplyr::across(dplyr::all_of(sample_cols), ~ sum(as.numeric(.x), na.rm = TRUE)),
      .groups = "drop"
    )

  list(df_agg = df_agg, sample_cols = sample_cols)
}

# 6) Build expression matrix
trem22_p <- prep_featurecounts(trem22)
trem22_mat <- as.matrix(trem22_p$df_agg[, trem22_p$sample_cols, drop = FALSE])
rownames(trem22_mat) <- trem22_p$df_agg$Gene_name
storage.mode(trem22_mat) <- "numeric"

# 7) Parse sample names
#    Internal parsing still uses WW/MM
#    Display uses WT/KO
make_group_id_internal <- function(samples) {
  s <- as.character(samples)

  is_lsd1ko <- grepl("LSD1\\s*KO|LSD1[_ ]?KO", s, ignore.case = TRUE)
  is_lsd1wt <- grepl("LSD1\\s*WT|LSD1[_ ]?WT", s, ignore.case = TRUE)

  lsd1_part <- dplyr::case_when(
    is_lsd1ko ~ "LSD1KO",
    is_lsd1wt ~ "LSD1WT",
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

# 8) Fixed group order
group_order_internal <- c(
  "LSD1WT_Trem2WW",
  "LSD1WT_Trem2MM",
  "LSD1KO_Trem2WW",
  "LSD1KO_Trem2MM"
)

group_order_display <- convert_internal_to_display(group_order_internal)

# 9) Collapse to 4-group mean
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

gid4_internal <- make_group_id_internal(colnames(trem22_mat))

cat("\n=== 4-group mapping counts ===\n")
print(table(gid4_internal, useNA = "ifany"))

if (any(is.na(gid4_internal))) {
  cat("\nUnmapped sample names (first 20):\n")
  print(head(colnames(trem22_mat)[is.na(gid4_internal)], 20))
}

# 10) Individual sample order + relabel
#     genotype_01 / 02 / 03 ...
make_individual_display_labels <- function(samples,
                                           group_order_internal = c("LSD1WT_Trem2WW",
                                                                    "LSD1WT_Trem2MM",
                                                                    "LSD1KO_Trem2WW",
                                                                    "LSD1KO_Trem2MM")) {
  gid_internal <- make_group_id_internal(samples)

  ord <- unlist(lapply(group_order_internal, function(g) which(gid_internal == g)))
  ordered_samples <- samples[ord]
  ordered_gid_internal <- gid_internal[ord]
  ordered_gid_display  <- convert_internal_to_display(ordered_gid_internal)

  labels <- character(length(ordered_samples))

  for (g in unique(ordered_gid_display)) {
    idx <- which(ordered_gid_display == g)
    labels[idx] <- paste0(g, "_", sprintf("%02d", seq_along(idx)))
  }

  list(
    ordered_samples = ordered_samples,
    ordered_gid_internal = ordered_gid_internal,
    ordered_gid_display = ordered_gid_display,
    display_labels = labels
  )
}

# 11) Top annotations
group_colors <- c(
  "LSD1WT" = "steelblue",
  "LSD1KO" = "firebrick"
)

trem2_colors <- c(
  "Trem2WT" = "#8E9AAF",
  "Trem2KO" = "#C77DFF"
)

make_sample_anno_4groups <- function(cols4_display) {
  s <- as.character(cols4_display)

  grp <- ifelse(grepl("^LSD1WT_", s), "LSD1WT",
                ifelse(grepl("^LSD1KO_", s), "LSD1KO", NA_character_))

  trem2 <- ifelse(grepl("_Trem2WT$", s), "Trem2WT",
                  ifelse(grepl("_Trem2KO$", s), "Trem2KO", NA_character_))

  HeatmapAnnotation(
    Group = factor(grp, levels = c("LSD1WT", "LSD1KO")),
    Trem2 = factor(trem2, levels = c("Trem2WT", "Trem2KO")),
    col = list(
      Group = group_colors,
      Trem2 = trem2_colors
    ),
    na_col = "grey90",
    show_annotation_name = TRUE,
    annotation_name_gp = grid::gpar(fontsize = 10)
  )
}

make_sample_anno_individual <- function(display_labels) {
  s <- as.character(display_labels)

  grp <- ifelse(grepl("^LSD1WT_", s), "LSD1WT",
                ifelse(grepl("^LSD1KO_", s), "LSD1KO", NA_character_))

  trem2 <- ifelse(grepl("_Trem2WT_", s), "Trem2WT",
                  ifelse(grepl("_Trem2KO_", s), "Trem2KO", NA_character_))

  HeatmapAnnotation(
    Group = factor(grp, levels = c("LSD1WT", "LSD1KO")),
    Trem2 = factor(trem2, levels = c("Trem2WT", "Trem2KO")),
    col = list(
      Group = group_colors,
      Trem2 = trem2_colors
    ),
    na_col = "grey90",
    show_annotation_name = TRUE,
    annotation_name_gp = grid::gpar(fontsize = 10)
  )
}

# 12) Left-side clean row index axis
#     Uses 0-based index:
#       first row = 0
#       last row  = n - 1
#
#     Aesthetic rule:
#       show clean rounded ticks only.
#       Do NOT force-label the final row if it is not a clean number.
#       The exact final row index is still saved in the gene-order CSV.
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

# 13) Title helper
make_module_title <- function(module_name, page_type = c("group mean", "individual sample")) {
  page_type <- match.arg(page_type)

  if (module_name == "Yellow") {
    main_title <- "Microglial and immune genes"
  } else {
    main_title <- paste0("Trem2 - Immune module: ", module_name)
  }

  suffix <- if (page_type == "group mean") {
    " (group mean)"
  } else {
    " (individual samples)"
  }

  paste0(main_title, suffix)
}

# 14) Choose modules and page order
priority_modules <- c("Yellow", "Light cyan")
remaining_modules <- setdiff(target_modules, priority_modules)
modules_to_plot <- c(intersect(priority_modules, target_modules), remaining_modules)

cat("\nModules to plot in order:\n")
print(modules_to_plot)

# 15) Heatmap color function
col_fun <- circlize::colorRamp2(
  c(-4, -2, 0, 2, 4),
  c("blue", "lightblue", "white", "salmon", "red")
)

# 16) ONE PDF ONLY
#     For each module:
#       page 1 = group mean
#       page 2 = individual sample
#     Shared gene order across both pages
#       group mean is calculated FROM individual z-scores
out_pdf <- "Trem22_Immune_Combined_OnePDF_GroupMean_Individual_with_LeftCleanRowIndex.pdf"

# Store exact gene order for every module
all_gene_orders <- list()

pdf(out_pdf, width = 11, height = 8)

for (m in modules_to_plot) {
  message("→ plotting module: ", m)

  # genes in this module
  genes_m <- gene2module %>%
    dplyr::filter(Module == m) %>%
    dplyr::pull(Gene_Symbol) %>%
    unique()

  if (length(genes_m) == 0) {
    message("Skipped module with no mapped genes: ", m)
    next
  }

  # Individual sample matrix
  keep_ind <- intersect(rownames(trem22_mat), genes_m)
  if (length(keep_ind) == 0) {
    message("Skipped individual page (no genes found): ", m)
    next
  }

  ind_mat <- trem22_mat[keep_ind, , drop = FALSE]
  ind_log <- log2(ind_mat + 1)
  ind_z <- t(scale(t(ind_log), center = TRUE, scale = TRUE))
  ind_z[!is.finite(ind_z)] <- 0

  # Group mean matrix
  # mean of individual z-scores
  gid_module <- make_group_id_internal(colnames(ind_z))

  mean_z <- collapse_mat_by_group_mean(
    ind_z,
    gid_module,
    group_order = group_order_internal
  )
  colnames(mean_z) <- group_order_display

  # Shared genes + shared row order
  common_genes <- intersect(rownames(ind_z), rownames(mean_z))
  if (length(common_genes) == 0) {
    message("Skipped module (no shared genes between individual and mean): ", m)
    next
  }

  ind_z2  <- ind_z[common_genes, , drop = FALSE]
  mean_z2 <- mean_z[common_genes, , drop = FALSE]

  mix_all <- cbind(mean_z2, ind_z2)
  hc_all  <- hclust(dist(mix_all), method = "ward.D2")
  row_ord <- rownames(mix_all)[hc_all$order]

  # Clean row-axis ticks for this module
  axis_ticks <- get_clean_axis_ticks(length(row_ord))
  axis_tick_step <- choose_clean_tick_step(length(row_ord))
  axis_last_labeled_index <- max(axis_ticks)

  # Print exact gene order to console
  cat("\n=== Gene order for module:", m, "===\n")
  cat("Total genes:", length(row_ord), "\n")
  cat("Clean row-index ticks shown on heatmap:", paste(axis_ticks, collapse = ", "), "\n")
  print(data.frame(
    Module = m,
    Total_genes_in_module = length(row_ord),
    Row_index_0_based = seq_along(row_ord) - 1,
    Row_order_1_based = seq_along(row_ord),
    Axis_tick_step = axis_tick_step,
    Axis_last_labeled_index = axis_last_labeled_index,
    Gene = row_ord
  ))

  # Store exact gene order for final combined CSV
  all_gene_orders[[m]] <- data.frame(
    Module = m,
    Total_genes_in_module = length(row_ord),
    Row_index_0_based = seq_along(row_ord) - 1,
    Row_order_1_based = seq_along(row_ord),
    Axis_tick_step = axis_tick_step,
    Axis_last_labeled_index = axis_last_labeled_index,
    Gene = row_ord
  )

  mean_z_shared <- mean_z2[row_ord, , drop = FALSE]
  ind_z_shared  <- ind_z2[row_ord, , drop = FALSE]

  # Column order
  mean_z_shared <- mean_z_shared[, group_order_display, drop = FALSE]

  label_info <- make_individual_display_labels(colnames(ind_z_shared), group_order_internal)
  ind_z_shared <- ind_z_shared[, label_info$ordered_samples, drop = FALSE]
  colnames(ind_z_shared) <- label_info$display_labels

  # Annotations
  ha_mean <- make_sample_anno_4groups(colnames(mean_z_shared))
  ha_ind  <- make_sample_anno_individual(colnames(ind_z_shared))

  # Same left-side row index logic for group mean and individual pages
  # Create separate annotation objects to avoid reusing the same ComplexHeatmap object twice.
  row_axis_mean <- make_row_index_axis(nrow(mean_z_shared))
  row_axis_ind  <- make_row_index_axis(nrow(ind_z_shared))

  # Page 1: group mean
  ht_mean <- Heatmap(
    mean_z_shared,
    name = "z-score",
    col = col_fun,
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    show_row_dend = FALSE,
    show_column_dend = FALSE,
    column_order = colnames(mean_z_shared),
    top_annotation = ha_mean,
    left_annotation = row_axis_mean,
    show_row_names = FALSE,
    show_column_names = TRUE,
    column_names_gp = grid::gpar(fontsize = 10),
    column_title = make_module_title(m, "group mean"),
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

  # Page 2: individual sample
  ht_ind <- Heatmap(
    ind_z_shared,
    name = "z-score",
    col = col_fun,
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    show_row_dend = FALSE,
    show_column_dend = FALSE,
    column_order = colnames(ind_z_shared),
    top_annotation = ha_ind,
    left_annotation = row_axis_ind,
    show_row_names = FALSE,
    show_column_names = TRUE,
    column_names_gp = grid::gpar(fontsize = 8),
    column_title = make_module_title(m, "individual sample"),
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

# 17) Export exact gene order
gene_order_df <- dplyr::bind_rows(all_gene_orders)

write.csv(
  gene_order_df,
  "ALL_HEATMAP_GENE_ORDER_with_LeftCleanRowIndex.csv",
  row.names = FALSE
)

cat("\nSaved:", out_pdf, "\n")
cat("Saved gene order file: ALL_HEATMAP_GENE_ORDER_with_LeftCleanRowIndex.csv\n")
cat("Done.\n")
