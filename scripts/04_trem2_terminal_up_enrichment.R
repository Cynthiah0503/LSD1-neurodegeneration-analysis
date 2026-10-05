# 04_trem2_terminal_up_enrichment.R
#
# Purpose: Trem2/LSD1 terminal enrichment analysis for UP-regulated gene lists.
# Inputs are expected under data/processed/ or data/external/ relative to this repository.
# Outputs are written under results/ or script-defined subfolders.

# 3组版本
# FINAL PIPELINE:
# 3 Trem2 UP gene groups only
#
# Analyses:
#   1) Cell marker enrichment
#   2) GO BP enrichment
#   3) KEGG enrichment
#   4) Reactome enrichment
#
# Correct cell marker file:
#   Cell_marker_Mouse.xlsx
#   sheet = "mouse"
#   gene column = Symbol
#   marker/category column = cell_name


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
library(clusterProfiler)
library(org.Mm.eg.db)
library(ReactomePA)
library(dplyr)
library(ggplot2)
library(tidyr)
library(stringr)
library(readxl)
library(grid)

select  <- dplyr::select
filter  <- dplyr::filter
mutate  <- dplyr::mutate
arrange <- dplyr::arrange

# 1) Input files: 3 Trem2 groups

files <- list(
  Trem2KO_LSD1WT_vs_Trem2WT_LSD1WT =
    "UP_Trem2KO_LSD1_WT_vs_Trem2_WT_LSD1_WT.csv",

  Trem2WT_LSD1KO_vs_Trem2WT_LSD1WT =
    "UP_Trem2_WT_LSD1_Del_vs_Trem2_WT_LSD1_WT.csv",

  Trem2KO_LSD1KO_vs_Trem2WT_LSD1WT =
    "UP_Trem2KO_LSD1_Del_vs_Trem2_WT_LSD1_WT.csv"
)

cell_marker_file <- "Cell_marker_Mouse.xlsx"

clean_genes <- function(x) {
  x <- as.character(x)
  x <- trimws(x)
  x <- x[!is.na(x) & x != ""]
  unique(x)
}

convert_entrez <- function(genes) {
  df <- suppressWarnings(
    bitr(
      genes,
      fromType = "SYMBOL",
      toType = "ENTREZID",
      OrgDb = org.Mm.eg.db
    )
  )

  df <- df %>%
    distinct(ENTREZID, .keep_all = TRUE)

  return(df$ENTREZID)
}

format_cluster_names_exact <- function(x) {
  map <- c(
    "Trem2WT_LSD1KO_vs_Trem2WT_LSD1WT" =
      "Trem2WT_LSD1KO\nvs\nTrem2WT_LSD1WT",

    "Trem2KO_LSD1KO_vs_Trem2WT_LSD1WT" =
      "Trem2KO_LSD1KO\nvs\nTrem2WT_LSD1WT",

    "Trem2KO_LSD1WT_vs_Trem2WT_LSD1WT" =
      "Trem2KO_LSD1WT\nvs\nTrem2WT_LSD1WT"
  )

  out <- unname(map[x])
  out[is.na(out)] <- x[is.na(out)]
  return(out)
}

cluster_order_raw <- names(files)
cluster_order_clean <- format_cluster_names_exact(cluster_order_raw)

# 3) Build gene lists
gene_list_symbol <- list()
gene_list_entrez <- list()

cat("\n========== BUILD 3-GROUP TREM2 GENE LIST ==========\n")

for (nm in names(files)) {

  cat("\nProcessing:", nm, "\n")

  df <- read.csv(files[[nm]], check.names = FALSE)

  if (!("Gene_name" %in% names(df))) {
    cat("Missing Gene_name column -> skip\n")
    next
  }

  genes <- clean_genes(df$Gene_name)
  cat("Total SYMBOL genes:", length(genes), "\n")

  if (length(genes) < 5) {
    cat("Too few genes -> skip enrichment\n")
    next
  }

  entrez <- convert_entrez(genes)

  cat("Mapped ENTREZ genes:", length(entrez), "\n")
  cat("Mapping rate:", round(length(entrez) / length(genes) * 100, 2), "%\n")

  if (length(entrez) < 5) {
    cat("Too few mapped genes -> skip enrichment\n")
    next
  }

  gene_list_symbol[[nm]] <- genes
  gene_list_entrez[[nm]] <- entrez
}

cat("\nFinal groups used for enrichment:\n")
print(names(gene_list_entrez))

# 4) CELL MARKER ENRICHMENT

cat("\n========== RUN CELL MARKER ENRICHMENT ==========\n")

cell_marker_raw <- readxl::read_excel(cell_marker_file, sheet = "mouse")

cell_marker_term2gene <- cell_marker_raw %>%
  dplyr::select(cell_name, tissue_type, Symbol) %>%
  dplyr::mutate(
    cell_name   = trimws(as.character(cell_name)),
    tissue_type = trimws(as.character(tissue_type)),
    Symbol      = trimws(as.character(Symbol)),
    Cell_marker = paste(cell_name, tissue_type, "Mouse", sep = " "),
    Gene        = Symbol
  ) %>%
  dplyr::filter(!is.na(cell_name), cell_name != "") %>%
  dplyr::filter(!is.na(tissue_type), tissue_type != "") %>%
  dplyr::filter(!is.na(Gene), Gene != "") %>%
  dplyr::distinct(Cell_marker, Gene)

cell_compare <- compareCluster(
  geneCluster = gene_list_symbol,
  fun = "enricher",
  TERM2GENE = cell_marker_term2gene,
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.2
)

write.csv(
  as.data.frame(cell_compare),
  "COMPARE_CELL_MARKER_3_TREM2_GROUPS.csv",
  row.names = FALSE
)

# 5) GO BP ENRICHMENT

cat("\n========== RUN GO BP COMPARE ==========\n")

ego_compare <- compareCluster(
  geneCluster = gene_list_entrez,
  fun = "enrichGO",
  OrgDb = org.Mm.eg.db,
  keyType = "ENTREZID",
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.2
)

write.csv(
  as.data.frame(ego_compare),
  "COMPARE_GO_BP_3_TREM2_GROUPS.csv",
  row.names = FALSE
)

# 6) KEGG ENRICHMENT

cat("\n========== RUN KEGG COMPARE ==========\n")

kegg_compare <- compareCluster(
  geneCluster = gene_list_entrez,
  fun = "enrichKEGG",
  organism = "mmu",
  pvalueCutoff = 0.05
)

write.csv(
  as.data.frame(kegg_compare),
  "COMPARE_KEGG_3_TREM2_GROUPS.csv",
  row.names = FALSE
)

# 7) REACTOME ENRICHMENT

cat("\n========== RUN REACTOME COMPARE ==========\n")

react_compare <- compareCluster(
  geneCluster = gene_list_entrez,
  fun = "enrichPathway",
  organism = "mouse",
  pvalueCutoff = 0.05
)

write.csv(
  as.data.frame(react_compare),
  "COMPARE_REACTOME_3_TREM2_GROUPS.csv",
  row.names = FALSE
)

# 8) HEATMAP FUNCTION
# Missing enrichment = 0

make_enrichment_heatmap <- function(input_csv,
                                    output_pdf,
                                    plot_title,
                                    top_n_pathways = 30) {

  df <- read.csv(input_csv, check.names = FALSE)

  required_cols <- c("Cluster", "Description", "p.adjust")
  missing_cols <- setdiff(required_cols, names(df))
  if (length(missing_cols) > 0) {
    stop(paste("Missing columns:", paste(missing_cols, collapse = ", ")))
  }

  df2 <- df %>%
    filter(!is.na(Cluster), !is.na(Description), !is.na(p.adjust)) %>%
    mutate(
      Cluster_clean = format_cluster_names_exact(Cluster),
      negLog10Padj = -log10(p.adjust)
    )

  if (nrow(df2) == 0) {
    cat("No enrichment result in:", input_csv, "\n")
    return(NULL)
  }

  top_pathways <- df2 %>%
    group_by(Description) %>%
    summarise(best_padj = min(p.adjust, na.rm = TRUE), .groups = "drop") %>%
    arrange(best_padj) %>%
    slice_head(n = top_n_pathways) %>%
    pull(Description)

  df_plot <- df2 %>%
    filter(Description %in% top_pathways) %>%
    select(Cluster_clean, Description, negLog10Padj)

  df_plot_complete <- expand.grid(
    Cluster_clean = cluster_order_clean,
    Description = top_pathways,
    stringsAsFactors = FALSE
  ) %>%
    left_join(df_plot, by = c("Cluster_clean", "Description")) %>%
    mutate(
      negLog10Padj = ifelse(is.na(negLog10Padj), 0, negLog10Padj)
    )

  pathway_order <- df_plot_complete %>%
    group_by(Description) %>%
    summarise(best_score = max(negLog10Padj, na.rm = TRUE), .groups = "drop") %>%
    arrange(best_score) %>%
    pull(Description)

  df_plot_complete <- df_plot_complete %>%
    mutate(
      Cluster_clean = factor(Cluster_clean, levels = cluster_order_clean),
      Description = factor(Description, levels = pathway_order)
    )

  p <- ggplot(df_plot_complete, aes(x = Cluster_clean, y = Description, fill = negLog10Padj)) +
    geom_tile(color = "white", linewidth = 0.3) +
    scale_fill_gradient(
      low = "white",
      high = "firebrick",
      name = "-log10\np.adjust"
    ) +
    labs(
      title = plot_title,
      x = "Group",
      y = "Pathway"
    ) +
    theme_bw(base_size = 11) +
    theme(
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      axis.text.x = element_text(size = 8.5, lineheight = 0.85),
      axis.text.y = element_text(size = 8),
      axis.title.x = element_text(size = 11),
      axis.title.y = element_text(size = 11),
      axis.ticks = element_line(color = "black"),
      axis.ticks.length = grid::unit(0.15, "cm"),
      panel.border = element_rect(color = "black", fill = NA),
      panel.grid = element_blank()
    )

  pdf(output_pdf, width = 12, height = 10)
  print(p)
  dev.off()

  cat("Saved:", output_pdf, "\n")
}

# 9) Generate heatmaps

make_enrichment_heatmap(
  input_csv = "COMPARE_CELL_MARKER_3_TREM2_GROUPS.csv",
  output_pdf = "HEATMAP_CELL_MARKER_3_TREM2_GROUPS.pdf",
  plot_title = "Pathway-level CellMarker enrichment across genotypes",
  top_n_pathways = 30
)

make_enrichment_heatmap(
  input_csv = "COMPARE_GO_BP_3_TREM2_GROUPS.csv",
  output_pdf = "HEATMAP_GO_BP_3_TREM2_GROUPS.pdf",
  plot_title = "Pathway-level GO biological process enrichment across genotypes",
  top_n_pathways = 30
)

make_enrichment_heatmap(
  input_csv = "COMPARE_KEGG_3_TREM2_GROUPS.csv",
  output_pdf = "HEATMAP_KEGG_3_TREM2_GROUPS.pdf",
  plot_title = "Pathway-level KEGG enrichment across genotypes",
  top_n_pathways = 30
)

make_enrichment_heatmap(
  input_csv = "COMPARE_REACTOME_3_TREM2_GROUPS.csv",
  output_pdf = "HEATMAP_REACTOME_3_TREM2_GROUPS.pdf",
  plot_title = "Pathway-level CellMarker enrichment across genotypes",
  top_n_pathways = 30
)

# 10) MODULE-LEVEL HEATMAPS
# Generate module-level heatmaps for:
#   1) CellMarker
#   2) GO Biological Process
#   3) KEGG
#   4) Reactome

assign_module <- function(desc) {

  desc <- tolower(desc)

  case_when(
    str_detect(desc, "microglia|macrophage|myeloid|monocyte|dendritic|neutrophil|granulocyte") ~
      "Myeloid / microglia",

    str_detect(desc, "t cell|b cell|lymphocyte|nk cell|plasma cell|adaptive immune|antigen") ~
      "Lymphoid / adaptive immune",

    str_detect(desc, "immune|cytokine|chemokine|inflammatory|inflammation|interferon|complement|toll-like|nod-like|tnf|nf-kappa|nfkb|jak-stat") ~
      "Immune / inflammatory signaling",

    str_detect(desc, "chemotaxis|migration|locomotion|adhesion|focal adhesion|integrin|actin cytoskeleton|cytoskeleton") ~
      "Migration / adhesion / cytoskeleton",

    str_detect(desc, "apoptotic|apoptosis|cell death|neuron death|oxidative stress|stress|reactive oxygen") ~
      "Cell death / stress response",

    str_detect(desc, "extracellular matrix|ecm|collagen|matrix|angiogenesis|vasculature|blood vessel|endothelial|atherosclerosis") ~
      "ECM / vascular remodeling",

    str_detect(desc, "synapse|neuron|axon|dendrite|neurogenesis|gliogenesis|myelination|glial|astrocyte|oligodendrocyte") ~
      "Neuronal / glial biology",

    str_detect(desc, "metabolic|metabolism|lipid|cholesterol|fatty acid|mitochondrial|oxidative phosphorylation|diabetic") ~
      "Metabolism / mitochondrial function",

    str_detect(desc, "cancer|carcinoma|papillomavirus|tumor|tumour|leukemia|melanoma") ~
      "Cancer-related pathways",

    str_detect(desc, "cardiomyopathy|cardiac|heart|muscle") ~
      "Muscle / cardiac pathways",

    TRUE ~ "Other"
  )
}

module_order <- c(
  "Myeloid / microglia",
  "Lymphoid / adaptive immune",
  "Immune / inflammatory signaling",
  "Migration / adhesion / cytoskeleton",
  "Cell death / stress response",
  "ECM / vascular remodeling",
  "Neuronal / glial biology",
  "Metabolism / mitochondrial function",
  "Cancer-related pathways",
  "Muscle / cardiac pathways",
  "Other"
)

make_module_heatmap <- function(input_csv,
                                output_pdf,
                                output_csv,
                                plot_title) {

  df <- read.csv(input_csv, check.names = FALSE)

  required_cols <- c("Cluster", "Description", "p.adjust")
  missing_cols <- setdiff(required_cols, names(df))
  if (length(missing_cols) > 0) {
    stop(paste("Missing columns:", paste(missing_cols, collapse = ", ")))
  }

  df2 <- df %>%
    filter(!is.na(Cluster), !is.na(Description), !is.na(p.adjust)) %>%
    mutate(
      Cluster_clean = format_cluster_names_exact(Cluster),
      Module = assign_module(Description),
      negLog10 = -log10(p.adjust)
    )

  if (nrow(df2) == 0) {
    cat("No enrichment result in:", input_csv, "\n")
    return(NULL)
  }

  df_module <- df2 %>%
    group_by(Cluster_clean, Module) %>%
    summarise(
      score = mean(negLog10, na.rm = TRUE),
      n_pathways = n(),
      .groups = "drop"
    )

  df_module_complete <- expand.grid(
    Cluster_clean = cluster_order_clean,
    Module = module_order,
    stringsAsFactors = FALSE
  ) %>%
    left_join(df_module, by = c("Cluster_clean", "Module")) %>%
    mutate(
      score = ifelse(is.na(score), 0, score),
      n_pathways = ifelse(is.na(n_pathways), 0, n_pathways),
      Cluster_clean = factor(Cluster_clean, levels = cluster_order_clean),
      Module = factor(Module, levels = rev(module_order))
    )

  p_module <- ggplot(df_module_complete, aes(x = Cluster_clean, y = Module, fill = score)) +
    geom_tile(color = "white", linewidth = 0.4) +
    geom_text(aes(label = ifelse(n_pathways > 0, n_pathways, "")), size = 3) +
    scale_fill_gradient(
      low = "white",
      high = "firebrick",
      name = "Mean\n-log10(p.adjust)"
    ) +
    labs(
      title = plot_title,
      x = "Group",
      y = "Module"
    ) +
    theme_bw(base_size = 12) +
    theme(
      axis.text.x = element_text(size = 8.5, lineheight = 0.85),
      axis.text.y = element_text(size = 10),
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      axis.ticks = element_line(color = "black"),
      axis.ticks.length = grid::unit(0.15, "cm"),
      panel.border = element_rect(color = "black", fill = NA),
      panel.grid = element_blank()
    )

  pdf(output_pdf, width = 12, height = 8)
  print(p_module)
  dev.off()

  write.csv(
    df_module_complete,
    output_csv,
    row.names = FALSE
  )

  cat("Saved:", output_pdf, "\n")
  cat("Saved:", output_csv, "\n")
}

# 10.1 CellMarker module-level heatmap

make_module_heatmap(
  input_csv = "COMPARE_CELL_MARKER_3_TREM2_GROUPS.csv",
  output_pdf = "HEATMAP_MODULE_LEVEL_CELL_MARKER_3_TREM2_GROUPS.pdf",
  output_csv = "MODULE_LEVEL_CELL_MARKER_SUMMARY_3_TREM2_GROUPS.csv",
  plot_title = "Module-level CellMarker enrichment across genotypes"
)

# 10.2 GO BP module-level heatmap

make_module_heatmap(
  input_csv = "COMPARE_GO_BP_3_TREM2_GROUPS.csv",
  output_pdf = "HEATMAP_MODULE_LEVEL_GO_BP_3_TREM2_GROUPS.pdf",
  output_csv = "MODULE_LEVEL_GO_BP_SUMMARY_3_TREM2_GROUPS.csv",
  plot_title = "Module-level biological process enrichment across genotypes"
)

# 10.3 KEGG module-level heatmap

make_module_heatmap(
  input_csv = "COMPARE_KEGG_3_TREM2_GROUPS.csv",
  output_pdf = "HEATMAP_MODULE_LEVEL_KEGG_3_TREM2_GROUPS.pdf",
  output_csv = "MODULE_LEVEL_KEGG_SUMMARY_3_TREM2_GROUPS.csv",
  plot_title = "Module-level KEGG enrichment across genotypes"
)

# 10.4 Reactome module-level heatmap

make_module_heatmap(
  input_csv = "COMPARE_REACTOME_3_TREM2_GROUPS.csv",
  output_pdf = "HEATMAP_MODULE_LEVEL_REACTOME_3_TREM2_GROUPS.pdf",
  output_csv = "MODULE_LEVEL_REACTOME_SUMMARY_3_TREM2_GROUPS.csv",
  plot_title = "Module-level Reactome enrichment across genotypes"
)

cat("\nALL DONE: pathway-level and module-level enrichment analyses completed for CellMarker, GO BP, KEGG, and Reactome.\n")

# FINAL PIPELINE:
# 4 Trem2 UP gene groups only
# HEATMAP-ONLY VERSION
#
# Analyses:
#   1) Cell marker enrichment
#   2) GO BP enrichment
#   3) KEGG enrichment
#   4) Reactome enrichment
#   5) pathway-level heatmaps
#   6) module-level heatmaps

library(clusterProfiler)
library(org.Mm.eg.db)
library(ReactomePA)
library(dplyr)
library(ggplot2)
library(tidyr)
library(stringr)
library(readxl)
library(grid)

select  <- dplyr::select
filter  <- dplyr::filter
mutate  <- dplyr::mutate
arrange <- dplyr::arrange

# 1) Input files: 4 Trem2 groups
files <- list(
  Trem2WT_LSD1KO_vs_Trem2WT_LSD1WT =
    "UP_Trem2_WT_LSD1_Del_vs_Trem2_WT_LSD1_WT.csv",

  Trem2KO_LSD1KO_vs_Trem2WT_LSD1KO =
    "UP_Trem2KO_LSD1_Del_vs_Trem2_WT_LSD1_Del.csv",

  Trem2KO_LSD1KO_vs_Trem2WT_LSD1WT =
    "UP_Trem2KO_LSD1_Del_vs_Trem2_WT_LSD1_WT.csv",

  Trem2KO_LSD1WT_vs_Trem2WT_LSD1WT =
    "UP_Trem2KO_LSD1_WT_vs_Trem2_WT_LSD1_WT.csv"
)

cell_marker_file <- "Cell_marker_Mouse.xlsx"

clean_genes <- function(x) {
  x <- as.character(x)
  x <- trimws(x)
  x <- x[!is.na(x) & x != ""]
  unique(x)
}

convert_entrez <- function(genes) {
  df <- suppressWarnings(
    bitr(
      genes,
      fromType = "SYMBOL",
      toType = "ENTREZID",
      OrgDb = org.Mm.eg.db
    )
  )

  df <- df %>%
    distinct(ENTREZID, .keep_all = TRUE)

  return(df$ENTREZID)
}

safe_neglog10 <- function(p) {
  p <- as.numeric(p)
  p[p == 0] <- .Machine$double.xmin
  -log10(p)
}

format_cluster_names_exact <- function(x) {
  map <- c(
    "Trem2WT_LSD1KO_vs_Trem2WT_LSD1WT" =
      "Trem2WT_LSD1KO\nvs\nTrem2WT_LSD1WT",

    "Trem2KO_LSD1KO_vs_Trem2WT_LSD1KO" =
      "Trem2KO_LSD1KO\nvs\nTrem2WT_LSD1KO",

    "Trem2KO_LSD1KO_vs_Trem2WT_LSD1WT" =
      "Trem2KO_LSD1KO\nvs\nTrem2WT_LSD1WT",

    "Trem2KO_LSD1WT_vs_Trem2WT_LSD1WT" =
      "Trem2KO_LSD1WT\nvs\nTrem2WT_LSD1WT"
  )

  out <- unname(map[x])
  out[is.na(out)] <- x[is.na(out)]
  return(out)
}

cluster_order_raw <- names(files)
cluster_order_clean <- format_cluster_names_exact(cluster_order_raw)

# 3) Build gene lists
gene_list_symbol <- list()
gene_list_entrez <- list()

cat("\n========== BUILD 4-GROUP TREM2 GENE LIST ==========\n")

for (nm in names(files)) {

  cat("\nProcessing:", nm, "\n")

  df <- read.csv(files[[nm]], check.names = FALSE)

  if (!("Gene_name" %in% names(df))) {
    cat("Missing Gene_name column -> skip\n")
    next
  }

  genes <- clean_genes(df$Gene_name)
  cat("Total SYMBOL genes:", length(genes), "\n")

  if (length(genes) < 5) {
    cat("Too few genes -> skip enrichment\n")
    next
  }

  entrez <- convert_entrez(genes)

  cat("Mapped ENTREZ genes:", length(entrez), "\n")
  cat("Mapping rate:", round(length(entrez) / length(genes) * 100, 2), "%\n")

  if (length(entrez) < 5) {
    cat("Too few mapped genes -> skip enrichment\n")
    next
  }

  gene_list_symbol[[nm]] <- genes
  gene_list_entrez[[nm]] <- entrez
}

cat("\nFinal groups used for enrichment:\n")
print(names(gene_list_entrez))

# 4) CELL MARKER ENRICHMENT

cat("\n========== RUN CELL MARKER ENRICHMENT ==========\n")

cell_marker_raw <- readxl::read_excel(cell_marker_file, sheet = "mouse")

cell_marker_term2gene <- cell_marker_raw %>%
  dplyr::select(cell_name, tissue_type, Symbol) %>%
  dplyr::mutate(
    cell_name   = as.character(cell_name),
    tissue_type = as.character(tissue_type),
    Symbol      = trimws(as.character(Symbol)),
    Cell_marker = paste(cell_name, tissue_type, "Mouse", sep = " "),
    Gene        = Symbol
  ) %>%
  dplyr::filter(!is.na(Cell_marker), Cell_marker != "") %>%
  dplyr::filter(!is.na(Gene), Gene != "") %>%
  dplyr::distinct(Cell_marker, Gene)

cell_compare <- compareCluster(
  geneCluster = gene_list_symbol,
  fun = "enricher",
  TERM2GENE = cell_marker_term2gene,
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.2
)

write.csv(
  as.data.frame(cell_compare),
  "COMPARE_CELL_MARKER_4_TREM2_GROUPS_HEATMAP.csv",
  row.names = FALSE
)

# 5) GO BP ENRICHMENT

cat("\n========== RUN GO BP COMPARE ==========\n")

ego_compare <- compareCluster(
  geneCluster = gene_list_entrez,
  fun = "enrichGO",
  OrgDb = org.Mm.eg.db,
  keyType = "ENTREZID",
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.2
)

write.csv(
  as.data.frame(ego_compare),
  "COMPARE_GO_BP_4_TREM2_GROUPS_HEATMAP.csv",
  row.names = FALSE
)

# 6) KEGG ENRICHMENT

cat("\n========== RUN KEGG COMPARE ==========\n")

kegg_compare <- compareCluster(
  geneCluster = gene_list_entrez,
  fun = "enrichKEGG",
  organism = "mmu",
  pvalueCutoff = 0.05
)

write.csv(
  as.data.frame(kegg_compare),
  "COMPARE_KEGG_4_TREM2_GROUPS_HEATMAP.csv",
  row.names = FALSE
)

# 7) REACTOME ENRICHMENT

cat("\n========== RUN REACTOME COMPARE ==========\n")

react_compare <- compareCluster(
  geneCluster = gene_list_entrez,
  fun = "enrichPathway",
  organism = "mouse",
  pvalueCutoff = 0.05
)

write.csv(
  as.data.frame(react_compare),
  "COMPARE_REACTOME_4_TREM2_GROUPS_HEATMAP.csv",
  row.names = FALSE
)

# 8) Module classifier

assign_module <- function(desc) {

  desc <- tolower(desc)

  case_when(
    str_detect(desc, "microglia|macrophage|myeloid|monocyte|dendritic|neutrophil|granulocyte") ~
      "Myeloid / microglia",

    str_detect(desc, "t cell|b cell|lymphocyte|nk cell|plasma cell|adaptive immune|antigen") ~
      "Lymphoid / adaptive immune",

    str_detect(desc, "immune|cytokine|chemokine|inflammatory|inflammation|interferon|complement|toll-like|nod-like|tnf|nf-kappa|nfkb|jak-stat") ~
      "Immune / inflammatory signaling",

    str_detect(desc, "chemotaxis|migration|locomotion|adhesion|focal adhesion|integrin|actin cytoskeleton|cytoskeleton") ~
      "Migration / adhesion / cytoskeleton",

    str_detect(desc, "apoptotic|apoptosis|cell death|neuron death|oxidative stress|stress|reactive oxygen") ~
      "Cell death / stress response",

    str_detect(desc, "extracellular matrix|ecm|collagen|matrix|angiogenesis|vasculature|blood vessel|endothelial|atherosclerosis") ~
      "ECM / vascular remodeling",

    str_detect(desc, "synapse|neuron|axon|dendrite|neurogenesis|gliogenesis|myelination|glial|astrocyte|oligodendrocyte") ~
      "Neuronal / glial biology",

    str_detect(desc, "metabolic|metabolism|lipid|cholesterol|fatty acid|mitochondrial|oxidative phosphorylation|diabetic") ~
      "Metabolism / mitochondrial function",

    str_detect(desc, "cancer|carcinoma|papillomavirus|tumor|tumour|leukemia|melanoma") ~
      "Cancer-related pathways",

    str_detect(desc, "cardiomyopathy|cardiac|heart|muscle") ~
      "Muscle / cardiac pathways",

    TRUE ~ "Other"
  )
}

module_order <- c(
  "Myeloid / microglia",
  "Lymphoid / adaptive immune",
  "Immune / inflammatory signaling",
  "Migration / adhesion / cytoskeleton",
  "Cell death / stress response",
  "ECM / vascular remodeling",
  "Neuronal / glial biology",
  "Metabolism / mitochondrial function",
  "Cancer-related pathways",
  "Muscle / cardiac pathways",
  "Other"
)

# 9) pathway-level heatmap function

make_enrichment_heatmap <- function(input_csv,
                                    output_pdf,
                                    top_n_pathways = 30) {

  df <- read.csv(input_csv, check.names = FALSE)

  required_cols <- c("Cluster", "Description", "p.adjust")
  missing_cols <- setdiff(required_cols, names(df))
  if (length(missing_cols) > 0) {
    stop(paste("Missing columns:", paste(missing_cols, collapse = ", ")))
  }

  df2 <- df %>%
    filter(!is.na(Cluster), !is.na(Description), !is.na(p.adjust)) %>%
    mutate(
      Cluster_clean = format_cluster_names_exact(Cluster),
      negLog10Padj = safe_neglog10(p.adjust)
    )

  if (nrow(df2) == 0) {
    cat("No enrichment result in:", input_csv, "\n")
    return(NULL)
  }

  top_pathways <- df2 %>%
    group_by(Description) %>%
    summarise(best_padj = min(p.adjust, na.rm = TRUE), .groups = "drop") %>%
    arrange(best_padj) %>%
    slice_head(n = top_n_pathways) %>%
    pull(Description)

  df_plot <- df2 %>%
    filter(Description %in% top_pathways) %>%
    select(Cluster_clean, Description, negLog10Padj)

  df_plot_complete <- expand.grid(
    Cluster_clean = cluster_order_clean,
    Description = top_pathways,
    stringsAsFactors = FALSE
  ) %>%
    left_join(df_plot, by = c("Cluster_clean", "Description")) %>%
    mutate(
      negLog10Padj = ifelse(is.na(negLog10Padj), 0, negLog10Padj)
    )

  pathway_order <- df_plot_complete %>%
    group_by(Description) %>%
    summarise(best_score = max(negLog10Padj, na.rm = TRUE), .groups = "drop") %>%
    arrange(best_score) %>%
    pull(Description)

  df_plot_complete <- df_plot_complete %>%
    mutate(
      Cluster_clean = factor(Cluster_clean, levels = cluster_order_clean),
      Description = factor(Description, levels = pathway_order)
    )

  final_title <- case_when(
    output_pdf == "HEATMAP_CELL_MARKER_4_TREM2_GROUPS_HEATMAP.pdf" ~ "Pathway-level CellMarker enrichment across genotypes",
    output_pdf == "HEATMAP_GO_BP_4_TREM2_GROUPS_HEATMAP.pdf" ~ "Pathway-level GO biological process enrichment across genotypes",
    output_pdf == "HEATMAP_KEGG_4_TREM2_GROUPS_HEATMAP.pdf" ~ "Pathway-level KEGG enrichment across genotypes",
    output_pdf == "HEATMAP_REACTOME_4_TREM2_GROUPS_HEATMAP.pdf" ~ "Pathway-level Reactome enrichment across genotypes",
    TRUE ~ "Pathway-level enrichment across genotypes"
  )

  p <- ggplot(df_plot_complete, aes(x = Cluster_clean, y = Description, fill = negLog10Padj)) +
    geom_tile(color = "white", linewidth = 0.3) +
    scale_fill_gradient(
      low = "white",
      high = "firebrick",
      name = "-log10\np.adjust"
    ) +
    labs(
      title = final_title,
      x = "Group",
      y = "Pathway"
    ) +
    theme_bw(base_size = 11) +
    theme(
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      axis.text.x = element_text(size = 8.5, lineheight = 0.85),
      axis.text.y = element_text(size = 8),
      axis.title.x = element_text(size = 11),
      axis.title.y = element_text(size = 11),
      axis.ticks = element_line(color = "black"),
      axis.ticks.length = grid::unit(0.15, "cm"),
      panel.border = element_rect(color = "black", fill = NA),
      panel.grid = element_blank()
    )

  pdf(output_pdf, width = 12, height = 10)
  print(p)
  dev.off()

  cat("Saved:", output_pdf, "\n")
}

# 10) module-level heatmap function

make_module_heatmap <- function(input_csv,
                                output_pdf,
                                output_csv) {

  df <- read.csv(input_csv, check.names = FALSE)

  required_cols <- c("Cluster", "Description", "p.adjust")
  missing_cols <- setdiff(required_cols, names(df))
  if (length(missing_cols) > 0) {
    stop(paste("Missing columns:", paste(missing_cols, collapse = ", ")))
  }

  df2 <- df %>%
    filter(!is.na(Cluster), !is.na(Description), !is.na(p.adjust)) %>%
    mutate(
      Cluster_clean = format_cluster_names_exact(Cluster),
      Module = assign_module(Description),
      negLog10 = safe_neglog10(p.adjust)
    )

  if (nrow(df2) == 0) {
    cat("No enrichment result in:", input_csv, "\n")
    return(NULL)
  }

  df_module <- df2 %>%
    group_by(Cluster_clean, Module) %>%
    summarise(
      score = mean(negLog10, na.rm = TRUE),
      n_pathways = n(),
      .groups = "drop"
    )

  df_module_complete <- expand.grid(
    Cluster_clean = cluster_order_clean,
    Module = module_order,
    stringsAsFactors = FALSE
  ) %>%
    left_join(df_module, by = c("Cluster_clean", "Module")) %>%
    mutate(
      score = ifelse(is.na(score), 0, score),
      n_pathways = ifelse(is.na(n_pathways), 0, n_pathways),
      Cluster_clean = factor(Cluster_clean, levels = cluster_order_clean),
      Module = factor(Module, levels = rev(module_order))
    )

  final_title <- case_when(
    output_pdf == "HEATMAP_MODULE_CELL_MARKER_4_TREM2_GROUPS_HEATMAP.pdf" ~ "Module-level CellMarker enrichment across genotypes",
    output_pdf == "HEATMAP_MODULE_GO_BP_4_TREM2_GROUPS_HEATMAP.pdf" ~ "Module-level GO biological process enrichment across genotypes",
    output_pdf == "HEATMAP_MODULE_KEGG_4_TREM2_GROUPS_HEATMAP.pdf" ~ "Module-level KEGG enrichment across genotypes",
    output_pdf == "HEATMAP_MODULE_REACTOME_4_TREM2_GROUPS_HEATMAP.pdf" ~ "Module-level Reactome enrichment across genotypes",
    TRUE ~ "Module-level enrichment across genotypes"
  )

  p_module <- ggplot(df_module_complete, aes(x = Cluster_clean, y = Module, fill = score)) +
    geom_tile(color = "white", linewidth = 0.4) +
    geom_text(aes(label = ifelse(n_pathways > 0, n_pathways, "")), size = 3) +
    scale_fill_gradient(
      low = "white",
      high = "firebrick",
      name = "Mean\n-log10(p.adjust)"
    ) +
    labs(
      title = final_title,
      x = "Group",
      y = "Module"
    ) +
    theme_bw(base_size = 12) +
    theme(
      axis.text.x = element_text(size = 8.5, lineheight = 0.85),
      axis.text.y = element_text(size = 10),
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      axis.ticks = element_line(color = "black"),
      axis.ticks.length = grid::unit(0.15, "cm"),
      panel.border = element_rect(color = "black", fill = NA),
      panel.grid = element_blank()
    )

  pdf(output_pdf, width = 12, height = 8)
  print(p_module)
  dev.off()

  write.csv(
    df_module_complete,
    output_csv,
    row.names = FALSE
  )

  cat("Saved:", output_pdf, "\n")
  cat("Saved:", output_csv, "\n")
}

# 11) Generate pathway-level heatmaps

make_enrichment_heatmap(
  input_csv = "COMPARE_CELL_MARKER_4_TREM2_GROUPS_HEATMAP.csv",
  output_pdf = "HEATMAP_CELL_MARKER_4_TREM2_GROUPS_HEATMAP.pdf",
  top_n_pathways = 30
)

make_enrichment_heatmap(
  input_csv = "COMPARE_GO_BP_4_TREM2_GROUPS_HEATMAP.csv",
  output_pdf = "HEATMAP_GO_BP_4_TREM2_GROUPS_HEATMAP.pdf",
  top_n_pathways = 30
)

make_enrichment_heatmap(
  input_csv = "COMPARE_KEGG_4_TREM2_GROUPS_HEATMAP.csv",
  output_pdf = "HEATMAP_KEGG_4_TREM2_GROUPS_HEATMAP.pdf",
  top_n_pathways = 30
)

make_enrichment_heatmap(
  input_csv = "COMPARE_REACTOME_4_TREM2_GROUPS_HEATMAP.csv",
  output_pdf = "HEATMAP_REACTOME_4_TREM2_GROUPS_HEATMAP.pdf",
  top_n_pathways = 30
)

# 12) Generate module-level heatmaps

make_module_heatmap(
  input_csv = "COMPARE_CELL_MARKER_4_TREM2_GROUPS_HEATMAP.csv",
  output_pdf = "HEATMAP_MODULE_CELL_MARKER_4_TREM2_GROUPS_HEATMAP.pdf",
  output_csv = "MODULE_LEVEL_CELL_MARKER_SUMMARY_4_TREM2_GROUPS_HEATMAP.csv"
)

make_module_heatmap(
  input_csv = "COMPARE_GO_BP_4_TREM2_GROUPS_HEATMAP.csv",
  output_pdf = "HEATMAP_MODULE_GO_BP_4_TREM2_GROUPS_HEATMAP.pdf",
  output_csv = "MODULE_LEVEL_GO_BP_SUMMARY_4_TREM2_GROUPS_HEATMAP.csv"
)

make_module_heatmap(
  input_csv = "COMPARE_KEGG_4_TREM2_GROUPS_HEATMAP.csv",
  output_pdf = "HEATMAP_MODULE_KEGG_4_TREM2_GROUPS_HEATMAP.pdf",
  output_csv = "MODULE_LEVEL_KEGG_SUMMARY_4_TREM2_GROUPS_HEATMAP.csv"
)

make_module_heatmap(
  input_csv = "COMPARE_REACTOME_4_TREM2_GROUPS_HEATMAP.csv",
  output_pdf = "HEATMAP_MODULE_REACTOME_4_TREM2_GROUPS_HEATMAP.pdf",
  output_csv = "MODULE_LEVEL_REACTOME_SUMMARY_4_TREM2_GROUPS_HEATMAP.csv"
)

cat("\nALL DONE: pathway-level and module-level Trem2 enrichment heatmaps completed for 4 groups.\n")
