# 12_lsd1_timecourse_corrected_terminal_enrichment_no_michael.R
#
# Purpose: Corrected LSD1 time-course enrichment analysis from no-Michael volcano outputs.
# Inputs are expected under data/processed/ or data/external/ relative to this repository.
# Outputs are written under results/ or script-defined subfolders.

# CORRECTED enrichment from corrected volcano output
# SINGLE multi-page PDF version
#
# Input:
#   Corrected volcano output files created from:
#   LSD1 terminal + Trem2WT-background terminal, NO Michael
#
#   - Does NOT create many separate PDFs.
#   - Creates ONE multi-page PDF containing:
#       UP CellMarker pathway heatmap
#       UP CellMarker module heatmap
#       UP GO BP pathway heatmap
#       UP GO BP module heatmap
#       UP KEGG pathway heatmap
#       UP KEGG module heatmap
#       UP Reactome pathway heatmap
#       UP Reactome module heatmap
#       DOWN CellMarker pathway heatmap
#       DOWN CellMarker module heatmap
#       DOWN GO BP pathway heatmap
#       DOWN GO BP module heatmap
#       DOWN KEGG pathway heatmap
#       DOWN KEGG module heatmap
#       DOWN Reactome pathway heatmap
#       DOWN Reactome module heatmap
#
# Still writes CSV tables/check files for troubleshooting.


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
  library(clusterProfiler)
  library(org.Mm.eg.db)
  library(ReactomePA)
  library(dplyr)
  library(ggplot2)
  library(readxl)
  library(tidyr)
  library(stringr)
})

select  <- dplyr::select
filter  <- dplyr::filter
mutate  <- dplyr::mutate
arrange <- dplyr::arrange

katz_root <- repo_root
project_dir <- file.path(katz_root, "timepoint_lsd1")

# Corrected volcano output folder from the corrected no-Michael volcano script.
base_dir <- file.path(
  project_dir,
  "正确combine terminal",
  "CORRECTED_volcano_LSD1terminal_plus_Trem2WTterminal_no_Michael_GALAXY_style_FC"
)

if (!dir.exists(base_dir)) {
  stop(
    "Input folder does not exist:\n",
    base_dir,
    "\n\nRun the corrected volcano script first, or check the folder name."
  )
}

# Support files
cellmarker_file <- file.path(project_dir, "Cell_marker_Mouse.xlsx")
allc_file       <- file.path(project_dir, "NIHMS472534-supplement-02.csv")
mod_file        <- file.path(project_dir, "modules.csv")

required_support_files <- c(cellmarker_file, allc_file, mod_file)
missing_support_files <- required_support_files[!file.exists(required_support_files)]

if (length(missing_support_files) > 0) {
  stop("Missing support files:\n", paste(missing_support_files, collapse = "\n"))
}

# One output folder, one output PDF.
out_dir <- file.path(base_dir, "enrichment_SINGLE_PDF_CORRECTED_GALAXY_STYLE_FC")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

single_pdf <- file.path(
  out_dir,
  "CORRECTED_LSD1terminal_plus_Trem2WTterminal_NO_MICHAEL_enrichment_UP_DOWN_SINGLE_MULTIPAGE.pdf"
)

# 1) Mouse symbol normalization
background_entrez <- keys(org.Mm.eg.db, keytype = "ENTREZID")
official_symbols  <- keys(org.Mm.eg.db, keytype = "SYMBOL")

symbol_case_lookup <- data.frame(
  SYMBOL_UPPER = toupper(official_symbols),
  Official_SYMBOL = official_symbols,
  stringsAsFactors = FALSE
) %>%
  distinct(SYMBOL_UPPER, .keep_all = TRUE)

clean_genes <- function(x) {
  x <- as.character(x)
  x <- trimws(x)
  x <- x[!is.na(x) & x != ""]
  unique(x)
}

normalize_mouse_symbols <- function(x) {
  x <- clean_genes(x)
  if (length(x) == 0) return(character(0))

  x_upper <- toupper(x)
  matched <- symbol_case_lookup$Official_SYMBOL[
    match(x_upper, symbol_case_lookup$SYMBOL_UPPER)
  ]

  matched <- matched[!is.na(matched) & matched != ""]
  unique(matched)
}

map_genes <- function(genes) {
  genes <- normalize_mouse_symbols(genes)

  if (length(genes) == 0) {
    return(data.frame(SYMBOL = character(0), ENTREZID = character(0)))
  }

  gene_df <- suppressWarnings(
    bitr(
      genes,
      fromType = "SYMBOL",
      toType = "ENTREZID",
      OrgDb = org.Mm.eg.db
    )
  )

  if (is.null(gene_df) || nrow(gene_df) == 0) {
    return(data.frame(SYMBOL = character(0), ENTREZID = character(0)))
  }

  gene_df %>%
    distinct(ENTREZID, .keep_all = TRUE)
}

safe_neglog10 <- function(p) {
  p <- as.numeric(p)
  p[p == 0] <- .Machine$double.xmin
  -log10(p)
}

# 2) Module classification
classify_module <- function(term) {
  term_lower <- tolower(term)

  case_when(
    str_detect(term_lower, "cytokine|chemokine|interleukin|interferon|tnf|nf-kappa|nfkb|toll|immune|inflammatory|inflammation|complement|leukocyte|lymphocyte|macrophage|microglia|myeloid|phagocyt|antigen|mhc|il-17|jak-stat") ~
      "Immune / inflammatory activation",

    str_detect(term_lower, "migration|chemotaxis|adhesion|integrin|focal adhesion|extravasation|motility") ~
      "Cell migration / chemotaxis",

    str_detect(term_lower, "mapk|pi3k|akt|ras|rap1|camp|calcium|signaling|receptor|cascade|transduction") ~
      "Signaling cascade",

    str_detect(term_lower, "apoptosis|cell death|necrosis|ferroptosis|stress|oxidative|hypoxia|p53|dna damage") ~
      "Cell death / stress",

    str_detect(term_lower, "extracellular matrix|ecm|collagen|vascular|angiogenesis|endothelial|basement membrane|matrix|atherosclerosis|shear stress") ~
      "ECM / vascular remodeling",

    str_detect(term_lower, "neuron|neuronal|synapse|axon|dendrite|glial|astrocyte|oligodendrocyte|myelin|brain|cerebrum|cerebellum") ~
      "Neuronal / glial biology",

    TRUE ~ "Other"
  )
}

make_blank_page <- function(title_text) {
  ggplot() +
    annotate("text", x = 0, y = 0, label = title_text, size = 6) +
    theme_void()
}

make_pathway_heatmap <- function(enrich_df,
                                 analysis_name,
                                 direction_label,
                                 cluster_levels,
                                 top_n = 30,
                                 high_color = "firebrick") {

  if (is.null(enrich_df) || nrow(enrich_df) == 0) {
    return(make_blank_page(
      paste0("No significant ", analysis_name, " enrichment for ", direction_label, " genes")
    ))
  }

  df <- enrich_df %>%
    mutate(
      score = safe_neglog10(p.adjust),
      Description = as.character(Description),
      Cluster = as.character(Cluster)
    )

  top_terms <- df %>%
    group_by(Description) %>%
    summarise(max_score = max(score, na.rm = TRUE), .groups = "drop") %>%
    arrange(desc(max_score)) %>%
    slice_head(n = top_n) %>%
    pull(Description)

  if (length(top_terms) == 0) {
    return(make_blank_page(
      paste0("No plottable ", analysis_name, " enrichment for ", direction_label, " genes")
    ))
  }

  plot_df <- df %>%
    filter(Description %in% top_terms) %>%
    mutate(
      Description = factor(Description, levels = rev(top_terms)),
      Cluster = factor(Cluster, levels = cluster_levels)
    )

  ggplot(plot_df, aes(x = Cluster, y = Description, fill = score)) +
    geom_tile(color = "white") +
    scale_fill_gradient(
      low = "white",
      high = high_color,
      name = expression(-log[10]("p.adjust"))
    ) +
    labs(
      title = paste0(direction_label, " genes: ", analysis_name, " pathway enrichment"),
      subtitle = "Corrected combined terminal = LSD1 terminal + Trem2WT-background terminal; no Michael",
      x = "Time point",
      y = "Pathway"
    ) +
    theme_bw(base_size = 11) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5, size = 9),
      axis.text.x = element_text(angle = 35, hjust = 1),
      axis.text.y = element_text(size = 8),
      panel.grid = element_blank()
    )
}

make_module_heatmap <- function(enrich_df,
                                analysis_name,
                                direction_label,
                                cluster_levels,
                                high_color = "firebrick") {

  if (is.null(enrich_df) || nrow(enrich_df) == 0) {
    return(make_blank_page(
      paste0("No significant ", analysis_name, " module enrichment for ", direction_label, " genes")
    ))
  }

  module_df <- enrich_df %>%
    mutate(
      score = safe_neglog10(p.adjust),
      Module = classify_module(Description),
      Cluster = as.character(Cluster)
    ) %>%
    group_by(Cluster, Module) %>%
    summarise(
      mean_score = mean(score, na.rm = TRUE),
      n_terms = n(),
      .groups = "drop"
    ) %>%
    mutate(
      Cluster = factor(Cluster, levels = cluster_levels)
    )

  if (nrow(module_df) == 0) {
    return(make_blank_page(
      paste0("No plottable ", analysis_name, " module enrichment for ", direction_label, " genes")
    ))
  }

  module_levels <- module_df %>%
    group_by(Module) %>%
    summarise(max_score = max(mean_score, na.rm = TRUE), .groups = "drop") %>%
    arrange(max_score) %>%
    pull(Module)

  module_df <- module_df %>%
    mutate(Module = factor(Module, levels = module_levels))

  ggplot(module_df, aes(x = Cluster, y = Module, fill = mean_score)) +
    geom_tile(color = "white") +
    geom_text(aes(label = n_terms), size = 3) +
    scale_fill_gradient(
      low = "white",
      high = high_color,
      name = expression("Mean " * -log[10]("p.adjust"))
    ) +
    labs(
      title = paste0(direction_label, " genes: ", analysis_name, " module summary"),
      subtitle = "Numbers inside tiles = number of enriched terms assigned to that module",
      x = "Time point",
      y = "Biological module"
    ) +
    theme_bw(base_size = 11) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5, size = 9),
      axis.text.x = element_text(angle = 35, hjust = 1),
      panel.grid = element_blank()
    )
}

# 4) CellMarker TERM2GENE
cellmarker <- readxl::read_excel(cellmarker_file, sheet = "mouse")

cellmarker_term2gene <- cellmarker %>%
  select(cell_name, tissue_type, Symbol) %>%
  mutate(
    cell_name = as.character(cell_name),
    tissue_type = as.character(tissue_type),
    Symbol_raw = as.character(Symbol),
    SYMBOL_UPPER = toupper(trimws(Symbol_raw)),
    Symbol = symbol_case_lookup$Official_SYMBOL[
      match(SYMBOL_UPPER, symbol_case_lookup$SYMBOL_UPPER)
    ],
    Term = paste(cell_name, tissue_type, "Mouse", sep = " ")
  ) %>%
  filter(!is.na(Term), !is.na(Symbol), Term != "", Symbol != "") %>%
  distinct(Term, Symbol)

# 5) Input files
get_input_files <- function(direction_label) {
  list(
    `3weeks` = file.path(base_dir, paste0(direction_label, "_genes_3weeks_CORRECTED_NO_MICHAEL.csv")),
    `4weeks` = file.path(base_dir, paste0(direction_label, "_genes_4weeks_CORRECTED_NO_MICHAEL.csv")),
    `Early_onset` = file.path(base_dir, paste0(direction_label, "_genes_Early_onset_CORRECTED_NO_MICHAEL.csv")),
    `Terminal` = file.path(base_dir, paste0(direction_label, "_genes_Terminal_CORRECTED_LSD1terminal_plus_Trem2WTterminal_NO_MICHAEL.csv"))
  )
}

# More permissive fallback in case file names differ slightly from the corrected volcano script.
find_input_file_fallback <- function(direction_label, timepoint_key) {
  pattern_map <- list(
    `3weeks` = paste0("^", direction_label, "_genes_3.*\\.csv$"),
    `4weeks` = paste0("^", direction_label, "_genes_4.*\\.csv$"),
    `Early_onset` = paste0("^", direction_label, "_genes_Early.*\\.csv$"),
    `Terminal` = paste0("^", direction_label, "_genes_Terminal.*NO_MICHAEL.*\\.csv$")
  )

  files <- list.files(base_dir, pattern = pattern_map[[timepoint_key]], full.names = TRUE, ignore.case = TRUE)
  if (length(files) == 0 && timepoint_key == "Terminal") {
    files <- list.files(base_dir, pattern = paste0("^", direction_label, "_genes_Terminal.*\\.csv$"), full.names = TRUE, ignore.case = TRUE)
  }
  if (length(files) == 0) return(NA_character_)
  files[1]
}

resolve_input_files <- function(direction_label) {
  input_files <- get_input_files(direction_label)

  for (nm in names(input_files)) {
    if (!file.exists(input_files[[nm]])) {
      fallback <- find_input_file_fallback(direction_label, nm)
      if (!is.na(fallback) && file.exists(fallback)) {
        input_files[[nm]] <- fallback
      }
    }
  }

  input_files
}

run_enrichment_for_direction <- function(direction = c("UP", "DOWN")) {

  direction <- match.arg(direction)
  direction_label <- direction
  high_color <- ifelse(direction == "UP", "firebrick", "steelblue")

  input_files <- resolve_input_files(direction_label)

  check_input_df <- data.frame(
    Direction = direction_label,
    Timepoint = names(input_files),
    File = unlist(input_files),
    Exists = file.exists(unlist(input_files)),
    stringsAsFactors = FALSE
  )

  write.csv(
    check_input_df,
    file.path(out_dir, paste0("CHECK_input_files_", direction_label, ".csv")),
    row.names = FALSE
  )

  if (!all(check_input_df$Exists)) {
    stop(
      "Missing ", direction_label, " gene input files. See:\n",
      file.path(out_dir, paste0("CHECK_input_files_", direction_label, ".csv"))
    )
  }

  cluster_levels <- names(input_files)

  cat("\n============================================================\n")
  cat("Running enrichment for ", direction_label, " genes\n", sep = "")
  cat("Input folder: ", base_dir, "\n", sep = "")
  cat("Output folder: ", out_dir, "\n", sep = "")
  cat("============================================================\n")

  gene_list_symbol <- list()
  gene_list_entrez <- list()
  gene_count_rows <- list()

  for (nm in names(input_files)) {

    df <- read.csv(input_files[[nm]], check.names = FALSE)

    if (!("Gene_name" %in% names(df))) {
      stop("Gene_name column missing in: ", input_files[[nm]])
    }

    genes_raw <- clean_genes(df$Gene_name)
    genes <- normalize_mouse_symbols(genes_raw)

    if (length(genes) == 0) {
      gene_df <- data.frame(SYMBOL = character(0), ENTREZID = character(0))
    } else {
      gene_df <- map_genes(genes)
    }

    gene_list_symbol[[nm]] <- genes
    gene_list_entrez[[nm]] <- gene_df$ENTREZID

    gene_count_rows[[nm]] <- data.frame(
      Direction = direction_label,
      Timepoint = nm,
      InputFile = input_files[[nm]],
      RawGenes = length(genes_raw),
      OfficialMouseSymbols = length(genes),
      EntrezMapped = nrow(gene_df),
      stringsAsFactors = FALSE
    )

    write.csv(
      gene_df,
      file.path(out_dir, paste0("SYMBOL_to_ENTREZ_", nm, "_", direction_label, "_CORRECTED.csv")),
      row.names = FALSE
    )
  }

  gene_count_df <- bind_rows(gene_count_rows)

  write.csv(
    gene_count_df,
    file.path(out_dir, paste0("CHECK_gene_counts_", direction_label, "_CORRECTED_GALAXY_STYLE_FC.csv")),
    row.names = FALSE
  )

  # Keep only usable gene sets for compareCluster.
  gene_list_symbol_usable <- gene_list_symbol[lengths(gene_list_symbol) >= 5]
  gene_list_entrez_usable <- gene_list_entrez[lengths(gene_list_entrez) >= 10]

  results <- list()

  # CellMarker
  if (length(gene_list_symbol_usable) > 0) {
    cellmarker_compare <- compareCluster(
      geneCluster   = gene_list_symbol_usable,
      fun           = "enricher",
      TERM2GENE     = cellmarker_term2gene,
      pAdjustMethod = "BH",
      pvalueCutoff  = 0.05,
      qvalueCutoff  = 0.2
    )
    cellmarker_df <- as.data.frame(cellmarker_compare)
  } else {
    cellmarker_df <- data.frame()
  }

  write.csv(
    cellmarker_df,
    file.path(out_dir, paste0("COMPARE_CELL_MARKER_", direction_label, "_TIMEPOINTS_CORRECTED.csv")),
    row.names = FALSE
  )

  results$CellMarker_pathway <- make_pathway_heatmap(
    cellmarker_df, "CellMarker", direction_label, cluster_levels,
    top_n = 30, high_color = high_color
  )
  results$CellMarker_module <- make_module_heatmap(
    cellmarker_df, "CellMarker", direction_label, cluster_levels,
    high_color = high_color
  )

  # GO BP
  if (length(gene_list_entrez_usable) > 0) {
    go_compare <- compareCluster(
      geneCluster   = gene_list_entrez_usable,
      fun           = "enrichGO",
      universe      = background_entrez,
      OrgDb         = org.Mm.eg.db,
      keyType       = "ENTREZID",
      ont           = "BP",
      pAdjustMethod = "BH",
      pvalueCutoff  = 0.05,
      qvalueCutoff  = 0.2,
      readable      = TRUE
    )
    go_df <- as.data.frame(go_compare)
  } else {
    go_df <- data.frame()
  }

  write.csv(
    go_df,
    file.path(out_dir, paste0("COMPARE_GO_BP_", direction_label, "_TIMEPOINTS_CORRECTED.csv")),
    row.names = FALSE
  )

  results$GO_BP_pathway <- make_pathway_heatmap(
    go_df, "GO biological process", direction_label, cluster_levels,
    top_n = 30, high_color = high_color
  )
  results$GO_BP_module <- make_module_heatmap(
    go_df, "GO biological process", direction_label, cluster_levels,
    high_color = high_color
  )

  # KEGG
  if (length(gene_list_entrez_usable) > 0) {
    kegg_compare <- compareCluster(
      geneCluster   = gene_list_entrez_usable,
      fun           = "enrichKEGG",
      organism      = "mmu",
      pvalueCutoff  = 0.05,
      pAdjustMethod = "BH"
    )
    kegg_df <- as.data.frame(kegg_compare)
  } else {
    kegg_df <- data.frame()
  }

  write.csv(
    kegg_df,
    file.path(out_dir, paste0("COMPARE_KEGG_", direction_label, "_TIMEPOINTS_CORRECTED.csv")),
    row.names = FALSE
  )

  results$KEGG_pathway <- make_pathway_heatmap(
    kegg_df, "KEGG", direction_label, cluster_levels,
    top_n = 30, high_color = high_color
  )
  results$KEGG_module <- make_module_heatmap(
    kegg_df, "KEGG", direction_label, cluster_levels,
    high_color = high_color
  )

  # Reactome
  if (length(gene_list_entrez_usable) > 0) {
    reactome_compare <- compareCluster(
      geneCluster   = gene_list_entrez_usable,
      fun           = "enrichPathway",
      organism      = "mouse",
      pvalueCutoff  = 0.05,
      pAdjustMethod = "BH",
      readable      = TRUE
    )
    reactome_df <- as.data.frame(reactome_compare)
  } else {
    reactome_df <- data.frame()
  }

  write.csv(
    reactome_df,
    file.path(out_dir, paste0("COMPARE_REACTOME_", direction_label, "_TIMEPOINTS_CORRECTED.csv")),
    row.names = FALSE
  )

  results$Reactome_pathway <- make_pathway_heatmap(
    reactome_df, "Reactome", direction_label, cluster_levels,
    top_n = 30, high_color = high_color
  )
  results$Reactome_module <- make_module_heatmap(
    reactome_df, "Reactome", direction_label, cluster_levels,
    high_color = high_color
  )

  results
}

up_plots <- run_enrichment_for_direction("UP")
down_plots <- run_enrichment_for_direction("DOWN")

pdf(single_pdf, width = 13, height = 9, onefile = TRUE)

for (nm in names(up_plots)) {
  print(up_plots[[nm]])
}

for (nm in names(down_plots)) {
  print(down_plots[[nm]])
}

dev.off()

cat("\nALL DONE\n")
cat("Single multi-page PDF saved to:\n")
cat(single_pdf, "\n\n")
cat("CSV/check files saved to:\n")
cat(out_dir, "\n")
