# 10_plx_terminal_combined_enrichment.R
#
# Purpose: PLX terminal combined-sex enrichment analysis for UP/DOWN gene lists.
# Inputs are expected under data/processed/ or data/external/ relative to this repository.
# Outputs are written under results/ or script-defined subfolders.

# PLX TERMINAL COMBINED-GENDER GENOME-WIDE UP/DOWN ENRICHMENT
# Gold + Forestgreen overlay volcano output folder
#
# Input:
#   Genome-wide UP_genes_* and DOWN_genes_* files created from:
#
# Goal:
#   Combined male/female version.
#   NO sex split.
#
#   For each volcano comparison, run enrichment separately for:
#     1) UP genes
#     2) DOWN genes
#
#   4 comparisons = 4 analysis pages.
#   Each page contains:
#     - CellMarker pathway heatmap + module heatmap
#     - GO BP pathway heatmap + module heatmap
#     - KEGG pathway heatmap + module heatmap
#     - Reactome pathway heatmap + module heatmap
#
# Output:
#   One 4-page PDF + CSV result tables.
#   Module heatmaps are summarized only from the exact pathway terms displayed in the paired pathway heatmap.
#
#   - This uses genome-wide UP/DOWN gene files.
#   - It does NOT use GOLD_FORESTGREEN_immune_module_genes_found_in_PLX_combined_male_female_volcano_tables.csv.
#   - It does NOT split by Female/Male because the volcano input is combined gender.


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
  library(patchwork)
  library(gridExtra)
  library(grid)
})

select  <- dplyr::select
filter  <- dplyr::filter
mutate  <- dplyr::mutate
arrange <- dplyr::arrange

# 0.5) Real folders
katz_root <- repo_root

project_dir <- file.path(katz_root, "PLX terminal experiment")

base_dir <- file.path(
  project_dir,
  "results",
  "volcano combined male and female",
  "volcano_plots_cleaned_2pages_gold_forestgreen_overlay"
)

if (!dir.exists(base_dir)) {
  stop(
    "Input volcano folder does not exist:\n",
    base_dir,
    "\n\nCheck whether the folder was moved or renamed."
  )
}

# Support files.
# Cell_marker_Mouse.xlsx exists in timepoint_lsd1 and is also duplicated in Trem2 terminal experiment.
cellmarker_candidates <- c(
  file.path(katz_root, "timepoint_lsd1", "Cell_marker_Mouse.xlsx"),
  file.path(katz_root, "Trem2 terminal experiment", "Cell_marker_Mouse.xlsx"),
  file.path(project_dir, "Cell_marker_Mouse.xlsx")
)

cellmarker_file <- cellmarker_candidates[file.exists(cellmarker_candidates)][1]

if (is.na(cellmarker_file) || !file.exists(cellmarker_file)) {
  stop(
    "Could not find Cell_marker_Mouse.xlsx. Checked:\n",
    paste(cellmarker_candidates, collapse = "\n")
  )
}

out_dir <- file.path(base_dir, "genome_wide_UP_DOWN_enrichment_4pages_combined_gender")
if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}

main_pdf <- file.path(out_dir, "PLX_terminal_combined_gender_genomewide_UP_DOWN_enrichment_4pages.pdf")

cat("\nInput volcano folder:\n", base_dir, "\n", sep = "")
cat("\nCellMarker file:\n", cellmarker_file, "\n", sep = "")
cat("\nOutput folder:\n", out_dir, "\n\n", sep = "")

# 1) Mouse symbol normalization
background_entrez <- keys(org.Mm.eg.db, keytype = "ENTREZID")
official_symbols <- keys(org.Mm.eg.db, keytype = "SYMBOL")

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

safe_file_name <- function(x) {
  x %>%
    str_replace_all("[[:space:]]+", "_") %>%
    str_replace_all("[^A-Za-z0-9_\\-]+", "_") %>%
    str_replace_all("_+", "_") %>%
    str_replace_all("_$", "")
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

make_empty_plot <- function(msg) {
  ggplot() +
    annotate("text", x = 0, y = 0, label = msg, size = 4.5) +
    theme_void() +
    xlim(-1, 1) +
    ylim(-1, 1)
}

select_pathway_terms_for_display <- function(enrich_df,
                                             top_n = 12) {

  if (is.null(enrich_df) || nrow(enrich_df) == 0) {
    return(data.frame())
  }

  df <- enrich_df %>%
    mutate(
      score = safe_neglog10(p.adjust),
      Description = as.character(Description),
      Direction = as.character(Cluster)
    )

  # The pathway heatmap displays only the top terms selected here.
  # The module heatmap MUST be summarized from this exact displayed pathway subset.
  # Otherwise the module plot can show modules that are not grounded in the pathway plot.
  top_terms <- df %>%
    group_by(Description) %>%
    summarise(max_score = max(score, na.rm = TRUE), .groups = "drop") %>%
    arrange(desc(max_score)) %>%
    slice_head(n = top_n) %>%
    pull(Description)

  df %>%
    filter(Description %in% top_terms)
}

make_pathway_heatmap <- function(display_df,
                                 analysis_name,
                                 direction_label,
                                 high_color = "firebrick") {

  if (is.null(display_df) || nrow(display_df) == 0) {
    return(make_empty_plot(paste0("No significant ", analysis_name, "\n", direction_label)))
  }

  top_terms <- display_df %>%
    group_by(Description) %>%
    summarise(max_score = max(score, na.rm = TRUE), .groups = "drop") %>%
    arrange(desc(max_score)) %>%
    pull(Description)

  plot_df <- display_df %>%
    mutate(
      Description = factor(Description, levels = rev(top_terms)),
      Direction = factor(Direction, levels = c("UP", "DOWN"))
    )

  ggplot(plot_df, aes(x = Direction, y = Description, fill = score)) +
    geom_tile(color = "white") +
    scale_fill_gradient(
      low = "white",
      high = high_color,
      name = expression(-log[10]("p.adjust"))
    ) +
    labs(
      title = paste0(analysis_name, " pathways"),
      x = NULL,
      y = NULL
    ) +
    theme_bw(base_size = 8.5) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 9),
      axis.text.x = element_text(angle = 0, hjust = 0.5),
      axis.text.y = element_text(size = 6.2),
      panel.grid = element_blank(),
      legend.position = "right"
    )
}

make_module_heatmap <- function(display_df,
                                analysis_name,
                                direction_label,
                                high_color = "firebrick") {

  if (is.null(display_df) || nrow(display_df) == 0) {
    return(make_empty_plot(paste0("No displayed ", analysis_name, " pathway\nfor module summary: ", direction_label)))
  }

  # Modules are NOT independently enriched here.
  # They are a summary/classification of the exact pathways shown in the paired pathway heatmap.
  # Therefore, no module tile can appear unless at least one displayed pathway term exists for that same direction.
  module_df <- display_df %>%
    mutate(
      Module = classify_module(Description),
      Direction = as.character(Direction)
    ) %>%
    group_by(Direction, Module) %>%
    summarise(
      mean_score = mean(score, na.rm = TRUE),
      n_terms = n(),
      .groups = "drop"
    ) %>%
    mutate(
      Direction = factor(Direction, levels = c("UP", "DOWN"))
    )

  if (nrow(module_df) == 0) {
    return(make_empty_plot(paste0("No displayed ", analysis_name, " pathway\nfor module summary: ", direction_label)))
  }

  module_levels <- module_df %>%
    group_by(Module) %>%
    summarise(max_score = max(mean_score, na.rm = TRUE), .groups = "drop") %>%
    arrange(max_score) %>%
    pull(Module)

  module_df <- module_df %>%
    mutate(Module = factor(Module, levels = module_levels))

  ggplot(module_df, aes(x = Direction, y = Module, fill = mean_score)) +
    geom_tile(color = "white") +
    geom_text(aes(label = n_terms), size = 2.6) +
    scale_fill_gradient(
      low = "white",
      high = high_color,
      name = expression("Mean " * -log[10]("p.adjust"))
    ) +
    labs(
      title = paste0(analysis_name, " modules"),
      x = NULL,
      y = NULL
    ) +
    theme_bw(base_size = 8.5) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 9),
      axis.text.x = element_text(angle = 0, hjust = 0.5),
      axis.text.y = element_text(size = 6.8),
      panel.grid = element_blank(),
      legend.position = "right"
    )
}

# 4) Shared CellMarker TERM2GENE
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

# 5) Define the 4 combined-gender volcano panels
panel_info <- tibble::tribble(
  ~comparison_old,                         ~comparison_new,
  "Ctrl_LSD1WT_vs_Ctrl_LSD1KO",            "Ctrl_LSD1KO_vs_Ctrl_LSD1WT",
  "Ctrl_LSD1WT_vs_PLX_LSD1WT",             "PLX_LSD1WT_vs_Ctrl_LSD1WT",
  "Ctrl_LSD1WT_vs_PLX_LSD1KO",             "PLX_LSD1KO_vs_Ctrl_LSD1WT",
  "Ctrl_LSD1KO_vs_PLX_LSD1KO",             "PLX_LSD1KO_vs_Ctrl_LSD1KO"
) %>%
  mutate(
    up_file = file.path(base_dir, paste0("UP_genes_", comparison_old, ".csv")),
    down_file = file.path(base_dir, paste0("DOWN_genes_", comparison_old, ".csv")),
    page_label = comparison_new,
    page_title = comparison_new
  )

missing_inputs <- c(panel_info$up_file, panel_info$down_file)
missing_inputs <- missing_inputs[!file.exists(missing_inputs)]

if (length(missing_inputs) > 0) {
  stop(
    "Missing UP/DOWN input files:\n",
    paste(missing_inputs, collapse = "\n"),
    "\n\nExpected these files inside:\n",
    base_dir
  )
}

read_gene_file <- function(file) {
  df <- read.csv(file, check.names = FALSE)

  gene_col <- dplyr::case_when(
    "Gene_name" %in% names(df) ~ "Gene_name",
    "Gene" %in% names(df) ~ "Gene",
    "gene" %in% names(df) ~ "gene",
    "Symbol" %in% names(df) ~ "Symbol",
    "SYMBOL" %in% names(df) ~ "SYMBOL",
    TRUE ~ NA_character_
  )

  if (is.na(gene_col)) {
    stop(
      "No recognizable gene column in:\n",
      file,
      "\nExpected one of: Gene_name, Gene, gene, Symbol, SYMBOL."
    )
  }

  clean_genes(df[[gene_col]])
}

run_one_panel <- function(up_file, down_file, page_label) {

  cat("\n============================================================\n")
  cat("Running panel: ", page_label, "\n", sep = "")
  cat("UP file: ", up_file, "\n", sep = "")
  cat("DOWN file: ", down_file, "\n", sep = "")
  cat("============================================================\n")

  raw_up <- read_gene_file(up_file)
  raw_down <- read_gene_file(down_file)

  up_symbols <- normalize_mouse_symbols(raw_up)
  down_symbols <- normalize_mouse_symbols(raw_down)

  up_mapped <- map_genes(up_symbols)
  down_mapped <- map_genes(down_symbols)

  gene_list_symbol <- list(
    UP = up_symbols,
    DOWN = down_symbols
  )

  gene_list_entrez <- list(
    UP = up_mapped$ENTREZID,
    DOWN = down_mapped$ENTREZID
  )

  gene_counts <- data.frame(
    Panel = page_label,
    Direction = c("UP", "DOWN"),
    Raw_gene_count = c(length(raw_up), length(raw_down)),
    Official_symbol_count = c(length(up_symbols), length(down_symbols)),
    Entrez_mapped_count = c(nrow(up_mapped), nrow(down_mapped))
  )

  panel_out_dir <- file.path(out_dir, safe_file_name(page_label))
  if (!dir.exists(panel_out_dir)) {
    dir.create(panel_out_dir, recursive = TRUE)
  }

  write.csv(gene_counts, file.path(panel_out_dir, paste0("GENE_COUNTS_", safe_file_name(page_label), ".csv")), row.names = FALSE)
  write.csv(up_mapped, file.path(panel_out_dir, paste0("SYMBOL_to_ENTREZ_", safe_file_name(page_label), "_UP.csv")), row.names = FALSE)
  write.csv(down_mapped, file.path(panel_out_dir, paste0("SYMBOL_to_ENTREZ_", safe_file_name(page_label), "_DOWN.csv")), row.names = FALSE)

  # Keep only usable gene sets.
  gene_list_symbol_usable <- gene_list_symbol[lengths(gene_list_symbol) >= 5]
  gene_list_entrez_usable <- gene_list_entrez[lengths(gene_list_entrez) >= 10]

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
    file.path(panel_out_dir, paste0("COMPARE_CELL_MARKER_", safe_file_name(page_label), "_UP_DOWN.csv")),
    row.names = FALSE
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
    file.path(panel_out_dir, paste0("COMPARE_GO_BP_", safe_file_name(page_label), "_UP_DOWN.csv")),
    row.names = FALSE
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
    file.path(panel_out_dir, paste0("COMPARE_KEGG_", safe_file_name(page_label), "_UP_DOWN.csv")),
    row.names = FALSE
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
    file.path(panel_out_dir, paste0("COMPARE_REACTOME_", safe_file_name(page_label), "_UP_DOWN.csv")),
    row.names = FALSE
  )

  # Build one PDF page
  # First select the exact pathway terms that will be displayed.
  # Then summarize modules ONLY from those displayed pathway terms.
  cellmarker_display_df <- select_pathway_terms_for_display(cellmarker_df, top_n = 12)
  go_display_df         <- select_pathway_terms_for_display(go_df,         top_n = 12)
  kegg_display_df       <- select_pathway_terms_for_display(kegg_df,       top_n = 12)
  reactome_display_df   <- select_pathway_terms_for_display(reactome_df,   top_n = 12)

  write.csv(
    cellmarker_display_df,
    file.path(panel_out_dir, paste0("DISPLAYED_CELL_MARKER_PATHWAYS_USED_FOR_MODULES_", safe_file_name(page_label), "_UP_DOWN.csv")),
    row.names = FALSE
  )
  write.csv(
    go_display_df,
    file.path(panel_out_dir, paste0("DISPLAYED_GO_BP_PATHWAYS_USED_FOR_MODULES_", safe_file_name(page_label), "_UP_DOWN.csv")),
    row.names = FALSE
  )
  write.csv(
    kegg_display_df,
    file.path(panel_out_dir, paste0("DISPLAYED_KEGG_PATHWAYS_USED_FOR_MODULES_", safe_file_name(page_label), "_UP_DOWN.csv")),
    row.names = FALSE
  )
  write.csv(
    reactome_display_df,
    file.path(panel_out_dir, paste0("DISPLAYED_REACTOME_PATHWAYS_USED_FOR_MODULES_", safe_file_name(page_label), "_UP_DOWN.csv")),
    row.names = FALSE
  )

  p_cell_path <- make_pathway_heatmap(cellmarker_display_df, "CellMarker", "UP/DOWN", high_color = "firebrick")
  p_cell_mod  <- make_module_heatmap(cellmarker_display_df, "CellMarker", "UP/DOWN", high_color = "firebrick")

  p_go_path <- make_pathway_heatmap(go_display_df, "GO BP", "UP/DOWN", high_color = "firebrick")
  p_go_mod  <- make_module_heatmap(go_display_df, "GO BP", "UP/DOWN", high_color = "firebrick")

  p_kegg_path <- make_pathway_heatmap(kegg_display_df, "KEGG", "UP/DOWN", high_color = "firebrick")
  p_kegg_mod  <- make_module_heatmap(kegg_display_df, "KEGG", "UP/DOWN", high_color = "firebrick")

  p_react_path <- make_pathway_heatmap(reactome_display_df, "Reactome", "UP/DOWN", high_color = "firebrick")
  p_react_mod  <- make_module_heatmap(reactome_display_df, "Reactome", "UP/DOWN", high_color = "firebrick")

  page_plot <- (
    (p_cell_path | p_cell_mod) /
      (p_go_path | p_go_mod) /
      (p_kegg_path | p_kegg_mod) /
      (p_react_path | p_react_mod)
  ) +
    plot_annotation(
      title = page_label,
      subtitle = paste0(
        "Combined-gender genome-wide enrichment from volcano UP/DOWN genes; ",
        "UP raw/official/mapped = ", length(raw_up), "/", length(up_symbols), "/", nrow(up_mapped),
        "; DOWN raw/official/mapped = ", length(raw_down), "/", length(down_symbols), "/", nrow(down_mapped)
      ),
      theme = theme(
        plot.title = element_text(hjust = 0.5, face = "bold", size = 15),
        plot.subtitle = element_text(hjust = 0.5, size = 10)
      )
    )

  return(list(
    plot = page_plot,
    counts = gene_counts,
    cellmarker = cellmarker_df,
    go = go_df,
    kegg = kegg_df,
    reactome = reactome_df
  ))
}

all_counts <- list()

pdf(main_pdf, width = 17, height = 22, onefile = TRUE)

for (i in seq_len(nrow(panel_info))) {
  res <- run_one_panel(
    up_file = panel_info$up_file[i],
    down_file = panel_info$down_file[i],
    page_label = panel_info$page_label[i]
  )

  print(res$plot)
  all_counts[[panel_info$page_label[i]]] <- res$counts
}

dev.off()

all_counts_df <- bind_rows(all_counts)
write.csv(
  all_counts_df,
  file.path(out_dir, "ALL_PANEL_UP_DOWN_GENE_COUNTS.csv"),
  row.names = FALSE
)

cat("\nALL DONE\n")
cat("Input files were read from:\n")
cat(base_dir, "\n\n")

cat("Main 4-page PDF saved to:\n")
cat(main_pdf, "\n\n")

cat("CSV result folders saved under:\n")
cat(out_dir, "\n\n")

cat("This script used combined-gender genome-wide UP/DOWN volcano files and ignored gold/forestgreen module-only files.\n")
