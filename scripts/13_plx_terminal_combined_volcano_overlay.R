# 13_plx_terminal_combined_volcano_overlay.R
#
# Purpose: PLX terminal combined-sex differential-expression volcano plots with immune-module overlay.
# Inputs are expected under data/processed/ or data/external/ relative to this repository.
# Outputs are written under results/ or script-defined subfolders.

# PLX terminal volcano plots
# Combined male + female version
# PDF only
# 2 pages in one PDF:
#   Page 1: Whole genome gene expression changes
#   Page 2: Whole genome with gold + forestgreen module genes overlaid
#
# Layout on every page:
#   4 volcano plots total, using the 4 files directly under:
#
# Input:
#   Galaxy-style DESeq2 Excel files with unnamed columns
#   Uses Sheet1 from each Excel file for the full volcano background
#
#   This script follows the original volcano logic.
#   Only the module overlay was changed from yellow to gold + forestgreen.
#   The yellow-only and whole-genome-without-yellow pages were removed.


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
results_dir <- file.path(repo_root, "results", "plx_terminal")
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
library(dplyr)
library(tidyr)
library(readxl)
library(ggplot2)
library(ggrepel)
library(patchwork)

base_dir <- file.path(repo_root, "data", "processed", "results")
combined_dir <- file.path(base_dir, "volcano combined male and female")

# Module annotation files
# These are the same module files used in the previous volcano-plot logic.
module_dir <- data_external_dir
allc_file  <- file.path(module_dir, "NIHMS472534-supplement-02.csv")
mod_file   <- file.path(module_dir, "modules.csv")

out_dir <- file.path(combined_dir, "volcano_plots_cleaned_2pages_gold_forestgreen_overlay")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# 1. Input files
# One file has an extra space before .xlsx:
# Ctrl_LSD1WT vs. PLX_LSD1WT .xlsx

plot_info <- tibble::tribble(
  ~comparison_id, ~file_path, ~plot_title,

  "Ctrl_LSD1WT_vs_Ctrl_LSD1KO",
  file.path(combined_dir, "Ctrl_LSD1WT vs. Ctrl_LSD1KO.xlsx"),
  "Ctrl_LSD1WT vs Ctrl_LSD1KO",

  "Ctrl_LSD1WT_vs_PLX_LSD1WT",
  file.path(combined_dir, "Ctrl_LSD1WT vs. PLX_LSD1WT .xlsx"),
  "Ctrl_LSD1WT vs PLX_LSD1WT",

  "Ctrl_LSD1KO_vs_PLX_LSD1KO",
  file.path(combined_dir, "Ctrl_LSD1KO vs. PLX_LSD1KO.xlsx"),
  "Ctrl_LSD1KO vs PLX_LSD1KO",

  "Ctrl_LSD1WT_vs_PLX_LSD1KO",
  file.path(combined_dir, "Ctrl_LSD1WT vs. PLX_LSD1KO.xlsx"),
  "Ctrl_LSD1WT vs PLX_LSD1KO"
)

missing_files <- plot_info$file_path[!file.exists(plot_info$file_path)]
if (length(missing_files) > 0) {
  stop(
    "These input files were not found:\n",
    paste(missing_files, collapse = "\n"),
    "\n\nCheck file names carefully. One file has an extra space before .xlsx: Ctrl_LSD1WT vs. PLX_LSD1WT .xlsx"
  )
}

if (!file.exists(allc_file)) stop("Cannot find NIHMS472534-supplement-02.csv at: ", allc_file)
if (!file.exists(mod_file))  stop("Cannot find modules.csv at: ", mod_file)

clean_volcano_table <- function(file_path, comparison_id) {

  raw <- readxl::read_excel(
    path = file_path,
    sheet = "Sheet1",
    col_names = TRUE,
    .name_repair = "minimal"
  )

  # Expected structure from the original volcano code:
  # blank Ensembl ID column + Base mean + log2(FC) + StdErr + Wald-Stats + P-value + P-adj
  # + chromosome + start + end + strand + blank/extra + gene symbol
  if (ncol(raw) < 13) {
    stop("Expected at least 13 columns in: ", basename(file_path), " but found ", ncol(raw))
  }

  raw <- raw[, 1:13]

  names(raw) <- c(
    "Gene_ID",
    "Base.mean",
    "log2.FC.",
    "StdErr",
    "Wald.Stats",
    "P.value",
    "P.adj",
    "Chr",
    "Start",
    "End",
    "Strand",
    "Extra",
    "Gene_name"
  )

  cleaned <- raw %>%
    mutate(
      Gene_ID    = as.character(Gene_ID),
      Gene_name  = as.character(Gene_name),
      Chr        = as.character(Chr),
      Strand     = as.character(Strand),
      Base.mean  = as.numeric(Base.mean),
      log2.FC.   = as.numeric(log2.FC.),
      StdErr     = as.numeric(StdErr),
      Wald.Stats = as.numeric(Wald.Stats),
      P.value    = as.numeric(P.value),
      P.adj      = as.numeric(P.adj),
      Start      = as.numeric(Start),
      End        = as.numeric(End),
      Comparison = comparison_id
    ) %>%
    filter(
      !is.na(Gene_name),
      Gene_name != "",
      !is.na(log2.FC.),
      !is.na(P.adj)
    ) %>%
    arrange(P.adj) %>%
    distinct(Gene_name, .keep_all = TRUE) %>%
    arrange(Gene_name) %>%
    mutate(
      P.adj.plot = ifelse(P.adj <= 0, .Machine$double.xmin, P.adj),
      neg_log10_Padj = -log10(P.adj.plot),
      Gene_UP = toupper(trimws(as.character(Gene_name))),
      Significance = case_when(
        P.adj < 0.05 & log2.FC. > 1  ~ "Up regulated",
        P.adj < 0.05 & log2.FC. < -1 ~ "Down regulated",
        TRUE ~ "Not Significant"
      )
    )

  return(cleaned)
}

volcano_list <- vector("list", nrow(plot_info))

for (i in seq_len(nrow(plot_info))) {
  df <- clean_volcano_table(
    file_path     = plot_info$file_path[i],
    comparison_id = plot_info$comparison_id[i]
  )

  volcano_list[[i]] <- df

  clean_name <- paste0("CLEANED_", plot_info$comparison_id[i], ".csv")
  write.csv(df, file.path(out_dir, clean_name), row.names = FALSE)
}

names(volcano_list) <- plot_info$comparison_id

all_volcano_df <- bind_rows(volcano_list)

summary_df <- all_volcano_df %>%
  group_by(Comparison, Significance) %>%
  summarise(n = n(), .groups = "drop") %>%
  tidyr::pivot_wider(names_from = Significance, values_from = n, values_fill = 0)

write.csv(summary_df, file.path(out_dir, "PLX_terminal_combined_male_female_volcano_significance_counts.csv"), row.names = FALSE)

# 4. Build gold + forestgreen immune-module gene sets

allc <- read.csv(allc_file, header = TRUE, check.names = FALSE)
mod  <- read.csv(mod_file, header = TRUE, check.names = FALSE)

# Clean column names robustly.
# Some annotation CSVs contain blank/NA header names, which breaks dplyr::mutate().
names(allc) <- sub("^\\ufeff", "", names(allc))
names(mod)  <- sub("^\\ufeff", "", names(mod))

bad_allc_names <- is.na(names(allc)) | names(allc) == ""
if (any(bad_allc_names)) names(allc)[bad_allc_names] <- paste0("V", seq_len(sum(bad_allc_names)))
names(allc) <- make.unique(names(allc))

bad_mod_names <- is.na(names(mod)) | names(mod) == ""
if (any(bad_mod_names)) names(mod)[bad_mod_names] <- paste0("V", seq_len(sum(bad_mod_names)))
names(mod) <- make.unique(names(mod))

if (!all(c("Gene_Symbol", "Module") %in% names(allc))) {
  stop("NIHMS472534-supplement-02.csv must contain Gene_Symbol and Module columns.")
}
if (!all(c("CategoryTerm", "Module") %in% names(mod))) {
  stop("modules.csv must contain CategoryTerm and Module columns.")
}

target_term <- "Immune functions"
target_modules <- c("gold", "forestgreen")

immune_modules <- mod %>%
  filter(trimws(CategoryTerm) == target_term) %>%
  pull(Module) %>%
  unique()

if (length(immune_modules) == 0) {
  stop("No modules found for CategoryTerm == 'Immune functions'.")
}

gene2module <- allc %>%
  mutate(
    Gene_Symbol = toupper(trimws(as.character(Gene_Symbol))),
    Module_clean = tolower(trimws(as.character(Module)))
  ) %>%
  select(Gene_Symbol, Module, Module_clean) %>%
  filter(!is.na(Gene_Symbol), Gene_Symbol != "", !is.na(Module)) %>%
  distinct()

all_genes_universe <- unique(all_volcano_df$Gene_UP)

gold_forestgreen_gene2module <- gene2module %>%
  filter(Module %in% immune_modules) %>%
  filter(Module_clean %in% target_modules) %>%
  filter(Gene_Symbol %in% all_genes_universe) %>%
  mutate(
    Overlay_Module = case_when(
      Module_clean == "gold" ~ "Gold module",
      Module_clean == "forestgreen" ~ "Forestgreen module",
      TRUE ~ Module
    )
  ) %>%
  arrange(Overlay_Module, Gene_Symbol) %>%
  distinct(Gene_Symbol, Overlay_Module, .keep_all = TRUE)

gold_set <- gold_forestgreen_gene2module %>%
  filter(Overlay_Module == "Gold module") %>%
  pull(Gene_Symbol) %>%
  unique()

forestgreen_set <- gold_forestgreen_gene2module %>%
  filter(Overlay_Module == "Forestgreen module") %>%
  pull(Gene_Symbol) %>%
  unique()

gold_forestgreen_set <- unique(c(gold_set, forestgreen_set))

if (length(gold_forestgreen_set) == 0) {
  stop("Gold + forestgreen gene set is empty after intersecting with volcano tables.")
}

write.csv(
  gold_forestgreen_gene2module %>% select(Gene_Symbol, Overlay_Module),
  file.path(out_dir, "GOLD_FORESTGREEN_immune_module_genes_found_in_PLX_combined_male_female_volcano_tables.csv"),
  row.names = FALSE
)

cat("Gold immune-module genes found in volcano tables:", length(gold_set), "\n")
cat("Forestgreen immune-module genes found in volcano tables:", length(forestgreen_set), "\n")
cat("Gold + forestgreen immune-module genes found in volcano tables:", length(gold_forestgreen_set), "\n")
cat("Whole-genome genes:", length(all_genes_universe), "\n")

# 5. Unified x/y axis limits across all pages

X_LIM <- c(-10, 10)

# Keep the old lower limit style from the previous volcano code.
# Top is dynamic but rounded to a clean 5-unit interval so no high-significance points are clipped.
Y_AXIS_TOP <- max(30, ceiling(max(all_volcano_df$neg_log10_Padj, na.rm = TRUE) / 5) * 5)
Y_LIM <- c(-2, Y_AXIS_TOP)

cat("Unified x-axis:", paste(X_LIM, collapse = " to "), "\n")
cat("Unified y-axis:", paste(Y_LIM, collapse = " to "), "\n")

legend_no_title <- list(
  guides(
    color = guide_legend(
      title = NULL,
      override.aes = list(label = "", shape = 16, size = 3, alpha = 1)
    ),
    fill  = guide_legend(
      title = NULL,
      override.aes = list(shape = 21, color = "black", size = 3, alpha = 1)
    ),
    shape = "none",
    size  = "none",
    alpha = "none"
  ),
  theme(legend.title = element_blank())
)

add_sig_class <- function(df) {
  df %>%
    mutate(
      Significance = case_when(
        P.adj < 0.05 & log2.FC. > 1  ~ "Up regulated",
        P.adj < 0.05 & log2.FC. < -1 ~ "Down regulated",
        TRUE ~ "Not Significant"
      )
    )
}

make_volcano_plot <- function(df, plot_title, subtitle_text = NULL, top_n = 15) {

  if (nrow(df) == 0) {
    # Empty gene-set plots should not crash the entire PDF.
    empty_df <- data.frame(x = 0, y = 0)
    return(
      ggplot(empty_df, aes(x = x, y = y)) +
        geom_blank() +
        coord_cartesian(xlim = X_LIM, ylim = Y_LIM) +
        labs(
          title = plot_title,
          subtitle = paste0(subtitle_text, "\nNo genes in this set"),
          x = "log2(Fold Change)",
          y = "-log10(P.adj)"
        ) +
        theme_minimal(base_size = 12) +
        theme(
          plot.title = element_text(hjust = 0.5, face = "bold", size = 10),
          plot.subtitle = element_text(hjust = 0.5, size = 8),
          axis.title = element_text(size = 10),
          axis.text = element_text(size = 8),
          legend.position = "none"
        )
    )
  }

  df <- add_sig_class(df)

  top_genes <- df %>%
    filter(P.adj < 0.05, abs(log2.FC.) > 1) %>%
    arrange(P.adj) %>%
    slice_head(n = top_n)

  ggplot(df, aes(x = log2.FC., y = neg_log10_Padj, color = Significance)) +
    geom_point(alpha = 0.65, size = 1.25, show.legend = TRUE) +
    geom_text_repel(
      data = top_genes,
      aes(label = Gene_name),
      size = 2.7,
      max.overlaps = Inf,
      box.padding = 0.35,
      point.padding = 0.25,
      min.segment.length = 0,
      show.legend = FALSE
    ) +
    scale_color_manual(
      name = NULL,
      values = c(
        "Up regulated" = "firebrick",
        "Down regulated" = "steelblue",
        "Not Significant" = "grey70"
      ),
      breaks = c("Down regulated", "Not Significant", "Up regulated")
    ) +
    geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
    coord_cartesian(xlim = X_LIM, ylim = Y_LIM) +
    labs(
      title = plot_title,
      subtitle = subtitle_text,
      x = "log2(Fold Change)",
      y = "-log10(P.adj)"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 9.5),
      plot.subtitle = element_text(hjust = 0.5, size = 8),
      axis.title = element_text(size = 9.5),
      axis.text = element_text(size = 8),
      legend.title = element_blank()
    ) +
    legend_no_title
}

make_volcano_with_gold_forestgreen_overlay <- function(df, plot_title, subtitle_text = NULL, top_n = 15,
                                                       gold_genes = gold_set,
                                                       forestgreen_genes = forestgreen_set) {

  df2 <- add_sig_class(df) %>%
    mutate(
      Overlay_Module = case_when(
        Gene_UP %in% gold_genes ~ "Gold module",
        Gene_UP %in% forestgreen_genes ~ "Forestgreen module",
        TRUE ~ NA_character_
      )
    )

  top_genes <- df2 %>%
    filter(P.adj < 0.05, abs(log2.FC.) > 1) %>%
    arrange(P.adj) %>%
    slice_head(n = top_n)

  ggplot(df2, aes(x = log2.FC., y = neg_log10_Padj)) +
    geom_point(aes(color = Significance), alpha = 0.65, size = 1.25, show.legend = TRUE) +
    geom_point(
      data = df2 %>% filter(!is.na(Overlay_Module)),
      aes(x = log2.FC., y = neg_log10_Padj, fill = Overlay_Module),
      shape = 21,
      color = "black",
      stroke = 0.35,
      size = 1.15,
      alpha = 0.95,
      show.legend = TRUE
    ) +
    geom_text_repel(
      data = top_genes,
      aes(label = Gene_name),
      size = 2.7,
      max.overlaps = Inf,
      box.padding = 0.35,
      point.padding = 0.25,
      min.segment.length = 0,
      show.legend = FALSE
    ) +
    scale_color_manual(
      name = NULL,
      values = c(
        "Up regulated" = "firebrick",
        "Down regulated" = "steelblue",
        "Not Significant" = "grey70"
      ),
      breaks = c("Down regulated", "Not Significant", "Up regulated")
    ) +
    scale_fill_manual(
      name = NULL,
      values = c(
        "Gold module" = "gold",
        "Forestgreen module" = "forestgreen"
      ),
      breaks = c("Gold module", "Forestgreen module")
    ) +
    geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
    coord_cartesian(xlim = X_LIM, ylim = Y_LIM) +
    labs(
      title = plot_title,
      subtitle = subtitle_text,
      x = "log2(Fold Change)",
      y = "-log10(P.adj)"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 9.5),
      plot.subtitle = element_text(hjust = 0.5, size = 8),
      axis.title = element_text(size = 9.5),
      axis.text = element_text(size = 8),
      legend.title = element_blank()
    ) +
    legend_no_title
}

build_4plot_page <- function(plot_type = c("whole", "gold_forestgreen_overlay")) {

  plot_type <- match.arg(plot_type)
  p_list <- vector("list", nrow(plot_info))

  for (i in seq_len(nrow(plot_info))) {
    key <- plot_info$comparison_id[i]
    df <- volcano_list[[key]]

    if (plot_type == "whole") {
      p_list[[i]] <- make_volcano_plot(
        df = df,
        plot_title = plot_info$plot_title[i],
        subtitle_text = "Whole genome",
        top_n = 15
      )
    }

    if (plot_type == "gold_forestgreen_overlay") {
      p_list[[i]] <- make_volcano_with_gold_forestgreen_overlay(
        df = df,
        plot_title = plot_info$plot_title[i],
        subtitle_text = "Gold + forestgreen module genes overlaid",
        top_n = 15,
        gold_genes = gold_set,
        forestgreen_genes = forestgreen_set
      )
    }
  }

  page_title <- switch(
    plot_type,
    whole = "PLX terminal combined male and female volcano plots: whole genome",
    gold_forestgreen_overlay = "PLX terminal combined male and female volcano plots: gold + forestgreen module overlay"
  )

  wrap_plots(p_list, ncol = 2, byrow = TRUE) +
    plot_layout(guides = "collect") +
    plot_annotation(
      title = page_title,
      subtitle = "Each page uses the same x/y axis limits; input files are combined male and female",
      theme = theme(
        plot.title = element_text(hjust = 0.5, face = "bold", size = 18),
        plot.subtitle = element_text(hjust = 0.5, size = 11)
      )
    ) &
    theme(
      legend.position = "right",
      legend.title = element_blank()
    ) &
    guides(
      color = guide_legend(
        title = NULL,
        override.aes = list(shape = 16, size = 3, alpha = 1, label = "")
      ),
      fill = guide_legend(
        title = NULL,
        override.aes = list(shape = 21, color = "black", size = 3, alpha = 1)
      )
    )
}

# 7. Build 2 pages and save ONE PDF only

combined_whole <- build_4plot_page("whole")
combined_gold_forestgreen_overlay <- build_4plot_page("gold_forestgreen_overlay")

graphics.off()

pdf_file <- file.path(out_dir, "PLX_terminal_COMBINED_MALE_FEMALE_VOLCANO_2pages_whole_goldForestgreen_overlay.pdf")

pdf(pdf_file, width = 16, height = 16, onefile = TRUE)
print(combined_whole)
print(combined_gold_forestgreen_overlay)
dev.off()

# 8. Export UP/DOWN gene lists for downstream enrichment

get_sig_table <- function(df, direction) {
  if (direction == "UP") {
    df %>% filter(P.adj < 0.05, log2.FC. > 1) %>% arrange(P.adj)
  } else if (direction == "DOWN") {
    df %>% filter(P.adj < 0.05, log2.FC. < -1) %>% arrange(P.adj)
  } else {
    stop("direction must be UP or DOWN")
  }
}

for (i in seq_len(nrow(plot_info))) {
  key <- plot_info$comparison_id[i]
  df <- volcano_list[[key]]

  up_file <- file.path(out_dir, paste0("UP_genes_", plot_info$comparison_id[i], ".csv"))
  down_file <- file.path(out_dir, paste0("DOWN_genes_", plot_info$comparison_id[i], ".csv"))

  write.csv(get_sig_table(df, "UP"), up_file, row.names = FALSE)
  write.csv(get_sig_table(df, "DOWN"), down_file, row.names = FALSE)
}

cat("Done. Files saved to:\n", out_dir, "\n")
cat("Main 2-page PDF only:\n", pdf_file, "\n")
cat("No PNG was created.\n")
cat("Cleaned CSVs, UP/DOWN gene lists, and gold + forestgreen gene list were also exported.\n")
