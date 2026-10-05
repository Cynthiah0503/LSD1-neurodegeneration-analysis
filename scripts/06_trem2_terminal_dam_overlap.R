# 06_trem2_terminal_dam_overlap.R
#
# Purpose: Overlap analysis between Trem2/LSD1 differential-expression gene lists and DAM-associated genes.
# Inputs are expected under data/processed/ or data/external/ relative to this repository.
# Outputs are written under results/ or script-defined subfolders.

# Trem2 Venn analysis:
# Trem2 UP/DOWN significant gene lists vs DAM gene lists
#
# This does NOT compare UP vs DOWN.
# It compares each Trem2 UP list vs each DAM list.
# It compares each Trem2 DOWN list vs each DAM list.


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
library(readr)
library(readxl)
library(dplyr)
library(stringr)
library(purrr)
library(VennDiagram)
library(grid)
library(tibble)

katz_dir <- repo_root

trem2_dir <- file.path(katz_dir, "Trem2 terminal experiment")

dam_dir <- file.path(
  katz_dir,
  "timepoint_lsd1/combine terminal/正确volcano enrichment：combine terminal batch aware"
)

dam_file <- file.path(dam_dir, "DAM.xlsx")

out_dir <- file.path(trem2_dir, "DAM_overlap_Trem2_UP_DOWN")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

venn_pdf <- file.path(out_dir, "VENN_Trem2_UP_DOWN_vs_DAM.pdf")
summary_csv <- file.path(out_dir, "SUMMARY_Trem2_UP_DOWN_vs_DAM.csv")
overlap_csv <- file.path(out_dir, "OVERLAP_GENES_Trem2_UP_DOWN_vs_DAM.csv")
trem2_check_csv <- file.path(out_dir, "TREM2_presence_check_Trem2_vs_DAM.csv")

# 2. Trem2 UP and DOWN files

up_files <- c(
  "UP_Trem2KO_LSD1WT_vs_Trem2WT_LSD1WT" =
    "UP_Trem2KO_LSD1_WT_vs_Trem2_WT_LSD1_WT.csv",

  "UP_Trem2KO_LSD1Del_vs_Trem2WT_LSD1WT" =
    "UP_Trem2KO_LSD1_Del_vs_Trem2_WT_LSD1_WT.csv",

  "UP_Trem2KO_LSD1Del_vs_Trem2WT_LSD1Del" =
    "UP_Trem2KO_LSD1_Del_vs_Trem2_WT_LSD1_Del.csv",

  "UP_Trem2WT_LSD1Del_vs_Trem2WT_LSD1WT" =
    "UP_Trem2_WT_LSD1_Del_vs_Trem2_WT_LSD1_WT.csv"
)

down_files <- c(
  "DOWN_Trem2KO_LSD1WT_vs_Trem2WT_LSD1WT" =
    "DOWN_genes_Trem2KO_LSD1WT_vs_Trem2WT_LSD1WT.csv",

  "DOWN_Trem2KO_LSD1KO_vs_Trem2WT_LSD1WT" =
    "DOWN_genes_Trem2KO_LSD1KO_vs_Trem2WT_LSD1WT.csv",

  "DOWN_Trem2KO_LSD1KO_vs_Trem2WT_LSD1KO" =
    "DOWN_genes_Trem2KO_LSD1KO_vs_Trem2WT_LSD1KO.csv",

  "DOWN_Trem2WT_LSD1KO_vs_Trem2WT_LSD1WT" =
    "DOWN_genes_Trem2WT_LSD1KO_vs_Trem2WT_LSD1WT.csv"
)

up_files <- file.path(trem2_dir, up_files)
down_files <- file.path(trem2_dir, down_files)

names(up_files) <- c(
  "UP_Trem2KO_LSD1WT_vs_Trem2WT_LSD1WT",
  "UP_Trem2KO_LSD1Del_vs_Trem2WT_LSD1WT",
  "UP_Trem2KO_LSD1Del_vs_Trem2WT_LSD1Del",
  "UP_Trem2WT_LSD1Del_vs_Trem2WT_LSD1WT"
)

names(down_files) <- c(
  "DOWN_Trem2KO_LSD1WT_vs_Trem2WT_LSD1WT",
  "DOWN_Trem2KO_LSD1KO_vs_Trem2WT_LSD1WT",
  "DOWN_Trem2KO_LSD1KO_vs_Trem2WT_LSD1KO",
  "DOWN_Trem2WT_LSD1KO_vs_Trem2WT_LSD1WT"
)

# 3. Check files exist

all_input_files <- c(up_files, down_files, DAM = dam_file)

missing_files <- all_input_files[!file.exists(all_input_files)]

if (length(missing_files) > 0) {
  stop(
    "These input files do not exist:\n",
    paste(missing_files, collapse = "\n")
  )
}

clean_gene_vector <- function(x) {
  x %>%
    as.character() %>%
    str_trim() %>%
    str_replace_all("[\r\n\t]", "") %>%
    na_if("") %>%
    discard(is.na) %>%
    toupper() %>%
    unique()
}

guess_gene_column <- function(df) {

  possible_gene_cols <- c(
    "Gene_name", "gene_name",
    "Gene", "gene",
    "GeneSymbol", "gene_symbol",
    "Symbol", "SYMBOL",
    "Mouse_gene", "mouse_gene",
    "geneName",
    "gene_name_final"
  )

  hit <- intersect(possible_gene_cols, colnames(df))

  if (length(hit) > 0) {
    return(hit[1])
  }

  non_numeric_cols <- colnames(df)[!sapply(df, is.numeric)]

  if (length(non_numeric_cols) == 0) {
    stop("No gene-like column found. Please check column names.")
  }

  return(non_numeric_cols[1])
}

read_gene_list_csv <- function(file) {
  df <- read_csv(file, show_col_types = FALSE)
  gene_col <- guess_gene_column(df)
  genes <- clean_gene_vector(df[[gene_col]])

  cat("\nReading file:\n", file, "\n")
  cat("Using gene column:", gene_col, "\n")
  cat("Number of genes:", length(genes), "\n")
  cat("Contains TREM2:", "TREM2" %in% genes, "\n")

  return(genes)
}

up_gene_lists <- map(up_files, read_gene_list_csv)
down_gene_lists <- map(down_files, read_gene_list_csv)

names(up_gene_lists) <- names(up_files)
names(down_gene_lists) <- names(down_files)

trem2_gene_lists <- c(up_gene_lists, down_gene_lists)

trem2_size_check <- tibble(
  list_name = names(trem2_gene_lists),
  direction = ifelse(str_starts(names(trem2_gene_lists), "UP_"), "UP", "DOWN"),
  n_genes = map_int(trem2_gene_lists, length),
  contains_TREM2 = map_lgl(trem2_gene_lists, ~ "TREM2" %in% .x)
)

cat("\nTrem2 gene list check:\n")
print(trem2_size_check)

# Handles:
# A) DAM.xlsx has 3 sheets, each sheet = one gene list
# B) DAM.xlsx has 1 sheet with 3 gene-list columns

dam_sheets <- excel_sheets(dam_file)

if (length(dam_sheets) >= 3) {

  dam_gene_lists <- map(dam_sheets, function(sh) {
    df <- read_excel(dam_file, sheet = sh)
    gene_col <- guess_gene_column(df)
    genes <- clean_gene_vector(df[[gene_col]])

    cat("\nReading DAM sheet:", sh, "\n")
    cat("Using gene column:", gene_col, "\n")
    cat("Number of genes:", length(genes), "\n")
    cat("Contains TREM2:", "TREM2" %in% genes, "\n")

    return(genes)
  })

  names(dam_gene_lists) <- make.names(dam_sheets)

} else {

  dam_df <- read_excel(dam_file, sheet = dam_sheets[1])

  gene_list_cols <- colnames(dam_df)[!sapply(dam_df, is.numeric)]

  if (length(gene_list_cols) < 3) {
    stop("DAM.xlsx does not look like 3 sheets or 3 gene-list columns. Please check the file.")
  }

  dam_gene_lists <- map(gene_list_cols, function(col) {
    genes <- clean_gene_vector(dam_df[[col]])

    cat("\nReading DAM column:", col, "\n")
    cat("Number of genes:", length(genes), "\n")
    cat("Contains TREM2:", "TREM2" %in% genes, "\n")

    return(genes)
  })

  names(dam_gene_lists) <- make.names(gene_list_cols)
}

dam_size_check <- tibble(
  DAM_list = names(dam_gene_lists),
  n_genes = map_int(dam_gene_lists, length),
  contains_TREM2 = map_lgl(dam_gene_lists, ~ "TREM2" %in% .x)
)

cat("\nDAM gene list check:\n")
print(dam_size_check)

# 7. TREM2 explicit check

all_lists_for_trem2 <- c(trem2_gene_lists, dam_gene_lists)

trem2_check <- tibble(
  list_name = names(all_lists_for_trem2),
  n_genes = map_int(all_lists_for_trem2, length),
  contains_TREM2 = map_lgl(all_lists_for_trem2, ~ "TREM2" %in% .x)
)

write_csv(trem2_check, trem2_check_csv)

cat("\nTREM2 presence check:\n")
print(trem2_check)

# 8. Pairwise overlap:
# Trem2 UP/DOWN lists vs DAM lists
# NO UP vs DOWN comparison

summary_results <- list()
overlap_results <- list()

counter <- 1

for (trem2_name in names(trem2_gene_lists)) {

  trem2_genes <- trem2_gene_lists[[trem2_name]]
  direction <- ifelse(str_starts(trem2_name, "UP_"), "UP", "DOWN")

  for (dam_name in names(dam_gene_lists)) {

    dam_genes <- dam_gene_lists[[dam_name]]
    overlap_genes <- intersect(trem2_genes, dam_genes)

    summary_results[[counter]] <- tibble(
      Trem2_list = trem2_name,
      direction = direction,
      DAM_list = dam_name,
      Trem2_n = length(trem2_genes),
      DAM_n = length(dam_genes),
      overlap_n = length(overlap_genes),
      Trem2_percent_overlap = ifelse(
        length(trem2_genes) > 0,
        length(overlap_genes) / length(trem2_genes) * 100,
        NA_real_
      ),
      DAM_percent_overlap = ifelse(
        length(dam_genes) > 0,
        length(overlap_genes) / length(dam_genes) * 100,
        NA_real_
      ),
      has_overlap = length(overlap_genes) > 0,
      contains_TREM2_overlap = "TREM2" %in% overlap_genes,
      overlap_genes = paste(overlap_genes, collapse = "; ")
    )

    if (length(overlap_genes) > 0) {
      overlap_results[[counter]] <- tibble(
        Trem2_list = trem2_name,
        direction = direction,
        DAM_list = dam_name,
        overlap_gene = overlap_genes
      )
    } else {
      overlap_results[[counter]] <- tibble(
        Trem2_list = trem2_name,
        direction = direction,
        DAM_list = dam_name,
        overlap_gene = character(0)
      )
    }

    counter <- counter + 1
  }
}

summary_df <- bind_rows(summary_results)
overlap_df <- bind_rows(overlap_results)

write_csv(summary_df, summary_csv)
write_csv(overlap_df, overlap_csv)

cat("\nOverlap summary:\n")
print(summary_df)

# 9. Venn PDF

pdf(venn_pdf, width = 8, height = 6.8)

for (trem2_name in names(trem2_gene_lists)) {

  trem2_genes <- trem2_gene_lists[[trem2_name]]

  for (dam_name in names(dam_gene_lists)) {

    dam_genes <- dam_gene_lists[[dam_name]]
    overlap_genes <- intersect(trem2_genes, dam_genes)
    overlap_n <- length(overlap_genes)

    grid.newpage()

    venn_plot <- draw.pairwise.venn(
      area1 = length(trem2_genes),
      area2 = length(dam_genes),
      cross.area = overlap_n,
      category = c(trem2_name, dam_name),
      fill = c("#9ecae1", "#fdae6b"),
      alpha = c(0.55, 0.55),
      lty = "blank",
      cex = 1.45,
      cat.cex = 0.72,
      cat.pos = c(-20, 20),
      cat.dist = c(0.05, 0.05),
      scaled = FALSE
    )

    grid.draw(venn_plot)

    grid.text(
      paste0(
        trem2_name, " vs ", dam_name,
        "\nOverlap genes = ", overlap_n
      ),
      x = 0.5,
      y = 0.96,
      gp = gpar(fontsize = 11.5, fontface = "bold")
    )

    if (overlap_n > 0) {
      grid.text(
        paste0("Overlap: ", paste(overlap_genes, collapse = ", ")),
        x = 0.5,
        y = 0.05,
        gp = gpar(fontsize = 7.5)
      )
    }
  }
}

dev.off()

# 10. Final output

cat("\nDone.\n")
cat("Venn PDF saved to:\n", venn_pdf, "\n\n")
cat("Summary saved to:\n", summary_csv, "\n\n")
cat("Overlap genes saved to:\n", overlap_csv, "\n\n")
cat("TREM2 check saved to:\n", trem2_check_csv, "\n\n")
