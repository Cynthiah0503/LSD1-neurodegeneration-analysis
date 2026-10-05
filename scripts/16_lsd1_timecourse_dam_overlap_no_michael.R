# 16_lsd1_timecourse_dam_overlap_no_michael.R
#
# Purpose: Corrected overlap analysis between LSD1 time-course gene lists and DAM-associated genes.
# Inputs are expected under data/processed/ or data/external/ relative to this repository.
# Outputs are written under results/ or script-defined subfolders.

# CORRECTED ALL-timepoint DAM overlap Venn
#
# FIXED:
#   The corrected volcano output uses these file names:
#     UP_genes_3weeks_CORRECTED_COMBINE_terminal.csv
#     DOWN_genes_3weeks_CORRECTED_COMBINE_terminal.csv
#     UP_genes_4weeks_CORRECTED_COMBINE_terminal.csv
#     DOWN_genes_4weeks_CORRECTED_COMBINE_terminal.csv
#     UP_genes_Early_onset_CORRECTED_COMBINE_terminal.csv
#     DOWN_genes_Early_onset_CORRECTED_COMBINE_terminal.csv
#     UP_genes_Terminal_CORRECTED_LSD1terminal_plus_Trem2WTterminal_NO_MICHAEL.csv
#     DOWN_genes_Terminal_CORRECTED_LSD1terminal_plus_Trem2WTterminal_NO_MICHAEL.csv
#
# Goal:
#   Compare ALL corrected volcano UP/DOWN gene lists against DAM gene lists:
#     3 weeks UP/DOWN vs DAM
#     4 weeks UP/DOWN vs DAM
#     Early onset UP/DOWN vs DAM
#     Terminal corrected UP/DOWN vs DAM
#
# Output:
#   One multi-page Venn PDF + summary CSV + overlap gene CSV + TREM2 check


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
  library(readr)
  library(readxl)
  library(dplyr)
  library(stringr)
  library(purrr)
  library(VennDiagram)
  library(grid)
  library(tibble)
  library(ggplot2)
})

katz_root <- repo_root
project_dir <- file.path(katz_root, "timepoint_lsd1")
corrected_root <- file.path(project_dir, "正确combine terminal")

base_dir <- file.path(
  corrected_root,
  "CORRECTED_volcano_LSD1terminal_plus_Trem2WTterminal_no_Michael_GALAXY_style_FC"
)

dam_file <- file.path(corrected_root, "DAM.xlsx")

out_dir <- file.path(
  base_dir,
  "DAM_overlap_ALL_timepoints_CORRECTED_LSD1terminal_plus_Trem2WTterminal_NO_MICHAEL"
)
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

venn_pdf <- file.path(out_dir, "VENN_ALL_timepoints_UP_DOWN_vs_DAM_CORRECTED_NO_MICHAEL.pdf")
summary_csv <- file.path(out_dir, "SUMMARY_ALL_timepoints_UP_DOWN_vs_DAM_CORRECTED_NO_MICHAEL.csv")
overlap_csv <- file.path(out_dir, "OVERLAP_GENES_ALL_timepoints_UP_DOWN_vs_DAM_CORRECTED_NO_MICHAEL.csv")
trem2_check_csv <- file.path(out_dir, "TREM2_presence_check_ALL_timepoints_vs_DAM_CORRECTED_NO_MICHAEL.csv")
input_check_csv <- file.path(out_dir, "CHECK_input_files_CORRECTED_DAM_overlap_ALL_timepoints.csv")
heatmap_pdf <- file.path(out_dir, "HEATMAP_overlap_counts_ALL_timepoints_UP_DOWN_vs_DAM_CORRECTED_NO_MICHAEL.pdf")

# 2. Corrected volcano UP/DOWN files

timepoint_files <- tribble(
  ~Timepoint,     ~Direction, ~List_name,          ~File,
  "3 weeks",     "UP",      "UP_3weeks",        file.path(base_dir, "UP_genes_3weeks_CORRECTED_COMBINE_terminal.csv"),
  "3 weeks",     "DOWN",    "DOWN_3weeks",      file.path(base_dir, "DOWN_genes_3weeks_CORRECTED_COMBINE_terminal.csv"),
  "4 weeks",     "UP",      "UP_4weeks",        file.path(base_dir, "UP_genes_4weeks_CORRECTED_COMBINE_terminal.csv"),
  "4 weeks",     "DOWN",    "DOWN_4weeks",      file.path(base_dir, "DOWN_genes_4weeks_CORRECTED_COMBINE_terminal.csv"),
  "Early onset", "UP",      "UP_Early_onset",   file.path(base_dir, "UP_genes_Early_onset_CORRECTED_COMBINE_terminal.csv"),
  "Early onset", "DOWN",    "DOWN_Early_onset", file.path(base_dir, "DOWN_genes_Early_onset_CORRECTED_COMBINE_terminal.csv"),
  "Terminal",    "UP",      "UP_Terminal",      file.path(base_dir, "UP_genes_Terminal_CORRECTED_LSD1terminal_plus_Trem2WTterminal_NO_MICHAEL.csv"),
  "Terminal",    "DOWN",    "DOWN_Terminal",    file.path(base_dir, "DOWN_genes_Terminal_CORRECTED_LSD1terminal_plus_Trem2WTterminal_NO_MICHAEL.csv")
)

input_check <- bind_rows(
  tibble(
    input_name = c("base_dir", "DAM.xlsx"),
    path = c(base_dir, dam_file),
    exists = c(dir.exists(base_dir), file.exists(dam_file))
  ),
  timepoint_files %>%
    transmute(
      input_name = paste(List_name, "corrected"),
      path = File,
      exists = file.exists(File)
    )
)

write_csv(input_check, input_check_csv)

cat("\nInput file check:\n")
print(input_check)

if (!dir.exists(base_dir)) {
  stop("Corrected volcano output folder does not exist:\n", base_dir)
}
if (!file.exists(dam_file)) {
  stop("DAM.xlsx not found:\n", dam_file)
}

missing_timepoint_files <- timepoint_files$File[!file.exists(timepoint_files$File)]
if (length(missing_timepoint_files) > 0) {
  cat("\nFiles available in base_dir:\n")
  print(list.files(base_dir, pattern = "genes_.*\\.csv$", full.names = FALSE))
  stop(
    "Missing corrected UP/DOWN gene files:\n",
    paste(missing_timepoint_files, collapse = "\n"),
    "\n\nI printed available gene-list CSV files above. Check file names in:\n",
    base_dir
  )
}

# Safety check:
# Allow NO_MICHAEL in corrected file names.
# Stop only for old explicit Michael or old batch-aware folder strings.
forbidden_patterns <- c(
  "Terminal_plus_Michael",
  "featurecountsLSD1del",
  "featurecountsLSD1del vs",
  "combine terminal batch aware",
  "正确volcano enrichment"
)

bad_paths <- c(base_dir, timepoint_files$File)
bad_hit <- unlist(lapply(forbidden_patterns, function(pat) {
  grepl(pat, bad_paths, ignore.case = TRUE)
}))

if (any(bad_hit)) {
  stop("A forbidden old/Michael/batch-aware path appears in the corrected DAM workflow.")
}

clean_gene_vector <- function(x) {
  x %>%
    as.character() %>%
    stringr::str_trim() %>%
    na_if("") %>%
    purrr::discard(is.na) %>%
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
    "gene_name_final",
    "GeneID", "Gene_ID"
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
  df <- readr::read_csv(file, show_col_types = FALSE)
  gene_col <- guess_gene_column(df)
  genes <- clean_gene_vector(df[[gene_col]])

  cat("\nReading gene list file:\n", file, "\n", sep = "")
  cat("Using gene column: ", gene_col, "\n", sep = "")
  cat("Number of genes: ", length(genes), "\n", sep = "")

  genes
}

timepoint_gene_lists <- list()

for (i in seq_len(nrow(timepoint_files))) {
  nm <- timepoint_files$List_name[i]
  timepoint_gene_lists[[nm]] <- read_gene_list_csv(timepoint_files$File[i])
}

timepoint_size_check <- timepoint_files %>%
  mutate(
    n_genes = purrr::map_int(List_name, ~ length(timepoint_gene_lists[[.x]])),
    contains_TREM2 = purrr::map_lgl(List_name, ~ "TREM2" %in% timepoint_gene_lists[[.x]])
  )

write_csv(
  timepoint_size_check,
  file.path(out_dir, "CHECK_corrected_timepoint_gene_list_sizes.csv")
)

cat("\nCorrected timepoint gene list check:\n")
print(timepoint_size_check)

# Handles:
#   A) DAM.xlsx has multiple sheets, each sheet = one gene list
#   B) DAM.xlsx has one sheet with multiple gene-list columns

dam_sheets <- readxl::excel_sheets(dam_file)

if (length(dam_sheets) >= 3) {

  dam_gene_lists <- purrr::map(dam_sheets, function(sh) {
    df <- readxl::read_excel(dam_file, sheet = sh)
    df <- as.data.frame(df)
    names(df) <- make.unique(names(df))

    gene_col <- guess_gene_column(df)
    genes <- clean_gene_vector(df[[gene_col]])

    cat("\nReading DAM sheet: ", sh, "\n", sep = "")
    cat("Using gene column: ", gene_col, "\n", sep = "")
    cat("Number of genes: ", length(genes), "\n", sep = "")

    genes
  })

  names(dam_gene_lists) <- make.names(dam_sheets)

} else {

  dam_df <- readxl::read_excel(dam_file, sheet = dam_sheets[1])
  dam_df <- as.data.frame(dam_df)
  names(dam_df) <- make.unique(names(dam_df))

  gene_list_cols <- colnames(dam_df)[!sapply(dam_df, is.numeric)]

  if (length(gene_list_cols) < 1) {
    stop("DAM.xlsx does not contain gene-list-like columns. Please check the file.")
  }

  dam_gene_lists <- purrr::map(gene_list_cols, function(col) {
    genes <- clean_gene_vector(dam_df[[col]])

    cat("\nReading DAM column: ", col, "\n", sep = "")
    cat("Number of genes: ", length(genes), "\n", sep = "")

    genes
  })

  names(dam_gene_lists) <- make.names(gene_list_cols)
}

dam_size_check <- tibble(
  DAM_list = names(dam_gene_lists),
  n_genes = purrr::map_int(dam_gene_lists, length),
  contains_TREM2 = purrr::map_lgl(dam_gene_lists, ~ "TREM2" %in% .x)
)

write_csv(
  dam_size_check,
  file.path(out_dir, "CHECK_DAM_gene_list_sizes.csv")
)

cat("\nDAM gene list check:\n")
print(dam_size_check)

# 6. TREM2 presence check

all_lists_for_trem2 <- c(timepoint_gene_lists, dam_gene_lists)

trem2_check <- tibble(
  list_name = names(all_lists_for_trem2),
  n_genes = purrr::map_int(all_lists_for_trem2, length),
  contains_TREM2 = purrr::map_lgl(all_lists_for_trem2, ~ "TREM2" %in% .x)
)

write_csv(trem2_check, trem2_check_csv)

cat("\nTREM2 presence check:\n")
print(trem2_check)

# 7. ALL timepoints UP/DOWN vs DAM overlap

summary_results <- list()
overlap_results <- list()
counter <- 1

for (i in seq_len(nrow(timepoint_files))) {

  list_name <- timepoint_files$List_name[i]
  timepoint_name <- timepoint_files$Timepoint[i]
  direction_name <- timepoint_files$Direction[i]
  query_genes <- timepoint_gene_lists[[list_name]]

  for (dam_name in names(dam_gene_lists)) {

    dam_genes <- dam_gene_lists[[dam_name]]
    overlap_genes <- intersect(query_genes, dam_genes)

    summary_results[[counter]] <- tibble(
      Timepoint = timepoint_name,
      Direction = direction_name,
      Query_list = list_name,
      DAM_list = dam_name,
      Query_n = length(query_genes),
      DAM_n = length(dam_genes),
      overlap_n = length(overlap_genes),
      Query_percent_overlap = ifelse(
        length(query_genes) > 0,
        length(overlap_genes) / length(query_genes) * 100,
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
        Timepoint = timepoint_name,
        Direction = direction_name,
        Query_list = list_name,
        DAM_list = dam_name,
        overlap_gene = overlap_genes
      )
    } else {
      overlap_results[[counter]] <- tibble(
        Timepoint = character(0),
        Direction = character(0),
        Query_list = character(0),
        DAM_list = character(0),
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

# 8. Venn PDF

pdf(venn_pdf, width = 7.8, height = 6.8)

for (i in seq_len(nrow(timepoint_files))) {

  list_name <- timepoint_files$List_name[i]
  timepoint_name <- timepoint_files$Timepoint[i]
  direction_name <- timepoint_files$Direction[i]
  query_genes <- timepoint_gene_lists[[list_name]]

  for (dam_name in names(dam_gene_lists)) {

    dam_genes <- dam_gene_lists[[dam_name]]
    overlap_genes <- intersect(query_genes, dam_genes)
    overlap_n <- length(overlap_genes)

    grid.newpage()

    venn_plot <- draw.pairwise.venn(
      area1 = length(query_genes),
      area2 = length(dam_genes),
      cross.area = overlap_n,
      category = c(list_name, dam_name),
      fill = c("#9ecae1", "#fdae6b"),
      alpha = c(0.55, 0.55),
      lty = "blank",
      cex = 1.45,
      cat.cex = 0.80,
      cat.pos = c(-20, 20),
      cat.dist = c(0.05, 0.05),
      scaled = FALSE
    )

    grid.draw(venn_plot)

    grid.text(
      paste0(
        timepoint_name, " ", direction_name, " vs ", dam_name,
        "\nOverlap genes = ", overlap_n,
        "\nQuery n = ", length(query_genes),
        " | DAM n = ", length(dam_genes)
      ),
      x = 0.5,
      y = 0.96,
      gp = gpar(fontsize = 11, fontface = "bold")
    )

    if (overlap_n > 0) {
      overlap_label <- paste(overlap_genes, collapse = ", ")
      if (nchar(overlap_label) > 800) {
        overlap_label <- paste0(substr(overlap_label, 1, 800), " ...")
      }

      grid.text(
        paste0("Overlap: ", overlap_label),
        x = 0.5,
        y = 0.05,
        gp = gpar(fontsize = 7.5)
      )
    }
  }
}

dev.off()

# 9. Compact overlap-count heatmap PDF

summary_plot_df <- summary_df %>%
  mutate(
    Timepoint = factor(Timepoint, levels = c("3 weeks", "4 weeks", "Early onset", "Terminal")),
    Direction = factor(Direction, levels = c("UP", "DOWN")),
    Plot_row = paste(Timepoint, Direction, sep = " | ")
  )

p_heat <- ggplot(summary_plot_df, aes(x = DAM_list, y = Plot_row, fill = overlap_n)) +
  geom_tile(color = "white") +
  geom_text(aes(label = overlap_n), size = 3) +
  scale_fill_gradient(low = "white", high = "firebrick", name = "Overlap n") +
  labs(
    title = "Corrected LSD1 timepoint UP/DOWN gene overlap with DAM lists",
    x = "DAM gene list",
    y = "Corrected LSD1 timepoint list"
  ) +
  theme_bw(base_size = 11) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.text.x = element_text(angle = 35, hjust = 1),
    panel.grid = element_blank()
  )

pdf(heatmap_pdf, width = 12, height = 7)
print(p_heat)
dev.off()

# 10. Final output

cat("\nDone.\n")
cat("Venn PDF saved to:\n", venn_pdf, "\n\n")
cat("Overlap-count heatmap saved to:\n", heatmap_pdf, "\n\n")
cat("Summary saved to:\n", summary_csv, "\n\n")
cat("Overlap genes saved to:\n", overlap_csv, "\n\n")
cat("TREM2 check saved to:\n", trem2_check_csv, "\n\n")
cat("Input check saved to:\n", input_check_csv, "\n\n")
