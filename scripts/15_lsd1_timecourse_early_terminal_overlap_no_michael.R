# 15_lsd1_timecourse_early_terminal_overlap_no_michael.R
#
# Purpose: Corrected early-onset versus terminal overlap analysis excluding Michael terminal samples.
# Inputs are expected under data/processed/ or data/external/ relative to this repository.
# Outputs are written under results/ or script-defined subfolders.

# Early onset vs Terminal
# UPREGULATED GENES:
# Venn diagram + Hypergeometric test
# Put hypergeometric result directly on the Venn diagram


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
library(dplyr)
library(VennDiagram)
library(grid)
library(readxl)

# 防冲突
select  <- dplyr::select
filter  <- dplyr::filter
mutate  <- dplyr::mutate
arrange <- dplyr::arrange

# 手动排除的基因（大小写不敏感）
excluded_genes <- c("ESR1")

# File paths
input_dir <- "./timepoint_lsd1/正确combine terminal/CORRECTED_volcano_LSD1terminal_plus_Trem2WTterminal_no_Michael_GALAXY_style_FC"
output_dir <- file.path(input_dir, "VENN_Early_onset_vs_Terminal_CORRECTED_COMBINE_terminal_NO_ESR1")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# 1) Input files
# UP genes files
up_early_file    <- file.path(input_dir, "UP_genes_Early_onset_CORRECTED_COMBINE_terminal.csv")
up_terminal_file <- file.path(input_dir, "UP_genes_Terminal_CORRECTED_LSD1terminal_plus_Trem2WTterminal_NO_MICHAEL.csv")

# 原始 full DEG files（直接用它们做 universe，不再依赖 CLEANED 文件）
full_early_file    <- file.path(input_dir, "Early_onset_GALAXY_STYLE_DESeq2_KO_vs_WT.csv")
full_terminal_file <- file.path(input_dir, "Terminal_LSD1terminal_plus_Trem2WTterminal_NO_MICHAEL_GALAXY_STYLE_DESeq2_KO_vs_WT.csv")

read_any_table <- function(path) {
  if (grepl("\\.xlsx$", path, ignore.case = TRUE)) {
    df <- readxl::read_excel(path, sheet = 1)
  } else if (grepl("\\.csv$", path, ignore.case = TRUE)) {
    df <- read.csv(path, header = TRUE, check.names = FALSE)
  } else {
    stop(paste("Unsupported file format:", path))
  }
  as.data.frame(df)
}

# 3) Standardize columns
standardize_volcano_cols <- function(df) {
  names(df) <- sub("^\\ufeff", "", names(df))

  bad <- is.na(names(df)) | names(df) == ""
  if (any(bad)) {
    names(df)[bad] <- paste0("V", seq_len(sum(bad)))
  }

  names(df) <- make.unique(names(df))

  if (!("Gene_name" %in% names(df))) {
    gene_candidates <- c("Gene_name", "Gene", "gene", "Gene.Symbol", "Gene_Symbol", "SYMBOL")
    hit <- gene_candidates[gene_candidates %in% names(df)]
    if (length(hit) > 0) names(df)[names(df) == hit[1]] <- "Gene_name"
  }

  if (!("log2.FC." %in% names(df))) {
    fc_candidates <- c(
      "log2.FC.", "log2FC", "log2_FC", "log2FoldChange",
      "log2 fold change", "log2(FC)"
    )
    hit <- fc_candidates[fc_candidates %in% names(df)]
    if (length(hit) > 0) names(df)[names(df) == hit[1]] <- "log2.FC."
  }

  if (!("P.adj" %in% names(df))) {
    padj_candidates <- c(
      "P.adj", "padj", "adj.P.Val", "FDR", "p_adj",
      "adjusted_pvalue", "P-adj"
    )
    hit <- padj_candidates[padj_candidates %in% names(df)]
    if (length(hit) > 0) names(df)[names(df) == hit[1]] <- "P.adj"
  }

  needed <- c("Gene_name", "log2.FC.", "P.adj")
  missing_cols <- setdiff(needed, names(df))
  if (length(missing_cols) > 0) {
    stop(paste("Missing required columns:", paste(missing_cols, collapse = ", ")))
  }

  df %>%
    select(-any_of(c("Delete", "Delete2", "Delete3", "Delete4"))) %>%
    mutate(
      Gene_name = as.character(Gene_name),
      log2.FC.  = as.numeric(log2.FC.),
      P.adj     = as.numeric(P.adj)
    ) %>%
    filter(!is.na(Gene_name), !is.na(log2.FC.), !is.na(P.adj)) %>%
    distinct(Gene_name, .keep_all = TRUE) %>%
    arrange(Gene_name)
}

up_early <- read.csv(up_early_file, header = TRUE, stringsAsFactors = FALSE)
up_terminal <- read.csv(up_terminal_file, header = TRUE, stringsAsFactors = FALSE)

full_early <- read_any_table(full_early_file)
full_terminal <- read_any_table(full_terminal_file)

full_early <- standardize_volcano_cols(full_early)
full_terminal <- standardize_volcano_cols(full_terminal)

# 明确记录 Esr1 是否在原始输入中；后续所有集合与 universe 均会排除它
cat("Manual exclusion: ESR1/Esr1 is removed from Early onset, Terminal, and universe.\n")

# 5) Extract gene sets
up_early_genes <- up_early %>%
  filter(
    !is.na(Gene_name),
    Gene_name != "",
    !(toupper(trimws(Gene_name)) %in% excluded_genes)
  ) %>%
  pull(Gene_name) %>%
  unique() %>%
  sort()

up_terminal_genes <- up_terminal %>%
  filter(
    !is.na(Gene_name),
    Gene_name != "",
    !(toupper(trimws(Gene_name)) %in% excluded_genes)
  ) %>%
  pull(Gene_name) %>%
  unique() %>%
  sort()

# universe：两个 full DEG 表中所有检测到基因的并集
universe_genes <- union(
  full_early %>%
    filter(
    !is.na(Gene_name),
    Gene_name != "",
    !(toupper(trimws(Gene_name)) %in% excluded_genes)
  ) %>%
    pull(Gene_name) %>%
    unique(),
  full_terminal %>%
    filter(
    !is.na(Gene_name),
    Gene_name != "",
    !(toupper(trimws(Gene_name)) %in% excluded_genes)
  ) %>%
    pull(Gene_name) %>%
    unique()
) %>%
  sort()

# 6) Overlap / unique genes
overlap_genes <- intersect(up_early_genes, up_terminal_genes) %>% sort()
early_only_genes <- setdiff(up_early_genes, up_terminal_genes) %>% sort()
terminal_only_genes <- setdiff(up_terminal_genes, up_early_genes) %>% sort()

# 7) Print counts
cat("Upregulated genes in Early onset:", length(up_early_genes), "\n")
cat("Upregulated genes in Terminal:", length(up_terminal_genes), "\n")
cat("Overlap upregulated genes:", length(overlap_genes), "\n")
cat("Early-only upregulated genes:", length(early_only_genes), "\n")
cat("Terminal-only upregulated genes:", length(terminal_only_genes), "\n")
cat("Universe genes:", length(universe_genes), "\n\n")

# 8) Hypergeometric test
N <- length(universe_genes)
m <- length(up_early_genes)
k <- length(up_terminal_genes)
x <- length(overlap_genes)

expected_overlap <- (m * k) / N
fold_enrichment  <- x / expected_overlap
p_value <- phyper(q = x - 1, m = m, n = N - m, k = k, lower.tail = FALSE)

hyper_result <- data.frame(
  Comparison = "Early_onset_UP_vs_Terminal_UP",
  Universe_size = N,
  Early_onset_UP = m,
  Terminal_UP = k,
  Overlap_UP = x,
  Expected_overlap = expected_overlap,
  Fold_enrichment = fold_enrichment,
  Hypergeometric_p_value = p_value,
  stringsAsFactors = FALSE
)

print(hyper_result)

# write.csv(
#   data.frame(Gene_name = overlap_genes),
#   "UP_OVERLAP_Early_onset_vs_Terminal.csv",
#   row.names = FALSE
# )
#
# write.csv(
#   data.frame(Gene_name = early_only_genes),
#   "UP_ONLY_Early_onset_vs_Terminal.csv",
#   row.names = FALSE
# )
#
# write.csv(
#   data.frame(Gene_name = terminal_only_genes),
#   "UP_ONLY_Terminal_vs_Early_onset.csv",
#   row.names = FALSE
# )
#
# write.csv(
#   hyper_result,
#   "HYPERGEOMETRIC_TEST_Early_onset_vs_Terminal_UP.csv",
#   row.names = FALSE
# )

# 10) Format p-value text
p_label <- if (p_value < 2.2e-16) {
  "p < 2.2e-16"
} else if (p_value < 0.001) {
  paste0("p = ", format(p_value, scientific = TRUE, digits = 2))
} else {
  paste0("p = ", signif(p_value, 3))
}

expected_label <- paste0("Expected overlap = ", round(expected_overlap, 2))
fold_label     <- paste0("Fold enrichment = ", round(fold_enrichment, 2))
universe_label <- paste0("Universe = ", N)

# 11) Venn diagram
venn_plot <- venn.diagram(
  x = list(
    Early_onset = up_early_genes,
    Terminal = up_terminal_genes
  ),
  filename = NULL,
  fill = c("#6BAED6", "#F4A6A6"),
  alpha = 0.6,
  cex = 1.8,

  # 关掉系统自动放的组名
  category.names = c("", ""),

  lwd = 2,
  col = "gray40",
  main = "Overlap of Upregulated Genes\nEarly onset vs Terminal",
  main.cex = 1.5
)

pdf(file.path(output_dir, "VENN_UP_Early_onset_vs_Terminal_CORRECTED_COMBINE_terminal.pdf"), width = 8.5, height = 8.5)
grid.newpage()
grid.draw(venn_plot)

# 手动把两个组名放到你红线的位置（圈外、同一水平线）
grid.text(
  "Terminal",
  x = 0.12, y = 0.57,
  just = c("right", "center"),
  gp = gpar(fontsize = 16)
)

grid.text(
  "Early_onset",
  x = 0.85, y = 0.57,
  just = c("left", "center"),
  gp = gpar(fontsize = 16)
)

# hypergeometric 结果放下方空白处
grid.text(
  label = paste(
    "Hypergeometric test",
    p_label,
    # expected_label,
    # fold_label,
    # universe_label,
    sep = "\n"
  ),
  x = 0.84, y = 0.88,
  gp = gpar(fontsize = 11)
)

dev.off()

sink(file.path(output_dir, "SUMMARY_UP_Early_onset_vs_Terminal_CORRECTED_COMBINE_terminal.txt"))
cat("Comparison: Early onset vs Terminal (UP genes)\n\n")
cat("Upregulated genes in Early onset:", m, "\n")
cat("Upregulated genes in Terminal:", k, "\n")
cat("Overlap upregulated genes:", x, "\n")
cat("Universe genes:", N, "\n")
cat("Expected overlap:", expected_overlap, "\n")
cat("Fold enrichment:", fold_enrichment, "\n")
cat("Hypergeometric p-value:", p_value, "\n")
sink()

cat("Done.\n")

# Early onset vs Terminal
# DOWNREGULATED GENES:
# Venn diagram + Hypergeometric test
# Put hypergeometric result directly on the Venn diagram

library(dplyr)
library(VennDiagram)
library(grid)
library(readxl)

# 防冲突
select  <- dplyr::select
filter  <- dplyr::filter
mutate  <- dplyr::mutate
arrange <- dplyr::arrange

# 手动排除的基因（大小写不敏感）
excluded_genes <- c("ESR1")

# 1) Input files
# DOWN genes files
down_early_file    <- file.path(input_dir, "DOWN_genes_Early_onset_CORRECTED_COMBINE_terminal.csv")
down_terminal_file <- file.path(input_dir, "DOWN_genes_Terminal_CORRECTED_LSD1terminal_plus_Trem2WTterminal_NO_MICHAEL.csv")

# 原始 full DEG files（直接用它们做 universe，不依赖 CLEANED 文件）
full_early_file    <- file.path(input_dir, "Early_onset_GALAXY_STYLE_DESeq2_KO_vs_WT.csv")
full_terminal_file <- file.path(input_dir, "Terminal_LSD1terminal_plus_Trem2WTterminal_NO_MICHAEL_GALAXY_STYLE_DESeq2_KO_vs_WT.csv")

read_any_table <- function(path) {
  if (grepl("\\.xlsx$", path, ignore.case = TRUE)) {
    df <- readxl::read_excel(path, sheet = 1)
  } else if (grepl("\\.csv$", path, ignore.case = TRUE)) {
    df <- read.csv(path, header = TRUE, check.names = FALSE)
  } else {
    stop(paste("Unsupported file format:", path))
  }
  as.data.frame(df)
}

# 3) Standardize columns
standardize_volcano_cols <- function(df) {
  names(df) <- sub("^\\ufeff", "", names(df))

  bad <- is.na(names(df)) | names(df) == ""
  if (any(bad)) {
    names(df)[bad] <- paste0("V", seq_len(sum(bad)))
  }

  names(df) <- make.unique(names(df))

  if (!("Gene_name" %in% names(df))) {
    gene_candidates <- c("Gene_name", "Gene", "gene", "Gene.Symbol", "Gene_Symbol", "SYMBOL")
    hit <- gene_candidates[gene_candidates %in% names(df)]
    if (length(hit) > 0) names(df)[names(df) == hit[1]] <- "Gene_name"
  }

  if (!("log2.FC." %in% names(df))) {
    fc_candidates <- c(
      "log2.FC.", "log2FC", "log2_FC", "log2FoldChange",
      "log2 fold change", "log2(FC)"
    )
    hit <- fc_candidates[fc_candidates %in% names(df)]
    if (length(hit) > 0) names(df)[names(df) == hit[1]] <- "log2.FC."
  }

  if (!("P.adj" %in% names(df))) {
    padj_candidates <- c(
      "P.adj", "padj", "adj.P.Val", "FDR", "p_adj",
      "adjusted_pvalue", "P-adj"
    )
    hit <- padj_candidates[padj_candidates %in% names(df)]
    if (length(hit) > 0) names(df)[names(df) == hit[1]] <- "P.adj"
  }

  needed <- c("Gene_name", "log2.FC.", "P.adj")
  missing_cols <- setdiff(needed, names(df))
  if (length(missing_cols) > 0) {
    stop(paste("Missing required columns:", paste(missing_cols, collapse = ", ")))
  }

  df %>%
    select(-any_of(c("Delete", "Delete2", "Delete3", "Delete4"))) %>%
    mutate(
      Gene_name = as.character(Gene_name),
      log2.FC.  = as.numeric(log2.FC.),
      P.adj     = as.numeric(P.adj)
    ) %>%
    filter(!is.na(Gene_name), !is.na(log2.FC.), !is.na(P.adj)) %>%
    distinct(Gene_name, .keep_all = TRUE) %>%
    arrange(Gene_name)
}

down_early <- read.csv(down_early_file, header = TRUE, stringsAsFactors = FALSE)
down_terminal <- read.csv(down_terminal_file, header = TRUE, stringsAsFactors = FALSE)

full_early <- read_any_table(full_early_file)
full_terminal <- read_any_table(full_terminal_file)

full_early <- standardize_volcano_cols(full_early)
full_terminal <- standardize_volcano_cols(full_terminal)

# 明确记录 Esr1 是否在原始输入中；后续所有集合与 universe 均会排除它
cat("Manual exclusion: ESR1/Esr1 is removed from Early onset, Terminal, and universe.\n")

# 5) Extract gene sets
down_early_genes <- down_early %>%
  filter(
    !is.na(Gene_name),
    Gene_name != "",
    !(toupper(trimws(Gene_name)) %in% excluded_genes)
  ) %>%
  pull(Gene_name) %>%
  unique() %>%
  sort()

down_terminal_genes <- down_terminal %>%
  filter(
    !is.na(Gene_name),
    Gene_name != "",
    !(toupper(trimws(Gene_name)) %in% excluded_genes)
  ) %>%
  pull(Gene_name) %>%
  unique() %>%
  sort()

# universe：两个 full DEG 表中所有检测到基因的并集
universe_genes <- union(
  full_early %>%
    filter(
    !is.na(Gene_name),
    Gene_name != "",
    !(toupper(trimws(Gene_name)) %in% excluded_genes)
  ) %>%
    pull(Gene_name) %>%
    unique(),
  full_terminal %>%
    filter(
    !is.na(Gene_name),
    Gene_name != "",
    !(toupper(trimws(Gene_name)) %in% excluded_genes)
  ) %>%
    pull(Gene_name) %>%
    unique()
) %>%
  sort()

# 6) Overlap / unique genes
overlap_genes <- intersect(down_early_genes, down_terminal_genes) %>% sort()
early_only_genes <- setdiff(down_early_genes, down_terminal_genes) %>% sort()
terminal_only_genes <- setdiff(down_terminal_genes, down_early_genes) %>% sort()

# 7) Print counts
cat("Downregulated genes in Early onset:", length(down_early_genes), "\n")
cat("Downregulated genes in Terminal:", length(down_terminal_genes), "\n")
cat("Overlap downregulated genes:", length(overlap_genes), "\n")
cat("Early-only downregulated genes:", length(early_only_genes), "\n")
cat("Terminal-only downregulated genes:", length(terminal_only_genes), "\n")
cat("Universe genes:", length(universe_genes), "\n\n")

# 8) Hypergeometric test
N <- length(universe_genes)
m <- length(down_early_genes)
k <- length(down_terminal_genes)
x <- length(overlap_genes)

expected_overlap <- (m * k) / N
fold_enrichment  <- x / expected_overlap
p_value <- phyper(q = x - 1, m = m, n = N - m, k = k, lower.tail = FALSE)

hyper_result <- data.frame(
  Comparison = "Early_onset_DOWN_vs_Terminal_DOWN",
  Universe_size = N,
  Early_onset_DOWN = m,
  Terminal_DOWN = k,
  Overlap_DOWN = x,
  Expected_overlap = expected_overlap,
  Fold_enrichment = fold_enrichment,
  Hypergeometric_p_value = p_value,
  stringsAsFactors = FALSE
)

print(hyper_result)

# write.csv(
#   data.frame(Gene_name = overlap_genes),
#   "DOWN_OVERLAP_Early_onset_vs_Terminal.csv",
#   row.names = FALSE
# )
#
# write.csv(
#   data.frame(Gene_name = early_only_genes),
#   "DOWN_ONLY_Early_onset_vs_Terminal.csv",
#   row.names = FALSE
# )
#
# write.csv(
#   data.frame(Gene_name = terminal_only_genes),
#   "DOWN_ONLY_Terminal_vs_Early_onset.csv",
#   row.names = FALSE
# )
#
# write.csv(
#   hyper_result,
#   "HYPERGEOMETRIC_TEST_Early_onset_vs_Terminal_DOWN.csv",
#   row.names = FALSE
# )

# 10) Format p-value text
p_label <- if (p_value < 2.2e-16) {
  "p < 2.2e-16"
} else if (p_value < 0.001) {
  paste0("p = ", format(p_value, scientific = TRUE, digits = 2))
} else {
  paste0("p = ", signif(p_value, 3))
}

expected_label <- paste0("Expected overlap = ", round(expected_overlap, 2))
fold_label     <- paste0("Fold enrichment = ", round(fold_enrichment, 2))
universe_label <- paste0("Universe = ", N)

# 11) Venn diagram
venn_plot <- venn.diagram(
  x = list(
    Early_onset = down_early_genes,
    Terminal = down_terminal_genes
  ),
  filename = NULL,
  fill = c("#6BAED6", "#F4A6A6"),
  alpha = 0.6,
  cex = 1.8,

  # 关掉系统自动放的组名
  category.names = c("", ""),

  # 关键：不要把 overlap 数字拉到圈外并画引导线
  ext.text = FALSE,

  lwd = 2,
  col = "gray40",
  main = "Overlap of Downregulated Genes\nEarly onset vs Terminal",
  main.cex = 1.5
)

pdf(file.path(output_dir, "VENN_DOWN_Early_onset_vs_Terminal_CORRECTED_COMBINE_terminal.pdf"), width = 8.5, height = 8.5)
grid.newpage()
grid.draw(venn_plot)

# 手动把两个组名放到同样的位置
grid.text(
  "Terminal",
  x = 0.12, y = 0.57,
  just = c("right", "center"),
  gp = gpar(fontsize = 16)
)

grid.text(
  "Early_onset",
  x = 0.85, y = 0.57,
  just = c("left", "center"),
  gp = gpar(fontsize = 16)
)

# hypergeometric 结果放同样的位置
grid.text(
  label = paste(
    "Hypergeometric test",
    p_label,
    # expected_label,
    # fold_label,
    # universe_label,
    sep = "\n"
  ),
  x = 0.84, y = 0.88,
  gp = gpar(fontsize = 11)
)

dev.off()

sink(file.path(output_dir, "SUMMARY_DOWN_Early_onset_vs_Terminal_CORRECTED_COMBINE_terminal.txt"))
cat("Comparison: Early onset vs Terminal (DOWN genes)\n\n")
cat("Downregulated genes in Early onset:", m, "\n")
cat("Downregulated genes in Terminal:", k, "\n")
cat("Overlap downregulated genes:", x, "\n")
cat("Universe genes:", N, "\n")
cat("Expected overlap:", expected_overlap, "\n")
cat("Fold enrichment:", fold_enrichment, "\n")
cat("Hypergeometric p-value:", p_value, "\n")
sink()

cat("Done.\n")
