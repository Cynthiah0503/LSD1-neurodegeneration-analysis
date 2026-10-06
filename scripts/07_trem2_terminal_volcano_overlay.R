# 07_trem2_terminal_volcano_overlay.R
#
# Purpose: Trem2/LSD1 terminal differential-expression volcano plots with immune-module overlay.
# Inputs are expected under data/processed/ or data/external/ relative to this repository.
# Outputs are written under results/ or script-defined subfolders.

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
library(dplyr)
library(ggplot2)
library(ggrepel)
library(patchwork)

legend_no_title <- list(
  guides(
    color = guide_legend(
      title = NULL,
      override.aes = list(label = "")
    ),
    fill  = "none",
    shape = "none",
    size  = "none",
    alpha = "none"
  ),
  theme(legend.title = element_blank())
)

# Load cleaned differential-expression tables.
f1 <- read.csv(file.path(data_processed_dir, "CLEANED_Trem2_WT; LSD1_Del vs. Trem2_WT; LSD1_WT.csv"), header = TRUE)
f2 <- read.csv(file.path(data_processed_dir, "CLEANED_Trem2KO; LDS1_Del vs. Trem2_WT; LSD1_Del.csv"), header = TRUE)
f3 <- read.csv(file.path(data_processed_dir, "CLEANED_Trem2KO; LDS1_Del vs. Trem2_WT; LSD1_WT.csv"), header = TRUE)
f4 <- read.csv(file.path(data_processed_dir, "CLEANED_Trem2KO; LSD1_WT vs. Trem2_WT; LSD1_WT.csv"), header = TRUE)

# Define comparison titles.
title_f1 <- "Trem2WT_LSD1KO_vs_Trem2WT_LSD1WT"
title_f2 <- "Trem2KO_LSD1KO_vs_Trem2WT_LSD1KO"
title_f3 <- "Trem2KO_LSD1KO_vs_Trem2WT_LSD1WT"
title_f4 <- "Trem2KO_LSD1WT_vs_Trem2WT_LSD1WT"

# Set the left-to-right plot order.
# f1, f4, f3, f2

# 1) Trem2WT_LSD1KO vs Trem2WT_LSD1WT
f1 <- f1 %>%
  mutate(Significance = case_when(
    P.adj < 0.05 & log2.FC. > 1  ~ "Up regulated",
    P.adj < 0.05 & log2.FC. < -1 ~ "Down regulated",
    TRUE ~ "Not Significant"
  ))

top_genes_f1 <- f1 %>%
  filter(P.adj < 0.05, abs(log2.FC.) > 1) %>%
  arrange(P.adj) %>%
  slice_head(n = 15)

p_f1 <- ggplot(f1, aes(x = log2.FC., y = -log10(P.adj), color = Significance)) +
  geom_point(alpha = 0.65, size = 1.6, show.legend = TRUE) +
  geom_text_repel(
    data = top_genes_f1,
    aes(label = Gene_name),
    size = 4,
    max.overlaps = Inf,
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
  coord_cartesian(xlim = c(-10, 10), ylim = c(-2, 30)) +
  labs(
    title = title_f1,
    x = "log2(Fold Change)",
    y = "-log10(P.adj)"
  ) +
  theme_minimal(base_size = 14) +
  legend_no_title

# 2) Trem2KO_LSD1KO vs Trem2WT_LSD1KO
f2 <- f2 %>%
  mutate(Significance = case_when(
    P.adj < 0.05 & log2.FC. > 1  ~ "Up regulated",
    P.adj < 0.05 & log2.FC. < -1 ~ "Down regulated",
    TRUE ~ "Not Significant"
  ))

top_genes_f2 <- f2 %>%
  filter(P.adj < 0.05, abs(log2.FC.) > 1) %>%
  arrange(P.adj) %>%
  slice_head(n = 15)

p_f2 <- ggplot(f2, aes(x = log2.FC., y = -log10(P.adj), color = Significance)) +
  geom_point(alpha = 0.65, size = 1.6, show.legend = TRUE) +
  geom_text_repel(
    data = top_genes_f2,
    aes(label = Gene_name),
    size = 4,
    max.overlaps = Inf,
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
  coord_cartesian(xlim = c(-10, 10), ylim = c(-2, 30)) +
  labs(
    title = title_f2,
    x = "log2(Fold Change)",
    y = "-log10(P.adj)"
  ) +
  theme_minimal(base_size = 14) +
  legend_no_title

# 3) Trem2KO_LSD1KO vs Trem2WT_LSD1WT
f3 <- f3 %>%
  mutate(Significance = case_when(
    P.adj < 0.05 & log2.FC. > 1  ~ "Up regulated",
    P.adj < 0.05 & log2.FC. < -1 ~ "Down regulated",
    TRUE ~ "Not Significant"
  ))

top_genes_f3 <- f3 %>%
  filter(P.adj < 0.05, abs(log2.FC.) > 1) %>%
  arrange(P.adj) %>%
  slice_head(n = 15)

p_f3 <- ggplot(f3, aes(x = log2.FC., y = -log10(P.adj), color = Significance)) +
  geom_point(alpha = 0.65, size = 1.6, show.legend = TRUE) +
  geom_text_repel(
    data = top_genes_f3,
    aes(label = Gene_name),
    size = 4,
    max.overlaps = Inf,
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
  coord_cartesian(xlim = c(-10, 10), ylim = c(-2, 30)) +
  labs(
    title = title_f3,
    x = "log2(Fold Change)",
    y = "-log10(P.adj)"
  ) +
  theme_minimal(base_size = 14) +
  legend_no_title

# 4) Trem2KO_LSD1WT vs Trem2WT_LSD1WT
f4 <- f4 %>%
  mutate(Significance = case_when(
    P.adj < 0.05 & log2.FC. > 1  ~ "Up regulated",
    P.adj < 0.05 & log2.FC. < -1 ~ "Down regulated",
    TRUE ~ "Not Significant"
  ))

top_genes_f4 <- f4 %>%
  filter(P.adj < 0.05, abs(log2.FC.) > 1) %>%
  arrange(P.adj) %>%
  slice_head(n = 15)

p_f4 <- ggplot(f4, aes(x = log2.FC., y = -log10(P.adj), color = Significance)) +
  geom_point(alpha = 0.65, size = 1.6, show.legend = TRUE) +
  geom_text_repel(
    data = top_genes_f4,
    aes(label = Gene_name),
    size = 4,
    max.overlaps = Inf,
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
  coord_cartesian(xlim = c(-10, 10), ylim = c(-2, 30)) +
  labs(
    title = title_f4,
    x = "log2(Fold Change)",
    y = "-log10(P.adj)"
  ) +
  theme_minimal(base_size = 14) +
  legend_no_title

# Page 1: genome-wide gene expression changes.
# 1×4 + 一个 legend
# 左到右顺序：
# f1, f4, f3, f2
combined_4 <- (p_f1 | p_f4 | p_f3 | p_f2) +
  plot_layout(guides = "collect", ncol = 4) &
  theme(
    legend.position = "right",
    legend.title = element_blank(),
    plot.title = element_text(hjust = 0.5, face = "bold")
  ) &
  guides(
    color = guide_legend(
      title = NULL,
      override.aes = list(shape = 16, size = 3, alpha = 1, label = "")
    )
  )

# Gold / Forestgreen module gene sets

allc <- read.csv(file.path(data_external_dir, "NIHMS472534-supplement-02.csv"), header = TRUE, check.names = FALSE)
mod  <- read.csv(file.path(data_external_dir, "modules.csv"), header = TRUE, check.names = FALSE)

names(mod) <- sub("^\\ufeff", "", names(mod))
bad <- is.na(names(mod)) | names(mod) == ""
if (any(bad)) names(mod)[bad] <- paste0("V", seq_len(sum(bad)))
names(mod) <- make.unique(names(mod))

names(allc) <- sub("^\\ufeff", "", names(allc))
allc$Gene_Symbol <- toupper(as.character(allc$Gene_Symbol))

gene2module <- allc %>%
  select(Gene_Symbol, Module) %>%
  filter(!is.na(Gene_Symbol), !is.na(Module)) %>%
  mutate(
    Gene_Symbol = toupper(as.character(Gene_Symbol)),
    Module_clean = tolower(trimws(as.character(Module)))
  ) %>%
  distinct()

all_genes_universe <- unique(toupper(c(f1$Gene_name, f2$Gene_name, f3$Gene_name, f4$Gene_name)))

gold_set <- gene2module %>%
  filter(Module_clean == "gold") %>%
  pull(Gene_Symbol) %>%
  unique()

forestgreen_set <- gene2module %>%
  filter(Module_clean == "forestgreen") %>%
  pull(Gene_Symbol) %>%
  unique()

gold_set <- intersect(gold_set, all_genes_universe)
forestgreen_set <- intersect(forestgreen_set, all_genes_universe)

if (length(gold_set) == 0) stop("Gold module gene set is empty after intersecting with volcano tables.")
if (length(forestgreen_set) == 0) stop("Forestgreen module gene set is empty after intersecting with volcano tables.")

cat("Gold module genes found in volcano tables:", length(gold_set), "\\n")
cat("Forestgreen module genes found in volcano tables:", length(forestgreen_set), "\\n")

# Page 2: immune-module overlay.

make_volcano_with_gold_forestgreen <- function(df, top_n = 15, plot_title, xlab_text, subtitle_text = NULL,
                                               gold_genes = NULL,
                                               forestgreen_genes = NULL,
                                               base_point_size = 1.6,
                                               overlay_point_size = 1.5) {

  stopifnot(all(c("Gene_name","P.adj","log2.FC.") %in% names(df)))

  df2 <- df %>%
    mutate(
      Gene_UP = toupper(as.character(Gene_name)),
      Significance = case_when(
        P.adj < 0.05 & log2.FC. > 1  ~ "Up regulated",
        P.adj < 0.05 & log2.FC. < -1 ~ "Down regulated",
        TRUE ~ "Not Significant"
      ),
      Overlay_Module = case_when(
        Gene_UP %in% toupper(gold_genes) ~ "Gold",
        Gene_UP %in% toupper(forestgreen_genes) ~ "Forestgreen",
        TRUE ~ NA_character_
      )
    )

  top_genes <- df2 %>%
    filter(P.adj < 0.05, abs(log2.FC.) > 1) %>%
    arrange(P.adj) %>%
    slice_head(n = top_n)

  ggplot(df2, aes(x = log2.FC., y = -log10(P.adj))) +
    geom_point(aes(color = Significance), alpha = 0.65, size = base_point_size, show.legend = TRUE) +
    geom_point(
      data = df2 %>% filter(Overlay_Module == "Gold"),
      aes(x = log2.FC., y = -log10(P.adj)),
      shape = 21, fill = "gold", color = "black",
      stroke = 0.35, size = overlay_point_size, alpha = 0.95,
      show.legend = FALSE
    ) +
    geom_point(
      data = df2 %>% filter(Overlay_Module == "Forestgreen"),
      aes(x = log2.FC., y = -log10(P.adj)),
      shape = 21, fill = "forestgreen", color = "black",
      stroke = 0.35, size = overlay_point_size, alpha = 0.95,
      show.legend = FALSE
    ) +
    geom_text_repel(
      data = top_genes,
      aes(label = Gene_name),
      size = 4,
      max.overlaps = Inf,
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
    coord_cartesian(xlim = c(-10, 10), ylim = c(-2, 30)) +
    labs(
      title = plot_title,
      subtitle = subtitle_text,
      x = "log2(Fold Change)",
      y = "-log10(P.adj)"
    ) +
    theme_minimal(base_size = 14) +
    legend_no_title
}

p_gf1 <- make_volcano_with_gold_forestgreen(
  f1, top_n = 15,
  plot_title = title_f1,
  xlab_text  = "log2(Fold Change)",
  subtitle_text = "Gold and forestgreen module overlay",
  gold_genes = gold_set,
  forestgreen_genes = forestgreen_set
)

p_gf2 <- make_volcano_with_gold_forestgreen(
  f2, top_n = 15,
  plot_title = title_f2,
  xlab_text  = "log2(Fold Change)",
  subtitle_text = "Gold and forestgreen module overlay",
  gold_genes = gold_set,
  forestgreen_genes = forestgreen_set
)

p_gf3 <- make_volcano_with_gold_forestgreen(
  f3, top_n = 15,
  plot_title = title_f3,
  xlab_text  = "log2(Fold Change)",
  subtitle_text = "Gold and forestgreen module overlay",
  gold_genes = gold_set,
  forestgreen_genes = forestgreen_set
)

p_gf4 <- make_volcano_with_gold_forestgreen(
  f4, top_n = 15,
  plot_title = title_f4,
  xlab_text  = "log2(Fold Change)",
  subtitle_text = "Gold and forestgreen module overlay",
  gold_genes = gold_set,
  forestgreen_genes = forestgreen_set
)

combined_gold_forestgreen_overlay <- (p_gf1 | p_gf4 | p_gf3 | p_gf2) +
  plot_layout(guides = "collect", ncol = 4) &
  theme(
    legend.position = "right",
    legend.title = element_blank(),
    plot.title = element_text(hjust = 0.5, face = "bold")
  ) &
  guides(
    color = guide_legend(
      title = NULL,
      override.aes = list(shape = 16, size = 3, alpha = 1, label = "")
    )
  )

# FINAL: 合并成一个 PDF
# 现在只输出 2 页：
#   Page 1 = whole genome gene expression changes
#   Page 2 = gold + forestgreen overlay

graphics.off()

pdf(file.path(results_dir, "trem2_terminal_volcano_overlay.pdf"), width = 28, height = 8, onefile = TRUE)

print(combined_4)
print(combined_gold_forestgreen_overlay)

dev.off()

# Export DOWN genes
# Export significant DOWN-regulated gene lists for downstream enrichment.
get_down_table <- function(df, comparison_name) {
  df %>%
    filter(P.adj < 0.05, log2.FC. < -1) %>%
    arrange(P.adj) %>%
    mutate(Comparison = comparison_name)
}

down_f1 <- get_down_table(f1, "Trem2WT_LSD1KO_vs_Trem2WT_LSD1WT")
down_f2 <- get_down_table(f2, "Trem2KO_LSD1KO_vs_Trem2WT_LSD1KO")
down_f3 <- get_down_table(f3, "Trem2KO_LSD1KO_vs_Trem2WT_LSD1WT")
down_f4 <- get_down_table(f4, "Trem2KO_LSD1WT_vs_Trem2WT_LSD1WT")

write.csv(down_f1, file.path(results_dir, "DOWN_genes_Trem2WT_LSD1KO_vs_Trem2WT_LSD1WT.csv"), row.names = FALSE)
write.csv(down_f2, file.path(results_dir, "DOWN_genes_Trem2KO_LSD1KO_vs_Trem2WT_LSD1KO.csv"), row.names = FALSE)
write.csv(down_f3, file.path(results_dir, "DOWN_genes_Trem2KO_LSD1KO_vs_Trem2WT_LSD1WT.csv"), row.names = FALSE)
write.csv(down_f4, file.path(results_dir, "DOWN_genes_Trem2KO_LSD1WT_vs_Trem2WT_LSD1WT.csv"), row.names = FALSE)

cat("Downregulated genes exported:\\n")
cat("f1:", nrow(down_f1), "genes\\n")
cat("f2:", nrow(down_f2), "genes\\n")
cat("f3:", nrow(down_f3), "genes\\n")
cat("f4:", nrow(down_f4), "genes\\n")
