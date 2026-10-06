# LSD1 neurodegeneration downstream analysis

This repository contains curated R scripts for downstream transcriptomic analyses used in an experimental neurodegeneration project involving LSD1, Trem2, PLX treatment, and disease-stage comparisons.

The goal of this repository is computational transparency for the downstream analysis, not one-click reproduction from raw FASTQ files. Raw sequencing data, alignment, and count generation were handled outside this repository. The scripts here start from processed count matrices, differential-expression tables, gene lists, and annotation files.

## Repository structure

```text
LSD1_neurodegeneration_analysis/
├── README.md
├── scripts/
│   ├── 01_lsd1_timecourse_terminal_volcano.R
│   ├── 02_lsd1_timecourse_terminal_immune_heatmap.R
│   ├── 03_lsd1_timecourse_terminal_reactive_astrocyte_heatmap.R
│   ├── 04_lsd1_timecourse_terminal_enrichment.R
│   ├── 05_lsd1_timecourse_early_terminal_overlap.R
│   ├── 06_lsd1_timecourse_dam_overlap.R
│   ├── 07_trem2_terminal_volcano_overlay.R
│   ├── 08_trem2_terminal_microglia_immune_heatmap.R
│   ├── 09_trem2_terminal_reactive_astrocyte_heatmap.R
│   ├── 10_trem2_terminal_up_enrichment.R
│   ├── 11_trem2_terminal_down_enrichment.R
│   ├── 12_trem2_terminal_dam_overlap.R
│   ├── 13_plx_terminal_combined_volcano_overlay.R
│   ├── 14_plx_terminal_microglia_immune_heatmap_by_sex.R
│   ├── 15_plx_terminal_reactive_astrocyte_heatmap_by_sex.R
│   └── 16_plx_terminal_combined_enrichment.R
├── data/
│   ├── README.md
│   ├── processed/
│   └── external/
├── results/
└── docs/
    └── source_script_manifest.csv
```

## Analysis scope

The curated scripts follow the analysis order used in the project narrative:

1. **LSD1 time-course analysis**
   - Terminal-stage volcano workflow across disease-stage datasets
   - Immune-response heatmaps across 3 weeks, 4 weeks, early-onset, and terminal stages
   - Reactive-astrocyte heatmaps across disease stages
   - Enrichment analysis from time-course gene lists
   - Early-onset versus terminal overlap analysis
   - DAM overlap across time points

2. **Trem2 terminal analysis**
   - Terminal-stage differential-expression volcano plots
   - Microglia/immune-response heatmaps
   - Reactive-astrocyte heatmaps
   - UP/DOWN gene-list enrichment
   - DAM gene-list overlap analysis

3. **PLX terminal analysis**
   - Combined-sex volcano plots
   - Sex-stratified immune heatmaps
   - Sex-stratified reactive-astrocyte heatmaps
   - UP/DOWN gene-list enrichment

## Reproducibility boundary

This repository is intended to reproduce the downstream statistical analysis and figure-generation logic when the required processed input files are available.

Expected inputs include:

- Processed featureCounts or count-matrix Excel files
- Cleaned DESeq2 comparison tables
- UP/DOWN gene-list CSV files
- Module annotation files such as `modules.csv`
- Gene annotation/module mapping files such as `NIHMS472534-supplement-02.csv`
- External gene-set resources such as CellMarker, DAM gene lists, or GMT files

Large raw sequencing files, BAM files, Galaxy histories, and collaborator-specific intermediate folders should not be committed to this repository. If this repository is used for publication, raw and processed sequencing data should be linked through the appropriate data accession, such as GEO or SRA.

## Running the scripts

Each script is designed to be run from the repository root or by `Rscript` using project-relative paths.

Example:

```bash
Rscript scripts/01_lsd1_timecourse_terminal_volcano.R
```

Before running a script, place the required processed input files into `data/processed/` or `data/external/` following the notes in `data/README.md`. Outputs should be written to `results/` or to script-defined analysis subfolders.

## R packages

The scripts use common R and Bioconductor packages including:

- `dplyr`
- `tidyr`
- `readr`
- `readxl`
- `ggplot2`
- `ggrepel`
- `patchwork`
- `ComplexHeatmap`
- `circlize`
- `DESeq2`
- `clusterProfiler`
- `org.Mm.eg.db`
- `ReactomePA`

For a formal publication archive, add a package-version file generated from `sessionInfo()` or `renv::snapshot()` after confirming the final execution environment.

## Code availability statement template

All custom R scripts used for downstream RNA-seq analysis, gene-list overlap analysis, enrichment analysis, and figure generation are available in this repository. Raw sequencing data and upstream alignment/count-generation workflows are documented separately and should be accessed through the associated data accession when available.
