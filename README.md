# LSD1 neurodegeneration downstream analysis

This repository contains curated R scripts for downstream transcriptomic analyses used in an experimental neurodegeneration project involving LSD1, Trem2, PLX treatment, and disease-stage comparisons.

The goal of this repository is computational transparency for the downstream analysis, not one-click reproduction from raw FASTQ files. Raw sequencing data, alignment, and count generation were handled outside this repository. The scripts here start from processed count matrices, differential-expression tables, gene lists, and annotation files.

## Repository structure

```text
LSD1_neurodegeneration_analysis/
├── README.md
├── scripts/
│   ├── 01_trem2_terminal_volcano_overlay.R
│   ├── 02_trem2_terminal_microglia_immune_heatmap.R
│   ├── 03_trem2_terminal_reactive_astrocyte_heatmap.R
│   ├── 04_trem2_terminal_up_enrichment.R
│   ├── 05_trem2_terminal_down_enrichment.R
│   ├── 06_trem2_terminal_dam_overlap.R
│   ├── 07_plx_terminal_combined_volcano_overlay.R
│   ├── 08_plx_terminal_microglia_immune_heatmap_by_sex.R
│   ├── 09_plx_terminal_reactive_astrocyte_heatmap_by_sex.R
│   ├── 10_plx_terminal_combined_enrichment.R
│   ├── 11_lsd1_timecourse_corrected_terminal_volcano_no_michael.R
│   ├── 12_lsd1_timecourse_corrected_terminal_enrichment_no_michael.R
│   ├── 13_lsd1_timecourse_corrected_terminal_immune_heatmap_no_michael.R
│   ├── 14_lsd1_timecourse_corrected_terminal_reactive_astrocyte_heatmap_no_michael.R
│   ├── 15_lsd1_timecourse_early_terminal_overlap_no_michael.R
│   └── 16_lsd1_timecourse_dam_overlap_no_michael.R
├── data/
│   ├── README.md
│   ├── processed/
│   └── external/
├── results/
└── docs/
    └── source_script_manifest.csv
```

## Analysis scope

The curated scripts cover three connected downstream analysis modules:

1. **Trem2 terminal experiment**
   - Differential-expression volcano plots
   - Microglia/immune-response heatmaps
   - Reactive-astrocyte heatmaps
   - UP/DOWN gene-list enrichment
   - DAM gene-list overlap analysis

2. **PLX terminal experiment**
   - Combined-sex volcano plots
   - Sex-stratified immune heatmaps
   - Sex-stratified reactive-astrocyte heatmaps
   - UP/DOWN gene-list enrichment

3. **Corrected LSD1 time-course analysis**
   - Corrected terminal-stage comparison excluding Michael terminal samples
   - DESeq2-based volcano workflow
   - Enrichment analysis
   - Immune and reactive-astrocyte heatmaps
   - Early-onset versus terminal overlap
   - DAM overlap across time points

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
Rscript scripts/11_lsd1_timecourse_corrected_terminal_volcano_no_michael.R
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
