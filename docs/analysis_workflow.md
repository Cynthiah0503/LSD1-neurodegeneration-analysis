# Analysis workflow

## Overview

The downstream analysis starts from processed RNA-seq count matrices or cleaned differential-expression tables. The scripts then generate statistical summaries, gene lists, enrichment results, overlap analyses, heatmaps, and volcano plots used for biological interpretation.

## Workflow modules

### 1. LSD1 time-course analysis

1. Analyze gene-expression changes across 3 weeks, 4 weeks, early-onset paralysis, and terminal paralysis.
2. Use terminal-stage integration for cross-dataset terminal comparisons.
3. Run DESeq2 contrasts for LSD1KO versus LSD1WT.
4. Export UP/DOWN gene lists and volcano plots.
5. Generate immune and reactive-astrocyte heatmaps.
6. Run enrichment analysis from time-course gene lists.
7. Compare early-onset and terminal gene lists.
8. Compare time-course gene lists with DAM-associated genes.

### 2. Trem2 terminal analysis

1. Generate terminal-stage volcano plots with immune-module overlay.
2. Build microglia/immune-response heatmaps from featureCounts matrices.
3. Build reactive-astrocyte marker heatmaps.
4. Run UP- and DOWN-regulated gene-list enrichment.
5. Compare differential-expression gene lists with DAM-associated genes.

### 3. PLX terminal analysis

1. Generate combined-sex volcano plots.
2. Build sex-stratified microglia/immune-response heatmaps.
3. Build sex-stratified reactive-astrocyte heatmaps.
4. Run combined-sex UP/DOWN enrichment analysis.

## Interpretation notes

Positive log2 fold change in the LSD1 time-course volcano workflow represents higher expression in LSD1KO relative to LSD1WT for the specified contrast. Heatmap group means are calculated from individual-sample row-z scores where specified in the scripts, so the heatmaps emphasize relative expression patterns across groups rather than raw count magnitude.
