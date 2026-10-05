# Data directory

This repository expects processed analysis inputs rather than raw sequencing files.

## `data/processed/`

Place project-derived processed inputs here, such as:

- FeatureCounts or count-matrix Excel files
- Cleaned DESeq2 result tables
- UP/DOWN gene-list CSV files
- Sample metadata tables
- Module/gene-order tables created from final analysis workflows

Example expected files include:

- `Trem2heatmap feature counts.xlsx`
- `3 weeks feature counts.xlsx`
- `4 weeks feature counts.xlsx`
- `early_onset_featurecounts.xlsx`
- `Terminal feature counts.xlsx`
- `modules.csv`
- `NIHMS472534-supplement-02.csv`
- `genotypes.xlsx`

## `data/external/`

Place external annotation or reference gene-set files here, such as:

- CellMarker tables
- DAM gene lists
- GMT files
- Public module annotation tables

## Files intentionally excluded from GitHub

Do not commit raw FASTQ files, BAM files, Galaxy histories, large intermediate alignment files, or private collaborator folders. For publication, link raw data through GEO/SRA or another approved public repository and include only the processed inputs required for downstream analysis reproduction.
