# sc-cellranger-workflow

## Workflow to perform single-cell RNA-seq analysis

This repo will contain a Cell Ranger workflow for going from BCL/FASTQ/10x files to fully analyzed single-cell data (incomplete). More info will be available soon.

- `daemon.sh`: Daemon file for controlling the whole workflow (incomplete)

- `config_daemon`: File with default and global variables (incomplete)

- `manage_input_files.sh`: For creating soft links to the input data. It must be edited manually

- `samples_to_process.lst`: File with samples sames, one per row

- `R/` -> R libraries
    - `preprocessing_library.R`: R functions for performing pre-processing of scRNA-seq feature-barcode matrices

- `autoflow/` -> [AutoFlow](https://github.com/seoanezonjic/autoflow) templates
    - `QC_template.txt`: Template for Quality Control with FastQC
    - `count_template.txt`: Template for converting FASTQ files to feature-barcode matrices using Cell Ranger `count`
    - `preprocessing_template.txt`: Template for pre-processing analysis with `Seurat`

- `aux_sh/` -> Bash scripts for each stage of the workflow
    - `autoflow_launcher.sh`: Generic [AutoFlow](https://github.com/seoanezonjic/autoflow) launcher for templates in `autoflow/`
    - `cellranger_count.sh`: **[Deprecated]** Bash script for stage 1: converting FASTQ files to a proper format for single-cell analysis using Cell Ranger `count`
    - `cellranger_mkfastq.sh`: Bash script for stage 0: converting BCL files to FASTQ files using Cell Ranger `mkfastq`
    - `compare_samples.sh`: Bash script for stage 2b: samples comparison
    - `preprocess.sh`: **[Deprecated]** Bash script for stage 2: preprocessing using R package `Seurat`

- `scripts/` -> Scripts for some stages of the workflow
    - `count_table.R`: **[On pause]** R script for generating a table with basic FASTQ metrics
    - `preprocessing.R`: R script for applying `Seurat` preprocessing and creating a report
    - `raw_preprocessing.R`: **[Deprecated]** Old file with some non organized pre-processing code

- `templates/` -> R Markdown templates for reports
    - `mapping_report.Rmd`: **[Imported from DEG workflow]** Template for mapping report after using [SeqTrimBB](https://github.com/rafnunser/seqtrimbb)
    - `preprocessing_report.Rmd`: Template for pre-processing report
    - `test.Rmd`: **[Deprecated]** Test template
