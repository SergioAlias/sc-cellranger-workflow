# sc-cellranger-workflow

## Workflow to perform single-cell RNA-seq analysis

This repo will contain a Cell Ranger workflow for going from BCL/FASTQ/10x files to fully analyzed single-cell data (incomplete). More info will be available soon.

### Workflow files and directories

- `daemon.sh`: Daemon file for controlling the whole workflow (incomplete)

- `config_daemon`: File with default and global variables (incomplete)

- `manage_input_files.sh`: For creating soft links to the input data. It must be edited manually

- `samples_to_process.lst`: File with samples sames, one per row

- `R/` -> R libraries
    - `preprocessing_library.R`: R functions for performing pre-processing of scRNA-seq feature-barcode matrices (stages 3a and 3b)

    - `qc_library.R`: R functions for generating the QC report (stage 2b)

- `autoflow/` -> [AutoFlow](https://github.com/seoanezonjic/autoflow) templates
    - `QC_template.txt`: Template for Quality Control with FastQC

    - `count_template.txt`: Template for converting FASTQ files to feature-barcode matrices using Cell Ranger `count`

    - `preprocessing_template.txt`: Template for pre-processing analysis with `Seurat`

- `aux_sh/` -> Bash scripts for each stage of the workflow
    - `autoflow_launcher.sh`: Generic [AutoFlow](https://github.com/seoanezonjic/autoflow) launcher for templates in `autoflow/`

    - `cellranger_mkfastq.sh`: Bash script for stage 0: converting BCL files to FASTQ files using Cell Ranger `mkfastq`

    - `compare_samples.sh`: Bash script for stage 2b: samples comparison

    - `general_report.sh`: Bash script for stage 3b: general preprocessing report

    - `qualimap.sh`: Bash script for stage 2c: [QualiMap](http://qualimap.conesalab.org/)

- `scripts/` -> Scripts for some stages of the workflow

    - `compare_samples.R`: R script for generating the QC report (stage 2b)

    - `general_report.R`: R script for generating the general preprocessing report (stage 3b)

    - `preprocessing.R`: R script for applying `Seurat` preprocessing and creating a report

- `templates/` -> R Markdown templates for reports

    - `clustering.Rmd`: Child template used in `preprocessing_report.Rmd` that performs clustering

    - `dimensionality_reduction.Rmd`: Child template used in `preprocessing_report.Rmd` that performs dimensionality reduction

    - `fastqc_report.Rmd`: Template for QC report (name will be changed soon to better represent the file function)

    - `feature_selection.Rmd`: Child template used in `preprocessing_report.Rmd` that performs feature selection

    - `preprocessing_report.Rmd`: Template for pre-processing report

    - `qc.Rmd`: Child template used in `preprocessing_report.Rmd` that performs cells quality control