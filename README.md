# sc-cellranger-workflow

## Workflow to perform single-cell RNA-seq analysis

*Last README update applies for v1.6.5*

This repo will contain a Cell Ranger workflow for going from BCL/FASTQ/10x files to fully analyzed single-cell data (incomplete). More info will be available soon.

---

### Usage

#### Before you start

![](https://github.com/SergioAlias/sc-cellranger-workflow/assets/96663968/88c5749e-ce56-4474-9f97-488843dafa18)

1. Edit `manage_input_files.sh` manually, from which all raw samples must be linked with a suitable name to another custom folder (e.g. `raw_data`). After editing the source folder, the output folder and the link commands, run the script and make sure everything is alright. This step must be performed for customizate samples according to the experiment and keep original sample names in order to facilitate sample backtracking (in case of failure) at same time. New sample names must follow [Illumina's naming convention](https://support.illumina.com/help/BaseSpace_OLH_009008/Content/Source/Informatics/BS/NamingConvention_FASTQ-files-swBS.htm) (`SampleName_S1_L001_R1_001.fastq.gz`). The first part of the name should be easy to interpret and it is recommended to use short and non redundant names.

2. Edit `samples_to_process.lst` adding sample names, one per line. **IMPORTANT:** Sample name must not include paired-end related information and file extension. If your sample name is `SampleName_S1_L001_R1_001.fastq.gz`, add a line with "SampleName" written in it. **MORE IMPORTANT:** Leave one empty line after all samples names. For instance, if you have 4 samples, your file sohuld have 5 lines, being the last one an empty line. Most text editors provides you with the line number so this should be easy to ensure.

3. Edit `config_daemon` and set the variables to the desired values. The config file already contains sensible values for the variables so you can do a first analysis with the default values and come back later to fine-tune.

---

#### Daemon execution

![imagen](https://github.com/SergioAlias/sc-cellranger-workflow/assets/96663968/6fd4492a-bd8c-48e2-93f7-996c5427ef13)

Next step is to execute a daemon stage. The execution is as simple as:

```
./daemon.sh N
```
Being N the stage number. You can use one of the following:

- **0**: BCL to FASTQ. The binary base call (BCL) sequence file format requires conversion to FASTQ format for use with user-developed or third-party data analysis tools. This stage asumes you have the full raw data folder with the BCL files

- **1**: FASTQ to counts. Aligns FASTQ reads to a reference transcriptome and generates feature-barcode matrices. This stage assumes you have the trasncriptome folder. You can get it from [here](https://support.10xgenomics.com/single-cell-gene-expression/software/downloads/latest/)

- **2a**: Launchs [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) on the FASTQ files
- **2b**: Generates a QC report with metrics from Cell Ranger count and FastQC
- **2c**: Launchs [QualiMap](http://qualimap.conesalab.org/) on the BAM file

- **3a**: Launchs single-cell analysis using the R package [Seurat](https://satijalab.org/seurat/) and generates one individual preprocessing report per sample
- **3b**: Generates a preprocessing report comparing all samples

---

### Developer section

Here you can find sub-sections useful for developers working on the workflow

#### List of files and directories

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
