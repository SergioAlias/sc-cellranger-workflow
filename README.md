# sc-cellranger-workflow

## Workflow to perform single-cell RNA-seq analysis

This repo will contain a Cell Ranger workflow for going from BCL/FASTQ/10x files to fully analyzed single-cell data (incomplete). More info will be available soon.

- `daemon.sh`: Daemon file for controlling the whole workflow (incomplete)

- `config_daemon`: File with default and global variables (incomplete)

- `manage_input_files.sh`: For creating soft links to the input data. It must be edited manually

- `aux_sh/` -> Bash scripts for each stage of the workflow
    - `cellranger_mkfastq.sh`: Bash script for stage 0: converting BCL files to FASTQ files using Cell Ranger `mkfastq`
    - `cellranger_count.sh`: Bash script for stage 1: converting FASTQ files to a proper format for single-cell analysis using Cell Ranger `count`
     - `preprocess.sh`: Bash script for stage 2: preprocessing using R package `Seurat` (incomplete)

- `aux_R/` -> R scripts for some stages of the workflow
    - `preprocessing.R`: R script implementing `Seurat` preprocessing (incomplete)

---

- `deprecated/` -> Deprecated scripts (will be removed soon)
    - `test-mkfastq.sh`: Bash script with some test on Cell Ranger ``mkfastq``
    - `sbatch-mkfastq.sh`: SBATCH script for sending ``test-mkfastq.sh`` to the Slurm queue system
    - `test-count.sh`: Bash script with some test on Cell Ranger `count`
    - `sbatch-count.sh`: SBATCH script for sending `test-count.sh` to the Slurm queue system
    - `AutoFlow_mkfastq`: Cell Ranger `mkfastq` wrapped in the [AutoFlow](https://github.com/seoanezonjic/autoflow) system
    - `sbatch-autoflow.sh`: SBATCH script fo sending `AutoFlow_mkfastq` to the Slurm queue system
