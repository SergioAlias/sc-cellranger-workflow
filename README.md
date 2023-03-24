# sc-cellranger-workflow

This repo will contain a Cell Ranger workflow for going from FASTQ files to 10x input data (incomplete). More info will be available soon.

- `daemon.sh`: Daemon file for controlling the whole workflow (incomplete)
- `test-mkfastq.sh`: Bash script with some test on Cell Ranger ``mkfastq``
- `sbatch-mkfastq.sh`: SBATCH script for sending ``test-mkfastq.sh`` to the Slurm queue system
- `test-count.sh`: Bash script with some test on Cell Ranger `count`
- `sbatch-count.sh`: SBATCH script for sending `test-count.sh` to the Slurm queue system
- `AutoFlow_mkfastq`: Cell Ranger `mkfastq` wrapped in the [AutoFlow](https://github.com/seoanezonjic/autoflow) system
- `sbatch-autoflow.sh`: SBATCH script fo sending `AutoFlow_mkfastq` to the Slurm queue system
