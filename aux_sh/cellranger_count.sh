#! /usr/bin/env bash

# Sergio Alías, 20230406
# Last modified 20230427

# STAGE 1 OBTAINING COUNTS FROM FASTQ FILES

#SBATCH --cpus-per-task=10
#SBATCH --mem=30gb
#SBATCH --constraint=cal
#SBATCH --error=job.counts.%J.err
#SBATCH --output=job.counts.%J.out


# Picasso modules

module load cellranger/7.0.0

# Main

/usr/bin/time cellranger count --id=$experiment_name \
                               --sample=$experiment_name \
                               --transcriptome=$count_transcriptome \
                               --fastqs=$read_path \
                               --localcores=10 \
                               --localmem=30


# time count_table.R --input $experiment_name