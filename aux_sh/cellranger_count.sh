#! /usr/bin/env bash

# Sergio Alías, 20230406
# Last modified 20230406

# STAGE 1 OBTAINING COUNTS FROM FASTQ FILES

#SBATCH --cpus-per-task=10
#SBATCH --mem=60gb
#SBATCH --constraint=cal
#SBATCH --error=job.counts.%J.err
#SBATCH --output=job.counts.%J.out


# Picasso modules

module load cellranger/7.0.0

# Vars

EXP_ID=$1

# Main

time cellranger count --id=$EXP_ID \
                      --transcriptome=refdata-gex-GRCh38-2020-A \
                      --fastqs=raw_data/ \
                      --sample=pbmc_1k_v3 \
                      --localcores=10 \
                      --localmem=60