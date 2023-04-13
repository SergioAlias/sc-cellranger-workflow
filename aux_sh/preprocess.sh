#! /usr/bin/env bash

# Sergio Alías, 20230411
# Last modified 20230414

# STAGE 2 PREPROCESSING

#SBATCH --cpus-per-task=2
#SBATCH --mem=20gb
#SBATCH --constraint=cal
#SBATCH --error=job.preproc.%J.err
#SBATCH --output=job.preproc.%J.out

# Main

time preprocessing.R --input $preproc_dir \
                     --output $preproc_outdir \
                     --name $experiment_name \
                     --assay $preproc_assay \
                     --mincells $preproc_min_cells \
                     --minfeats $preproc_min_feats
