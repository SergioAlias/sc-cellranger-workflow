#! /usr/bin/env bash

# Sergio Alías, 20230627
# Last modified 20230627

# STAGE 3 PREPROCESSING

#SBATCH --cpus-per-task=2
#SBATCH --mem=30gb
#SBATCH --constraint=cal
#SBATCH --error=job.genrep.%J.err
#SBATCH --output=job.genrep.%J.out


mkdir -p $PREPROC_RESULTS_FOLDER
. ~soft_bio_267/initializes/init_R

# Main

/usr/bin/time general_report.R
