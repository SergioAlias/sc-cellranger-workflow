#! /usr/bin/env bash

# Sergio Alías, 20230411
# Last modified 20230411

# STAGE 2 PREPROCESSING

#SBATCH --cpus-per-task=10
#SBATCH --mem=30gb
#SBATCH --constraint=cal
#SBATCH --error=job.counts.%J.err
#SBATCH --output=job.counts.%J.out

# Picasso modules

