#! /usr/bin/env bash

# Sergio Alías, 20230401
# Last modified 20230401

# STAGE 0 CONVERTING BCL FILES INTO FASTQ

#SBATCH --ntasks=1
#SBATCH --mem=4gb
#SBATCH --error=bcl.%J.err
#SBATCH --output=bcl.%J.out


# Picasso modules

module load bcl2fastq/2.20
module load cellranger/7.0.0

# Vars

EXP_ID=$1

# Main

cellranger mkfastq --id=$EXP_ID \
                   --run=cellranger-tiny-bcl-1.2.0 \
                   --csv=cellranger-tiny-bcl-simple-1.2.0.csv #\
                   # --rc-i2-override=true
                   # --samplesheet=cellranger-tiny-bcl-samplesheet-1.2.0.csv