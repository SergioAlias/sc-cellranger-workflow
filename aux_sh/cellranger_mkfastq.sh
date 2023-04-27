#! /usr/bin/env bash

# Sergio Alías, 20230401
# Last modified 20230427

# STAGE 0 CONVERTING BCL FILES INTO FASTQ

#SBATCH --cpus-per-task=1
#SBATCH --mem=4gb
#SBATCH --constraint=cal
#SBATCH --error=job.bcl.%J.err
#SBATCH --output=job.bcl.%J.out


# Picasso modules

module load bcl2fastq/2.20
module load cellranger/7.0.0


# Main



if [ "$bcl_use_sample_sheet" == "TRUE" ] ; then
    /usr/bin/time cellranger mkfastq --id=$experiment_name \
                                     --run=$bcl_folder \
                                     --samplesheet=$bcl_sample_sheet
else
    /usr/bin/time cellranger mkfastq --id=$experiment_name \
                                     --run=$bcl_folder \
                                     --csv=$bcl_simple_csv
fi