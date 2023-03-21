#!/usr/bin/env bash
# ######################
# Sergio Alías, 20230315
# Last modified 20230321
# ######################
# Leave only one comment symbol on selected options
# Those with two commets will be ignored:
# The name to show in queue lists for this job:
##SBATCH -J sbatch-mkfast.sh

# Number of desired cpus (can be in any node):
#SBATCH --ntasks=1

# Number of desired cpus (all in same node):
##SBATCH --cpus-per-task=10

# Amount of RAM needed for this job:
#SBATCH --mem=60gb

# The time the job will be running:
#SBATCH --time=10:00:00

# To use GPUs you have to request them:
##SBATCH --gres=gpu:1

# If you need nodes with special features leave only one # in the desired SBATCH constraint line. cal is selected by default:
# * to request any machine without GPU - DEFAULT
#SBATCH --constraint=cal
# * to request only the machines with 128 cores and 1800GB of usable RAM
##SBATCH --constraint=bigmem
# * to request only the machines with 128 cores and 450GB of usable RAM (
##SBATCH --constraint=sr
# * to request only the machines with 52 cores and 187GB of usable RAM (
##SBATCH --constraint=sd

# Set output and error files
#SBATCH --error=job.%J.err
#SBATCH --output=job.%J.out

# Leave one comment in following line to make an array job. Then N jobs will be launched. In each one SLURM_ARRAY_TASK_ID will take one value from 1 to 100
##SBATCH --array=1-100

# To load some software (you can show the list with 'module avail'):
# module load software
# module load bcl2fastq/2.20
module load cellranger/7.0.0

# the program to execute with its parameters:
time ./test-count.sh
