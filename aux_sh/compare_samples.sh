#! /usr/bin/env bash

# Sergio Alias, 20230530
# Last modified 20230707

# STAGE 2 SAMPLES COMPARISON

#SBATCH -J compare_samples.sh
#SBATCH --cpus-per-task=3
#SBATCH --mem='5gb'
#SBATCH --constraint=cal
#SBATCH --time=0-01:00:00
#SBATCH --error=job.comp.%J.err
#SBATCH --output=job.comp.%J.out

# Setup

hostname

mkdir -p $report_folder

. ~soft_bio_267/initializes/init_ruby
. ~soft_bio_267/initializes/init_R
create_metric_table.rb $experiment_folder'/metrics' sample $experiment_folder'/metric_table'
all_report_files=$experiment_folder/metric_table
headers="t"

# Main

/usr/bin/time $LAB_SCRIPTS/create_report.R -t $REPORT_TEMPLATES_FOLDER/mapping_report.Rmd \
                                           -o $report_folder/mapping_report.html \
                                           -d $all_report_files \
                                           -H $headers