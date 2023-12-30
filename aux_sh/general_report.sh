#! /usr/bin/env bash

# Sergio Alías, 20230627
# Last modified 20231204

# STAGE 3 PREPROCESSING

#SBATCH --cpus-per-task=2
#SBATCH --mem=500gb
#SBATCH --constraint=cal
#SBATCH --error=job.genrep.%J.err
#SBATCH --output=job.genrep.%J.out


mkdir -p $PREPROC_RESULTS_FOLDER
. ~soft_bio_267/initializes/init_R

if [ "$integrative_analysis" == "TRUE" ] ; then
    export SAMPLES_FILE=$integration_file
fi

# Main

/usr/bin/time general_report.R --input $SAMPLES_FILE \
                               --output $report_folder \
                               --filter $preproc_filter \
                               --mincells $preproc_init_min_cells \
                               --minfeats $preproc_init_min_feats \
                               --minqcfeats $preproc_qc_min_feats \
                               --percentmt $preproc_max_percent_mt \
                               --normalmethod $preproc_norm_method \
                               --scalefactor $preproc_scale_factor \
                               --hvgs $preproc_select_hvgs \
                               --ndims $preproc_pca_n_dims \
                               --dimheatmapcells $preproc_pca_n_cells \
                               --experiment_name $experiment_name \
                               --results_folder $PREPROC_RESULTS_FOLDER \
                               --resolution $preproc_resolution \
                               --integrative_analysis $integrative_analysis \
                               --int_sec_cond $int_sec_cond