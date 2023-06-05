#! /usr/bin/env bash

# Sergio Alías, 20230411
# Last modified 20230427

# STAGE 3 PREPROCESSING

#SBATCH --cpus-per-task=2
#SBATCH --mem=20gb
#SBATCH --constraint=cal
#SBATCH --error=job.preproc.%J.err
#SBATCH --output=job.preproc.%J.out


mkdir -p $preproc_outdir
. ~soft_bio_267/initializes/init_R

# Main

/usr/bin/time preprocessing.R --input $preproc_dir \
                              --output $preproc_outdir \
                              --name "CRR073022" \
                              --assay $preproc_assay \
                              --mincells $preproc_min_cells \
                              --minfeats $preproc_min_feats \
                              --maxfeats $preproc_max_feats \
                              --maxcounts $preproc_max_counts \
                              --percentmt $preproc_max_percent_mt \
                              --normalmethod $preproc_norm_method \
                              --scalefactor $preproc_scale_factor \
                              --hvgs $preproc_select_hvgs \
                              --ndims $preproc_pca_n_dims \
                              --dimheatmapcells $preproc_pca_n_cells
