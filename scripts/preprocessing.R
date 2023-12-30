#! /usr/bin/env Rscript

# Sergio Alías, 20230606
# Last modified 20231226


#################################
###   STAGE 3 PREPROCESSING   ###
###   preprocessing.R         ###
#################################

#################
### Libraries ###
#################

library(optparse)
library(Seurat)
library(scCustomize)

###################
### Custom libs ###
###################

root_path <- Sys.getenv("CODE_PATH") # daemon (TODO decide what to to with this)
source(file.path(root_path, "R", "preprocessing_library.R"))

############
### Args ###
############

option_list <- list(
  make_option(c("-i", "--input"), type = "character",
              help="Input folder with 10X data"),
  make_option(c("-o", "--output"), type = "character",
              help="Output folder"),
  make_option(c("-n", "--name"), type = "character",
              help="Sample name"),
  make_option(c("--filter"), type = "character",
              help="TRUE for using only detected cell-associated barcodes, FALSE for using all detected barcodes"),
  make_option(c("--mincells"), type = "integer",
              help="Min number of cells for which a feature was recorded"),
  make_option(c("--minfeats"), type = "integer",
              help="Min number of features for which a cell was recorded"),
  make_option(c("--minqcfeats"), type = "integer",
              help="Min number of features for which a cell was selected in QC"),
  make_option(c("--percentmt"), type = "integer",
              help="Max percentage of reads mapped to mitochondrial genes for which a cell is recorded"),
  make_option(c("--normalmethod"), type = "character",
              help="Method for normalization. LogNormalize, CLR or RC"),
  make_option(c("--scalefactor"), type = "integer",
              help="Scale factor for cell-level normalization"),
  make_option(c("--hvgs"), type = "integer",
              help="Number of HVG to be selected"),
  make_option(c("--ndims"), type = "integer",
              help="Number of PC to be used for clustering / UMAP / tSNE"),
  make_option(c("--dimheatmapcells"), type = "integer",
              help="Heatmap plots the 'extreme' cells on both ends of the spectrum"),
  make_option(c("--report_folder"), type = "character",
              help="Folder where the report is written"),
  make_option(c("--experiment_name"), type = "character",
              help="Experiment name"),
  make_option(c("--resolution"), type = "double",
              help="Granularity of the clustering"),
  make_option(c("--integrate"), type = "character",
              help="TRUE if integrative analysis, FALSE otherwise")
)  


opt <- parse_args(OptionParser(option_list = option_list))


############
### Main ###
############

main_preprocessing_analysis(name = opt$name,
                            experiment = opt$experiment_name,
                            input = opt$input,
                            output = opt$output,
                            filter = opt$filter,
                            mincells = opt$mincells,
                            minfeats = opt$minfeats,
                            minqcfeats = opt$minqcfeats,
                            percentmt = opt$percentmt,
                            normalmethod = opt$normalmethod,
                            scalefactor = opt$scalefactor,
                            hvgs = opt$hvgs,
                            ndims = opt$ndims,
                            resolution = opt$resolution,
                            integrate = as.logical(opt$integrate))

write_preprocessing_report(name = opt$name,
                           experiment = opt$experiment_name,
                           template = file.path(root_path, 
                                                "templates",
                                                "preprocessing_report.Rmd"),
                           outdir = opt$report_folder,
                           intermediate_files = "int_files",
                           minqcfeats = opt$minqcfeats,
                           percentmt = opt$percentmt,
                           hvgs = opt$hvgs,
                           resolution = opt$resolution)