#! /usr/bin/env Rscript

# Sergio Alías, 20230606
# Last modified 20230613

#################################
###   STAGE 3 PREPROCESSING   ###
###   preprocessing.R         ###
#################################

#################
### Libraries ###
#################

library(optparse)
library(Seurat)


###################
### Custom libs ###
###################

root_path <- Sys.getenv("CODE_PATH") # daemon
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
  make_option(c("--mincells"), type = "integer",
              help="Min number of cells for which a feature was recorded"),
  make_option(c("--minfeats"), type = "integer",
              help="Min number of features for which a cell was recorded"),
  make_option(c("--maxfeats"), type = "integer",
              help="Max number of features for which a cell was recorded"),
  make_option(c("--maxcounts"), type = "integer",
              help="Max number of features for which a cell was recorded"),
  make_option(c("--percentmt"), type = "integer",
              help="Max percentage of reads mapped to mitochondrial genes for which a cell is recorded"),
  make_option(c("--normalmethod"), type = "character",
              help="Method for normalization. LogNormalize, CLR or RC"),
  make_option(c("--scalefactor"), type = "integer",
              help="Scale factor for cell-level normalization"),
  make_option(c("--hvgs"), type = "integer",
              help="Number of HVG to be selected"),
  make_option(c("--ndims"), type = "integer",
              help="Number of PC to be plotted on the heatmap"),
  make_option(c("--dimheatmapcells"), type = "integer",
              help="Heatmap plots the 'extreme' cells on both ends of the spectrum")  
)  

opt <- parse_args(OptionParser(option_list = option_list))


############
### Main ###
############

samples <- readLines(file.path(root_path, "samples_to_process.lst"))

report_folder <- Sys.getenv("report_folder") # config_daemon

aux_plots <- file.path(report_folder, "aux_plots")

experiment_name <- Sys.getenv("experiment_name") # config_daemon

results <- main_preprocessing_analysis(aux_plots = aux_plots,
                                       name = opt$name,
                                       experiment = experiment_name,
                                       input = opt$input,
                                       mincells = opt$mincells,
                                       minfeats = opt$minfeats,
                                       maxfeats = opt$maxfeats,
                                       maxcounts = opt$maxcounts,
                                       percentmt = opt$percentmt)

write_preprocessing_report(name = opt$name,
                           experiment = experiment_name,
                           template = file.path(root_path, 
                                                "templates",
                                                "preprocessing_report.Rmd"),
                           outdir = report_folder,
                           intermediate_files = "int_files")