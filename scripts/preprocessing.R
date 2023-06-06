#! /usr/bin/env Rscript

# Sergio Alías, 20230606
# Last modified 20230606

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
              help="Project name"),
  make_option(c("-a", "--assay"), type = "character",
              help="Assay name"),
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

results <- append_message(samples)

outdir <- Sys.getenv("report_folder") # config_daemon
int_files <- file.path(outdir, "int_files")

if (!file.exists(int_files)){
  dir.create(int_files)
}

rmarkdown::render(file.path(root_path, 
                            "templates",
                            "preprocessing_report.Rmd"),
                  output_file = file.path(outdir,
                                          paste0(Sys.getenv("experiment_name"), # config_daemon
                                                 "_preprocessing_report.html")), 
                  clean = TRUE,
                  intermediates_dir = int_files)