#! /usr/bin/env Rscript

# Sergio Alías, 20230627
# Last modified 20230929


################################################
###   STAGE 3 GENERAL PREPROCESSING REPORT   ###
###   general_report.R                       ###
################################################

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
              help="Input file with samples names"),
  make_option(c("-o", "--output"), type = "character",
              help="Folder where the report is written"),
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
  make_option(c("--experiment_name"), type = "character",
              help="Experiment name"),
  make_option(c("--results_folder"), type = "character",
              help="Folder with the preprocessing results"),
  make_option(c("--resolution"), type = "double",
              help="Granularity of the clustering"),
  make_option(c("--integrative_analysis"), type = "character",
              help="TRUE if we want to integrate samples, FALSE otherwise")
)  

opt <- parse_args(OptionParser(option_list = option_list))


############
### Main ###
############

seu_list <- list()
raw_seu_list <- list()
marker_gene_list <- list()
rds_files <- scan(opt$input,
                  what = character())
main_folder <- opt$results_folder
experiment_name <- opt$experiment_name

for (name in rds_files) {
  if (opt$integrative_analysis == "FALSE") { # We don't need preprocessed objects if we do integrative analysis
      rds_object <- readRDS(file.path(main_folder,
                                      name,
                                      "preprocessing.R_0000",
                                      paste0(experiment_name,
                                             ".",
                                             name,
                                             ".seu.RDS")))
      seu_list <- append(seu_list, rds_object)
  }
  
  raw_object <- readRDS(file.path(main_folder,
                                  name,
                                  "preprocessing.R_0000",
                                  paste0(experiment_name,
                                         ".",
                                         name,
                                         ".before.seu.RDS")))
  raw_seu_list <- append(raw_seu_list, raw_object)

  marker_genes <- readRDS(file.path(main_folder,
                                    name,
                                    "preprocessing.R_0000",
                                    paste0(experiment_name,
                                           ".",
                                           name,
                                           ".markers.RDS")))
  marker_gene_list <- append(marker_gene_list, marker_genes)
}

if (opt$integrative_analysis == "FALSE") { # Same reason as line 82
  seu <- seu_list
  before.seu <- raw_seu_list
  report_name <- "All samples"
} else { # "If we want an integrative analysis"
  report_name <- "Integrated samples"
  main_preprocessing_analysis(name = report_name,
                              experiment = experiment_name,
                              input = opt$input,
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
                              integrate = opt$integrative_analysis)
  seu <- readRDS(file.path(Sys.getenv(PREPROC_RESULTS_FOLDER), paste0(experiment_name, ".", report_name, ".seu.RDS")))
}


write_preprocessing_report(name = report_name,
                           experiment = experiment_name,
                           template = file.path(root_path, 
                                                "templates",
                                                "preprocessing_report.Rmd"),
                           outdir = opt$output,
                           intermediate_files = "int_files",
                           minqcfeats = opt$minqcfeats,
                           percentmt = opt$percentmt,
                           hvgs = opt$hvgs,
                           resolution = opt$resolution,
                           all_seu = list(seu, before.seu, marker_gene_list),
                           integrate = opt$integrative_analysis)