#! /usr/bin/env Rscript

# Sergio Alías, 20230707
# Last modified 20230721


######################################
###   STAGE 2 SAMPLES COMPARISON   ###
###   compare_samples.R            ###
######################################

#################
### Libraries ###
#################

library(optparse)

###################
### Custom libs ###
###################

root_path <- Sys.getenv("CODE_PATH") # daemon
source(file.path(root_path, "R", "qc_library.R"))


############
### Args ###
############

option_list <- list(
  make_option(c("-m", "--metrics"), type = "character",
              help="Metrics file in wide format"),
  make_option(c("-l", "--long_metrics"), type = "character",
              help="Metrics file in long format"),
  make_option(c("--cellranger_metrics"), type = "character",
              help="Cell Ranger metrics file in wide format"),
  make_option(c("--cellranger_long_metrics"), type = "character",
              help="Cell Ranger Metrics file in long format"),
  make_option(c("-o", "--output"), type = "character",
              help="Output folder"),
  make_option(c("-e", "--experiment_name"), type = "character",
              help="Experiment name")
)  

opt <- parse_args(OptionParser(option_list = option_list))


############
### Main ###
############

write_qc_report(name = "All samples",
                experiment = opt$experiment_name,
                template = file.path(root_path,
                                     "templates",
                                     "fastqc_report.Rmd"),
                outdir = opt$output,
                intermediate_files = "int_files",
                metrics = opt$metrics,
                long_metrics = opt$long_metrics,
                cellranger_metrics = opt$cellranger_metrics,
                cellranger_long_metrics = opt$cellranger_long_metrics)