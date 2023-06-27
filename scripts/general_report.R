#! /usr/bin/env Rscript

# Sergio Alías, 20230627
# Last modified 20230627


################################################
###   STAGE 3 GENERAL PREPROCESSING REPORT   ###
###   preprocessing.R                        ###
################################################

#################
### Libraries ###
#################

library(Seurat)
library(Matrix)


###################
### Custom libs ###
###################

root_path <- Sys.getenv("CODE_PATH") # daemon
source(file.path(root_path, "R", "preprocessing_library.R"))

############
### Main ###
############

seu_list <- list()
rds_files <- scan(Sys.getenv("SAMPLES_FILE"),
                  what = character())
main_folder <- Sys.getenv("PREPROC_RESULTS_FOLDER")
experiment <- Sys.getenv("experiment_name")

for (name in rds_files) {
  rds_object <- readRDS(file.path(main_folder,
                                  name,
                                  "preprocessing.R_0000",
                                  paste0(experiment,
                                         ".",
                                         name,
                                         ".seu.RDS")))
  
  seu_list <- append(seu_list, rds_object)
}

names(seu_list) <- rds_files

while (length(seu_list) > 1){
  seu_1 <- seu_list[1]
  seu_2 <- seu_list[2]
  name_1 <-names(seu_list)[1]
  name_2 <-names(seu_list)[2]
  prov_seu <- merge(seu_1, seu_2, add.cell.ids = c(name_1, name_2))
  seu_list <- seu_list[-c(1, 2)]
  seu_list <- c(prov_seu, seu_list)
}

seu <- seu_list[1]

write_preprocessing_report(name = "All samples",
                           experiment = experiment_name,
                           template = file.path(root_path, 
                                                "templates",
                                                "preprocessing_report.Rmd"),
                           outdir = report_folder,
                           intermediate_files = "int_files",
                           minqcfeats = Sys.getenv("preproc_qc_min_feats"),
                           percentmt = Sys.getenv("preproc_max_percent_mt"),
                           all_seu = seu)