#! /usr/bin/env Rscript

# Sergio Alías, 20231011
# Last modified 20231226


#################################
###   STAGE 3 PREPROCESSING   ###
###   prior_integration.R     ###
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
  make_option(c("-d", "--exp_design"), type = "character",
              help="Input file with the experiment design"),
  make_option(c("-o", "--output"), type = "character",
              help="Output folder"),
  make_option(c("-c", "--condition"), type = "character",
              help="Name of the column in the experimental design file we want to use for integration"),
  make_option(c("-i", "--integration_file"), type = "character",
              help="Name of file that will contain integration subset names"),
  make_option(c("-e", "--experiment_name"), type = "character",
              help="Experiment name"),
  make_option(c("--count_folder"), type = "character",
              help="Count results folder")
)  


opt <- parse_args(OptionParser(option_list = option_list))



############
### Main ###
############

# Divide sample names by the condition

exp_subsets <- do_subsetting(exp_design = opt$exp_design,
	                           column = opt$condition)


# Add them to a file so we pass it to the SAMPLES_FILE variable in autoflow_launcher.sh

if (file.exists(opt$integration_file)) { # Integration subsets file creation (must be empty if already exists from previous runs)
  file_conn <- file(opt$integration_file, open = "w")
  close(file_conn)
} else {
  file.create(opt$integration_file)
}

cat(names(exp_subsets), file = opt$integration_file, sep = "\n", append = FALSE) # the second quotes are for having an extra empty line at the end (see README in Github for details)


# Create and save the "before" (with no preprocessing) Seurat object

for (cond in names(exp_subsets)){
  folder_name <- file.path(opt$output, cond)
  if (!file.exists(folder_name)) {
    dir.create(folder_name)
  }
  seu <- merge_condition(exp_cond = cond,
                         samples = exp_subsets[[cond]],
                         exp_design = opt$exp_design,
                         count_path = opt$count_folder)
  saveRDS(seu, file = file.path(folder_name, paste0(opt$experiment_name, ".", cond, ".before.seu.RDS")))
}