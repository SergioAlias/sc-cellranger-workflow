#! /usr/bin/env Rscript

# Sergio Alías, 20231011
# Last modified 20231011


#################################
###   STAGE 3 PREPROCESSING   ###
###   subset.R                ###
#################################

#################
### Libraries ###
#################

library(optparse)


###################
### Custom libs ###
###################

root_path <- Sys.getenv("CODE_PATH") # daemon (TODO decide what to to with this)
source(file.path(root_path, "R", "preprocessing_library.R"))


############
### Args ###
############

option_list <- list(
  make_option(c("-e", "--exp_design"), type = "character",
              help="Input file with the experiment design"),
  make_option(c("-o", "--output"), type = "character",
              help="Output folder"),
  make_option(c("-c", "--column"), type = "character",
              help="Sample name")
)  


opt <- parse_args(OptionParser(option_list = option_list))



############
### Main ###
############

subsets <- do_subsetting(exp_design = opt$exp_design,
	                     column = opt$column)