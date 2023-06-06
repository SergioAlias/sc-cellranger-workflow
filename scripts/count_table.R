#!/usr/bin/env Rscript

# Sergio Alías, 20230421
# Last modified 20230421

#####################################################
###   STAGE 1 OBTAINING COUNTS FROM FASTQ FILES   ###
###   count_table.R                               ###
#####################################################

#################
### Libraries ###
#################

library(optparse)
library(data.table)

############
### Args ###
############

option_list <- list(
  make_option(c("-i", "--input"), type = "character",
              help="Input folder with Cellranger count results")
)
  
opt <- parse_args(OptionParser(option_list = option_list))


############
### Main ###
############


#### Read input data ####

metrics <- fread(file.path(opt$input,
                           "outs",
                           "metrics_summary.csv"))

count_table <- data.table(Sample = opt$input,
                          Number_of_cells = metrics$`Estimated Number of Cells`,
                          Number_of_reads = metrics$`Number of Reads`,
                          Mean_reads_per_cell = metrics$`Mean Reads per Cell`,
                          Median_genes_per_cell = metrics$`Median Genes per Cell`
                          )