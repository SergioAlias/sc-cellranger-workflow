#!/usr/bin/env Rscript

# Sergio Alías, 20230411
# Last modified 20230417

#################################
###   STAGE 2 PREPROCESSING   ###
###   preprocessing.R         ###
#################################

################  DISCLAIMER  ###################
#                                               #
#          This script is incomplete            #
#                                               #
# This script has not been tried in Picasso yet #
#                                               #
#################################################

#################
### Libraries ###
#################

library(optparse)
library(Seurat)


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
  make_option(c("-mc", "--mincells"), type = "int",
              help="Min number of cells for which a feature was recorded"),
  make_option(c("-mf", "--minfeats"), type = "int",
              help="Min number of features for which a cell was recorded"),
  make_option(c("-xf", "--maxfeats"), type = "int",
              help="Max number of features for which a cell was recorded"),
  make_option(c("-xc", "--maxcounts"), type = "int",
              help="Max number of features for which a cell was recorded"),
  make_option(c("-mt", "--percentmt"), type = "int",
              help="Max percentage of reads mapped to mitochondrial genes for which a cell is recorded"),
  make_option(c("-nm", "--normalmethod"), type = "character",
              help="Method for normalization. LogNormalize, CLR or RC"),
  make_option(c("-sf", "--scalefactor"), type = "int",
              help="Scale factor for cell-level normalization")
)  

opt <- parse_args(OptionParser(option_list = option_list))


############
### Main ###
############


#### Read input data ####

mtx <- Read10X(opt$input)

seu <- CreateSeuratObject(counts = mtx,
                          project = opt$name, 
                          assay = opt$assay, 
                          min.cells = opt$mincells, 
                          min.features = opt$minfeats
                         )


#### QC ####

##### Reads mapped to mitochondrial genes #####

# In human ENSEMBL gene symbol annotations, mitochondrial genes are annotated starting with a MT- string

seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = "^MT-")

# We can visualize metrics:

VlnPlot(seu, features = c('nFeature_RNA', 'nCount_RNA', 'percent.mt'))

# TODO think about best way of putting all plots in a future HTML report

##### Filtering out cells #####

seu <- subset(seu, 
              nFeature_RNA < opt$maxfeats &
              nCount_RNA < opt$maxcounts & 
              percent.mt < opt$percentmt)


#### Normalization ####

seu <- NormalizeData(seu,
                     normalization.method = opt$normalmethod,
                     scale.factor = opt$scalefactor)

seu <- ScaleData(seu)

######################################
######################################
######################################




# Caution!!!!!!!!!!!!! #
# The following code is not from Seurat and will be deleted after finding the best Seurat alternative to it #####
# Take care!!!!!!!!!!! #


## Feature selection

dec <- scran::modelGeneVar(sce)
sce.original <- sce
chosen <- scran::getTopHVGs(dec, prop=0.1)
sce <- sce[chosen,]
altExp(sce, "original") <- sce.original
rm(sce.original)


## Dimensionality reduction (PCA)

set.seed(100)
sce <- scran::fixedPCA(sce, subset.row = NULL)


## Clustering

nn.clusters <- scran::clusterCells(sce, use.dimred = "PCA")
colLabels(sce) <- nn.clusters


## Saving sce as RDS

saveRDS(sce, paste0(name, ".sce.RDS"))