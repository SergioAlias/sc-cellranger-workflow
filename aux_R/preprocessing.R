#!/usr/bin/env Rscript

# Sergio Alías, 20230411
# Last modified 20230414

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
              help="Threshold of the number of cells for which a feature was recorded"),
  make_option(c("-mf", "--minfeats"), type = "int",
              help="Threshold of the number of features for which a cell was recorded")
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

pbmc.filtered <- subset(pbmc.seurat, 
                        nFeature_RNA < 1250 &
                        nCount_RNA < 4000 & 
                        percent.mt < 5)

# TODO add variables to config_daemon and to preprocess.sh

######################################
######################################
######################################

# Caution!!!!!!!!!!!!! #
# The following code is not from Seurat and will be deleted after finding the best Seurat alternative to it #####
# Take care!!!!!!!!!!! #

## QC

sce <- scuttle::addPerCellQCMetrics(sce) # adds QC to colData
reasons <- scuttle::perCellQCFilters(colData(sce))
# colSums(as.matrix(reasons)) # for visualization
sce <- sce[, !reasons$discard]


## Normalization

set.seed(100)
clust <- scran::quickCluster(sce)
# table(clust) # for visualization
# deconv.sf <- scran::calculateSumFactors(sce, cluster = clust) # for visualization
# summary(deconv.sf) # for visualization
sce <- scran::computeSumFactors(sce, cluster = clust)
sce <- scuttle::logNormCounts(sce)


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