#!/usr/bin/env Rscript

# Sergio Alías, 20230411
# Last modified 20230518

#################################
###   STAGE 2 PREPROCESSING   ###
###   preprocessing.R         ###
#################################

################  DISCLAIMER  ###################
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
              help="Scale factor for cell-level normalization"),
  make_option(c("-fs", "--hvgs"), type = "int",
              help="Number of HVG to be selected"),
  make_option(c("-nd", "--ndims"), type = "int",
              help="Number of PC to be plotted on the heatmap"),
  make_option(c("-ph", "--dimheatmapcells"), type = "int",
              help="Heatmap plots the 'extreme' cells on both ends of the spectrum")  
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
# Maybe plotting before and after (?)

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


#### Feature selection (HVG) ####

seu <- FindVariableFeatures(seu,
                            nfeatures = opt$hvgs)

fs_plot <- VariableFeaturePlot(pbmc.filtered)
fs_plot <- LabelPoints(plot = fs_plot, 
                     points = head(VariableFeatures(pbmc.filtered),
                                   10), 
                     repel = TRUE)

# TODO plot to the report
# TODO set an option for selecting a % of features, not an absolute number


#### PCA ####

seu <- RunPCA(seu,
              features = VariableFeatures(object = seu))

VizDimLoadings(seu, dims = 1:2, reduction = "pca")

DimPlot(seu, reduction = "pca")

DimHeatmap(seu,
           dims = 1:opt$ndims,
           cells = opt$dimheatmapcells,
           balanced = TRUE)

# TODO plots to the report

ElbowPlot(seu)

# TODO reconsirer adding JackStraw method


#### Clustering ####

seu <- FindNeighbors(seu, dims = 1:10)
seu <- FindClusters(seu, resolution = 0.5)


#### Non-linear dimensional reduction: UMAP, tSNE ####

# Maybe we need to install UMAP
# reticulate::py_install(packages = 'umap-learn')

seu <- RunUMAP(seu, dims = 1:10)

DimPlot(seu, reduction = "umap")

# TODO Add tSNE


#### Saving object as RDS ####

saveRDS(seu, paste0(name, ".seu.RDS"))