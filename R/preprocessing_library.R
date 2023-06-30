# Sergio Alías, 20230606
# Last modified 20230630

##########################################################################
########################## PRE-PROCESSING LIBRARY ########################
##########################################################################


#' append_message
#' Append experiments names with a Hello World message. Useful for testing :D
#'
#' @param name: expermient name
#' @param message: message to append. Default: "Hello World"
#' 
#' @keywords test
#' 
#' @return character vector
append_message <- function(name, message = "Hello World") {paste(name, message)}


##########################################################################


#' read_input
#' Read Cellranger input into Seurat
#'
#' @param name: sample name
#' @param input: directory with the single-cell data
#' @param mincells: min number of cells for which a feature is recorded
#' @param minfeats: min number of features for which a cell is recorded
#' 
#' @keywords preprocessing, input
#' 
#' @return Seurat object
read_input <- function(name, input, mincells, minfeats){
  mtx <- Read10X(input)
  seu <- CreateSeuratObject(counts = mtx,
                            project = name, 
                            assay = "scRNAseq", 
                            min.cells = mincells, 
                            min.features = minfeats
                            )
  mtx <- NULL
  return(seu)
}


##########################################################################


#' do_qc
#' Perform Quality Control
#'
#' @param name: sample name
#' @param expermient: experiment name
#' @param seu: Seurat object
#' @param minqcfeats: min number of features for which a cell is selected
#' @param percentmt: max percentage of reads mapped to mitochondrial genes for which a cell is selected
#' 
#' @keywords preprocessing, qc
#' 
#' @return Seurat object
do_qc <- function(name, experiment, seu, minqcfeats, percentmt){
  #### QC ####
  
  ##### Reads mapped to mitochondrial genes #####
  
  # In human ENSEMBL gene symbol annotations, mitochondrial genes are annotated starting with a MT- string
  
  seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = "^MT-")
  
  # We can define ribosomal proteins (their names start with RPS or RPL)
  
  seu[["percent.rb"]] <- PercentageFeatureSet(seu, pattern = "^RP[SL]")
  
  # Save before version
  
  saveRDS(seu, paste0(experiment, ".", name, ".before.seu.RDS"))

  ##### Filtering out cells #####

  # seu[['QC']] <- ifelse(seu@meta.data$Is_doublet == 'True','Doublet','Pass')
  seu[['QC']] <- ifelse(TRUE,'Pass','This should not happen') # provisional until I code the dublet detection stuff (see previous line)
  seu[['QC']] <- ifelse(seu@meta.data$nFeature_scRNAseq < minqcfeats & seu@meta.data$QC == 'Pass','Low_nFeature',seu@meta.data$QC)
  seu[['QC']] <- ifelse(seu@meta.data$nFeature_scRNAseq < minqcfeats & seu@meta.data$QC != 'Pass' & seu@meta.data$QC != 'Low_nFeature',paste('Low_nFeature',seu@meta.data$QC,sep = ','),seu@meta.data$QC)
  seu[['QC']] <- ifelse(seu@meta.data$percent.mt > percentmt & seu@meta.data$QC == 'Pass','High_MT',seu@meta.data$QC)
  seu[['QC']] <- ifelse(seu@meta.data$nFeature_scRNAseq < minqcfeats & seu@meta.data$QC != 'Pass' & seu@meta.data$QC != 'High_MT',paste('High_MT',seu@meta.data$QC,sep = ','),seu@meta.data$QC)
  table(seu[['QC']])
  
  seu <- subset(seu, subset = QC == 'Pass')
  
  return(seu)
}


##########################################################################


#' do_dimred
#' Perform linear (PCA) and non-linear (UMAP/tSNE) dimensionality reduction
#'
#' @param seu: Seurat object
#' @param ndims: Number of PC to be used for UMAP / tSNE
#' 
#' @keywords preprocessing, dimensionality, reduction, PCA, UMAP, tSNE
#' 
#' @return Seurat object
do_dimred <- function(seu, ndims){
  seu <- RunPCA(seu, features = VariableFeatures(object = seu))
  seu <- RunUMAP(seu, dims = 1:ndims)
  seu <- RunTSNE(seu, dims = 1:ndims)
  return(seu)
}


##########################################################################


#' do_clustering
#' Perform clustering of cells
#'
#' @param seu: Seurat object
#' @param ndims: Number of PC to be used for clustering
#' @param: resolution: Granularity of the downstream clustering (higher values -> greater number of clusters)
#' 
#' @keywords preprocessing, clustering
#' 
#' @return Seurat object
do_clustering <- function(seu, ndims, resolution){
  seu <- FindNeighbors(seu, dims = 1:ndims)
  seu <- FindClusters(seu, resolution = resolution)
  return(seu)
}


##########################################################################



#' main_preprocessing_analysis
#' Main preprocessing function that performs all the analyses
#'
#' @param name: sample name
#' @param expermient: experiment name
#' @param input: directory with the single-cell data
#' @param filter: TRUE for using only detected cell-associated barcodes, FALSE for using all detected barcodes
#' @param mincells: min number of cells for which a feature is recorded
#' @param minfeats: min number of features for which a cell is recorded
#' @param minqcfeats: min number of features for which a cell is selected
#' @param percentmt: max percentage of reads mapped to mitochondrial genes for which a cell is selected
#' @param normalmethod: Normalization method
#' @param hvgs: Number of HVGs to be selected
#' @param ndims: Number of PC to be used for clustering / UMAP / tSNE
#' @param: resolution: Granularity of the downstream clustering (higher values -> greater number of clusters)
#' 
#' @keywords preprocessing, main
#' 
#' @return Seurat object
main_preprocessing_analysis <- function(name, experiment, input, filter, mincells, minfeats, minqcfeats,
                                        percentmt, normalmethod, scalefactor, hvgs, ndims, resolution){

  # Input selection
  
  if (filter == "TRUE"){
    input <- file.path(input, "filtered_feature_bc_matrix")
  } else {
    input <- file.path(input, "raw_feature_bc_matrix")
  }
  
  # Input reading
  
  seu <- read_input(name = name, 
                    input = input,
                    mincells = mincells,
                    minfeats = minfeats)
  
  # QC
  
  seu <- do_qc(name = name,
               experiment = experiment,
               seu = seu,
               minqcfeats = minqcfeats, 
               percentmt = percentmt)
  
  # Normalization
  
  seu <- NormalizeData(seu, normalization.method = normalmethod, scale.factor = scalefactor)
  
  # Feature selection
  
  seu <- FindVariableFeatures(seu, nfeatures = hvgs)
  
  # Data scaling
  
  seu <- ScaleData(seu, features = rownames(seu))
  
  # Dimensionality reduction (PCA/UMAP/tSNE)
  
  seu <- do_dimred(seu = seu,
                   ndims = ndims)
  
  # Clustering
  
  seu <- do_clustering(seu = seu,
                       ndims = ndims,
                       resolution = resolution)
  
  # Save final Seurat object
  
  saveRDS(seu, paste0(experiment, ".", name, ".seu.RDS"))
}


##########################################################################


#' write_preprocessing_report
#' Write preprocessing HTML report
#' 
#' @param name: sample name
#' @param expermient: experiment name
#' @param template: Rmd template
#' @param outdir: output directory
#' @param intermediate_files: directory for saving intermediate files in case pandoc fails
#' @param minqcfeats: min number of features for which a cell is selected
#' @param percentmt: max percentage of reads mapped to mitochondrial genes for which a cell is selected
#' @param hvgs: Number of HVGs to be selected
#' @param: resolution: Granularity of the downstream clustering (higher values -> greater number of clusters)
#' @param all_seu: NULL if creating an indiviual report (daemon 3a). A list of 2 Seurat objects if creating a general report (daemon 3b)
#' 
#' @keywords preprocessing, write, report
#' 
#' @return nothing
write_preprocessing_report <- function(name, experiment, template, outdir, intermediate_files, minqcfeats, percentmt, hvgs, resolution, all_seu = NULL){
  int_files <- file.path(outdir, intermediate_files)
  if (!file.exists(int_files)){
    dir.create(int_files)
  }
  if (is.null(all_seu)){
    main_folder <- Sys.getenv("PREPROC_RESULTS_FOLDER")
    
    seu <- readRDS(file.path(main_folder,
                             name,
                             "preprocessing.R_0000",
                             paste0(experiment,
                                    ".",
                                    name,
                                    ".seu.RDS")))
    before.seu <- readRDS(file.path(main_folder,
                                    name,
                                    "preprocessing.R_0000",
                                    paste0(experiment,
                                           ".",
                                           name,
                                           ".before.seu.RDS")))
    
  } else {
    seu <- all_seu[[1]]
    before.seu <- all_seu[[2]]
  }  
  rmarkdown::render(template,
                    output_file = file.path(outdir,
                                            paste0(experiment,
                                                   "_",
                                                   name,
                                                   "_preprocessing_report.html")), 
                    clean = TRUE,
                    intermediates_dir = int_files)
}