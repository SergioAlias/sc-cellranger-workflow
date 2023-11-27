# Sergio Alías, 20230606
# Last modified 20231127

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
#' @param minqcfeats: Min number of features for which a cell is selected
#' @param percentmt: Max percentage of reads mapped to mitochondrial genes for which a cell is selected
#' @param save_before_seu: Whether to save Seurat object before QC
#'
#' @keywords preprocessing, qc
#' 
#' @return Seurat object
do_qc <- function(name, experiment, seu, minqcfeats, percentmt, save_before_seu = TRUE){
  #### QC ####
  
  ##### Reads mapped to mitochondrial genes #####
  
  # In human ENSEMBL gene symbol annotations, mitochondrial genes are annotated starting with a MT- string
  
  seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = "^MT-")
  
  # We can define ribosomal proteins (their names start with RPS or RPL)
  
  seu[["percent.rb"]] <- PercentageFeatureSet(seu, pattern = "^RP[SL]")
  
  ##### Filtering out cells #####

  # seu[['QC']] <- ifelse(seu@meta.data$Is_doublet == 'True','Doublet','Pass')
  seu[['QC']] <- ifelse(TRUE,'Pass','This should not happen') # provisional until I code the dublet detection stuff (see previous line)
  seu[['QC']] <- ifelse(seu@meta.data$nFeature_scRNAseq < minqcfeats & seu@meta.data$QC == 'Pass','Low_nFeature',seu@meta.data$QC)
  seu[['QC']] <- ifelse(seu@meta.data$nFeature_scRNAseq < minqcfeats & seu@meta.data$QC != 'Pass' & seu@meta.data$QC != 'Low_nFeature',paste('Low_nFeature',seu@meta.data$QC,sep = ','),seu@meta.data$QC)
  seu[['QC']] <- ifelse(seu@meta.data$percent.mt > percentmt & seu@meta.data$QC == 'Pass','High_MT',seu@meta.data$QC)
  seu[['QC']] <- ifelse(seu@meta.data$nFeature_scRNAseq < minqcfeats & seu@meta.data$QC != 'Pass' & seu@meta.data$QC != 'High_MT',paste('High_MT',seu@meta.data$QC,sep = ','),seu@meta.data$QC)
  table(seu[['QC']])
  
  # Save before version
  
  if (save_before_seu){
    saveRDS(seu, paste0(experiment, ".", name, ".before.seu.RDS"))
  }
  
  seu <- subset(seu, subset = QC != 'High_MT,Low_nFeature')
  
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
#' @param resolution: Granularity of the downstream clustering (higher values -> greater number of clusters)
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


#' do_marker_gene_selection
#' Perform marker gene selection
#' TODO this function is harcoded - make proper variables
#'
#' @param seu: Seurat object
#' @param name: sample name
#' @param expermient: experiment name
#' 
#' @keywords preprocessing, marker, gene
#' 
#' @return Nothing
do_marker_gene_selection <- function(seu, name, experiment){
  markers <- FindAllMarkers(seu, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
  saveRDS(markers, paste0(experiment, ".", name, ".markers.RDS"))
}


##########################################################################


#' do_subsetting
#' Subset samples according to experimental condition
#' TODO this function is harcoded - make proper variables
#'
#' @param exp_design: Experiment design table in CSV format
#' @param column: Column with the condition used for subsetting
#' 
#' @keywords preprocessing, subsetting, integration
#' 
#' @return List of vectors with sample names
do_subsetting <- function(exp_design, column){
  exp_design <- read.csv(exp_design)
  subsets <- split(exp_design$code, exp_design[[column]])
  return(subsets)
}


##########################################################################


##########################################################################


#' merge_condition
#' Merge samples sharing an experimental condition
#'
#' @param exp_cond: Experimental condition
#' @param samples: Vector of samples with that experimental condition
#' 
#' @keywords preprocessing, merging, integration
#' 
#' @return Merged Seurat object
do_subsetting <- function(exp_cond, samples){ # (TODO Finsih this and add secondary condition)
  seu.list.48 <- sapply(ids.48, function(i){ # Loading
    d10x <- Read10X(file.path(dataset_loc,i,"cellranger_0000",i,"outs/filtered_feature_bc_matrix"))
    colnames(d10x) <- paste(sapply(strsplit(colnames(d10x),split="-"),'[[',1L),i,sep="-")
    seu <- CreateSeuratObject(counts = d10x, project = "neurorg", min.cells = 1, min.features = 1)
    seu@meta.data$samplename <- c(rep(i, nrow(seu@meta.data)))
    seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = "^MT-")
    seu <- subset(seu, subset = nFeature_RNA > 500 & nCount_RNA > 1000 & percent.mt < 20)
    seu
  })
  return(seu)
}


##########################################################################


#' do_seu_integration (DEPRECATED, USE do_harmony INSTEAD)
#' Perform integration of multiple samples
#'
#' @param names: List of sample names
#' @param experiment: Experiment name
#' @param normalmethod: Normalization method
#' @param scalefactor: Scale factor for cell-level normalization
#' @param hvgs: Number of HVGs to be selected
#' 
#' @keywords preprocessing, integration
#' 
#' @return Multi-sample integrated Seurat object
do_seu_integration <- function(names, experiment, normalmethod, scalefactor, hvgs){
  main_folder <- Sys.getenv("PREPROC_RESULTS_FOLDER")
  seu.list <- lapply(X = names, FUN = function(x) { # Load, normalize and identify variable features for each sample independently
    x <- readRDS(file.path(main_folder,
                                    x,
                                    "preprocessing.R_0000",
                                    paste0(experiment,
                                           ".",
                                           x,
                                           ".before.seu.RDS")))
    x <- NormalizeData(x, normalization.method = normalmethod, scale.factor = scalefactor)
    x <- FindVariableFeatures(x, nfeatures = hvgs)
    })
  features <- SelectIntegrationFeatures(object.list = seu.list) # Select features that are repeatedly variable across samples for integration
  anchors <- FindIntegrationAnchors(object.list = seu.list, anchor.features = features, reference = 1) # Identify anchors
  seu <- IntegrateData(anchorset = anchors) # Create an integrated data assay
  DefaultAssay(seu) <- "integrated" # Set integrated data as default assay
  return(seu)
}


##########################################################################


#' main_preprocessing_analysis
#' Main preprocessing function that performs all the analyses (individually or combining all samples with and without integration)
#'
#' @param name: sample name, or condition if integrate is TRUE
#' @param expermient: experiment name
#' @param input: directory with the single-cell data
#' @param filter: TRUE for using only detected cell-associated barcodes, FALSE for using all detected barcodes
#' @param mincells: min number of cells for which a feature is recorded
#' @param minfeats: min number of features for which a cell is recorded
#' @param minqcfeats: min number of features for which a cell is selected
#' @param percentmt: max percentage of reads mapped to mitochondrial genes for which a cell is selected
#' @param normalmethod: Normalization method
#' @param scalefactor: Scale factor for cell-level normalization
#' @param hvgs: Number of HVGs to be selected
#' @param ndims: Number of PC to be used for clustering / UMAP / tSNE
#' @param resolution: Granularity of the downstream clustering (higher values -> greater number of clusters)
#' @param integrate: FALSE if we don't run integrative analysis, TRUE otherwise
#' @param column: TODO remove var ---- Column with the condition used for subsetting, default NULL just in case integrate is FALSE
#' @param subset: vector of sample names, default NULL just in case integrate is FALSE
#'
#' @keywords preprocessing, main
#' 
#' @return Seurat object
main_preprocessing_analysis <- function(name, experiment, input, filter, mincells, minfeats, minqcfeats,
                                        percentmt, normalmethod, scalefactor, hvgs, ndims, resolution, integrate = FALSE, column = NULL, subset = NULL){

  # Input selection
  
  if (filter == "TRUE"){
    input <- file.path(input, "filtered_feature_bc_matrix")
  } else {
    input <- file.path(input, "raw_feature_bc_matrix")
  }
  
  # Input reading and integration 
  if (identical(integrate, FALSE)) {
    seu <- read_input(name = name, 
                      input = input,
                      mincells = mincells,
                      minfeats = minfeats)
    save_before_seu = TRUE
  } else { # TODO this is supposed to be executed before so it's not needed here
    subsets <- do_subsetting(exp_design = Sys.getenv("exp_design"),
                             column = column)
    seu <- do_integration(names = sample_names,
                          experiment = experiment,
                          normalmethod = normalmethod,
                          scalefactor = scalefactor,
                          hvgs = hvgs)
    save_before_seu = FALSE # For QC (we should already have the Seurat objects before QC)
  }
  
  # QC
  
  seu <- do_qc(name = name,
               experiment = experiment,
               seu = seu,
               minqcfeats = minqcfeats, 
               percentmt = percentmt,
               save_before_seu = save_before_seu)
  
  # Normalization

  if (integrate == FALSE) {
    seu <- NormalizeData(seu, normalization.method = normalmethod, scale.factor = scalefactor)
  }

  # Feature selection
  
  if (integrate == FALSE) {
    seu <- FindVariableFeatures(seu, nfeatures = hvgs)
  }
  
  # Data scaling
  
  seu <- ScaleData(seu, features = rownames(seu))
  
  # Dimensionality reduction (PCA/UMAP/tSNE)
  
  seu <- do_dimred(seu = seu,
                   ndims = ndims)
  
  # Clustering
  
  seu <- do_clustering(seu = seu,
                       ndims = ndims,
                       resolution = resolution)
                       
  # Marker gene selection
  
  do_marker_gene_selection(seu = seu,
                           name = name,
                           experiment = experiment)
  
  # Save final Seurat object
  
  if (identical(integrate, FALSE)){
    saveRDS(seu, paste0(experiment, ".", name, ".seu.RDS"))
  } else {
    saveRDS(seu, file.path(Sys.getenv(PREPROC_RESULTS_FOLDER), paste0(experiment, ".", name, ".seu.RDS")))
  }
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
#' @param resolution: Granularity of the downstream clustering (higher values -> greater number of clusters)
#' @param all_seu: NULL if creating an individual report (daemon 3a). A list of 2 Seurat objects and a matrix of markers if creating a general report (daemon 3b)
#' 
#' @keywords preprocessing, write, report
#' 
#' @return nothing
write_preprocessing_report <- function(name, experiment, template, outdir, intermediate_files, minqcfeats, percentmt, hvgs, resolution, all_seu = NULL, integrate = FALSE){
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
    markers <- readRDS(file.path(main_folder,
                                 name,
                                 "preprocessing.R_0000",
                                 paste0(experiment,
                                        ".",
                                        name,
                                        ".markers.RDS")))
  } else {
    seu <- all_seu[[1]]
    before.seu <- all_seu[[2]]
    markers <- all_seu[[3]]
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


##########################################################################


#' extract_metadata
#' Extract metadata dataframe from Seurat objects
#' 
#' @param seu: Seurat object / list of Seurat objects
#' 
#' @keywords preprocessing, report, metadata
#' 
#' @return Dataframe with metadata
extract_metadata <- function(seu){
  if (!is.list(seu)){
    seu <- seu[[]]
  } else {
    seu <- lapply(seu, "[[")
    seu <- do.call(rbind, seu)
    }
return(seu)
}

##########################################################################

#' make_vln
#' Make Violin plot
#' 
#' @param seu: Seurat object / list of Seurat objects
#' @param feature: metadata feature to plot
#' 
#' @keywords preprocessing, report, plot, violin
#' 
#' @return nothing
make_vln <- function(seu, feature){
  seu <- extract_metadata(seu)
  seu <- seu[, c("orig.ident", feature)]
  colnames(seu)[2] <- "values"
  ggplot() + 
    geom_point(seu,
               mapping = aes(orig.ident, values, fill = orig.ident),
               shape = 21,
               colour = "white",
               size = 2,
               stroke = 0.5,
               position = "jitter",
               alpha = 0.3) +
    geom_violin(seu,
                mapping = aes(orig.ident, values, fill = orig.ident),
                width = 0.5,
                color = "black") +
    xlab(NULL) +
    ylab(feature) +
    labs(fill = NULL) +
    theme_bw()
}

##########################################################################


#' ensure_list
#' Makes sure you have list of Seurat objects (even if you have only one)
#' 
#' @param seu: Seurat object / list of Seurat objects
#' 
#' @keywords preprocessing, report, list
#' 
#' @return List of Seurat objects
ensure_list <- function(seu){
  if (!is.list(seu)){
    seu <- list(seu)
  }
  return(seu)
}

##########################################################################

