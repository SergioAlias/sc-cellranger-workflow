# Sergio Alías, 20230606
# Last modified 20230628

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
  
  # seu[["orig.ident"]] <- rep.int(name, times = length(Cells(seu)))
  mtx <- NULL
  return(seu)
}


##########################################################################


#' do_qc
#' Perform Quality Control
#'
#' @param seu: Seurat object
#' @param minqcfeats: min number of features for which a cell is selected
#' @param percentmt: max percentage of reads mapped to mitochondrial genes for which a cell is selected
#' 
#' @keywords preprocessing, qc
#' 
#' @return Seurat object
do_qc <- function(seu, minqcfeats, percentmt){
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
  
  return(seu)
}


##########################################################################


#' do_fsel
#' Perform Feature Selection
#'
#' @param seu: Seurat object
#' @param minqcfeats: min number of features for which a cell is selected
#' @param percentmt: max percentage of reads mapped to mitochondrial genes for which a cell is selected
#' 
#' @keywords preprocessing, qc
#' 
#' @return Seurat object
do_fsel <- function(seu, minqcfeats, percentmt){
  #### Feature selection ####
  

  
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
#' 
#' @keywords preprocessing, main
#' 
#' @return Seurat object
main_preprocessing_analysis <- function(name, experiment, input, filter, mincells, minfeats, minqcfeats, percentmt, normalmethod, scalefactor){

  if (filter == "TRUE"){
    input <- file.path(input, "filtered_feature_bc_matrix")
  } else {
    input <- file.path(input, "raw_feature_bc_matrix")
  }
  
  seu <- read_input(name = name, 
                    input = input,
                    mincells = mincells,
                    minfeats = minfeats)
  
  seu <- do_qc(seu = seu,
               minqcfeats = minqcfeats, 
               percentmt = percentmt)
  
  seu <- NormalizeData(seu, normalization.method = normalmethod, scale.factor = scalefactor)
  
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
#' @param all_seu: NULL if creating an indiviual report (daemon 3a). A Seurat object if creating a general report (daemon 3b)
#' 
#' @keywords preprocessing, write, report
#' 
#' @return nothing
write_preprocessing_report <- function(name, experiment, template, outdir, intermediate_files, minqcfeats, percentmt, all_seu = NULL){
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
  } else {
    seu <- all_seu
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