# Sergio Alías, 20230606
# Last modified 20230608

##########################################################################
########################## PRE-PROCESSING LIBRARY ########################
##########################################################################


#' append_message
#' Append experiments names with a Hello World message
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
#' @param aux_plots: directory for plots implemented by Seurat
#' @param pdf_prefix: prefix for PDF files containing Seurat plots
#' @param seu: Seurat object
#' @param maxfeats: max number of features for which a cell is recorded
#' @param maxcounts: max number of counts for which a cell is recorded
#' @param percentmt: 
#' 
#' @keywords preprocessing, qc
#' 
#' @return Seurat object
do_qc <- function(aux_plots, pdf_prefix, seu, maxfeats, maxcounts, percentmt){
  #### QC ####
  
  ##### Reads mapped to mitochondrial genes #####
  
  # In human ENSEMBL gene symbol annotations, mitochondrial genes are annotated starting with a MT- string
  
  seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = "^MT-")
  
  # We can define ribosomal proteins (their names start with RPS or RPL)
  
  seu[["percent.rb"]] <- PercentageFeatureSet(seu, pattern = "^RP[SL]")
  
  # Violin plot before filtering
  
  pdf(file.path(aux_plots, paste0(pdf_prefix, "VlnPlot_before.pdf")))
  VlnPlot(seu,
          features = c('nFeature_scRNAseq', 'nCount_scRNAseq', 'percent.mt', 'percent.rb'),
          ncol = 4)
  dev.off()

  ##### Filtering out cells #####
  
  # seu[['QC']] <- ifelse(seu@meta.data$Is_doublet == 'True','Doublet','Pass')
  seu[['QC']] <- ifelse(seu@meta.data$nFeature_scRNAseq < maxfeats & seu@meta.data$QC == 'Pass','Low_nFeature',seu@meta.data$QC)
  seu[['QC']] <- ifelse(seu@meta.data$nFeature_scRNAseq < maxfeats & seu@meta.data$QC != 'Pass' & seu@meta.data$QC != 'Low_nFeature',paste('Low_nFeature',seu@meta.data$QC,sep = ','),seu@meta.data$QC)
  seu[['QC']] <- ifelse(seu@meta.data$percent.mt > percentmt & seu@meta.data$QC == 'Pass','High_MT',seu@meta.data$QC)
  seu[['QC']] <- ifelse(seu@meta.data$nFeature_scRNAseq < maxfeats & seu@meta.data$QC != 'Pass' & seu@meta.data$QC != 'High_MT',paste('High_MT',seu@meta.data$QC,sep = ','),seu@meta.data$QC)
  # table(seu[['QC']])
  
  pdf(file.path(aux_plots, paste0(pdf_prefix, "VlnPlot_after.pdf")))
  VlnPlot(seu,
          features = c('nFeature_scRNAseq', 'nCount_scRNAseq', 'percent.mt', 'percent.rb'),
          ncol = 4)
  dev.off()
  
  return(seu)
}


##########################################################################


#' main_preprocessing_analysis
#' Main preprocessing function that performs all the analyses
#'
#' @param aux_plots: directory for plots implemented by Seurat
#' @param name: sample name
#' @param expermient: experiment name
#' @param input: directory with the single-cell data
#' @param mincells: min number of cells for which a feature is recorded
#' @param minfeats: min number of features for which a cell is recorded
#' 
#' @keywords preprocessing, main
#' 
#' @return Seurat object
main_preprocessing_analysis <- function(aux_plots, name, experiment, input, mincells, minfeats, maxfeats, maxcounts, percentmt){
  
  if (!file.exists(aux_plots)){
    dir.create(aux_plots)
  }
  
  pdf_prefix <- paste0(experiment, "_", name, "_")
  
  seu <- read_input(name = name, 
                    input = input,
                    mincells = mincells,
                    minfeats = minfeats)
  
  seu <- do_qc(aux_plots = aux_plots,
               pdf_prefix = pdf_prefix,
               seu = seu,
               maxfeats = maxfeats, 
               maxcounts = maxcounts, 
               percentmt = percentmt)
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
#' 
#' @keywords preprocessing, write, report
#' 
#' @return nothing
write_preprocessing_report <- function(name, experiment, template, outdir, intermediate_files){
  int_files <- file.path(outdir, intermediate_files)
  if (!file.exists(int_files)){
    dir.create(int_files)
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