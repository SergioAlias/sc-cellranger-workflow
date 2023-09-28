# Sergio Alías, 20230710
# Last modified 20230922
##############################################################
########################## QC LIBRARY ########################
##############################################################


#' write_qc_report
#' Write QC HTML report
#' 
#' @param name: sample name
#' @param expermient: experiment name
#' @param template: Rmd template
#' @param outdir: output directory
#' @param intermediate_files: directory for saving intermediate files in case pandoc fails
#' @param metrics: metrics file in wide format
#' @param long_metrics: metrics file in long format
#' 
#' @keywords QC, write, report
#' 
#' @return nothing
write_qc_report <- function(name, experiment, template, outdir, intermediate_files, metrics, long_metrics, cellranger_metrics, cellranger_long_metrics){
  int_files <- file.path(outdir, intermediate_files)
  if (!file.exists(int_files)){
    dir.create(int_files)
  }
  rmarkdown::render(template,
                    output_file = file.path(outdir,
                                            paste0(experiment,
                                                   "_",
                                                   name,
                                                   "_QC_report.html")), 
                    clean = TRUE,
                    intermediates_dir = int_files)
}


##########################################################################


#' make_barplot
#' Make Barplot
#' 
#' @param metric_df: DataFrame with the metric table in wide format
#' @param feature: metric to plot
#' 
#' @keywords preprocessing, report, plot, barplot
#' 
#' @return nothing
make_barplot <- function(metric_df, feature){
  metric_df <- metric_df[, c("sample", feature)]
  colnames(metric_df)[2] <- "values"
  ggplot() +
    geom_bar(metric_df,
             mapping = aes(sample, values, fill = sample),
             width = 0.5,
             color = "black",
             stat = "identity") +
    xlab(NULL) +
    ylab(feature) +
    labs(fill = NULL) +
    scale_x_discrete(guide = guide_axis(n.dodge = 2)) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
}


##########################################################################


#' make_hplot
#' Make horizontal plot
#' 
#' @param df: DataFrame with the metric table in long format
#' @param features: metrics to plot
#' 
#' @keywords preprocessing, report, plot, horizontal
#' 
#' @return nothing
make_hplot <- function(metric_df, features){
  metric_df <- metric_df[metric_df[[2]] %in% features,]
  metric_df[[3]] <- as.numeric(metric_df[[3]])
  ggplot() +
    geom_line(metric_df,
             mapping = aes(V1, V3, group = V2, colour = V2)) +
    geom_point(metric_df,
              mapping = aes(V1, V3, group = V2, colour = V2)) +
    xlab(NULL) +
    ylab(NULL) +
    labs(fill = NULL) +
    theme_bw() +
    scale_x_discrete(guide = guide_axis(n.dodge = 2)) +
    scale_y_continuous(limits = c(0, NA)) +
    theme(legend.position="bottom",
          legend.title=element_blank(),
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
}

##########################################################################

