# Sergio Alías, 20230606
# Last modified 20230606

##########################################################################
########################## PRE-PROCESSING LIBRARY ########################
##########################################################################


#' Append experiments names with a Hello World message
#' @param name: expermient name
#' @param message: message to append. Default: "Hello World"
#' @keywords method
#' @return character vector
append_message <- function(name, message = "Hello World") {paste(name, message)}