#' @useDynLib housepi12345, .registration = TRUE
#' @import methods Rcpp
#' @exportPattern ^[[:alpha:]]+
NULL

#' housepi12345
#'
#' Household Epidemiological Modeling and Inference Rcpp Toolkit
#'
#' A first draft for what could become a package. 
"_PACKAGE"

Rcpp::loadModule("modelhousepi12345", TRUE)
