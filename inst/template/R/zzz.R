#' @useDynLib hemirt12345, .registration = TRUE
#' @import methods Rcpp
#' @exportPattern ^[[:alpha:]]+
NULL

#' hemirt12345
#'
#' Household Epidemiological Modeling and Inference Rcpp Toolkit
#'
#' A first draft for what could become a package. 
"_PACKAGE"

Rcpp::loadModule("model_hemirt12345", TRUE)
