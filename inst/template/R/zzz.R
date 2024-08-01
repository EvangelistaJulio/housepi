#' @useDynLib housepi12345, .registration = TRUE
#' @import methods Rcpp
#' @exportPattern ^[[:alpha:]]+
NULL

#' Modeling environment
#'
#' @description
#' Type ?Modelhousepi12345 to retrieve all information about the modeling class.
"_PACKAGE"

Rcpp::loadModule("modelhousepi12345_cpp", TRUE)
