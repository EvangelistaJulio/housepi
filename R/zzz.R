#' housepi
#'
#' Household Epidemiological Modeling and Inference Rcpp Toolkit
#'
#' Rcpp package for the epidemiological modeling and indference of household survey datasets.
"_PACKAGE"

.onLoad <- function(libname, pkgname) {
  re_compile_library(pkgname, force = FALSE)
}