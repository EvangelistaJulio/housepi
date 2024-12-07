#' housepi
#'
#' Household Epidemiological Modeling and Inference Rcpp Toolkit
#'
#' Rcpp package for epidemiological modeling and inference of household survey datasets.
"_PACKAGE"

.onLoad <- function(libname, pkgname) {
  re_compile_library(pkgname, force = FALSE)
}