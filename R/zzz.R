#' housepi
#'
#' Household Epidemiological Modeling and Inference Rcpp Toolkit
#'
#' Rcpp package for the epidemiological modeling and indference of household survey datasets.
"_PACKAGE"

.onLoad <- function(libname, pkgname) {
    template_dir <- file.path(system.file(package = pkgname), "template")
    pkgbuild::compile_dll(template_dir, debug = FALSE)
}