## housepi - Epidemiological Analysis of Household Outbreaks.
##
## Copyright © 2024-2026 Institut Pasteur (Paris), CNRS and INSERM.
## Distributed under the terms of the GNU General Public License (GPLv3).
## See the LICENSE file for details.
##

#' @keywords internal
#' @details
#' We introduce "housepi", an open source R package for the statistical analysis
#' and modeling of epidemiological data from household outbreaks.
#' Mathematical modeling has proven to be a valuable tool for estimating
#' key epidemiological parameters from household transmission surveys.
#' The relative ease of sampling household data makes it a readily available
#' source of information during the early phases of an epidemic
#' (e.g. longitudinal studies).
#' Additionally, retrospective serological surveys can provide critical insights
#' into the history of transmission within a population (final size studies).
#' The package aims to combine ease of use with flexibility, supporting both
#' straightforward analyses and advanced epidemiological investigations.
"_PACKAGE"

#' @useDynLib housepi, .registration = TRUE
#' @importFrom RcppParallel RcppParallelLibs
#' @importFrom Rdpack reprompt
NULL

.engine <- NULL

# nocov start
.onLoad <- function(libname, pkgname) {
  pkgbuild::check_build_tools(quiet = TRUE)
  .engine <<- Rcpp::Module(pkgname, PACKAGE = pkgname)
}

.onUnload <- function(libpath) {
  library.dynam.unload("housepi", libpath)
}
# nocov end
