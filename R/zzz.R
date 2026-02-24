## housepi - Epidemiological Analysis of Household Outbreaks.
##
## Copyright © 2026 Institut Pasteur (Paris), CNRS and INSERM.
## Distributed under the terms of the GNU General Public License (GPLv3).
## See the LICENSE file for details.
##

#' housepi
#'
#' Epidemiological Analysis of Household Outbreaks.
#'
#' We introduce "housepi", an open source R package for the statistical analysis
#' and modeling of epidemiological data from household outbreaks.
#' Mathematical modeling has proven to be a valuable tool for estimating
#' key epidemiological parameters from household transmission surveys.
#' The relative ease of sampling household data makes it a readily available
#' source of information during the early phases of an epidemic
#' (e. g. longitudinal studies).
#' Additionally, retrospective serological surveys can provide critical insights
#' into the history of transmission within a population (final size studies).
#' The package aims to combine ease of use with flexibility, supporting both
#' straightforward analyses and advanced epidemiological investigations.
"_PACKAGE"

#' @useDynLib housepi, .registration = TRUE
NULL

.engine <- NULL

.onLoad <- function(libname, pkgname) {
  dll_path <- getLoadedDLLs()[[pkgname]][["path"]]
  if (Sys.info()["sysname"] %in% c("Darwin", "Linux")) {
    dyn.unload(dll_path)
    dyn.load(dll_path, local = FALSE, now = TRUE)
  }
  .engine <<- Rcpp::Module(pkgname, PACKAGE = pkgname)
}
