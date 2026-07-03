## housepi - Epidemiological Analysis of Household Outbreaks.
##
## Copyright © 2024-2026 Institut Pasteur (Paris), CNRS and INSERM.
## Distributed under the terms of the GNU General Public License (GPLv3).
## See the LICENSE file for details.
##

#' @keywords internal
#' @details
#' Mathematical modeling is a valuable tool for estimating
#' key epidemiological parameters from household transmission surveys.
#' The relative ease of collecting household data makes it a readily available
#' source of information during an epidemic early phases
#' (e.g. longitudinal studies).
#' Additionally, retrospective serological surveys can provide critical insights
#' into the history of transmission within a population (final size studies).
#' The package "housepi" aims at supporting straightforward analysis
#' and advanced epidemiological investigations from household outbreaks.
"_PACKAGE"

#' @useDynLib housepi, .registration = TRUE
#' @importFrom RcppParallel RcppParallelLibs
#' @importFrom Rdpack reprompt
NULL

# Thin proxy over the Rcpp module.
# Its `$` members are wrapped alltogether (see `with_cpp_conditions`).
.engine <- NULL

# Fetch a raw C++ member via `$`, bypassing any condition interceptor.
cpp_member_ <- function(x, name) {
  eval(as.call(list(as.name("$"), x, as.name(name))))
}
# Fetch a C++ member via `$`, wrapping its *call* (not just the fetch) so
# conditions raised while the member runs route through the interceptor.
wrap_cpp_member_ <- function(x, name) {
  fn <- cpp_member_(x, name)
  function(...) with_cpp_conditions(fn(...))
}
#' Define the wrapper for to the .engine Rcpp class.
#' @noRd
`$.housepi_engine` <- function(x, name) {
  wrap_cpp_member_(.subset2(x, "module"), name)
}

# nocov start'
.onLoad <- function(libname, pkgname) {
  pkgbuild::check_build_tools(quiet = TRUE)
  module <- Rcpp::Module(pkgname, PACKAGE = pkgname)
  .engine <<- structure(list(module = module), class = "housepi_engine")
  registerS3method("$", "housepi_engine", `$.housepi_engine`,
                   envir = asNamespace(pkgname))
  registerS3method("$", "housepi_exposed", `$.housepi_exposed`,
                   envir = asNamespace(pkgname))
  registerS3method("print", "housepi_exposed", print.housepi_exposed,
                   envir = asNamespace(pkgname))
  # default label for the susceptible status (code 0) in character status output
  op <- options()
  if (is.null(op[["housepi.susceptible_label"]])) {
    options(housepi.susceptible_label = "susceptible")
  }
}

.onUnload <- function(libpath) {
  library.dynam.unload("housepi", libpath)
}
# nocov end
