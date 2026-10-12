## housepi - Epidemiological Analysis of Household Outbreaks.
## Copyright 2024-2026 Institut Pasteur (Paris), France.
## SPDX-License-Identifier: GPL-3.0-or-later
##

#' @keywords internal
#' @details
#' To create a model call [new_housepi_model()]
#' from one of the templated transmission models or a custom model.
#' Provide household data at construction or later with [add_data_wide()].
#' Specify the model features with the decorators [add_parameter()],
#' [add_group_risk()], [add_contact_risk()], [add_infectivity_profile()],
#' [add_latent()] and [add_observation()].
#' Run the inference with [run_inference()], extract draws with
#' [get_parameter_draws()] and summarise them with [summarise_draws()].
#' Run simulations with [run_simulation()],
#' predictive checks with [check_predictive()].
#' Explore available building blocks using [browse_templates()],
#' [browse_distributions()], [browse_functions()]
#' and [browse_infectivity_profiles()].
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
  \(...) with_cpp_conditions(fn(...))
}
#' Define the wrapper for to the .engine Rcpp class.
#' @exportS3Method "$" housepi_engine
#' @noRd
`$.housepi_engine` <- function(x, name) {
  wrap_cpp_member_(.subset2(x, "module"), name)
}

# nocov start
.onLoad <- function(libname, pkgname) {
  pkgbuild::check_build_tools(quiet = TRUE)
  module <- Rcpp::Module(pkgname, PACKAGE = pkgname)
  .engine <<- structure(list(module = module), class = "housepi_engine")
}

.onUnload <- function(libpath) {
  library.dynam.unload("housepi", libpath)
}
# nocov end
