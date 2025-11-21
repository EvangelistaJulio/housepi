#
#        __
#       / /_  ____  __  __ ____ __  ____  __
#      / __ \/ __ \/ / / / ___/ _ \/ __ \/ /
#     / / / / /_/ / /_/ (__  )  __/ /_/ / /
#    /_/ /_/\____/\____/____/\___/ ____/_/
#                               /_/
#
# housepi - Epidemiological Analysis of Household Outbreaks.
#
# Author: Francesco Camaglia
# Copyright © 2025 Institut Pasteur (Paris), CNRS and INSERM.
#
# Distributed under the terms of the GNU General Public License (GPLv3).
# See the LICENSE file for details.
#
# expose_internal.R
#

# Get template name from DLL file
get_template_from_dll <- function(dll) {
  if (is.null(dll) || is.na(dll) || !file.exists(dll)) NA_character_
  else tools::file_path_sans_ext(basename(dll))
}

# Get Rcpp module from DLL file
get_rcpp_module_from_dll <- function(module_name, dll) {
  dyn.load(dll)
  Rcpp::Module(module_name, PACKAGE = get_template_from_dll(dll))
}

# Get Rcpp module from model object
get_rcpp_module_from_model <- function(module_name, model = NULL) {
  dll <- ifelse(is.null(model), NA_character_,
                model$.__enclos_env__$private$.dll)
  if (is.na(dll)) dll <- blank_model_dll_file()
  get_rcpp_module_from_dll(module_name, dll)
}

rd_test_name_ <- "`character(1)` with the name of the %s."
rd_test_x_ <- "`numeric()` with values."
rd_test_args_ <- "named arguments with corresponding values of the %s"
rd_test_discr_ <- "`logical(1)` to treat definition as discrete."
rd_test_model_ <- "`HousepiModel` object containing the custom registred %s.\\cr
If `NULL` (default) the built-in registry is used."

#' @title
#' New Housepi Function
#' @description
#' Expose an object from housepi registered function.
#' Used to access user defined functions.
#'
#' @param name `r sprintf(rd_test_name_, "function")`
#' @param ... `r sprintf(rd_test_args_, "function")`
#' @param min `r rd_var_min_`
#' @param max `r rd_var_max_`
#' @param discr `r rd_test_discr_`
#' @param model `r rd_test_model_`\cr
#'
#' @return
#' An object of class `FunctionHousepi`.
#'
#' profile <- expose_function("weibull", shape = 2, scale = 1)
#' profile$pdf(seq(0, 10))
#' profile$cdf(seq(0, 10))
#'
#' @export
expose_function <- function(name, ...,
                            min = NA_real_, max = NA_real_,
                            discr = FALSE, model = NULL) {
  args <- unlist(list(...))
  m <- get_rcpp_module_from_model("housepi_function", model)
  m$HousepiRcppFunction$new(name, args = args, min, max, discr)
}

#' @title
#' Test Housepi Distribution
#' @description
#' Expose an object from housepi registered distribution.
#' Used to access user defined distributions.
#'
#' @param name `r sprintf(rd_test_name_, "distribution")`
#' @param ... `r sprintf(rd_test_args_, "distribution")`
#' @param min `r rd_var_min_`
#' @param max `r rd_var_max_`
#' @param discr `r rd_test_discr_`
#' @param model `r rd_test_model_`\cr
#'
#' @return
#' An object of class `DistributionHousepi`.
#'
#' profile <- expose_distribution("weibull", shape = 2, scale = 1)
#' profile$pdf(seq(0, 10))
#' profile$cdf(seq(0, 10))
#'
#' @export
expose_distribution <- function(name, ...,
                                min = NA_real_, max = NA_real_,
                                discr = FALSE, model = NULL) {
  args <- unlist(list(...))
  m <- get_rcpp_module_from_model("housepi_distribution", model)
  m$HousepiRcppDistribution$new(name, args = args, min, max, discr)
}

#' @title
#' Expose Housepi Infectivity Profile
#' @description
#' Expose an object from housepi registered infectivity profile.
#' Used to access user defined infectivity profiles.
#'
#' @param name `r sprintf(rd_test_name_, "infectivity profile")`
#' @param ... `r sprintf(rd_test_args_, "infectivity profile")`
#' @param min `r rd_var_min_`
#' @param max `r rd_var_max_`
#' @param discr `r rd_test_discr_`
#' @param model `r rd_test_model_`\cr
#'
#' @return
#' An object of class `InfecProfileHousepi`.
#'
#' @examples
#' profile <- expose_infec_profile("gamma", shape = 2, scale = 1)
#' profile$pdf(seq(0, 10))
#' profile$cdf(seq(0, 10))
#'
#' @export
expose_infec_profile <- function(name, ...,
                                 min = NA_real_, max = NA_real_,
                                 discr = FALSE, model = NULL) {
  args <- unlist(list(...))
  default_min <- ifelse(discr, 1, 0)
  if (is.na(min)) min <- default_min
  else if (min != default_min)  warning_h(
    "Overwriting default 'min' for an infectivity profile."
  )
  m <- get_rcpp_module_from_model("housepi_infec_profile", model)
  m$HousepiRcppInfecProfile$new(name, args = args, default_min, max, discr)
}
