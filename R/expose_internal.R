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

#' @importFrom Rcpp Module

# Get template name from DLL file
get_template_from_dll <- function(dll) {
  if (is.null(dll) || is.na(dll) || !file.exists(dll)) NA_character_
  else tools::file_path_sans_ext(basename(dll))
}

# Get Rcpp module from dll
get_rcpp_module_from_dll <- function(module_name, dll) {
  dyn.load(dll)
  Rcpp::Module(module_name, PACKAGE = get_template_from_dll(dll))
}

# Get Rcpp module from model
get_rcpp_module_from_model <- function(module_name, hhmodel) {
  if (is.null(hhmodel)) {
    Rcpp::Module(module_name, PACKAGE = "housepi")
  } else {
    get_rcpp_module_from_dll(module_name, hhmodel$.__enclos_env__$private$.dll)
  }
}

rd_test_name_ <- "`character(1)` with the name of the %s."
rd_test_x_ <- "`numeric()` with values."
rd_test_args_ <- "named arguments with corresponding values of the %s"
rd_test_discr_ <- "`logical(1)` to treat definition as discrete."
rd_test_model_ <- "`HousepiModel` object containing the custom registred %s.\\cr
If `NULL` (default) the built-in registry is used."

#' @title
#' Expose Housepi Function
#' @description
#' Expose an object from housepi registered function.
#' Used to access user defined functions.
#'
#' @param name `r sprintf(rd_test_name_, "function")`
#' @param ... `r sprintf(rd_test_args_, "function")`
#' @param min `r rd_var_min_`
#' @param max `r rd_var_max_`
#' @param discr `r rd_test_discr_`
#' @param hhmodel `r rd_test_model_`
#'
#' @return
#' An object of class `FunctionHousepi`.
#'
#' fnc <- expose_function("pow", x2 = 3)
#' fnc$func(seq(0, 5))
#' fnc$prim(seq(0, 5))
#'
#' @export
expose_function <- function(name, ...,
                            min = NA_real_, max = NA_real_,
                            discr = FALSE, hhmodel = NULL) {
  m <- get_rcpp_module_from_model("housepi_function", hhmodel)
  m$HousepiRcppFunction$new(name, args = unlist(list(...)),
                            min = min, max = max, discr = discr)
}

#' @title
#' Expose Housepi Distribution
#' @description
#' Expose an object from housepi registered distribution.
#' Used to access user defined distributions.
#'
#' @param name `r sprintf(rd_test_name_, "distribution")`
#' @param ... `r sprintf(rd_test_args_, "distribution")`
#' @param min `r rd_var_min_`
#' @param max `r rd_var_max_`
#' @param discr `r rd_test_discr_`
#' @param hhmodel `r rd_test_model_`
#'
#' @return
#' An object of class `DistributionHousepi`.
#'
#' dst <- expose_distribution("weibull", shape = 2, scale = 1)
#' dst$pdf(seq(0, 5))
#' dst$cdf(seq(0, 5))
#'
#' @export
expose_distribution <- function(name, ...,
                                min = NA_real_, max = NA_real_,
                                discr = FALSE, hhmodel = NULL) {
  m <- get_rcpp_module_from_model("housepi_distribution", hhmodel)
  m$HousepiRcppDistribution$new(name, args = unlist(list(...)),
                                min = min, max = max, discr = discr)
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
#' @param hhmodel `r rd_test_model_`
#'
#' @return
#' An object of class `InfecProfileHousepi`.
#'
#' @examples
#' ipf <- expose_infec_profile("gamma", shape = 2, scale = 1)
#' ipf$pdf(seq(0, 5))
#' ipf$cdf(seq(0, 5))
#'
#' @export
expose_infec_profile <- function(name, ...,
                                 min = NA_real_, max = NA_real_,
                                 discr = FALSE, hhmodel = NULL) {
  default_min <- ifelse(discr, 1, 0)
  if (is.na(min)) min <- default_min
  else if (min != default_min) warning_h(
    "Overwriting default 'min' for an infectivity profile."
  )
  m <- get_rcpp_module_from_model("housepi_infec_profile", hhmodel)
  m$HousepiRcppInfecProfile$new(name, args = unlist(list(...)),
                                min = min, max = max, discr = discr)
}
