#      __
#     / /_  ____  __  __ ____ __  ____  __
#    / __ \/ __ \/ / / / ___/ _ \/ __ \/ /
#   / / / / /_/ / /_/ (__  )  __/ /_/ / /
#  /_/ /_/\____/\____/____/\___/ ____/_/
#                             /_/
#
# Copyright (C) 2025 CAMAGLIA Francesco - Institut Pasteur, Paris (France).
#
# expose_internal.R
#

# Get template name from DLL file
get_template_from_dll <- function(dll) {
  if (is.null(dll) || is.na(dll) || !file.exists(dll)) NA
  else tools::file_path_sans_ext(basename(dll))
}

# Get Rcpp module from DLL file
get_rcpp_module_from_dll <- function(module_name, dll) {
  dyn.load(dll)
  Rcpp::Module(module_name, PACKAGE = get_template_from_dll(dll))
}

# Get Rcpp module from model object
get_rcpp_module_from_model <- function(module_name, model = NULL) {
  dll <- ifelse(is.null(model), NA, model$path_to_dll)
  if (is.na(dll)) {
    housepi_dir <- system.file("template", package = "housepi")
    dll <- file.path(housepi_dir, "src", dll_file_ext("housepienvironment"))
  }
  get_rcpp_module_from_dll(module_name, dll)
}

#' @title
#' Test Housepi Function
#' @description
#' Call housepi registered function.
#' This function internally calls a large C++ class,
#' thus it is not optimized for repeated  operations.
#'
#' @param name `character(1)`\cr
#' Name of the function.
#' @param x `numeric()`\cr
#' Input vector.
#' @param args `numeric()`\cr
#' Additional arguments.
#' @param min `double(1)`\cr
#' Value of the minimum.
#' @param max `double(1)`\cr
#' Value of the maximum.
#' @param model `HouseholdModel`\cr
#' Custom household model to load from.
#'
#' @return
#' `numeric()`
#'
#' @export
test_function <- function(name, x, args = NULL,
                          min = NA, max = NA, model = NULL) {
  # note: for simple function no reshaping -> discr = FALSE
  get_rcpp_module_from_model(
    "housepi_function", model
  )$FunctionHousepi$new(
    name, FALSE, args, c(min = min, max = max)
  )$call(x)
}

#' @title
#' Test Housepi Distribution
#' @description
#' Call housepi registered distribution.
#' This function internally calls a large C++ class,
#' thus it is not optimized for repeated  operations.
#'
#' @param name `character(1)`\cr
#' Name of the distribution.
#' @param x `numeric()`\cr
#' Input vector.
#' @param discr `logical(1)`\cr
#' If `TRUE` discretize the distribution.
#' @param args `numeric()`\cr
#' Additional arguments.
#' @param min `double(1)`\cr
#' Value of the minimum.
#' @param max `double(1)`\cr
#' Value of the maximum.
#' @param model `HouseholdModel`\cr
#' Custom household model to load from.
#'
#' @return
#' `numeric()`
#'
#' @export
test_distribution <- function(name, x, discr = FALSE, args = NULL,
                              min = NA, max = NA, model = NULL) {
  get_rcpp_module_from_model(
    "housepi_distribution", model
  )$DistributionHousepi$new(
    name, discr, args, c(min = min, max = max)
  )$call(x)
}

#' @title
#' Test Housepi Infectivity Profile
#' @description
#' Call housepi registered infectivity profile.
#' This function internally calls a large C++ class,
#' thus it is not optimized for repeated  operations.
#'
#' @param name `character(1)`\cr
#' Name of the infectivity profile
#' @param x `numeric()`\cr
#' Input vector.
#' @param discr `logical(1)`\cr
#' If `TRUE` discretize the infectivity profile.
#' @param args `numeric()`\cr
#' Additional arguments.
#' @param max `double(1)`\cr
#' Value of the maximum.
#' @param model `HouseholdModel`\cr
#' Custom household model to load from.
#'
#' @return
#' `numeric()`
#'
#' @export
test_infec_profile <- function(name, x, discr = FALSE, args = NULL,
                               max = NA, model = NULL) {
  # UPDATE: the minimum should be imposed from within
  get_rcpp_module_from_model(
    "housepi_infec_profile", model
  )$InfecProfileHousepi$new(
    name, discr, args, c(min = ifelse(discr, 1, 0), max = max)
  )$call(x)
}
