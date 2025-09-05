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

rd_test_disclaimer_ <- "This function is meant for testing only and not for repeated operations."
rd_test_name_ <- "`character(1)` with the name of the %s."
rd_test_x_ <- "`numeric()` with values."
rd_test_prim_ <- "`logical(1)` to consider %s instead.\\cr
Default is `FALSE`."
rd_test_args_ <- "named `numeric()` or `list` with custom arguments."
rd_test_discr_ <- "`logical(1)` to treat definition as discrete."
rd_test_model_ <- "`HousepiModel` object eventually containing the custom registred %s.\cr
If `NULL` (default) the built-in registry is used."

#' @title
#' Test Housepi Function
#' @description
#' Call housepi registered function.
#' `r rd_test_disclaimer_`
#'
#' @param name `r sprintf(rd_test_name_, "function")`
#' @param x `r rd_test_x_`
#' @param prim `r sprintf(rd_test_prim_, "primitive")`
#' @param args `r rd_test_args_`
#' @param min `r rd_var_min_`
#' @param max `r rd_var_max_`
#' @param model `r rd_test_model_`\cr
#'
#' @return
#' `numeric()`
#'
#' @export
test_function <- function(name, x, prim = FALSE, args = NULL,
                          min = NA, max = NA, discr = FALSE, model = NULL) {
  # note: for simple function no reshaping -> discr = FALSE
  get_rcpp_module_from_model(
    "housepi_function", model
  )$FunctionHousepi$new(
    name, prim, args, c(min = min, max = max), discr
  )$call(x)
}

#' @title
#' Test Housepi Distribution
#' @description
#' Call housepi registered distribution.
#' `r rd_test_disclaimer_`
#'
#' @param name `r sprintf(rd_test_name_, "distribution")`
#' @param x `r rd_test_x_`
#' @param cumul `r sprintf(rd_test_prim_, "cumulative")`
#' @param args `r rd_test_args_`
#' @param min `r rd_var_min_`
#' @param max `r rd_var_max_`
#' @param discr `r rd_test_discr_`
#' @param model `r rd_test_model_`\cr
#'
#' @return
#' `numeric()`
#'
#' @export
test_distribution <- function(name, x, cumul = FALSE, args = NULL,
                              min = NA, max = NA, discr = FALSE, model = NULL) {
  get_rcpp_module_from_model(
    "housepi_distribution", model
  )$DistributionHousepi$new(
    name, cumul, args, c(min = min, max = max), discr
  )$call(x)
}

#' @title
#' Test Housepi Infectivity Profile
#' @description
#' Call housepi registered infectivity profile.
#' `r rd_test_disclaimer_`
#'
#' @param name `r sprintf(rd_test_name_, "infectivity profile")`
#' @param x `r rd_test_x_`
#' @param cumul `r sprintf(rd_test_prim_, "cumulative")`
#' @param args `r rd_test_args_`
#' @param max `r rd_var_max_`
#' @param discr `r rd_test_discr_`
#' @param model `r rd_test_model_`\cr
#'
#' @return
#' `numeric()`
#'
#' @export
test_infec_profile <- function(name, x, cumul = FALSE, args = NULL,
                               max = NA, discr = FALSE, model = NULL) {
  # UPDATE: the minimum should be imposed from within
  get_rcpp_module_from_model(
    "housepi_infec_profile", model
  )$InfecProfileHousepi$new(
    name, cumul, args, c(min = ifelse(discr, 1, 0), max = max), discr
  )$call(x)
}
