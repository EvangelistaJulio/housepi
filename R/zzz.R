#
# zzz.R
#
# Copyright (C) 2025 CAMAGLIA Francesco - Institut Pasteur, Paris (France).
#

#' housepi
#'
#' Household Epidemiological Modeling and Inference Rcpp Toolkit
#'
#' Rcpp package for epidemiological modeling and inference
#' of household survey datasets.
"_PACKAGE"

.onLoad <- function(libname, pkgname) {
  compile_template_library(pkgname)
}

# @description
# Compile all template objects of the housepi library.
#
# @param pkgname `character(1)`\cr
# Name of the package.
compile_template_library <- function(pkgname) {
  # compile the core
  housepi_dir <- file.path(system.file(package = pkgname), "template")
  pkgbuild::compile_dll(housepi_dir, debug = FALSE, quiet = FALSE)

  # compile each template model
  for (m in list.files(template_models_dir(), full.names = TRUE)) {
    dll <- template_dll_file(m)
    if (isFALSE(file.exists(dll)) || untracked_changes(m, dll)) {
      HousepiGenerator$new(m,
                           as_template = TRUE,
                           force = FALSE,
                           quiet_pkgbuild = FALSE)$generate()
    }
  }
}
