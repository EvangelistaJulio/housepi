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
# See the COPYRIGHT file for details.
#
# Distributed under the terms of the GNU General Public License (GPLv3).
# See the COPYING file for details.
#
# zzz.R
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
  # UPDATE: this should be droppable
  # # compile the core
  # withr::with_envvar(
  #   get_cpp_environment_flags(),
  #   pkgbuild::compile_dll(system.file("template", package = pkgname),
  #                         quiet = FALSE, debug = FALSE)
  # )

  # # compile each template model
  # for (m in list.files(template_models_dir(), full.names = TRUE)) {
  #   dll <- template_dll_file(m)
  #   if (isFALSE(file.exists(dll)) || untracked_changes(m, dll)) {
  #     msg <- paste0("COMPILING ", toupper(spaced_template_choice(basename(m))))
  #     cat("\n")
  #     cli::cat_rule(msg, col = "cyan")
  #     HousepiGenerator$new(m, as_template = TRUE, quiet = FALSE)
  #   }
  # }
}
