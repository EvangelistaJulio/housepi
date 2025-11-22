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
# zzz.R
#

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
#' into the history of transmission within a population (finite size studies).
#' The package aims to combine ease of use with flexibility, supporting both
#' straightforward analyses and advanced epidemiological investigations.
"_PACKAGE"

#' @import BH
#' @useDynLib housepi, .registration = TRUE
NULL

# UPDATE: the goal is to compile all template models upon loading the package
.onLoad <- function(libname, pkgname) {
  compile_template_library(pkgname)
}

# @description
# Compile all template objects of the housepi library.
#
# @param pkgname `character(1)` with the name of the package.
compile_template_library <- function(pkgname) {
  # compile each template model
  for (m in list.files(template_models_dir(), full.names = TRUE)) {
    dll <- template_dll_file(m)
    if (isFALSE(file.exists(dll)) || untracked_changes(m, dll)) {
      msg <- paste0("COMPILING ", toupper(spaced_template_choice(basename(m))))
      cat("\n")
      cli::cat_rule(msg, col = "cyan")
      HousepiGenerator$new(m, as_template = TRUE, quiet = FALSE)
    }
  }
}
