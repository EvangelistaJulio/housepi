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

#' @import Rcpp
#' @import methods
#' @useDynLib housepienvironment, .registration = TRUE
NULL

Rcpp::loadModule("housepi_flags", TRUE)
Rcpp::loadModule("housepi_function", TRUE)
Rcpp::loadModule("housepi_distribution", TRUE)
Rcpp::loadModule("housepi_infec_profile", TRUE)
Rcpp::loadModule("housepi_model_core", TRUE)
