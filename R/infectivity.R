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
# infectivity.R
#

############################################################################
############################################################################
########################### Housepi Infec Profile ##########################
############################################################################
############################################################################

#' @title
#' Housepi Infectivity Profile
#' @description
#' Auxiliary class to deal with infectivity profile.
#'
#' @keywords internal
#' @noRd
HousepiInfecProfile <- R6::R6Class(
  classname = "HousepiInfecProfile",
  cloneable = TRUE,

  ############################################################################
  ################################## PRIVATE #################################
  ############################################################################

  private = list(
    deep_clone = function(name, value) helper_clone(value)
  ),

  ############################################################################
  ################################## ACTIVE ##################################
  ############################################################################

  active = list(
    #' @field augm_required `character()`\cr
    #' Augmented variables required by infectivity profile definition.
    augm_required = function(value) {
      if (!missing(value)) stop_h("'augm_required' is read only.")
      infec_prof_find_required_def(self$defs, "A")
    },

    #' @field data_required `character()`\cr
    #' Dataset entries required by infectivity profile definition.
    data_required = function(value) {
      if (!missing(value)) stop_h("'data_required' is read only.")
      infec_prof_find_required_def(self$defs, "D")
    },

    #' @field group_rates_required `character()`\cr
    #' Find group rates required by infectivity profile definition.
    group_rates_required = function(value) {
      if (!missing(value)) stop_h("'group_rates_required' is read only.")
      infec_prof_find_required_def(self$defs, "G")
    },

    #' @field param_required `character()`\cr
    #' Parameters required by infectivity profile definition.
    param_required = function(value) {
      if (!missing(value)) stop_h("'param_required' is read only.")
      infec_prof_find_required_def(self$defs, "P")
    }
  ),

  ############################################################################
  ################################## PUBLIC ##################################
  ############################################################################

  public = list(
    #' @field defs `data.frame`\cr
    #' Infectivity profile definition data.frame.
    defs = NULL,

    #' @description
    #' New object of class
    #' [HousepiInfecProfile][housepi::HousepiInfecProfile].
    initialize = function() {
      self$defs <- NULL
      invisible(self)
    },

    #
    add_quo = function(name, status = NA_integer_,
                       discr = FALSE, ..., max = NULL) {
      # load non-quosures
      name <- housepi_smart_unwrap(name, "name")
      status <- housepi_smart_unwrap(status, "status")
      discr <- housepi_smart_unwrap(discr, "discr")
      # load quosures
      max <- housepi_extractor(!!max)
      # load dots # UPDATE: not super elegant...
      dots <- list(...)
      # flatten if single quosures list (from enquos)
      if (length(dots) == 1 && rlang::is_quosures(dots[[1]])) dots <- dots[[1]]
      args <- housepi_extractor(!!!dots)

      # check infectivity profile name
      if (is.na(name)) stop_h("Incomplete definition: missing name.")
      assert_valid_status(status)

      self$defs <- rbind(self$defs, data.frame(list(
        name = ifelse(is.null(args), name, housepi_combine_to_call(name, args)),
        status = status,
        min = ifelse(isTRUE(discr), 1, 0),
        max = ifelse(is.null(max), NA_character_, max),
        discr = isTRUE(discr)
      )))

      invisible(self)
    },

    #' @inherit add_infec_profile
    add = function(name, status = NA_integer_,
                   discr = FALSE, ..., max = NULL) {
      self$add_quo(name = name, status = status,
                   discr = discr, rlang::enquos(...),
                   max = rlang::enquo(max))
    }
  )
)
