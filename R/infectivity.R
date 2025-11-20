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
    deep_clone = function(name, value) helper_clone(value),

    # Definitions required by augmented variable definitions
    # according to the "type" of the definition (A, D, P).
    .find_def_required = function(type) {
      output <- c()
      if (length(self$defs) > 0) {
        for (typed_col in c("name", "max")) {
          for (x in self$defs[[typed_col]]) {
            output <- c(output, DefinitionWalker$new(x, type)$x)
          }
        }
      }
      # remove duplicates and add source name
      set_source_name_(unique(output), "infectivity profile")
    }
  ),

  ############################################################################
  ################################## ACTIVE ##################################
  ############################################################################

  active = list(
    #' @field data_required `character()`\cr
    #' Dataset entries required by infectivity profile definition.
    data_required = function(value) {
      if (!missing(value)) stop_h("'data_required' is read only.")
      private$.find_def_required("D")
    },

    #' @field param_required `character()`\cr
    #' Parameters required by infectivity profile definition.
    param_required = function(value) {
      if (!missing(value)) stop_h("'param_required' is read only.")
      private$.find_def_required("P")
    },

    #' @field augm_required `character()`\cr
    #' Augmented variables required by infectivity profile definition.
    augm_required = function(value) {
      if (!missing(value)) stop_h("'augm_required' is read only.")
      private$.find_def_required("A")
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

    #' @inherit add_infec_profile
    add = function(name, status = NA, discr = FALSE, ..., max = NULL) {
      # check infectivity profile name
      if (is.na(name)) stop_h("Incomplete definition: missing name.")
      args <- unlist(list(...))
      self$defs <- rbind(self$defs, data.frame(list(
        # load name/arguments
        name = ifelse(
          is.null(args), name,
          DefinitionEncoder$new(as.list(args))$combine_to_call(name)
        ),
        status = status,
        min = ifelse(isTRUE(discr), 1, 0),
        # load maximum value
        max = DefinitionEncoder$new(max)$x,
        discr = isTRUE(discr)
      )))

      invisible(self)
    }
  )
)
