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

    #' @inherit add_infec_profile
    add_quo = function(func, status = NA_integer_,
                       discr = FALSE, ..., max = NULL) {
      # load non-quosures
      func <- housepi_smart_unwrap(func, "func")
      status <- housepi_smart_unwrap(status, "status")
      discr <- housepi_smart_unwrap(discr, "discr")
      # load quosures
      max <- housepi_extractor(!!max)
      # load dots # UPDATE: not super elegant...
      dots <- list(...)
      # flatten if single quosures list (from enquos)
      if (length(dots) == 1 && rlang::is_quosures(dots[[1]])) dots <- dots[[1]]
      func_args <- housepi_extractor(!!!dots)
      # ensure the argument name doesn't mistakenly be interpreted as argument
      func_args <- func_args[names(func_args) != "name"]

      # check infectivity profile function
      if (is.na(func)) stop_h("Incomplete definition: missing function name.")
      assert_valid_status(status)

      self$defs <- rbind(self$defs, data.frame(list(
        func = ifelse(is.null(func_args), func,
                      housepi_combine_to_call(func, func_args)),
        status = status,
        min = ifelse(isTRUE(discr), 1, 0),
        max = ifelse(is.null(max), NA_character_, max),
        discr = isTRUE(discr)
      )))

      # assert no multiple definitions have been provided
      if (length(unique(self$defs$status)) != nrow(self$defs)) {
        self$initialize()
        stop_h(
          "Found multiple infectivity profile definitions for the same status.",
          "Falling back to template definition."
        )
      }

      invisible(self)
    },

    #' @inherit load_infec_profile
    load_quosures = function(xquos, append = FALSE) {
      # UPDATE: append requires checking existing definitions.
      # I can store old ones, compare them and then drop at the end
      # remove previous definitions if append is false
      if (isFALSE(append)) self$initialize()
      self <- housepi_load_defin_quos(self, "func", xquos)

      invisible(self)
    },

    #' @description
    #' Returns a readable definition object.
    showcase = function() {
      if (length(self$defs) > 0) {
        output <- housepi_defin_showcaser(self$defs)

        # append func_args if not present
        i <- which(names(output) == "func")
        if (names(output)[i + 1] != "func_args") {
          output <- append(output, list(func_args = NA), after = i)
        }
        return(output)
      }
      NULL
    }
  )
)
