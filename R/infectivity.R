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
    #' @field family `character(1)` with the module family.
    family = "infectivity profile",
    #' @field store_warn_over `character()` with cached overwritten names.
    store_warn_over = character(),
    #' @field defs `data.frame` with formatted infectivity profile definitions.
    defs = NULL,

    #' @description
    #' New object of class
    #' [HousepiInfecProfile][housepi::HousepiInfecProfile].
    initialize = function() {
      self$defs <- NULL
      invisible(self)
    },

    #' @description
    #' Remove a given infectivity profile definition.
    #'
    #' @param status
    #' An `integer(1)` with the associated status of the definition.
    #' @param warn_removal
    #' `logical(1)`to show warnings for failed removal.
    remove = function(status, warn_removal = FALSE) {
      mask <- sapply(self$defs[["status"]], function(x) {
        identical(x, as.integer(status))
      })
      if (any(mask)) {
        self$defs <- self$defs[!mask, ]
        return(status)
      } else if (isTRUE(warn_removal)) {
        warning_h(sprintf(
          "Attempt to remove unknown infectivity profile for status '%s'.",
          status
        ))
      }
      c()
    },

    #' @inherit add_infec_profile
    add_quo = function(func, status = NA_integer_,
                       discr = FALSE, ..., max = NULL, .warn_over = TRUE) {
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

      # drop previous definition eventually
      rmv_name <- self$remove(status)
      if (.warn_over) template_overwrite_warn(rmv_name, self$family)
      else self$store_warn_over <- c(self$store_warn_over, rmv_name)

      self$defs <- rbind(self$defs, data.frame(list(
        func = ifelse(is.null(func_args), func,
                      housepi_combine_to_call(func, func_args)),
        status = as.integer(status),
        min = ifelse(isTRUE(discr), 1, 0),
        max = ifelse(is.null(max), NA_character_, max),
        discr = isTRUE(discr)
      )))

      invisible(self)
    },

    #' @inherit load_infec_profile
    load_quosures = function(xquos, append = FALSE) {
      if (isFALSE(append)) self$initialize()
      housepi_load_defin_quos(self, "func", xquos)
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
