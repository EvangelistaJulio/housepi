#      __
#     / /_  ____  __  __ ____ __  ____  __
#    / __ \/ __ \/ / / / ___/ _ \/ __ \/ /
#   / / / / /_/ / /_/ (__  )  __/ /_/ / /
#  /_/ /_/\____/\____/____/\___/ ____/_/
#                             /_/
#
# Copyright (C) 2025 CAMAGLIA Francesco - Institut Pasteur, Paris (France).
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
        for (typed_col in c("func", "max")) {
          for (x in self$defs[[typed_col]]) {
            output <- c(output, walk_encoded_func(x, type))
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
      if (missing(value)) return(private$.find_def_required("D"))
      stop_h("'data_required' is read only.")
    },

    #' @field param_required `character()`\cr
    #' Parameters required by infectivity profile definition.
    param_required = function(value) {
      if (missing(value)) return(private$.find_def_required("P"))
      stop_h("'param_required' is read only.")
    },

    #' @field augm_required `character()`\cr
    #' Augmented variables required by infectivity profile definition.
    augm_required = function(value) {
      if (missing(value)) return(private$.find_def_required("A"))
      stop_h("'augm_required' is read only.")
    }
  ),

  ############################################################################
  ################################## PUBLIC ##################################
  ############################################################################

  public = list(
    #' @field defs `list`\cr
    #' Infectivity profile definition list.
    defs = NULL,

    #' @description
    #' New object of class
    #' [HousepiInfecProfile][housepi::HousepiInfecProfile].
    initialize = function() {
      self$defs <- NULL
      invisible(self)
    },

    #' @description
    #' Add definition of infectivity profile.
    #'
    #' @param name `character(1)`\cr
    #' Name of the infectivity profile.
    #' @param discr `logical(1)`\cr
    #' IF `TRUE` the profile is discrete.
    #' @param args named `list`\cr
    #' Arguments (eventually with types) for the profile.
    #' @param max `double(1)` | `character(1)`\cr
    #' One between:
    #' - Scalar value;
    #' - Parameter name (followed by `@P`).
    add = function(name,
                   discr = FALSE,
                   args = NULL,
                   max = NULL) {
      # check infectivity profile name
      if (is.na(name)) stop_h("Incomplete definition: missing name.")
      self$defs <- list(func = name, discr = isTRUE(discr))

      # load arguments
      if (!is.null(args)) {
        self$defs$func <- DefEncoder$new(as.list(args))$combine_to_call(name)
      }

      # load maximum value
      self$defs$max <- DefEncoder$new(max)$x

      invisible(self)
    }
  )
)
