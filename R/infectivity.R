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
    # according to the "at" key of the definition (A, D, P).
    .find_def_required = function(at) {
      # tuples (value, type) amongst definitions
      tuple_at_cols <- list(c("arg_val", "arg_at"),
                            c("max", "max_at"))

      required <- c()
      if (length(self$defs) > 0) {
        for (tuple in tuple_at_cols) {
          xdef <- self$defs[tuple]
          rw <- unlist(xdef)
          if (!any(is.na(rw))) {
            arg2par <- data.frame(
              val = strsplit(rw[1], "/")[[1]],
              at = strsplit(rw[2], "/")[[1]]
            )
            required <- c(required, arg2par[arg2par$at == at, "val"])
          }
        }
      }
      # remove duplicates and add source name
      set_source_name_(unique(required), "infectivity profile")
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
                   max = NA) {
      # check infectivity profile name
      if (is.na(name)) stop_h(
        "Incomplete definition: missing name."
      )
      self$defs <- list(name = name, discr = isTRUE(discr))

      # load arguments
      parse_args <- data.frame(t(sapply(args, extract_func_)))
      if (any(.is_func(parse_args$func))) stop_h(
        "Calls in the argument definitions of infectivity profiles
        are not supported."
      )
      self$defs <- modifyList(
        self$defs,
        list(arg_key = type_join_(row.names(parse_args)),
             arg_val = type_join_(parse_args$val),
             arg_at = type_join_(parse_args$at))
      )

      # load maximum value
      parse_max <- data.frame(t(sapply(args, extract_func_)))
      if (.is_func(parse_max$func)) stop_h(
        "Calls in the maximum definition of infectivity profiles
        are not supported."
      )
      self$defs <- modifyList(
        self$defs,
        list(max = parse_max$val, max_at = parse_max$at)
      )

      invisible(self)
    }
  )
)
