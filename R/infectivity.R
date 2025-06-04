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
#' Housepi Infectivity Profile R6 Class
#' @description
#' Auxiliary class to deal with infectivity profile.
#'
#' @keywords internal
HousepiInfecProfile <- R6::R6Class(
  classname = "HousepiInfecProfile",
  cloneable = TRUE,

  ############################################################################
  ################################## ACTIVE ##################################
  ############################################################################

  active = list(
    #' @field param_required `character()`\cr
    #' Parameters required by infectivity profile definition.
    param_required = function(value) {
      if (missing(value)) {
        requirements <- c()
        if (is.list(self$defs) && (length(self$defs) > 0)) {
          if (!any(is.na((self$defs[c("arg_val", "arg_at")])))) {
            aux <- data.frame(
              vals = strsplit(self$defs[["arg_val"]], "/")[[1]],
              keys = strsplit(self$defs[["arg_at"]], "/")[[1]]
            )
            aux <- aux[aux[["keys"]] == "P", "vals"]
            requirements <- c(requirements, aux[!is.na(aux)])
          }
        }
        return(requirements)
      }
      stop_h("'param_required' is read only.")
    },

    #' @field augm_required `character()`\cr
    #' Augmented variables required by infectivity profile definition.
    augm_required = function(value) {
      if (missing(value)) {
        requirements <- c()
        if (is.list(self$defs) && (length(self$defs) > 0)) {
          if (!any(is.na((self$defs[c("arg_val", "arg_at")])))) {
            aux <- data.frame(
              vals = strsplit(self$defs[["arg_val"]], "/")[[1]],
              keys = strsplit(self$defs[["arg_at"]], "/")[[1]]
            )
            aux <- aux[aux[["keys"]] == "A", "vals"]
            requirements <- c(requirements, aux[!is.na(aux)])
          }
        }
        return(requirements)
      }
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
    #' New instance of class
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
    #' @param args (named `list()`\cr
    #' Arguments (eventually with types) for the profile.
    #' @param min `double(1)` | `character(1)`\cr
    #' One between:
    #' - Scalar value;
    #' - Parameter name (followed by `@P`).
    #' @param max `double(1)` | `character(1)`\cr
    #' One between:
    #' - Scalar value;
    #' - Parameter name (followed by `@P`).
    add = function(name,
                   discr = FALSE,
                   args = list(),
                   min = NA,
                   max = NA) {
      #----- name -----#
      if (is.na(name)) {
        stop_h("Incomplete definition: missing name.")
      }

      #----- arguments -----#
      arg_key_ <- NA
      arg_val_ <- NA
      arg_at_ <- NA
      # assign method only if arguments are provided
      if (length(args) > 0) {
        arg_key_ <- type_join_(names(args))
        tmp <- type_split_(args)
        arg_val_ <- type_join_(tmp$val)
        arg_at_ <- type_join_(tmp$at)
      }

      #----- minimum value -----#
      min_ <- NA
      min_at_ <- NA
      if (!is.na(min)) {
        min_ <- type_split_(min)$val
        min_at_ <- type_split_(min)$at
      }

      #----- maximum value -----#
      max_ <- NA
      max_at_ <- NA
      if (!is.na(max)) {
        max_ <- type_split_(max)$val
        max_at_ <- type_split_(max)$at
      }

      # finally store definition
      self$defs <- list(
        name = name, discr = isTRUE(discr),
        arg_key = arg_key_, arg_val = arg_val_, arg_at = arg_at_,
        min = min_, min_at = min_at_, max = max_, max_at = max_at_
      )

      invisible(self)
    }
  )
)
