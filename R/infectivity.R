#
# infectivity
#
# Copyright (C) 2025 CAMAGLIA Francesco - Institut Pasteur, Paris (France).
#

############################################################################
############################################################################
############################################################################
############################################################################
############################################################################

#' @title
#' Housepi Infectivity Profile R6 Class
#' @description
#' Auxiliary class to deal with infectivity profile.
#'
#' @noRd
HousepiInfecProfile <- R6::R6Class(
  classname = "HousepiInfecProfile",
  cloneable = TRUE,

  ############################################################################
  ############################################################################
  ################################## ACTIVE ##################################
  ############################################################################
  ############################################################################

  active = list(
    #' @description
    #' Find parameters required by infectivity profile variable definition.
    param_require = function(value) {
      if (missing(value)) {
        required_param_set <- c()
        if (is.list(self$defs) && (length(self$defs) > 0)) {
          if (!any(is.na((self$defs[c("arg_value", "arg_type")])))) {
            aux <- data.frame(
              vals = strsplit(self$defs[["arg_value"]], "/")[[1]],
              keys = strsplit(self$defs[["arg_type"]], "/")[[1]]
            )
            aux <- aux[aux[["keys"]] == "P", "vals"]
            required_param_set <- c(required_param_set, aux[!is.na(aux)])
          }
        }
        return(required_param_set)
      }
      stop_h("'param_require' is read only.")
    }
  ),

  ###################
  ##### PUBLIC ######
  ###################

  public = list(
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
    #' @param name (`character(1)`)\cr
    #' Name of the infectivity profile.
    #' @param discrete (`logical(1)`)\cr
    #' IF `TRUE` the profile is discrete.
    #' @param args (named `list()`)\cr
    #' Arguments (eventually with types) for the profile.
    #' @param min (`double(1)` | `character(1)`)\cr
    #' One between:
    #' - Scalar value;
    #' - Parameter name (followed by `@P`).
    #' @param max (`double(1)` | `character(1)`)\cr
    #' One between:
    #' - Scalar value;
    #' - Parameter name (followed by `@P`).
    add = function(name,
                   discrete = FALSE,
                   args = list(),
                   min = NA,
                   max = NA) {
      #----- name -----#
      if (is.na(name)) {
        stop_h("Incomplete definition: missing name.")
      }

      #----- arguments -----#
      arg_name_ <- NA
      arg_value_ <- NA
      arg_type_ <- NA
      # assign method only if arguments are provided
      if (length(args) > 0) {
        arg_name_ <- str_join_(names(args))
        tmp <- str_split_(args)
        arg_value_ <- str_join_(tmp$var)
        arg_type_ <- str_join_(tmp$typ)
      }

      #----- minimum value -----#
      min_ <- NA
      min_type_ <- NA
      if (!is.na(min)) {
        min_ <- str_split_(min)$var
        min_type_ <- str_split_(min)$typ
      }

      #----- maximum value -----#
      max_ <- NA
      max_type_ <- NA
      if (!is.na(max)) {
        max_ <- str_split_(max)$var
        max_type_ <- str_split_(max)$typ
      }

      # finally store definition
      self$defs <- list(
        name = name, discrete = isTRUE(discrete),
        arg_name = arg_name_, arg_value = arg_value_, arg_type = arg_type_,
        min = min_, min_type = min_type_, max = max_, max_type = max_type_
      )

      invisible(self)
    }
  )
)
