#
# model_class.R
#
# Copyright (C) 2025 CAMAGLIA Francesco - Institut Pasteur, Paris (France).
#

#' @title Housepi Augmented R6 Class
#' @description
#' Help with augmented variable definitions of model Rcpp object.
#'
#' @noRd
HousepiAugmented <- R6::R6Class(
  classname = "HousepiAugmented",
  cloneable = TRUE,

  ############################################################################
  ############################################################################
  ################################# PRIVATE ##################################
  ############################################################################
  ############################################################################

  private = list(
    # from : augmented variable model argument definitions
    .find_def_require = function(key) {
      required_set <- c()
      if (length(self$defs) > 0) {
        x <- self$defs[c("arg_value", "arg_at")]
        for (i in seq_along(x)) {
          r <- unlist(x[i, ])
          if (!any(is.na(r))) {
            arg2par <- data.frame(
              vals = strsplit(r[1], "/")[[1]],
              keys = strsplit(r[2], "/")[[1]]
            )
            required_set <- c(required_set,
                              arg2par[arg2par$keys == key]$vals)
          }
        }
      }
      return(unique(required_set))
    }
  ),

  ############################################################################
  ############################################################################
  ################################## ACTIVE ##################################
  ############################################################################
  ############################################################################

  active = list(
    #' @description
    #' Find parameters required by augmented variable definitions.
    param_require = function(value) {
      if (missing(value)) {
        return(private$.find_def_require("P"))
      }
      stop_h("'param_require' is read only.")
    },

    #' @description
    #' Find data required by augmented variable definitions.
    data_require = function(value) {
      if (missing(value)) {
        return(private$.find_def_require("D"))
      }
      stop_h("'data_require' is read only.")
    }
  ),

  ############################################################################
  ############################################################################
  ################################## PUBLIC ##################################
  ############################################################################
  ############################################################################

  public = list(
    defs = NULL,

    #' @description
    #' New instance of class
    #' [HousepiAugmented][housepi::HousepiAugmented].
    initialize = function() {
      self$defs <- NULL
      invisible(self)
    },

    #' @description
    #' Find augmented variables required by the model.
    #'
    #' @param refer_file
    #' Path/to/.refer source file.
    require = function(refer_file) {
      # vector of regular expression to use for assertion
      regular_expr <- c(
        # from : .augm(, "string" or augm(, 'string'
        '\\.augm\\s*\\(\\s*[^,]*\\s*,\\s*["\'](.*?)["\']',
        # from : .set_augm_value(, "string"
        '\\.set_augm_value\\s*\\(\\s*[^,]*\\s*,\\s*["\'](.*?)["\']'
      )
      # find required augmented variables from regular expressions
      return(template_find_requirement(refer_file, regular_expr))
    },

    #' @description
    #' Remove a given augmented variable
    #'
    #' @param name `character(1)`\cr
    #' Name of the augmented variable.
    #' @param status `integer(1)`\cr
    #' The associated status of the definition.
    #' Default NA is for all infectious.
    #' @param show_warnings `logical(1)`\cr
    #' Wheter to show warnings or not. Default is `TRUE`.
    remove = function(name,
                      status = NA,
                      show_warnings = TRUE) {
      mask <- (self$defs[["name"]] == name)
      if (any(mask)) {
        if (is.na(status)) {
          # drop previous definition
          self$defs <- self$defs[!mask, ]
          return(TRUE)
        } else {
          # assert if same name and status
          mask <- mask & (self$defs[["status"]] == status)
          mask[is.na(mask)] <- FALSE
          if (any(mask)) {
            # drop previous definition
            self$defs <- self$defs[!mask, ]
            return(TRUE)
          }
        }

      } else {
        if (isTRUE(show_warnings)) {
          warning_h(
            sprintf(
              "Attempt to remove unknown augmented variable \"%s\".", name
            )
          )
        }
      }

      return(FALSE)
    },

    #' @description
    #' Add a new augmented variable definition.
    #'
    #' @param name `character(1)`\cr
    #' Name of the augmented variable.
    #' @param status `integer(1)`\cr
    #' The associated status of the definition.
    #' Default NA is for all infectious.
    #' @param value `double(1)` | `character(1)`\cr
    #' One between:
    #' - Scalar value;
    #' - Parameter name (followed by `@P`);
    #' - Data name (followed by `@D`).
    #' @param init_rand `character(1)`\cr
    #' If `TRUE` draw initial vale from model distribution.
    #' @param constant `logical(1)`\cr
    #' If `FALSE` the augmented variable is treated as a constant.
    #' @param discrete `logical(1)`\cr
    #' IF `TRUE` the augmenter variable is discrete.
    #' @param model `character(1)`\cr
    #' Name of the model distribution.
    #' @param args (named `list()`\cr
    #' Arguments for the model distribution.
    #' @param min `double(1)` | `character(1)`\cr
    #' Minimum value. One between:
    #' - Scalar value;
    #' - Parameter name (followed by `@P`);
    #' - Data name (followed by `@D`).
    #' @param max `double(1)` | `character(1)`\cr
    #' Maxmimum value. One between:
    #' - Scalar value;
    #' - Parameter name (followed by `@P`);
    #' - Data name (followed by `@D`).
    #' @param rule `character(1)`\cr
    #' Name of the inference proposal rule.
    #' @param scale `double(1)`\cr
    #' Scale parameter for the proposal rule.
    #' @param from `double(1)` | `character(1)`\cr
    #' Advanced options for non random walk moves.
    #' `FIXME: missing`
    #' @noRd
    add = function(name = NA,
                   status = NA,
                   value = NA,
                   init_rand = FALSE,
                   constant = NA,
                   discrete = FALSE,
                   model = NA,
                   args = list(),
                   min = NA,
                   max = NA,
                   rule = NA,
                   scale = NA,
                   from = NA) {
      #----- name -----#
      if (is.na(name)) {
        stop_h("Incomplete definition: missing name.")
      }

      # remove previous definition if it exists
      if (self$remove(name, status, FALSE)) {
        if (is.na(status)) {
          stop_h(
            sprintf("Overwriting existing \"%s\" definition.", name)
          )
        } else {
          stop_h(
            sprintf(
              "Overwriting existing \"%s\" definition for status %d.",
              name, status
            )
          )
        }
      }

      # assume it's a constant if model and rule are not provided
      if (is.na(constant)) {
        constant <- is.na(model) && is.na(rule)
      }

      #----- value/constant/init_rand -----#
      value_ <- NA
      value_at_ <- NA
      if (isTRUE(constant)) {
        if (is.na(value)) {
          stop_h("Missing value for constant augmented variable.")
        } else {
          value_ <- str_split_(value)$var
          value_at_ <- str_split_(value)$typ
        }
      } else if (isTRUE(init_rand)) {
        if (is.na(model)) {
          stop_h("Unknown model for random initialization.")
        }
      } else {
        if (is.na(value)) {
          if (is.na(model)) {
            if (is.na(from)) {
              stop_h(
                "Incomplete definition: missing value.",
                "Alternative: provide a model and set 'init_rand = TRUE'."
              )
            } else {
              value_ <- str_split_(from)$var
              value_at_ <- str_split_(from)$typ
            }
          } else {
            init_rand <- TRUE
          }
        } else {
          value_ <- str_split_(value)$var
          value_at_ <- str_split_(value)$typ
        }
      }

      #----- method -----#
      arg_name_ <- NA
      arg_value_ <- NA
      arg_at_ <- NA
      # assign method only if arguments are provided
      if ((!is.na(model)) && (length(args) > 0)) {
        arg_name_ <- str_join_(names(args))
        tmp <- str_split_(args)
        arg_value_ <- str_join_(tmp$var)
        arg_at_ <- str_join_(tmp$typ)
      }

      #----- rule -----#
      scale_ <- NA
      scale_at_ <- NA
      from_ <- NA
      from_at_ <- NA

      if (isFALSE(constant)) {
        if (is.na(rule)) {
          if (is.na(model)) {
            constant <- TRUE
          } else {
            rule <- "gibbs"
          }
        } else if (rule != "gibbs") {
          if (is.na(scale)) {
            scale_ <- 1
            scale_at_ <- "S"
          } else {
            scale_ <- str_split_(scale)$var
            scale_at_ <- str_split_(scale)$typ
          }
          from_ <- str_split_(from)$var
          from_at_ <- str_split_(from)$typ
        }
      }

      #----- minimum value -----#
      min_ <- NA
      min_at_ <- NA
      if (!is.na(min)) {
        min_ <- str_split_(min)$var
        min_at_ <- str_split_(min)$typ
      }

      #----- maximum value -----#
      max_ <- NA
      max_at_ <- NA
      if (!is.na(max)) {
        max_ <- str_split_(max)$var
        max_at_ <- str_split_(max)$typ
      }

      #----- append the new definitions -----#
      self$defs <- rbind(
        self$defs,
        data.frame(
          name = name,
          status = status,
          value = value_,
          value_at = value_at_,
          init_rand = init_rand,
          constant = isTRUE(constant),
          discrete = isTRUE(discrete),
          model = model,
          arg_name = arg_name_,
          arg_value = arg_value_,
          arg_at = arg_at_,
          min = min_,
          min_at = min_at_,
          max = max_,
          max_at = max_at_,
          rule = rule,
          scale = scale_,
          scale_at = scale_at_,
          from = from_,
          from_at = from_at_,
          stringsAsFactors = FALSE
        )
      )
      invisible(self)
    },

    #' @description
    #' Associate augmented variables from a list.
    #'
    #' @param ...
    #' Augmented variable names with respective list of definitions.
    #' @param from (named `list()`\cr
    #' Named list of augmented variable definitions.
    #' @param append `logical(1)`\cr
    #' If `TRUE` appends the list to current definition (default).
    #' If `FALSE` drop all other augmented variable definitions.
    load_list = function(...,
                         from = NULL,
                         append = TRUE) {
      definitions <- list(...)

      if (length(definitions) == 0) {
        if (length(from) == 0) {
          stop_h("No definitions provided.")
        } else if (is.list(from)) {
          definitions <- from
        } else {
          stop_h("Argument 'from' should be a list.")
        }
      } else if (any(names(definitions) == "")) {
        stop_h("Found definition without name.")
      }

      if (isFALSE(append)) {
        self$defs <- NULL
      }

      for (i in seq_along(definitions)){
        this_def <- as.list(definitions[[i]])
        def_name <- names(definitions[i])
        if (length(names(this_def)) == 0) {
          if (length(this_def) == 1) {
            names(this_def) <- "value"
          } else {
            stop_h(
              sprintf(
                "Wrong definition for augmented variable \"%s\".", def_name
              ),
              "Provide a definition list.",
              "Alternative: if constant, a scalar or a typed value."
            )
          }
        }
        this_def[["name"]] <- def_name
        do.call(self$add, this_def)
      }

      invisible(self)
    },

    #' @description
    #' Load augmented variables from template definitions.
    #'
    #' @param templ_file
    #' Path/to/.refer template file.
    #' @param requirements `character()`\cr
    #' Vector containing required parameters.
    #' @param append `logical(1)`\cr
    #' If `TRUE` append to previous definitions (default).
    #' If `FALSE` drop previous definitions.
    load_template = function(templ_file,
                             requirements = NULL,
                             append = TRUE) {
      # check if template file exists
      if (file.exists((templ_file))) {
        # template augmented variable definitions
        templ_defs <- dget(templ_file)

        if (isTRUE(append)) {
          # check set of required augmented variables
          if (length(requirements) == 0) {
            warning_h(
              "Appending template augmented variable definitions failed.",
              "Alternatively use option 'append = FALSE'."
            )
          } else {
            missing_set <- template_find_missing(requirements,
                                                 self$defs[["name"]])

            if (length(missing_set) == 0) {
              warning_h(
                "No template augmented variable to append.",
                "Definitions are complete."
              )
            } else if (!all(missing_set %in% names(templ_defs))) {
              stop_h(
                "No template found for the missing augmented variables."
              )
            } else {
              self$load_list(from = templ_defs[missing_set], append = TRUE)
            }
          }
        } else {
          self$load_list(from = templ_defs, append = FALSE)
        }

      } else {
        stop_h(
          "Impossible to automatically initialize augmented variables.",
          "A file with template definitions was not found."
        )
      }

      invisible(self)
    },

    #' @description
    #' Assert augmented definitions from the model implementation script.
    #'
    #' @param requirements `character()`\cr
    #' Vector containing required parameters.
    #'
    #' @return `logical(1)`\cr
    assert = function(requirements) {
      # logical state for complete definitions
      is_def_complete <- TRUE

      if (length(requirements) == 0) {
        # if requirement set is NULL it means nothing can be said
        is_def_complete <- NA

      } else {
        # find missing augmented variables
        missing_set <- template_find_missing(requirements,
                                             self$defs[["name"]])

        if (length(missing_set) > 0) {
          is_def_complete <- FALSE
          warning_h(
            "Missing augmented variable(s):",
            character_to_string_list(missing_set)
          )
        }

        # point out augmented variables that are not used
        unused_set <- template_find_unused(requirements,
                                           self$defs[["name"]])

        # drop unused augmented variables
        if (length(unused_set) > 0) {
          warning_h(
            "Dropping not used augmented variable(s):",
            character_to_string_list(unused_set)
          )

          for (name in unused_set) {
            self$remove(name)
          }
        }
      }

      return(is_def_complete)
    }
  )
)
