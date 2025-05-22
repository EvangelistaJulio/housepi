#      __
#     / /_  ____  __  __ ____ __  ____  __
#    / __ \/ __ \/ / / / ___/ _ \/ __ \/ /
#   / / / / /_/ / /_/ (__  )  __/ /_/ / /
#  /_/ /_/\____/\____/____/\___/ ____/_/
#                             /_/
#
# Copyright (C) 2025 CAMAGLIA Francesco - Institut Pasteur, Paris (France).
#
# augmented.R
#

#' @title Housepi Augmented R6 Class
#' @description
#' Help with augmented variable definitions of model Rcpp object.
#'
#' @keywords internal
HousepiAugmented <- R6::R6Class(
  classname = "HousepiAugmented",
  cloneable = TRUE,

  ############################################################################
  ############################################################################
  ################################# PRIVATE ##################################
  ############################################################################
  ############################################################################

  private = list(
    # from : augmented variable definitions
    .find_def_required = function(key) {

      # tuples (value, type) amongst definitions
      tuple_at_cols <- list(c("farg_val", "farg_at"),
                            c("value", "value_at"), c("arg_val", "arg_at"),
                            c("min", "min_at"), c("max", "max_at"),
                            c("scale", "scale_at"), c("from", "from_at"))

      required_set <- c()
      if (length(self$defs) > 0) {
        for (tuple in tuple_at_cols) {
          xdef <- self$defs[tuple]
          for (i in seq_len(dim(xdef)[[1]])) {
            rw <- unlist(xdef[i, ])
            if (!any(is.na(rw))) {
              arg2par <- data.frame(
                vals = strsplit(rw[1], "/")[[1]],
                keys = strsplit(rw[2], "/")[[1]]
              )
              required_set <- c(required_set,
                                arg2par[arg2par$keys == key]$vals)
            }
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
    #' @field param_required `character()`\cr
    #' Find parameters required by augmented variable definitions.
    param_required = function(value) {
      if (missing(value)) return(private$.find_def_required("P"))
      stop_h("'param_required' is read only.")
    },

    #' @field augm_required `character()`\cr
    #' Find augmented required by augmented variable definitions.
    augm_required = function(value) {
      if (missing(value)) return(private$.find_def_required("A"))
      stop_h("'augm_required' is read only.")
    },

    #' @field data_required `character()`\cr
    #' Find data required by augmented variable definitions.
    data_required = function(value) {
      if (missing(value)) return(private$.find_def_required("D"))
      stop_h("'data_required' is read only.")
    }
  ),

  ############################################################################
  ############################################################################
  ################################## PUBLIC ##################################
  ############################################################################
  ############################################################################

  public = list(
    #' @field augm_defs `data.frame`\cr
    #' Augmented variable formatted definitions.
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
      return(c(
        template_find_requirement(refer_file, regular_expr),
        self$augm_required
      ))
    },

    #' @description
    #' Remove a given augmented variable
    #'
    #' @param this `character(1)`\cr
    #' Name of the augmented variable to remove.
    #' @param status `integer(1)`\cr
    #' The associated status of the definition.
    #' Default NA is for all infectious.
    #' @param show_warnings `logical(1)`\cr
    #' Wheter to show warnings or not. Default is `TRUE`.
    remove = function(this,
                      status = NA,
                      show_warnings = TRUE) {
      mask <- (self$defs[["name"]] == this)
      if (any(mask)) {
        if (is.na(status)) {
          # drop previous definition
          self$defs <- self$defs[!mask, ]
          return(TRUE)
        } else {
          # check if same name and status
          mask <- mask & (self$defs[["status"]] == status)
          mask[is.na(mask)] <- FALSE
          if (any(mask)) {
            # drop previous definition
            self$defs <- self$defs[!mask, ]
            return(TRUE)
          }
        }

      } else if (isTRUE(show_warnings)) {
        warning_h(
          sprintf("Unknown augmented variable \"%s\" to remove.", this)
        )
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
    #' @param func `character()`\cr
    #' Name of the augmented function.
    #' @param value `double(1)` | `character(1)`\cr
    #' One between:
    #' - Scalar value;
    #' - Parameter name (followed by `@P`);
    #' - Data name (followed by `@D`).
    #' @param init_rand `character(1)`\cr
    #' If `TRUE` draw initial vale from model distribution.
    #' @param const `logical(1)`\cr
    #' If `FALSE` the augmented variable is treated as a constant.
    #' @param discr `logical(1)`\cr
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
    #' 
    #' @keywords internal
    add = function(name = NA,
                   status = NA,
                   func = NA,
                   value = NA,
                   init_rand = FALSE,
                   const = NA,
                   discr = FALSE,
                   model = NA,
                   args = NULL,
                   min = NA,
                   max = NA,
                   rule = NA,
                   scale = NA,
                   from = NA) {
      # auxiliary typed options for Rcpp build method
      .opts <- list(
        farg_key = NA, farg_val = NA, farg_at = NA,
        value = NA, value_at = NA, 
        arg_key = NA, arg_val = NA, arg_at = NA,
        min = NA, min_at = NA,
        max = NA, max_at = NA,
        scale = NA, scale_at = NA,
        from = NA, from_at = NA
      )

      # check for missing name
      if (is.na(name)) {
        stop_h("Incomplete definition: missing name.")
      }

      # remove previous definition if it exists
      if (self$remove(name, status, FALSE)) {
        msg <- if (is.na(status)) 
          sprintf("Overwriting existing \"%s\" definition.", name) 
        else 
          sprintf("Overwriting existing \"%s\" definition for status %d.", name, status)
        stop_h(msg)
      }

      # check if a function is provided
      if (is.na(func)) {
        # check for augmented function string definition
        afunc_str_def <- extract_augm_func_(value)
        
        if (!is.null(afunc_str_def)) {
          func <- afunc_str_def$func
          .opts$farg_key <- afunc_str_def$key
          .opts$farg_val <- afunc_str_def$val
          .opts$farg_at <- afunc_str_def$at
          init_rand <- FALSE
          const <- FALSE
          discr <- FALSE

        } else if (is.na(const)) {
          # assume constant if func or model/rule not provided
          const <- is.na(model) && is.na(rule)
        }

      } else {
        # FIXME // UPDATE
        #const <- FALSE     
        stop_h(
          "Augmented function is not supported yet.",
        )
      }

      if (is.na(func)) {
        #----- value/const/init_rand -----#
        if (isTRUE(const)) {
          if (is.na(value)) {
            stop_h("Missing value for constant augmented variable.")
          } else {
            .opts$value <- type_split_(value)$val
            .opts$value_at <- type_split_(value)$at
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
                .opts$value <- type_split_(from)$val
                .opts$value_at <- type_split_(from)$at
              }
            } else {
              init_rand <- TRUE
            }
          } else {
            .opts$value <- type_split_(value)$val
            .opts$value_at <- type_split_(value)$at
          }
        }

        #----- method -----#
        # assign method only if arguments are provided
        if ((!is.na(model)) && (length(args) > 0)) {
          .opts$arg_key <- type_join_(names(args))
          .opts$arg_val <- type_join_(type_split_(args)$val)
          .opts$arg_at <- type_join_(type_split_(args)$at)
        }

        #----- rule -----#

        if (isFALSE(const)) {
          if (is.na(rule)) {
            if (is.na(model)) {
              const <- TRUE
            } else {
              rule <- "draw"
            }
          } else if (rule != "draw") {
            if (is.na(scale)) {
              .opts$scale <- 1
              .opts$scale_at <- "S"
            } else {
              .opts$scale <- type_split_(scale)$val
              .opts$scale_at <- type_split_(scale)$at
            }
            .opts$from <- type_split_(from)$val
            .opts$from_at <- type_split_(from)$at
          }
        }

        #----- minimum value -----#
        if (!is.na(min)) {
          .opts$min <- type_split_(min)$val
          .opts$min_at <- type_split_(min)$at
        }

        #----- maximum value -----#
        if (!is.na(max)) {
          .opts$max <- type_split_(max)$val
          .opts$max_at <- type_split_(max)$at
        }
      }

      #----- append the new definitions -----#
      self$defs <- rbind(
        self$defs,
        data.frame(
          name = name, status = status,
          func = func, 
          farg_key = .opts$farg_key, farg_val = .opts$farg_val, farg_at = .opts$farg_at,
          value = .opts$value, value_at = .opts$value_at, init_rand = init_rand,
          const = isTRUE(const), discr = isTRUE(discr),
          model = model,
          arg_key = .opts$arg_key, arg_val = .opts$arg_val, arg_at = .opts$arg_at,
          min = .opts$min, min_at = .opts$min_at, max = .opts$max, max_at = .opts$max_at,
          rule = rule, scale = .opts$scale, scale_at = .opts$scale_at,
          from = .opts$from, from_at = .opts$from_at,
          stringsAsFactors = FALSE
        )
      )

      invisible(self)
    },

    #' @description
    #' Load augmented variable definitions either as direct arguments
    #' or from predefined list.
    #'
    #' @param ...
    #' Named arguments with respective list of definition details.
    #' @param from (named `list()`\cr
    #' Named list of parameter definitions.
    #' @param append `logical(1)`\cr
    #' If `TRUE` appends to current definitions.
    #' If `FALSE` drop all pre-existing definitions (default).
    load_list = function(...,
                         from = NULL,
                         append = FALSE) {
      definitions <- load_definition_list(..., from = from)

      # remove previos definitions if not append
      if (isFALSE(append)) self$defs <- NULL

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
        this_def$name <- def_name
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
    #' If `TRUE` append missing augmented variable definitions.
    #' If `FALSE` drop all pre-existing definitions (default).
    load_template = function(templ_file,
                             requirements = NULL,
                             append = FALSE) {
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
    assert = function(requirements) {
      defined <- self$defs[["name"]]

      # point out augmented variables that are not used
      unused_set <- template_find_unused(requirements, defined)
      # drop unused augmented variables

      if (length(unused_set) > 0) {
        warning_h(
          sprintf("Dropping unused augmented variable%s %s",
                  ifelse(length(unused_set) > 1, "s:", ""),
                  character_to_string_list(unused_set)),
          find_closest_match(unused_set, requirements)
        )
        for (name in unused_set) {
          self$remove(name)
        }
      }

      # find missing augmented variables
      missing_set <- template_find_missing(requirements, defined)
      if (length(missing_set) > 0) {
        stop_h(
          sprintf("Missing augmented variable definition%s for:",
                  ifelse(length(missing_set) > 1, "s", "")),
          character_to_string_list(missing_set)
        )
      }

      invisible(self)
    }
  )
)
