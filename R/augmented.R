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
    # handle existing definition
    .handle_existing_definition = function(name, status) {
      # check for missing name
      if (is.na(name)) {
        stop_h("Incomplete definition: missing name.")
      }

      # remove previous definition # FIXME: why stop?
      if (self$remove(name, status, FALSE)) {
        warning_h(sprintf(
          "Overwriting existing \"%s\" definition%s",
          name, ifelse(is.na(status), ".", sprintf(" for status %d.", status))
        ))
      }
      invisible(self)
    },

    # parse augmented implicit function
    .handle_augmented_implicit_func = function(opts, value) {
      parsed <- extract_func_(value)
      if (!is.na(parsed$func)) {
        if (is.na(opts$func)) {
          opts$func <- parsed$func
          opts$farg_key <- parsed$key
          opts$farg_val <- parsed$val
          opts$farg_at  <- parsed$at
        } else {
          stop_h(
            "Augmented function multiple definitions.",
            "Use either 'value' or 'func' argument."
          )
        }
      }
      return(opts)
    },

    # parse augmented implicit function
    .handle_augmented_explicit_func = function(opts, func_args) {
      parsed <- extract_func_(opts$func)

      if (is.null(func_args)) {
        if (!is.na(parsed$func)) {
          opts$func <- parsed$func
          opts$farg_key <- parsed$key
          opts$farg_val <- parsed$val
          opts$farg_at  <- parsed$at
        }

      } else {
        if (is.na(parsed$func)) {
          if (is.na(opts$func)) {
            stop_h("Augmented function not defined.")
          } else {
            opts$farg_key <- type_join_(names(func_args))
            opts$farg_val <- type_join_(type_split_(func_args)$val)
            opts$farg_at  <- type_join_(type_split_(func_args)$at)
          }
        } else {
          stop_h(
            "Augmented function defined both implicitly and explicitly.",
            "Use either implicit 'func' or 'func_args' argument."
          )
        }
      }

      return(opts)
    },

    # process model arguments
    .handle_model_arguments = function(opts, args) {
      if (is.na(opts$model))  {
        message <- sprintf(
          "Spare arguments for augmented variable \"%s\".", opts$name
        )
        if ((!is.na(opts$func)) && is.na(opts$farg_key)) {
          message <- c(message, "Did you mean 'func_args'?")
        }
        stop_h(message)

      } else {
        opts$arg_key <- type_join_(names(args))
        opts$arg_val <- type_join_(type_split_(args)$val)
        opts$arg_at  <- type_join_(type_split_(args)$at)
      }

      return(opts)
    },

    # process constant or random
    .handle_constant_or_random = function(opts, value, from) {
      if (isTRUE(opts$const)) {
        if (is.na(value)) {
          stop_h("Missing value for constant augmented variable.")
        }
        opts$value <- type_split_(value)$val
        opts$value_at <- type_split_(value)$at
      } else if (isTRUE(opts$init_rand)) {
        if (is.na(opts$model)) {
          stop_h("Unknown model for random initialization.")
        }
      } else {
        if (is.na(value)) {
          if (is.na(opts$model)) {
            if (is.na(from)) {
              stop_h("Incomplete definition: missing value.",
                     "Alternative: provide a model and set 'init_rand = TRUE'.")
            } else {
              opts$value <- type_split_(from)$val
              opts$value_at <- type_split_(from)$at
            }
          } else {
            opts$init_rand <- TRUE
          }
        } else {
          opts$value <- type_split_(value)$val
          opts$value_at <- type_split_(value)$at
        }
      }
      return(opts)
    },

    # interpret user choice for proposer rule
    .handle_proposer_rule = function(opts, scale, from) {
      if (isFALSE(opts$const)) {
        if (is.na(opts$rule)) {
          if (is.na(opts$model)) {
            opts$const <- TRUE
          } else {
            opts$rule <- "draw"
          }
        } else if (opts$rule != "draw") {
          if (!is.na(scale)) {
            opts$scale <- type_split_(scale)$val
            opts$scale_at <- type_split_(scale)$at
          }
          opts$from <- type_split_(from)$val
          opts$from_at <- type_split_(from)$at
        }
      }
      return(opts)
    },

    # interpret user choice for bounds
    .handle_bounds = function(opts, min, max) {
      opts$min <- type_split_(min)$val
      opts$min_at <- type_split_(min)$at
      opts$max <- type_split_(max)$val
      opts$max_at <- type_split_(max)$at
      return(opts)
    },

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
                                arg2par[arg2par$keys == key, "vals"])
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
    #' @field defs `data.frame`\cr
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
        # from : .augm("string"
        '.augm\\s*\\(\\s*["\'](.*?)["\']',
        # from : .get_augm_variable("string"
        '.get_augm_variable\\s*\\(\\s*["\'](.*?)["\']',
        # from : .set_augm_value("string"
        '.set_augm_value\\s*\\(\\s*["\'](.*?)["\']'
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
    #' @param x `character(1)`\cr
    #' Name of the augmented variable to remove.
    #' @param status `integer(1)`\cr
    #' The associated status of the definition.
    #' Default NA is for all infectious.
    #' @param show_warnings `logical(1)`\cr
    #' Wheter to show warnings or not. Default is `TRUE`.
    #'
    #' @return `logical(1)`
    remove = function(x, status = NA, show_warnings = TRUE) {
      mask <- (self$defs[["name"]] == x)
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

          } else if (isTRUE(show_warnings)) {
            warning_h(sprintf(
              "Unknown augmented variable \"%s\" for status %d to remove.",
              x, status
            ))
          }
        }

      } else if (isTRUE(show_warnings)) {
        warning_h(sprintf("Unknown augmented variable \"%s\" to remove.", x))
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
    #' @param const `logical(1)`\cr
    #' If `FALSE` the augmented variable is treated as a constant.
    #' @param discr `logical(1)`\cr
    #' IF `TRUE` the augmenter variable is discrete.
    #' @param func `character()`\cr
    #' Name of the augmented function.
    #' @param func_args (named `list()`\cr
    #' Arguments of the augmented function.
    #' @param value `double(1)` | `character(1)`\cr
    #' One between:
    #' - Scalar value;
    #' - Parameter name (followed by `@P`);
    #' - Data name (followed by `@D`).
    #' @param init_rand `character(1)`\cr
    #' If `TRUE` draw initial vale from model distribution.
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
                   const = NA,
                   discr = FALSE,
                   func = NA,
                   func_args = NULL,
                   value = NA,
                   init_rand = FALSE,
                   model = NA,
                   args = NULL,
                   min = NA,
                   max = NA,
                   rule = NA,
                   scale = NA,
                   from = NA) {
      # hande (name, status) and check for existing definitions
      private$.handle_existing_definition(name, status)

      # auxiliary typed options for Rcpp build method
      opts <- list(
        name = as.character(name),
        status = as.integer(status),
        const = isTRUE(const),
        discr = isTRUE(discr),
        func = as.character(func),
        farg_key = NA, farg_val = NA, farg_at = NA,
        value = NA, value_at = NA,
        init_rand = isTRUE(init_rand),
        model = as.character(model),
        arg_key = NA, arg_val = NA, arg_at = NA,
        min = NA, min_at = NA,
        max = NA, max_at = NA,
        rule = as.character(rule),
        scale = 1, scale_at = NA,
        from = NA, from_at = NA
      )

      if (is.na(opts$func)) {
        opts <- private$.handle_augmented_implicit_func(opts, value)
        if (!is.na(opts$func)) {
          opts$init_rand <- opts$const <- opts$discr <- FALSE
        } else if (is.na(opts$const)) {
          opts$const <- is.na(opts$model) && is.na(opts$rule)
        }
      } else {
        opts <- private$.handle_augmented_explicit_func(opts, func_args)
        opts$init_rand <- opts$const <- opts$discr <- FALSE
      }

      if (length(args) > 0) {
        opts <- private$.handle_model_arguments(opts, args)
      }

      if (is.na(opts$func)) {
        opts <- private$.handle_constant_or_random(opts, value, from)
        opts <- private$.handle_proposer_rule(opts, scale, from)
      }

      opts <- private$.handle_bounds(opts, min, max)
      self$defs <- rbind(self$defs, data.frame(opts))

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
          sprintf("Augmented variable%s declared but not used: %s",
                  ifelse(length(unused_set) > 1, "s", ""),
                  character_to_string_list(unused_set)),
          find_closest_match(unused_set, requirements),
          "Unused augmented variables should be dropped."
        )
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
