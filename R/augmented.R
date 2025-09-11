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

extract_at_values <- function(defs, at) {
  # tuples (value, type) amongst augmented definitions
  at_tuples <- list(c("farg_val", "farg_at"),
                    c("value", "value_at"), c("arg_val", "arg_at"),
                    c("min", "min_at"), c("max", "max_at"),
                    c("from", "from_at"))

  output <- c()
  if (length(defs) > 0) {
    for (tuple in at_tuples) {
      xdef <- defs[tuple]
      for (i in seq_len(dim(xdef)[[1]])) {
        rw <- unlist(xdef[i, ])
        if (!any(is.na(rw))) {
          arg2par <- data.frame(val = strsplit(rw[1], "/")[[1]],
                                at = strsplit(rw[2], "/")[[1]])
          output <- c(output, arg2par[arg2par$at == at, "val"])
        }
      }
    }
  }
  unique(output)
}

#' @title Housepi Augmented
#' @description
#' Help with augmented variable definitions of model Rcpp object.
#'
#' @keywords internal
#' @noRd
HousepiAugmented <- R6::R6Class(
  classname = "HousepiAugmented",
  cloneable = TRUE,

  ############################################################################
  ################################# PRIVATE ##################################
  ############################################################################

  private = list(
    deep_clone = function(name, value) helper_clone(value),

    # handle existing definition
    .helper_existing_definition = function(name, status) {
      # check for missing name
      if (is.na(name)) stop_h("Incomplete definition: missing name.")
      # wrong status
      else if (!is.na(status) && !is.numeric(status)) stop_h(
        "Wrong status: must be NA or a positive integer alue."
      )
      else if (!is.na(status) && (status <= 0)) stop_h(
        "Wrong status: must be NA or a positive integer value."
      )
      # remove previous definition
      else if (self$remove(name, status, FALSE)) warning_h(sprintf(
        "Overwriting existing \"%s\" definition%s",
        name, ifelse(is.na(status), ".", sprintf(" for status %d.", status))
      ))
      invisible(self)
    },

    # parse augmented implicit function
    .helper_implicit_func = function(opts, expr) {
      parsed <- extract_func_(expr)
      if (!is.na(parsed$func)) {
        if (is.na(opts$func)) opts <- modifyList(
          opts,
          list(func = parsed$func, farg_key = parsed$key,
               farg_val = parsed$val, farg_at  = parsed$at)
        )
        else stop_h(
          sprintf("Mutiple function definitions for \"%s\".", opts$name),
          "Use either 'expr' or 'func' argument."
        )
      }

      if (!is.na(opts$func)) {
        opts$init_rand <- opts$const <- opts$discr <- FALSE
      } else if (is.na(opts$const)) {
        opts$const <- is.na(opts$model) && is.na(opts$rule)
      }

      opts
    },

    # process function arguments
    .helper_function_arguments = function(opts, func_args) {

      if (!is.null(func_args)) {
        parsed <- data.frame(t(sapply(func_args, extract_func_)))
        if (any(.is_func(parsed$func))) stop_h(sprintf(
          "Calls in the function arguments definition of augmented \"%s\" 
          are not supported.", opts$name
        ))

        opts <- modifyList(opts, list(farg_key = type_join_(row.names(parsed)),
                                      farg_val = type_join_(parsed$val),
                                      farg_at = type_join_(parsed$at)))
      }
      opts
    },

    # process model arguments
    .helper_model_arguments = function(opts, args) {
      if (is.na(opts$model))  {
        msg <- sprintf("Spare arguments for augmented \"%s\".", opts$name)
        if ((!is.na(opts$func)) && is.na(opts$farg_key)) {
          msg <- c(msg, "Did you mean 'func_args'?")
        }
        stop_h(msg)

      } else {
        parsed <- data.frame(t(sapply(args, extract_func_)))
        if (any(.is_func(parsed$func))) stop_h(sprintf(
          "Calls in the model arguments definition of augmented \"%s\" 
          are not supported.", opts$name
        ))
        opts <- modifyList(opts, list(arg_key = type_join_(row.names(parsed)),
                                      arg_val = type_join_(parsed$val),
                                      arg_at = type_join_(parsed$at)))
      }

      opts
    },

    # process constant or random
    .helper_constant_or_random = function(opts, value, from) {
      parsed <- extract_func_(value)

      if (.is_func(parsed$func)) stop_h(sprintf(
        "Calls in the initial value definition of augmented \"%s\" 
        are not supported.", opts$name
      ))

      if (isTRUE(opts$const)) {
        if (is.na(parsed$val)) stop_h(sprintf(
          "Missing initial value for constant augmented \"%s\".", opts$name
        ))
        opts <- modifyList(opts, list(value = parsed$val, value_at = parsed$at))

      } else if (isTRUE(opts$init_rand)) {
        if (is.na(opts$model)) stop_h(sprintf(
          "Unknown model for random initialization of augmented \"%s\".",
          opts$name
        ))

      } else {
        if (is.na(parsed$val)) {
          if (is.na(opts$model)) {
            parse_from <- extract_func_(from)
            if (.is_func(parse_from$func)) stop_h(sprintf(
              "Calls in the from value definition of augmented \"%s\" 
              are not supported.", opts$name
            ))
            else if (!is.na(parse_from$val)) opts <- modifyList(
              opts, list(value = parse_from$val, value_at = parse_from$at)
            )
            # else: constant NA augmented variable
          } else {
            opts$init_rand <- TRUE
          }
        } else {
          opts <- modifyList(opts,
                             list(value = parsed$val, value_at = parsed$at))
        }
      }

      opts
    },

    # interpret user choice for proposer rule
    .helper_proposer_rule = function(opts, scale, from) {
      if (isFALSE(opts$const)) {
        if (is.na(opts$rule)) {
          if (is.na(opts$model)) opts$const <- TRUE
          else opts$rule <- "sample"
        } else if (opts$rule != "sample") {
          if (!is.na(scale)) {
            if (is.numeric(scale)) opts$scale <- scale
            else stop_h("Rule scale should be NA or a numeric value.")
          }

          # extract proposer rule 'from' value
          parsed <- extract_func_(from)
          if (.is_func(parsed$func)) stop_h(sprintf(
            "Calls in the from value definition of augmented \"%s\" 
              are not supported.", opts$name
          ))
          opts <- modifyList(opts, list(from = parsed$val, from_at = parsed$at))
        }
      }

      opts
    },

    # interpret user choice for bounds
    .helper_bounds = function(opts, min, max) {
      parse_min <- extract_func_(min)
      parse_max <- extract_func_(max)
      # UPDATE: develop smarter bounds
      if (.is_func(parse_min$func)) stop_h(sprintf(
        "Call in the minimum definition of augmented \"%s\" is not supported.",
        opts$name
      ))
      if (.is_func(parse_max$func)) stop_h(sprintf(
        "Call in the maximum definition of augmented \"%s\" is not supported.",
        opts$name
      ))

      modifyList(opts, list(min = parse_min$val, min_at = parse_min$at,
                            max = parse_max$val, max_at = parse_max$at))
    },

    # Definitions required by augmented variable definitions
    # according to the "at" key of the definition (A, D, P).
    .find_def_required = function(at) {
      required <- extract_at_values(self$defs, at)
      # add source name
      set_source_name_(required, "augmented")
    },

    # recursive explore of the dependencies tree to point out errors
    .explore_dependencies = function(origin, root_def, status_set) {
      root_name <- unique(root_def$name) # name of the closest root definition
      root_status <- root_def$status     # status of the closest root definition
      is_generic <- is.na(root_status)   # whether generic or status-specific

      if (any(is_generic) && length(root_status) > 1) {
        # assert for unique definition
        stop_h(
          sprintf("Equivocal definitions for augmented \"%s\"", root_name),
          ifelse(is.all(is_generic),
                 "Provide a single generic definitions.",
                 "Provide a generic or a status-specific definition.")
        )

      } else if (length(root_status) == 1) {
        # check all depedendecies
        for (r in extract_at_values(root_def, "A")) {
          # check for nested definitions
          if (r == origin) {
            stop_h(sprintf(
              "Found nested definitions for \"%s\" and \"%s\". Please correct.",
              origin, r
            ))
          }

          branch_def <- subset(self$defs, name == r)
          branch_status <- branch_def$status   # status of the branch definition

          # check if at a status-specific dead end
          is_dead_end <- FALSE
          if (length(branch_status) == 1) if (!is.na(branch_status)) {
            # if branch is generic, than it's not a dead end
            is_dead_end <- !identical(branch_status, root_status)
          }

          if (is_dead_end) {
            stop_h(sprintf(
              "Missing definition for augmented \"%s\"%s.",
              r, ifelse(is_generic, "", sprintf(" status %s", root_status))
            ))
          } else {
            if (!is_generic) status_set <- root_status
            private$.explore_dependencies(origin, branch_def, status_set)
          }
        }

      } else {
        # explore dependencies of status-specific, if present
        miss_status <- setdiff(status_set, root_status)
        if (length(miss_status) == 0) {
          for (i in seq_len(dim(root_def)[1])) {
            private$.explore_dependencies(origin, root_def[i, ], status_set)
          }
        } else {
          stop_h(
            sprintf("Missing definition for augmented \"%s\"", root_name),
            sprintf("Provide a definition for status %s.",
                    paste(miss_status, collapse = ", "))
          )
        }
      }
    }
  ),

  ############################################################################
  ################################## ACTIVE ##################################
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
  ################################## PUBLIC ##################################
  ############################################################################

  public = list(
    #' @field defs `data.frame`\cr
    #' Augmented variable formatted definitions.
    defs = NULL,

    #' @description
    #' New object of class
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
      regexp <- c(
        # from : .augm("string"
        '.augm\\s*\\(\\s*["\'](.*?)["\']',
        # from : get_augmented("string"
        '.get_augm_variable\\s*\\(\\s*["\'](.*?)["\']',
        # from : .set_augm_value("string"
        '.set_augm_value\\s*\\(\\s*["\'](.*?)["\']'
      )
      # find required augmented variables from regular expressions
      output <- template_required_find(refer_file, regexp)

      c(set_source_name_(output, "model"), self$augm_required)
    },

    #' @description
    #' Remove a given augmented variable.
    #'
    #' @param x `character(1)`\cr
    #' Name of the augmented variable to remove.
    #' @param status `integer(1)`\cr
    #' The associated status of the definition.
    #' Default NA is for all infected status.
    #' @param show_warnings `logical(1)`\cr
    #' Whether to show warnings or not. Default is `TRUE`.
    #'
    #' @return `logical(1)`
    remove = function(x, status = NA, show_warnings = TRUE) {
      mask <- (self$defs[["name"]] == x)
      if (any(mask)) {
        if (is.na(status)) {
          self$defs <- self$defs[!mask, ]
        } else {
          # check if same name and status
          mask <- mask & (self$defs[["status"]] == status)
          mask[is.na(mask)] <- FALSE

          if (any(mask)) self$defs <- self$defs[!mask, ]
          else if (isTRUE(show_warnings)) warning_h(sprintf(
            "Unknown augmented variable \"%s\" for status %d to remove.",
            x, status
          ))
        }
      } else if (isTRUE(show_warnings)) {
        warning_h(sprintf("Unknown augmented variable \"%s\" to remove.", x))
      }
      return(any(mask))
    },

    #' @inherit add_augmented
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
                   min = NULL,
                   max = NULL,
                   rule = NA,
                   scale = NA,
                   from = NULL,
                   expr = NULL) {
      # handle (name, status) and check for existing definitions
      private$.helper_existing_definition(name, status)

      # auxiliary typed options for Rcpp build method
      opts <- list(
        name = as.character(name),
        status = as.integer(status),
        const = isTRUE(const),
        discr = isTRUE(discr),
        func = as.character(func),
        fcumul = FALSE, # UPDATE: to be implemented
        farg_key = NA, farg_val = NA, farg_at = NA,
        value = NA, value_at = NA,
        init_rand = isTRUE(init_rand),
        model = as.character(model),
        arg_key = NA, arg_val = NA, arg_at = NA,
        min = NA, min_at = NA,
        max = NA, max_at = NA,
        rule = as.character(rule),
        scale = 1,
        from = NA, from_at = NA
      )

      # augmented function
      if (!is.null(expr)) {
        opts <- private$.helper_implicit_func(opts, expr)
      } else if (!is.na(opts$func)) {
        opts <- private$.helper_function_arguments(opts, func_args)
        opts$init_rand <- opts$const <- opts$discr <- FALSE
      }

      # model arguments
      if (length(args) > 0) {
        opts <- private$.helper_model_arguments(opts, args)
      }

      if (is.na(opts$func)) {
        opts <- private$.helper_constant_or_random(opts, value, from)
        opts <- private$.helper_proposer_rule(opts, scale, from)
      }

      opts <- private$.helper_bounds(opts, min, max)
      self$defs <- rbind(self$defs, data.frame(opts))

      invisible(self)
    },

    #' @description
    #' Load augmented variable definitions either as direct arguments
    #' or from predefined list.
    #'
    #' @param ...
    #' Named arguments with respective list of definition details.
    #' @param from named `list`\cr
    #' Named list of parameter definitions.
    #' @param append `logical(1)`\cr
    #' If `TRUE` appends to current definitions.
    #' If `FALSE` drop all pre-existing definitions (default).
    load_list = function(...,
                         from = NULL,
                         append = FALSE) {
      definitions <- load_definition_list(..., from = from)

      # remove previos definitions if not append
      if (isFALSE(append)) self$initialize()

      for (i in seq_along(definitions)){
        def_name <- names(definitions[i])
        expr <- definitions[[i]]
        if (rlang::is_call(expr)) this_def <- list("expr" = expr)
        else this_def <- as.list(expr)

        # raw definition
        if (length(names(this_def)) == 0) {
          if (length(this_def) == 1) {
            # distinguish between constant value definition and expression
            names(this_def) <- ifelse(is.numeric(expr, "value", "expr"))
          } else {
            stop_h(
              sprintf("Wrong definition for augmented \"%s\".", def_name),
              "Provide scalar, expression or named list (see 'add_augmented')."
            )
          }
        }

        # associate name
        this_def$name <- def_name

        # multi-status definitions
        if (is.null(this_def$status)) do.call(self$add, this_def)
        else for (s in this_def$status) {
          do.call(self$add, modifyList(this_def, list(status = s)))
        }
      }

      invisible(self)
    },

    #' @description
    #' Load augmented variables from template definitions.
    #'
    #' @param templ_file
    #' Path/to/.refer template file.
    #' @param required `character()`\cr
    #' Vector containing required parameters.
    #' @param append `logical(1)`\cr
    #' If `TRUE` append missing augmented variable definitions.
    #' If `FALSE` drop all pre-existing definitions (default).
    load_template = function(templ_file, required = NULL, append = FALSE) {
      # check if template file exists
      if (file.exists((templ_file))) {
        # template augmented variable definitions
        templ_defs <- dget(templ_file)

        if (isTRUE(append)) {
          # check set of required augmented variables
          if (length(required) == 0) {
            warning_h(
              "Appending template augmented variable definitions failed.",
              "Alternatively use option 'append = FALSE'."
            )
          } else {
            miss <- template_missing_find(required, self$defs[["name"]])

            if (length(miss) == 0) warning_h(
              "No template augmented variable to append.",
              "Definitions are complete."
            )
            else if (!all(miss %in% names(templ_defs))) stop_h(
              "No template found for the missing augmented variables."
            )
            else self$load_list(from = templ_defs[miss], append = TRUE)
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
    #' @param required `character()`\cr
    #' Vector containing required parameters.
    #' @param status_set `integer()`\cr
    #' All statuses present in current data.
    assert = function(required, status_set) {
      status_set <- sort(status_set)
      defined <- self$defs[["name"]]
      # point out augmented variables that are not used
      template_not_used_warn(required, defined, "augmented")

      # stop if augmented variables are missing
      template_missing_stop(required, defined, "augmented")

      # explore status-specific definitions tree of dependencies
      if ("model" %in% names(required)) {
        for (r in required[["model"]]) private$.explore_dependencies(
          r, subset(self$defs, name == r), status_set
        )
      }

      invisible(self)
    }
  )
)
