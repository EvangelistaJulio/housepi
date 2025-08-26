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
  at_tuples = list(c("farg_val", "farg_at"),
                  c("value", "value_at"), c("arg_val", "arg_at"),
                  c("min", "min_at"), c("max", "max_at"),
                  c("scale", "scale_at"), c("from", "from_at"))

  output <- c()
  if (length(defs) > 0) {
    for (tuple in at_tuples) {
      xdef <- defs[tuple]
      for (i in seq_len(dim(xdef)[[1]])) {
        rw <- unlist(xdef[i, ])
        if (!any(is.na(rw))) {
          arg2par <- data.frame(
            val = strsplit(rw[1], "/")[[1]],
            at = strsplit(rw[2], "/")[[1]]
          )
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
    .handle_existing_definition = function(name, status) {
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
    .handle_implicit_func = function(opts, expr) {
      parsed <- extract_func_(expr)
      if (!is.na(parsed$func)) {
        if (is.na(opts$func)) {
          opts$func <- parsed$func
          opts$farg_key <- parsed$key
          opts$farg_val <- parsed$val
          opts$farg_at  <- parsed$at
        } else {
          stop_h(
            sprintf("Redundant function definitions for \"%s\".", opts$name),
            "Use either 'expr' or 'func' argument."
          )
        }
      }

      if (!is.na(opts$func)) {
        opts$init_rand <- opts$const <- opts$discr <- FALSE
      } else if (is.na(opts$const)) {
        opts$const <- is.na(opts$model) && is.na(opts$rule)
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
        if (is.na(value)) stop_h(
          "Missing value for constant augmented variable."
        )
        opts$value <- type_split_(value)$val
        opts$value_at <- type_split_(value)$at
      } else if (isTRUE(opts$init_rand)) {
        if (is.na(opts$model)) stop_h(
          "Unknown model for random initialization."
        )
      } else {
        if (is.na(value)) {
          if (is.na(opts$model)) {
            if (!is.na(from)) {
              opts$value <- type_split_(from)$val
              opts$value_at <- type_split_(from)$at
            }
            # else: constant NA augmetned variable
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
          if (is.na(opts$model)) opts$const <- TRUE
          else opts$rule <- "sample"
        } else if (opts$rule != "sample") {
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

    # Definitions required by augmented variable definitions
    # according to the "at" key of the definition (A, D, P).
    .find_def_required = function(at) {
      required <- extract_at_values(self$defs, at)
      # add source name
      return(set_source_name_(required, "augmented"))
    },

    # recursive explore of the dependencies tree to point out errors
    .explore_dependencies = function(root_def, status_set) {
      root_name <- unique(root_def$name) # name of the root definition
      root_status <- root_def$status     # status of the root definition
      is_generic <- is.na(root_status)   # wheter generic or status-specific
      
      if (any(is_generic) && length(root_status) > 1) {
        # assert for unique definition
        stop_h(
          sprintf("Equivocal definitions for augmented \"%s\"", root_name),
          ifelse(is.all(is_generic),
                 "Impossible to provide multiple generic definitions.",
                 "Impossible to provide both generic and status-specific definitions.")
        )

      } else if (length(root_status) == 1) {
        # check all depedendecies
        for (r in extract_at_values(root_def, "A")){
          branch_def <- subset(self$defs, name == r)
          branch_status <- branch_def$status   # status of the branch definition

          # check if at status specific tips
          is_tip <- (length(branch_status) == 1) && 
                    !identical(branch_status, root_status)

          if (is_tip) stop_h(sprintf(
            "Missing definition for augmented \"%s\"%s.",
            r, ifelse(is_generic, "", sprintf(" status %s", root_status))
          ))
          else {
            if (!is_generic) status_set <- root_status
            private$.explore_dependencies(branch_def, status_set)
          }
        }
        
      } else {
        # explore dependencies of status-specific, if present
        miss_status <- setdiff(status_set, root_def$status)
        if (length(miss_status) == 0) {
          for (i in seq_len(dim(root_def)[1])) {
            private$.explore_dependencies(root_def[i, ], status_set)
          }
        } else stop_h(
          sprintf("Missing definition for augmented \"%s\"", root_name),
          sprintf("Provide a definition for status %s.", 
                  paste(miss_status, collapse = ", "))
        )
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

    #' @field dataset_required `character()`\cr
    #' Find data required by augmented variable definitions.
    dataset_required = function(value) {
      if (missing(value)) return(private$.find_def_required("D"))
      stop_h("'dataset_required' is read only.")
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
      regular_expr <- c(
        # from : .augm("string"
        '.augm\\s*\\(\\s*["\'](.*?)["\']',
        # from : .get_augm_variable("string"
        '.get_augm_variable\\s*\\(\\s*["\'](.*?)["\']',
        # from : .set_augm_value("string"
        '.set_augm_value\\s*\\(\\s*["\'](.*?)["\']'
      )
      # find required augmented variables from regular expressions
      output <- template_required_find(refer_file, regular_expr)
      return(c(set_source_name_(output, "model"), self$augm_required))
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
        if (is.na(status)) self$defs <- self$defs[!mask, ]
        else {
          # check if same name and status
          mask <- mask & (self$defs[["status"]] == status)
          mask[is.na(mask)] <- FALSE

          if (any(mask)) self$defs <- self$defs[!mask, ]
          else if (isTRUE(show_warnings)) warning_h(sprintf(
            "Unknown augmented variable \"%s\" for status %d to remove.",
            x, status
          ))
        }
      } else if (isTRUE(show_warnings)) warning_h(
        sprintf("Unknown augmented variable \"%s\" to remove.", x)
      )
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
                   min = NA,
                   max = NA,
                   rule = NA,
                   scale = NA,
                   from = NA,
                   expr = NA) {
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

      # agumented function 
      if (!is.na(expr)) {
        opts <- private$.handle_implicit_func(opts, expr)

      } else if (!is.na(opts$func)){
        if (!is.null(func_args)) {
          opts$farg_key <- type_join_(names(func_args))
          opts$farg_val <- type_join_(type_split_(func_args)$val)
          opts$farg_at  <- type_join_(type_split_(func_args)$at)
        }
        opts$init_rand <- opts$const <- opts$discr <- FALSE
      }

      # model arguments
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
        this_def <- as.list(definitions[[i]])
        def_name <- names(definitions[i])
        if (length(names(this_def)) == 0) {
          if (length(this_def) == 1) names(this_def) <- "expr"
          else stop_h(
            sprintf("Wrong definition for augmented variable \"%s\".", 
                    def_name),
            "Provide a scalar, an expression or a named list to pass to 'add_augmented'."
          )
        }
        # associate name
        this_def$name <- def_name
        
        # multi status definitions
        if (is.null(this_def$status)) do.call(self$add, this_def)
        else {
          status_def <- this_def
          for (s in this_def$status) {
            status_def$status <- s
            do.call(self$add, status_def)
          }
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
    #' All statuses present in current dataset.
    assert = function(required, status_set) {
      status_set <- sort(status_set)
      defined <- self$defs[["name"]]
      # point out augmented variables that are not used
      template_not_used_warn(required, defined, "augmented variable")

      # stop if augmented variables are missing
      template_missing_stop(required, defined, "augmented variable")

      # explore status-specific definitions tree of depedencies
      if ("model" %in% names(required)) {
        for (r in required[["model"]]) {
          private$.explore_dependencies(subset(self$defs, name == r),
                                        status_set)
        }
      }

      invisible(self)
    }
  )
)
