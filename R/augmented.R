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
      parsed <- DefEncoder$new(expr)
      if (!is.na(parsed$x)) {
        if (is.na(opts$func)) opts$func <- parsed$x
        else stop_h(
          sprintf("Mutiple function definitions for \"%s\".", opts$name),
          "Use either 'expr' or 'func' argument."
        )
      }

      if (!is.na(opts$func)) {
        opts$const <- FALSE
        opts$discr <- FALSE
      } else if (is.na(opts$const)) {
        opts$const <- is.na(opts$model) && is.na(opts$rule)
      }

      opts
    },

    # process function arguments
    .helper_function_arguments = function(opts, func_args) {
      if (!is.null(func_args)) {
        parsed <- DefEncoder$new(as.list(func_args))
        opts$func <- parsed$combine_to_call(opts$func)
      }
      opts
    },

    # process model arguments
    .helper_model_args = function(opts, args) {
      if (is.na(opts$model))  {
        msg <- sprintf("Spare arguments for augmented \"%s\".", opts$name)
        if (!is.na(opts$func)) msg <- c(msg, "Did you mean 'func_args'?")
        stop_h(msg)

      } else if (!is.null(args)) {
        parsed <- DefEncoder$new(as.list(args))
        opts$model <- parsed$combine_to_call(opts$model)
      }

      opts
    },

    # Process constant or random.
    .helper_const_vs_rand = function(opts, value) {
      parsed <- DefEncoder$new(value)
      opts$value <- parsed$x

      if (isFALSE(opts$const) && is.na(opts$value) && is.na(opts$model)) stop_h(
        sprintf("Impossible to randomly initialize undefined augmented \"%s\".",
                opts$name),
        "Provide one between value and model."
      )
      else if (is.na(opts$const)) opts$const <- TRUE

      opts
    },

    # Interpret user choice for proposer rule.
    .helper_proposer_rule = function(opts, scale) {
      # sampling rule but no constant parameter provided
      if (is.na(opts$const)) {
        if (!(is.na(opts$rule) && is.na(opts$model))) {
          opts$const <- FALSE
        }
      }

      if (isFALSE(opts$const)) {
        if (is.na(opts$rule)) {
          if (is.na(opts$model)) opts$const <- TRUE
          else opts$rule <- "sample"
        } else if (opts$rule == "sample") {
          if (is.na(opts$model)) stop_h(
            sprintf("No model to sample from for augmented \"%s\".",
                    opts$name)
          )
        } else {
          if (is.na(scale)) opts$scale <- 1
          else if (is.numeric(scale) && (scale > 0)) opts$scale <- scale
          else stop_h("Rule scale should be a numeric value greater than 0.")
        }
      } else {
        opts$rule <- NA
      }

      opts
    },

    # Extract value of given type from augmented definitions.
    .extract_types = function(type) {
      # columns with types amongst augmented definitions
      output <- c()
      if (length(self$defs) > 0) {
        for (typed_col in c("value", "func", "model", "min", "max")) {
          for (x in self$defs[[typed_col]]) {
            output <- c(output, walk_encoded_func(x, type))
          }
        }
      }
      unique(output)
    },

    # Definitions required by augmented variable definitions
    # according to the "at" key of the definition (A, D, P).
    .find_def_required = function(type) {
      required <- private$.extract_types(type)
      # add source name
      set_source_name_(required, "augmented")
    }
  ),

  ############################################################################
  ################################## ACTIVE ##################################
  ############################################################################

  active = list(

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
    },

    #' @field group_rate_required `character()`\cr
    #' Find group rates required by augmented variable definitions.
    group_rate_required = function(value) {
      if (missing(value)) return(private$.find_def_required("G"))
      stop_h("'group_rate_required' is read only.")
    },

    #' @field param_required `character()`\cr
    #' Find parameters required by augmented variable definitions.
    param_required = function(value) {
      if (missing(value)) return(private$.find_def_required("P"))
      stop_h("'param_required' is read only.")
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
    #' @param x
    #' A `character(1)` with the name of the augmented definition to remove.
    #' @param status
    #' An `integer(1)` with the associated status of the definition.
    #' @param show_warnings
    #' `logical(1)`to show warnings or not.
    #'
    #' @return `logical(1)`
    remove = function(x, status, show_warnings = TRUE) {
      mask <- self$defs[["name"]] == x
      if (any(mask)) {
        # check if same name and status
        mask <- mask & sapply(self$defs[["status"]], function(x) {
          identical(x, as.integer(status))
        })

        if (any(mask)) self$defs <- self$defs[!mask, ]
        else if (isTRUE(show_warnings)) warning_h(sprintf(
          "Unknown augmented variable \"%s\" for status %d to remove.",
          x, status
        ))
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
                   model = NA,
                   args = NULL,
                   min = NULL,
                   max = NULL,
                   rule = NA,
                   scale = NA,
                   expr = NULL) {
      # handle (name, status) and check for existing definitions
      private$.helper_existing_definition(name, status)

      # auxiliary typed options for Rcpp build method
      opts <- list(
        name = as.character(name),
        status = as.integer(status),
        const = const,
        discr = isTRUE(discr),
        func = as.character(func),
        fcumul = FALSE,
        value = NA,
        model = as.character(model),
        min = NA,
        max = NA,
        rule = as.character(rule),
        scale = NA
      )

      # augmented function
      if (!is.null(expr)) {
        opts <- private$.helper_implicit_func(opts, expr)
      } else if (!is.na(opts$func)) {
        opts <- private$.helper_function_arguments(opts, func_args)
        opts$const <- FALSE
        opts$discr <- FALSE
      }

      # model arguments
      if (length(args) > 0) opts <- private$.helper_model_args(opts, args)

      if (is.na(opts$func)) {
        opts <- private$.helper_proposer_rule(opts, scale)
        opts <- private$.helper_const_vs_rand(opts, value)
      }

      opts$min <- DefEncoder$new(min)$x
      opts$max <- DefEncoder$new(max)$x
      self$defs <- rbind(self$defs, data.frame(opts))

      invisible(self)
    },

    #' @description
    #' Load augmented variable definitions either as direct arguments
    #' or from predefined list.
    #'
    #' @param ...
    #' Named arguments with respective list of definition details.
    #' @param from
    #' A named `list` of parameter definitions.
    #' @param append
    #' A `logical(1)` to append to current definitions.
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
          # distinguish between constant value definition and expression
          if (length(this_def) == 1) {
            if (is.na(expr)) names(this_def) <- "value"
            else if (is.numeric(expr)) names(this_def) <- "value"
            else names(this_def) <- "expr"
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
          do.call(self$add, utils::modifyList(this_def, list(status = s)))
        }
      }

      invisible(self)
    },

    #' @description
    #' Load augmented variables from template definitions.
    #'
    #' @param templ_file
    #' Path/to/.refer template file.
    #' @param required
    #' A `character()` containing required parameters.
    #' @param append
    #' `logical(1)` to append missing augmented variable definitions.
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
    #' @param required
    #' A `character()` containing required parameters.
    #' @param status_set
    #' A `integer()` with all statuses present in the current data.
    assert = function(required, status_set) {
      status_set <- sort(status_set)
      defined <- self$defs[["name"]]
      # point out augmented variables that are not used
      template_not_used_warn(required, defined, "augmented")

      # stop if augmented variables are missing
      template_missing_stop(required, defined, "augmented")

      invisible(self)
    }
  )
)
