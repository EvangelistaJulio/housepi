#      __
#     / /_  ____  __  __ ____ __  ____  __
#    / __ \/ __ \/ / / / ___/ _ \/ __ \/ /
#   / / / / /_/ / /_/ (__  )  __/ /_/ / /
#  /_/ /_/\____/\____/____/\___/ ____/_/
#                             /_/
#
# Copyright (C) 2025 CAMAGLIA Francesco - Institut Pasteur, Paris (France).
#
# parameter.R
#

############################################################################
############################################################################
############################# Housepi Parameter ############################
############################################################################
############################################################################

#' @title Housepi Parameter
#' @description
#' Help with parameter definitions of model Rcpp object.
#'
#' @keywords internal
#' @noRd
HousepiParameter <- R6::R6Class(
  classname = "HousepiParameter",
  cloneable = TRUE,

  private = list(
    deep_clone = function(name, value) {
      helper_clone(value)
    },

    # Helper to associate prior definition.
    .helper_prior = function(opts, prior, args) {
      # prior to tell if constant
      if (is.na(prior)) stop_h(
        sprintf("No prior provided for parameter \"%s\".", opts$name)
      )

      if (!is.null(args)) {
        # prior arguments
        parsed <- DefEncoder$new(as.list(args))
        prior <- parsed$combine_to_call(prior)
      }

      utils::modifyList(opts, list(const = FALSE, prior = prior))
    },

    # Helper with mcmc proposer rule.
    .helper_proposer = function(opts, discr, rule, scale) {
      opts$discr <- isTRUE(discr)

      # default proposer rule
      opts$rule <- "sample"
      if (!is.na(rule)) {
        opts$rule <- rule
        opts$scale <- ifelse(is.na(scale), 1, scale)
      }
      opts
    }
  ),

  ############################################################################
  ################################## PUBLIC ##################################
  ############################################################################

  public = list(
    #' @field defs `data.frame`\cr
    #' Parameter formatted definitions.
    defs = NULL,

    #' @description
    #' New object of class
    #' [HousepiParameter][housepi::HousepiParameter].
    initialize = function() {
      self$defs <- NULL
      invisible(self)
    },

    #' @description
    #' Find parameters required by the model.
    #'
    #' @param refer_file
    #' Path/to/.refer source file.
    require = function(refer_file) {
      # vector of regular expression to use for assertion
      regexp <- c(
        # from : param("string" or get_param('string'
        'param\\s*\\(\\s*["\'](.*?)["\']',
        # from : param_pointer("string" or param_pointer('string'
        'param_pointer\\s*\\(\\s*["\'](.*?)["\']',
        # from : emplace_contact_rate("string"
        'emplace_contact_rate\\s*\\(\\s*["\'](.*?)["\']',
        # from : emplace_group_rate(..., "string"
        'emplace_group_rate\\s*\\(\\s*[^,]*\\s*,\\s*["\'](.*?)["\']'
      )
      # retrieve required parameters from regular expressions
      set_source_name_(
        template_required_find(refer_file, regexp),
        "model"
      )
    },

    #' @description
    #' Remove a given parameter.
    #'
    #' @param name `character(1)`\cr
    #' Name of the parameter to remove.
    #' @param show_warnings `logical(1)`\cr
    #' If `TRUE` show warnings if parameter not found (default).
    #' If `FALSE` do not show warnings.
    remove = function(name,
                      show_warnings = TRUE) {
      mask <- (self$defs[["name"]] == name)
      if (any(mask)) self$defs <- self$defs[!mask, ]
      else if (isTRUE(show_warnings)) warning_h(
        sprintf("Attempt to remove unknown parameter \"%s\".", name)
      )
      return(any(mask))
    },

    #' @inherit add_parameter
    add = function(name,
                   const = NA,
                   discr = FALSE,
                   value = NULL,
                   prior = NA_character_,
                   args = NULL,
                   min = NULL,
                   max = NULL,
                   rule = NA_character_,
                   scale = NA_real_) {
      if (is.na(name)) stop_h(
        "Incomplete parameter definition: missing name."
      )
      # remove previous definition if it exists
      if (self$remove(name, FALSE)) warning_h(
        sprintf("Overwriting definition of parameter \"%s\".", name)
      )

      opts <- list(name = name,
                   const = TRUE, discr = FALSE,
                   value = NA_real_, prior = NA_character_,
                   min = NA_real_, max = NA_real_,
                   rule = NA_character_, scale = NA_real_)

      # assume it's a constant if a prior is not provided
      if (is.na(const)) const <- is.na(prior)

      # parse initial value
      if (is.na(DefEncoder$new(value)$x)) opts$value <- NA_real_
      else opts$value <- DefEncoder$new(value)$x

      if (isTRUE(const)) {
        # check that a value is provided
        if (is.na(opts$value)) stop_h(
          sprintf("Missing value for constant parameter \"%s\".", opts$name)
        )
      } else {
        opts <- private$.helper_prior(opts, prior, args)
        opts$min <- DefEncoder$new(min)$x
        opts$max <- DefEncoder$new(max)$x
        opts <- private$.helper_proposer(opts, discr, rule, scale)
      }

      # append as a row
      self$defs <- rbind(self$defs, data.frame(opts))

      invisible(self)
    },

    #' @description
    #' Load parameter definitions either as direct arguments
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
      # remove previous definitions if not append
      if (isFALSE(append)) self$initialize()

      # load definition list
      definitions <- load_definition_list(..., from = from)

      for (i in seq_along(definitions)){
        this_def <- as.list(definitions[[i]])
        def_name <- names(definitions[i])

        if (length(names(this_def)) == 0) {
          if (length(this_def) == 1) names(this_def) <- "value"
          else stop_h(
            sprintf("Wrong definition for parameter \"%s\"", def_name),
            "Provide a definition list or a scalar (if constant)."
          )
        }
        this_def$name <- def_name
        do.call(self$add, this_def)
      }

      invisible(self)
    },

    #' @description
    #' Load missing parameter from template definitions.
    #'
    #' @param templ_file
    #' Path/to/.refer template file.
    #' @param required `character()`\cr
    #' Vector containing required parameters.
    #' @param append `logical(1)`\cr
    #' If `TRUE` append missing parameter definitions.
    #' If `FALSE` drop all pre-existing definitions (default).
    load_template = function(templ_file,
                             required = NULL,
                             append = FALSE) {
      # check if template file exists
      if (file.exists((templ_file))) {
        # template parameter definitions
        templ_defs <- dget(templ_file)

        if (isTRUE(append)) {
          # check the set of required parameters
          if (length(required) == 0) {
            warning_h("No parameter definitions required by the model.")
          } else {
            miss <- template_missing_find(required, self$defs[["name"]])
            if (length(miss) == 0) warning_h(
              "No template parameters to append.",
              "Definitions are complete."
            )
            else if (!all(miss %in% names(templ_defs))) stop_h(
              "No template found for the missing parameters."
            )
            else self$load_list(from = templ_defs[miss], append = TRUE)
          }
        } else {
          self$load_list(from = templ_defs, append = FALSE)
        }

      } else {
        stop_h(
          "Impossible to automatically initialize parameters.",
          "A file with template definitionse was not found."
        )
      }

      invisible(self)
    },

    #' @description
    #' Assert parameter definitions match the requirement list.
    #'
    #' @param required `character()`\cr
    #' Vector containing required parameters.
    assert = function(required) {
      defined <- self$defs[["name"]]
      # point out parameters that are not used
      template_not_used_warn(required, defined, "parameter")
      # stop if parameters are missing
      template_missing_stop(required, defined, "parameter")
      invisible(self)
    }
  )
)
