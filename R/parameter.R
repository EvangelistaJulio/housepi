#
# parameter.R
#
# Copyright (C) 2025 Francesco Camaglia - Institut Pasteur, Paris (France).
#

#' @title Housepi Parameter R6 Class
#' @description
#' Help with parameter definitions of model Rcpp object.
#' @export HousepiParameterClass
#' @docType class
#' @noRd
HousepiParameterClass <- R6::R6Class(
  classname = "HousepiParameterClass",
  cloneable = FALSE,

  ###################
  ##### PUBLIC ######
  ###################

  public = list(
    defs = NULL,

    #' @description
    #' New instance of class [HousepiParameterClass][housepi::HousepiParameterClass].
    #' @noRd
    initialize = function() {
      self$drop_all()
      invisible(self)
    },

    #' @description
    #' Find parameters required by the model.
    #' @param refer_file
    #' Path/to/.refer source file.
    #' @noRd
    require = function(refer_file) {
      # vector of regular expression to use for assertion
      regular_expr <- c(
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
      return(template_find_requirement(refer_file, regular_expr))                                  
    },

    # @description
    # Remove a given parameter.
    remove = function(name,
                      show_warnings = TRUE) {
      mask <- (self$defs[["name"]] == name)
      if (any(mask)) {
        # drop previous definition
        self$defs <- self$defs[!mask, ]
        return(TRUE)

      } else {
        if (isTRUE(show_warnings)) {
          cli::cli_alert_warning(cli::format_warning(sprintf(
            "Unknown parameter \"%s\" to remove.", name
          )))
        }
      }

      return(FALSE)
    },

    #' @description
    #' Add a new parameter definition.
    #' @noRd
    add = function(name,
                   value = NA,
                   constant = NA,
                   discrete = FALSE,
                   prior = NA,
                   args = list(),
                   rule = NA,
                   rule_scale = NA,
                   min = NA,
                   max = NA) {
      if (is.na(name)) {
        stop(cli::format_error("Incomplete definition: missing name."))
      }

      # remove previous definition if it exists
      if (self$remove(name, FALSE)) {
        cli::cli_alert_warning(cli::format_warning(sprintf(
          "Existing \"%s\" definition has been overwritten.", name
        )))
      }

      # assume it's a constant if a prior is not provided
      if (is.na(constant)) {
        constant <- is.na(prior)
      }

      if (isTRUE(constant)) {
        # check that value is provided
        if (is.na(value)) {
          stop(cli::format_error(
            sprintf("Missing value for constant parameter \"%s\".", name)
          ))
        }
        # append constant parameter
        self$defs <- rbind(
          self$defs,
          data.frame(
            name = name,
            value = value,
            constant = TRUE,
            discrete = is.integer(value),
            prior = NA,
            arg_name = NA,
            arg_value = NA,
            rule = NA,
            rule_scale = NA,
            min = NA,
            max = NA,
            stringsAsFactors = FALSE
          )
        )

      } else {
        # assign prior
        if (is.na(prior)) {
          stop(cli::format_error(sprintf(
            "Incomplete definition: provide prior name for \"%s\".", name
          )))
        }
        # prior arguments
        arg_name_ <- NA
        arg_value_ <- NA
        if (length(args) > 0) {
          arg_name_ <- str_join_(names(args))
          arg_value_ <- str_join_(unname(unlist(args)))
        }

        # mcmc proposer rule
        if (is.na(rule)) {
          # default proposer rule
          rule <- "gibbs"
          rule_scale <- NA
        }

        # append parameter to be inferred
        self$defs <- rbind(
          self$defs,
          data.frame(
            name = name,
            value = value,
            constant = FALSE,
            discrete = isTRUE(discrete),
            prior = prior,
            arg_name = arg_name_,
            arg_value = arg_value_,
            rule = rule,
            rule_scale = rule_scale,
            min = min,
            max = max,
            stringsAsFactors = FALSE
          )
        )
      }

      invisible(self)
    },

    #' @description
    #' Drop all current parameter definitions.
    #' @noRd
    drop_all = function() {
      self$defs <- data.frame(
        name = character(),
        value = numeric(),
        constant = logical(),
        discrete = logical(),
        prior = character(),
        arg_name = character(),
        arg_value = character(),
        rule = character(),
        rule_scale = numeric(),
        min = numeric(),
        max = numeric(),
        stringsAsFactors = FALSE
      )
      invisible(self)
    },

    #' @description
    #' Associate parameters from a list. Previous definitions are dropped.
    #' @param def_list (named `list()`)\cr
    #' A list containg parameter names with respective list of definitions.
    #' @noRd
    from_list = function(definitions,
                         append = TRUE) {
      if (!is.list(definitions)) {
        stop(cli::format_error("Argument {.code def_list} should be a list."))
      }
      if (isFALSE(append)) {
        self$drop_all()
      }

      for (i in seq_along(definitions)){
        this_def <- as.list(definitions[[i]])
        def_name <- names(definitions[i])
        if (is.null(names(this_def))) {
          if (length(this_def) == 1) {
            names(this_def) <- "value"
          } else {
            stop(cli::format_error(c(
              sprintf("Wrong definition for parameter \"%s\"", def_name),
              i = "Provide a definition list or a scalar (if constant)."
            )))
          }
        }
        this_def[["name"]] <- def_name
        do.call(self$add, this_def)
      }

      invisible(self)
    },

    #' @description
    #' Add missing parameter from template definitions.
    #' @param templ_file
    #' Path/to/.refer template file.
    #' @param requirements (`character()`)\cr
    #' Vector containing required parameters.
    #' @param append (`logical(1)`)\cr
    #' If `TRUE` append to previous definitions (default).
    #' If `FALSE` drop previous definitions.
    #' @noRd
    template = function(templ_file,
                        requirements = NULL,
                        append = TRUE) {
      # check if template file exists
      if (file.exists((templ_file))) {
        # template parameter definitions
        templ_defs <- dget(templ_file)

        if (isTRUE(append)) {
          # check the set of required parameters
          if (is.null(requirements)) {
            cli::cli_alert_warning(cli::format_warning(c(
              "Appending template parameter definitions failed.",
              i = "Alternatively use option {.code append = FALSE}."
            )))
          } else {
            missing_set <- template_find_missing(requirements,
                                                 self$defs[["name"]])
            self$from_list(templ_defs[missing_set], append = TRUE)
          }
        } else {
          self$from_list(templ_defs, append = FALSE)
        }

      } else {
        stop(cli::format_error(c(
          "Impossible to automatically initialize parameters.",
          i = "A file with template definitionse was not found."
        )))
      }

      invisible(self)
    },

    #' @description
    #' Assert parameter definitions match the requirement list.
    #' @param requirements (`character()`)\cr
    #' Vector containing required parameters.
    #' @return (`logical(1)`)\cr
    #' @noRd
    assert = function(requirements) {
      # logical state for complete definitions
      is_def_complete <- TRUE

      if (is.null(requirements)) {
        # if requirement set is NULL it means nothing can be said
        is_def_complete <- NA

      } else {
        # unique and ordered set of parameters
        requirements <- unique(requirements[order(requirements)])

        # find missing parameters
        missing_set <- template_find_missing(requirements,
                                             self$defs[["name"]])
        if (length(missing_set) > 0) {
          is_def_complete <- FALSE
          oss <- paste0("Missing parameter(s): ",
                        character_to_string_list(missing_set), ".")
          cli::cli_alert_warning(cli::format_warning(oss))
        }

        # point out parameters that are not used
        unused_set <- template_find_unused(requirements,
                                           self$defs[["name"]])
        # drop unused parameters
        if (length(unused_set) > 0) {
          oss <- paste0("Dropping not used parameter(s): ",
                        character_to_string_list(unused_set), ".")
          cli::cli_alert_warning(cli::format_warning(oss))

          for (name in unused_set) {
            self$remove(name)
          }
        }
      }

      return(is_def_complete)
    }
  )
)
