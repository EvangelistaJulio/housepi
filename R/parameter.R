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

#' @title Housepi Parameter R6 Class
#' @description
#' Help with parameter definitions of model Rcpp object.
#'
#' @keywords internal
HousepiParameter <- R6::R6Class(
  classname = "HousepiParameter",
  cloneable = TRUE,

  ############################################################################
  ################################## PUBLIC ##################################
  ############################################################################

  public = list(
    #' @field defs `data.frame`\cr
    #' Parameter formatted definitions.
    defs = NULL,

    #' @description
    #' New instance of class
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
      if (any(mask)) {
        # drop previous definition
        self$defs <- self$defs[!mask, ]
        return(TRUE)

      } else {
        if (isTRUE(show_warnings)) {
          warning_h(
            sprintf("Attempt to remove unknown parameter \"%s\".", name)
          )
        }
      }

      return(FALSE)
    },

    #' @description
    #' Add a new parameter definition.
    #'
    #' @param name `character(1)`\cr
    #' Name of the parameter.
    #' @param const `logical(1)`\cr
    #' If `FALSE` the parameter is treated as a constant.
    #' @param discr `logical(1)`\cr
    #' If `TRUE` the parameter is discrete.
    #' @param value `double(1)`\cr
    #' Initial value of the parameter.
    #' @param prior `character(1)`\cr
    #' Name of the prior distribution.
    #' @param args named `list()`\cr
    #' Arguments for the prior distribution.
    #' @param rule `character(1)`\cr
    #' Name of the inference proposal rule.
    #' @param scale `double(1)`\cr
    #' Scale parameter for the proposal rule.
    #' Changed automatically during inference if `n_adpat`>0.
    #' @param min `double(1)`\cr
    #' Value of the minimum.
    #' @param max `double(1)`\cr
    #' Value of the maximum.
    add = function(name,
                   const = NA,
                   discr = FALSE,
                   value = NA,
                   prior = NA,
                   args = list(),
                   rule = NA,
                   scale = NA,
                   min = NA,
                   max = NA) {
      if (is.na(name)) {
        stop_h("Incomplete definition: missing name.")
      }

      # remove previous definition if it exists
      if (self$remove(name, FALSE)) {
        warning_h(
          sprintf("Overwriting \"%s\" definition.", name)
        )
      }

      # assume it's a constant if a prior is not provided
      if (is.na(const)) {
        const <- is.na(prior)
      }

      if (isTRUE(const)) {
        # check that value is provided
        if (is.na(value)) {
          stop_h(
            sprintf("Missing value for constant parameter \"%s\".", name)
          )
        }
        # append constant parameter
        self$defs <- rbind(
          self$defs,
          data.frame(
            name = name,
            const = TRUE,
            discr = is.integer(value),
            value = value,
            prior = NA,
            arg_key = NA,
            arg_val = NA,
            rule = NA,
            scale = NA,
            min = NA,
            max = NA,
            stringsAsFactors = FALSE
          )
        )

      } else {
        # prior to tell if constant
        if (is.na(prior)) {
          stop_h(
            sprintf("No prior provided for parameter \"%s\".", name)
          )
        }

        # prior arguments
        arg_key_ <- NA
        arg_val_ <- NA
        if (length(args) > 0) {
          arg_key_ <- type_join_(names(args))
          arg_val_ <- type_join_(unname(unlist(args)))
        }

        # mcmc proposer rule
        if (is.na(rule)) {
          # default proposer rule
          rule <- "draw"
          scale <- NA
        }

        # append parameter to be inferred
        self$defs <- rbind(
          self$defs,
          data.frame(
            name = name,
            const = FALSE,
            discr = isTRUE(discr),
            value = value,
            prior = prior,
            arg_key = arg_key_,
            arg_val = arg_val_,
            min = min,
            max = max,
            rule = rule,
            scale = scale,
            stringsAsFactors = FALSE
          )
        )
      }

      invisible(self)
    },

    #' @description
    #' Load parameter definitions either as direct arguments
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
      # load definition list
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
              sprintf("Wrong definition for parameter \"%s\"", def_name),
              "Provide a definition list or a scalar (if constant)."
            )
          }
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
    #' @param requirements `character()`\cr
    #' Vector containing required parameters.
    #' @param append `logical(1)`\cr
    #' If `TRUE` append missing parameter definitions.
    #' If `FALSE` drop all pre-existing definitions (default).
    load_template = function(templ_file,
                             requirements = NULL,
                             append = FALSE) {
      # check if template file exists
      if (file.exists((templ_file))) {
        # template parameter definitions
        templ_defs <- dget(templ_file)

        if (isTRUE(append)) {
          # check the set of required parameters
          if (length(requirements) == 0) {
            warning_h("No parameter definitions required by the model.")
          } else {
            missing_set <- template_find_missing(requirements,
                                                 self$defs[["name"]])
            if (length(missing_set) == 0) {
              warning_h(
                "No template parameters to append.",
                "Definitions are complete."
              )
            } else if (!all(missing_set %in% names(templ_defs))) {
              stop_h(
                "No template found for the missing parameters."
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
          "Impossible to automatically initialize parameters.",
          "A file with template definitionse was not found."
        )
      }

      invisible(self)
    },

    #' @description
    #' List containing current parameter definitions.
    #'
    #' @return `list()`
    as_list = function() {
      output <- list()
      for (r in seq_len(dim(self$defs)[1])) {
        this_par <- as.list(self$defs[r, ])
        this_par$name <- NULL

        if (!is.na(this_par$arg_val)) {
          args <- stringr::str_split(this_par$arg_val, "/")[[1]]
          names(args) <- stringr::str_split(this_par$arg_key, "/")[[1]]
          this_par$args <- args
        }
        this_par$arg_key <- NULL
        this_par$arg_val <- NULL

        this_par[is.na(this_par)] <- NULL

        output[[self$defs[r, ]$name]] <- this_par
      }

      return(output)
    },

    #' @description
    #' Assert parameter definitions match the requirement list.
    #'
    #' @param requirements `character()`\cr
    #' Vector containing required parameters.
    assert = function(requirements) {
      defined <- self$defs[["name"]]

      # point out parameters that are not used
      unused_set <- template_find_unused(requirements, defined)
      # drop unused parameters
      if (length(unused_set) > 0) {
        warning_h(
          sprintf("Dropping unused parameter%s %s",
                  ifelse(length(unused_set) > 1, "s:", ""),
                  character_to_string_list(unused_set)),
          find_closest_match(unused_set, requirements)
        )
        for (name in unused_set) self$remove(name)
      }

      # find missing parameters
      missing_set <- template_find_missing(requirements, defined)
      if (length(missing_set) > 0) {
        stop_h(
          sprintf("Missing parameter definition%s for:",
                  ifelse(length(missing_set) > 1, "s", "")),
          character_to_string_list(missing_set)
        )
      }

      invisible(self)
    }
  )
)
