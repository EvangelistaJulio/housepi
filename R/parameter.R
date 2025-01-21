#
# parameter.R
#
# Copyright (C) 2025 Francesco Camaglia - Institut Pasteur, Paris (France).
#

# @title Housepi Parameter R6 Class
# @description
# Help with parameter definitions of model Rcpp object.
# @export HousepiParameterClass
# @docType class
HousepiParameterClass <- R6::R6Class(
  classname = "HousepiParameterClass",
  cloneable = FALSE,

  ###################
  ##### PUBLIC ######
  ###################

  public = list(
    defs = NULL,

    # @description
    # Assign the empty parameter data.frame
    initialize = function() {
      self$drop_all()
      invisible(self)
    },

    # @description
    # Vector of required parameters.
    # @param refer_file
    # Path/to/.refer source file.
    # @param augm_frame
    # Data.frane containing 
    requisites = function(refer_file,
                          augm_frame) {
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
      require_set <- template_find_requirement(refer_file,
                                          regular_expr)                                  

      # from : augmented variable model argument definitions
      augm_require_set <- c()
      if (length(augm_frame) > 0) {
        x <- augm_frame[c("arg_value", "arg_type")]
        for (i in seq_along(x)) {
          aux <- data.frame(vals = strsplit(x[i, 1], "/")[[1]],
                            keys = strsplit(x[i, 2], "/")[[1]])
          aux <- aux[aux[["keys"]] == "P", "vals"]
          augm_require_set <- c(augm_require_set, aux[!is.na(aux)])
        }
      }

      # format output
      if (is.null(require_set)) {
        if (length(augm_require_set) == 0) {
          require_set <- NULL
        } else {
          require_set <- augm_require_set
          require_set <- unique(augm_require_set[order(augm_require_set)])
        }
      } else {
        require_set <- c(require_set, augm_require_set)
        require_set <- unique(require_set[order(require_set)])
      }

      return(require_set)
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
          cli::cli_alert_warning(cli::format_warning(c(
            sprintf("Unknown parameter \"%s\" to remove.", name
          ))))
        }
      }

      return(FALSE)
    },

    # @description
    # Add a new parameter definition.
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
        cli::cli_alert_warning(cli::format_warning(c(
          sprintf(
            "Existing \"%s\" definition has been overwritten.",
            name
        ))))
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
          stop(cli::format_error(
            sprintf("Incomplete definition: provide prior name for \"%s\".", name)
          ))
        }
        # prior arguments
        if (length(args) == 0) {
          arg_name_ <- NA
          arg_value_ <- NA
        } else {
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

    # @description
    # Drop all current parameter definitions.
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

    # @description
    # Associate parameters from a list. Previous definitions are dropped.
    # @param def_list (named `list()`)\cr
    # A list containg parameter names with respective list of definitions.
    from_list = function(definitions,
                         append = TRUE) {
      if (!is.list(definitions)) {
        stop(cli::format_error("Argument {.code def_list} should be a list."))
      }
      if (isFALSE(append)) {
        self$drop_all()
      }
      for (i in seq_along(definitions)){
        this_def <- definitions[[i]]
        this_def[["name"]] <- names(definitions[i])
        do.call(self$add, this_def)
      }

      invisible(self)
    },

    # @description
    # Add missing parameter from template definitions.
    # @param refer_file
    # Path/to/.refer source file.
    # @param augm_frame
    # @param append
    template = function(refer_file,
                        augm_frame = NULL,
                        append = TRUE) {
      # retrieve the template parameters file
      templ_file <- file.path(dirname(refer_file), "template_params.R")

      # check if template file exists
      if (file.exists((templ_file))) {
        # template parameter definitions
        templ_defs <- dget(templ_file)

        if (isTRUE(append)) {
          # set of required parameters
          require_set <- self$requisites(refer_file, augm_frame)
          if (is.null(require_set)) {
            cli::cli_alert_warning(cli::format_warning(c(
              "Appending template parameter definitions failed.",
              i = "Alternatively use option {.code append = FALSE}."
            )))
          } else {
            missing_set <- template_find_missing(require_set,
                                            self$defs[["name"]])
            self$from_list(templ_defs[intersect(names(templ_defs), missing_set)],
                           append = TRUE)
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

    # @description
    # Assert parameter definitions from the implementation script.
    # @param refer_file
    # @param augm_frame
    # @return (`logical(1)`)\cr
    assert = function(refer_file,
                      augm_frame) {
      # logical state for complete definitions
      is_def_complete <- TRUE

      # set of required parameters
      require_set <- self$requisites(refer_file, augm_frame)

      if (is.null(require_set)) {
        # if requirement set is NULL it means nothing can be said
        is_def_complete <- NA

      } else {
        
        # find missing parameters
        missing_set <- template_find_missing(require_set,
                                        self$defs[["name"]])
        if (length(missing_set) > 0) {
          is_def_complete <- FALSE
          oss <- paste0("Missing parameter(s): ",
                        character_to_string_list(missing_set), ".")
          cli::cli_alert_warning(cli::format_warning(oss))
        }

        # point out parameters that are not used
        unused_set <- template_find_unused(require_set, 
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
