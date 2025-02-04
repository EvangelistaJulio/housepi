#
# model_class.R
#
# Copyright (C) 2025 Francesco Camaglia - Institut Pasteur, Paris (France).
#

# @title Housepi Augmented R6 Class
# @description
# Help with augmented variable definitions of model Rcpp object.
# @export HousepiAugmentedClass
# @docType class
HousepiAugmentedClass <- R6::R6Class(
  classname = "HousepiAugmentedClass",
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
    # Vector of required augmented variables.
    # @param refer_file
    # Path/to/.refer source file.
    requisites = function(refer_file) {
      # vector of regular expression to use for assertion
      regular_expr <- c(
        # from : .augm(, "string" or augm(, 'string'
        '\\.augm\\s*\\(\\s*[^,]*\\s*,\\s*["\'](.*?)["\']',
        # from : .set_augm_value(, "string"
        '\\.set_augm_value\\s*\\(\\s*[^,]*\\s*,\\s*["\'](.*?)["\']'
      )
      # retrieve require augmented variables from regular expressions
      require_set <- template_find_requirement(refer_file, regular_expr)

      return(require_set)
    },

    # @description
    # Remove a given augmented variable
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
          cli::cli_alert_warning(cli::format_warning(sprintf(
            "Unknown augmented variable \"%s\" to remove.", name
          )))
        }
      }

      return(FALSE)
    },

    # @description
    # Add a new augmented variable definition.
    add = function(name = NA,
                   value = NA,
                   discrete = FALSE,
                   model = NA,
                   args = list(),
                   rule = NA,
                   rule_scale = NA,
                   rule_from = NA,
                   min = NA,
                   max = NA,
                   init_rand = FALSE,
                   status = NA) {
      #----- name -----#
      if (is.na(name)) {
        stop(cli::format_error("Incomplete definition: missing name."))
      }

      # remove previous definition if it exists
      if (self$remove(name, status, FALSE)) {
        if (is.na(status)) {
          cli::cli_alert_warning(cli::format_warning(sprintf(
            "Existing \"%s\" definition has been overwritten.",
            name
          )))
        } else {
          cli::cli_alert_warning(cli::format_warning(sprintf(
            "Existing \"%s\" with status %d definition has been overwritten.",
            name, status
          )))
        }
      }

      #----- value/init_rand -----#
      value_ <- NA
      value_type_ <- NA
      if (isTRUE(init_rand)) {
        if (is.na(model)) {
          stop(cli::format_error("Unknown model for random initialization."))
        }
      } else {
        if (is.na(value)) {
          if (is.na(model)) {
            if (is.na(rule_from)) {
              stop(cli::format_error(c(
                "Incomplete definition: missing value.",
                i = "Alternative: provide a model and set {.code init_rand = TRUE}."
              )))
            } else {
              value_ <- str_split_(rule_from)$var
              value_type_ <- str_split_(rule_from)$typ
            }
          } else {
            init_rand <- TRUE
          }
        } else {
          value_ <- str_split_(value)$var
          value_type_ <- str_split_(value)$typ
        }
      }

      #----- method -----#
      arg_name_ <- NA
      arg_value_ <- NA
      arg_type_ <- NA
      # assign method only if arguments are provided
      if ((!is.na(model)) && (length(args) > 0)) {
        arg_name_ <- str_join_(names(args))
        tmp <- str_split_(args)
        arg_value_ <- str_join_(tmp$var)
        arg_type_ <- str_join_(tmp$typ)
      }

      #----- rule -----#
      rule_scale_ <- NA
      rule_scale_type_ <- NA
      rule_from_ <- NA
      rule_from_type_ <- NA

      if (is.na(rule)) {
        rule <- "gibbs"
      } else if (rule != "gibbs") {
        if (is.na(rule_scale)) {
          rule_scale_ <- 1
          rule_scale_type_ <- "S"
        } else {
          rule_scale_ <- str_split_(rule_scale)$var
          rule_scale_type_ <- str_split_(rule_scale)$typ
        }

        rule_from_ <- str_split_(rule_from)$var
        rule_from_type_ <- str_split_(rule_from)$typ
      }

      #----- minimum value -----#
      min_ <- NA
      min_type_ <- NA
      if (!is.na(min)) {
        min_ <- str_split_(min)$var
        min_type_ <- str_split_(min)$typ
      }

      #----- maximum value -----#
      max_ <- NA
      max_type_ <- NA
      if (!is.na(max)) {
        max_ <- str_split_(max)$var
        max_type_ <- str_split_(max)$typ
      }

      #----- append the new definitions -----#
      self$defs <- rbind(
        self$defs,
        data.frame(
          name = name,
          status = status,
          value = value_,
          value_type = value_type_,
          discrete = isTRUE(discrete),
          model = model,
          arg_name = arg_name_,
          arg_value = arg_value_,
          arg_type = arg_type_,
          rule = rule,
          rule_scale = rule_scale_,
          rule_scale_type = rule_scale_type_,
          rule_from = rule_from_,
          rule_from_type = rule_from_type_,
          min = min_,
          min_type = min_type_,
          max = max_,
          max_type = max_type_,
          init_rand = init_rand,
          stringsAsFactors = FALSE
        )
      )
      invisible(self)
    },

    # @description
    # Associate augmented variables from a list.
    from_list = function(definitions,
                         append = TRUE) {
      if (!is.list(definitions)) {
        stop("The argument definitions should be a list.")
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
              sprintf("Wrong definition for augmented variable \"%s\"", def_name),
              i = "Provide a definition list; if constant, a scalar or a typed value."
            )))
          }
        }
        this_def[["name"]] <- def_name
        do.call(self$add, this_def)
      }

      invisible(self)
    },

    # @description
    # Drop all current augmented variable definitions.
    drop_all = function() {
      self$defs <- data.frame(
        name = character(),
        status = character(),
        value = character(),
        value_type = character(),
        discrete = logical(),
        model = character(),
        arg_name = character(),
        arg_value = character(),
        arg_type = character(),
        rule = character(),
        rule_scale = character(),
        rule_scale_type = character(),
        rule_from = character(),
        rule_from_type = character(),
        min = character(),
        min_type = character(),
        max = character(),
        max_type = character(),
        init_rand = logical(),
        stringsAsFactors = FALSE
      )
      invisible(self)
    },

    # @description
    # Add missing augmented variables from template definitions.
    # @param refer_file
    # Path/to/.refer source file.
    # @param augm_frame
    # @param append
    template = function(refer_file,
                        append = TRUE) {
      # retrieve the template augmented file
      templ_file <- file.path(dirname(refer_file), "template_augm.R")

      # check if template file exists
      if (file.exists((templ_file))) {
        # template augmented variable definitions
        templ_defs <- dget(templ_file)

        if (isTRUE(append)) {
          # set of required augmented variables
          require_set <- self$requisites(refer_file)
          if (is.null(require_set)) {
            cli::cli_alert_warning(cli::format_warning(c(
              "Appending template augmented variable definitions failed.",
              i = "Alternatively use option {.code append = FALSE}."
            )))
          } else {
            missing_set <- template_find_missing(require_set,
                                                 self$defs[["name"]])
            self$from_list(templ_defs[missing_set],
                           append = TRUE)
          }
        } else {
          self$from_list(templ_defs, append = FALSE)
        }

      } else {
        stop(cli::format_error(c(
          "Impossible to automatically initialize augmented variables.",
          i = "A file with template definitionse was not found."
        )))
      }

      invisible(self)
    },

    # @description
    # Check augmented definitions from the model implementation script.
    # @param refer_file
    # @return (`logical(1)`)\cr
    assert = function(refer_file) {
      # logical state for complete definitions
      is_def_complete <- TRUE

      # set of required augmented variables
      require_set <- self$requisites(refer_file)

      if (is.null(require_set)) {
        # if requirement set is NULL it means nothing can be said
        is_def_complete <- NA

      } else {

        # find missing augmented variables
        missing_set <- template_find_missing(require_set,
                                             self$defs[["name"]])
        if (length(missing_set) > 0) {
          is_def_complete <- FALSE
          oss <- paste0("Missing parameter(s): ",
                        character_to_string_list(missing_set), ".")
          cli::cli_alert_warning(cli::format_warning(oss))
        }

        # point out augmented variables that are not used
        unused_set <- template_find_unused(require_set,
                                           self$defs[["name"]])

        # drop unused augmented variables
        if (length(unused_set) > 0) {
          oss <- paste0("Dropping not used augmented variable(s): ",
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
