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
          cli::cli_alert_warning(cli::format_warning(c(
            sprintf("Unknown augmented variable \"%s\" to remove.", name
          ))))
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
      if(self$remove(name, status, FALSE))
      {
        if (is.na(status)) {
          cli::cli_alert_warning(cli::format_warning(c(
            sprintf(
              "Existing \"%s\" definition has been overwritten.",
              name
          ))))
        } else {
          cli::cli_alert_warning(cli::format_warning(c(
            sprintf(
              "Existing \"%s\" with status %d definition has been overwritten.", 
              name, status
          ))))
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
              value_ = str_split_(rule_from)$var
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
        this_def <- definitions[[i]]
        this_def[["name"]] <- names(definitions[i])
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
    # Overwrite augmented variable definitions with template definitions.
    template = function(envir_path) {
      # retrieve the template augmented file
      templ_augm_file <- file.path(envir_path, "template_augm.R")
      if (!file.exists((templ_augm_file))) {
        stop(cli::format_error(c(
          "Impossible to automatically initialize augmented variables.",
          i = "A file with template definitionse was not found."
        )))
      }
      # remove previous definitions
      self$from_list(dget(templ_augm_file), append = FALSE)
      invisible(self)
    },

    # @description
    # Check augmented definitions from the model implementation script.
    # @param refer_file
    # @param fill (`logical(1)`)\cr
    # If `TRUE` fills missing augmented variables with template definitions.
    # @return (`logical(1)`)\cr
    assert = function(refer_file,
                      fill = FALSE) {
      # logical state for complete definitions
      is_def_complete <- TRUE

      # check if reference file exists
      if (file.exists(refer_file)) {

        # read the script and extract augmented variables
        script <- readChar(refer_file, file.info(refer_file)$size)

        # from : augm(, "string" or augm(, 'string'
        regular_expr <- 'augm\\s*\\(\\s*[^,]*\\s*,\\s*["\'](.*?)["\']'
        matches <- stringr::str_match_all(script, regular_expr)
        model_augm_set <- unique(unlist(matches[[1]][, 2]))

        # from : set_augm(, "string"
        regular_expr <- 'set_augm\\s*\\(\\s*[^,]*\\s*,\\s*["\'](.*?)["\']'
        matches <- stringr::str_match_all(script, regular_expr)
        model_augm_set <- unique(c(model_augm_set, unique(unlist(matches[[1]][, 2]))))

        # re-order alphabetically
        model_augm_set <- model_augm_set[order(model_augm_set)]

        ####################################
        # find missing augmented variables #
        ####################################
        missing <- c()
        for (name in model_augm_set) {
          if (!(name %in% self$defs[["name"]])) {
            missing <- c(missing, name)
          }
        }

        # eventually fill with template definitions
        if (length(missing) > 0) {
          if (isTRUE(fill)) {
            # retrieve the template augmented file
            templ_augm_file <- file.path(dirname(refer_file), "template_augm.R")
            if (file.exists((templ_file))) {
              # UPDATE : I should check if the fill is complete
              self$from_list(dget(templ_augm_file)[missing], append = TRUE)
              oss <- sprintf(
                "Added missing augmented variable(s) from template: \"%s\"",
                missing[[1]]
              )
            } else {
              stop(cli::format_error(c(
                "Impossible to fill missing augmented without a template.",
                i = "Set {.code fill = FALSE} to silence this warning."
              )))              
            }

          } else {
            oss <- sprintf("Missing augmented variable(s): \"%s\"", missing[[1]])
            is_def_complete <- FALSE
          }

          # complete warning messages
          for (name in missing[-1]) {
            oss <- sprintf("%s, \"%s\"", oss, name)
          }
          cli::cli_alert_warning(cli::format_warning(c(paste0(oss, "."))))
        }

        # point out augmented variables that are not used
        unused <- c()
        for (name in self$defs[["name"]]) {
          if (!(name %in% model_augm_set)) {
            unused <- c(unused, name)
          }
        }
        if (length(unused) > 0) {
          oss <- sprintf("Augmented variable(s) not used in model: \"%s\"", unused[[1]])
          for (name in unused[-1]) {
            oss <- sprintf("%s, \"%s\"", oss, name)
          }
          cli::cli_alert_warning(cli::format_warning(c(paste0(oss, "."))))
        }

      } else {
        cli::cli_alert_warning(cli::format_warning(c(
          sprintf("Impossible to find reference file \"%s\".", refer_file),
          i = "Try to generate a new model object with {.code force = TRUE}."
        )))
        is_def_complete <- NA
      } 

      return(is_def_complete)
    }
  )
)
