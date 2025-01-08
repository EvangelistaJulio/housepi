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
      self$clear()
      invisible(self)
    },

    # @description
    # Add a new parameter definition.
    add = function(name = NA,
                   value = NA,
                   constant = NA,
                   discrete = FALSE,
                   prior = NA,
                   args = list(),
                   rule = NA,
                   rule_scale = NA,
                   range = list()) {
      if (is.na(name)) {
        stop(cli::format_error("Incomplete definition: missing name."))
      }

      if (name %in% self$defs[["name"]]) {
        cli::cli_alert_warning(cli::format_warning(c(
          sprintf("Overwriting existing '%s' definition.", name)
        )))
        mask <- self$defs[["name"]] != name
        self$defs <- self$defs[mask, ]
      }

      # assume it's a constant if a prior is not provided
      if (is.na(constant)) {
        constant <- is.na(prior)
      }

      if (isTRUE(constant)) {
        # check that value is provided
        if (is.na(value)) {
          stop(cli::format_error(
            sprintf("Missing value for constant parameter '%s'.", name)
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
        #----- prior -----#
        if (is.na(prior)) {
          stop(cli::format_error(
            sprintf("Incomplete definition: provide prior name for '%s'.", name)
          ))
        }

        if (length(args) == 0) {
          arg_name_ <- NA
          arg_value_ <- NA
        } else {
          arg_name_ <- str_join_(names(args))
          arg_value_ <- str_join_(unname(unlist(args)))
        }

        #----- rule -----#
        if (is.na(rule)) {
          # default proposer rule
          rule <- "gibbs"
          rule_scale <- NA
        }

        #----- range -----#
        min_ <- get_val_(names(range), "min")
        max_ <- get_val_(names(range), "max")

        #----- self$defs : add row -----#
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
            min = min_,
            max = max_,
            stringsAsFactors = FALSE
          )
        )
      }

      invisible(self)
    },

    # @description
    # Clear all current parameter definitions.
    clear = function() {
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
    from_list = function(definitions, append = FALSE) {
      if (!is.list(definitions)) {
        stop("The argument def_list should be a list.")
      }
      if (!append) {
        self$clear()
      }
      for (i in seq_along(definitions)){
        this_def <- definitions[[i]]
        this_def[["name"]] <- names(definitions[i])
        do.call(self$add, this_def)
      }

      invisible(self)
    },

    # @description
    # Overwrite current parameter definitions with template definitions.
    template = function(envir_path) {
      # retrieve the template parameters file
      templ_params_file <- file.path(envir_path, "template_params.R")
      if (!file.exists((templ_params_file))) {
        stop(cli::format_error(c(
          "Impossible to automatically initialize parameters.",
          i = "A file with template definitionse was not found."
        )))
      }
      # remove previous definitions
      self$from_list(dget(templ_params_file), append = FALSE)
      invisible(self)
    },

    # @description
    # Assert parameter definitions from the implementation script.
    # @param refer_file
    # @param augm_frame
    # @param fill (`logical(1)`)\cr
    # If `TRUE` fills the missing parameter with template definitions.
    assert = function(refer_file,
                      augm_frame,
                      fill = FALSE) {

      # read the script and extract parameters
      script <- readChar(refer_file, file.info(refer_file)$size)

      # from : param("string" or get_param('string'
      regular_expr <- 'param\\s*\\(\\s*["\'](.*?)["\']'
      matches <- stringr::str_match_all(script, regular_expr)
      model_params_set <- unlist(matches[[1]][, 2])

      # param from : param_pointer("string" or param_pointer('string'
      regular_expr <- 'param_pointer\\s*\\(\\s*["\'](.*?)["\']'
      matches <- stringr::str_match_all(script, regular_expr)
      model_params_set <- c(model_params_set, unlist(matches[[1]][, 2]))

      # param from : assign_contact_rate("string"
      regular_expr <- 'emplace_contact_rate\\s*\\(\\s*["\'](.*?)["\']'
      matches <- stringr::str_match_all(script, regular_expr)
      model_params_set <- c(model_params_set, unlist(matches[[1]][, 2]))

      # param from : vec_parID.emplace_back(..., "string"
      # UPDATE: need for smtg more solid here
      regular_expr <- 'emplace_group_rate\\s*\\(\\s*[^,]*\\s*,\\s*["\'](.*?)["\']'
      matches <- stringr::str_match_all(script, regular_expr)
      model_params_set <- c(model_params_set, unlist(matches[[1]][, 2]))

      # param from : augmented variable model argument definitions
      # UPDATE : to be controlled
      if (length(augm_frame) > 0) {
        x <- augm_frame[c("arg_value", "arg_type")]
        for (i in seq_along(x)) {
          aux <- data.frame(vals = strsplit(x[i, 1], "/")[[1]],
                            keys = strsplit(x[i, 2], "/")[[1]])
          aux <- aux[aux[["keys"]] == "P", "vals"]
          model_params_set <- c(model_params_set, aux[!is.na(aux)])
        }
      }

      # re-order alphabetically
      model_params_set <- unique(model_params_set[order(model_params_set)])

      # add missing parameters with standard options
      missing <- c()
      for (name in model_params_set) {
        if (!(name %in% self$defs[["name"]])) {
          missing <- c(missing, name)
        }
      }

      # eventually fill with template definitions
      if (length(missing) > 0) {
        if (isTRUE(fill)) {
          # retrieve the template parameters file
          templ_augm_file <- file.path(dirname(refer_file), "template_augm.R")
          if (!file.exists((templ_file))) {
            stop(cli::format_error(c(
              "Impossible to fill missing parameters without a template.",
              i = "Set {.code fill = FALSE} to silence this warning."
            )))
          }
          self$from_list(dget(templ_file)[missing], append = TRUE)
          oss <- sprintf(
            "Added missing parameter(s) from template: '%s'",
            missing[[1]]
          )

        } else {
          oss <- sprintf("Missing parameter(s): '%s'", missing[[1]])
        }

        for (name in missing[-1]) {
          oss <- sprintf("%s, '%s'", oss, name)
        }
        cli::cli_alert_warning(cli::format_warning(c(paste0(oss, "."))))
      }

      # point out parameters that are not used
      unused <- c()
      for (name in self$defs[["name"]]) {
        if (!(name %in% model_params_set)) {
          unused <- c(unused, name)
        }
      }
      if (length(unused) > 0) {
        oss <- paste0("Parameter(s) not used in model: '", unused[[1]], "'")
        for (name in unused[-1]) {
          oss <- paste0(oss, ", '", name, "'")
        }
        cli::cli_alert_warning(cli::format_warning(c(paste0(oss, "."))))
      }

      invisible(self)
    }
  )
)
