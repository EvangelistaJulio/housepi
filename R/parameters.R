#
# parameters.R
#
# Copyright (C) 2024 Francesco Camaglia - Institut Pasteur, Paris (France).
#

HousepiParameterClass <- R6::R6Class(
  classname = "HousepiParameterClass",
  cloneable = FALSE,
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
    # @param name (`character(1)`)\cr
    # Name of the parameter.
    # @param value (`double(1)`)\cr
    # Initial value of the parameter.
    # @param infer (`logical(1)`)\cr
    # If `FALSE` the parameter is treated as a constant.
    # @param prior (`character(1)`)\cr
    # Name of the prior distribution.
    # @param args (named `list()` of `double(1)`)\cr
    # Arguments for the prior distribution.
    # @param rule (`character(1)`)\cr
    # Name of the inference rule.
    # @param rule_args (named `list()` of `double(1)`)\cr
    # Arguments for the inference rule.
    # They can be changed automatically during iterations with adaptation.
    # @param range (named `list()` of `double(1)`)\cr
    # The prior truncation 'min', 'max' or no bound for the parameter value.
    add = function(name = NA,
                   value = NA,
                   infer = NA,
                   prior = NA,
                   args = list(),
                   rule = NA,
                   rule_args = list(),
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

      #----- extrapolate infer value -----#
      if (is.na(infer)) {
        infer <- !is.na(rule)
        if (is.na(rule) && (!is.na(prior))) {
          cli::cli_alert_warning(cli::format_warning(c(
            sprintf("'%s' interpreted as if infer = FALSE.", name),
            i = "Provide one between rule/infer to silence this warning."
          )))
        }
      }

      if (isTRUE(infer)) {
        #----- prior -----#
        if (is.na(prior)) {
          stop(cli::format_error(
            sprintf("Incomplete definition: provide prior name for '%s'.", name)
          ))
        }

        if (length(args) == 0) {
          args_ <- NA
          argsval_ <- NA
        } else {
          args_ <- str_join_(names(args))
          argsval_ <- str_join_(unname(unlist(args)))
        }

        #----- rule -----#
        if (is.na(rule)) {
          stop(cli::format_error(
            sprintf("Incomplete definition: provide rule name for '%s'.", name)
          ))
        } else if (length(rule_args) == 0) {
          rule_args_ <- NA
          rule_argsval_ <- NA
        } else {
          rule_args_ <- str_join_(names(rule_args))
          rule_argsval_ <- str_join_(unname(unlist(rule_args)))
        }

        #----- range -----#
        min_ <- get_val_(names(range), "min")
        max_ <- get_val_(names(range), "max")

        #----- self$defs : append -----#
        self$defs <- rbind(
          self$defs,
          data.frame(
            name = name,
            value = value,
            infer = TRUE,
            prior = prior,
            args = args_,
            argsval = argsval_,
            rule = rule,
            rule_args = rule_args_,
            rule_argsval = rule_argsval_,
            min = min_,
            max = max_,
            stringsAsFactors = FALSE
          )
        )

      } else {
        if (is.na(value)) {
          stop(cli::format_error(
            sprintf("Incomplete definition: missing value for constant parameter '%s'.", name)
          ))
        }
        self$defs <- rbind(
          self$defs,
          data.frame(
            name = name,
            value = value,
            infer = FALSE,
            prior = NA,
            args = NA,
            argsval = NA,
            rule = NA,
            rule_args = NA,
            rule_argsval = NA,
            min = NA,
            max = NA,
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
        infer = logical(),
        prior = character(),
        args = character(),
        argsval = character(),
        rule = character(),
        rule_args = character(),
        rule_argsval = character(),
        min = numeric(),
        max = numeric(),
        stringsAsFactors = FALSE
      )
      invisible(self)
    },

    # @description
    # Append a list of parameter definitions.
    # @param ... (named `list()`s)\cr
    # Parameter names with respective list of definitions.
    append = function(...) {
      # convert definitions to list
      defs <- list(...)

      if (length(defs) == 0) {
        # reinitialize according to templated parameters
        stop(cli::format_error(c(
          x = "No definitions were provided.",
          i = "To use template definitions use method param_template()."
        )))

      } else {
        # user definitions
        for (p in names(defs)) {
          this <- defs[[p]]
          self$add(
            name = p,
            value = get_val_(this, "value"),
            infer = get_val_(this, "infer"),
            prior = get_val_(this, "prior"),
            args = get_val_(this, "args", default = list()),
            rule = get_val_(this, "rule"),
            rule_args = get_val_(this, "rule_args", default = list()),
            range = get_val_(this, "range", default = list())
          )
        }
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
          x = "Impossible to automatically initialize parameters.",
          i = "A file with template definitionse was not found."
        )))
      }
      # remove previous definitions
      self$clear()
      # append definitions
      do.call(self$append, dget(templ_params_file))
      invisible(self)
    },

    # @description
    # Check parameter definitions from the implementation script.
    # @param fill_template (`logical(1)`)\cr
    # If `TRUE` fills the missing parameter with template definitions.
    # @param fill_template (`logical(1)`)\cr
    # If `TRUE` fills the missing parameter with template definitions.
    check = function(envir_path, augm_frame, fill_template = FALSE) {

      # retrieve the model source code
      model_path <- file.path(envir_path, ".source.txt")

      if (!file.exists(file.path(envir_path, ".source.txt"))) {
        cli::cli_alert_warning(cli::format_warning(c(
         x = "Impossible to find '.source.txt': file has been moved or deleted.",
         i = "Try to generate a new model object with {.code force = TRUE}."
        )))

      } else {

        # read the script and extract parameters
        script <- readChar(model_path, file.info(model_path)$size)

        # from : param("string" or get_param('string'
        regular_expr <- 'param\\s*\\(\\s*["\'](.*?)["\']'
        matches <- stringr::str_match_all(script, regular_expr)
        model_params_set <- unique(unlist(matches[[1]][, 2]))
        
        # param from : param_pointer("string" or param_pointer('string'
        regular_expr <- 'param_pointer\\s*\\(\\s*["\'](.*?)["\']'
        matches <- stringr::str_match_all(script, regular_expr)
        model_params_set <- unique(c(model_params_set, unlist(matches[[1]][, 2])))

        # param from : assign_contact_rate(, , "string"
        regular_expr <- 'assign_contact_rate\\s*\\(\\s*[^,]*\\s*,\\s*[^,]*\\s*,\\s*["\'](.*?)["\']'
        matches <- stringr::str_match_all(script, regular_expr)
        model_params_set <- unique(c(model_params_set, unlist(matches[[1]][, 2])))

        # param from : assign_group_rate(, , "string"
        regular_expr <- 'assign_group_rate\\s*\\(\\s*[^,]*\\s*,\\s*[^,]*\\s*,\\s*["\'](.*?)["\']'
        matches <- stringr::str_match_all(script, regular_expr)
        model_params_set <- unique(c(model_params_set, unlist(matches[[1]][, 2])))

        # param from : augmented variable model argument definitions
        # UPDATE : to be controlled
        if (length(augm_frame) > 0) {
          x <- augm_frame[c("argsval", "argstype")]
          for (i in seq_along(x)) {
            aux <- data.frame(vals = strsplit(x[i, 1], "/")[[1]],
                              keys = strsplit(x[i, 2], "/")[[1]])
            aux <- aux[aux[["keys"]] == "P", "vals"]
            model_params_set <- unique(c(model_params_set, aux[!is.na(aux)]))
          }
        }

        # re-order alphabetically
        model_params_set <- model_params_set[order(model_params_set)]

        # add missing parameters with standard options
        missing <- c()
        for (name in model_params_set) {
          if (!(name %in% self$defs[["name"]])) {
            missing <- c(missing, name)
          }
        }

        # eventually fill with template definitions
        if (length(missing) > 0) {
          if (fill_template) {
            # retrieve the template parameters file
            templ_file <- file.path(envir_path, "template_params.R")
            if (!file.exists((templ_file))) {
              stop(cli::format_error(c(
                x = "Impossible to fill missing parameters without a template.",
                i = "If you just want to check, set {.code fill_template = FALSE}."
              )))
            }
            templ_params_list <- dget(templ_file)[missing]
            do.call(self$append, templ_params_list)
            oss <- sprintf("Added missing parameter(s) from template: '%s'", missing[[1]])
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
      }

      invisible(self)
    }
  )
)