#
# model_class.R
#
# Copyright (C) 2024 Francesco Camaglia - Institut Pasteur, Paris (France).
#

HousepiAugmentedClass <- R6::R6Class(
  classname = "HousepiAugmentedClass",
  cloneable = FALSE,

  public = list(
    defs = NULL,

    # @description
    # Assign the empty parameter data.frame
    initialize = function(){
      self$clear()
      invisible(self)
    },

    # @description Add a new augmented variable definition.
    # @param name (`character(1)` )\cr
    # Name of the augmented variable.
    # @param value (`double(1)` | `character(1)`)\cr
    # Scalar value
    # or parameter name (followed by `@P`)
    # or data name (followed by `@D`).
    # @param model (`character(1)` )\cr
    # Name of the model distribution.
    # @param args - List of arguments for the prior distribution.
    # @param rule (`character(1)` )\cr
    # Name of the inference rule.
    # @param rule_args (named `list()` of `double(1)` | `character(1)`)\cr
    # Arguments for the inference rule.
    # @param rule_from (`double(1)` | `character(1)`)\cr
    # Advanced options for non random walk moves.
    # @param range - List with 'min' and/or 'max' values (or none).
    # @param init_rand (`character(1)` )\cr
    # If `TRUE` draw initial vale from model distribution.
    # @param status (`integer(1)` )\cr
    # The associated status of the definition. Default NA is for all infectious.
    add = function(name = NA,
                   value = NA,
                   model = NA,
                   args = list(),
                   rule = NA,
                   rule_args = list(),
                   rule_from = NA,
                   range = list(),
                   init_rand = FALSE,
                   status = NA) {
      #----- name -----#
      if (is.na(name)) {
        stop(cli::format_error("Incomplete definition: missing name."))

      } else {
        if (name %in% self$defs[["name"]]) {
          cli::cli_alert_warning(cli::format_warning(c(
            sprintf("Existing '%s' definition has been overwritten.", name)
          )))
          mask <- self$defs[["name"]] != name
          self$defs <- self$defs[mask, ]
        }
      }

      #----- value/init_rand -----#
      if (isTRUE(init_rand)) {
        if (is.na(model)) {
          stop(cli::format_error("Unknown model for random initialization."))
        } else {
          value_ <- NA
          valuetype_ <- NA
        }
      } else {
        if (is.na(value)) {
          stop(cli::format_error(c(
              x = "Incomplete definition: missing value.",
              i = "Alternatively provide a model and set {.code init_rand = TRUE}."
          )))
        } else {
          value_ <- str_split_(value)[["values"]]
          valuetype_ <- str_split_(value)[["types"]]
        }
      }

      #----- method -----#
      args_ <- NA
      argsval_ <- NA
      argstype_ <- NA
      # assign method only if arguments are provided
      if (!is.na(model))  {
        if (length(args) == 0) {
          cli::cli_alert_warning(cli::format_warning(c(
            sprintf("Missing args for augmented variable '%s'.", name)
          )))
        } else {
          args_ <- str_join_(names(args))
          tmp <- str_split_(args)
          argsval_ <- str_join_(tmp[["values"]])
          argstype_ <- str_join_(tmp[["types"]])
        }
      }

      #----- rule -----#
      rule_args_ <- NA
      rule_argsval_ <- NA
      rule_argstype_ <- NA
      rule_from_ <- NA
      rule_fromtype_ <- NA
      if (!is.na(rule)) {
        if ((rule != "gibbs") && (length(rule_args) == 0)) {
          cli::cli_alert_warning(cli::format_warning(c(
            sprintf("Missing rule_args for augmented variable '%s'.", name)
          )))
        } else {
          rule_args_ <- str_join_(names(rule_args))
          tmp <- str_split_(rule_args)
          rule_argsval_ <- str_join_(tmp[["values"]])
          rule_argstype_ <- str_join_(tmp[["types"]])
          rule_from_ <- str_split_(rule_from)[["values"]]
          rule_fromtype_ <- str_split_(rule_from)[["types"]]
        }
      }

      #----- range -----#
      if (length(range) == 0) {
        # empty range
        min_ <- NA
        mintype_ <- NA
        max_ <- NA
        maxtype_ <- NA
      } else {
        # minimum value
        if ("min" %in% names(range)) {
          tmp <- str_split_(range[names(range) == "min"])
          min_ <- tmp[["values"]]
          mintype_ <- tmp[["types"]]
        } else {
          min_ <- NA
          mintype_ <- NA
        }
        # maximum value
        if ("max" %in% names(range)) {
          tmp <- str_split_(range[names(range) == "max"])
          max_ <- tmp[["values"]]
          maxtype_ <- tmp[["types"]]
        } else {
          max_ <- NA
          maxtype_ <- NA
        }
      }

      #----- append the new definitions -----#
      self$defs <- rbind(
        self$defs,
        data.frame(
          name = name,
          status = status,
          value = value_,
          valuetype = valuetype_,
          model = model,
          args = args_,
          argsval = argsval_,
          argstype = argstype_,
          rule = rule,
          rule_args = rule_args_,
          rule_argsval = rule_argsval_,
          rule_argstype = rule_argstype_,
          rule_from = rule_from_,
          rule_fromtype = rule_fromtype_,
          min = min_,
          mintype = mintype_,
          max = max_,
          maxtype = maxtype_,
          init_rand = init_rand,
          stringsAsFactors = FALSE
        )
      )
      invisible(self)
    },

    # @description
    # Remove all current augmented variable definitions.
    clear = function() {
      self$defs <- data.frame(
        name = character(),
        status = character(),
        value = character(),
        valuetype = character(),
        model = character(),
        args = character(),
        argsval = character(),
        argstype = character(),
        rule = character(),
        rule_args = character(),
        rule_argsval = character(),
        rule_argstype = character(),
        rule_from = character(),
        rule_fromtype = character(),
        min = character(),
        mintype = character(),
        max = character(),
        maxtype = character(),
        init_rand = logical(),
        stringsAsFactors = FALSE
      )
      invisible(self)
    },

    # @description
    # Append a data frame of augmented variable definitions.
    # @param ... (named `list()`s)\cr
    # Augmented variable names with respective list of definitions.
    append = function(...) {
      # convert definitions to list
      definitions <- list(...)

      if (length(definitions) == 0) {
        # reinitialize according to templated parameters
        stop(cli::format_error(c(
          x = "No definitions were provided.",
          i = "To use template definitions use method augm_template()."
        )))

      } else {
        # add each entry of the list
        for (a in names(definitions)) {
          this <- definitions[[a]]
          self$add(
            name = a,
            status = get_val_(this, "status"),
            value = get_val_(this, "value"),
            model = get_val_(this, "model"),
            args = get_val_(this, "args", default = list()),
            rule = get_val_(this, "rule"),
            rule_args = get_val_(this, "rule_args", default = list()),
            rule_from = get_val_(this, "rule_from"),
            range = get_val_(this, "range", default = list()),
            init_rand = get_val_(this, "init_rand", default = FALSE)
          )
        }
      }
      invisible(self)
    },

    # @description
    # Overwrite augmented variable definitions with template definitions.
    template = function(envir_path) {
      # retrieve the template augmented file
      templ_augm_file <- file.path(envir_path, "template_augm.R")
      if (!file.exists((templ_augm_file))) {
        stop(cli::format_error(c(
          x = "Impossible to automatically initialize augmented variables.",
          i = "A file with template definitionse was not found."
        )))
      }
      # remove previous definitions
      self$clear()
      # append definitions
      do.call(self$append, dget(templ_augm_file))
      invisible(self)
    },

    # @description
    # Check augmented definitions from the model implementation script.
    # @param fill_template (`logical(1)`)\cr
    # If `TRUE` fills missing augmented variables with template definitions.
    check = function(envir_path, fill_template = FALSE) {

      # retrieve the model script
      if (file.exists(file.path(envir_path, "model.cpp"))) {
        model_path <- file.path(envir_path, "model.cpp")
      } else if (file.exists(file.path(envir_path, "src", "model.cpp"))) {
        model_path <- file.path(envir_path, "src", "model.cpp")
      } else {
        stop(cli::format_error(c(
         x = "Impossible to find the source code 'model.cpp'.",
         i = "Check if the original file has been moved or deleted."
        )))
      }

      # read the script and extract augmented variables
      script <- readChar(model_path, file.info(model_path)$size)

      # from : augm(, "string" or augm(, 'string'
      regular_expr <- 'augm\\s*\\(\\s*[^,]*\\s*,\\s*["\'](.*?)["\']'
      matches <- stringr::str_match_all(script, regular_expr)
      model_augm_set <- unique(unlist(matches[[1]][, 2]))

      # from : set_augm(, "string"
      regular_expr <- 'set_augm\\s*\\(\\s*[^,]*\\s*,\\s*["\'](.*?)["\']'
      matches <- stringr::str_match_all(script, regular_expr)
      model_augm_set <- unique(c(model_augm_set, unique(unlist(matches[[1]][, 2]))))

      # UPDATE : import from infectivity profile?

      # re-order alphabetically
      model_augm_set <- model_augm_set[order(model_augm_set)]

      # find missing augmented variables
      missing <- c()
      for (name in model_augm_set) {
        if (!(name %in% self$defs[["name"]])) {
          missing <- c(missing, name)
        }
      }

      # eventually fill with template definitions
      if (length(missing) > 0) {
        if (isTRUE(fill_template)) {
          # retrieve the template augmented file
          templ_augm_file <- file.path(envir_path, "template_augm.R")
          if (!file.exists((templ_augm_file))) {
            stop(cli::format_error(c(
              x = "Impossible to fill missing augmented variables without a template.",
              i = "If you just want to check, set {.code fill_template = FALSE}."
            )))
          }
          templ_augm_list <- dget(templ_augm_file)[missing]
          do.call(self$append, templ_augm_list)
          oss <- sprintf("Added missing augmented variable(s) from template: '%s'", missing[[1]])

        } else {
          oss <- sprintf("Missing augmented variable(s): '%s'", missing[[1]])
        }

        for (name in missing[-1]) {
          oss <- sprintf("%s, '%s'", oss, name)
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
        oss <- sprintf("Augmented variable(s) not used in model: '%s'", unused[[1]])
        for (name in unused[-1]) {
          oss <- paste0("%s, '%s'", oss, name)
        }
        cli::cli_alert_warning(cli::format_warning(c(paste0(oss, "."))))
      }

      invisible(self)
    }
  )
)