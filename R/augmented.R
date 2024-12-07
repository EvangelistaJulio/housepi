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
    initialize = function() {
      self$clear()
      invisible(self)
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
                   range = list(),
                   init_rand = FALSE,
                   status = NA) {
      #----- name -----#
      if (is.na(name)) {
        stop(cli::format_error("Incomplete definition: missing name."))

      } else {

        mask <- (self$defs[["name"]] == name)
        if (any(mask)) {
          if (is.na(status)) {
            cli::cli_alert_warning(cli::format_warning(c(
              sprintf("Existing '%s' definition has been overwritten.", name)
            )))
            # drop previous definitions
            self$defs <- self$defs[!mask, ]

          } else {
            # check if same name and status
            mask <- mask & (self$defs[["status"]] == status)
            mask[is.na(mask)] <- FALSE
            if (any(mask)) {
              cli::cli_alert_warning(cli::format_warning(c(
                sprintf("Existing '%s' with status %d definition has been overwritten.", name, status)
              )))
              # drop previous definitions
              self$defs <- self$defs[!mask, ]
            }
          }
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
            "Incomplete definition: missing value.",
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
      rule_scale_ <- NA
      rule_scaletype_ <- NA
      rule_from_ <- NA
      rule_fromtype_ <- NA
      if ((!is.na(rule)) & (rule != "gibbs")) {
        rule_scale_ <- str_split_(rule_scale)[["values"]]
        rule_scaletype_ <- str_split_(rule_scale)[["types"]]

        rule_from_ <- str_split_(rule_from)[["values"]]
        rule_fromtype_ <- str_split_(rule_from)[["types"]]
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
          discrete = isTRUE(discrete),
          model = model,
          args = args_,
          argsval = argsval_,
          argstype = argstype_,
          rule = rule,
          rule_scale = rule_scale_,
          rule_scaletype = rule_scaletype_,
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
    # Associate augmented variables from a list.
    from_list = function(definitions, append = FALSE) {
      if (!is.list(definitions)) {
        stop("The argument definitions should be a list.")
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
    # Remove all current augmented variable definitions.
    clear = function() {
      self$defs <- data.frame(
        name = character(),
        status = character(),
        value = character(),
        valuetype = character(),
        discrete = logical(),
        model = character(),
        args = character(),
        argsval = character(),
        argstype = character(),
        rule = character(),
        rule_scale = character(),
        rule_scaletype = character(),
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
    # @param fill_template (`logical(1)`)\cr
    # If `TRUE` fills missing augmented variables with template definitions.
    check = function(envir_path, fill_template = FALSE) {
      # retrieve the model source code
      model_path <- file.path(envir_path, ".source")

      if (!file.exists(file.path(envir_path, ".source"))) {
        cli::cli_alert_warning(cli::format_warning(c(
          "Impossible to find '.source': file moved or deleted.",
          i = "Try to generate a new model object with {.code force = TRUE}."
        )))

      } else {

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
                "Impossible to fill missing augmented variables without a template.",
                i = "If you just want to check, set {.code fill_template = FALSE}."
              )))
            }
            self$from_list(dget(templ_augm_file)[missing], append = TRUE)
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
            oss <- sprintf("%s, '%s'", oss, name)
          }
          cli::cli_alert_warning(cli::format_warning(c(paste0(oss, "."))))
        }
      }

      invisible(self)
    }
  )
)
