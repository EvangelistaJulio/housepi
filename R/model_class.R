#
# model_class.R
#
# Copyright (C) 2024 Francesco Camaglia - Institut Pasteur, Paris (France).
#

#' @title Household Model Class
#' @description
#' Compile household model, manage your data and parameters from a single class.
#' @export HouseholdModelClass
#' @docType class
HouseholdModelClass <- R6::R6Class(
  classname = "HouseholdModelClass",
  cloneable = FALSE,
  private = list(
    setup = NULL,
    param_frame = NULL,
    augm_frame = NULL,
    traces_frame = NULL,
    obj = NULL,
    envir_path = NA,
    model_path = NA,

    # @description
    # Check that traces exists, otherwise raise error
    traces_check = function() {
      if (is.null(private$traces_frame)) {
        stop(cli::format_error(c(
          "No MCMC traces were found.",
          i = "Be sure to run your inference first."
        )))
      }
    },

    # @description
    # Associate data and definitions to the model.
    # @param data (`character(1)` | `data.frame`)\cr
    # Path to a csv data file or a data.frame.
    # @param verbose (`logical(1)`)\cr
    # If `TRUE` prints automatic build information.
    build = function(data = NULL, verbose = FALSE) {
      # change scientific notation
      store_scipen <- options()$scipen
      options(scipen = 500)

      #----- 1) SETUP -----#
      char_setup <- lapply(private$setup[!is.na(private$setup)], as.character)
      private$obj$import_setup(char_setup)

      #----- 2) PARAMETERS -----#
      char_params <- data.frame(lapply(private$param_frame, as.character))
      private$obj$import_params(char_params)

      #----- 3) DATA -----#
      if (is.null(data)) {
        private$obj$read_data("<NA>")
      } else if (is.character(data)) {
        private$obj$read_data(data)
      } else if (is.data.frame(data)) {
        char_data <- data.frame(lapply(data, as.character))
        private$obj$import_data(char_data)
      } else {
        stop(cli::format_error(c(
          "Wrong type for data",
          i = "It should be data.frame or path/to/file."
        )))
      }

      #----- 4) AUGMENTED VARIABLES -----#
      char_augm <- data.frame(lapply(private$augm_frame, as.character))
      private$obj$import_augm(char_augm)

      # restore scientific notation option
      options(scipen = store_scipen)
    },

    # @description
    # Setup the list of configurations for the different analysis.
    # @param seed (`integer(1)`)\cr
    # Seed of the random number generator.
    # @param num_iter (`integer(1)`)\cr
    # Number of iterations.
    # @param burn_in (`integer(1)`)\cr
    # Number of iterations to discard.
    # @param thinning (`integer(1)`)\cr
    # Name of the parameter.
    # @param num_chains (`integer(1)`)\cr
    # Name of the parameter.
    # @param num_adapt (`integer(1)`)\cr
    # Name of the parameter.
    # @param num_threads (`integer(1)`)\cr
    # Name of the parameter.
    # @param rnd_init_params (`logical(1)`)\cr
    # Name of the parameter.
    # @param inv_temper (`double(1)`)\cr
    # Name of the parameter.
    setup_list = function(seed = NA,
                          num_iter = 1000,
                          burn_in = 0,
                          thinning = 1,
                          num_chains = 1,
                          num_adapt = 0,
                          num_threads = 1,
                          rnd_init_params = TRUE,
                          inv_temper = 1.0) {
      private$setup[["seed"]] <- seed
      private$setup[["num_iter"]] <- num_iter
      private$setup[["burn_in"]] <- burn_in
      private$setup[["thinning"]] <- thinning
      private$setup[["num_chains"]] <- num_chains
      private$setup[["num_adapt"]] <- num_adapt
      private$setup[["num_threads"]] <- num_threads
      private$setup[["rnd_init_params"]] <- rnd_init_params
      private$setup[["inv_temper"]] <- inv_temper
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
    param_add = function(name = NA,
                         value = NA,
                         infer = NA,
                         prior = NA,
                         args = list(),
                         rule = NA,
                         rule_args = list(),
                         range = list()) {
      if (is.na(name)) {
        stop(cli::format_error(c("Missing parameter name.")))
      }

      if (name %in% private$param_frame[["name"]]) {
        cli::cli_alert_warning(cli::format_warning(c(
          paste0("Overwriting existing '", name, "' definition.")
        )))
        mask <- private$param_frame$name != name
        private$param_frame <- private$param_frame[mask, ]
      }

      #----- extrapolate infer value -----#
      if (is.na(infer)) {
        infer <- !is.na(rule)
        if (is.na(rule) && (!is.na(prior))) {
          cli::cli_alert_warning(cli::format_warning(
            c(paste0(
              "'", name, "' interpreted as if infer = FALSE."
            ), i = "Provide one between rule/infer to silence this warning.")
          ))
        }
      }

      if (isTRUE(infer)) {
        #----- prior -----#
        if (is.na(prior)) {
          stop(cli::format_error(c(
            paste0("Provide a prior name for '", name, "'.")
          )))
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
          stop(cli::format_error(c(
            paste0("Provide a rule name for '", name, "'.")
          )))
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

        #----- private$param_frame : append -----#
        private$param_frame <- rbind(
          private$param_frame,
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
            paste0("Missing value for constant parameter '", name, "'.")
          ))
        }
        private$param_frame <- rbind(
          private$param_frame,
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
    augm_add = function(name = NA,
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
        if (name %in% private$augm_frame[["name"]]) {
          cli::cli_alert_warning(cli::format_warning(c(
            paste0("Existing '", name, "' definition has been overwritten.")
          )))
          mask <- private$augm_frame$name != name
          private$augm_frame <- private$augm_frame[mask, ]
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
          stop(cli::format_error(
            c(
              "Incomplete definition: missing value.",
              i = "Alternatively provide a model and set init_rand = TRUE"
            )
          ))
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
            paste0("Missing args for augmented variable '", name, "'.")
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
            paste0("Missing rule_args for augmented variable '", name, "'.")
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
      private$augm_frame <- rbind(
        private$augm_frame,
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
    }
  ),

  ############################################################################
  ############################################################################
  ############################################################################
  ############################################################################
  ############################################################################

  public = list(
    #' @description
    #' Creates a new instance of class [HouseholdModelClass][housepi::HouseholdModelClass].
    #' If not specified otherwise, it create a one per session model environment.
    #' @param model (`character(1)`)\cr
    #' Template name, path to folder with implementation or with previously compiled model.
    #' @param store_path (`character(1)`)\cr
    #' Path to folder where to store the compiled model environment.
    #' @param quiet (`logical(1)`)\cr
    #' If `TRUE` hide compilation commands.
    #' @param force (`logical(1)`)\cr
    #' If `TRUE` compile even if a cached model already exists.
    initialize = function(model = ".",
                          store_path = NA,
                          quiet = TRUE,
                          force = FALSE) {

      dll <- generate_model(model, store_path = store_path, 
                            quiet = quiet, force = force)
      # load shared library
      dyn.load(dll)
      # expose internally the module
      m <- Rcpp::Module("model_housepi_rcpp", PACKAGE = "housepienvironment")

      #----- assign fields -----#
      # assign the path to model environment source
      private$envir_path <- dirname(dll)
      # assign the model object
      private$obj <- m$ModelHousepi$new()

      # assign the empty setup list
      private$setup <- list()
      # assign the empty parameter data.frame
      self$param_clear()
      # assign the empty augmented data.frame
      self$augm_clear()

      invisible(self)
    },

    #' @description
    #' Show model description in RStudio Help window.
    help = function() {
      templ_info_file <- file.path(private$envir_path, "template_info.Rd")
      if (!file.exists((templ_info_file))) {
        stop(cli::format_error(c(
          "Impossible to show model description.", 
          i = "No template info were found."
        )))
      }
      rstudioapi::previewRd(templ_info_file)
      invisible(self)
    },

    ############################################################################
    ############################################################################
    ############################################################################
    ############################################################################
    ############################################################################

    #' @description
    #' Estimate the loglikelihood for the given parameters.
    #' @param data (`character(1)` | `data.frame`)\cr
    #' Path to a csv data file or a data.frame.
    #' @param ...
    #' further arguments passed to \code{setup}.
    #' @return (`double(1)`)\cr
    #' The value of the loglikelihood.
    loglikelihood = function(data = NULL, ...) {
      # assign setup
      do.call(private$setup_list, list(...))
      # internal build function
      private$build(data = data)
      # evaluate marginal from the internal Rcpp object
      return(private$obj$evaluate())
    },

    #' @description
    #' Run Bayesian MCMC inference and store internally the chain traces.
    #' @param data (`character(1)` | `data.frame`)\cr
    #' Path to a csv data file or a data.frame.
    #' @param progbar (`logical(1)`)\cr
    #' If `FALSE` suppress the inference progress bar.
    #' @param ... further arguments passed to \code{setup}.
    infer = function(data = NULL, progbar = TRUE, ...) {
      # assign setup
      do.call(private$setup_list, list(...))
      # internal build function
      private$build(data = data)
      # infer method from cpp model object
      private$obj$infer(progbar)
      # load traces from last inference
      df <- private$obj$export_traces()
      # extract setup information for coda
      load_setup <- private$obj$export_setup()
      # Check if multiple chains
      if ("CHAIN__" %in% colnames(df)) {
        df_list <- lapply(split(df, df[["CHAIN__"]]), function(x) {
          coda::mcmc(
            data = x[, !names(x) %in% "CHAIN__"],
            start = load_setup[["burn_in"]] + load_setup[["thinning"]],
            end = load_setup[["num_iter"]],
            thin = load_setup[["thinning"]]
          )
        })
        private$traces_frame <- coda::as.mcmc.list(df_list)
      } else {
        private$traces_frame <- coda::mcmc(
          data = df,
          start = load_setup[["burn_in"]] + load_setup[["thinning"]],
          end = load_setup[["num_iter"]],
          thin = load_setup[["thinning"]]
        )
      }
      invisible(self)
    },

    #' @description
    #' Export the traces of last inference call as a coda object.
    #' @return
    #' (`coda::mcmc` | `coda::mcmc.list`)
    traces_export = function() {
      private$traces_check()
      return(private$traces_frame)
    },

    #' @description
    #' Export the correlations between traces values of last inference call.
    #' @return
    #' A (`matrix` | 3d `array`)
    traces_crosscorr = function() {
      private$traces_check()
      return(coda::crosscorr(private$traces_frame))
    },

    #' @description
    #' Export the traces of last inference call as a coda object.
    #' @return
    #' (`data.frame`)
    traces_summary = function() {
      private$traces_check()
      res <- summary(private$traces_frame)
      imask <- private$param_frame["infer"] == TRUE
      inames <- private$param_frame[imask, "name"]
      output <- data.frame(
        Mean = res[["statistics"]][inames, "Mean"],
        SD = res[["statistics"]][inames, "SD"],
        Q2.5 = res[["quantiles"]][inames, "2.5%"],
        Q97.5 = res[["quantiles"]][inames, "97.5%"]
      )
      return(output)
    },

    #' @description
    #' Calls `coda` traceplot function to plot traces of a single chain.
    #' @param chain (`integer(1)`)
    #' Option to specifiy which chain to plot if multiple chains.
    #' @param col (`character(1)`)
    #' Name of the color or RGB hex value.
    #' @param ...
    #' Further arguments pass to [coda][coda::traceplot].
    traces_traceplot = function(chain = NA,
                                col = "dodgerblue",
                                ...) {
      private$traces_check()
      if (is.null(private$traces_frame)) {
        stop("No traces found.")
      } else if (coda::is.mcmc.list(private$traces_frame)) {
        if (is.na(chain)) {
          coda::traceplot(private$traces_frame, col = col, ...)
        } else {
          coda::traceplot(private$traces_frame[chain], col = col, ...)
        }
      } else {
        coda::traceplot(private$traces_frame, col = col, ...)
      }

      invisible(self)
    },

    #' @description
    #' Run the simulation.
    #' @param data (`character(1)` | `data.frame`)\cr
    #' Path to a csv data file or a data.frame.
    #' @param progbar (`logical(1)`)\cr
    #' If `FALSE` suppress the simulation progress bar.
    #' @param ... further arguments passed to \code{setup}.
    #' Named build options.
    simulate = function(data = NULL, progbar = TRUE, ...) {
      # assign setup
      do.call(private$setup_list, list(...))
      # internal build function
      private$build(data = data)
      # simulate method from internal Rcpp object
      private$obj$simulate(progbar)
      invisible(self)
    },

    #' @description
    #' Export the loaded database at current state.
    #' If the case, it also return the associated augmented variables.
    #' @return (`data.frame`)\cr
    #' A data.frame with individual data.
    data_export = function() {
      return(private$obj$export_data())
    },

    #' @description
    #' Print a summary of the loaded database.
    data_summary = function() {
      return(private$obj$show_summary())
    },

    ############################################################################
    ############################################################################
    ############################################################################
    ############################################################################
    ############################################################################

    #' @description
    #' Use loo package to compute the leave-one-out cross-validation.
    #' @param ... arguments passed to [loo][loo::loo].
    #' For arguments and further details type ?loo::loo.
    loo = function(...) {
      return(loo::loo(private$obj$export_loglike_matrix(), ...))
    },

    #' @description
    #' Use loo package to compute the widely applicable information criterion.
    #' @param ... arguments passed to [loo][loo::waic].
    #' #' For arguments and further details type ?loo::waic.
    waic = function(...) {
      return(loo::waic(private$obj$export_loglike_matrix(), ...))
    },

    #' @description
    #' Estimate the model evidence (or marginal likelihood).
    #' It uses the method of power posteriors.
    #' @param data (`character(1)` | `data.frame`)\cr
    #' Path to a csv data file or a data.frame.
    #' @param num_rungs (`integer`)\cr
    #' The number of rungs in the thermodynamic ladder.
    #' @param adaptive (`logical(1)`)\cr
    #' If `TRUE` uses adaptive scheme for thermodynamic ladder.
    #' Else use power fractions.
    #' @param progbar (`logical(1)`)\cr
    #' If `FALSE` suppress the simulation progress bar.
    #' @return (`double(1)`)\cr
    #' The value of the model evidence.
    #' @param ... further arguments passed to  \code{setup}.
    evidence = function(data = NULL,
                        num_rungs = 40,
                        adaptive = TRUE,
                        progbar = TRUE,
                        ...) {
      # assign setup
      do.call(private$setup_list, list(...))
      # internal build function
      private$build(data = data)
      # evaluate marginal from the internal Rcpp object
      output <- private$obj$evaluate_marginal(num_rungs, adaptive, progbar)
      # drop traces
      private$traces_frame <- NULL
      return(output)
    },

    ############################################################################
    ############################################################################
    ############################################################################
    ############################################################################
    ############################################################################

    #' @description
    #' Export the parameter definitions at current state.
    #' @return (`data.frame`)\cr
    #' A data.frame with formatted parameter definitions.
    param_export = function() {
      return(private$param_frame)
    },

    #' @description
    #' Clear all current parameter definitions.
    param_clear = function() {
      private$param_frame <- data.frame(
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

    #' @description
    #' Associate parameters from a list. Previous definitions are dropped.
    #' @param def_list (named `list()`)\cr
    #' A list containg parameter names with respective list of definitions.
    param_from_list = function(def_list) {
      if (!is.list(def_list)) {
        stop("The argument def_list should be a list.")
      }
      self$param_clear()
      do.call(self$param_append, def_list)
    },

    #' @description
    #' Append a list of parameter definitions.
    #' @param ... (named `list()`s)\cr
    #' Parameter names with respective list of definitions.
    param_append = function(...) {
      # convert definitions to list
      definitions <- list(...)

      if (length(definitions) == 0) {
        # reinitialize according to templated parameters
        stop(cli::format_error(c(
          "No definitions provided.",
          i = "To use template definitions use method param_template()."
        )))

      } else {
        # user definitions
        for (p in names(definitions)) {
          this <- definitions[[p]]
          private$param_add(
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

        # check if definitions match the model
        self$param_check()
      }
      invisible(self)
    },

    #' @description
    #' Overwrite current parameter definitions with template definitions.
    param_template = function() {
      # retrieve the template parameters file
      templ_params_file <- file.path(private$envir_path, "template_params.R")
      if (!file.exists((templ_params_file))) {
        stop(cli::format_error(c(
          "Impossible to automatically initialize parameters from template.",
          i = "No template file was found."
        )))
      }
      # remove previous definitions
      self$param_clear()
      # append definitions
      do.call(self$param_append, dget(templ_params_file))
      invisible(self)
    },

    #' @description
    #' Check parameter definitions from the implementation script.
    #' @param fill_template (`logical(1)`)\cr
    #' If `TRUE` fills the missing parameter with template definitions.
    param_check = function(fill_template = FALSE) {
      # retrieve the model script
      model_path <- file.path(private$envir_path, "src", "model.cpp")

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
      # FIXME : to be checked
      if (length(private$augm_frame) > 0) {
        x <- private$augm_frame[c("argsval", "argstype")]
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
        if (!(name %in% private$param_frame[["name"]])) {
          missing <- c(missing, name)
        }
      }

      # eventually fill with template definitions
      if (length(missing) > 0) {
        if (fill_template) {
          # retrieve the template parameters file
          templ_file <- file.path(private$envir_path, "template_params.R")
          if (!file.exists((templ_file))) {
            stop(cli::format_error(c(
              "Impossible to automatically fill missing parameters.",
              i = "No template file was found."
            )))
          }
          templ_params_list <- dget(templ_file)[missing]
          do.call(self$param_append, templ_params_list)
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
      for (name in private$param_frame[["name"]]) {
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
    },

    ############################################################################
    ############################################################################
    ############################################################################
    ############################################################################
    ############################################################################

    #' @description
    #' Export the augmented variable definitions at current state.
    #' @return (`data.frame`)\cr
    #' A data.frame with formatted augmented variable definitions.
    augm_export = function() {
      return(private$augm_frame)
    },

    #' @description
    #' Remove all current augmented variable definitions.
    augm_clear = function() {
      private$augm_frame <- data.frame(
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

    #' @description
    #' Associate augmented variables from a list.
    #' Previous definitions are dropped.
    #' @param def_list (named `list()`)\cr
    #' Augmented variable names with respective list of definitions.
    augm_from_list = function(def_list) {
      if (!is.list(def_list)) {
        stop("The argument def_list should be a list.")
      }
      self$augm_clear()
      do.call(self$augm_append, def_list)
    },

    #' @description
    #' Append a data frame of augmented variable definitions.
    #' @param ... (named `list()`s)\cr
    #' Augmented variable names with respective list of definitions.
    augm_append = function(...) {
      # convert definitions to list
      definitions <- list(...)

      if (length(definitions) == 0) {
        # reinitialize according to templated parameters
        stop(cli::format_error(c(
          "No definitions provided.",
          i = "To use template definitions use method param_template()."
        )))

      } else {
        # add each entry of the list
        for (a in names(definitions)) {
          this <- definitions[[a]]
          private$augm_add(
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

        # check if definitions match the model
        self$augm_check()
      }
      invisible(self)
    },


    #' @description
    #' Overwrite augmented variable definitions with template definitions.
    augm_template = function() {
      # retrieve the template augmented file
      templ_augm_file <- file.path(private$envir_path, "template_augm.R")
      if (!file.exists((templ_augm_file))) {
        stop(cli::format_error(c(
          "Impossible to automatically initialize augmented variables.",
          i = "No template file was found."
        )))
      }
      # remove previous definitions
      self$augm_clear()
      # append definitions
      do.call(self$augm_append, dget(templ_augm_file))
      invisible(self)
    },

    #' @description
    #' Check augmented definitions from the model implementation script.
    #' @param fill_template (`logical(1)`)\cr
    #' If `TRUE` fills missing augmented variables with template definitions.
    augm_check = function(fill_template = FALSE) {
      # retrieve the model script
      model_path <- file.path(private$envir_path, "src", "model.cpp")

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

      # UPDATE : check from infectivity profile?

      # re-order alphabetically
      model_augm_set <- model_augm_set[order(model_augm_set)]

      # find missing augmented variables
      missing <- c()
      for (name in model_augm_set) {
        if (!(name %in% private$augm_frame[["name"]])) {
          missing <- c(missing, name)
        }
      }

      # eventually fill with template definitions
      if (length(missing) > 0) {
        if (isTRUE(fill_template)) {
          # retrieve the template augmented file
          templ_augm_file <- file.path(private$envir_path, "template_augm.R")
          if (!file.exists((templ_augm_file))) {
            stop(cli::format_error(c(
              "Impossible to automatically initialize augmented variables.",
              i = "No template file was found."
            )))
          }
          templ_augm_list <- dget(templ_augm_file)[missing]
          do.call(self$augm_append, templ_augm_list)
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
      for (name in private$augm_frame[["name"]]) {
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