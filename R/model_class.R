#
# model_class.R
#
# Copyright (C) 2025 Francesco Camaglia - Institut Pasteur, Paris (France).
#

#' @title Household Model Class
#' @description
#' Compile household model, manage your data and parameters from a single class.
#' @export HouseholdModelClass
#' @docType class
HouseholdModelClass <- R6::R6Class(
  classname = "HouseholdModelClass",
  cloneable = FALSE,

  ###################
  ##### PRIVATE #####
  ###################

  private = list(
    model = NULL,
    setup = NULL,
    param = NULL,
    augm = NULL,
    group_rate = NULL,
    traces_frame = NULL,

    # @description
    # Associate dataset to the model and options for the setup.
    build = function(dataset = NULL,
                     n_steps = 1000,
                     n_burnin = 0,
                     thin = 1,
                     n_adapt = 0,
                     n_chains = 1,
                     n_threads = 1,
                     rnd_init_params = TRUE,
                     inv_temper = 1.0,
                     seed = NA,
                     verbose = FALSE) {
      #----- LOAD MODEL OBJECT -----#
      # UPDATE : I should assign group_rates here

      # if not developer mode, check definitions
      if (!private$model$as_developer) {
        self$param_assert()
        self$augm_assert()
      }

      # compile
      private$model$generate()
      private$model$load()

      #----- ASSOCIATE DEFINITIONS -----#

      # change scientific notation
      store_scipen <- options()$scipen
      options(scipen = 100)

      #----- 1) SETUP -----#
      private$setup[["n_steps"]] <- n_steps
      private$setup[["n_burnin"]] <- n_burnin
      private$setup[["thin"]] <- thin
      private$setup[["n_adapt"]] <- n_adapt
      private$setup[["n_chains"]] <- n_chains
      private$setup[["n_threads"]] <- n_threads
      private$setup[["rnd_init_params"]] <- rnd_init_params
      private$setup[["inv_temper"]] <- inv_temper
      private$setup[["seed"]] <- seed

      setup <- lapply(private$setup[!is.na(private$setup)], as.character)
      private$model$rcpp$import_setup(setup)
      if (isTRUE(verbose)) {
        cli::cli_alert_success("Imported setup definitions.")
      }

      #----- 2) PARAMETERS -----#
      params <- data.frame(lapply(private$param$defs, as.character))
      private$model$rcpp$import_params(params)
      if (isTRUE(verbose)) {
        cli::cli_alert_success("Imported parameter definitions.")
      }

      #----- 3) DATA -----#
      if (is.null(dataset)) {
        private$model$rcpp$read_data("<NA>")
        if (isTRUE(verbose)) {
          cli::cli_alert_success("Using previously loaded dataset.")
        }

      } else if (is.character(dataset)) {
        private$model$rcpp$read_data(dataset)
        if (isTRUE(verbose)) {
          cli::cli_alert_success("Imported dataset from file.")
        }

      } else if (is.data.frame(dataset)) {
        dataset <- data.frame(lapply(dataset, as.character))
        private$model$rcpp$import_data(dataset)
        if (isTRUE(verbose)) {
          cli::cli_alert_success("Imported dataset from data.frame")
        }

      } else {
        stop(cli::format_error(c(
          "Unsupported type for dataset.",
          i = "It should be data.frame or path/to/csv."
        )))
      }

      #----- 4) AUGMENTED -----#
      augm <- data.frame(lapply(private$augm$defs, as.character))
      private$model$rcpp$import_augm(augm)
      if (isTRUE(verbose)) {
        cli::cli_alert_success("Imported augmented variable definitions.")
      }

      # restore scientific notation option
      options(scipen = store_scipen)
      invisible(self)
    },

    # @description
    # Check that traces exists, otherwise raise error
    traces_check = function() {
      if (is.null(private$traces_frame)) {
        stop(cli::format_error(c(
          "No MCMC traces were found.",
          i = "Be sure to run your inference first."
        )))
      }
      invisible(self)
    }
  ),

  ############################################################################
  ############################################################################
  ############################################################################
  ############################################################################
  ############################################################################

  ###################
  ##### PUBLIC ######
  ###################

  public = list(
    #' @description
    #' Create a new instance of class [HouseholdModelClass][housepi::HouseholdModelClass].
    #' @param model (`character(1)`)\cr
    #' One between:
    #' - Template name;
    #' - path/to/folder with implementation;
    #' - path/to/folder with previously compiled model;
    #' - "developer".
    #' @param quiet (`logical(1)`)\cr
    #' If `TRUE` hide compilation commands.
    #' @param force (`logical(1)`)\cr
    #' If `TRUE` compile even if a cached model already exists.
    #' @param as_template (`logical(1)`)\cr
    #' If `TRUE` store the dynamic library as a local template.
    #' IF `FALSE` recompile the model at every new session.
    initialize = function(model,
                          quiet = TRUE,
                          force = FALSE,
                          as_template = FALSE) {
      # assign model definitions
      private$model <- HousepiGeneratorClass$new(model = model,
                                                 quiet = quiet,
                                                 force = force,
                                                 as_template = as_template)
      # assign the empty setup list
      private$setup <- list()
      # assign the empty parameter data.frame
      private$param <- HousepiParameterClass$new()
      # assign the empty augmented data.frame
      private$augm <- HousepiAugmentedClass$new()
      # assign the empty group rate list
      private$group_rate <- HousepiGroupRateClass$new()

      invisible(self)
    },

    #' @description
    #' Show model description in RStudio Help window.
    help = function() {
      templ_info_file <- file.path(private$model$path, "template_info.Rd")
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
    #' @param dataset (`character(1)` | `data.frame`)\cr
    #' Path to a csv file or a data.frame containing the dataset.
    #' @param n_threads (`integer(1)`)\cr
    #' Number of threads for the parallelization of the loglikelihood.
    #' @param rnd_init_params (`logical(1)`)\cr
    #' If `TRUE` the initial value of parameters is drawn from prior.
    #' Default is `FALSE`.
    #' @param inv_temper (`double(1)`)\cr
    #' Value of the inverse temperature in the range [0, 1].
    #' For standard usage this is 1.
    #' @param verbose (`logical(1)`)\cr
    #' If `TRUE` prints automatic build information.
    #' @return (`double(1)`)\cr
    #' The value of the loglikelihood.
    loglikelihood = function(dataset = NULL,
                             n_threads = 1,
                             rnd_init_params = FALSE,
                             inv_temper = 1.0,
                             verbose = FALSE) {
      # internal build function
      private$build(dataset = dataset, n_threads = n_threads,
                    rnd_init_params = rnd_init_params,
                    inv_temper = inv_temper, verbose = verbose)
      # evaluate loglikelihood from the internal Rcpp object
      return(private$model$rcpp$evaluate())
    },

    #' @description
    #' Run Bayesian MCMC inference and store internally the chain traces.
    #' @param dataset (`character(1)` | `data.frame`)\cr
    #' Path to a csv file or a data.frame containing the dataset.
    #' @param n_steps (`integer(1)`)\cr
    #' Number of MCMC steps.
    #' @param n_burnin (`integer(1)`)\cr
    #' Number of steps to discard from the sample.
    #' @param thin (`integer(1)`)\cr
    #' Steps discarded between each sample.
    #' @param n_adapt (`integer(1)`)\cr
    #' Number of steps with finite adptation before sampling.
    #' @param n_chains (`integer(1)`)\cr
    #' Number of independent MCMC chains.
    #' @param n_threads (`integer(1)`)\cr
    #' Number of threads for the parallelization of the loglikelihood.
    #' @param rnd_init_params (`logical(1)`)\cr
    #' If `TRUE` the initial value of the parameters is drawn from prior.
    #' @param inv_temper (`double(1)`)\cr
    #' Value of the inverse temperature.
    #' For standard usage this is 1.
    #' @param seed (`integer(1)`)\cr
    #' Seed of the random number generator.
    #' If `NA` uses the internal clock (default).
    #' @param progress (`logical(1)`)\cr
    #' If `FALSE` suppress the inference progress bar.
    #' @param verbose (`logical(1)`)\cr
    #' If `TRUE` prints automatic build information.
    infer = function(dataset = NULL,
                     n_steps = 1000,
                     n_burnin = 0,
                     thin = 1,
                     n_adapt = 0,
                     n_chains = 1,
                     n_threads = 1,
                     rnd_init_params = TRUE,
                     inv_temper = 1.0,
                     seed = NA,
                     progress = TRUE,
                     verbose = FALSE) {
      # internal build function
      private$build(dataset, n_steps = n_steps, n_burnin = n_burnin,
                    thin = thin, n_adapt = n_adapt, n_chains = n_chains,
                    n_threads = n_threads, rnd_init_params = rnd_init_params,
                    inv_temper = inv_temper, seed = seed, verbose = verbose)

      # infer method from cpp model object
      private$model$rcpp$infer(progress)
      # load traces from last inference
      df <- private$model$rcpp$export_traces()
      # extract setup information for coda
      load_setup <- private$model$rcpp$export_setup()

      # Check if multiple chains
      if ("CHAIN__" %in% colnames(df)) {
        df_list <- lapply(split(df, df[["CHAIN__"]]), function(x) {
          coda::mcmc(
            data = x[, !names(x) %in% "CHAIN__"],
            start = load_setup[["n_burnin"]] + load_setup[["thin"]],
            end = load_setup[["n_steps"]],
            thin = load_setup[["thin"]]
          )
        })
        private$traces_frame <- coda::as.mcmc.list(df_list)

      } else {
        private$traces_frame <- coda::mcmc(
          data = df,
          start = load_setup[["n_burnin"]] + load_setup[["thin"]],
          end = load_setup[["n_steps"]],
          thin = load_setup[["thin"]]
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
      imask <- private$param$defs["infer"] == TRUE
      inames <- private$param$defs[imask, "name"]
      output <- data.frame(
        Mean = res[["statistics"]][inames, "Mean"],
        SD = res[["statistics"]][inames, "SD"],
        Q2.5 = res[["quantiles"]][inames, "2.5%"],
        Q97.5 = res[["quantiles"]][inames, "97.5%"]
      )
      return(output)
    },

    #' @description
    #' Calls coda wrapper of lattice:xyplot for the traces traces.
    #' @param chain (`integer(1)`)
    #' Option to specifiy which chain to plot if multiple chains.
    #' @param col (`character()`)
    #' Color names.
    #' @param ...
    #' Further arguments pass to [xyplot][lattice::xyplot].
    traces_xyplot = function(chain = NA,
                             col = "dodgerblue",
                             ...) {
      private$traces_check()
      xlab <- "iteration"

      if (coda::is.mcmc.list(private$traces_frame)) {
        if (is.na(chain)) {
          lattice::xyplot(
            private$traces_frame, auto.key = TRUE, xlab = xlab, ...
          )
        } else {
          lattice::xyplot(
            private$traces_frame[chain], col = col, xlab = xlab, ...
          )
        }
      } else {
        lattice::xyplot(private$traces_frame, col = col, xlab = xlab, ...)
      }
      invisible(self)
    },

    #' @description
    #' Run the simulation.
    #' @param dataset (`character(1)` | `data.frame`)\cr
    #' Path to a csv file or a data.frame containing the dataset.
    #' @param n_threads (`integer(1)`)\cr
    #' Number of threads for the parallelization of the loglikelihood.
    #' @param rnd_init_params (`logical(1)`)\cr
    #' If `TRUE` the initial value of the parameters is drawn from prior.
    #' Default is `FALSE`.
    #' @param seed (`integer(1)`)\cr
    #' Seed of the random number generator.
    #' If `NA` uses the internal clock (default).
    #' @param progress (`logical(1)`)\cr
    #' If `FALSE` suppress the simulation progress bar.
    #' @param verbose (`logical(1)`)\cr
    #' If `TRUE` prints automatic build information.
    simulate = function(dataset = NULL,
                        n_threads = 1,
                        rnd_init_params = FALSE,
                        seed = NA,
                        progress = TRUE,
                        verbose = FALSE) {
      # internal build function
      private$build(dataset, n_threads = n_threads,
                    rnd_init_params = rnd_init_params,
                    seed = seed, verbose = verbose)
      # call simulate method from internal Rcpp object
      private$model$rcpp$simulate(progress)
      invisible(self)
    },

    #' @description
    #' Export the current dataset at current state.
    #' If the case, it also return the associated augmented variables.
    #' @return (`data.frame`)\cr
    #' A data.frame with individual level information.
    data_export = function() {
      if (is.null(private$model$rcpp)) {
        stop(cli::format_error(c(
          "Impossible to export data before model build.",
        )))
      }
      return(private$model$rcpp$export_data())
    },

    #' @description
    #' Print a summary of the current dataset.
    data_summary = function() {
      if (is.null(private$model$rcpp)) {
        stop(cli::format_error(c(
          "Impossible to show data summary before model build.",
        )))
      }
      return(private$model$rcpp$show_summary())
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
      if (is.null(private$model$rcpp)) {
        stop(cli::format_error(c(
          "Impossible to run loo command before model build.",
        )))
      }
      return(loo::loo(private$model$rcpp$export_household_loglike_traces(), ...))
    },

    #' @description
    #' Use loo package to compute the widely applicable information criterion.
    #' @param ... arguments passed to [loo][loo::waic].
    #' #' For arguments and further details type ?loo::waic.
    waic = function(...) {
      if (is.null(private$model$rcpp)) {
        stop(cli::format_error(c(
          "Impossible to run waic command before model build.",
        )))
      }
      return(loo::waic(private$model$rcpp$export_household_loglike_traces(), ...))
    },

    #' @description
    #' Estimate the model evidence (or marginal likelihood).
    #' It uses the method of power posteriors.
    #' @param dataset (`character(1)` | `data.frame`)\cr
    #' Path to a csv file or a data.frame containing the dataset.
    #' @param n_rungs (`integer`)\cr
    #' The number of rungs in the thermodynamic ladder.
    #' @param scheme (`character(1)`)\cr
    #' If "adaptive" uses adaptive scheme for thermodynamic ladder (default).
    #' Other available scheme is "power-fractions".
    #' @param n_steps (`integer(1)`)\cr
    #' Number of MCMC steps.
    #' @param n_burnin (`integer(1)`)\cr
    #' Number of steps to discard from the sample.
    #' @param thin (`integer(1)`)\cr
    #' Steps discarded between each sample.
    #' @param n_adapt (`integer(1)`)\cr
    #' Number of steps with finite adptation before sampling.
    #' @param n_chains (`integer(1)`)\cr
    #' Number of independent MCMC chains.
    #' @param n_threads (`integer(1)`)\cr
    #' Number of threads for the parallelization of the loglikelihood.
    #' @param seed (`integer(1)`)\cr
    #' Seed of the random number generator.
    #' If `NA` uses the internal clock (default).
    #' @param progress (`logical(1)`)\cr
    #' If `FALSE` suppress the simulation progress bar.
    #' @param verbose (`logical(1)`)\cr
    #' If `TRUE` prints automatic build information.
    #' @return (`double(1)`)\cr
    #' The value of the model evidence.
    evidence = function(dataset = NULL,
                        n_rungs = 40,
                        scheme = "adaptive",
                        n_steps = 1000,
                        n_burnin = 0,
                        thin = 1,
                        n_adapt = 0,
                        n_chains = 1,
                        n_threads = 1,
                        seed = NA,
                        progress = TRUE,
                        verbose = FALSE) {
      # internal build function
      private$build(dataset, n_steps = n_steps, n_burnin = n_burnin,
                    thin = thin, n_adapt = n_adapt, n_chains = n_chains,
                    n_threads = n_threads, rnd_init_params = TRUE,
                    seed = seed, verbose = verbose)
      # evaluate marginal from the internal Rcpp object
      output <- private$model$rcpp$evaluate_marginal(n_rungs, scheme, progress)
      # drop traces eventually
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
      return(private$param$defs)
    },

    #' @description
    #' Clear all current parameter definitions.
    param_clear = function() {
      private$param$clear()
      invisible(self)
    },

    #' @description
    #' Add a new parameter definition.
    #' @param name (`character(1)`)\cr
    #' Name of the parameter.
    #' @param value (`double(1)`)\cr
    #' Initial value of the parameter.
    #' @param infer (`logical(1)`)\cr
    #' If `FALSE` the parameter is treated as a constant.
    #' @param discrete (`logical(1)`)\cr
    #' If `TRUE` the parameter is discrete.
    #' @param prior (`character(1)`)\cr
    #' Name of the prior distribution.
    #' @param args (named `list()`)\cr
    #' Arguments for the prior distribution.
    #' @param rule (`character(1)`)\cr
    #' Name of the inference proposal rule.
    #' @param rule_scale (`double(1)`)\cr
    #' Scale parameter for the proposal rule.
    #' Changed automatically during inference if `n_adpat`>0.
    #' @param range (named `list()`)\cr
    #' The prior truncation 'min', 'max' or no bound.
    param_add = function(name = NA,
                         value = NA,
                         infer = NA,
                         discrete = FALSE,
                         prior = NA,
                         args = list(),
                         rule = NA,
                         rule_scale = NA,
                         range = list()) {
      private$param$add(name = name, value = value, infer = infer,
                        discrete = discrete, prior = prior, args = args,
                        rule = rule, rule_scale = rule_scale, range = range)
      invisible(self)
    },

    #' @description
    #' Associate parameters from a list. Previous definitions are dropped.
    #' @param definitions (named `list()`)\cr
    #' A list containg parameter names with respective list of definitions.
    #' @param append (`logical(1)`)\cr
    #' If `TRUE` appends the list to current definition.
    #' If `FALSE` clears previous definitions (default).
    param_from_list = function(definitions, append = FALSE) {
      private$param$from_list(definitions, append = append)
      invisible(self)
    },

    #' @description
    #' Check parameter definitions from the implementation script.
    #' @param fill (`logical(1)`)\cr
    #' If `TRUE` fills the missing parameter with template definitions.
    param_assert = function(fill = FALSE) {
      if (private$model$as_developer) {
        cli::cli_alert_warning(cli::format_warning(
          "No check available in developer mode."
        ))
      } else {
        private$param$check(private$model$path, private$augm$defs, fill)
      }
      invisible(self)
    },

    #' @description
    #' Overwrite current parameter definitions with template definitions.
    param_template = function() {
      if (private$model$as_developer) {
        cli::cli_alert_warning(cli::format_warning(
          "No template available in developer mode."
        ))
      } else {
        private$param$template(private$model$path)
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
      return(private$augm$defs)
    },

    #' @description
    #' Remove all current augmented variable definitions.
    augm_clear = function() {
      private$augm$clear()
      invisible(self)
    },

    #' @description
    #' Add a new augmented variable definition.
    #' @param name (`character(1)`)\cr
    #' Name of the augmented variable.
    #' @param value (`double(1)` | `character(1)`)\cr
    #' Scalar value
    #' or parameter name (followed by `@P`)
    #' or data name (followed by `@D`).
    #' @param discrete (`logical(1)`)\cr
    #' IF `TRUE` the parameter is discrete.
    #' @param model (`character(1)`)\cr
    #' Name of the model distribution.
    #' @param args (named `list()`)\cr
    #' Arguments for the model distribution.
    #' @param rule (`character(1)`)\cr
    #' Name of the inference proposal rule.
    #' @param rule_scale (`double(1)`)\cr
    #' Scale parameter for the proposal rule.
    #' @param rule_from (`double(1)` | `character(1)`)\cr
    #' Advanced options for non random walk moves.
    #' @param range  (named `list()`)\cr
    #' List with 'min' and/or 'max' values (or none).
    #' @param init_rand (`character(1)`)\cr
    #' If `TRUE` draw initial vale from model distribution.
    #' @param status (`integer(1)`)\cr
    #' The associated status of the definition.
    #' Default NA is for all infectious.
    augm_add = function(name = NA,
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
      private$augm$add(name = name, value = value, discrete = discrete,
                       model = model, args = args, rule = rule,
                       rule_scale = rule_scale, rule_from = rule_from,
                       range = range, init_rand = init_rand, status = status)
      invisible(self)
    },

    #' @description
    #' Associate augmented variables from a list.
    #' Previous definitions are dropped.
    #' @param definitions (named `list()`)\cr
    #' Augmented variable names with respective list of definitions.
    #' @param append (`logical(1)`)\cr
    #' If `TRUE` appends the list to current definition.
    #' If `FALSE` clears previous definitions (default).
    augm_from_list = function(definitions, append = FALSE) {
      private$augm$from_list(definitions, append = append)
      invisible(self)
    },

    #' @description
    #' Check augmented definitions from the model implementation script.
    #' @param fill (`logical(1)`)\cr
    #' If `TRUE` fills missing augmented variables with template definitions.
    augm_ = function(fill = FALSE) {
      if (private$model$as_developer) {
        cli::cli_alert_warning(cli::format_warning(
          "No check available in developer mode."
        ))
      } else {
        private$augm$check(private$model$path, fill)
      }
      invisible(self)
    },

    #' @description
    #' Overwrite augmented variable definitions with template definitions.
    augm_template = function() {
      if (private$model$as_developer) {
        cli::cli_alert_warning(cli::format_warning(
          "No template available in developer mode."
        ))
      } else {
        private$augm$template(private$model$path)
      }
      invisible(self)
    },

    ############################################################################
    ############################################################################
    ############################################################################
    ############################################################################
    ############################################################################

    #' @description
    #' Remove all current augmented variable definitions.
    group_rate_clear = function() {
      private$group_rate$clear()
      invisible(self)
    },

    #' @description
    #' Add a new group rate definition.
    #' @param name (`character(1)`)\cr
    #' Name of the group rate.
    #' @param features (`list()`)\cr
    #' Feature(s) associated to group rate as dataset headers.
    #' Corresponding dataset column have to be categorical.
    #' @param prefix (`character(1)`)\cr
    #' Prefix to parameter names associated to the group rate.
    #' @param reference (named `list()`)\cr
    #' Choose reference group providing corresponding features.
    #' Default is the group with the largest sample size.
    #' @param independent (`logical(1)`)\cr
    #' If `TRUE` features are considered separately.
    group_rate_add = function(name,
                              features,
                              prefix = NA,
                              reference = NULL,
                              independent = TRUE) {
      private$group_rate$add(name = name, features = features, prefix = prefix,
                             reference = reference, independent = independent)
      invisible(self)
    }
  )
)
