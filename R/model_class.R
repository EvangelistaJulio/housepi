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
    rcppmodel = NULL,
    envir_path = NA,
    setup = NULL,
    param = NULL,
    augm = NULL,
    traces_frame = NULL,

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
      setup <- lapply(private$setup[!is.na(private$setup)], as.character)
      private$rcppmodel$import_setup(setup)

      #----- 2) PARAMETERS -----#
      params <- data.frame(lapply(private$param$defs, as.character))
      private$rcppmodel$import_params(params)

      #----- 3) DATA -----#
      if (is.null(data)) {
        private$rcppmodel$read_data("<NA>")
      } else if (is.character(data)) {
        private$rcppmodel$read_data(data)
      } else if (is.data.frame(data)) {
        data <- data.frame(lapply(data, as.character))
        private$rcppmodel$import_data(data)
      } else {
        stop(cli::format_error(c(
          "Wrong type for data",
          i = "It should be data.frame or path/to/file."
        )))
      }

      #----- 4) AUGMENTED VARIABLES -----#
      augm <- data.frame(lapply(private$augm$defs, as.character))
      private$rcppmodel$import_augm(augm)

      # restore scientific notation option
      options(scipen = store_scipen)
    },

    # @description
    # Setup the list of configurations for the different analysis.
    # @param seed (`integer(1)`)\cr
    # Seed of the random number generator.
    # @param n_steps (`integer(1)`)\cr
    # Number of iterations.
    # @param n_burnin (`integer(1)`)\cr
    # Number of iterations to discard.
    # @param thin (`integer(1)`)\cr
    # Name of the parameter.
    # @param n_chains (`integer(1)`)\cr
    # Name of the parameter.
    # @param n_adapt (`integer(1)`)\cr
    # Name of the parameter.
    # @param n_threads (`integer(1)`)\cr
    # Name of the parameter.
    # @param rnd_init_params (`logical(1)`)\cr
    # Name of the parameter.
    # @param inv_temper (`double(1)`)\cr
    # Name of the parameter.
    setup_list = function(seed = NA,
                          n_steps = 1000,
                          n_burnin = 0,
                          thin = 1,
                          n_chains = 1,
                          n_adapt = 0,
                          n_threads = 1,
                          rnd_init_params = TRUE,
                          inv_temper = 1.0) {
      private$setup[["seed"]] <- seed
      private$setup[["n_steps"]] <- n_steps
      private$setup[["n_burnin"]] <- n_burnin
      private$setup[["thin"]] <- thin
      private$setup[["n_chains"]] <- n_chains
      private$setup[["n_adapt"]] <- n_adapt
      private$setup[["n_threads"]] <- n_threads
      private$setup[["rnd_init_params"]] <- rnd_init_params
      private$setup[["inv_temper"]] <- inv_temper
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
      private$rcppmodel <- m$ModelHousepi$new()

      # assign the empty setup list
      private$setup <- list()
      # assign the empty parameter data.frame
      private$param <- HousepiParameterClass$new()
      # assign the empty augmented data.frame
      private$augm <- HousepiAugmentedClass$new()

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
      return(private$rcppmodel$evaluate())
    },

    #' @description
    #' Run Bayesian MCMC inference and store internally the chain traces.
    #' @param data (`character(1)` | `data.frame`)\cr
    #' Path to a csv data file or a data.frame.
    #' @param progress (`logical(1)`)\cr
    #' If `FALSE` suppress the inference progress bar.
    #' @param ... further arguments passed to \code{setup}.
    infer = function(data = NULL, progress = TRUE, ...) {
      # assign setup
      do.call(private$setup_list, list(...))
      # internal build function
      private$build(data = data)
      # infer method from cpp model object
      private$rcppmodel$infer(progress)
      # load traces from last inference
      df <- private$rcppmodel$export_traces()
      # extract setup information for coda
      load_setup <- private$rcppmodel$export_setup()
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
    #' Calls coda version of lattice:xyplot for the traces traces.
    #' @param chain (`integer(1)`)
    #' Option to specifiy which chain to plot if multiple chains.
    #' @param col (`character()`)
    #' Color names.
    #' @param ...
    #' Further arguments pass to [xyplot][lattice::xyplot].
    traces_plot = function(chain = NA, col = "dodgerblue", ...) {
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
    },

    #' @description
    #' Run the simulation.
    #' @param data (`character(1)` | `data.frame`)\cr
    #' Path to a csv data file or a data.frame.
    #' @param progress (`logical(1)`)\cr
    #' If `FALSE` suppress the simulation progress bar.
    #' @param ... further arguments passed to \code{setup}.
    #' Named build options.
    simulate = function(data = NULL, progress = TRUE, ...) {
      # assign setup
      do.call(private$setup_list, list(...))
      # internal build function
      private$build(data = data)
      # simulate method from internal Rcpp object
      private$rcppmodel$simulate(progress)
      invisible(self)
    },

    #' @description
    #' Export the loaded database at current state.
    #' If the case, it also return the associated augmented variables.
    #' @return (`data.frame`)\cr
    #' A data.frame with individual data.
    data_export = function() {
      return(private$rcppmodel$export_data())
    },

    #' @description
    #' Print a summary of the loaded database.
    data_summary = function() {
      return(private$rcppmodel$show_summary())
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
      return(loo::loo(private$rcppmodel$export_loglike_matrix(), ...))
    },

    #' @description
    #' Use loo package to compute the widely applicable information criterion.
    #' @param ... arguments passed to [loo][loo::waic].
    #' #' For arguments and further details type ?loo::waic.
    waic = function(...) {
      return(loo::waic(private$rcppmodel$export_loglike_matrix(), ...))
    },

    #' @description
    #' Estimate the model evidence (or marginal likelihood).
    #' It uses the method of power posteriors.
    #' @param data (`character(1)` | `data.frame`)\cr
    #' Path to a csv data file or a data.frame.
    #' @param n_rungs (`integer`)\cr
    #' The number of rungs in the thermodynamic ladder.
    #' @param adaptive (`logical(1)`)\cr
    #' If `TRUE` uses adaptive scheme for thermodynamic ladder.
    #' Else use power fractions.
    #' @param progress (`logical(1)`)\cr
    #' If `FALSE` suppress the simulation progress bar.
    #' @return (`double(1)`)\cr
    #' The value of the model evidence.
    #' @param ... further arguments passed to  \code{setup}.
    evidence = function(data = NULL,
                        n_rungs = 40,
                        adaptive = TRUE,
                        progress = TRUE,
                        ...) {
      # assign setup
      do.call(private$setup_list, list(...))
      # internal build function
      private$build(data = data)
      # evaluate marginal from the internal Rcpp object
      output <- private$rcppmodel$evaluate_marginal(n_rungs, adaptive, progress)
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
      return(private$param$defs)
    },

    #' @description
    #' Clear all current parameter definitions.
    param_clear = function() {
      private$param$clear()
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
      invisible(self)
    },

    #' @description
    #' Append a list of parameter definitions.
    #' @param ... (named `list()`s)\cr
    #' Parameter names with respective list of definitions.
    param_append = function(...) {
      private$param$append(...)
      self$param_check(FALSE)
      invisible(self)
    },

    #' @description
    #' Check parameter definitions from the implementation script.
    #' @param fill_template (`logical(1)`)\cr
    #' If `TRUE` fills the missing parameter with template definitions.
    param_check = function(fill_template = FALSE) {
      private$param$check(private$envir_path, private$augm$defs, fill_template)
      invisible(self)
    },

    #' @description
    #' Overwrite current parameter definitions with template definitions.
    param_template = function() {
      private$param$template(private$envir_path)
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
      invisible(self)
    },

    #' @description
    #' Append a data frame of augmented variable definitions.
    #' @param ... (named `list()`s)\cr
    #' Augmented variable names with respective list of definitions.
    augm_append = function(...) {
      private$augm$append(...)
      self$augm_check(FALSE)
      invisible(self)
    },

    #' @description
    #' Check augmented definitions from the model implementation script.
    #' @param fill_template (`logical(1)`)\cr
    #' If `TRUE` fills missing augmented variables with template definitions.
    augm_check = function(fill_template = FALSE) {
      private$augm$check(private$envir_path, fill_template)
      invisible(self)
    },

    #' @description
    #' Overwrite augmented variable definitions with template definitions.
    augm_template = function() {
      private$augm$template(private$envir_path)
      invisible(self)
    }
  )
)