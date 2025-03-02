#
# model_class.R
#
# Copyright (C) 2025 CAMAGLIA Francesco - Institut Pasteur, Paris (France).
#

#' @title
#' Household Model R6 Class
#' @description
#' Compile household model, manage your data and parameters from a single class.
#'
#' @export HouseholdModel
#' @docType class
HouseholdModel <- R6::R6Class(
  classname = "HouseholdModel",
  cloneable = TRUE,

  ############################################################################
  ############################################################################
  ################################## PRIVATE #################################
  ############################################################################
  ############################################################################

  private = list(
    .model = NULL,
    .param = NULL,
    .augm = NULL,
    .group_rate = NULL,
    .infec_profile = NULL,

    # @description
    # Assert if parameter definitions match the implementation script.
    #
    # @return (`logical(1)`)
    .assert_param = function() {
      # extract requirements from reference file
      # note is.null(c()) -> TRUE
      if (private$.model$is_blank) {
        requirements <- c()
      } else {
        refer_file <- private$.model$reference
        requirements <- private$.param$require(refer_file)
      }
      requirements <- c(requirements,
                        private$.augm$param_require,
                        private$.infec_profile$param_require)
      return(private$.param$assert(requirements))
    },

    # @description
    # Assert if augmented definitions match the model implementation script.
    #
    # @return (`logical(1)`)
    .assert_augm = function() {
      if (private$.model$is_blank) {
        return(TRUE)
      } else {
        refer_file <- private$.model$reference
        requirements <- private$.augm$require(refer_file)
        return(private$.augm$assert(requirements))
      }
    },

    # @description
    # Assert group rate definitions from the model implementation script.
    #
    # @return (`logical(1)`)\cr
    .assert_group_rate = function() {
      if (private$.model$is_blank) {
        return(TRUE)
      } else {
        refer_file <- private$.model$reference
        requirements <- private$.group_rate$require(refer_file)
        return(private$.group_rate$assert(requirements))
      }
    },

    # @description
    # Assert that dataset headers match implementation file.
    #
    # @param refer_file (`character()`)\cr
    # Path/to/.refer source file.
    # @param data_headers (`character()`)\cr
    # Vector with the dataset headers.
    #
    # @return (`logical(1)`)\cr
    .assert_data = function(refer_file, data_headers) {
      # logical state for complete definitions
      is_def_complete <- TRUE

      # vector of regular expression to use for assertion
      regular_expr <- c(
        # from : .data(, "string" or augm(, 'string'
        '.data\\s*\\(\\s*[^,]*\\s*,\\s*["\'](.*?)["\']'
      )
      # retrieve required data headers from regular expressions
      require_set <- template_find_requirement(refer_file, regular_expr)

      if (length(require_set) == 0) {
        # if requirement set is NULL it means nothing can be said
        is_def_complete <- NA

      } else {
        # find missing data headers
        missing_set <- template_find_missing(require_set,
                                             data_headers)
        if (length(missing_set) > 0) {
          is_def_complete <- FALSE
          warning_h(
            "Missing data headers(s):",
            character_to_string_list(missing_set)
          )
        }
      }

      return(is_def_complete)
    },

    ############################################################################
    ############################################################################
    ############################################################################

    # @description
    # Pre-build to deal with group rates.
    #
    # @param dataset (`data.frame`)\cr
    # The dataset from which to extract group rates.
    .build_group_rate = function(dataset) {
      # proceed only if group rate definitions are new and complete
      if (isFALSE(
        private$.group_rate$cpp_assert(dataset, private$.model$reference)
      )) {
        # NOTE: at the moment I just cache group rate changes
        # search for the cache file
        source_file <- file.path(private$.model$cache_dir,
                                 "src", "model.cpp")
        # deal with the case cache has been deleted
        if (!file.exists(source_file)) {
          # check if it's a template
          model_choice <- private$.model$path
          if (file.exists(file.path(model_choice, "model.cpp"))) {
            # re-assign model object forcing copy of source to cache
            private$.model <- HousepiGenerator$new(model_choice,
                                                   force = TRUE)
            # assign new cache file
            source_file <- file.path(private$.model$cache_dir,
                                     "src", "model.cpp")
          } else {
            stop_h("Impossible to find reference file.")
          }
        }

        # write the cpp script and save new reference
        private$.group_rate$cpp_write(dataset, source_file)
        private$.model$save_reference()
      }

      # append new parameter definitions
      new_defs <- list()
      for (rate_id in names(private$.group_rate$logic_table)) {
        for (par_id in private$.group_rate$logic_table[[rate_id]]$.param_name) {
          if (par_id != "REFERENCE") {
            new_defs[[par_id]] <- private$.group_rate$param_defs[[rate_id]]
          }
        }
      }
      if (length(new_defs) > 0) {
        self$load_param_list(from = new_defs, append = TRUE)
      }

      invisible(self)
    },

    # @description
    # Associate dataset and parameters for a minimal build.
    # @param dataset (`data.frame`)\cr
    # A data.frame containing the dataset.
    .pre_build = function(dataset = NULL) {
      if (private$.model$is_blank) {
        # blank dataset
        # note: no assert definitions in such case
        id_trace <- 1
        dataset <- data.frame(id = c(id_trace),
                              house_id = c("HH_1"),
                              status = c(1))

      } else if (length(dataset) > 0) {
        # check that dataset type is correct
        if (!is.data.frame(dataset)) {
          stop_h(
            "Unsupported type for dataset.",
            "Dataset should be a data.frame."
          )
        }

        # deal with group rates
        private$.build_group_rate(dataset)

        # assert dataset headers
        if (isFALSE(private$.assert_data(private$.model$reference,
                                         colnames(dataset)))) {
          stop_h("Incomplete dataset headers.")
        }
      }

      # assert parameter definitions
      if (isFALSE(private$.assert_param())) {
        stop_h("Incomplete parameter definitions.")
      }

      # 0) model generation
      private$.model$generate()

      # 1) parameters
      private$.model$import_parameter(self$param_defs)

      # 2) dataset
      private$.model$import_dataset(dataset)

      invisible(self)
    },

    # @description
    # Full model build with complete list of definitions.
    #
    # @param dataset (`data.frame`)\cr
    # A data.frame containing the dataset.
    # @param sample_size (`integer(1)`)\cr
    # Sample size of a single chain according to posterior distribution.
    # @param n_steps (`integer(1)`)\cr
    # Number of MCMC steps.
    # @param n_burnin (`integer(1)`)\cr
    # Number of steps to discard at the beginning of the sample.
    # @param thin (`integer(1)`)\cr
    # Number of steps to discard between each sample.
    # @param n_adapt (`integer(1)`)\cr
    # Number of steps with finite adptation before sampling.
    # @param n_chains (`integer(1)`)\cr
    # Number of independent MCMC chains.
    # @param n_threads (`integer(1)`)\cr
    # Number of threads for the parallelization of the loglikelihood.
    # @param rand_init (`logical(1)`)\cr
    # If `TRUE` the initial value of the parameters is drawn from prior.
    # @param id_trace (`character(1)`)\cr
    # Individual identifier to show augmented traces. Default is NA.
    # @param inv_temper (`double(1)`)\cr
    # Value of the inverse temperature.
    # For standard usage this is 1.
    # @param seed (`integer(1)`)\cr
    # Seed of the random number generator.
    # If `NA` uses the internal clock (default).
    # @param progress (`logical(1)`)\cr
    # If `TRUE` prints the progress bar.
    .build = function(dataset = NULL,
                      sample_size = NA,
                      n_steps = 0,
                      n_burnin = 0,
                      thin = 1,
                      n_adapt = 0,
                      n_chains = 1,
                      n_threads = 1,
                      rand_init = TRUE,
                      id_trace = NA,
                      inv_temper = 1.0,
                      seed = NA,
                      progress = FALSE) {
      # run model pre-build
      private$.pre_build(dataset)

      # assert augmented variable definitions
      if (isFALSE(private$.assert_augm())) {
        stop_h("Incomplete augmented variable definitions.")
      }

      # 3) setup
      if (private$.model$is_blank) id_trace <- 1
      private$.model$import_setup(list(
        sample_size = sample_size, n_steps = n_steps, thin = thin,
        n_burnin = n_burnin, n_adapt = n_adapt, n_chains = n_chains,
        n_threads = n_threads, rand_init = rand_init,
        id_trace = id_trace, inv_temper = inv_temper,
        seed = seed, progress = progress
      ))

      # 4) group & contact rates
      # UPDATE: rebuild group rates <==> something has changed
      private$.model$build_all_rates()

      # 5) augmented variables
      if (length(self$augm_defs) > 0) {
        private$.model$import_augmented(self$augm_defs)
      }

      # 6) infectivity profile
      if (length(private$.infec_profile$defs) > 0) {
        private$.model$import_infectivity_profile(private$.infec_profile$defs)
      }

      invisible(self)
    }
  ),

  ############################################################################
  ############################################################################
  ################################## ACTIVE ##################################
  ############################################################################
  ############################################################################

  active = list(
    #' @field dataset (`data.frame`)\cr
    #' Dataset with individual level information of current dataset.
    #' If the case, it also return the associated augmented variables.
    dataset = function(value) {
      if (missing(value)) return(private$.model$dataset)
      stop_h("'dataset' is read only.")
    },

    #' @field setup (`list()`)\cr
    #' Last setup used by the build process.
    setup = function(value) {
      if (missing(value)) return(private$.model$setup)
      stop_h("'setup' is read only.")
    },

    #' @field param_defs (`data.frame`)\cr
    #' Parameter formatted definitions.
    param_defs = function(value) {
      if (missing(value)) return(private$.param$defs)
      stop_h("'param_defs' is read only.")
    },

    #' @field augm_defs (`data.frame`)\cr
    #' Augmented variable formatted definitions.
    augm_defs = function(value) {
      if (missing(value)) return(private$.augm$defs)
      stop_h("'augm_defs' is read only.")
    },

    #' @field traces (`matrix()` | `list()` of `matrix()`)\cr
    #' Traces of last inference.
    traces = function(value) {
      if (missing(value)) return(private$.model$traces)
      stop_h("'traces' is read only.")
    },

    #' @field traces_no_loglike (`matrix()` | `list()` of `matrix()`)\cr
    #' Traces of last inference without loglikelihood.
    traces_no_loglike = function(value) {
      if (missing(value)) return(private$.model$traces_no_loglike)
      stop_h("'traces_no_loglike' is read only.")
    },

    #' @field traces_as_coda (`coda::mcmc` | `coda::mcmc.list`)\cr
    #' Traces of last inference as a `coda` object.
    traces_as_coda = function(value) {
      if (missing(value)) return(coerce_traces_to_coda(private$.model$traces))
      stop_h("'traces' is read only.")
    },

    #' @field traces_results (`data.frame`)\cr
    #' Results extracted from the traces of last inference.
    traces_results = function(value) {
      if (missing(value)) {
        return(retrieve_traces_results(private$.model$traces_no_loglike))
      }
      stop_h("'traces_results' is read only.")
    }
  ),

  ############################################################################
  ############################################################################
  ################################## PUBLIC ##################################
  ############################################################################
  ############################################################################

  public = list(
    #' @description
    #' New instance of class
    #' [HouseholdModel][housepi::HouseholdModel].
    #'
    #' @param model (`character(1)`)\cr
    #' One between:
    #' - template name;
    #' - path/to/folder with implementation;
    #' - path/to/folder with previously compiled model;
    #' - "blank": template name for blank model (default).
    #' @param force (`logical(1)`)\cr
    #' If `TRUE` compile even if a cached model already exists.
    #' @param quiet_pkgbuild (`logical(1)`)\cr
    #' If `TRUE` hide `pkgbuild` compilation commands.
    #' @param as_template (`logical(1)`)\cr
    #' If `TRUE` store the dynamic library as a local template.
    #' IF `FALSE` recompile the model at every new session.
    initialize = function(model = NA,
                          force = FALSE,
                          quiet_pkgbuild = TRUE,
                          as_template = FALSE) {
      if (is.na(model)) {
        model <- "blank"
        warning_h("No model specified. Falling back to \"blank\" model.")
      }

      # assign model definitions
      private$.model <- HousepiGenerator$new(model = model,
                                             force = force,
                                             quiet_pkgbuild = quiet_pkgbuild,
                                             as_template = as_template)
      # initialize internal fields
      self$drop_all()
      invisible(self)
    },

    #' @description
    #' Override the default printing behavior.
    print = function() {
      cli::cli_alert_info("Household model {.pkg {private$.model$model_name}}.")
      invisible(self)
    },

    #' @description
    #' Show model description in RStudio Help window.
    help = function() {
      if (requireNamespace("rstudioapi", quietly = TRUE)) {
        templ_info_file <- file.path(private$.model$path, "template_info.Rd")
        if (!file.exists((templ_info_file))) {
          stop_h(
            "Impossible to show model description.",
            "No template info was found."
          )
        }
        rstudioapi::previewRd(templ_info_file)
      } else {
        warning_h(
          "Help method not yet supported outside of R studio.",
          "Sorry for the inconvenience."
        )
      }

      invisible(self)
    },

    ############################################################################
    ############################################################################
    ############################################################################

    #' @description
    #' Drop all current parameter definitions.
    drop_param_defs = function() {
      private$.param <- HousepiParameter$new()
      invisible(self)
    },

    #' @description
    #' Drop all current augmented variable definitions.
    drop_augm_defs = function() {
      private$.augm <- HousepiAugmented$new()
      invisible(self)
    },

    #' @description
    #' Drop all current augmented variable definitions.
    drop_group_rate_defs = function() {
      private$.group_rate <- HousepiGroupRate$new()
      invisible(self)
    },

    #' @description
    #' Drop all current definitions.
    drop_all = function() {
      self$drop_param_defs()
      self$drop_augm_defs()
      self$drop_group_rate_defs()
      private$.infec_profile <- HousepiInfecProfile$new()
      invisible(self)
    },

    ############################################################################
    ############################################################################
    ############################################################################

    #' @description
    #' Add a new parameter definition.
    #'
    #' @param name (`character(1)`)\cr
    #' Name of the parameter.
    #' @param value (`double(1)`)\cr
    #' Initial value of the parameter.
    #' @param constant (`logical(1)`)\cr
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
    #' @param min (`double(1)`)\cr
    #' Value of the minimum.
    #' @param max (`double(1)`)\cr
    #' Value of the maximum.
    add_param = function(name = NA,
                         value = NA,
                         constant = NA,
                         discrete = FALSE,
                         prior = NA,
                         args = list(),
                         rule = NA,
                         rule_scale = NA,
                         min = NA,
                         max = NA) {
      private$.param$add(name = name, value = value, constant = constant,
                         discrete = discrete, prior = prior, args = args,
                         rule = rule, rule_scale = rule_scale,
                         min = min, max = max)
      invisible(self)
    },

    #' @description
    #' Associate parameters from a list.
    #'
    #' @param ...
    #' Parameter names with respective list of definitions.
    #' @param from (named `list()`)\cr
    #' Named list of parameter definitions.
    #' @param append (`logical(1)`)\cr
    #' If `TRUE` appends the list to current definition (default).
    #' If `FALSE` drop all other parameter definitions.
    load_param_list = function(...,
                               from = NA,
                               append = TRUE) {
      private$.param$load_list(..., from = from, append = append)
      invisible(self)
    },

    #' @description
    #' Load parameter definitions from model template.
    #'
    #' @param append
    #' If `TRUE` append only missing definitions (default).
    #' Else drop all previous definitions.
    load_param_template = function(append = TRUE) {
      if (private$.model$is_blank) {
        warning_h("No template parameters available for blank model.")

      } else {
        # retrieve the template parameters file
        private$.param$load_template(
          file.path(private$.model$path, "template_params.R"),
          private$.param$require(private$.model$reference),
          isTRUE(append)
        )
      }
      invisible(self)
    },

    ############################################################################
    ############################################################################
    ############################################################################

    #' @description
    #' Add a new augmented variable definition.
    #'
    #' @param name (`character(1)`)\cr
    #' Name of the augmented variable.
    #' @param status (`integer(1)`)\cr
    #' The associated status of the definition.
    #' Default NA is for all infectious.
    #' @param value (`double(1)` | `character(1)`)\cr
    #' One between:
    #' - Scalar value;
    #' - Parameter name (optionally followed by `@P`);
    #' - Data name (followed by `@D`).
    #' @param constant (`logical(1)`)\cr
    #' If `FALSE` the augmented variable is treated as a constant.
    #' @param discrete (`logical(1)`)\cr
    #' IF `TRUE` the augmenter variable is discrete.
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
    #' `FIXME: missing`
    #' @param min (`double(1)` | `character(1)`)\cr
    #' Minimum value. One between:
    #' - Scalar value;
    #' - Parameter name (optionally followed by `@P`);
    #' - Data name (followed by `@D`).
    #' @param max (`double(1)` | `character(1)`)\cr
    #' Maxmimum tvalue. One between:
    #' - Scalar value;
    #' - Parameter name (optionally followed by `@P`);
    #' - Data name (followed by `@D`).
    #' @param init_rand (`character(1)`)\cr
    #' If `TRUE` draw initial vale from model distribution.
    add_augm = function(name = NA,
                        status = NA,
                        value = NA,
                        constant = NA,
                        discrete = FALSE,
                        model = NA,
                        args = list(),
                        rule = NA,
                        rule_scale = NA,
                        rule_from = NA,
                        min = NA,
                        max = NA,
                        init_rand = FALSE) {
      private$.augm$add(name = name, status = status,
                        value = value, constant = constant, discrete = discrete,
                        model = model, args = args, rule = rule,
                        rule_scale = rule_scale, rule_from = rule_from,
                        min = min, max = max, init_rand = init_rand)
      invisible(self)
    },

    #' @description
    #' Associate augmented variables from a list.
    #'
    #' @param ...
    #' Augmented variable names with respective list of definitions.
    #' @param from (named `list()`)\cr
    #' Named list of augmented variable definitions.
    #' @param append (`logical(1)`)\cr
    #' If `TRUE` appends the list to current definition (default).
    #' Else drop all other augmented variable definitions.
    load_augm_list = function(...,
                              from = NA,
                              append = TRUE) {
      private$.augm$load_list(..., from = from, append = append)
      invisible(self)
    },

    #' @description
    #' Load augmented variables from template definitions.
    #'
    #' @param append
    #' If `TRUE` append only missing definitions (default).
    #' Else drop all previous definitions.
    load_augm_template = function(append = TRUE) {
      if (private$.model$is_blank) {
        warning_h("No template augmented variables available for blank model.")

      } else {
        # retrieve the template augmented file
        private$.augm$load_template(
          file.path(private$.model$path, "template_augm.R"),
          private$.augm$require(private$.model$reference),
          isTRUE(append)
        )
      }

      invisible(self)
    },

    ############################################################################
    ############################################################################
    ############################################################################

    #' @description
    #' Add definitions of a group rate by given data feature(s).
    #' Reference value has unitary rate and can be specified
    #' (default is most represented group).
    #'
    #' @param name (`character(1)`)\cr
    #' Name of the group rate.
    #' @param by (`list()`)\cr
    #' Feature(s) associated to group rate as dataset headers.
    #' Corresponding dataset column have to be categorical.
    #' @param prefix (`character(1)`)\cr
    #' Prefix to parameter names associated to the group rate.
    #' @param reference (named `list()`)\cr
    #' Choose reference group providing corresponding features.
    #' Default is the group with the largest sample size.
    #' @param independent (`logical(1)`)\cr
    #' If `TRUE` then features are considered independently.
    #' @param prior (`character(1)`)\cr
    #' Name of the prior distribution.
    #' Default associated parameter definitions
    #' have prior "lnorm" and `args = list(meanlog = 0, sdlog = 3)`.
    #' @param args (named `list()`)\cr
    #' Arguments for the prior distribution.
    #' @param rule (`character(1)`)\cr
    #' Name of the inference proposal rule.
    #' @param rule_scale (`double(1)`)\cr
    #' Scale parameter for the proposal rule.
    #' Changed automatically during inference if `n_adpat` > 0.
    #' @param min (`double(1)`)\cr
    #' Value of the minimum.
    #' @param max (`double(1)`)\cr
    #' Value of the maximum.
    #'
    #' @references
    #' Cauchemez, et al. (2009)
    #' Household Transmission of 2009 Pandemic Influenza A (H1N1) Virus in the United States.
    #' New England Journal of Medicine. 361: 2619–2627. doi:10.1056/NEJMoa0905498
    add_group_rate = function(name,
                              by,
                              prefix = NA,
                              reference = NULL,
                              independent = TRUE,
                              prior = "lnorm",
                              args = list(meanlog = 0, sdlog = 3),
                              rule = "lnorm",
                              rule_scale = 1,
                              min = NA,
                              max = 10) {
      private$.group_rate$add(name = name, by = by, prefix = prefix,
                              reference = reference, independent = independent,
                              prior = prior, args = args, rule = rule,
                              rule_scale = rule_scale, min = min, max = max)
      invisible(self)
    },

    #' @description
    #' Add definition of infectivity profile.
    #'
    #' @param name (`character(1)`)\cr
    #' Name of the infectivity profile.
    #' @param discrete (`logical(1)`)\cr
    #' IF `TRUE` the profile is discrete.
    #' @param args (named `list()`)\cr
    #' Arguments (eventually with types) for the profile.
    #' @param min (`double(1)` | `character(1)`)\cr
    #' One between:
    #' - Scalar value;
    #' - Parameter name (optionally followed by `@P`).
    #' @param max (`double(1)` | `character(1)`)\cr
    #' One between:
    #' - Scalar value;
    #' - Parameter name (optionally followed by `@P`).
    add_infectivity_profile = function(name,
                                       discrete = FALSE,
                                       args = list(),
                                       min = NA,
                                       max = NA) {
      private$.infec_profile$add(name = name, discrete = discrete,
                                 args = args, min = min, max = max)
      invisible(self)
    },

    #' @description
    #' Export instantaneous and cumulative values of the infectivity profile.
    #'
    #' @param opts (`double()`)\cr
    #' Additional values for the computation.
    #' @param xlim (`double(2)`)\cr
    #' Customize the left and/or right limits for the x axis.
    #' To keep the default value of one limit, set it to `NA`.
    #' @param n_bins (`integer(1)`)\cr
    #' Number of bins for continuous profile (ignored if discrete).
    #'
    #' @return
    #' (`data.frame`)\cr
    export_infectivity_profile = function(opts = list(),
                                          xlim = double(),
                                          n_bins = 2000) {
      # internal build function
      private$.pre_build()
      return(private$.model$export_infectivity_profile(opts = opts,
                                                       xlim = xlim,
                                                       n_bins = n_bins))
    },

    ############################################################################
    ############################################################################
    ############################################################################

    #' @description
    #' Run Bayesian MCMC inference and store internally the chain traces.
    #'
    #' @param dataset (`data.frame`)\cr
    #' A data.frame containing the dataset.
    #' @param sample_size (`integer(1)`)\cr
    #' Sample size of a single chain according to posterior distribution.
    #' @param n_steps (`integer(1)`)\cr
    #' Number of MCMC steps.
    #' @param n_burnin (`integer(1)`)\cr
    #' Number of steps to discard at the beginning of the sample.
    #' @param thin (`integer(1)`)\cr
    #' Number of steps to discard between each sample.
    #' @param n_adapt (`integer(1)`)\cr
    #' Number of steps with finite adaptation before sampling.
    #' @param n_chains (`integer(1)`)\cr
    #' Number of independent MCMC chains.
    #' @param n_threads (`integer(1)`)\cr
    #' Number of threads for the parallel computations.
    #' @param rand_init (`logical(1)`)\cr
    #' If `TRUE` the initial value of the parameters is drawn from prior.
    #' @param id_trace (`character(1)`)\cr
    #' Individual identifier to show augmented traces. Default is NA.
    #' @param inv_temper (`double(1)`)\cr
    #' Value of the inverse temperature.
    #' For standard usage this is 1.
    #' @param seed (`integer(1)`)\cr
    #' Seed of the random number generator.
    #' If `NA` uses the internal clock (default).
    #' @param progress (`logical(1)`)\cr
    #' If `FALSE` suppress the inference progress bar.
    infer = function(dataset = NULL,
                     sample_size = NA,
                     n_steps = 10000,
                     n_burnin = 0,
                     thin = 1,
                     n_adapt = 0,
                     n_chains = 1,
                     n_threads = 1,
                     rand_init = TRUE,
                     id_trace = NA,
                     inv_temper = 1.0,
                     seed = NA,
                     progress = TRUE) {
      # internal build function
      private$.build(dataset, sample_size = sample_size,
                     n_steps = n_steps, n_burnin = n_burnin, thin = thin,
                     n_adapt = n_adapt, n_chains = n_chains,
                     n_threads = n_threads, rand_init = rand_init,
                     id_trace = id_trace, inv_temper = inv_temper,
                     seed = seed, progress = progress)

      # infer method from cpp model object
      private$.model$infer()
      invisible(self)
    },

    #' @description
    #' After `$infer` has been run, use the sample to perform simulations.
    #' Extract the expected number of household cases per size
    #' and the credible interval.
    #'
    #' @param cred_level (`double(1)`)\cr
    #' Level for the credible interval.
    #' Value between 0 and 1 (default is 0.95).
    #' @param progress (`logical(1)`)\cr
    #' Progress of simulations. Default is inherited from the model.
    #' @return
    #' (`HousepiOutExpCases`)
    #' An instance of class
    #' [HousepiOutExpCases][housepi::HousepiOutExpCases].
    expec_cases_by_size = function(cred_level = 0.95,
                                   progress = NA) {
      return(HousepiOutExpCases$new(self,
                                    cred_level = cred_level,
                                    progress = progress))
    },

    #' @description
    #' Run the simulation.
    #'
    #' @param dataset (`data.frame`)\cr
    #' A data.frame containing the dataset.
    #' @param n_threads (`integer(1)`)\cr
    #' Number of threads for the parallelization of the loglikelihood.
    #' @param seed (`integer(1)`)\cr
    #' Seed of the random number generator.
    #' If `NA` uses the internal clock (default).
    #' @param progress (`logical(1)`)\cr
    #' If `FALSE` suppress the simulation progress bar.
    simulate = function(dataset = NULL,
                        n_threads = 1,
                        seed = NA,
                        progress = FALSE) {
      # internal build function
      private$.build(dataset = dataset, n_threads = n_threads,
                     seed = seed, progress = progress)

      # call simulate method from internal Rcpp object
      private$.model$simulate()
      invisible(self)
    },

    #' @description
    #' Evaluate the loglikelihood for the given parameters.
    #'
    #' @param dataset (`data.frame`)\cr
    #' A data.frame containing the dataset.
    #' @param n_threads (`integer(1)`)\cr
    #' Number of threads for the parallelization of the loglikelihood.
    #' @param rand_init (`logical(1)`)\cr
    #' If `TRUE` the initial value of parameters is drawn from prior.
    #' Default is `FALSE`.
    #'
    #' @return (`double(1)`)\cr
    #' The value of the loglikelihood.
    evaluate_loglikelihood = function(dataset = NULL,
                                      n_threads = 1,
                                      rand_init = FALSE) {
      # internal build function
      private$.build(dataset = dataset, n_threads = n_threads,
                     rand_init = rand_init)
      # evaluate loglikelihood from the internal Rcpp object
      return(private$.model$evaluate())
    },

    #' @description
    #' Estimate the model evidence (or marginal likelihood).
    #' It uses the method of power posteriors with different available schemes.
    #'
    #' @param dataset (`data.frame`)\cr
    #' A data.frame containing the dataset.
    #' @param scheme (`character(1)`)\cr
    #' If "adaptive" uses adaptive scheme for thermodynamic ladder (default).
    #' Other available scheme is "power-fractions".
    #' @param n_rungs (`integer`)\cr
    #' Number of rungs in the thermodynamic ladder.
    #' @param sample_size (`integer(1)`)\cr
    #' Sample size of a single chain according to posterior distribution.
    #' @param n_steps (`integer(1)`)\cr
    #' Number of MCMC steps.
    #' @param n_burnin (`integer(1)`)\cr
    #' Number of steps to discard at the beginning of the sample.
    #' @param thin (`integer(1)`)\cr
    #' Number of steps to discard between each sample.
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
    #' @return (`double(1)`)\cr
    #' The value of the model evidence.
    #'
    #' @references
    #' Friel N, Hurn M, Wyse J.
    #' Improving power posterior estimation of statistical evidence.
    #' Stat Comput. 2014;24: 709–723. doi:10.1007/s11222-013-9397-1
    estimate_evidence = function(dataset = NULL,
                                 scheme = "adaptive",
                                 n_rungs = 40,
                                 sample_size = NA,
                                 n_steps = 1000,
                                 n_burnin = 0,
                                 thin = 1,
                                 n_adapt = 0,
                                 n_chains = 1,
                                 n_threads = 1,
                                 seed = NA,
                                 progress = TRUE) {
      # internal build function
      private$.build(dataset, sample_size = sample_size, n_steps = n_steps,
                     n_burnin = n_burnin, thin = thin, n_adapt = n_adapt,
                     n_chains = n_chains, n_threads = n_threads,
                     rand_init = TRUE, seed = seed, progress = progress)
      # evaluate marginal from the internal Rcpp object
      return(private$.model$evaluate_marginal(scheme = scheme,
                                              n_rungs = n_rungs))
    },

    ############################################################################
    ############################################################################
    ############################################################################

    #' @description
    #' Use loo package to compute the leave-one-out cross-validation.
    #'
    #' @param ... arguments passed to [loo][loo::loo].
    #' For arguments and further details type ?loo::loo.
    loo = function(...) {
      return(loo::loo(private$.model$household_loglike_matrix, ...))
    },

    #' @description
    #' Use loo package to compute the widely applicable information criterion.
    #'
    #' @param ... arguments passed to [waic][loo::waic].
    #' For arguments and further details type ?loo::waic.
    waic = function(...) {
      return(loo::waic(private$.model$household_loglike_matrix, ...))
    }
  )
)
