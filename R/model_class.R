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
    infec_profile = NULL,

    ############################################################################
    ############################################################################
    ############################################################################
    ############################################################################
    ############################################################################

    #' @description
    #' Assert if parameter definitions match the implementation script.
    #' @return (`logical(1)`)
    #' @noRd
    assert_param = function() {
      # extract requirements from reference file
      # note is.null(c()) -> TRUE
      if (private$model$is_blank) {
        requirements <- c()
      } else {
        refer_file <- private$model$reference()
        requirements <- private$param$require(refer_file)
      }
      requirements <- c(requirements,
                        private$augm$param_require(),
                        private$infec_profile$param_require())
      return(private$param$assert(requirements))
    },

    #' @description
    #' Assert if augmented definitions match the model implementation script.
    #' @return (`logical(1)`)
    #' @noRd
    assert_augm = function() {
      if (private$model$is_blank) {
        return(TRUE)
      } else {
        refer_file <- private$model$reference()
        requirements <- private$augm$require(refer_file)
        return(private$augm$assert(requirements))
      }
    },

    #' @description
    #' Assert group rate definitions from the model implementation script.
    #' @return (`logical(1)`)\cr
    #' @noRd
    assert_group_rate = function() {
      if (private$model$is_blank) {
        return(TRUE)
      } else {
        refer_file <- private$model$reference()
        requirements <- private$group_rate$require(refer_file)
        return(private$group_rate$assert(requirements))
      }
    },

    #' @description
    #' Assert that dataset headers match implementation file.
    #' @param refer_file (`character()`)\cr
    #' Path/to/.refer source file.
    #' @param data_headers (`character()`)\cr
    #' Vector with the dataset headers.
    #' @return (`logical(1)`)\cr
    #' @noRd
    assert_data = function(refer_file, data_headers) {
      # logical state for complete definitions
      is_def_complete <- TRUE

      # vector of regular expression to use for assertion
      regular_expr <- c(
        # from : .data(, "string" or augm(, 'string'
        '.data\\s*\\(\\s*[^,]*\\s*,\\s*["\'](.*?)["\']'
      )
      # retrieve required data headers from regular expressions
      require_set <- template_find_requirement(refer_file, regular_expr)

      if (is.null(require_set)) {
        # if requirement set is NULL it means nothing can be said
        is_def_complete <- NA

      } else {
        # find missing data headers
        missing_set <- template_find_missing(require_set,
                                             data_headers)
        if (length(missing_set) > 0) {
          is_def_complete <- FALSE
          oss <- paste0("Missing data headers(s): ",
                        character_to_string_list(unused_set), ".")
          cli::cli_alert_warning(cli::format_warning(oss))
        }
      }

      return(is_def_complete)
    },

    ############################################################################
    ############################################################################
    ############################################################################
    ############################################################################
    ############################################################################

    #' @description
    #' Pre-build to deal with group rates.
    #' @param dataset (`data.frame`)\cr
    #' The dataset from which to extract group rates.
    #' @noRd
    build_group_rate = function(dataset) {
      # proceed only if group rate definitions are new and complete
      if (isFALSE(
        private$group_rate$cpp_assert(dataset, private$model$reference())
      )) {
        # NOTE: at the moment I just cache group rate changes
        # search for the cache file
        source_file <- file.path(private$model$cache_dir,
                                 "src", "model.cpp")
        # deal with the case cache has been deleted
        if (!file.exists(source_file)) {
          # check if it's a template
          model_choice <- private$model$path()
          if (file.exists(file.path(model_choice, "model.cpp"))) {
            # re-assign model object forcing copy of source to cache
            private$model <- HousepiGeneratorClass$new(model_choice,
                                                       force = TRUE)
            # assign new cache file
            source_file <- file.path(private$model$cache_dir,
                                     "src", "model.cpp")
          } else {
            stop(cli::format_error(
              "Unknown error: impossible to find source file."
            ))
          }
        }

        # write the cpp script and save new reference
        private$group_rate$cpp_write(dataset, source_file)
        private$model$save_reference()
      }

      # append new parameter definitions
      new_param_defs <- list()
      for (rate_id in names(private$group_rate$logic_table)) {
        for (par_id in private$group_rate$logic_table[[rate_id]]$.param_name) {
          if (par_id != "REFERENCE") {
            new_param_defs[[par_id]] <- private$group_rate$param_defs[[rate_id]]
          }
        }
      }
      if (length(new_param_defs) > 0) {
        self$load_param_list(new_param_defs, append = TRUE)
      }

      return(self)
    },

    #' @description
    #' Associate dataset to the model and options for the setup.
    #' @param dataset (`character(1)` | `data.frame`)\cr
    #' Path to a csv file or a data.frame containing the dataset.
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
    #' @noRd
    build = function(dataset = NULL,
                     n_steps = 1000,
                     n_burnin = 0,
                     thin = 1,
                     n_adapt = 0,
                     n_chains = 1,
                     n_threads = 1,
                     rand_init = TRUE,
                     id_trace = NA,
                     inv_temper = 1.0,
                     seed = NA) {

      #----- MODEL PRE-BUILD -----#

      if (private$model$is_blank) {
        # blank dataset
        # note: no assert definitions in such case
        id_trace <- "HH_1_1"
        dataset <- data.frame(id = c(id_trace),
                              house_id = c("HH_1"),
                              status = c(1))

      } else {
        # operations on non-null dataset
        if (!is.null(dataset)) {
          # check that dataset type is correct
          if (!is.data.frame(dataset)) {
            stop(cli::format_error(c(
              "Unsupported type for dataset.",
              i = "Dataset should be a data.frame."
            )))
          }

          # deal with group rates
          private$build_group_rate(dataset)

          # assert dataset headers
          if (isFALSE(private$assert_data(private$model$reference(),
                                          colnames(dataset)))) {
            stop("Incomplete dataset headers.")
          }
        }
      }

      # assert parameter definitions
      if (isFALSE(private$assert_param())) {
        stop("Incomplete parameter definitions.")
      }
      # assert augmented variable definitions
      if (isFALSE(private$assert_augm())) {
        stop("Incomplete augmented variable definitions.")
      }

      #----- MODEL GENERATION -----#
      private$model$generate()

      #----- ASSOCIATE DEFINITIONS -----#

      # change scientific notation
      store_scipen <- options()$scipen
      options(scipen = 100)

      #----- 1) SETUP -----#
      load_setup <- list()
      load_setup[["n_steps"]] <- n_steps
      load_setup[["n_burnin"]] <- n_burnin
      load_setup[["thin"]] <- thin
      load_setup[["n_adapt"]] <- n_adapt
      load_setup[["n_chains"]] <- n_chains
      load_setup[["n_threads"]] <- n_threads
      load_setup[["rand_init"]] <- rand_init
      load_setup[["id_trace"]] <- id_trace
      load_setup[["inv_temper"]] <- inv_temper
      load_setup[["seed"]] <- seed

      private$model$rcpp$import_setup(
        lapply(load_setup, as.character)
      )

      #----- 2) PARAMETERS -----#
      private$model$rcpp$import_params(
        data.frame(lapply(private$param$defs, as.character))
      )

      #----- 3) DATA -----#
      if (!is.null(dataset)) {
        private$model$rcpp$import_data(
          data.frame(lapply(dataset, as.character))
        )
      }

      #----- 4) AUGMENTED -----#
      if (!is.null(private$augm$defs)) {
        private$model$rcpp$import_augm(
          data.frame(lapply(private$augm$defs, as.character))
        )
      }

      #----- 5) INFECTIVITY PROFILE -----#
      if (!is.null(private$infec_profile$defs)) {
        private$model$rcpp$import_infec_profile(
          lapply(private$infec_profile$defs, as.character)
        )
      }

      # restore scientific notation option
      options(scipen = store_scipen)
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
    #' New instance of class [HouseholdModelClass][housepi::HouseholdModelClass].
    #' @param model (`character(1)`)\cr
    #' One between:
    # - template name;
    # - path/to/folder with implementation;
    # - path/to/folder with previously compiled model;
    # - "blank": template name for blank model (default).
    #' @param quiet (`logical(1)`)\cr
    #' If `TRUE` hide compilation commands.
    #' @param force (`logical(1)`)\cr
    #' If `TRUE` compile even if a cached model already exists.
    #' @param as_template (`logical(1)`)\cr
    #' If `TRUE` store the dynamic library as a local template.
    #' IF `FALSE` recompile the model at every new session.
    initialize = function(model = NA,
                          quiet = TRUE,
                          force = FALSE,
                          as_template = FALSE) {
      if (is.na(model)) {
        model <- "blank"
        cli::cli_alert_warning(cli::format_warning(
          "No model specified. Falling back to default \"blank\" model."
        ))
      }

      # assign model definitions
      private$model <- HousepiGeneratorClass$new(model = model,
                                                 quiet = quiet,
                                                 force = force,
                                                 as_template = as_template)
      # initialize internal fields
      self$drop_all()

      invisible(self)
    },

    #' @description
    #' Show model description in RStudio Help window.
    help = function() {
      templ_info_file <- file.path(private$model$path(), "template_info.Rd")
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
    #' Run Bayesian MCMC inference and store internally the chain traces.
    #' @param dataset (`character(1)` | `data.frame`)\cr
    #' Path to a csv file or a data.frame containing the dataset.
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
      private$build(dataset,
                    n_steps = n_steps, n_burnin = n_burnin, thin = thin,
                    n_adapt = n_adapt, n_chains = n_chains,
                    n_threads = n_threads, rand_init = rand_init,
                    id_trace = id_trace, inv_temper = inv_temper,
                    seed = seed)

      # infer method from cpp model object
      private$model$rcpp$infer(progress)
      invisible(self)
    },

    #' @description
    #' Run the simulation.
    #' @param dataset (`character(1)` | `data.frame`)\cr
    #' Path to a csv file or a data.frame containing the dataset.
    #' @param n_threads (`integer(1)`)\cr
    #' Number of threads for the parallelization of the loglikelihood.
    #' @param rand_init (`logical(1)`)\cr
    #' If `TRUE` the initial value of the parameters is drawn from prior.
    #' Default is `FALSE`.
    #' @param seed (`integer(1)`)\cr
    #' Seed of the random number generator.
    #' If `NA` uses the internal clock (default).
    #' @param progress (`logical(1)`)\cr
    #' If `FALSE` suppress the simulation progress bar.
    simulate = function(dataset = NULL,
                        n_threads = 1,
                        rand_init = FALSE,
                        seed = NA,
                        progress = FALSE) {
      # internal build function
      private$build(dataset = dataset, n_threads = n_threads,
                    rand_init = rand_init, seed = seed)
      # call simulate method from internal Rcpp object
      private$model$rcpp$simulate(progress)
      invisible(self)
    },

    #' @description
    #' Estimate the loglikelihood for the given parameters.
    #' @param dataset (`character(1)` | `data.frame`)\cr
    #' Path to a csv file or a data.frame containing the dataset.
    #' @param n_threads (`integer(1)`)\cr
    #' Number of threads for the parallelization of the loglikelihood.
    #' @param rand_init (`logical(1)`)\cr
    #' If `TRUE` the initial value of parameters is drawn from prior.
    #' Default is `FALSE`.
    #' @param inv_temper (`double(1)`)\cr
    #' Value of the inverse temperature in the range [0, 1].
    #' For standard usage this is 1.
    #' @return (`double(1)`)\cr
    #' The value of the loglikelihood.
    loglikelihood = function(dataset = NULL,
                             n_threads = 1,
                             rand_init = FALSE,
                             inv_temper = 1.0) {
      # internal build function
      private$build(dataset = dataset, n_threads = n_threads,
                    rand_init = rand_init, inv_temper = inv_temper)
      # evaluate loglikelihood from the internal Rcpp object
      return(private$model$rcpp$evaluate())
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
                        progress = TRUE) {
      # internal build function
      private$build(dataset, n_steps = n_steps, n_burnin = n_burnin,
                    thin = thin, n_adapt = n_adapt, n_chains = n_chains,
                    n_threads = n_threads, rand_init = TRUE, seed = seed)
      # evaluate marginal from the internal Rcpp object
      output <- private$model$rcpp$evaluate_marginal(n_rungs, scheme, progress)
      return(output)
    },

    ############################################################################
    ############################################################################
    ############################################################################
    ############################################################################
    ############################################################################

    #' @description
    #' Add a new parameter definition.
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
      private$param$add(name = name, value = value, constant = constant,
                        discrete = discrete, prior = prior, args = args,
                        rule = rule, rule_scale = rule_scale,
                        min = min, max = max)
      invisible(self)
    },

    #' @description
    #' Add a new augmented variable definition.
    #' @param name (`character(1)`)\cr
    #' Name of the augmented variable.
    #' @param status (`integer(1)`)\cr
    #' The associated status of the definition.
    #' Default NA is for all infectious.
    #' @param value (`double(1)` | `character(1)`)\cr
    #' One between:
    #' - Scalar value;
    #' - Parameter name (followed by `@P`);
    #' - Data name (followed by `@D`).
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
    #' - Parameter name (followed by `@P`);
    #' - Data name (followed by `@D`).
    #' @param max (`double(1)` | `character(1)`)\cr
    #' Maxmimum tvalue. One between:
    #' - Scalar value;
    #' - Parameter name (followed by `@P`);
    #' - Data name (followed by `@D`).
    #' @param init_rand (`character(1)`)\cr
    #' If `TRUE` draw initial vale from model distribution.
    add_augm = function(name = NA,
                        status = NA,
                        value = NA,
                        discrete = FALSE,
                        model = NA,
                        args = list(),
                        rule = NA,
                        rule_scale = NA,
                        rule_from = NA,
                        min = NA,
                        max = NA,
                        init_rand = FALSE) {
      private$augm$add(name = name, status = status,
                       value = value, discrete = discrete,
                       model = model, args = args, rule = rule,
                       rule_scale = rule_scale, rule_from = rule_from,
                       min = min, max = max, init_rand = init_rand)
      invisible(self)
    },

    #' @description
    #' Add definitions of a group rate by given data feature(s). 
    #' Reference value has unitary rate and can be specified (default is most represented group).
    #' Associated parameter definition(s) with default prior "lnorm" (meanlog = 0, sdlog = 3).
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
      private$group_rate$add(name = name, by = by, prefix = prefix,
                             reference = reference, independent = independent,
                             prior = prior, args = args, rule = rule,
                             rule_scale = rule_scale, min = min, max = max)
      invisible(self)
    },

    #' @description
    #' Add definition of infectivity profile.
    #' @param name (`character(1)`)\cr
    #' Name of the infectivity profile.
    #' @param discrete (`logical(1)`)\cr
    #' IF `TRUE` the profile is discrete.
    #' @param args (named `list()`)\cr
    #' Arguments (eventually with types) for the profile.
    #' @param min (`double(1)` | `character(1)`)\cr
    #' One between:
    #' - Scalar value;
    #' - Parameter name (followed by `@P`).
    #' @param max (`double(1)` | `character(1)`)\cr
    #' One between:
    #' - Scalar value;
    #' - Parameter name (followed by `@P`).
    add_infectivity_profile = function(name,
                                       discrete = FALSE,
                                       args = list(),
                                       min = NA,
                                       max = NA) {
      private$infec_profile$add(
        name = name, discrete = discrete,
        args = args, min = min, max = max)
      invisible(self)
    },

    ############################################################################
    ############################################################################
    ############################################################################
    ############################################################################
    ############################################################################

    #' @description
    #' Associate parameters from a list.
    #' @param definitions (named `list()`)\cr
    #' A list containg parameter names with respective list of definitions.
    #' @param append (`logical(1)`)\cr
    #' If `TRUE` appends the list to current definition (default).
    #' If `FALSE` drop all other parameter definitions.
    load_param_list = function(definitions, append = TRUE) {
      private$param$from_list(definitions, append = append)
      invisible(self)
    },

    #' @description
    #' Associate augmented variables from a list.
    #' @param definitions (named `list()`)\cr
    #' Augmented variable names with respective list of definitions.
    #' @param append (`logical(1)`)\cr
    #' If `TRUE` appends the list to current definition (default).
    #' If `FALSE` drop all other augmented variable definitions.
    load_augm_list = function(definitions, append = TRUE) {
      private$augm$from_list(definitions, append = append)
      invisible(self)
    },

    #' @description
    #' Add missing paramters from template definitions.
    #' @param append
    #' If `TRUE` append only missing definitions (default).
    #' Else drop all previous definitions
    template_param = function(append = TRUE) {
      if (private$model$is_blank) {
        cli::cli_alert_warning(cli::format_warning(
          "No template parameters available for blank model."
        ))
      } else {
        # retrieve the template parameters file
        private$param$template(
          file.path(private$model$path(), "template_params.R"),
          private$param$require(private$model$reference()), 
          isTRUE(append)
        )
      }
      invisible(self)
    },

    #' @description
    #' Add missing augmented variables from template definitions.
    #' @param append
    #' If `TRUE` append only missing definitions (default).
    #' Else drop all previous definitions
    template_augm = function(append = TRUE) {
      if (private$model$is_blank) {
        cli::cli_alert_warning(cli::format_warning(
          "No template augmented variables available for blank model."
        ))
      } else {
        # retrieve the template augmented file
        private$augm$template(
          file.path(private$model$path(), "template_augm.R"),
          private$augm$require(private$model$reference()), 
          isTRUE(append)
        )
      }

      invisible(self)
    },

    ############################################################################
    ############################################################################
    ############################################################################
    ############################################################################
    ############################################################################

    #' @description
    #' Export the current dataset at current state.
    #' If the case, it also return the associated augmented variables.
    #' @return (`data.frame`)\cr
    #' A data.frame with individual level information.
    export_data = function() {
      if (is.null(private$model$rcpp)) {
        stop(cli::format_error(c(
          "Impossible to export data before model build.",
        )))
      }
      return(private$model$rcpp$export_data())
    },

    #' @description
    #' Export the traces of last inference call as a coda object.
    #' @param keep_loglike (`logical(1)`)\cr
    #' If `FALSE` don't export the trace of the loglikelihood.
    #' Default is `TRUE`.
    #' @return
    #' (`coda::mcmc` | `coda::mcmc.list`)
    export_traces = function(keep_loglike = TRUE) {

      # load traces from last inference
      df <- NULL
      is_wrong_usage <- is.null(private$model$rcpp)
      if (!is_wrong_usage) {
        df <- private$model$rcpp$export_traces()
        is_wrong_usage <- (dim(df)[1] == 0)
      }

      if (is_wrong_usage) {
        stop(cli::format_error(c(
          "MCMC traces weren't found.",
          i = "Be sure to have run your inference first (method {.code infer()})."
        )))
      }

      # drop loglikelihood if flag is true
      if (isFALSE(keep_loglike)) {
        df <- df[, names(df) != "loglikelihood"]
      }

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
        return(coda::as.mcmc.list(df_list))

      } else {
        return(coda::mcmc(
          data = df,
          start = load_setup[["n_burnin"]] + load_setup[["thin"]],
          end = load_setup[["n_steps"]],
          thin = load_setup[["thin"]]
        ))
      }
    },

    #' @description
    #' Export the parameter definitions at current state.
    #' @return (`data.frame`)\cr
    #' A data.frame with formatted parameter definitions.
    export_param = function() {
      return(private$param$defs)
    },

    #' @description
    #' Export the augmented variable definitions at current state.
    #' @return (`data.frame`)\cr
    #' A data.frame with formatted augmented variable definitions.
    export_augm = function() {
      return(private$augm$defs)
    },

    #' @description
    #' Export instantaneoud and cumulative values of the infectivity profile.
    #' @param opts (`double()`)\cr
    #' Additional values for the computation.
    #' @param xlim (`double(2)`)\cr
    #' Customize the left and/or right limits for the x axis.
    #' To keep the default value of one limit, set it to `NA`.
    #' @return @return (`data.frame`)\cr
    export_infectivity_profile = function(opts = list(),
                                          xlim = double()) {
      # internal build function
      private$build()
      
      # load options
      opts_ <- as.list(as.numeric(opts))

      # load x-limits
      xlim_ <- list()
      xlim <- unlist(xlim)
      if (length(xlim) == 2) {
        if (!all(is.na(xlim))) {
          # check 
          if (!any(is.na(xlim))) {
            if (xlim[[1]] == xlim[[2]]) {
              stop(cli::format_error(
                "Left x limit larger than right limit."
              ))
            } else if (xlim[[1]] > xlim[[2]]) {
              stop(cli::format_error(
                "Left x limit equal to right limit."
              ))
            }
          }
          xlim_ <- lapply(xlim, as.character)
          names(xlim_) <- c("min", "max")
        }
      } else if (length(xlim) > 0) {
        stop(cli::format_error(
          "Provide exactly 2 limits for the x axis."
        ))
      }

      return(private$model$rcpp$export_infec_profile(opts_, xlim_))
    },

    #' @description
    #' Print a summary of the current dataset.
    summarise_data = function() {
      if (is.null(private$model$rcpp)) {
        stop(cli::format_error(c(
          "Impossible to show data summary before model build.",
        )))
      }
      return(private$model$rcpp$show_summary())
    },

    #' @description
    #' Export the traces of last inference call as a coda object.
    #' @return
    #' (`data.frame`)
    summarise_traces = function() {
      res <- summary(self$export_traces())
      imask <- private$param$defs["constant"] == FALSE
      inames <- private$param$defs[imask, "name"]
      output <- data.frame(
        Mean = res[["statistics"]][inames, "Mean"],
        SD = res[["statistics"]][inames, "SD"],
        Q2.5 = res[["quantiles"]][inames, "2.5%"],
        Q97.5 = res[["quantiles"]][inames, "97.5%"]
      )
      return(output)
    },

    ############################################################################
    ############################################################################
    ############################################################################
    ############################################################################
    ############################################################################

    #' @description
    #' Drop all current parameter definitions.
    drop_all_param = function() {
      private$param$drop_all()
      invisible(self)
    },

    #' @description
    #' Drop all current augmented variable definitions.
    drop_all_augm = function() {
      private$augm$drop_all()
      invisible(self)
    },

    #' @description
    #' Drop all current augmented variable definitions.
    drop_all_group_rate = function() {
      private$group_rate$drop_all()
      invisible(self)
    },

    #' @description
    #' Drop all current definitions.
    drop_all = function() {
      # initialize parameter helper
      private$param <- HousepiParameterClass$new()
      # initialize augmented helper
      private$augm <- HousepiAugmentedClass$new()
      # initialize group rate helper
      private$group_rate <- HousepiGroupRateClass$new()
      # initialize infectivity profile helper
      private$infec_profile <- HousepiInfecProfileClass$new()
      invisible(self)
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
    }
  )
)
