#      __
#     / /_  ____  __  __ ____ __  ____  __
#    / __ \/ __ \/ / / / ___/ _ \/ __ \/ /
#   / / / / /_/ / /_/ (__  )  __/ /_/ / /
#  /_/ /_/\____/\____/____/\___/ ____/_/
#                             /_/
#
# Copyright (C) 2025 CAMAGLIA Francesco - Institut Pasteur, Paris (France).
#
# model_class.R
#

############################################################################
############################################################################
############################## Household Model #############################
############################################################################
############################################################################

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
  lock_class = TRUE,

  ############################################################################
  ################################## PRIVATE #################################
  ############################################################################

  private = list(
    .model = NULL,
    .param = NULL,
    .augm = NULL,
    .group_rate = NULL,
    .infec_profile = NULL,
    .cri_level = 0.95,

    # @description
    # Dataset headers required by implementation file and other definitions.
    #
    # @return `character()`
    .require_dataset = function() {
      requirements <- c()
      if (!private$.model$is_blank) {
        requirements <- private$.model$data_required
      }
      return(c(
        requirements, # by the model
        private$.augm$data_required # by augmented variable
      ))
    },

    # @description
    # Assert if dataset headers match requirements.
    #
    # @param data_headers `character()`\cr
    # Vector with the dataset headers.
    .assert_dataset = function(data_headers) {
      if (length(data_headers) == 0) {
        stop_h("Impossible to proceed. No database was provided.")
      }
      private$.model$assert(data_headers, private$.require_dataset())
      invisible(self)
    },

    # @description
    # Parameters required by implementation file and other definitions.
    #
    # @return `character()`
    .require_param = function() {
      if (private$.model$is_blank) {
        requirements <- self$param_defs$name
      } else {
        requirements <- private$.param$require(private$.model$reference)
      }
      return(c(
        requirements, # by the model
        private$.augm$param_required, # by augmented variable
        private$.infec_profile$param_required # by infectivity profile
      ))
    },

    # @description
    # Assert if parameter definitions match requirements
    .assert_param = function() {
      private$.param$assert(private$.require_param())
      invisible(self)
    },

    # @description
    # Augmented variables required by implementation file and other definitions.
    #
    # @return `character()`
    .require_augm = function() {
      return(c(
        private$.augm$require(private$.model$reference), # by the model
        private$.infec_profile$augm_required # by infectivity profile
      ))
    },

    # @description
    # Assert if augmented definitions match the requirements.
    .assert_augm = function() {
      if (!private$.model$is_blank) {
        private$.augm$assert(private$.require_augm())
      }
      invisible(self)
    },

    # @description
    # Group rates required by implementation file and other definitions.
    #
    # @return `character()`
    .require_group_rate = function() {
      return(private$.group_rate$require(private$.model$reference))
    },

    # @description
    # Assert group rate definitions from the model implementation script.
    .assert_group_rate = function() {
      if (!private$.model$is_blank) {
        private$.group_rate$assert(private$.require_group_rate())
      }
      invisible(self)
    },

    ############################################################################
    ############################################################################
    ############################################################################

    # @description
    # Enable multiplicity engine.
    #
    # @param dataset `data.frame`\cr
    # A data.frame containing the dataset.
    .multiplicity_engine = function(dataset = NULL) {
      # check if model allows for multiplicity
      refer_file <- private$.model$reference
      regular_expr <- ".multiplicity\\s*\\(\\s*\\)"

      if (template_assert_expression(refer_file, regular_expr)) {

        # note: this might be a modeler error
        if (length(private$.require_augm()) > 0) {
          warning_h(
            "Model with data augmentation called multiplicity engine."
          )
        }

        # if dataset not provided, load previous
        # UPDATE: this "back & forth" of dataset is not optimal!
        if ((length(dataset) == 0) && (length(private$.model$last_data_cols) > 0)) {
          dataset <- self$dataset
        }

        headers <- c(
          # dataset headers required by implementation
          private$.require_dataset(),
          # dataset headers required by group rates
          unlist(lapply(
            private$.group_rate$logic_tables, 
            function(x) unlist(lapply(x, names))
          ))
        )
        # remove logic tables auxiliary column names
        headers <- setdiff(unique(headers), c(".counts", ".param.name"))

        # list of splitted dataset by house_id
        split_data <- split(dataset[c("status", headers)], dataset$house_id)
        # sort individuals in each household by all columns
        sorted_data <- lapply(split_data, function(x) {
            x[do.call(order, x), , drop = FALSE]
        })
        # create hash strings from sorted rows
        hashes <- sapply(sorted_data, function(x) {
          paste(apply(x, 1, paste, collapse = "-"), collapse = "|")
        })
        # regroup households by hash
        hash_groups <- split(names(hashes), hashes)

        # create summary table with class_id and multiplicity
        multi_hh_table <- data.frame(class_id = sapply(hash_groups, `[`, 1),
                                     multiplicity = lengths(hash_groups),
                                     stringsAsFactors = FALSE)

        # use hashes to produce multiplicity table
        output <- data.frame(house_id = unique(dataset$house_id),
                             multiplicity = 0,
                             class_id = NA,
                             stringsAsFactors = FALSE)
        rownames(output) <- output$house_id
        output[multi_hh_table$class_id, "multiplicity"] <- multi_hh_table$multiplicity

        for (i in seq_len(nrow(multi_hh_table))) {
          output[hash_groups[[i]], "class_id"] <- multi_hh_table$class_id[i]
        }
        # call rcpp method
        private$.model$import_multiplicity(output)
      }

      invisible(self)
    },

    ############################################################################
    ############################################################################
    ############################################################################

    # @description
    # Pre-build to deal with group rates.
    #
    # @param dataset `data.frame`\cr
    # The dataset from which to extract group rates.
    # @param .hold `logical(1)`\cr
    # If `FALSE` does not ultimate model build.
    # For internal use.
    .build_group_rate = function(dataset = NULL, .hold = TRUE) {
      # if dataset not provided, load previous
      # UPDATE: this "back & forth" of dataset is not optimal!
      if ((length(dataset) == 0) && (length(private$.model$last_data_cols) > 0)) {
        dataset <- self$dataset
      }

      # assert and implement group rates
      is_updated <- private$.group_rate$cpp_assert(
        dataset, private$.model$reference, private$.model$factor_info
      )

      # proceed only if group rate definitions have been updated
      if (isFALSE(.hold)) {
        if (isFALSE(is_updated)) {
          # force caching a new model
          # UPDATE: can't I change in-place?

          # new model name
          new_name <- private$.model$name
          gr_label <- "+group+rates"
          if (!stringr::str_detect(new_name, stringr::fixed(gr_label))) {
            new_name <- paste0(new_name, gr_label)
          }
          # re-assign model object forcing copy of source to cache
          private$.model <- HousepiWrapper$new(
            private$.model$source_dir, as_template = FALSE,
            quiet = private$.model$quiet, force = TRUE,
            name = new_name
          )

          # write the cpp script for new group rates in the new cache file
          private$.group_rate$cpp_write(
            file.path(private$.model$cache_dir, "src", "model.cpp")
          )
          # save the new reference
          private$.model$save_reference()
        }
      }

      # append parameter definitions if missing
      new_param_defs <- private$.group_rate$param_group_rate_defs
      # find missing definitions
      mask <- setdiff(names(new_param_defs), self$param_defs$name)
      new_param_defs <- new_param_defs[mask]
      if (length(new_param_defs) > 0) {
        self$load_param_list(from = new_param_defs, append = TRUE)
      }

      # perform standard generation and data importation
      private$.model$generate()
      private$.model$import_dataset(dataset)

      invisible(self)
    },

    # @description
    # Import only parameters and infectivity profile for a minimal build.
    .build_minimal = function() {
      # assert parameter definitions before importing
      private$.assert_param()
      private$.model$import_param(self$param_defs)

      # infectivity profile
      if (length(private$.infec_profile$defs) > 0) {
        private$.model$import_infec_profile(private$.infec_profile$defs)
      }

      invisible(self)
    },

    # @description
    # Full model build with complete list of definitions.
    #
    # @param dataset `data.frame`\cr
    # A data.frame containing the dataset.
    # @param sample_size `integer(1)`\cr
    # Sample size of a single chain according to posterior distribution.
    # @param n_steps `integer(1)`\cr
    # Number of MCMC steps.
    # @param n_burnin `integer(1)`\cr
    # Number of steps to discard at the beginning of the sample.
    # @param thin `integer(1)`\cr
    # Number of steps to discard between each sample.
    # @param n_adapt `integer(1)`\cr
    # Number of steps with finite adptation before sampling.
    # @param n_chains `integer(1)`\cr
    # Number of independent MCMC chains.
    # @param n_threads `integer(1)`\cr
    # Number of threads for the parallelization of the loglikelihood.
    # @param rand_init `logical(1)`\cr
    # If `TRUE` the initial value of the parameters is drawn from prior.
    # @param id_trace `character(1)`\cr
    # Individual identifier to show augmented traces. Default is NA.
    # @param inv_temper `double(1)`\cr
    # Value of the inverse temperature.
    # For standard usage this is 1.
    # @param seed `integer(1)`\cr
    # Seed of the random number generator.
    # If `NA` uses the internal clock (default).
    # @param progress `logical(1) | integer(1)`\cr
    # If `TRUE` prints the progress bar.
    .build_full = function(dataset = NULL,
                           sample_size = NA, n_steps = 0,
                           n_burnin = 0, thin = 1, n_adapt = 0,
                           n_chains = 1, n_threads = 1,
                           rand_init = TRUE, id_trace = NA, inv_temper = 1.0,
                           seed = NA, progress = FALSE) {
      # import dataset
      self$import_dataset(dataset, .hold = FALSE)
      # run model pre-build
      private$.build_minimal()
      # assert augmented variable definitions
      private$.assert_augm()

      # individual id for blank model augmented traces
      if (private$.model$is_blank && is.na(id_trace)) id_trace <- 1
      # control that seed is 32bit
      if (isTRUE(abs(seed) > .Machine$integer.max)) {
        stop_h(
          "Supplied seed is not supported.",
          sprintf("Maximum value for a integer in R is %.0f.",
                  .Machine$integer.max)
        )
      }
      # import setup
      private$.model$import_setup(list(
        sample_size = ifelse(is.na(sample_size), 0, sample_size),
        n_steps = n_steps, thin = thin,
        n_burnin = n_burnin, n_adapt = n_adapt, n_chains = n_chains,
        n_threads = n_threads, rand_init = rand_init,
        id_trace = id_trace, inv_temper = inv_temper,
        seed = seed, progress = as.integer(progress)
      ))

      # build group & contact rates
      # UPDATE: rebuild group rates <==> something has changed
      private$.model$build_all_rates()

      # build augmented variables
      if (length(self$augm_defs) > 0) {
        private$.model$import_augm(self$augm_defs)
      }

      invisible(self)
    }
  ),

  ############################################################################
  ################################## ACTIVE ##################################
  ############################################################################

  active = list(
    #' @field name `character(1)`\cr
    #' Model name.
    name = function(value) {
      if (missing(value)) return(private$.model$name)
      stop_h("'name' is read only.")
    },

    #' @field dataset `data.frame`\cr
    #' Dataset with individual level information of current dataset.
    #' If the case, it also return the associated augmented variables.
    dataset = function(value) {
      if (missing(value)) return(private$.model$dataset)
      else self$import_dataset(value, .hold = TRUE)
    },

    #' @field summary `data.frame`\cr
    #' A summary of statuses in the dataset.
    summary = function(value) {
      if (missing(value)) return(private$.model$summary)
      stop_h("'summary' is read only.")
    },

    #' @field summary_house_size `data.frame`\cr
    #' A summary of house sizes in the dataset.
    summary_house_size = function(value) {
      if (missing(value)) return(private$.model$summary_house_size)
      stop_h("'summary_house_size' is read only.")
    },

    #' @field setup `list()`\cr
    #' Last setup used by the build process.
    setup = function(value) {
      if (missing(value)) return(private$.model$setup)
      stop_h("'setup' is read only.")
    },

    #' @field param_defs `data.frame`\cr
    #' Parameter formatted definitions.
    param_defs = function(value) {
      if (missing(value)) return(private$.param$defs)
      else self$load_param_list(from = value, append = FALSE)
    },

    #' @field augm_defs `data.frame`\cr
    #' Augmented variable formatted definitions.
    augm_defs = function(value) {
      if (missing(value)) return(private$.augm$defs)
      else self$load_augm_list(from = value, append = FALSE)
    },

    #' @field infec_profile_defs `list`\cr
    #' Infectivity profile definition list.
    infec_profile_defs = function(value) {
      if (missing(value)) return(private$.infec_profile$defs)
      else do.call(self$add_infec_profile, value)
    },

    #' @field group_rate_defs `list`\cr
    #' List containing current group rate definitions.
    group_rate_defs = function(value) {
      if (missing(value)) return(private$.group_rate$defs)
      else self$load_group_rate_list(from = value, append = FALSE)
    },

    #' @field group_rate_tables `list`\cr
    #' List containing table of parameters associated to each group rate.
    group_rate_tables = function(value) {
      if (missing(value)) {
        private$.build_group_rate(NULL, .hold = TRUE)
        return(private$.group_rate$logic_tables)
      }
      stop_h("'group_rate_tables' is read only.")
    },

    #' @field house_loglike_mtx `matrix()`\cr
    #' Traces of household loglikelihoods.
    house_loglike_mtx = function(value) {
      if (missing(value)) return(private$.model$house_loglike_mtx)
      stop_h("'house_loglike_mtx' is read only.")
    },

    #' @field traces `matrix()` | `list()` of `matrix()`\cr
    #' Traces of last inference.
    traces = function(value) {
      if (missing(value)) return(private$.model$traces)
      stop_h("'traces' is read only.")
    },

    #' @field traces_no_loglike `matrix()` | `list()` of `matrix()`\cr
    #' Traces of last inference without loglikelihood.
    traces_no_loglike = function(value) {
      if (missing(value)) return(private$.model$traces_no_loglike)
      stop_h("'traces_no_loglike' is read only.")
    },

    #' @field traces_results `data.frame`\cr
    #' Results extracted from the traces of last inference.
    traces_results = function(value) {
      if (missing(value)) {
        return(retrieve_traces_results(
          private$.model$traces_no_loglike,
          cri_level = private$.cri_level
        ))
      }
      stop_h("'traces_results' is read only.")
    },

    #' @field cri_level `double(1)`\cr
    #' Level for the credible interval.
    #' Value between 0.5 and 1 (default is 0.95).
    cri_level = function(value) {
      if (missing(value)) {
        return(private$.cri_level)
      } else if ((value <= 0.5) || (value >= 1)) {
        stop_h(
          "Credible interval level must be a value between 0.5 and 1."
        )
      } else {
        private$.cri_level <- value
      }
    }
  ),

  ############################################################################
  ################################## PUBLIC ##################################
  ############################################################################

  public = list(
    #' @description
    #' New instance of class
    #' [HouseholdModel][housepi::HouseholdModel].
    #'
    #' @param choice `character(1)`\cr
    #' Choose a model providing one between:
    #' - template name;
    #' - path/to/folder with implementation;
    #' - path/to/folder with previously compiled model;
    #' - "blank": template name for blank model (default).
    #' @param as_template `logical(1)`\cr
    #' If `TRUE` store the dynamic library as a local template.
    #' IF `FALSE` recompile the model at every new session.
    #' @param force `logical(1)`\cr
    #' If `TRUE` compile even if a cached model already exists.
    #' @param name `chracter(1)`\cr
    #' Custom name of the model (default is inferred from `choice`).
    #' @param cri_level `double(1)`\cr
    #' Level for the credible interval.
    #' Value between 0 and 1 (default is 0.95).
    #' @param quiet `logical(1)`\cr
    #' If `TRUE` hide `pkgbuild` compilation commands (default).
    initialize = function(choice = NA, as_template = FALSE, force = FALSE,
                          name = NA, cri_level = 0.95, quiet = TRUE) {
      # generate model
      private$.model <- HousepiWrapper$new(
        choice = choice, as_template = as_template, force = force,
        name = name, quiet = quiet
      )$generate()

      # credible interval level
      self$cri_level <- cri_level

      # initialize internal fields
      self$drop_all()

      invisible(self)
    },

    #' @description
    #' Override the default printing behavior.
    print = function() {
      cli::cli_alert_info("Household model {.pkg {self$name}}.")
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

      } else if (interactive()) {
        # print as plain text in the console
        txt <- capture.output(tools::Rd2txt(rd_file))
        cat(txt, sep = "\n")

      } else {
        warning_h("Help method not supported in non-interactive mode.")
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
      for (p in names(private$.group_rate$param_group_rate_defs)) {
        private$.param$remove(p, show_warnings = FALSE)
      }
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
    #' Import a new dataset.
    #'
    #' @param dataset `data.frame`\cr
    #' A data.frame containing the dataset with individual information by row.
    #' @param .hold `logical(1)`\cr
    #' If `FALSE` does not ultimate model build. Default is `TRUE`.
    #' For internal use.
    import_dataset = function(dataset = NULL, .hold = TRUE) {

      if (private$.model$is_blank) {
        if (length(dataset) > 0) {
          # UPDATE: actually it could be nice to allow for this
          warning_h(
            "Impossible to provide a dataset with blank model.",
            "Falling back to default dataset."
          )
        }
        # blank dataset (id_trace <- 1)
        dataset <- data.frame(id = c(1),
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

        # assert dataset headers
        private$.assert_dataset(colnames(dataset))

      } else {
        # assert headers of last imported dataset
        private$.assert_dataset(private$.model$last_data_cols)
      }

      if (length(private$.group_rate$defs) > 0) {
        # special generation dealing with group rates
        private$.build_group_rate(dataset, .hold = .hold)

      } else {
        # default model generation and dataset import
        private$.model$generate()
        private$.model$import_dataset(dataset)
      }

      invisible(self)
    },

    ############################################################################
    ############################################################################
    ############################################################################

    #' @description
    #' Add a new parameter definition.
    #'
    #' @param name `character(1)`\cr
    #' Name of the parameter.
    #' @param value `double(1)`\cr
    #' Initial value of the parameter.
    #' @param const `logical(1)`\cr
    #' If `FALSE` the parameter is treated as a constant.
    #' @param discr `logical(1)`\cr
    #' If `TRUE` the parameter is discrete.
    #' @param prior `character(1)`\cr
    #' Name of the prior distribution.
    #' @param args named `list()`\cr
    #' Arguments for the prior distribution.
    #' @param rule `character(1)`\cr
    #' Name of the inference proposal rule.
    #' @param scale `double(1)`\cr
    #' Scale parameter for the proposal rule.
    #' Changed automatically during inference if `n_adpat`>0.
    #' @param min `double(1)`\cr
    #' Value of the minimum.
    #' @param max `double(1)`\cr
    #' Value of the maximum.
    add_param = function(name = NA, value = NA,
                         const = NA, discr = FALSE,
                         prior = NA, args = list(),
                         rule = NA, scale = NA,
                         min = NA, max = NA) {
      private$.param$add(name = name, value = value, const = const,
                         discr = discr, prior = prior, args = args,
                         rule = rule, scale = scale,
                         min = min, max = max)
      invisible(self)
    },

    #' @description
    #' Load parameter definitions either as direct arguments
    #' or from predefined list.
    #'
    #' @param ...
    #' Named arguments with respective list of definition details.
    #' @param from (named `list()`\cr
    #' Named list of parameter definitions.
    #' @param append `logical(1)`\cr
    #' If `TRUE` appends to current definitions.
    #' If `FALSE` drop all pre-existing definitions (default).
    load_param_list = function(..., from = NA, append = FALSE) {
      private$.param$load_list(..., from = from, append = append)
      invisible(self)
    },

    #' @description
    #' Load parameter definitions from model template.
    #'
    #' @param append `logical(1)`\cr
    #' If `TRUE` append missing parameter definitions.
    #' If `FALSE` drop all pre-existing definitions (default).
    load_param_template = function(append = FALSE) {
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
    #' @param name `character(1)`\cr
    #' Name of the augmented variable.
    #' @param status `integer(1)`\cr
    #' The associated status of the definition.
    #' Default NA is for all infectious.
    #' @param value `double(1)` | `character(1)`\cr
    #' One between:
    #' - Scalar value;
    #' - Parameter name (optionally followed by `@P`);
    #' - Data name (followed by `@D`).
    #' @param init_rand `character(1)`\cr
    #' If `TRUE` draw initial vale from model distribution.
    #' @param const `logical(1)`\cr
    #' If `FALSE` the augmented variable is treated as a constant.
    #' @param discr `logical(1)`\cr
    #' IF `TRUE` the augmenter variable is discrete.
    #' @param func `character()`\cr
    #' Name of the augmented function.
    #' @param func_args (named `list()`\cr
    #' Arguments of the augmented function.
    #' @param model `character(1)`\cr
    #' Name of the model distribution.
    #' @param args named `list()`\cr
    #' Arguments for the model distribution.
    #' @param min `double(1)` | `character(1)`\cr
    #' Minimum value. One between:
    #' - Scalar value;
    #' - Parameter name (optionally followed by `@P`);
    #' - Data name (followed by `@D`).
    #' @param max `double(1)` | `character(1)`\cr
    #' Maxmimum tvalue. One between:
    #' - Scalar value;
    #' - Parameter name (optionally followed by `@P`);
    #' - Data name (followed by `@D`).
    #' @param rule `character(1)`\cr
    #' Name of the inference proposal rule.
    #' @param scale `double(1)`\cr
    #' Scale parameter for the proposal rule.
    #' @param from `double(1)` | `character(1)`\cr
    #' Advanced options for non random walk moves.
    #' `FIXME: missing`
    add_augm = function(name = NA, status = NA,
                        const = NA, discr = FALSE,
                        func = NA, func_args = NULL,
                        value = NA, init_rand = FALSE,
                        model = NA, args = list(),
                        min = NA, max = NA,
                        rule = NA, scale = NA, from = NA) {
      private$.augm$add(name = name, status = status,
                        const = const, discr = discr,
                        func = func, func_args = func_args,
                        value = value, init_rand = init_rand,
                        model = model, args = args,
                        min = min, max = max,
                        rule = rule, scale = scale, from = from)
      invisible(self)
    },

    #' @description
    #' Load augmented variable definitions either as direct arguments
    #' or from predefined list.
    #'
    #' @param ...
    #' Named arguments with respective list of definition details.
    #' @param from (named `list()`\cr
    #' Named list of parameter definitions.
    #' @param append `logical(1)`\cr
    #' If `TRUE` appends to current definitions.
    #' If `FALSE` drop all pre-existing definitions (default).
    load_augm_list = function(..., from = NA, append = FALSE) {
      private$.augm$load_list(..., from = from, append = append)
      invisible(self)
    },

    #' @description
    #' Load augmented variables from template definitions.
    #'
    #' @param append
    #' If `TRUE` append missing augmented variable definitions.
    #' If `FALSE` drop all pre-existing definitions (default).
    load_augm_template = function(append = FALSE) {
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
    #' @param name `character(1)`\cr
    #' Name of the group rate.
    #' @param by `list()`\cr
    #' Feature(s) associated to group rate as dataset headers.
    #' Corresponding dataset column have to be categorical.
    #' @param prefix `character(1)`\cr
    #' Prefix to parameter names associated to the group rate.
    #' @param reference named `list()`\cr
    #' Choose reference group providing corresponding features.
    #' Default is the group with the largest sample size.
    #' @param independent `logical(1)`\cr
    #' If `TRUE` then features are considered independently.
    #' @param prior `character(1)`\cr
    #' Name of the prior distribution.
    #' Default associated parameter definitions
    #' have prior "lnorm" and `args = list(meanlog = 0, sdlog = 3)`.
    #' @param args named `list()`\cr
    #' Arguments for the prior distribution.
    #' @param min `double(1)`\cr
    #' Value of the minimum.
    #' @param max `double(1)`\cr
    #' Value of the maximum.
    #' @param rule `character(1)`\cr
    #' Name of the inference proposal rule.
    #' @param scale `double(1)`\cr
    #' Scale parameter for the proposal rule.
    #' Changed automatically during inference if `n_adpat` > 0.
    #'
    #' @references
    #' Cauchemez, et al. (2009)
    #' Household Transmission of 2009 Pandemic Influenza A (H1N1) Virus in the United States.
    #' New England Journal of Medicine. 361: 2619-2627. doi:10.1056/NEJMoa0905498
    add_group_rate = function(name, by,
                              prefix = NA, reference = NULL, independent = TRUE,
                              prior = "lnorm",
                              args = list(meanlog = 0, sdlog = 3),
                              min = NA, max = 10,
                              rule = "lnorm", scale = 1) {
      private$.group_rate$add(name = name, by = by, prefix = prefix,
                              reference = reference, independent = independent,
                              prior = prior, args = args,
                              min = min, max = max,
                              rule = rule, scale = scale)
      invisible(self)
    },

    #' @description
    #' Load group rate definitions either as direct arguments
    #' or from predefined list.
    #'
    #' @param ...
    #' Named arguments with respective list of definition details.
    #' @param from (named `list()`\cr
    #' Named list of parameter definitions.
    #' @param append `logical(1)`\cr
    #' If `TRUE` appends to current definitions.
    #' If `FALSE` drop all pre-existing definitions (default).
    load_group_rate_list = function(..., from = NA, append = FALSE) {
      private$.group_rate$load_list(..., from = from, append = append)
      invisible(self)
    },

    #' @description
    #' Add definition of infectivity profile.
    #'
    #' @param name `character(1)`\cr
    #' Name of the infectivity profile.
    #' @param discr `logical(1)`\cr
    #' If `TRUE` the profile is discrete.
    #' @param args named `list()`\cr
    #' Arguments (eventually with types) for the profile.
    #' @param min `double(1)` | `character(1)`\cr
    #' One between:
    #' - Scalar value;
    #' - Parameter name (optionally followed by `@P`).
    #' @param max `double(1)` | `character(1)`\cr
    #' One between:
    #' - Scalar value;
    #' - Parameter name (optionally followed by `@P`).
    add_infec_profile = function(name, discr = FALSE,
                                 args = list(), min = NA, max = NA) {
      private$.infec_profile$add(name = name, discr = discr,
                                 args = args, min = min, max = max)
      invisible(self)
    },

    #' @description
    #' Showcase instantaneous and cumulative values of infectivity profile.
    #'
    #' @param bins `integer(1)`\cr
    #' Number of bins for continuous profile (ignored if discrete).
    #'
    #' @return
    #' `list()`
    showcase_infec_profile = function(bins = 2000) {
      # internal build function
      private$.build_minimal()
      return(private$.model$showcase_infec_profile(bins = bins))
    },

    ############################################################################
    ############################################################################
    ############################################################################

    #' @description
    #' Run Bayesian Markov Chain Monte Carlo (MCMC) inference
    #' and store internally the chain traces.
    #'
    #' @param dataset `data.frame` \cr
    #' Input dataset for analysis.
    #' @param sample_size `integer(1)` \cr
    #' Number of samples to draw from the posterior distribution in each chain.
    #' @param n_steps `integer(1)` \cr
    #' Total number of MCMC steps to perform.
    #' @param n_burnin `integer(1)` \cr
    #' Number of initial steps to discard as burn-in period.
    #' @param thin `integer(1)` \cr
    #' Thinning interval: number of steps to skip between each saved sample.
    #' @param n_adapt `integer(1)` \cr
    #' Number of adaptation steps before starting the main sampling.
    #' @param n_chains `integer(1)` \cr
    #' Number of independent MCMC chains to run.
    #' @param n_threads `integer(1)` \cr
    #' Number of parallel threads to use for computations.
    #' @param rand_init `logical(1)` \cr
    #' Controls random initialization of parameters:
    #' - If `TRUE` (default), initial values are ndomly drawn from the prior distribution;
    #' - If `FALSE` initial value definitions are used.
    #' @param id_trace `character(1)` \cr
    #' Individual identifier for augmented trace plots. Default is `NA`.
    #' @param inv_temper `double(1)` \cr
    #' Inverse temperature parameter.
    #' Typically set to 1 for standard MCMC sampling.
    #' @param seed `integer(1)` \cr
    #' Seed for the random number generator.
    #' - If a specific integer is provided, sets the random number generator seed.
    #' - If `NA` (default), uses the system clock to generate a seed.
    #' @param progress `logical(1)` | `numeric(1)`\cr
    #' Controls progress visualization:
    #' - If `TRUE` (default), prints progress every second;
    #' - If an integer `n`, prints progress every `n` seconds;
    #' - If `FALSE` or `0`, suppresses progress display.
    infer = function(dataset = NULL,
                     sample_size = NA, n_steps = 10000,
                     n_burnin = 0, thin = 1, n_adapt = 0,
                     n_chains = 1, n_threads = 1,
                     rand_init = TRUE, id_trace = NA, inv_temper = 1.0,
                     seed = NA, progress = TRUE) {
      # internal build function
      private$.build_full(dataset = dataset, sample_size = sample_size,
                          n_steps = n_steps, n_burnin = n_burnin, thin = thin,
                          n_adapt = n_adapt, n_chains = n_chains,
                          n_threads = n_threads, rand_init = rand_init,
                          id_trace = id_trace, inv_temper = inv_temper,
                          seed = seed, progress = progress)
      private$.multiplicity_engine(dataset = dataset)

      # infer method from cpp model object
      private$.model$infer()
      invisible(self)
    },

    #' @description
    #' Run simulation and return data.frame with simulated dataset.
    #'
    #' @param dataset `data.frame` \cr
    #' Input dataset for simulation.
    #' @param seed `integer(1)` \cr
    #' Seed for the random number generator.
    #' - If a specific integer is provided, sets the random number generator seed.
    #' - If `NA` (default), uses the system clock to generate a seed.
    #' @param progress `logical(1)` | `numeric(1)`\cr
    #' Controls progress visualization:
    #' - If `TRUE` (default), prints progress every second;
    #' - If an integer `n`, prints progress every `n` seconds;
    #' - If `FALSE` or `0`, suppresses progress display.
    #'
    #' @return
    #' `data.frame`
    simulate = function(dataset = NULL, seed = NA, progress = FALSE) {
      # internal build function
      private$.build_full(dataset = dataset, rand_init = FALSE,
                          seed = seed, progress = progress)

      # call method from internal Rcpp object
      return(private$.model$simulate())
    },

    #' @description
    #' Extract the posterior expectation of a given measure.
    #' Perform simulations over posterior sample generated by `$infer`.
    #'
    #' @param FUN `closure()` | `list()` of `closure()`\cr
    #' Function(s) to extract the desired measure from a dataset.
    #' @param progress `logical(1)` | `numeric(1)`\cr
    #' Controls progress visualization. Default is inherited from the model.
    #' @param ...
    #' Named arguments with respective functions.
    #' If `FUN` is provided, additional arguments are ignored.
    #'
    #' @return
    #' `HousepiPosteriorExpec`\cr
    #' An instance of class
    #' [HousepiPosteriorExpec][housepi::HousepiPosteriorExpec].
    expec_posterior_measure = function(FUN = NULL, progress = NA, ...) {
      return(
        HousepiPosteriorExpec$new(self, FUN = FUN, progress = progress, ...)
      )
    },

    #' @description
    #' Extract the expected number of household cases per size.
    #' It can be used to estimate the secondary attack rate.
    #' Perform simulations over posterior sample generated by `$infer`.
    #'
    #' @param progress `logical(1)` | `numeric(1)`\cr
    #' Controls progress visualization. Default is inherited from the model.
    #'
    #' @return
    #' `ExpectedCasesBySize`\cr
    #' An instance of class
    #' [ExpectedCasesBySize][housepi::ExpectedCasesBySize]
    #' containing method to compare observed and expected cases.
    expec_cases_by_house_size = function(progress = NA) {
      return(
        ExpectedCasesBySize$new(self, progress = progress)
      )
    },

    #' @description
    #' Evaluate the loglikelihood for the given parameters.
    #'
    #' @param dataset `data.frame` \cr
    #' Input dataset for simulation.
    #' @param seed `integer(1)` \cr
    #' Seed for the random number generator.
    #' - If a specific integer is provided, sets the random number generator seed.
    #' - If `NA` (default), uses the system clock to generate a seed.
    #'
    #' @return
    #' `double(1)`
    loglikelihood = function(dataset = NULL, seed = NA) {
      # internal build function
      private$.build_full(dataset = dataset, seed = seed, rand_init = FALSE)
      private$.multiplicity_engine(dataset = dataset)
      return(private$.model$evaluate())
    },

    #' @description
    #' Estimate the model evidence (or marginal likelihood).
    #' It uses the method of power posteriors with different available schemes.
    #'
    #' @param dataset `data.frame` \cr
    #' Input dataset for simulation.
    #' @param scheme `character(1)` \cr
    #' Scheme used for evidence estimation.
    #' - If "power-fractions" (default) uses inverse temperatures (r / n_rungs)^5.
    #' - If "adaptive" uses adaptive thermodynamic ladder.
    #' @param n_rungs `integer(1)` \cr
    #' Number of rungs in the thermodynamic ladder.
    #' @param sample_size `integer(1)` \cr
    #' Number of samples to draw from the posterior distribution in each chain.
    #' @param n_steps `integer(1)` \cr
    #' Total number of MCMC steps to perform.
    #' @param n_burnin `integer(1)` \cr
    #' Number of initial steps to discard as burn-in period.
    #' @param thin `integer(1)` \cr
    #' Thinning interval: number of steps to skip between each saved sample.
    #' @param n_adapt `integer(1)` \cr
    #' Number of adaptation steps before starting the main sampling.
    #' @param n_chains `integer(1)` \cr
    #' Number of independent MCMC chains to run.
    #' @param n_threads `integer(1)` \cr
    #' Number of parallel threads to use for computations.
    #' @param seed `integer(1)` \cr
    #' Seed for the random number generator.
    #' - If a specific integer is provided, sets the random number generator seed.
    #' - If `NA` (default), uses the system clock to generate a seed.
    #' @param progress `logical(1)` | `numeric(1)`\cr
    #' Controls progress visualization:
    #' - If `TRUE` (default), prints progress every second;
    #' - If an integer `n`, prints progress every `n` seconds;
    #' - If `FALSE` or `0`, suppresses progress display.
    #'
    #' @return
    #' `double(1)`
    #'
    #' @references
    #' Friel N, Hurn M, Wyse J.
    #' Improving power posterior estimation of statistical evidence.
    #' Stat Comput. 2014;24: 709-723. doi:10.1007/s11222-013-9397-1
    evidence = function(dataset = NULL,
                        scheme = "power-fractions", n_rungs = 40,
                        sample_size = NA, n_steps = 10000,
                        n_burnin = 0, thin = 1, n_adapt = 0,
                        n_chains = 1, n_threads = 1,
                        seed = NA, progress = TRUE) {
      # internal build function
      private$.build_full(dataset, sample_size = sample_size, n_steps = n_steps,
                          n_burnin = n_burnin, thin = thin, n_adapt = n_adapt,
                          n_chains = n_chains, n_threads = n_threads,
                          rand_init = TRUE, seed = seed, progress = progress)
      private$.multiplicity_engine(dataset = dataset)

      # evaluate marginal from the internal Rcpp object
      return(private$.model$evidence(scheme = scheme, n_rungs = n_rungs))
    }
  )
)
