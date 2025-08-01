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
#' Compile household model, manage your dataset and parameters from a single class.
#'
#' @keywords internal
#' @noRd
HousepiModel <- R6::R6Class(
  classname = "HousepiModel",
  inherit = HousepiGenerator,
  cloneable = TRUE,
  lock_class = TRUE,

  ############################################################################
  ################################## PRIVATE #################################
  ############################################################################

  private = list(
    deep_clone = function(name, value) helper_clone(value),
    .param = NULL,
    .augm = NULL,
    .group_rates = NULL,
    .infec_profile = NULL,

    # @description
    # Dataset headers required by implementation file and other definitions.
    #
    # @return `character()`
    .get_all_required_dataset = function() {
      required <- c()
      if (!private$.is_blank()) required <- private$.dataset_required()
      return(c(
        required, # by the model
        private$.augm$dataset_required, # by augmented variable
        private$.infec_profile$dataset_required # by infectivity profile
      ))
    },

    # @description 
    # Assert if dataset headers match implementation file.
    # @param dataset_headers `character()`\cr
    # Vector with the dataset headers.
    .assert_dataset = function(dataset_headers) {
      if (length(dataset_headers) == 0) stop_h(
        "Impossible to proceed. No database was provided."
      )
      required <- private$.get_all_required_dataset()
      template_missing_stop(required, dataset_headers, "dataset")
      invisible(self)
    },

    # @description
    # Parameters required by implementation file and other definitions.
    #
    # @return `character()`
    .get_all_param = function() {
      if (private$.is_blank()) required <- self$parameters$name
      else required <- private$.param$require(private$.reference())
      return(c(
        required, # by the model
        private$.group_rates$param_required, # by group rates
        private$.augm$param_required, # by augmented variable
        private$.infec_profile$param_required # by infectivity profile
      ))
    },

    # @description
    # Assert if parameter definitions match required
    .assert_param = function() {
      private$.param$assert(private$.get_all_param())
      invisible(self)
    },

    # @description
    # Augmented variables required by implementation file and other definitions.
    #
    # @return `character()`
    .get_all_augm = function() {
      return(c(
        private$.augm$require(private$.reference()), # by the model
        private$.infec_profile$augm_required # by infectivity profile
      ))
    },

    # @description
    # Assert if augmented definitions match the required.
    .assert_augm = function() {
      if (!private$.is_blank()) private$.augm$assert(private$.get_all_augm())
      invisible(self)
    },

    # @description
    # Group rates required by implementation file and other definitions.
    #
    # @return `character()`
    .get_all_group_rate = function() {
      return(private$.group_rates$require(private$.reference()))
    },

    # @description
    # Assert group rate definitions from the model implementation script.
    .assert_group_rate = function() {
      if (!private$.is_blank()) {
        private$.group_rates$assert(private$.get_all_group_rate())
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
      refer_file <- private$.reference()
      regular_expr <- ".multiplicity\\s*\\(\\s*\\)"

      if (template_assert_expression(refer_file, regular_expr)) {

        # note: this might be a modeler error
        if (length(private$.get_all_augm()) > 0) warning_h(
          "Model with augmentated variable called multiplicity engine."
        )

        # if dataset not provided, load previous
        # UPDATE: this "back & forth" of dataset is not optimal!
        if ((length(dataset) == 0) && (length(private$.cols_dataset) > 0)) {
          dataset <- self$dataset
        }

        headers <- c(
          # dataset headers required by implementation
          private$.get_all_required_dataset(),
          # dataset headers required by group rates
          unlist(lapply(
            private$.group_rates$logic_tables,
            function(x) unlist(lapply(x, names))
          ))
        )
        # remove logic tables auxiliary column names
        headers <- setdiff(unique(headers), c("n", "param_name"))

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
        defs <- data.frame(house_id = unique(dataset$house_id),
                           multiplicity = 0,
                           class_id = NA,
                           stringsAsFactors = FALSE)
        rownames(defs) <- defs$house_id
        defs[multi_hh_table$class_id, "multiplicity"] <- multi_hh_table$multiplicity

        for (i in seq_len(nrow(multi_hh_table))) {
          defs[hash_groups[[i]], "class_id"] <- multi_hh_table$class_id[i]
        }
        # call rcpp method
        private$.rcpp$import_multiplicity(defs)
      }

      invisible(self)
    },

    ############################################################################
    ############################################################################
    ############################################################################

    # @description
    # Pre-build to deal with group rates.
    #
    # @param dataset `data.frame` \cr
    # Input dataset for simulation.
    .build_group_rates = function() {
      # assert definitions are complete
      private$.assert_group_rate()

      # proceed only if definitions exist
      if (length(private$.group_rates$defs) > 0) {
        # if dataset not provided, load previous
        # UPDATE: count in the wrapper while storing dataset_levels
        if (length(private$.cols_dataset) == 0) dataset <- NULL
        else dataset <- self$dataset

        # implement group rates
        private$.group_rates$implement(dataset, private$.info_dataset)

        # append parameter definitions if missing
        new_parameters <- private$.group_rates$parameters
        # find missing parameter definitions
        mask <- setdiff(names(new_parameters), self$parameters$name)
        new_parameters <- new_parameters[mask]
        if (length(new_parameters) > 0) {
          self$load_parameter_list(from = new_parameters, append = TRUE)
        }
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
    # @param id_draws `character()`\cr
    # Individual identifiers to show augmented traces. 
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
                           rand_init = TRUE, id_draws = NA,
                           seed = NA, progress = FALSE) {
      # build dataset
      self$add_dataset(dataset)

      # individual id for blank model augmented traces
      # if no default columns, assume debugging mode
      if (private$.is_blank() && all(is.na(id_draws))) {
        # FIXME: not robust in case of dataset with only default columns
        if (length(setdiff(private$.cols_dataset,
                           c("house_id", "id", "status"))) == 0) {
          id_draws <- template_blank_dataset_()$id[1]
        }
      }

      # build setup
      private$.rcpp$import_setup(char_setup_input_(list(
        sample_size = ifelse(is.na(sample_size), 0, sample_size),
        n_steps = n_steps, thin = thin,
        n_burnin = n_burnin, n_adapt = n_adapt, n_chains = n_chains,
        n_threads = n_threads, rand_init = rand_init,
        id_draws = id_draws, 
        seed = private$.adjust_seed(seed), 
        progress = as.integer(progress)
      )))

      # assert parameter definitions before importing
      private$.assert_param()
      private$.rcpp$import_param(self$parameters)

      # group rates (import depends on param)
      private$.build_group_rates()
      defs <- private$.group_rates$tidy_table
      defs$n <- NULL
      info <- private$.info_dataset
      for (k in intersect(colnames(defs), names(info))) {
        defs[[k]] <- as.integer(factor(defs[[k]], levels = info[[k]]$levels))
      }
      private$.rcpp$import_group_rates(defs)

      # build contact rates (FIXME: temporary)
      private$.rcpp$build_all_rates()

      # assert augmented variable definitions
      private$.assert_augm()

      # build infectivity profile
      private$.rcpp$import_infec_profile(private$.infec_profile$defs)

      # build augmented variables
      private$.rcpp$import_augm(self$augmented)

      invisible(self)
    },

    # UPDATE:
    # this could be designed to be more helpful
    # e. g. create a standalone folder somewhere with correct links to headers
    # FIXME:
    # setup is broken, follow method in copy to fix it
    .export_to_cpp_standalone = function(outpath = "../housepi/inst/standalone/") {
      write.csv(private$.rcpp$export_dataset(),
                file.path(outpath, "dataset.tmp"),
                quote = FALSE, row.names = FALSE)
      write.csv(private$.rcpp$export_datalist(),
                file.path(outpath, "datalist.tmp"),
                quote = FALSE, row.names = FALSE)

      defs <- self$group_rates
      info <- private$.info_dataset
      for (k in intersect(colnames(defs), names(info))) {
        defs[[k]] <- as.integer(factor(defs[[k]], 
                                levels = info[[k]]$levels))
      }
      write.csv(defs,
                file.path(outpath, "group_rates.tmp"),
                quote = FALSE, row.names = FALSE)

      write.csv(self$augmented,
                file.path(outpath, "augm.tmp"),
                quote = FALSE, row.names = FALSE)
      write.csv(self$parameters,
                file.path(outpath, "param.tmp"),
                quote = FALSE, row.names = FALSE)
      write.table(t(as.data.frame(char_setup_input_(self$setup))),
                  file.path(outpath, "setup.tmp"),
                  row.names = TRUE, col.names = FALSE, sep = ",", quote = FALSE)
      write.table(t(as.data.frame(private$.infec_profile$defs)),
                  file.path(outpath, "infec.tmp"),
                  row.names = TRUE, col.names = FALSE, sep = ",", quote = FALSE)
      invisible(self)
    }
  ),

  ############################################################################
  ################################## ACTIVE ##################################
  ############################################################################

  active = list(
    #' @field dataset `data.frame`\cr
    #' Dataset with individual level information of current dataset.
    dataset = function(value) {
      if (missing(value)) {
        if (is.null(private$.rcpp)) stop_h(
          "Impossible to export dataset before model build."
        )
        return(private$.refactor_dataset(private$.rcpp$export_dataset()))
      }
      else self$add_dataset(value)
    },

    #' @field datalist `data.frame`\cr
    #' Datalist with multiple rows indexed by individual.
    datalist = function(value) {
      if (missing(value)) {
        if (is.null(private$.rcpp)) stop_h(
          "Impossible to export datalist before model build."
        )
        return(private$.refactor_datalist(private$.rcpp$export_datalist()))
      }
      else self$add_datalist(value)
    },

    #' @field parameters `data.frame`\cr
    #' Parameter formatted definitions.
    parameters = function(value) {
      if (missing(value)) return(private$.param$defs)
      else self$load_parameter_list(from = value, append = FALSE)
    },

    #' @field augmented `data.frame`\cr
    #' Augmented variable formatted definitions.
    augmented = function(value) {
      if (missing(value)) return(private$.augm$defs)
      else self$load_augmented_list(from = value, append = FALSE)
    },

    #' @field infec_profile_defs `list`\cr
    #' Infectivity profile definition list.
    infec_profile_defs = function(value) {
      if (missing(value)) {
        defs <- private$.infec_profile$defs
        if (!is.null(defs)) return(defs)
        else return(private$.rcpp$export_infec_profile())
      }
      else do.call(self$add_infec_profile, value)
    },

    #' @field group_rates `list`\cr
    #' List containing current group rate definitions.
    group_rates = function(value) {
      if (missing(value)) {
        private$.build_group_rates()
        return(private$.group_rates$tidy_table)
      }
      else self$load_group_rate_list(from = value, append = FALSE)
    }
  ),

  ############################################################################
  ################################## PUBLIC ##################################
  ############################################################################

  public = list(
    #' @inherit new_housepi_model
    initialize = function(choice = NA, as_template = FALSE, force = FALSE,
                          name = NA, cri_level = 0.95, quiet = TRUE) {
      # initialize internal fields
      self$drop_all()

      # generate model wrapper
      super$initialize(choice = choice, as_template = as_template, 
                       force = force, name = name, quiet = quiet)

      # import blank template dataset
      if (private$.is_blank()) self$add_dataset(template_blank_dataset_())

      # credible interval level
      self$cri_level <- cri_level

      invisible(self)
    },

    ############################################################################
    ############################################################################
    ############################################################################

    #' @inherit drop_parameters
    drop_parameters = function() {
      private$.param <- HousepiParameter$new()
      invisible(self)
    },

    #' @inherit drop_augmented
    drop_augmented = function() {
      private$.augm <- HousepiAugmented$new()
      invisible(self)
    },

    #' @inherit drop_group_rates
    drop_group_rates = function() {
      # remove associated parameters
      if (!is.null(private$.group_rates)) {
        for (p in names(private$.group_rates$parameters)) {
          private$.param$remove(p, show_warnings = FALSE)
        }
      }
      private$.group_rates <- HousepiGroupRate$new()
      invisible(self)
    },

    #' @inherit drop_all
    drop_all = function() {
      self$drop_parameters()
      self$drop_augmented()
      self$drop_group_rates()
      private$.infec_profile <- HousepiInfecProfile$new()
      invisible(self)
    },

    ############################################################################
    ############################################################################
    ############################################################################

    #' @inherit add_dataset
    add_dataset = function(dataset) {
      if (length(dataset) > 0) {
        # check that dataset type is correct
        if (!is.data.frame(dataset)) stop_h(
          "Unsupported type for dataset.",
          "Dataset should be a data.frame."
        )

        # assert dataset headers
        private$.assert_dataset(colnames(dataset))
        # convert factors
        embed <- cpp_data_factorize(dataset)
        # store levels
        private$.info_dataset <- embed$info
        # store column names
        private$.cols_dataset <- embed$cols
        # import data.frame in rcpp
        private$.rcpp$import_dataset(embed$data)

      } else {
        # assert headers of last imported dataset
        private$.assert_dataset(private$.cols_dataset)
      }

      # implement group rates and build
      private$.build_group_rates()

      invisible(self)
    },

    #' @inherit add_datalist
    add_datalist = function(datalist) {
      if (!is.data.frame(datalist)) stop_h(
        "Unsupported type for datalist.",
        "Datalist should be a data.frame."
      )

      # convert factors
      embed <- cpp_data_factorize(datalist)
      # store levels
      private$.info_datalist <- embed$info
      # store column names
      private$.cols_datalist <- embed$cols
      # import data.frame in rcpp
      private$.rcpp$import_datalist(embed$data)
    
      invisible(self)
    },

    ############################################################################
    ############################################################################
    ############################################################################

    #' @inherit add_parameter
    add_parameter = function(name = NA, value = NA,
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

    #' @inherit load_parameter_list
    load_parameter_list = function(..., from = NA, append = FALSE) {
      private$.param$load_list(..., from = from, append = append)
      invisible(self)
    },

    #' @inherit load_parameter_template
    load_parameter_template = function(append = FALSE) {
      if (private$.is_blank()) warning_h(
        "No template parameters available for blank model."
      )
      # retrieve the template parameters file
      else private$.param$load_template(
        file.path(private$.path(), "template_param.R"),
        private$.param$require(private$.reference()),
        isTRUE(append)
      )
      invisible(self)
    },

    ############################################################################
    ############################################################################
    ############################################################################

    #' @inherit add_augmented
    add_augmented = function(name = NA, status = NA,
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

    #' @inherit load_augmented_list
    load_augmented_list = function(..., from = NA, append = FALSE) {
      private$.augm$load_list(..., from = from, append = append)
      invisible(self)
    },

    #' @inherit load_augmented_template
    load_augmented_template = function(append = FALSE) {
      if (private$.is_blank()) warning_h(
        "No template augmented variables available for blank model."
      )
      # retrieve the template augmented file
      else private$.augm$load_template(
        file.path(private$.path(), "template_augm.R"),
        private$.augm$require(private$.reference()),
        isTRUE(append)
      )
      invisible(self)
    },

    ############################################################################
    ############################################################################
    ############################################################################

    #' @inherit add_group_rate
    add_group_rate = function(name, by,
                              prefix = NA, 
                              reference = NULL, independent = TRUE,
                              prior = "lnorm", args = list(meanlog = 0, sdlog = 3),
                              min = NA, max = 10,
                              rule = "lnorm", scale = 1) {
      private$.group_rates$add(name = name, by = by, prefix = prefix,
                               reference = reference, independent = independent,
                               prior = prior, args = args,
                               min = min, max = max,
                               rule = rule, scale = scale)
      invisible(self)
    },

    #' @inherit load_group_rate_list
    load_group_rate_list = function(..., from = NA, append = FALSE) {
      if (isFALSE(append)) self$drop_group_rates()
      private$.group_rates$load_list(..., from = from, append = append)
      invisible(self)
    },

    #' @inherit add_infec_profile
    add_infec_profile = function(name, discr = FALSE,
                                 args = list(), max = NA) {
      private$.infec_profile$add(name = name, discr = discr,
                                 args = args, max = max)
      invisible(self)
    },

    ############################################################################
    ############################################################################
    ############################################################################

    #' @inherit infer
    infer = function(dataset = NULL,
                     sample_size = NA, n_steps = 10000,
                     n_burnin = 0, thin = 1, n_adapt = 0,
                     n_chains = 1, n_threads = 1,
                     rand_init = TRUE, id_draws = NA,
                     seed = NA, progress = TRUE) {
      # internal build function
      private$.build_full(dataset = dataset, sample_size = sample_size,
                          n_steps = n_steps, n_burnin = n_burnin, thin = thin,
                          n_adapt = n_adapt, n_chains = n_chains,
                          n_threads = n_threads, rand_init = rand_init,
                          id_draws = id_draws, seed = seed, progress = progress)
      private$.multiplicity_engine(dataset = dataset)

      # infer method from cpp model object
      private$.rcpp$infer()
      invisible(self)
    },

    #' @inherit simulate
    simulate = function(dataset = NULL, seed = NA, progress = FALSE) {
      # internal build function
      private$.build_full(dataset = dataset, rand_init = FALSE,
                          seed = seed, progress = progress)

      # refactor static dataset of simulation
      static_df <- private$.rcpp$dataset_for_simulation()

      # bind simulated entries with static part
      # UPDATE : n_rep = 1 for development
      return(private$.refactor_dataset(
        cbind(static_df, private$.rcpp$simulate(1)))
      )
    },

    #' @description
    #' Extract the posterior expectation of a given measure.
    #' Perform simulations over posterior sample generated by `$infer`.
    #'
    #' @param FUN `closure` | `list` of `closure`\cr
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
      HousepiPosteriorExpec$new(self, FUN = FUN, progress = progress, ...)
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
      ExpectedCasesBySize$new(self, progress = progress)
    },

    #' @inherit loglikelihood
    loglikelihood = function(dataset = NULL, seed = NA) {
      # internal build function
      private$.build_full(dataset = dataset, seed = seed, rand_init = FALSE)
      private$.multiplicity_engine(dataset = dataset)
      return(private$.rcpp$evaluate())
    },

    #' @inherit evidence
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

      # load options
      scheme <- as.character(scheme)
      n_rungs <- as.integer(n_rungs)

      # evaluate marginal from the internal Rcpp object
      output <- private$.rcpp$estimate_evidence(scheme, n_rungs)

      # reformat dataframe
      mask <- c("inv_temper", "mean_loglike", "var_loglike")
      output$detail <- as.data.frame(output[mask])
      output[mask] <- NULL

      return(output)
    }
  )
)
