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
#' @export HouseholdModel
#' @docType class
HouseholdModel <- R6::R6Class(
  classname = "HouseholdModel",
  inherit = HousepiGenerator,
  cloneable = TRUE,
  lock_class = TRUE,

  ############################################################################
  ################################## PRIVATE #################################
  ############################################################################

  private = list(
    .param = NULL,
    .augm = NULL,
    .group_rates = NULL,
    .infec_profile = NULL,

    # @description
    # Dataset headers required by implementation file and other definitions.
    #
    # @return `character()`
    .require_dataset = function() {
      requirements <- c()
      if (!private$.is_blank()) {
        requirements <- private$.dataset_required()
      }
      return(c(
        requirements, # by the model
        private$.augm$dataset_required # by augmented variable
      ))
    },

    # @description
    # Parameters required by implementation file and other definitions.
    #
    # @return `character()`
    .require_param = function() {
      if (private$.is_blank()) {
        requirements <- self$param_defs$name
      } else {
        requirements <- private$.param$require(private$.reference())
      }
      return(c(
        requirements, # by the model
        private$.group_rates$param_required, # by group rates
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
        private$.augm$require(private$.reference()), # by the model
        private$.infec_profile$augm_required # by infectivity profile
      ))
    },

    # @description
    # Assert if augmented definitions match the requirements.
    .assert_augm = function() {
      if (!private$.is_blank()) {
        private$.augm$assert(private$.require_augm())
      }
      invisible(self)
    },

    # @description
    # Group rates required by implementation file and other definitions.
    #
    # @return `character()`
    .require_group_rate = function() {
      return(private$.group_rates$require(private$.reference()))
    },

    # @description
    # Assert group rate definitions from the model implementation script.
    .assert_group_rate = function() {
      if (!private$.is_blank()) {
        private$.group_rates$assert(private$.require_group_rate())
      }
      invisible(self)
    },

    ############################################################################
    ############################################################################
    ############################################################################

    # @description
    # Import a new dataset and drop previous.
    #
    # @param dataset `data.frame`\cr
    # A data.frame containing the dataset with individual information by row.
    .import_dataset = function(dataset = NULL) {
      if (length(dataset) > 0) {
        # check that dataset type is correct
        if (!is.data.frame(dataset)) stop_h(
          "Unsupported type for dataset.",
          "Dataset should be a data.frame."
        )

        # assert dataset headers
        private$.assert_dataset(colnames(dataset))

        # convert factors to numeric and store levels
        private$.dataset_levels <- list()
        for (col in setdiff(colnames(dataset), c("id", "house_id", "status"))) {
          x <- dataset[[col]]
          if (is.character(x) || is.factor(x) || is.logical(x)) {
            x <- factor(x)
            private$.dataset_levels[[col]] <- levels(x)
            dataset[[col]] <- as.integer(x)
          }
        }

        # import data.frame in rcpp # note: it calls
        private$.rcpp$import_dataset(fmt_data_input_(dataset))
        # store column names
        private$.last_dataset_cols <- colnames(dataset)

      } else {
        # assert headers of last imported dataset
        private$.assert_dataset(private$.last_dataset_cols)
      }

      # implement group rates and build
      private$.build_group_rates()

      invisible(self)
    },

    # @description
    # Import a new datalist and drop previous.
    #
    # @param datalist `data.frame`\cr
    # A data.frame containing the dataset with individual information by row.
    .import_datalist = function(datalist = NULL) {
      if (!is.data.frame(datalist)) stop_h(
        "Unsupported type for datalist.",
        "Datalist should be a data.frame."
      )

      # convert factors to numeric and store levels
      private$.datalist_levels <- list()
      for (col in setdiff(colnames(datalist), c("id", "house_id", "status"))) {
        x <- datalist[[col]]
        if (is.character(x) || is.factor(x) || is.logical(x)) {
          x <- factor(x)
          private$.datalist_levels[[col]] <- levels(x)
          datalist[[col]] <- as.integer(x)
        }
      }

      # import data.frame in rcpp
      private$.rcpp$import_datalist(fmt_data_input_(datalist))
      # store column names
      private$.last_datalist_cols <- colnames(datalist)
    
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
        if (length(private$.require_augm()) > 0) warning_h(
          "Model with augmentated variable called multiplicity engine."
        )

        # if dataset not provided, load previous
        # UPDATE: this "back & forth" of dataset is not optimal!
        if ((length(dataset) == 0) && (length(private$.last_dataset_cols) > 0)) {
          dataset <- self$dataset
        }

        headers <- c(
          # dataset headers required by implementation
          private$.require_dataset(),
          # dataset headers required by group rates
          unlist(lapply(
            private$.group_rates$logic_tables,
            function(x) unlist(lapply(x, names))
          ))
        )
        # remove logic tables auxiliary column names
        headers <- setdiff(unique(headers), c("counts", "param_name"))

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
      refer_file <- private$.reference()
      requirements <- private$.group_rates$require(refer_file)
      private$.group_rates$assert(requirements)

      # proceed only if definitions exist
      if (length(private$.group_rates$defs) > 0) {
        # if dataset not provided, load previous
        # UPDATE: count in the wrapper while storing dataset_levels
        if (length(private$.last_dataset_cols) == 0) {
          dataset <- NULL
        } else {
          dataset <- self$dataset
        }

        # implement group rates
        private$.group_rates$implement(dataset, private$.dataset_levels)

        # append parameter definitions if missing
        new_param_defs <- private$.group_rates$param_defs
        # find missing parameter definitions
        mask <- setdiff(names(new_param_defs), self$param_defs$name)
        new_param_defs <- new_param_defs[mask]
        if (length(new_param_defs) > 0) {
          self$load_param_list(from = new_param_defs, append = TRUE)
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
    # @param id_trace `character(1)`\cr
    # Individual identifier to show augmented traces. Default is NA.
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
                           rand_init = TRUE, id_trace = NA,
                           seed = NA, progress = FALSE) {
      # import dataset
      private$.import_dataset(dataset)

      # assert parameter definitions before importing
      private$.assert_param()
      private$.rcpp$import_param(self$param_defs)

      # group rates (import depends on param)
      defs <- self$group_rate_defs
      lvls <- private$.dataset_levels
      for (k in intersect(colnames(defs), names(lvls))) {
        defs[[k]] <- as.integer(factor(defs[[k]], levels = lvls[[k]]))
      }
      private$.rcpp$import_group_rates(defs)

      # build contact rates (UPDATE: in development)
      private$.rcpp$build_all_rates()

      # assert augmented variable definitions
      private$.assert_augm()

      # infectivity profile
      private$.rcpp$import_infec_profile(private$.infec_profile$defs)

      # individual id for blank model augmented traces
      # if no default columns, assume debugging mode
      if (private$.is_blank() && is.na(id_trace)) {
        # UPDATE: not robust in case of dataset with only default columns
        if (length(setdiff(private$.last_dataset_cols,
                           c("house_id", "id", "status"))) == 0) {
          id_trace <- template_blank_dataset_()$id[1]
        }
      }

      # control that seed is 32bit
      if (isTRUE(abs(seed) > .Machine$integer.max)) stop_h(
        "Supplied seed is not supported.",
        sprintf("Maximum value for a integer in R is %.0f.",
                .Machine$integer.max)
      )

      # import setup
      private$.rcpp$import_setup(fmt_setup_input_(list(
        sample_size = ifelse(is.na(sample_size), 0, sample_size),
        n_steps = n_steps, thin = thin,
        n_burnin = n_burnin, n_adapt = n_adapt, n_chains = n_chains,
        n_threads = n_threads, rand_init = rand_init,
        id_trace = id_trace, seed = seed, progress = as.integer(progress)
      )))

      # build augmented variables
      if (length(self$augm_defs) > 0) {
        private$.rcpp$import_augm(self$augm_defs)
      }

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

      defs <- self$group_rate_defs
      lvls <- private$.dataset_levels
      for (k in intersect(colnames(defs), names(lvls))) {
        defs[[k]] <- as.integer(factor(defs[[k]], levels = lvls[[k]]))
      }
      write.csv(defs,
                file.path(outpath, "group_rates.tmp"),
                quote = FALSE, row.names = FALSE)

      write.csv(self$augm_defs,
                file.path(outpath, "augm.tmp"),
                quote = FALSE, row.names = FALSE)
      write.csv(self$param_defs,
                file.path(outpath, "param.tmp"),
                quote = FALSE, row.names = FALSE)
      write.table(t(as.data.frame(fmt_setup_input_(self$setup))),
                  file.path(outpath, "setup.tmp"),
                  row.names = TRUE, col.names = FALSE, sep = ",", quote = FALSE)
      write.table(t(as.data.frame(self$infec_profile_defs)),
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
      else private$.import_dataset(value)
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
      else private$.import_datalist(value)
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
      # initialize internal fields
      self$drop_all()

      # generate model wrapper
      super$initialize(choice = choice, as_template = as_template, 
                       force = force, name = name, quiet = quiet)

      # import blank template dataset
      if (private$.is_blank()) private$.import_dataset(template_blank_dataset_())

      # credible interval level
      self$cri_level <- cri_level

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
      # remove associated parameters
      if (!is.null(private$.group_rates)) {
        for (p in names(private$.group_rates$param_defs)) {
          private$.param$remove(p, show_warnings = FALSE)
        }
      }
      private$.group_rates <- HousepiGroupRate$new()
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

    #' @description
    #' Add definitions of a group rate by given dataset feature(s).
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
      private$.group_rates$add(name = name, by = by, prefix = prefix,
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
      private$.group_rates$load_list(..., from = from, append = append)
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
    #' @param max `double(1)` | `character(1)`\cr
    #' One between:
    #' - Scalar value;
    #' - Parameter name (optionally followed by `@P`).
    add_infec_profile = function(name, discr = FALSE,
                                 args = list(), max = NA) {
      private$.infec_profile$add(name = name, discr = discr,
                                 args = args, max = max)
      invisible(self)
    },

    #' @description
    #' Showcase instantaneous and cumulative values of infectivity profile.
    #' UPDATE: this should be dropped and just parameters should be exported
    #' to then be passed to housepi_infec_profile
    #'
    #' @param bins `integer(1)`\cr
    #' Number of bins for continuous profile (ignored if discrete).
    #'
    #' @return
    #' `list()`
    showcase_infec_profile = function(bins = 2000) {
      # load options
      bins <- as.integer(bins)

      # assert parameter definitions before importing
      requirements <- private$.infec_profile$param_required
      if (length(requirements) > 0) {
        private$.param$assert(requirements)
        private$.rcpp$import_param(self$param_defs)
      }

      # infectivity profile
      private$.rcpp$import_infec_profile(private$.infec_profile$defs)

      # call Rcpp method
      output <- private$.rcpp$showcase_infec_profile(bins)
      # reformat dataframe
      mask <- c("gentime", "inst", "cumul")
      output$profile <- as.data.frame(output[mask])
      output[mask] <- NULL

      return(output)
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
                     rand_init = TRUE, id_trace = NA,
                     seed = NA, progress = TRUE) {
      # internal build function
      private$.build_full(dataset = dataset, sample_size = sample_size,
                          n_steps = n_steps, n_burnin = n_burnin, thin = thin,
                          n_adapt = n_adapt, n_chains = n_chains,
                          n_threads = n_threads, rand_init = rand_init,
                          id_trace = id_trace, seed = seed, progress = progress)
      private$.multiplicity_engine(dataset = dataset)

      # infer method from cpp model object
      private$.rcpp$infer()
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
      return(private$.refactor_dataset(private$.rcpp$simulate()))
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
      return(private$.rcpp$evaluate())
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
