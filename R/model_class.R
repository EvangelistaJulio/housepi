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
#' Housepi Model
#' @description
#' Compile household model, manage your data
#' and variable defintions from a single class.
#'
#' @docType class
#' @keywords internal
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
    .get_all_required_data = function() {
      required <- c()
      if (!private$.is_blank()) required <- private$.data_required()
      c(
        required, # by the model
        private$.augm$data_required, # by augmented variable
        private$.infec_profile$data_required # by infectivity profile
      )
    },

    # @description 
    # Assert if data_wide headers match implementation file.
    .assert_data_wide = function(headers) {
      if (length(headers) == 0) stop_h(
        "Impossible to proceed. No data was provided."
      )
      required <- private$.get_all_required_data()
      template_missing_stop(required, headers, "data")
      invisible(self)
    },

    # @description
    # Parameters required by implementation file and other definitions.
    .get_all_param = function() {
      if (private$.is_blank()) required <- self$parameters$name
      else required <- private$.param$require(private$.reference())
      c(
        required, # by the model
        private$.group_rates$param_required, # by group rates
        private$.augm$param_required, # by augmented variable
        private$.infec_profile$param_required # by infectivity profile
      )
    },

    # @description
    # Assert if parameter definitions match required
    .assert_param = function() {
      private$.param$assert(private$.get_all_param())
      invisible(self)
    },

    # @description
    # Augmented variables required by implementation file and other definitions.
    .get_all_augm = function() {
      c(
        private$.augm$require(private$.reference()), # by the model
        private$.infec_profile$augm_required # by infectivity profile
      )
    },

    # @description
    # Assert if augmented definitions match the required.
    .assert_augm = function() {
      if (!private$.is_blank()) {
        private$.augm$assert(private$.get_all_augm(),
                             private$.infec_status_set)
      }
      invisible(self)
    },

    # @description
    # Group rates required by implementation file and other definitions.
    .get_all_group_rate = function() {
      private$.group_rates$require(private$.reference())
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
    .multiplicity_engine = function(data_wide = NULL) {
      # check if model allows for multiplicity
      refer_file <- private$.reference()
      regexp <- ".multiplicity\\s*\\(\\s*\\)"

      if (template_assert_expression(refer_file, regexp)) {

        # note: this might be a modeler error
        if (length(private$.get_all_augm()) > 0) warning_h(
          "Model with augmentated variable called multiplicity engine."
        )

        # if data_wide not provided, load previous
        # UPDATE: this "back & forth" of data_wide is not optimal!
        if ((length(data_wide) == 0) && (length(private$.cols_data_wide) > 0)) {
          data_wide <- self$data_wide
        }

        headers <- c(
          # data_wide headers required by implementation
          private$.get_all_required_data(),
          # data_wide headers required by group rates
          unlist(lapply(
            private$.group_rates$logic_tables,
            function(x) unlist(lapply(x, names))
          ))
        )
        # remove logic tables auxiliary column names
        headers <- setdiff(unique(headers), c("n", "param_name"))

        # list of splitted data_wide by house_id
        split_data <- split(data_wide[c("status", headers)], data_wide$house_id)
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
        defs <- data.frame(house_id = unique(data_wide$house_id),
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
    .build_group_rates = function() {
      # assert definitions are complete
      private$.assert_group_rate()

      # proceed only if definitions exist
      if (length(private$.group_rates$defs) > 0) {
        # if data_wide not provided, load previous
        # UPDATE: count in the wrapper while storing data_wide_levels
        if (length(private$.cols_data_wide) == 0) data_wide <- NULL
        else data_wide <- self$data_wide

        # implement group rates
        private$.group_rates$implement(data_wide, private$.info_data_wide)

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
    .build_full = function(data_wide = NULL,
                           sample_size = NA, n_steps = 0,
                           n_burnin = 0, thin = 1, n_adapt = 0,
                           n_chains = 1, n_threads = 1,
                           rand_init = TRUE, id_draws = NA,
                           seed = NA, progress = FALSE) {
      # build data_wide
      self$add_data_wide(data_wide)

      # individual id for blank model augmented draws
      # if no default columns, assume debugging mode
      if (private$.is_blank() && all(is.na(id_draws))) {
        # UPDATE: not robust in case of data_wide with only default columns
        if (length(setdiff(private$.cols_data_wide,
                           c("house_id", "id", "status"))) == 0) {
          id_draws <- template_blank_data_wide_()$id[1]
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
      private$.rcpp$import_group_rates(
        private$.group_rates$expose_to_cpp(private$.info_data_wide)
      )

      # build contact rates (UPDATE: temporary)
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
    .export_to_cpp_standalone = function(outpath = "../housepi/inst/standalone/") {
      write.csv(private$.rcpp$export_data_wide(),
                file.path(outpath, "data_wide.tmp"),
                quote = FALSE, row.names = FALSE)
      write.csv(private$.rcpp$export_data_long(),
                file.path(outpath, "data_long.tmp"),
                quote = FALSE, row.names = FALSE)

      write.csv(private$.group_rates$expose_to_cpp(private$.info_data_wide),
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
    #' @field data_wide `data.frame`\cr
    #' Dataset with individual level information of current data_wide.
    data_wide = function(value) {
      if (missing(value)) {
        if (is.null(private$.rcpp)) stop_h(
          "Impossible to export data_wide before model build."
        )
        return(private$.refactor_data_wide(private$.rcpp$export_data_wide()))
      } else {
        self$add_data_wide(value)
      }
    },

    #' @field data_long `data.frame`\cr
    #' Datalist with multiple rows indexed by individual.
    data_long = function(value) {
      if (missing(value)) {
        if (is.null(private$.rcpp)) stop_h(
          "Impossible to export data_long before model build."
        )
        return(private$.refactor_data_long(private$.rcpp$export_data_long()))
      } else {
        self$add_data_long(value)
      }
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
      else {
        do.call(self$add_infec_profile, value)
      }
    },

    #' @field group_rates `list`\cr
    #' List containing current group rate definitions.
    group_rates = function(value) {
      if (missing(value)) {
        private$.build_group_rates()
        return(private$.group_rates$tidy_table)
      } else {
        self$load_group_rate_list(from = value, append = FALSE)
      }
    }
  ),

  ############################################################################
  ################################## PUBLIC ##################################
  ############################################################################

  public = list(
    #' @description
    #' Create a new object of class 'HousepiModel'
    #' used to run inference and simulation routines.
    #'
    #' @param choice `r rd_mdl_choice_`
    #' @param as_template `r rd_mdl_astmpl_`
    #' @param force `r rd_mdl_force_`
    #' @param quiet `r rd_mdl_quiet_`
    initialize = function(choice = NA, as_template = FALSE,
                          force = FALSE, quiet = TRUE) {
      # initialize internal fields
      self$drop_all()
      # generate model wrapper
      super$initialize(choice = choice, as_template = as_template,
                       force = force, quiet = quiet)
      # import blank template data_wide
      if (private$.is_blank()) self$add_data_wide(template_blank_data_wide_())

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
      open_model_helper_(private$.path())
      invisible(self)
    },

    ############################################################################
    ############################################################################
    ############################################################################

    #' @description
    #' Make an internal copy of the provided household data_wide.
    #'
    #' @param data_wide `r rd_add_data_wide_`
    add_data_wide = function(data_wide) {
      if (length(data_wide) > 0) {
        # check that data_wide type is correct
        if (!is.data.frame(data_wide)) stop_h(
          "Unsupported type for data_wide.",
          "Dataset should be a data.frame."
        )

        # assert data_wide headers
        private$.assert_data_wide(colnames(data_wide))
        # convert factors
        embed <- cpp_data_factorize(data_wide)
        # store levels
        private$.info_data_wide <- embed$info
        # store column names
        private$.cols_data_wide <- embed$cols
        # store all postive status
        private$.infec_status_set <- embed$status

        # import data.frame in rcpp
        private$.rcpp$import_data_wide(embed$data)

      } else {
        # assert headers of last imported data_wide
        private$.assert_data_wide(private$.cols_data_wide)
      }

      # implement group rates and build
      private$.build_group_rates()

      invisible(self)
    },

    #' @description Summarise statuses in the data_wide.
    #' @return`data.frame`
    summarise_cases = function() {
      private$.rcpp$summarise_cases()
    },

    #' @description Summary house sizes in the data_wide.
    #' @return `data.frame`
    summarise_house_sizes = function() {
      private$.rcpp$summarise_house_sizes()
    },

    #' @description
    #' Maje an internal copy of the provided data_long.
    #'
    #' @param data_long `r rd_add_data_long_`
    add_data_long = function(data_long) {
      if (!is.data.frame(data_long)) stop_h(
        "Unsupported type for data_long.",
        "Datalist should be a data.frame."
      )

      # convert factors
      embed <- cpp_data_factorize(data_long)
      # store levels
      private$.info_data_long <- embed$info
      # store column names
      private$.cols_data_long <- embed$cols
      # import data.frame in rcpp
      private$.rcpp$import_data_long(embed$data)

      invisible(self)
    },

    ############################################################################
    ############################################################################
    ############################################################################

    #' @description
    #' Add a new parameter to houspepi model definitions.
    #' If the parameter already exists, it is updated.
    #'
    #' @param object `r rd_model_in_`
    #' @param name `r rd_var_name_`
    #' @param value `r rd_var_value_`
    #' @param const `r rd_var_const_`
    #' @param discr `r rd_par_discr_`
    #' @param prior `r rd_par_model_`
    #' @param args `r rd_par_args_`
    #' @param min `r rd_var_min_`
    #' @param max `r rd_var_max_`
    #' @param rule `r rd_par_rule_`
    #' @param scale `r rd_var_scale_`
    add_parameter = function(name = NA, value = NULL,
                             const = NA, discr = FALSE,
                             prior = NA, args = NULL,
                             min = NULL, max = NULL,
                             rule = NA, scale = NA) {
      private$.param$add(name = name, value = value, const = const,
                         discr = discr, prior = prior, args = args,
                         min = min, max = max,
                         rule = rule, scale = scale)
      invisible(self)
    },

    #' @description
    #' Load parameter definitions either as arguments or from predefined list.
    #'
    #' @param ... `r rd_list_args_`
    #' @param from `r rd_list_from_`
    #' @param append `r rd_list_append_`
    load_parameter_list = function(..., from = NA, append = FALSE) {
      private$.param$load_list(..., from = from, append = append)
      invisible(self)
    },

    #' @description
    #' Load parameter definitions from model template.
    #'
    #' @param append `r rd_tmpl_append_`
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

    #' @description
    #' Add a new augmented variable/function to houspepi model definitions.
    #' If the augmented variablw already exists, it is updated.
    #'
    #' @param name `r rd_var_name_`
    #' @param status `r rd_aug_status_`
    #' @param value `r rd_aug_val_`
    #' @param init_rand `r rd_aug_rndin_`
    #' @param const `r rd_var_const_`
    #' @param discr `r rd_aug_discr_`
    #' @param func `r rd_aug_func_`
    #' @param func_args `r rd_aug_funargs_`
    #' @param model `r rd_aug_model_`
    #' @param args `r rd_aug_args_`
    #' @param min `r rd_aug_min_`
    #' @param max `r rd_aug_max_`
    #' @param rule `r rd_var_rule_`
    #' @param scale `r rd_var_scale_`
    #' @param from `r rd_aug_from_`
    #' @param expr `r rd_aug_expr_`
    add_augmented = function(name = NA, status = NA,
                             const = NA, discr = FALSE,
                             func = NA, func_args = NULL,
                             value = NULL, init_rand = FALSE,
                             model = NA, args = list(),
                             min = NULL, max = NULL,
                             rule = NA, scale = NA, from = NULL,
                             expr = NULL) {
      private$.augm$add(name = name, status = status,
                        const = const, discr = discr,
                        func = func, func_args = func_args,
                        value = value, init_rand = init_rand,
                        model = model, args = args,
                        min = min, max = max,
                        rule = rule, scale = scale, from = from,
                        expr = expr)
      invisible(self)
    },

    #' @description
    #' Load augmented variable/function definitions
    #' either as arguments or from predefined list.
    #'
    #' @param ... `r rd_list_args_`
    #' @param from `r rd_list_from_`
    #' @param append `r rd_list_append_`
    load_augmented_list = function(..., from = NA, append = FALSE) {
      private$.augm$load_list(..., from = from, append = append)
      invisible(self)
    },

    #' @description
    #' Load augmented variable/function definitions from template.
    #'
    #' @param append `r rd_tmpl_append_`
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

    #' @description
    #' Add definitions of a group rate by given data_wide feature(s).
    #' Reference value has unitary rate and can be specified
    #' (default is most represented group).
    #'
    #' @param name `r rd_var_name_`
    #' @param by `r rd_grr_by_`
    #' @param prefix `r rd_grr_prefix_`
    #' @param reference `r rd_grr_refer_`
    #' @param independent `r rd_grr_indep_`
    #' @param prior `r rd_par_model_`
    #' Default is "lnorm".
    #' @param args `r rd_par_args_`
    #' Default is `list(meanlog = 0, sdlog = 3)`.
    #' @param min `r rd_var_min_`
    #' Default is 0.
    #' @param max `r rd_var_max_`
    #' Default is 10.
    #' @param rule `r rd_var_rule_`
    #' @param scale `r rd_var_scale_`
    add_group_rate = function(name, by,
                              prefix = NA,
                              reference = NULL, independent = TRUE,
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
    #' @param ... ` r rd_list_args_`
    #' @param from `r rd_list_from_`
    #' @param append `r rd_list_append_`
    load_group_rate_list = function(..., from = NA, append = FALSE) {
      if (isFALSE(append)) self$drop_group_rates()
      private$.group_rates$load_list(..., from = from, append = append)
      invisible(self)
    },

    ############################################################################
    ############################################################################
    ############################################################################

    #' @description
    #' Add definition of infectivity profile.
    #'
    #' @param object `r rd_model_in_`
    #' @param name `r rd_var_name_`
    #' @param discr `r rd_inf_discr_`
    #' @param args `r rd_inf_args_`
    #' @param max `r rd_aug_max_`
    add_infec_profile = function(name, discr = FALSE,
                                 args = NULL, max = NULL) {
      private$.infec_profile$add(name = name, discr = discr,
                                 args = args, max = max)
      invisible(self)
    },

    ############################################################################
    ############################################################################
    ############################################################################

    #' @description
    #' Compute the loglikelihood.
    #'
    #' @param object `r rd_model_in_`
    #' @param data_wide `r rd_smpl_data_wide_`
    #' @param seed `r rd_seed_`
    #'
    #' @return `numeric(1)`
    loglikelihood = function(data_wide = NULL, seed = NA) {
      # internal build function
      private$.build_full(data_wide = data_wide, seed = seed, rand_init = FALSE)
      private$.multiplicity_engine(data_wide = data_wide)
      return(private$.rcpp$evaluate())
    },

    #' @description
    #' Run inference with Bayesian Markov Chain Monte Carlo.
    #'
    #' @param data_wide `r rd_smpl_data_wide_`
    #' @param sample_size `r rd_smpl_size_`
    #' @param n_steps `r rd_smpl_steps_`
    #' @param n_burnin `r rd_smpl_burn_`
    #' @param thin `r rd_smpl_thin_`
    #' @param n_adapt `r rd_smpl_adapt_`
    #' @param n_chains `r rd_smpl_chain_`
    #' @param n_threads `r rd_smpl_thread_`
    #' @param rand_init `r rd_smpl_rndin_`
    #' @param id_draws `r rd_smpl_iddrw_`
    #' @param seed `r rd_seed_`
    #' @param progress `r rd_progress_`
    infer = function(data_wide = NULL,
                     sample_size = NA, n_steps = 10000,
                     n_burnin = 0, thin = 1, n_adapt = 0,
                     n_chains = 1, n_threads = 1,
                     rand_init = TRUE, id_draws = NA,
                     seed = NA, progress = TRUE) {
      # internal build function
      private$.build_full(data_wide = data_wide, sample_size = sample_size,
                          n_steps = n_steps, n_burnin = n_burnin, thin = thin,
                          n_adapt = n_adapt, n_chains = n_chains,
                          n_threads = n_threads, rand_init = rand_init,
                          id_draws = id_draws, seed = seed, progress = progress)
      private$.multiplicity_engine(data_wide = data_wide)

      # infer method from cpp model object
      private$.rcpp$infer()

      return(HousepiPosteriorDraws$new(private$.rcpp))
    },

    #' @description
    #' Estimate the model evidence (or marginal likelihood).
    #' It uses the method of power posteriors with different available schemes.
    #'
    #' @param data_wide `r rd_smpl_data_wide_`
    #' @param scheme `r rd_evi_scheme_`
    #' @param n_rungs `r rd_evi_rungs_`
    #' @param sample_size `r rd_smpl_size_`
    #' @param n_steps `r rd_smpl_steps_`
    #' @param n_burnin `r rd_smpl_burn_`
    #' @param thin `r rd_smpl_thin_`
    #' @param n_adapt `r rd_smpl_adapt_`
    #' @param n_chains `r rd_smpl_chain_`
    #' @param n_threads `r rd_smpl_thread_`
    #' @param seed `r rd_seed_`
    #' @param progress `r rd_progress_`
    #'
    #' @return `numeric(1)`
    evidence = function(data_wide = NULL,
                        scheme = "power-fractions", n_rungs = 40,
                        sample_size = NA, n_steps = 10000,
                        n_burnin = 0, thin = 1, n_adapt = 0,
                        n_chains = 1, n_threads = 1,
                        seed = NA, progress = TRUE) {
      # internal build function
      private$.build_full(data_wide, sample_size = sample_size,
                          n_steps = n_steps, thin = thin,
                          n_burnin = n_burnin, n_adapt = n_adapt,
                          n_chains = n_chains, n_threads = n_threads,
                          rand_init = TRUE, seed = seed, progress = progress)
      private$.multiplicity_engine(data_wide = data_wide)

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
    },

    #' @description
    #' Dataset with associated augmented variables at current state.
    #'
    #' @return
    #' `data.frame`
    export_augmented_data = function() {
      private$.refactor_data_wide(
        cbind(private$.rcpp$data_wide_for_simulation(),
              private$.rcpp$export_augmented_data_wide())
      )
    },

    #' @description
    #' A snapshot of all current internal state of a given household.
    #'
    #' @param house_id `character(1)`\cr
    #' Id of the household to export.
    #'
    #' @return `list`
    export_house_snapshot = function(house_id) {
      house_id <- as.character(house_id)
      output <- private$.rcpp$export_house_snapshot(house_id)
      output$data_wide <- private$.refactor_data_wide(output$data_wide)
      output$data_long <- private$.refactor_data_long(output$data_long)

      # first draft to reformat infectivity profile output
      lst <- output$infec_profile
      df <- data.frame(id = lst$id, discr = lst$discr, name = lst$name,
                       stringsAsFactors = FALSE)
      df$args <- lst$args
      df$range <- lst$range
      output$infec_profile <- df

      output
    },

    #' @description
    #' Generate a simulated outbreak with the same community structure of the loaded data.
    #'
    #' @param nsim `r rd_smpl_nsim_`
    #' @param seed `r rd_seed_`
    #' @param data_wide `r rd_smpl_data_wide_`
    #' @param progress `r rd_progress_`
    #'
    #' @return `data.frame`
    simulate = function(nsim = 1,
                        seed = NA,
                        data_wide = NULL,
                        progress = FALSE) {
      # internal build function
      private$.build_full(data_wide = data_wide, rand_init = FALSE,
                          seed = seed, progress = progress)

      # refactor static data_wide of simulation
      static_df <- private$.rcpp$data_wide_for_simulation()

      # bind simulated entries with static part
      # UPDATE : n_rep = 1 for development
      private$.refactor_data_wide(
        cbind(static_df, private$.rcpp$simulate(nsim))
      )
    },

    ############################################################################
    ############################################################################
    ############################################################################

    #' @description
    #' Drop parameter definitions.
    drop_parameters = function() {
      private$.param <- HousepiParameter$new()
      invisible(self)
    },

    #' @description
    #' Drop augmented variable definitions.
    drop_augmented = function() {
      private$.augm <- HousepiAugmented$new()
      invisible(self)
    },

    #' @description
    #' Drop group rates definitions.
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

    #' @description
    #' Drop all definitions.
    drop_all = function() {
      self$drop_parameters()
      self$drop_augmented()
      self$drop_group_rates()
      private$.infec_profile <- HousepiInfecProfile$new()
      invisible(self)
    }
  )
)
