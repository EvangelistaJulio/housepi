## housepi - Epidemiological Analysis of Household Outbreaks.
##
## Copyright © 2026 Institut Pasteur (Paris), CNRS and INSERM.
## Distributed under the terms of the GNU General Public License (GPLv3).
## See the LICENSE file for details.
##

############################################################################
############################################################################
############################## Household Model #############################
############################################################################
############################################################################

# auxiliary environment
.housepi_env <- new.env(parent = emptyenv())
.housepi_env$hollow <- FALSE

#' @title
#' R6 class with the household model
#' @description
#' Compile the household model, manage your data
#' and variable definitions from a single class.
#' @docType class
HousepiModel <- R6::R6Class(
  classname = "HousepiModel",
  inherit = HousepiGenerator,
  cloneable = FALSE,
  lock_class = TRUE,

  ############################################################################
  ################################## PRIVATE #################################
  ############################################################################

  private = list(
    param_ = NULL,
    augm_ = NULL,
    group_rates_ = NULL,
    contact_rates_ = NULL,
    infec_profile_ = NULL,

    ############################################################################
    ############################################################################
    ############################################################################

    # @description
    # Augmented variables required by implementation file and other definitions.
    get_all_required_augm_ = function() {
      c(private$augm_$require(private$get_reference_()), # by the model
        private$infec_profile_$augm_required)            # by infec profile
    },

    # @description
    # Dataset headers required by implementation file and other definitions.
    get_all_required_data_wide_ = function() {
      refer_file <- private$get_reference_()
      output <- template_extract_keys(refer_file, "DataKey")
      required <- set_source_name_(output, "model")
      c(required,                               # by the model
        private$group_rates_$data_required,     # by group rates
        private$contact_rates_$data_required,   # by contact rates
        private$augm_$data_required,            # by augmented
        private$infec_profile_$data_required)   # by infec profile
    },

    # @description
    # Group rates required by implementation file and other definitions.
    get_all_required_group_rates_ = function() {
      c(private$group_rates_$require(private$get_reference_()), # by the model
        private$augm_$group_rates_required,                 # by augmented
        private$infec_profile_$group_rates_required)        # by infec profile
    },

    # @description
    # Parameters required by implementation file and other definitions.
    get_all_required_param_ = function() {
      required <- private$param_$require(private$get_reference_())
      c(required,                                # by the model
        private$param_$param_required,           # by parameter functions
        private$group_rates_$param_required,     # by group rates
        private$contact_rates_$param_required,   # by contact rates
        private$augm_$param_required,            # by augmented
        private$infec_profile_$param_required)   # by infec profile
    },

    ############################################################################
    ############################################################################
    ############################################################################

    # @description
    # Assert if augmented definitions match the required.
    # @param capital `logical(1)` to raise error.
    assert_augm_ = function(capital = TRUE, warn_unused = TRUE) {
      private$augm_$assert(private$get_all_required_augm_(),
                           private$infec_status_set_,
                           capital = capital,
                           warn_unused = warn_unused && !private$is_blank_())
    },

    # @description
    # Assert if data_wide headers match implementation file.
    # @param capital `logical(1)` to raise error.
    assert_data_wide_ = function(capital = TRUE) {
      if (isTRUE(capital) && length(private$cols_data_wide_) == 0) stop_h(
        "Impossible to proceed. No data was provided."
      )
      required <- private$get_all_required_data_wide_()
      sms <- template_missing_stop(required, private$cols_data_wide_, "data")
      if (isTRUE(capital) && length(sms) > 0) stop_h(sms)
      else sms
    },

    # @description
    # Assert group rate definitions from the model implementation script.
    # @param capital `logical(1)` to raise error.
    assert_group_rates_ = function(capital = TRUE) {
      private$group_rates_$assert(private$get_all_required_group_rates_(),
                                  capital = capital)
    },

    # @description
    # Assert contact rate definitions are compatible with the model.
    assert_contact_rates_ = function() {
      has_defs <- private$contact_rates_$assert()
      if (has_defs && !.engine$has_contact_rates(private$mcmc_ptr_)) {
        self$drop_contact_rates()
        stop_h(
          sprintf("Contact rates are not supported for model '%s'.",
                  self$name),
          "Definitions will be dropped."
        )
      }
    },

    # @description
    # Assert if parameter definitions match required
    # @param capital `logical(1)` to raise error.
    # @param warn_unused `logical(1)` warn if unused.
    assert_param_ = function(capital = TRUE, warn_unused = TRUE) {
      private$param_$assert(private$get_all_required_param_(),
                            capital = capital,
                            warn_unused = warn_unused && !private$is_blank_())
    },

    ############################################################################
    ############################################################################
    ############################################################################

    # @description
    # Get template or previously implemented infectiviyt profile.
    get_internal_infec_profile_ = function() {
      .engine$export_infec_profile(private$mcmc_ptr_)
    },

    ############################################################################
    ############################################################################
    ############################################################################

    # @description
    # Helper to assert and re-build group and contact rates.
    # @param family
    helper_build_rates_ = function(module) {
      # proceed only if definitions exist
      if (length(module$defs) == 0) return(invisible(self))

      # UPDATE: this "back & forth" of data_wide is sub-optimal...
      module$implement(self$data_wide, private$info_data_wide_)

      # append parameter definitions if missing
      new_parameters <- module$parameters
      # find missing parameter definitions
      mask <- setdiff(names(new_parameters), private$param_$defs$name)
      new_parameters <- new_parameters[mask]
      if (length(new_parameters) > 0) {
        xquos <- convert_list_to_quosures(new_parameters)
        private$param_$load_quosures(xquos, append = TRUE)
      }

      invisible(self)
    },

    # @description
    # Assert and pre-build group rates.
    build_group_rates_ = function() {
      private$assert_group_rates_()
      private$helper_build_rates_(private$group_rates_)
    },

    # @description
    # Assert and pre-build contact rates.
    build_contact_rates_ = function() {
      private$assert_contact_rates_()
      private$helper_build_rates_(private$contact_rates_)
    },

    # @description
    # Get individual for blank model augmented draws.
    # UPDATE: not robust if data has only default columns
    get_blank_id_draws_ = function(id_draws) {
      # if no default columns, assume debugging mode
      if (private$is_blank_() &&
            is.null(id_draws) &&
            !is.null(private$augm_$defs) &&
            length(setdiff(private$cols_data_wide_,
                           c("house_id", "id", "status"))) == 0) {
        env <- new.env(parent = emptyenv())
        utils::data("blank_house", package = "housepi", envir = env)
        id_draws <- env$blank_house$id[1]
      }
      id_draws
    },

    # @description
    # Validate, factorize, and push `data` (wide format) into the C++ engine.
    import_data_wide_ = function(data) {
      if (!is.data.frame(data)) stop_h(
        "Data should be an object of class data.frame.",
        format_type_error(data)
      )
      if (dim(data)[[1]] == 0) stop_h("Provided is empty.")

      # convert factors
      embed <- cpp_data_factorize(data)
      # store levels
      private$info_data_wide_ <- embed$info
      # store column names
      private$cols_data_wide_ <- embed$cols
      # store all positive status
      private$infec_status_set_ <- embed$status

      # assert data_wide headers
      private$assert_data_wide_()

      # import data.frame in rcpp
      .engine$import_data_wide(private$mcmc_ptr_, embed$data)

      invisible(self)
    },

    # @description
    # Tolerabt re-push of R-side definitions to the C++ engine.
    # Hard validation lives in `build_full_`.
    rebuild_ = function() {
      # guard: pre-init (before super$initialize creates mcmc_ptr_)
      if (is.null(private$mcmc_ptr_)) return(invisible(self))
      # guard: data_wide not yet imported (during initialize drop_all phase)
      if (length(private$cols_data_wide_) == 0) return(invisible(self))

      # soft data_wide check (headers vs. implementation file)
      private$assert_data_wide_(FALSE)

      # build rates first (they append to param_$defs)
      if (length(private$group_rates_$defs) > 0) {
        private$build_group_rates_()
      }
      if (length(private$contact_rates_$defs) > 0) {
        private$build_contact_rates_()
      }

      #> parameters
      if (length(private$param_$defs) > 0) private$assert_param_(FALSE, FALSE)
      .engine$import_param(private$mcmc_ptr_, private$param_$defs)

      #> group rates
      .engine$import_group_rate(
        private$mcmc_ptr_,
        private$group_rates_$expose_to_cpp(private$info_data_wide_)
      )

      #> contact rates
      .engine$import_contact_rate(
        private$mcmc_ptr_,
        private$contact_rates_$expose_to_cpp(private$info_data_wide_)
      )

      #> infectivity profile
      .engine$import_infec_profile(private$mcmc_ptr_,
                                   private$infec_profile_$defs)

      #> augmented
      if (length(private$augm_$defs) > 0) private$assert_augm_(FALSE, FALSE)
      .engine$import_augm(private$mcmc_ptr_, private$augm_$defs)

      invisible(self)
    },

    # @description
    # Full model build. Refresh runtime setup, enforce hard validation,
    # then push the final state to the engine.
    build_full_ = function(steps = 0, adapt = 0, burn = 0,
                           thin = 1, chains = 1, threads = 1,
                           id_draws = NULL, progress = FALSE) {
      .engine$import_setup(private$mcmc_ptr_, char_setup_input_(list(
        steps = as.integer(steps),
        adapt = as.integer(adapt),
        burn = as.integer(burn),
        thin = as.integer(thin),
        chains = as.integer(chains),
        threads = as.integer(threads),
        id_draws = private$get_blank_id_draws_(id_draws),
        # generate a seed for Cpp random number generator
        seed = sample.int(.Machine$integer.max, 1L),
        progress = as.integer(progress)
      )))
      # hard validation before MCMC
      private$assert_data_wide_(capital = TRUE)
      private$assert_param_(capital = TRUE)
      private$assert_group_rates_(capital = TRUE)
      private$assert_contact_rates_()
      private$assert_augm_(capital = TRUE)
      # push final state (picks up any stage deferred by dep gating)
      private$rebuild_()
      invisible(self)
    },

    # Helper with parameter table import.
    # @param ptab Table with parameter values
    # param ess `logical(1)`
    helper_read_param_table_ = function(ptab, ess = TRUE) {
      # extract parameter draws from table
      if (inherits(ptab, "HousepiPosteriorDraws")) ptab <- ptab$param_draws

      if (posterior::is_draws(ptab)) {
        ptab <- posterior::as_draws_matrix(ptab)
        # reduce sample to effective sample size (ESS)
        if (isTRUE(ess)) {
          ess_bulk <- sapply(posterior::variables(ptab),
                             function(v) posterior::ess_bulk(ptab[, v]))
          # check no NA values in ESS
          if (any(is.na(ess_bulk))) stop_h(
            "ESS values return NA. Check your draws for degenerate parameters."
          )
          # suppress posterior message
          target_ndraws <- min(ess_bulk)
          thin <- floor(posterior::ndraws(ptab) / target_ndraws)
          ptab <- posterior::thin_draws(ptab, thin = thin)
        }
      }
      # UPDATE:  warn if sample is too small!

      tryCatch({
        return(as.data.frame(ptab))
      }, error = function(e) stop_h("Type of param_table not supported.")
      )
    },

    # Helper with simulation build.
    # FIXME: If there's a param table, build should change
    helper_simulation_build_ = function(param_table = NULL,
                                        nsim = NA_integer_,
                                        ess = TRUE) {
      # read param table
      if (is.null(param_table)) {
        if (is.na(nsim)) nsim <- 1
        else if (!is.numeric(nsim) || nsim < 0) stop_h(
          "Number of simulations should be a positive integer."
        )
      } else {
        param_table <- private$helper_read_param_table_(param_table, ess = ess)
        nsim <- dim(param_table)[1]
        # don't overwrite
        if (dim(self$parameters)[1] == 0) {
          param_dummy_list <- as.list(utils::head(param_table, 1))
          do.call(self$load_parameter_list, param_dummy_list)
        }
      }

      # internal build function
      private$build_full_(steps = nsim, progress = FALSE)

      # NOTE: there might be not parameters at all
      if (dim(self$parameters)[1] == 0) return(nsim)

      if (is.null(param_table)) {
        # replicates the single row nsim times
        param_table <- as.data.frame(stats::setNames(
          as.list(self$parameters$value), self$parameters$name
        ))
        param_table <- param_table[rep(1L, nsim), , drop = FALSE]
        rownames(param_table) <- NULL
      }
      # attach table with parameters
      .engine$import_param_table(private$mcmc_ptr_, param_table)

      nsim
    },

    # Helper with observation process import.
    # Combine observation process to measures and refactor.
    # UPDATE: I shall drop necessity for refactoring with variant data.
    helper_obs_process_ = function(obs_process) {
      if (is.null(obs_process)) {
        private$refactor_data_wide_
      } else if (is.function(obs_process)) {
        function(x) obs_process(private$refactor_data_wide_(x))
      } else {
        stop_h(
          "Observation process must be a function.",
          format_type_error(obs_process)
        )
      }
    },

    # @description
    # Run measure function on posterior simulations.
    # UPDATE: this has to be improved with loops in C++
    run_measure_ = function(param_table,
                            apply_wide,
                            apply_long = NULL,
                            obs_proc_wide = NULL,
                            obs_proc_long = NULL,
                            ess = TRUE,
                            threads = 1,
                            progress = TRUE) {
      param_table <- private$helper_read_param_table_(param_table, ess = ess)
      nsim <- dim(param_table)[1]

      # internal build function
      private$build_full_(steps = nsim, threads = threads, progress = progress)

      # attach param_table
      .engine$import_param_table(private$mcmc_ptr_, param_table)

      if (!is.null(apply_long) &&
            (is.null(obs_proc_long) || !is.function(obs_proc_long))) stop_h(
        "'obs_process_long' must be a function."
      )
      .engine$run_simulation(
        private$mcmc_ptr_,
        private$helper_obs_process_(obs_proc_wide),
        apply_wide,
        if (is.null(apply_long)) NULL
        else function(x) apply_long(obs_proc_long(x))
      )
    },

    # @description
    # Estimate the model logevidence (or marginal likelihood).
    # It uses the method of power posteriors with different available schemes.
    #
    # @param scheme `r rd_evi_scheme_`
    # @param rungs `r rd_evi_rungs_`
    # @param steps `r rd_smpl_steps_`
    # @param burn `r rd_smpl_burn_`
    # @param thin `r rd_smpl_thin_`
    # @param adapt `r rd_smpl_adapt_`
    # @param chains `r rd_smpl_chain_`
    # @param threads `r rd_thread_`
    # @param progress `r rd_progress_`
    #
    # @return `numeric(1)`
    logevidence_ = function(scheme = "power-fractions", rungs = 40,
                            steps = 10000, burn = 0,
                            thin = 1, adapt = 0,
                            chains = 1, threads = 1,
                            progress = TRUE) {
      warning_h(
        "[Internal] This method is currently under development.",
        "Its usage is deprecated."
      )

      # internal build function
      private$build_full_(steps = steps, burn = burn,
                          thin = thin, adapt = adapt,
                          chains = chains, threads = threads,
                          progress = progress)
      # load options
      scheme <- as.character(scheme)
      rungs <- as.integer(rungs)
      # evaluate marginal from the internal Rcpp object
      output <- .engine$logevidence(private$mcmc_ptr_, scheme, rungs)
      # reformat dataframe
      mask <- c("inv_temper", "mean_log_lik", "var_log_lik")
      output$detail <- as.data.frame(output[mask])
      output[mask] <- NULL
      output
    }
  ),

  ############################################################################
  ################################## ACTIVE ##################################
  ############################################################################

  active = list(
    #' @field data_wide `data.frame`. Individual level data (wide format).
    data_wide = function(value) {
      if (!missing(value)) stop_read_only("data_wide")
      private$refactor_data_wide_(.engine$export_data_wide(private$mcmc_ptr_))
    },

    #' @field data_long
    #' `data.frame`. Additional data rows indexed by individual (long format).
    data_long = function(value) {
      if (!missing(value)) stop_read_only("data_long")
      private$refactor_data_long_(.engine$export_data_long(private$mcmc_ptr_))
    },

    #' @field parameters `data.frame`. Parameter definitions.
    parameters = function(value) {
      if (!missing(value)) stop_read_only("parameters")
      housepi_defin_showcaser(private$param_$defs)
    },

    #' @field augmented `data.frame`. Augmented definitions.
    augmented = function(value) {
      if (!missing(value)) stop_read_only("augmented")
      housepi_defin_showcaser(private$augm_$defs)
    },

    #' @field infec_profile `data.frame`. Infectivity profile definitions.
    infec_profile = function(value) {
      if (!missing(value)) stop_read_only("infec_profile")
      # get custom definition
      custom <- private$infec_profile_$showcase()
      if (length(custom) > 0) return(custom)
      # get template/previously implemented
      stored <- private$get_internal_infec_profile_()
      if (length(stored) == 0) return(NULL)
      # format exported implementation
      output <- data.frame(func = stored$func,
                           func_args = base::I(list(stored$func_args)),
                           min = stored$min, max = stored$max,
                           discr = stored$discr,
                           stringsAsFactors = FALSE)
      class(output$func_args) <- "list"
      output
    },

    #' @field group_rates `data.frame`. Group rate definitions.
    group_rates = function(value) {
      if (!missing(value)) stop_read_only("group_rates")
      private$group_rates_$tidy_table
    },

    #' @field contact_rates `data.frame`. Contact rate definitions.
    contact_rates = function(value) {
      if (!missing(value)) stop_read_only("contact_rates")
      private$contact_rates_$tidy_table
    }
  ),

  ############################################################################
  ################################## PUBLIC ##################################
  ############################################################################

  public = list(
    #' @description
    #' Create a new object of class 'HousepiModel'
    #' used to run inference and simulation routines.
    #' @param data `r rd_add_wide_`
    #' @param model `r rd_mdl_choice_`
    #' @param as_template `r rd_mdl_astmpl_`
    #' @param force `r rd_mdl_force_`
    #' @param quiet `r rd_mdl_quiet_`
    #' @param ... `r rd_mdl_flags_`
    initialize = function(data = NULL, model = NA_character_,
                          as_template = FALSE,
                          force = FALSE, quiet = FALSE, ...) {
      # helper construction option fort copying
      if (isTRUE(.housepi_env$hollow)) return(invisible(self))

      # initialize internal fields
      self$drop_all()
      # generate model wrapper
      super$initialize(model = model, as_template = as_template,
                       force = force, quiet = quiet, ...)
      # import data_wide (invariant: always populated after construction)
      if (!is.null(data)) {
        private$import_data_wide_(data)
      } else if (private$is_blank_()) {
        env <- new.env(parent = emptyenv())
        utils::data("blank_house", package = "housepi", envir = env)
        private$import_data_wide_(env$blank_house)
      } else {
        stop_h(sprintf(
          "Missing household 'data' argument for model {.pkg %s}.",
          self$name
        ))
      }

      invisible(self)
    },

    #' @description Deep copy this object.
    #' @return A new independent object.
    copy = function() {
      .housepi_env$hollow <- TRUE
      on.exit(.housepi_env$hollow <- FALSE)
      obj <- HousepiModel$new()

      # clone default fields
      obj_encl <- obj$.__enclos_env__$private
      for (field in names(obj_encl)) {
        val <- private[[field]]
        if (!field %in% c("mcmc_ptr_") && !is.function(val)) {
          obj_encl[[field]] <- helper_clone(val)
        }
      }

      # use cloned paths to create new Rcpp pointers
      obj_encl$create_new_pointers_()

      # transfer setup
      .engine$import_setup(
        obj_encl$mcmc_ptr_,
        char_setup_input_(.engine$export_setup(private$mcmc_ptr_))
      )
      # share data store (just copy-on-write)
      .engine$share_data_store(obj_encl$mcmc_ptr_, private$mcmc_ptr_)

      # fresh pointer needs definitions pushed to the engine
      obj_encl$rebuild_()

      obj
    },

    #' @description
    #' Override the default printing behavior.
    print = function() {
      sms <- c(i = sprintf("Household model {.pkg %s}%s.", self$name,
                           ifelse(private$is_template_(), " (template)", "")))
      # UPDATE: to be rationalize in the future if more options are included
      if ("_HOUSEPI_DISCRETE_TIME" %in% private$last_flags_) sms <- c(
        sms, "Model is set to {.emph discrete time} mode."
      )
      if (length(private$get_all_required_data_wide_()) > 0) sms <- c(sms,
        fmt_assert_cli("Data", private$assert_data_wide_(FALSE))
      )
      if (length(private$get_all_required_param_()) > 0) sms <- c(sms,
        fmt_assert_cli("Parameters", private$assert_param_(FALSE))
      )
      if (length(private$get_all_required_augm_()) > 0) sms <- c(sms,
        fmt_assert_cli("Augmented", private$assert_augm_(FALSE))
      )
      if (length(private$get_all_required_group_rates_()) > 0) sms <- c(sms,
        fmt_assert_cli("Group rates", private$assert_group_rates_(FALSE),
                       minor = TRUE)
      )
      if (length(private$get_internal_infec_profile_()) > 0) sms <- c(sms,
        fmt_assert_cli(
          "Infectivity profile",
          if (length(private$infec_profile_$showcase()) == 0) c(
            "Using templated profile. A custom choice can be provided."
          )
          else NULL,
          minor = TRUE
        )
      )
      cli::cli_inform(sms)
      invisible(self)
    },

    #' @description
    #' Show model description in RStudio Help window.
    help = function() {
      open_model_helper_(private$source_dir_)
      invisible(self)
    },

    ############################################################################
    ############################################################################
    ############################################################################

    #' @description Summarise statuses in loaded data.
    #' @return`data.frame`
    summarise_cases = function() {
      .engine$summarise_cases(private$mcmc_ptr_)
    },

    #' @description Summary house sizes in loaded data.
    #' @return `data.frame`
    summarise_house_sizes = function() {
      .engine$summarise_house_sizes(private$mcmc_ptr_)
    },

    #' @description Store a copy of provided individual data.
    #' @param data_long `r rd_add_data_long_`
    add_data_long = function(data_long) {
      if (!is.data.frame(data_long)) stop_h(
        "Data-long should be an object of class data.frame.",
        format_type_error(data_long)
      )

      # convert factors
      embed <- cpp_data_factorize(data_long)
      # store levels
      private$info_data_long_ <- embed$info
      # store column names
      private$cols_data_long_ <- embed$cols
      # import data.frame in rcpp
      .engine$import_data_long(private$mcmc_ptr_, embed$data)

      private$rebuild_()
      invisible(self)
    },

    ############################################################################
    ############################################################################
    ############################################################################

    #' @description
    #' Add a new parameter to housepi model definitions.
    #' If the parameter already exists, it is updated.
    #' @param name `r rd_var_name_`
    #' @param value `r rd_var_value_`
    #' @param const `r rd_var_const_`
    #' @param discr `r rd_par_discr_`
    #' @param prior `r rd_par_model_`
    #' @param ... `r rd_par_args_`
    #' @param min `r rd_var_min_`
    #' @param max `r rd_var_max_`
    #' @param proposal `r rd_par_proposal_`
    #' @param proposal_scale `r rd_proposal_scale_`
    add_parameter = function(name = NA_character_, value = NULL,
                             const = NA, discr = FALSE,
                             prior = NA_character_, ...,
                             min = NULL, max = NULL,
                             proposal = NA_character_,
                             proposal_scale = NA_real_) {
      private$param_$add_quo(name = name, value = rlang::enquo(value),
                             const = const, discr = discr,
                             prior = prior, rlang::enquos(...),
                             min = rlang::enquo(min), max = rlang::enquo(max),
                             proposal = proposal,
                             proposal_scale = proposal_scale,
                             .assert = TRUE)
      private$rebuild_()
      invisible(self)
    },

    #' @description
    #' Load parameter definitions as named arguments.
    #' @param ... `r rd_list_args_`
    #' @param append `r rd_list_append_`
    load_parameter_list = function(..., append = FALSE) {
      private$param_$load_quosures(rlang::enquos(...), append = append)
      private$rebuild_()
      invisible(self)
    },

    #' @description Load parameter definitions from template.
    #' @param append `r rd_tmpl_append_`
    load_parameter_template = function(append = FALSE) {
      if (private$is_blank_()) warning_h(
        "No template parameters available for blank model."
      )
      # retrieve the template parameters file
      else self$load_parameter_list(
        !!!get_template_def_list(private$source_dir_, "parameter"),
        append = isTRUE(append)
      )
      invisible(self)
    },

    ############################################################################
    ############################################################################
    ############################################################################

    #' @description
    #' Add a new augmented variable/function to housepi model definitions.
    #' If the augmented variable already exists, it is updated.
    #' @param name `r rd_var_name_`
    #' @param value `r rd_aug_val_`
    #' @param status `r rd_aug_status_`
    #' @param const `r rd_var_const_`
    #' @param discr `r rd_aug_discr_`
    #' @param func `r rd_aug_func_`
    #' @param func_args `r rd_aug_funargs_`
    #' @param model `r rd_aug_model_`
    #' @param ... `r rd_aug_args_`
    #' @param min `r rd_aug_min_`
    #' @param max `r rd_aug_max_`
    #' @param proposal `r rd_var_proposal_`
    #' @param proposal_scale `r rd_proposal_scale_`
    add_augmented = function(name = NA_character_,  value = NULL,
                             status = NA_integer_, const = NA, discr = FALSE,
                             func = NULL, func_args = NULL,
                             model = NA_character_, ...,
                             min = NULL, max = NULL,
                             proposal = NA_character_,
                             proposal_scale = NA_real_) {
      private$augm_$add_quo(name = name, value = rlang::enquo(value),
                            status = status, const = const, discr = discr,
                            func = rlang::enquo(func),
                            func_args = rlang::enquos(func_args),
                            model = model, rlang::enquos(...),
                            min = rlang::enquo(min), max = rlang::enquo(max),
                            proposal = proposal,
                            proposal_scale = proposal_scale,
                            .assert = TRUE)
      private$rebuild_()
      invisible(self)
    },

    #' @description Load augmented definitions as named arguments.
    #' @param ... `r rd_list_args_`
    #' @param append `r rd_list_append_`
    load_augmented_list = function(..., append = FALSE) {
      private$augm_$load_quosures(rlang::enquos(...), append = append)
      private$rebuild_()
      invisible(self)
    },

    #' @description Load augmented definitions from template.
    #' @param append `r rd_tmpl_append_`
    load_augmented_template = function(append = FALSE) {
      if (private$is_blank_()) warning_h(
        "No template augmented variables available for blank model."
      )
      # retrieve the template augmented file
      else self$load_augmented_list(
        !!!get_template_def_list(private$source_dir_, "augmented"),
        append = isTRUE(append)
      )
      invisible(self)
    },

    ############################################################################
    ############################################################################
    ############################################################################

    #' @description `r rd_grr_desc_`
    #' @param name `r rd_var_name_`
    #' @param by `r rd_grr_by_`
    #' @param prefix `r rd_grr_prefix_`
    #' @param reference `r rd_grr_refer_`
    #' @param joint `r rd_grr_joint_`
    #' @param strata `r rd_grr_strata_`
    #' @param params `r rd_grr_params_`
    #' @param prior `r rd_par_model_`
    #' @param ... `r rd_par_args_`
    #' @param min `r rd_var_min_`
    #' @param max `r rd_var_max_`
    #' @param proposal `r rd_var_proposal_`
    #' @param proposal_scale `r rd_proposal_scale_`
    add_group_rate = function(name,
                              by = NA_character_, prefix = NA_character_,
                              reference = NULL, joint = FALSE,
                              strata = NULL, params = NULL,
                              prior = NA_character_, ...,
                              min = NA_real_, max = NA_real_,
                              proposal = NA_character_,
                              proposal_scale = NA_real_) {
      private$group_rates_$add_quo(name = name, by = by, prefix = prefix,
                                   reference = reference, joint = joint,
                                   strata = strata, params = params,
                                   prior = prior, rlang::enquos(...),
                                   min = rlang::enquo(min),
                                   max = rlang::enquo(max),
                                   proposal = proposal,
                                   proposal_scale = proposal_scale)
      # UPDATE: I didn't find a way to assert prior arguments here
      private$rebuild_()
      invisible(self)
    },

    #' @description
    #' Load group rate definitions as direct arguments.
    #' @param ... ` r rd_list_args_`
    #' @param append `r rd_list_append_`
    load_group_rate_list = function(..., append = FALSE) {
      if (isFALSE(append)) self$drop_group_rates()
      private$group_rates_$load_quosures(rlang::enquos(...), append = append)
      # UPDATE: I didn't find a way to assert prior arguments here
      private$rebuild_()
      invisible(self)
    },

    #' @description `r rd_ctr_desc_`
    #' @param by `r rd_ctr_by_`
    #' @param reference `r rd_ctr_refer_`
    #' @param prefix `r rd_ctr_prefix_`
    #' @param joint `r rd_ctr_joint_`
    #' @param symmetric `r rd_ctr_symmetric_`
    #' @param strata `r rd_ctr_strata_`
    #' @param params `r rd_ctr_params_`
    #' @param prior `r rd_par_model_`
    #' @param ... `r rd_par_args_`
    #' @param min `r rd_var_min_`
    #' @param max `r rd_var_max_`
    #' @param proposal `r rd_var_proposal_`
    #' @param proposal_scale `r rd_proposal_scale_`
    add_contact_rate = function(by = NA_character_,
                                reference = NULL, prefix = NA_character_,
                                joint = TRUE, symmetric = TRUE,
                                strata = NULL, params = NULL,
                                prior = NA_character_, ...,
                                min = NA_real_, max = NA_real_,
                                proposal = NA_character_,
                                proposal_scale = NA_real_) {
      private$contact_rates_$add_quo(
        by = by, reference = reference, prefix = prefix,
        joint = joint, symmetric = symmetric,
        strata = strata, params = params,
        prior = prior, rlang::enquos(...),
        min = rlang::enquo(min),
        max = rlang::enquo(max),
        proposal = proposal,
        proposal_scale = proposal_scale
      )
      private$rebuild_()
      invisible(self)
    },

    ############################################################################
    ############################################################################
    ############################################################################

    #' @description Add definition of infectivity profile.
    #' @param func `r rd_var_name_`
    #' @param status `r rd_aug_status_`
    #' @param discr `r rd_inf_discr_`
    #' @param ... `r rd_inf_args_`
    #' @param max `r rd_aug_max_`
    add_infec_profile = function(func, status = NA_integer_,
                                 discr = FALSE, ..., max = NULL) {
      private$infec_profile_$add_quo(func = func,
                                     status = status, discr = discr,
                                     rlang::enquos(...),
                                     max = rlang::enquo(max),
                                     .assert = TRUE)
      private$rebuild_()
      invisible(self)
    },

    #' @description Load infectivity profile definitions as named arguments.
    #' @param ... `r rd_list_args_`
    #' @param append `r rd_list_append_`
    load_infec_profile_list = function(..., append = FALSE) {
      private$infec_profile_$load_quosures(rlang::enquos(...), append = append)
      private$rebuild_()
      invisible(self)
    },

    ############################################################################
    ############################################################################
    ############################################################################

    #' @description Evaluate the log-likelihood.
    #' @param house_id `r rd_house_id_`
    #' If `NA` (default) it computes for the whole data.
    #' @return `numeric(1)`
    evaluate_log_lik = function(house_id = NA_character_) {
      private$build_full_()
      .engine$evaluate_log_lik(private$mcmc_ptr_, house_id)
    },

    #' @description Run inference with Bayesian Markov chain Monte Carlo.
    #' @param steps `r rd_smpl_steps_`
    #' @param burn `r rd_smpl_burn_`
    #' @param thin `r rd_smpl_thin_`
    #' @param adapt `r rd_smpl_adapt_`
    #' @param chains `r rd_smpl_chain_`
    #' @param threads `r rd_thread_`
    #' @param id_draws `r rd_smpl_iddrw_`
    #' @param progress `r rd_progress_`
    #' @return An object of class `HousepiPosteriorDraws`.
    run_inference = function(steps = 200, burn = 0,
                             thin = 1, adapt = 0,
                             chains = 1, threads = 1,
                             id_draws = NULL,
                             progress = TRUE) {
      private$build_full_(steps = steps, burn = burn,
                          thin = thin, adapt = adapt,
                          chains = chains, threads = threads,
                          id_draws = id_draws, progress = progress)
      .engine$run_inference(private$mcmc_ptr_)
      HousepiPosteriorDraws$new(private$name_, private$mcmc_ptr_)
    },

    #' @description
    #' Dataset with associated augmented variables for last state.
    #' @return `data.frame`
    export_augmented_data = function() {
      private$refactor_data_wide_(
        .engine$export_augmented_data_wide(private$mcmc_ptr_)
      )
    },

    #' @description
    #' All variables associated to given household for last state.
    #' @param house_id `r rd_house_id_`
    #' @return `list`
    export_house_snapshot = function(house_id) {
      house_id <- as.character(house_id)
      # export household snapshot
      output <- .engine$export_house_snapshot(private$mcmc_ptr_, house_id)
      output$data_wide <- private$refactor_data_wide_(output$data_wide)
      # > data long
      if (all(dim(output$data_long) == 0)) output$data_long <- NULL
      else output$data_long <- private$refactor_data_long_(output$data_long)
      # > group rates
      if (all(dim(output$group_rates) == 0)) output$group_rates <- NULL
      # > infectivity profile # UPDATE
      if (length(output$infec_profile) == 0) {
        output$infec_profile <- NULL
      } else {
        lst <- output$infec_profile
        df <- data.frame(id = lst$id, discr = lst$discr, func = lst$func,
                         min = lst$min, max = lst$max)
        df$func_args <- lst$func_args
        output$infec_profile <- df
      }
      output
    },

    #' @description
    #' Run simulation to generate a household survey.
    #' It assumes the same community structure of the loaded data.
    #' @param obs_process `r rd_simul_observ_`
    #' @param obs_process_long `r rd_simul_obs_long_`
    #' @param param_table `r rd_param_table_`
    #' @param nsim `r rd_simul_nsim_`
    #' @return `data.frame` or `list` (when `obs_process_long` is provided)
    run_simulation = function(obs_process = NULL,
                              obs_process_long = NULL,
                              param_table = NULL,
                              nsim = 1) {
      # simulation specific build
      nsim <- private$helper_simulation_build_(param_table = param_table,
                                               nsim = nsim)

      # run simulations
      output <- .engine$run_simulation(
        private$mcmc_ptr_,
        # combine observation process to measures and refactor
        private$helper_obs_process_(obs_process),
        NULL, # wide measure
        obs_process_long
      )

      # drop empty wide/long
      if (length(output$wide) == 0) output$wide <- NULL
      if (length(output$long) == 0) output$long <- NULL
      if (length(output) == 1) {
        output <- output[[1]]
        if (nsim == 1) output <- output[[1]]
      } else if (nsim == 1) {
        output$wide <- output$wide[[1]]
        output$long <- output$long[[1]]
      }

      output
    },

    ############################################################################
    ############################################################################
    ############################################################################

    #' @description Drop parameter definitions.
    drop_parameters = function() {
      private$param_ <- HousepiParameter$new()
      private$rebuild_()
      invisible(self)
    },

    #' @description Drop augmented variable definitions.
    drop_augmented = function() {
      private$augm_ <- HousepiAugmented$new()
      private$rebuild_()
      invisible(self)
    },

    #' @description Drop group rates definitions.
    drop_group_rates = function() {
      # remove associated parameters
      if (!is.null(private$group_rates_)) {
        for (p in names(private$group_rates_$parameters)) {
          private$param_$remove(p)
        }
        if (!is.null(private$param_$defs) && dim(private$param_$defs)[1] == 0) {
          self$drop_parameters()
        }
      }
      private$group_rates_ <- HousepiGroupRate$new()
      private$rebuild_()
      invisible(self)
    },

    #' @description Drop contact rate definitions.
    drop_contact_rates = function() {
      # remove associated parameters
      if (!is.null(private$contact_rates_)) {
        for (p in names(private$contact_rates_$parameters)) {
          private$param_$remove(p)
        }
        if (!is.null(private$param_$defs) && dim(private$param_$defs)[1] == 0) {
          self$drop_parameters()
        }
      }
      private$contact_rates_ <- HousepiContactRate$new()
      private$rebuild_()
      invisible(self)
    },

    #' @description Drop infectivity profile rates definitions.
    drop_infectivity_profile = function() {
      private$infec_profile_ <- HousepiInfecProfile$new()
      private$rebuild_()
      invisible(self)
    },

    #' @description Drop all definitions.
    drop_all = function() {
      self$drop_parameters()
      self$drop_augmented()
      self$drop_group_rates()
      self$drop_contact_rates()
      self$drop_infectivity_profile()
      invisible(self)
    }
  )
)

#' @title Export to Cpp Standalone
#' @description
#' Export definitions in standalone folder in order to run debuggin in Cpp.
#' @param hhmodel `r rd_model_in_`
#' @details
#' This could be designed to be more helpful, e. g.
#' by creating a standalone folder somewhere with correct links to headers
#' @keywords internal
#' @noRd
export_to_cpp_standalone <- function(hhmodel, outpath = "inst/standalone/") {
  stopifnothousepi(hhmodel)
  private <- hhmodel$.__enclos_env__$private
  private$build_group_rates_()
  private$build_contact_rates_()

  info_wide <- private$info_data_wide_
  utils::write.csv(private$augm_$defs,
                   file.path(outpath, "augm.csv"),
                   quote = FALSE, row.names = FALSE)
  utils::write.csv(private$contact_rates_$expose_to_cpp(info_wide),
                   file.path(outpath, "contact_rates.csv"),
                   quote = FALSE, row.names = FALSE)
  utils::write.csv(.engine$export_data_wide(private$mcmc_ptr_),
                   file.path(outpath, "data_wide.csv"),
                   quote = FALSE, row.names = FALSE)
  utils::write.csv(.engine$export_data_long(private$mcmc_ptr_),
                   file.path(outpath, "data_long.csv"),
                   quote = FALSE, row.names = FALSE)
  utils::write.csv(private$group_rates_$expose_to_cpp(info_wide),
                   file.path(outpath, "group_rates.csv"),
                   quote = FALSE, row.names = FALSE)
  utils::write.csv(private$infec_profile_$defs,
                   file.path(outpath, "infec.csv"),
                   quote = FALSE, row.names = FALSE)
  utils::write.csv(private$param_$defs,
                   file.path(outpath, "param.csv"),
                   quote = FALSE, row.names = FALSE)
  utils::write.csv(char_setup_input_(hhmodel$setup),
                   file.path(outpath, "setup.csv"),
                   quote = FALSE, row.names = FALSE)
}
