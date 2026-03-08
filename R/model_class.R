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
      c(required,                             # by the model
        private$augm_$data_required,          # by augmented
        private$infec_profile_$data_required) # by infec profile
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
      c(required,                              # by the model
        private$group_rates_$param_required,   # by group rates
        private$augm_$param_required,          # by augmented
        private$infec_profile_$param_required) # by infec profile
    },

    ############################################################################
    ############################################################################
    ############################################################################

    # @description
    # Assert if augmented definitions match the required.
    # @param capital `logical(1)` to raise error.
    assert_augm_ = function(capital = TRUE) {
      private$augm_$assert(private$get_all_required_augm_(),
                           private$infec_status_set_,
                           capital, !private$is_blank_())
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
                                  capital)
    },

    # @description
    # Assert if parameter definitions match required
    # @param capital `logical(1)` to raise error.
    assert_param_ = function(capital = TRUE) {
      private$param_$assert(private$get_all_required_param_(),
                            capital, !private$is_blank_())
    },

    ############################################################################
    ############################################################################
    ############################################################################

    # @description
    # Get template or previously implemented infectiviyt profile.
    get_internal_infec_profile_ = function() {
      .engine$export_infec_profile(private$mcmc_ptr_)
    },

    # @description
    # Generate multiplicity definitions
    # UPDATE: this will be moved inside C++ standalone in a future version.
    #
    # @return `data.frame`
    get_multiplicity_ = function() {
      # skip if augmented definitions are required
      if (length(private$get_all_required_augm_()) > 0) return(NULL)
      # skip if augmented definitions are provided
      # UPDATE: there could be the case of all constants...
      if (length(dim(private$augm_$defs)[1]) > 0) return(NULL)

      # load data_wide
      # UPDATE: this "back & forth" of data_wide is not optimal!
      data_wide <- self$data_wide

      headers <- c(
        # data headers required by implementation
        private$get_all_required_data_wide_(),
        # data headers required by group rates
        unlist(lapply(
          private$group_rates_$logic_tables,
          function(x) unlist(lapply(x, names))
        ))
      )
      # remove logic tables auxiliary column names
      headers <- setdiff(unique(headers), c("n", "param_name"))

      # list of splitted data by house_id
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
      mult_hhtable <- data.frame(class_id = sapply(hash_groups, `[`, 1),
                                 multiplicity = lengths(hash_groups))

      # use hashes to produce multiplicity table
      defs <- data.frame(house_id = unique(data_wide$house_id),
                         multiplicity = 0,
                         class_id = NA_integer_)
      rownames(defs) <- defs$house_id
      defs[mult_hhtable$class_id, "multiplicity"] <- mult_hhtable$multiplicity

      for (i in seq_len(nrow(mult_hhtable))) {
        defs[hash_groups[[i]], "class_id"] <- mult_hhtable$class_id[i]
      }
      defs
    },

    ############################################################################
    ############################################################################
    ############################################################################

    # @description
    # Pre-build to deal with group rates.
    build_group_rates_ = function() {
      # proceed only if definitions exist and data loaded
      if (length(private$group_rates_$defs) == 0) return(invisible(self))
      if (length(private$cols_data_wide_) == 0) {
        warning_h(
          "Missing data needed for building group rates.",
          "Ensure to provide all necessary data before proceeding."
        )
        return(invisible(self))
      }
      # load data_wide
      # UPDATE: this "back & forth" of data_wide is not optimal!
      # UPDATE: count in the wrapper while storing data_wide_levels
      data_wide <- self$data_wide

      # implement group rates
      private$group_rates_$implement(data_wide, private$info_data_wide_)

      # append parameter definitions if missing
      new_parameters <- private$group_rates_$parameters
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
    # Get individual for blank model augmented draws.
    # UPDATE: not robust if data has only default columns
    get_blank_id_draws_ = function(id_draws) {
      # if no default columns, assume debugging mode
      if (private$is_blank_() &&
            is.null(id_draws) &&
            !is.null(private$augm_$defs) &&
            length(setdiff(private$cols_data_wide_,
                           c("house_id", "id", "status"))) == 0) {
        id_draws <- template_blank_data_wide_()$id[1]
      }
      id_draws
    },

    # @description
    # Full model build with complete list of definitions.
    build_full_ = function(steps = 0, adapt = 0,
                           burnin = 0, thin = 1,
                           chains = 1, threads = 1,
                           id_draws = NULL,
                           seed = NA_integer_, progress = FALSE) {
      # import setup
      .engine$import_setup(private$mcmc_ptr_, char_setup_input_(list(
        steps = as.integer(steps),
        adapt = as.integer(adapt),
        burnin = as.integer(burnin),
        thin = as.integer(thin),
        chains = as.integer(chains),
        threads = as.integer(threads),
        id_draws = private$get_blank_id_draws_(id_draws),
        seed = private$adjust_seed_(seed),
        progress = as.integer(progress)
      )))

      # assert headers of last imported data_wide
      private$assert_data_wide_()

      # assert group rate definitions
      private$assert_group_rates_()
      # build group rate definitions
      private$build_group_rates_()

      # build contact rates
      # UPDATE: waiting for table definitions
      .engine$build_contact_rates(private$mcmc_ptr_)

      # assert parameter definitions
      private$assert_param_()
      # import parameters
      .engine$import_param(private$mcmc_ptr_,
                           private$param_$defs)

      # import group rates
      .engine$import_group_rate(
        private$mcmc_ptr_,
        private$group_rates_$expose_to_cpp(private$info_data_wide_)
      )

      # assert augmented variable definitions
      private$assert_augm_()

      # import infectivity profile
      .engine$import_infec_profile(private$mcmc_ptr_,
                                   private$infec_profile_$defs)

      # import augmented variables
      .engine$import_augm(private$mcmc_ptr_,
                          private$augm_$defs)

      # associate multiplicity engine
      .engine$import_multiplicity(private$mcmc_ptr_,
                                  private$get_multiplicity_())
      invisible(self)
    },

    # UPDATE:
    # this could be designed to be more helpful
    # e. g. create a standalone folder somewhere with correct links to headers
    export_to_standalone_ = function(outpath = "inst/standalone/") {
      private$build_group_rates_()

      write.csv(.engine$export_data_wide(private$mcmc_ptr_),
                file.path(outpath, "data_wide.csv"),
                quote = FALSE, row.names = FALSE)
      write.csv(private$get_multiplicity_(),
                file.path(outpath, "multiplicity.csv"),
                quote = FALSE, row.names = FALSE)
      write.csv(.engine$export_data_long(private$mcmc_ptr_),
                file.path(outpath, "data_long.csv"),
                quote = FALSE, row.names = FALSE)
      write.csv(private$group_rates_$expose_to_cpp(private$info_data_wide_),
                file.path(outpath, "group_rates.csv"),
                quote = FALSE, row.names = FALSE)
      write.csv(private$augm_$defs,
                file.path(outpath, "augm.csv"),
                quote = FALSE, row.names = FALSE)
      write.csv(private$param_$defs,
                file.path(outpath, "param.csv"),
                quote = FALSE, row.names = FALSE)
      write.csv(char_setup_input_(self$setup),
                file.path(outpath, "setup.csv"),
                quote = FALSE, row.names = FALSE)
      write.csv(private$infec_profile_$defs,
                file.path(outpath, "infec.csv"),
                quote = FALSE, row.names = FALSE)
      invisible(self)
    },

    # Helper with parameter table import.
    helper_import_table_ = function(param_table) {
      # extract parameter draws from table
      if (inherits(param_table, "HousepiPosteriorDraws")) {
        return(posterior::as_draws_matrix(param_table$param_draws))
      } else if (posterior::is_draws(param_table)) {
        return(posterior::as_draws_matrix(param_table))
      }
      tryCatch({
        return(as.data.frame(param_table))
      }, error = function(e) {
        stop_h("Type of param_table not supported.")
      })
    },

    # Helper with observation process import.
    # Combine observation process to measures and refactor.
    helper_obs_process_ = function(obs_process) {
      if (is.null(obs_process)) {
        private$refactor_data_wide_
      } else if (is.function(obs_process)) {
        function(x) obs_process(private$refactor_data_wide_(x))
      } else {
        stop_h("Unrecognized function passed as observation process.")
      }
    },

    # @description
    # Run measure function on posterior simulations.
    # UPDATE: this has to be improved with loops in C++
    run_measure_ = function(param_table,
                            apply_measure,
                            obs_process = NULL,
                            threads = 1,
                            seed = NA_integer_,
                            progress = TRUE) {
      standard_apply <- private$helper_obs_process_(obs_process)
      final_apply <- function(x) apply_measure(standard_apply(x))

      param_table <- private$helper_import_table_(param_table)
      nsim <- dim(param_table)[1]

      # internal build function
      private$build_full_(steps = nsim, threads = threads,
                          seed = seed, progress = progress)

      # attach param_table
      .engine$import_param_for_prediction(private$mcmc_ptr_, param_table)

      # run simulation and apply all measures
      .engine$run_simulation(private$mcmc_ptr_, final_apply)
    },

    # @description
    # Estimate the model logevidence (or marginal likelihood).
    # It uses the method of power posteriors with different available schemes.
    #
    # @param scheme `r rd_evi_scheme_`
    # @param rungs `r rd_evi_rungs_`
    # @param steps `r rd_smpl_steps_`
    # @param burnin `r rd_smpl_burn_`
    # @param thin `r rd_smpl_thin_`
    # @param adapt `r rd_smpl_adapt_`
    # @param chains `r rd_smpl_chain_`
    # @param threads `r rd_thread_`
    # @param seed `r rd_seed_`
    # @param progress `r rd_progress_`
    #
    # @return `numeric(1)`
    logevidence_ = function(scheme = "power-fractions", rungs = 40,
                            steps = 10000, burnin = 0,
                            thin = 1, adapt = 0,
                            chains = 1, threads = 1,
                            seed = NA_integer_, progress = TRUE) {
      warning_h(
        "[Internal] This method is currently under development.",
        "Its usage is deprecated."
      )

      # internal build function
      private$build_full_(steps = steps, burnin = burnin,
                          thin = thin, adapt = adapt,
                          chains = chains, threads = threads,
                          seed = seed, progress = progress)
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
      if (!missing(value)) stop_h("'data_wide' is read only.")
      private$refactor_data_wide_(.engine$export_data_wide(private$mcmc_ptr_))
    },

    #' @field data_long
    #' `data.frame`. Additional data rows indexed by individual (long format).
    data_long = function(value) {
      if (!missing(value)) stop_h("'data_long' is read only.")
      private$refactor_data_long_(.engine$export_data_long(private$mcmc_ptr_))
    },

    #' @field parameters `data.frame`. Parameter definitions.
    parameters = function(value) {
      if (!missing(value)) stop_h("'parameters' is read only.")
      private$build_group_rates_()
      housepi_defin_showcaser(private$param_$defs)
    },

    #' @field augmented `data.frame`. Augmented definitions.
    augmented = function(value) {
      if (!missing(value)) stop_h("'augmented' is read only.")
      housepi_defin_showcaser(private$augm_$defs)
    },

    #' @field infec_profile `data.frame`. Infectivity profile definitions.
    infec_profile = function(value) {
      if (!missing(value)) stop_h("'infec_profile' is read only.")
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
      if (!missing(value)) stop_h("'group_rates' is read only.")
      private$assert_group_rates_()
      private$build_group_rates_()
      private$group_rates_$tidy_table
    }
  ),

  ############################################################################
  ################################## PUBLIC ##################################
  ############################################################################

  public = list(
    #' @description
    #' Create a new object of class 'HousepiModel'
    #' used to run inference and simulation routines.
    #' @param choice `r rd_mdl_choice_`
    #' @param as_template `r rd_mdl_astmpl_`
    #' @param force `r rd_mdl_force_`
    #' @param quiet `r rd_mdl_quiet_`
    #' @param ... `r rd_mdl_flags_`
    initialize = function(choice = NA_character_, as_template = FALSE,
                          force = FALSE, quiet = FALSE, ...) {
      # helper construction option fort copying
      if (isTRUE(.housepi_env$hollow)) return(invisible(self))

      # initialize internal fields
      self$drop_all()
      # generate model wrapper
      super$initialize(choice = choice, as_template = as_template,
                       force = force, quiet = quiet, ...)
      # import blank template data_wide
      if (private$is_blank_()) self$add_data_wide(template_blank_data_wide_())

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
        if (!field %in% c("model_ptr_", "mcmc_ptr_") &&
              !is.function(val)) {
          obj_encl[[field]] <- helper_clone(val)
        }
      }

      # use cloned paths to create new Rcpp pointers
      obj_encl$create_new_pointers_()

      # UPDATE: this transfer shall be cleaned
      .engine$import_setup(
        obj_encl$mcmc_ptr_,
        char_setup_input_(.engine$export_setup(private$mcmc_ptr_))
      )
      # UPDATE: avoid continuos re-importing data
      .engine$import_data_wide(
        obj_encl$mcmc_ptr_,
        char_data_input_(.engine$export_data_wide(private$mcmc_ptr_))
      )
      .engine$import_data_long(
        obj_encl$mcmc_ptr_,
        char_data_input_(.engine$export_data_long(private$mcmc_ptr_))
      )

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
        fmt_assert_cli("Data", private$assert_data_wide_(FALSE),
                       minor = length(private$cols_data_wide_) == 0)
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

    #' @description Store a copy of provided household data.
    #' @param data_wide `r rd_add_data_wide_`
    add_data_wide = function(data_wide) {
      if (length(data_wide) == 0) return(invisible(self))

      # check that data_wide type is correct
      if (!is.data.frame(data_wide)) stop_h(
        "Unsupported type for data_wide.",
        "Dataset should be a data.frame."
      )

      # convert factors
      embed <- cpp_data_factorize(data_wide)
      # store levels
      private$info_data_wide_ <- embed$info
      # store column names
      private$cols_data_wide_ <- embed$cols
      # store all positive status
      private$infec_status_set_ <- embed$status

      # assert data_wide headers
      tryCatch({
        private$assert_data_wide_()
      }, error = function(e) {
        stop_h(e)
      })

      # import data.frame in rcpp
      .engine$import_data_wide(private$mcmc_ptr_, embed$data)

      invisible(self)
    },

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
        "Unsupported type for data_long.",
        "Datalist should be a data.frame."
      )

      # convert factors
      embed <- cpp_data_factorize(data_long)
      # store levels
      private$info_data_long_ <- embed$info
      # store column names
      private$cols_data_long_ <- embed$cols
      # import data.frame in rcpp
      .engine$import_data_long(private$mcmc_ptr_, embed$data)

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
      invisible(self)
    },

    #' @description
    #' Load parameter definitions as named arguments.
    #' @param ... `r rd_list_args_`
    #' @param append `r rd_list_append_`
    load_parameter_list = function(..., append = FALSE) {
      private$param_$load_quosures(rlang::enquos(...), append = append)
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
      invisible(self)
    },

    #' @description Load augmented definitions as named arguments.
    #' @param ... `r rd_list_args_`
    #' @param append `r rd_list_append_`
    load_augmented_list = function(..., append = FALSE) {
      private$augm_$load_quosures(rlang::enquos(...), append = append)
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
    #' @param independent `r rd_grr_indep_`
    #' @param prior `r rd_par_model_`
    #' @param ... `r rd_par_args_`
    #' @param min `r rd_var_min_`
    #' @param max `r rd_var_max_`
    #' @param proposal `r rd_var_proposal_`
    #' @param proposal_scale `r rd_proposal_scale_`
    add_group_rate = function(name, by, prefix = NA_character_,
                              reference = NULL, independent = TRUE,
                              prior = NA_character_, ...,
                              min = NULL, max = NULL,
                              proposal = NA_character_,
                              proposal_scale = NA_real_) {
      private$group_rates_$add_quo(name = name, by = by, prefix = prefix,
                                   reference = reference,
                                   independent = independent,
                                   prior = prior, rlang::enquos(...),
                                   min = rlang::enquo(min),
                                   max = rlang::enquo(max),
                                   proposal = proposal,
                                   proposal_scale = proposal_scale)
      # UPDATE: I didn't find a way to assert prior arguments here
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
                                 discr = FALSE,
                                 ..., max = NULL) {
      private$infec_profile_$add_quo(func = func,
                                     status = status, discr = discr,
                                     rlang::enquos(...),
                                     max = rlang::enquo(max),
                                     .assert = TRUE)
      invisible(self)
    },

    #' @description Load infectivity profile definitions as named arguments.
    #' @param ... `r rd_list_args_`
    #' @param append `r rd_list_append_`
    load_infec_profile_list = function(..., append = FALSE) {
      private$infec_profile_$load_quosures(rlang::enquos(...), append = append)
      invisible(self)
    },

    ############################################################################
    ############################################################################
    ############################################################################

    #' @description Evaluate the log-likelihood.
    #' @param house_id `r rd_house_id_`
    #' If `NA` (default) it computes for the whole data.
    #' @param seed `r rd_seed_`
    #' @return `numeric(1)`
    evaluate_log_lik = function(house_id = NA_character_, seed = NA_integer_) {
      private$build_full_(seed = seed)
      .engine$evaluate_log_lik(private$mcmc_ptr_, house_id)
    },

    #' @description Run inference with Bayesian Markov chain Monte Carlo.
    #' @param steps `r rd_smpl_steps_`
    #' @param burnin `r rd_smpl_burn_`
    #' @param thin `r rd_smpl_thin_`
    #' @param adapt `r rd_smpl_adapt_`
    #' @param chains `r rd_smpl_chain_`
    #' @param threads `r rd_thread_`
    #' @param id_draws `r rd_smpl_iddrw_`
    #' @param seed `r rd_seed_`
    #' @param progress `r rd_progress_`
    #' @return An object of class `HousepiPosteriorDraws`.
    run_inference = function(steps = 200, burnin = 0,
                             thin = 1, adapt = 0,
                             chains = 1, threads = 1,
                             id_draws = NULL,
                             seed = NA_integer_,
                             progress = TRUE) {
      private$build_full_(steps = steps, burnin = burnin,
                          thin = thin, adapt = adapt,
                          chains = chains, threads = threads,
                          id_draws = id_draws, seed = seed, progress = progress)
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
    #' @param nsim `r rd_simul_nsim_`
    #' @param seed `r rd_seed_`
    #' @return `data.frame`
    run_simulation = function(obs_process = NULL,
                              nsim = 1,
                              seed = NA_integer_) {
      # combine observation process to measures and refactor
      standard_apply <- private$helper_obs_process_(obs_process)

      # internal build function
      private$build_full_(steps = nsim, seed = seed, progress = FALSE)

      # attach param_table
      if (dim(self$parameters)[1] > 0) {
        param_table <- self$parameters |>
          dplyr::select(name, value) |>
          tidyr::pivot_wider(names_from = "name", values_from = "value") |>
          dplyr::mutate(n = nsim) |>
          tidyr::uncount(n) |>
          private$helper_import_table_()
        .engine$import_param_for_prediction(private$mcmc_ptr_, param_table)
      }

      # run simulations
      output <- .engine$run_simulation(private$mcmc_ptr_, standard_apply)
      # return a single data.frame if single simulation
      if (nsim == 1) output <- output[[1]]

      output
    },

    ############################################################################
    ############################################################################
    ############################################################################

    #' @description Drop parameter definitions.
    drop_parameters = function() {
      private$param_ <- HousepiParameter$new()
      invisible(self)
    },

    #' @description Drop augmented variable definitions.
    drop_augmented = function() {
      private$augm_ <- HousepiAugmented$new()
      invisible(self)
    },

    #' @description Drop group rates definitions.
    drop_group_rates = function() {
      # remove associated parameters
      if (!is.null(private$group_rates_)) {
        for (p in names(private$group_rates_$parameters)) {
          private$param_$remove(p)
        }
      }
      private$group_rates_ <- HousepiGroupRate$new()
      invisible(self)
    },

    #' @description Drop infectivity profile rates definitions.
    drop_infectivity_profile = function() {
      private$infec_profile_ <- HousepiInfecProfile$new()
      invisible(self)
    },

    #' @description Drop all definitions.
    drop_all = function() {
      self$drop_parameters()
      self$drop_augmented()
      self$drop_group_rates()
      self$drop_infectivity_profile()
      invisible(self)
    }
  )
)
