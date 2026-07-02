## housepi - Epidemiological Analysis of Household Outbreaks.
##
## Copyright © 2024-2026 Institut Pasteur (Paris), CNRS and INSERM.
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
#' R6 class for compiling and running a household transmission model
#' @description
#' Compile the household transmission model, manage your data, associate
#' definitions, run inference and simulations from a single class.
#' @docType class
HousepiModel <- R6::R6Class(
  classname = "HousepiModel",
  cloneable = FALSE,
  lock_class = TRUE,

  ############################################################################
  ################################## PRIVATE #################################
  ############################################################################

  private = list(
    generator_ = NULL,
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
      refer_file <- private$generator_$get_reference()
      c(private$augm_$require(refer_file),    # by the model
        private$infec_profile_$augm_required) # by infec profile
    },

    # @description
    # Dataset headers required by implementation file and other definitions.
    get_all_required_data_wide_ = function() {
      c(private$generator_$data_required,     # by the model
        private$group_rates_$data_required,   # by group rates
        private$contact_rates_$data_required, # by contact rates
        private$augm_$data_required,          # by augmented
        private$infec_profile_$data_required) # by infec profile
    },

    # @description
    # Group rates required by implementation file and other definitions.
    get_all_required_group_rates_ = function() {
      refer_file <- private$generator_$get_reference()
      c(private$group_rates_$require(refer_file),    # by the model
        private$augm_$group_rates_required,          # by augmented
        private$infec_profile_$group_rates_required) # by infec profile
    },

    # @description
    # Parameters required by implementation file and other definitions.
    get_all_required_param_ = function() {
      refer_file <- private$generator_$get_reference()
      c(private$param_$require(refer_file),    # by the model
        private$param_$param_required,         # by parameter functions
        private$group_rates_$param_required,   # by group rates
        private$contact_rates_$param_required, # by contact rates
        private$augm_$param_required,          # by augmented
        private$infec_profile_$param_required) # by infec profile
    },

    ############################################################################
    ############################################################################
    ############################################################################

    # @description
    # Assert if augmented definitions match the required.
    # @param capital `logical(1)` to raise error.
    assert_augm_ = function(capital = TRUE, warn_unused = TRUE) {
      private$augm_$assert(
        private$get_all_required_augm_(),
        private$generator_$infec_status_set,
        capital = capital,
        warn_unused = warn_unused && !private$generator_$is_blank()
      )
    },

    # @description
    # Assert if data_wide headers match implementation file.
    # @param capital `logical(1)` to raise error.
    assert_data_wide_ = function(capital = TRUE, simulating = FALSE) {
      cols <- private$generator_$cols_data_wide
      required <- private$get_all_required_data_wide_()

      # observations are generated during simulation: not required
      if (isTRUE(simulating)) {
        required <- setdiff(required, private$augm_$observation_names)
      }

      # household observations derived purely from data satisfy requirement
      house_obs <- private$augm_$house_observ_names
      clash <- intersect(house_obs, cols)
      if (isTRUE(capital) && length(clash) > 0) stop_h(
        "Household observation {.field {clash}} cannot stand in for required data.",
        "It is already provided as a data column."
      )
      cols <- c(cols, intersect(house_obs, required))

      sms <- template_missing_stop(required, cols, "data")
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
      if (has_defs &&
            !.engine$has_contact_rates(private$generator_$mcmc_ptr)) {
        private$drop_contact_rates_()
        stop_h(
          "Contact rates are not supported for model {.val {self$name}}.",
          `!` = "Definitions will be dropped."
        )
      }
    },

    # @description
    # Assert if parameter definitions match required.
    # @param capital `logical(1)` to raise error.
    # @param warn_unused `logical(1)` warn if unused.
    assert_param_ = function(capital = TRUE, warn_unused = TRUE) {
      private$param_$assert(
        private$get_all_required_param_(),
        capital = capital,
        warn_unused = warn_unused && !private$generator_$is_blank()
      )
    },

    ############################################################################
    ############################################################################
    ############################################################################

    # @description
    # Get template or previously implemented infectiviyt profile.
    get_internal_infec_profile_ = function() {
      .engine$export_infec_profile(private$generator_$mcmc_ptr)
    },

    ############################################################################
    ############################################################################
    ############################################################################

    # @description
    # Helper to assert and re-build group and contact rates.
    helper_build_rates_ = function(module) {
      # proceed only if definitions exist
      if (length(module$defs) == 0) return(invisible(self))

      # UPDATE: this "back & forth" of data_wide is sub-optimal...
      module$implement(self$data_wide, private$generator_$info_data_wide)

      # append parameter if missing
      new_parameters <- module$parameters
      # find missing parameters
      mask <- setdiff(names(new_parameters), private$param_$defs$name)
      new_parameters <- new_parameters[mask]
      if (length(new_parameters) > 0) {
        xquos <- convert_list_to_quosures(new_parameters)
        private$param_$add_quo_list(xquos, append = TRUE)
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
    # Reset every definition holder to an empty object.
    reset_definitions_ = function() {
      private$param_ <- HousepiParameter$new()
      private$augm_ <- HousepiAugmented$new()
      private$group_rates_ <- HousepiGroupRate$new()
      private$contact_rates_ <- HousepiContactRate$new()
      private$infec_profile_ <- HousepiInfecProfile$new()
      invisible(self)
    },

    # @description
    # Drop parameter definitions.
    drop_parameters_ = function() {
      private$param_ <- HousepiParameter$new()
      private$rebuild_()
      invisible(self)
    },

    # @description
    # Drop group rates definitions, removing any associated parameters.
    drop_group_rates_ = function() {
      if (!is.null(private$group_rates_)) {
        for (p in names(private$group_rates_$parameters)) {
          private$param_$remove(p)
        }
        if (!is.null(private$param_$defs) && nrow(private$param_$defs) == 0) {
          private$drop_parameters_()
        }
      }
      private$group_rates_ <- HousepiGroupRate$new()
      private$rebuild_()
      invisible(self)
    },

    # @description
    # Drop contact rate definitions, removing any associated parameters.
    drop_contact_rates_ = function() {
      if (!is.null(private$contact_rates_)) {
        for (p in names(private$contact_rates_$parameters)) {
          private$param_$remove(p)
        }
        if (!is.null(private$param_$defs) && nrow(private$param_$defs) == 0) {
          private$drop_parameters_()
        }
      }
      private$contact_rates_ <- HousepiContactRate$new()
      private$rebuild_()
      invisible(self)
    },

    # @description
    # Single, model-owned ordered infected-status dictionary. Data labels win
    # (the data is already coded), else the status process order, else literal
    # character labels declared on augmented definitions.
    status_dictionary_ = function() {
      data_lab <- private$generator_$status_data_labels
      if (length(data_lab) > 0) return(data_lab)
      private$augm_$status_dictionary
    },

    # @description
    # Soft re-push of R-side definitions to the C++ engine.
    # NOTE: Hard validation lives in `build_full_`.
    rebuild_ = function() {
      # guard: pre-init (before generator is created)
      if (is.null(private$generator_) ||
            is.null(private$generator_$mcmc_ptr)) {
        return(invisible(self))
      }
      # guard: data_wide not yet imported (during construction)
      if (length(private$generator_$cols_data_wide) == 0) {
        return(invisible(self))
      }

      # soft data_wide check (headers vs. implementation file)
      private$assert_data_wide_(FALSE)

      # build rates first (they append to param_$defs)
      if (length(private$group_rates_$defs) > 0) {
        private$build_group_rates_()
      }
      if (length(private$contact_rates_$defs) > 0) {
        private$build_contact_rates_()
      }

      mcmc_ptr <- private$generator_$mcmc_ptr
      info_wide <- private$generator_$info_data_wide

      #> parameters
      if (length(private$param_$defs) > 0) private$assert_param_(FALSE, FALSE)
      .engine$import_param(mcmc_ptr, private$param_$defs)

      #> group rates
      .engine$import_group_rate(
        mcmc_ptr,
        private$group_rates_$expose_to_cpp(info_wide)
      )

      #> contact rates
      .engine$import_contact_rate(
        mcmc_ptr,
        private$contact_rates_$expose_to_cpp(info_wide)
      )

      #> infectivity profile
      .engine$import_infec_profile(mcmc_ptr, private$infec_profile_$defs)

      #> augmented (resolve character status labels to integer codes for C++)
      status_dict <- private$status_dictionary_()
      if (private$generator_$escaped %in% status_dict) stop_h(
        paste("Infected status {.val {private$generator_$escaped}} collides",
              "with the escaped-status label ({.val {private$generator_$escaped}})."),
        "Rename the status or define a different {.code escaped =} label."
      )
      if (length(private$augm_$defs) > 0) private$assert_augm_(FALSE, FALSE)
      .engine$import_augm(
        mcmc_ptr,
        resolve_status_codes_(private$augm_$defs, status_dict)
      )

      # status dictionary drives the factor refactoring of the status column
      private$generator_$set_status_levels(status_dict)

      invisible(self)
    },

    # @description
    # Full model build. Refresh runtime setup, enforce hard validation,
    # then push the final state to the engine.
    build_full_ = function(steps = 0L, adapt = 0L, burn = 0L,
                           thin = 1L, chains = 1L, threads = 1L,
                           id_draws = NULL, progress = FALSE,
                           simulating = FALSE, ...) {
      # if no default columns on a blank model, assume debugging mode
      if (is.null(id_draws) &&
            private$generator_$is_blank() &&
            !is.null(private$augm_$defs) &&
            length(setdiff(private$generator_$cols_data_wide,
                           c("house_id", "id", "status"))) == 0) {
        env <- new.env(parent = emptyenv())
        utils::data("blank_house", package = "housepi", envir = env)
        id_draws <- env$blank_house$id[1]
      }
      .engine$import_setup(
        private$generator_$mcmc_ptr,
        char_setup_input_(c(list(
          steps = steps,
          adapt = adapt,
          burn = burn,
          thin = thin,
          chains = chains,
          threads = threads,
          id_draws = id_draws,
          # generate a seed for Cpp random number generator
          seed = sample.int(.Machine$integer.max, 1L),
          progress = progress
        ), list(...)))
      )
      # hard validation before MCMC
      private$assert_data_wide_(capital = TRUE, simulating = simulating)
      private$assert_param_(capital = TRUE)
      private$assert_group_rates_(capital = TRUE)
      private$assert_contact_rates_()
      private$assert_augm_(capital = TRUE)
      # push final state (picks up any stage deferred by dep gating)
      private$rebuild_()
      invisible(self)
    },

    # Assert peak memory needed to materialise draws is available.
    assert_posterior_memory_ = function(steps, burn, thin, chains, id_draws) {
      n_sample <- setup_sample_size_(steps = steps, burn = burn,
                                     thin = thin, chains = chains)$total

      # UPDATE: n. of augmented/parameters might be smaller if functions/private
      n_params_defs <- ifelse(is.null(private$param_$defs), 0L,
                              nrow(private$param_$defs))
      n_param_cols <- 2L * n_params_defs
      n_augm_defs <- ifelse(is.null(private$augm_$defs), 0L,
                            nrow(private$augm_$defs))
      n_augm_cols <- length(unique(id_draws)) * n_augm_defs
      n_house_total <- length(unique(self$data_wide$house_id))
      # total count of columns
      total_cols <- n_param_cols + n_augm_cols + n_house_total + 1L
      bytes_final <- 8L * n_sample * total_cols

      # byte budget for posterior draws materialisation
      # set 80% cap under R vsize, otherwise fall back to 16 GiB
      budget <- 16 * 2^30
      vsize_mb <- tryCatch(mem.maxVSize(), error = function(e) Inf)
      if (is.finite(vsize_mb) && vsize_mb > 0) budget <- 0.8 * vsize_mb * 2^20

      if (bytes_final <= budget) return(invisible())
      stop_h(c(
        sprintf(
          "Posterior draws would exceed the memory budget (%.f GB).",
          budget / 2^30
        ),
        "Consider increasing {.arg thin}, or reducing {.arg {c('steps', 'chains')}}."
      ))
    },

    # Helper with parameter table import.
    # @param ptab Table with parameter values
    # param thin_ess `logical(1)`
    helper_read_param_table_ = function(ptab, thin_ess = TRUE) {
      # extract parameter draws from table
      if (inherits(ptab, "HousepiPosteriorDraws")) ptab <- ptab$parameter_draws

      if (posterior::is_draws(ptab)) {
        ptab <- posterior::as_draws_matrix(ptab)
        # reduce sample to effective sample size (ESS)
        if (isTRUE(thin_ess)) {
          ess_bulk <- sapply(posterior::variables(ptab),
                             function(v) posterior::ess_bulk(ptab[, v]))
          # check no NA values in ESS
          if (any(is.na(ess_bulk))) stop_h(
            "ESS values return NA.",
            "Check your draws for degenerate parameters."
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
      }, error = function(e) stop_h("Type of {.arg param_table} not supported.")
      )
    },

    # Helper with simulation build.
    # FIXME: If there's a param table, build should change
    helper_simulation_build_ = function(param_table = NULL,
                                        nrep = NA_integer_,
                                        thin_ess = TRUE,
                                        ...) {
      # assert number of replicates
      if (!is.na(nrep) && (!is.numeric(nrep) || nrep < 1)) stop_h(
        "Number of replicates should be a positive integer."
      )

      # read param table
      if (is.null(param_table)) {
        if (is.na(nrep)) nrep <- 1
      } else {
        param_table <- private$helper_read_param_table_(
          param_table, thin_ess = thin_ess
        )
        if (!is.na(nrep)) param_table <- helper_thin_table_(param_table, nrep)
        nrep <- nrow(param_table)
        # don't overwrite
        if (nrow(self$parameters) == 0) {
          param_dummy_list <- as.list(utils::head(param_table, 1))
          private$param_$add_quo_list(
            convert_list_to_quosures(param_dummy_list), append = TRUE
          )
          private$rebuild_()
        }
      }

      # internal build function
      private$build_full_(
        steps = nrep, progress = FALSE, simulating = TRUE, ...
      )

      # NOTE: there might be not parameters at all
      if (nrow(self$parameters) == 0) return(nrep)

      if (is.null(param_table)) {
        # replicates the single row nrep times
        param_table <- as.data.frame(stats::setNames(
          as.list(self$parameters$value), self$parameters$name
        ))
        param_table <- param_table[rep(1L, nrep), , drop = FALSE]
        rownames(param_table) <- NULL
      }
      # attach table with parameters
      .engine$import_param_table(private$generator_$mcmc_ptr, param_table)

      nrep
    },

    # Helper with observation process import.
    # Combine observation process to measures and refactor.
    # UPDATE: I shall drop necessity for refactoring with variant data.
    helper_obs_process_ = function(obs_process) {
      gen <- private$generator_
      if (is.null(obs_process)) {
        function(x) gen$refactor_data_wide(x)
      } else if (is.function(obs_process)) {
        function(x) obs_process(gen$refactor_data_wide(x))
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
                            thin_ess = TRUE,
                            max_nrep = NA_integer_,
                            threads = 1L,
                            progress = TRUE) {
      param_table <- private$helper_read_param_table_(param_table,
                                                      thin_ess = thin_ess)
      param_table <- helper_thin_table_(param_table, max_nrep)
      nsim <- nrow(param_table)

      # internal build function
      private$build_full_(steps = nsim, threads = threads, progress = progress)

      # attach param_table
      .engine$import_param_table(private$generator_$mcmc_ptr, param_table)

      if (!is.null(apply_long) &&
            (is.null(obs_proc_long) || !is.function(obs_proc_long))) stop_h(
        "{.arg obs_process_long} must be a function."
      )

      obs_func <- private$helper_obs_process_(obs_proc_wide)
      .engine$run_simulation(
        private$generator_$mcmc_ptr,
        obs_func, apply_wide,
        if (is.null(apply_long)) NULL
        else function(x) apply_long(obs_proc_long(x))
      )
    }
  ),

  ############################################################################
  ################################## ACTIVE ##################################
  ############################################################################

  active = list(
    #' @field name `character(1)`. Model name.
    name = function(value) {
      if (!missing(value)) stop_read_only("name")
      private$generator_$name
    },

    #' @field setup `list`. Last setup used by the build process.
    setup = function(value) {
      if (!missing(value)) stop_read_only("setup")
      private$generator_$setup
    },

    #' @field data_wide `data.frame`. Individual level data.
    #' One row per individual index by corresponding ID.
    data_wide = function(value) {
      if (!missing(value)) stop_read_only("data_wide")
      private$generator_$refactor_data_wide(
        .engine$export_data_wide(private$generator_$mcmc_ptr)
      )
    },

    #' @field data_long
    #' `data.frame`. Additional data rows indexed by individual (long format).
    data_long = function(value) {
      if (!missing(value)) stop_read_only("data_long")
      private$generator_$refactor_data_long(
        .engine$export_data_long(private$generator_$mcmc_ptr)
      )
    },

    #' @field parameters `data.frame`. Parameter definitions.
    parameters = function(value) {
      if (!missing(value)) stop_read_only("parameters")
      housepi_defin_showcaser(private$param_$defs)
    },

    #' @field augmented `data.frame`. Augmented definitions.
    augmented = function(value) {
      if (!missing(value)) stop_read_only("augmented")
      defs <- private$augm_$defs_of_kind("augm")
      housepi_defin_showcaser(defs)
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
    #' @param escaped `r rd_mdl_escaped_`
    #' @param ... `r rd_mdl_flags_`
    initialize = function(data = NULL, model = NA_character_,
                          as_template = FALSE,
                          force = FALSE, quiet = FALSE,
                          escaped = "escaped", ...) {
      # helper construction option fort copying
      if (isTRUE(.housepi_env$hollow)) return(invisible(self))

      # check data and enable blank_house behavior
      data <- assert_household_data(data, model)

      # initialize internal fields
      private$reset_definitions_()
      # generate model wrapper
      private$generator_ <- HousepiGenerator$new(
        model = model, as_template = as_template,
        force = force, quiet = quiet, escaped = escaped, ...
      )
      # factorize data in order to pass it to Cpp
      data_factorized <- private$generator_$factorize_data_wide(data)

      # soft check: data may be intentionally incomplete at construction
      # NOTE: hard check lives at build_full_
      private$assert_data_wide_(capital = FALSE)
      .engine$import_data_wide(
        private$generator_$mcmc_ptr, data_factorized
      )

      invisible(self)
    },

    #' @description Deep copy this object.
    #' @return A new independent object.
    copy = function() {
      .housepi_env$hollow <- TRUE
      on.exit(.housepi_env$hollow <- FALSE)
      obj <- HousepiModel$new()
      obj_priv <- obj$.__enclos_env__$private

      # clone definition fields
      for (field in c("param_", "augm_", "group_rates_",
                      "contact_rates_", "infec_profile_")) {
        obj_priv[[field]] <- helper_clone(private[[field]])
      }

      # clone the generator (private fields except mcmc_ptr_)
      # then, re-init pointers
      new_gen <- HousepiGenerator$new()
      new_gen_priv <- new_gen$.__enclos_env__$private
      src_gen_priv <- private$generator_$.__enclos_env__$private
      for (field in names(new_gen_priv)) {
        val <- src_gen_priv[[field]]
        if (field != "mcmc_ptr_" && !is.function(val)) {
          new_gen_priv[[field]] <- helper_clone(val)
        }
      }
      new_gen$init_pointers()
      obj_priv$generator_ <- new_gen

      # NOTE: setup is run-scoped, fully resupplied by build_full_ at every
      # run/evaluate, so it is intentionally not carried across a copy.
      # share data store (just copy-on-write)
      .engine$share_data_store(new_gen$mcmc_ptr, private$generator_$mcmc_ptr)

      # fresh pointer needs definitions pushed to the engine
      obj_priv$rebuild_()

      obj
    },

    #' @description
    #' Override the default printing behavior.
    print = function() {
      is_tmpl <- !is.na(get_template_dir(private$generator_$raw_name))
      tmpl_tag <- if (is_tmpl) sprintf(
        " ({.run [template](housepi::browse_housepi_templates(\"%s\"))})",
        self$name
      ) else ""
      sms <- sprintf("Household model {.val %s}%s.", self$name, tmpl_tag)
      # UPDATE: to be rationalize in the future if more options are included
      if ("_HOUSEPI_DISCRETE_TIME" %in% private$generator_$last_flags) sms <- c(
        sms, "Model is set to {.emph discrete time} mode."
      )
      # data
      if (length(private$get_all_required_data_wide_()) > 0) sms <- c(sms,
        fmt_assert_cli("Data", private$assert_data_wide_(FALSE))
      )
      # parameters
      if (length(private$get_all_required_param_()) > 0) sms <- c(sms,
        fmt_assert_cli("Parameters", private$assert_param_(FALSE))
      )
      # augmented
      if (length(private$get_all_required_augm_()) > 0) sms <- c(sms,
        fmt_assert_cli("Augmented", private$assert_augm_(FALSE))
      )

      # inform about rates if dataset is suitable
      cols_extra <- setdiff(private$generator_$cols_data_wide,
                            data_reserved_cols)
      if (length(cols_extra) > 0) {
        # group rates
        if (length(private$get_all_required_group_rates_()) > 0) sms <- c(sms,
          fmt_assert_cli("Group rates", private$assert_group_rates_(FALSE),
                         minor = TRUE)
        )
        # contact rates
        if (.engine$has_contact_rates(private$generator_$mcmc_ptr)) {
          sms <- c(sms, fmt_assert_cli(
            "Contact rates",
            if (!private$contact_rates_$assert()) c(
              "Using default contact rates.",
              "A custom choice can be provided."
            )
            else NULL,
            minor = TRUE
          ))
        }
      }

      # infectivity profile
      if (length(private$get_internal_infec_profile_()) > 0) sms <- c(sms,
        fmt_assert_cli(
          "Infectivity profile",
          if (length(private$infec_profile_$showcase()) == 0) c(
            "Using templated profile.",
            "A custom choice can be provided."
          )
          else NULL,
          minor = TRUE
        )
      )

      fmt_cli_info(sms)
      invisible(self)
    },

    #' @description
    #' Show model description in RStudio Help window.
    help = function() {
      open_model_helper_(private$generator_$source_dir)
      invisible(self)
    },

    ############################################################################
    ############################################################################
    ############################################################################

    #' @description
    #' Dataset with associated augmented variables for last state.
    #' @return `data.frame`
    export_augmented_data = function() {
      private$generator_$refactor_data_wide(
        .engine$export_augmented_data_wide(private$generator_$mcmc_ptr)
      )
    },

    #' @description
    #' All variables associated to given household for last state.
    #' @param house_id `r rd_house_id_`
    #' @return `list`
    export_house_snapshot = function(house_id) {
      house_id <- as.character(house_id)
      # export household snapshot
      out <- .engine$export_house_snapshot(private$generator_$mcmc_ptr, house_id)
      out$data_wide <- private$generator_$refactor_data_wide(out$data_wide)
      # > data long
      if (all(dim(out$data_long) == 0)) out$data_long <- NULL
      else out$data_long <- private$generator_$refactor_data_long(out$data_long)
      # > group rates
      if (all(dim(out$group_rates) == 0)) out$group_rates <- NULL
      # > infectivity profile # UPDATE
      if (length(out$infec_profile) == 0) {
        out$infec_profile <- NULL
      } else {
        lst <- out$infec_profile
        df <- data.frame(id = lst$id, discr = lst$discr, func = lst$func,
                         min = lst$min, max = lst$max)
        df$func_args <- lst$func_args
        out$infec_profile <- df
      }
      out
    }
  )
)

# nocov start
#' @title Export to Cpp standalone (developer)
#' @description
#' Export definitions in standalone folder in order to run debuggin in Cpp.
#' Used for development.
#' @param hhmodel `r rd_model_in_`
#' @details
#' This could be designed to be more helpful, e. g.
#' by creating a standalone folder somewhere with correct links to headers
#' @keywords internal
#' @noRd
export_to_cpp_standalone <- function(hhmodel, outpath = "inst/standalone/") {
  stopifnothousepi(hhmodel)
  pvt <- hhmodel$.__enclos_env__$private
  pvt$build_group_rates_()
  pvt$build_contact_rates_()

  gen <- pvt$generator_
  info_wide <- gen$info_data_wide
  utils::write.csv(pvt$augm_$defs,
                   file.path(outpath, "augm.csv"),
                   quote = FALSE, row.names = FALSE)
  utils::write.csv(pvt$contact_rates_$expose_to_cpp(info_wide),
                   file.path(outpath, "contact_rates.csv"),
                   quote = FALSE, row.names = FALSE)
  utils::write.csv(.engine$export_data_wide(gen$mcmc_ptr),
                   file.path(outpath, "data_wide.csv"),
                   quote = FALSE, row.names = FALSE)
  utils::write.csv(.engine$export_data_long(gen$mcmc_ptr),
                   file.path(outpath, "data_long.csv"),
                   quote = FALSE, row.names = FALSE)
  utils::write.csv(pvt$group_rates_$expose_to_cpp(info_wide),
                   file.path(outpath, "group_rates.csv"),
                   quote = FALSE, row.names = FALSE)
  utils::write.csv(pvt$infec_profile_$defs,
                   file.path(outpath, "infec.csv"),
                   quote = FALSE, row.names = FALSE)
  utils::write.csv(pvt$param_$defs,
                   file.path(outpath, "param.csv"),
                   quote = FALSE, row.names = FALSE)
  utils::write.csv(char_setup_input_(hhmodel$setup),
                   file.path(outpath, "setup.csv"),
                   quote = FALSE, row.names = FALSE)
}
# nocov end
