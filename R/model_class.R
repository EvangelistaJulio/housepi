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
    group_risks_ = NULL,
    contact_risks_ = NULL,
    infec_profile_ = NULL,

    ############################################################################
    ############################################################################
    ############################################################################

    # @description
    # Latent variables required by implementation file and other definitions.
    get_all_required_latent_ = function() {
      c(private$generator_$latent_required,     # by the model
        private$augm_$latent_required,          # by latent/observation
        private$infec_profile_$latent_required) # by infectivity profile
    },

    # @description
    # Dataset headers required by implementation file and other definitions.
    get_all_required_data_wide_ = function() {
      c(private$generator_$data_required,     # by the model
        private$group_risks_$data_required,   # by group risks
        private$contact_risks_$data_required, # by contact risks
        private$augm_$data_required,          # by latent/observation
        private$infec_profile_$data_required) # by infectivity profile
    },

    # @description
    # Group risks required by implementation file and other definitions.
    get_all_required_group_risks_ = function() {
      c(private$generator_$group_risks_required,     # by the model
        private$augm_$group_risks_required,          # by latent/observation
        private$infec_profile_$group_risks_required) # by infectivity profile
    },

    # @description
    # Parameters required by implementation file and other definitions.
    get_all_required_param_ = function() {
      c(private$generator_$param_required,     # by the model
        private$param_$param_required,         # by parameter functions
        private$group_risks_$param_required,   # by group risks
        private$contact_risks_$param_required, # by contact risks
        private$augm_$param_required,          # by latent/observation
        private$infec_profile_$param_required) # by infectivity profile
    },

    ############################################################################
    ############################################################################
    ############################################################################

    # @description
    # Whether household data is attached.
    has_data_wide_ = function() {
      length(private$generator_$cols_data_wide) > 0
    },

    # @description
    # Assert if latent definitions match the required.
    # @param wall `logical(1)` to raise error.
    assert_latent_ = function(wall = TRUE, warn_unused = TRUE) {
      private$augm_$assert(
        private$get_all_required_latent_(),
        wall = wall,
        warn_unused = warn_unused && !private$generator_$is_test_model()
      )
    },

    # @description
    # Assert if data_wide headers match implementation file.
    # @param wall `logical(1)` to raise error.
    assert_data_wide_ = function(wall = TRUE, simulating = FALSE) {
      # no data attached at all: capital error at the wall
      if (isTRUE(wall) && !private$has_data_wide_()) stop_h(
        "No household data attached.",
        "Attach it with {.fn add_data_wide}."
      )
      required <- private$get_all_required_data_wide_()

      # observations are generated during simulation: not required
      if (isTRUE(simulating)) {
        required <- setdiff(required, private$augm_$observation_names)
      }

      # a deterministic observation produces a surrogate data column
      cols <- private$generator_$cols_data_wide
      surrogate <- private$augm_$surrogate_names
      cols <- c(cols, intersect(surrogate, required))

      sms <- template_missing_stop(required, cols, "data",
                                   noun = "column", style = "field")
      if (isTRUE(wall) && length(sms) > 0) stop_h(sms)
      else sms
    },

    # @description
    # Assert group risk definitions from the model implementation script.
    # @param wall `logical(1)` to raise error.
    assert_group_risks_ = function(wall = TRUE) {
      private$group_risks_$assert(private$get_all_required_group_risks_(),
                                  wall = wall)
    },

    # @description
    # Assert contact risk definitions are compatible with the model.
    assert_contact_risks_supported_ = function() {
      if (private$generator_$has_contact_risk)
        return(invisible(self))
      stop_h(
        "Contact risks are not supported for model {.val {self$name}}.",
        `!` = "Definitions will be dropped."
      )
    },

    # @description
    # Pin to 1 the baselines whose anchored risks all carry absolute values.
    # An absolute risk leaves no scale for its baseline to take.
    autofill_baseline_params_ = function() {
      anchors <- private$generator_$risk_anchors
      for (name in names(anchors)) {
        if (name %in% private$param_$defs[["name"]]) next
        absolute <- vapply(anchors[[name]]$risks,
                           private$group_risks_$is_absolute, logical(1))
        if (isTRUE(anchors[[name]]$contacts)) {
          absolute <- c(absolute, private$contact_risks_$is_absolute())
        }
        if (any(absolute)) private$param_$add_quo(
          name, value = 1, const = TRUE, .warn_over = FALSE
        )
      }
      invisible(self)
    },

    # @description
    # Assert if parameter definitions match required.
    # @param wall `logical(1)` to raise error.
    # @param warn_unused `logical(1)` warn if unused.
    assert_param_ = function(wall = TRUE, warn_unused = TRUE) {
      private$autofill_baseline_params_()
      private$param_$assert(
        private$get_all_required_param_(),
        wall = wall,
        warn_unused = warn_unused && !private$generator_$is_test_model()
      )
    },

    ############################################################################
    ############################################################################
    ############################################################################

    # @description
    # Assert infectivity profile definitions are compatible with the model.
    assert_infec_profile_supported_ = function() {
      if (private$generator_$has_infec_profile) return(invisible(self))
      stop_h(
        "Infectivity profiles are not supported for model {.val {self$name}}.",
        `!` = "Definitions will be dropped."
      )
    },

    # @description
    # Assert an infectivity profile is defined when the model needs one.
    # @param wall `logical(1)` to raise error.
    assert_infec_profile_ = function(wall = TRUE) {
      if (!private$generator_$has_infec_profile ||
            length(private$infec_profile_$defs) > 0) return(character())
      sms <- c("Missing infectivity profile definition.",
               "Define it with {.fn add_infectivity_profile}.")
      if (isTRUE(wall)) stop_h(sms) else sms
    },

    ############################################################################
    ############################################################################
    ############################################################################

    # @description
    # Reject a risk keyed on an observation column, at either level.
    # A risk feeds the transmission model the observation is generated from, so
    # its strata would depend on the order the definitions were declared in.
    # @param module `R6`. Group or contact risk holder.
    assert_no_obs_risk_collision_ = function(module) {
      obs <- union(private$augm_$observation_names,
                   private$augm_$surrogate_names)
      clash <- intersect(module$feature_names, obs)
      if (length(clash) > 0) stop_h(
        "A {module$family} is keyed on observation{?s} {.val {clash}}.",
        paste("An observation is generated by the model, which the risk feeds:",
              "it cannot stratify the risk it depends on."),
        "Key the {module$family} on a data column."
      )
      invisible(self)
    },

    # @description
    # Helper to assert and re-build group and contact risks.
    helper_build_risks_ = function(module) {
      # proceed only if definitions exist
      if (length(module$defs) == 0) return(invisible(self))

      private$assert_no_obs_risk_collision_(module)

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
    # Assert and pre-build group risks.
    build_group_risks_ = function() {
      private$assert_group_risks_()
      private$helper_build_risks_(private$group_risks_)
    },

    # @description
    # Pre-build contact risks.
    build_contact_risks_ = function() {
      private$helper_build_risks_(private$contact_risks_)
    },

    # @description
    # Factorize and push wide-format data to the C++ engine.
    # Data is assumed validated (non-null data.frame).
    # Definitions are (re)built by the caller via `rebuild_()`.
    set_data_wide_ = function(data) {
      data_factorized <- private$generator_$factorize_data_wide(data)
      # soft check: data may be intentionally incomplete
      # hard check at the wall
      private$assert_data_wide_(wall = FALSE)
      .engine$import_data_wide(private$generator_$mcmc_ptr, data_factorized)
      invisible(self)
    },

    # @description
    # Drop current risk-derived parameters so data change re-derives them.
    # NOTE: reads the pre-change logic tables.
    # NOTE: must run before the next `rebuild_()` risks re-implementation.
    drop_stale_risk_params_ = function() {
      stale <- c(names(private$group_risks_$parameters),
                 names(private$contact_risks_$parameters))
      for (p in stale) private$param_$remove(p)
      invisible(self)
    },

    # @description
    # Reset every definition holder to an empty object.
    reset_definitions_ = function() {
      private$param_ <- HousepiParameter$new()
      private$augm_ <- HousepiAugmented$new()
      private$group_risks_ <- HousepiGroupRisk$new()
      private$contact_risks_ <- HousepiContactRisk$new()
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
    # Drop group risks definitions, removing any associated parameters.
    drop_group_risks_ = function() {
      if (!is.null(private$group_risks_)) {
        for (p in names(private$group_risks_$parameters)) {
          private$param_$remove(p)
        }
        if (!is.null(private$param_$defs) && nrow(private$param_$defs) == 0) {
          private$drop_parameters_()
        }
      }
      private$group_risks_ <- HousepiGroupRisk$new()
      private$rebuild_()
      invisible(self)
    },

    # @description
    # Single, model-owned ordered infected-status dictionary. Data labels win
    # (the data is already coded), else the status process order, else literal
    # character labels declared on latent definitions.
    status_dictionary_ = function() {
      data_lab <- private$generator_$status_data_labels
      if (length(data_lab) > 0) return(data_lab)
      private$augm_$status_dictionary
    },

    # @description
    # Parameter definitions with `discr` resolved against the model on demand.
    resolved_param_ = function() {
      resolve_discr_flags_(private$param_$defs,
                           private$generator_$param_discr_required,
                           is_prior = TRUE)
    },

    # @description
    # Latent definitions with `discr` resolved against the model on demand.
    resolved_latent_ = function() {
      resolve_discr_flags_(
        private$augm_$defs,
        private$generator_$latent_discr_required,
        discr_col_types = private$generator_$discr_data_wide,
        discr_req_cols = private$generator_$data_wide_discr_required
      )
    },

    # @description
    # Whether the latent declaration can be built as it stands. Producers
    # land in the engine on their own decorator, and a status block needs a
    # definition to fill it; either gap waits for the next decorator.
    # NOTE: what the model requires stays out, being supplied by its own
    # decorator too (see `rebuild_`).
    # @param defs `data.frame()`. Declarations with status codes resolved.
    latent_declaration_complete_ = function(defs) {
      # the infectivity profile folds into the latent build
      need_latent <- c(private$augm_$latent_required,
                       private$infec_profile_$latent_required)
      need_param <- c(private$augm_$param_required,
                      private$infec_profile_$param_required)
      all(need_latent %in% private$augm_$defs[["name"]]) &&
        all(need_param %in% private$param_$defs[["name"]]) &&
        private$augm_$covers_statuses(defs,
                                      private$generator_$infec_status_set)
    },

    # @description
    # Soft re-push of R-side definitions to the C++ engine.
    # NOTE: Hard validation lives in `build_full_`.
    # @param wall `logical(1)` Whether to raise error.
    rebuild_ = function(wall = FALSE) {
      # guard: pre-init (before generator is created)
      if (is.null(private$generator_) ||
            is.null(private$generator_$mcmc_ptr)) {
        return(invisible(self))
      }
      # risks stratify on data: nothing to derive while none is attached
      if (length(private$generator_$cols_data_wide) == 0) {
        return(invisible(self))
      }

      # risks derive from their own covariates (they append to param_$defs)
      if (length(private$group_risks_$defs) > 0) {
        private$build_group_risks_()
      }
      if (length(private$contact_risks_$defs) > 0) {
        private$build_contact_risks_()
      }

      # defer the engine push while required data is incomplete
      if (length(private$assert_data_wide_(FALSE, simulating = TRUE)) > 0) {
        return(invisible(self))
      }

      mcmc_ptr <- private$generator_$mcmc_ptr
      info_wide <- private$generator_$info_data_wide

      #> parameters (resolve unspecified discretization against the model)
      private$autofill_baseline_params_()
      .engine$import_param(mcmc_ptr, private$resolved_param_())

      #> group risks
      .engine$import_group_risk(
        mcmc_ptr,
        private$group_risks_$expose_to_cpp(info_wide)
      )

      #> contact risks
      .engine$import_contact_risk(
        mcmc_ptr,
        private$contact_risks_$expose_to_cpp(info_wide)
      )

      #> infectivity profile
      .engine$import_infec_profile(mcmc_ptr, private$infec_profile_$defs)

      #> latent (resolve character status labels to integer codes for C++,
      #> and unspecified discretization against the model and the data)
      status_dict <- private$status_dictionary_()
      if (private$generator_$escaped %in% status_dict) stop_h(
        paste("Infected status {.val {private$generator_$escaped}} collides",
              "with the escaped-status label",
              "({.val {private$generator_$escaped}})."),
        "Rename the status or define a different {.code escaped =} label."
      )
      latent_defs <- resolve_status_codes_(private$resolved_latent_(),
                                           status_dict)
      # NOTE: unresolved definitions wait for completeness or for the wall
      if (wall || private$latent_declaration_complete_(latent_defs)) {
        .engine$import_latent(mcmc_ptr, latent_defs)
      }

      # status dictionary drives the factor refactoring of the status column
      private$generator_$set_status_levels(status_dict)

      invisible(self)
    },

    # @description
    # Normalise the `id_draws` selection into individual and household
    # identifiers plus the default-follow signal. Status filtering and the
    # default household pick are left to the engine.
    resolve_id_draws_ = function(id_draws) {
      out <- list(id = character(), house_id = character())
      # unspecified: let the engine follow one random household
      if (is.null(id_draws) || length(id_draws) == 0L) {
        return(out)
      }
      # explicit opt-out: forward the sentinel so the engine stores nothing
      if (identical(as.character(id_draws), ".none")) {
        out$id <- ".none"
        return(out)
      }
      # named list: split individual and household identifiers
      if (is.list(id_draws)) {
        extra <- setdiff(names(id_draws), c("id", "house_id"))
        if (length(extra) > 0) stop_h(
          "Unknown {.arg id_draws} name{?s}: {.val {extra}}."
        )
        if (!is.null(id_draws$id)) out$id <- as.character(id_draws$id)
        if (!is.null(id_draws$house_id))
          out$house_id <- as.character(id_draws$house_id)
        return(out)
      }
      # bare vector: individual identifiers
      out$id <- as.character(id_draws)
      out
    },

    # @description
    # Full model build. Refresh runtime setup, enforce hard validation,
    # then push the final state to the engine.
    build_full_ = function(steps = 0L, adapt = 0L, burn = 0L,
                           thin = 1L, chains = 1L, threads = 1L,
                           id_draws = NULL, progress = FALSE,
                           simulating = FALSE, ...) {
      # only run_inference resolves a selection; other builds store nothing
      if (!is.list(id_draws))
        id_draws <- list(id = character(), house_id = character())
      .engine$import_setup(
        private$generator_$mcmc_ptr,
        char_setup_input_(c(list(
          steps = steps,
          adapt = adapt,
          burn = burn,
          thin = thin,
          chains = chains,
          threads = threads,
          id_draws = id_draws$id,
          house_draws = id_draws$house_id,
          # generate a seed for Cpp random number generator
          seed = sample.int(.Machine$integer.max, 1L),
          progress = progress
        ), list(...)))
      )
      # hard validation before MCMC
      private$assert_data_wide_(wall = TRUE, simulating = simulating)
      private$assert_param_(wall = TRUE)
      private$assert_group_risks_(wall = TRUE)
      private$assert_infec_profile_(wall = TRUE)
      # engine build reports unresolved read first
      private$rebuild_(wall = TRUE)
      # latent assertion reports plainly missing definitions
      private$assert_latent_(wall = TRUE)
      invisible(self)
    },

    # Assert peak memory needed to materialise draws is available.
    assert_posterior_memory_ = function(steps, burn, thin, chains, id_draws) {
      n_sample <- setup_sample_size_(steps = steps, burn = burn,
                                     thin = thin, chains = chains)$total

      # UPDATE: n. of latent/parameters might be smaller if functions/private
      n_params_defs <- ifelse(is.null(private$param_$defs), 0L,
                              nrow(private$param_$defs))
      n_param_cols <- 2L * n_params_defs
      n_latent_defs <- ifelse(is.null(private$augm_$defs), 0L,
                              nrow(private$augm_$defs))
      # structural upper bound: every named/default household counts as its
      # largest possible size (status filtering happens in the engine)
      house_id <- self$data_wide$house_id
      max_hhsize <- if (length(house_id) == 0) 1L else max(table(house_id))

      n_house_pick <- length(unique(id_draws$house_id))
      n_members_upper <- length(unique(id_draws$id)) + n_house_pick * max_hhsize
      n_latent_cols <- n_members_upper * n_latent_defs
      n_house_total <- length(unique(house_id))
      # total count of columns
      total_cols <- n_param_cols + n_latent_cols + n_house_total + 1L
      # the column-by-sample product exceeds the 32-bit range on large runs
      bytes_final <- 8 * as.double(n_sample) * as.double(total_cols)

      # byte budget for posterior draws materialisation
      # set 80% cap under R vsize, otherwise fall back to 16 GiB
      budget <- 16 * 2^30
      vsize_mb <- tryCatch(mem.maxVSize(), error = function(e) Inf)
      if (is.finite(vsize_mb) && vsize_mb > 0) budget <- 0.8 * vsize_mb * 2^20

      if (bytes_final <= budget) return(invisible())
      stop_h(c(
        paste("Posterior draws would exceed the memory budget",
              "({round(budget / 2^30)} GB)."),
        paste("Consider increasing {.arg thin},",
              "or reducing {.arg {c('steps', 'chains')}}.")
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
        # replicates the single row nrep times; unpinned parameters (NA value)
        # carry no column, so a prior-only parameter is sampled instead
        pinned <- !is.na(self$parameters$value)
        param_table <- as.data.frame(stats::setNames(
          as.list(self$parameters$value[pinned]), self$parameters$name[pinned]
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
    #' @field name `character(1)`
    #'  Model name.
    name = function(value) {
      if (!missing(value)) stop_read_only("name")
      private$generator_$name
    },

    #' @field setup `list`
    #' Last setup used by the build process.
    setup = function(value) {
      if (!missing(value)) stop_read_only("setup")
      .engine$export_setup(private$generator_$mcmc_ptr)
    },

    #' @field data_wide `data.frame`
    #' Individual level data.
    #' One row per individual index by corresponding ID.
    data_wide = function(value) {
      if (!missing(value)) stop_read_only("data_wide")
      # no data attached yet (deferred until add_data_wide() or the wall)
      if (length(private$generator_$cols_data_wide) == 0) return(NULL)
      private$generator_$refactor_data_wide(
        .engine$export_data_wide(private$generator_$mcmc_ptr)
      )
    },

    #' @field data_long `data.frame`
    #' Additional data rows indexed by individual (long format).
    data_long = function(value) {
      if (!missing(value)) stop_read_only("data_long")
      private$generator_$refactor_data_long(
        .engine$export_data_long(private$generator_$mcmc_ptr)
      )
    },

    #' @field parameters `data.frame`
    #' Parameter definitions.
    parameters = function(value) {
      if (!missing(value)) stop_read_only("parameters")
      housepi_defin_showcaser(private$resolved_param_())
    },

    #' @field group_risks `data.frame`
    #' Group risk definitions.
    group_risks = function(value) {
      if (!missing(value)) stop_read_only("group_risks")
      helper_show_risk_table(private$group_risks_,
                             risk_reserved_cols,
                             private$has_data_wide_())
    },

    #' @field contact_risks `data.frame`
    #' Contact risk definitions.
    contact_risks = function(value) {
      if (!missing(value)) stop_read_only("contact_risks")
      helper_show_risk_table(private$contact_risks_,
                             setdiff(risk_reserved_cols, "group_risk"),
                             private$has_data_wide_())
    },

    #' @field latent `data.frame`
    #' Latent definitions.
    latent = function(value) {
      if (!missing(value)) stop_read_only("latent")
      defs <- private$augm_$defs_of_kind(private$resolved_latent_(), "latent")
      housepi_defin_showcaser(defs)
    },

    #' @field observations `data.frame`
    #' Observation definitions.
    #' Column 'level' tells member ("member") from household ("house") entries.
    observations = function(value) {
      if (!missing(value)) stop_read_only("observations")
      resolved <- private$resolved_latent_()
      defs <- do.call(rbind, lapply(c("member", "house"), function(level) {
        out <- private$augm_$defs_of_kind(resolved, level)
        if (is.null(out)) return(NULL)
        cbind(out, level = level, stringsAsFactors = FALSE)
      }))
      housepi_defin_showcaser(defs)
    },

    #' @field infectivity_profiles `data.frame`
    #' Infectivity profile definitions.
    infectivity_profiles = function(value) {
      if (!missing(value)) stop_read_only("infectivity_profiles")
      private$infec_profile_$showcase()
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

      # check data and enable dummy_household behavior
      # NOTE: NULL data is allowed; attach later with add_data_wide()
      data <- assert_household_data(data, model)

      # initialize internal fields
      private$reset_definitions_()
      # generate model wrapper
      private$generator_ <- HousepiGenerator$new(
        model = model, as_template = as_template,
        force = force, quiet = quiet, escaped = escaped, ...
      )
      # import data, if provided (no user definitions to rebuild yet)
      if (!is.null(data)) private$set_data_wide_(data)

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
      for (field in c("param_", "augm_", "group_risks_",
                      "contact_risks_", "infec_profile_")) {
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
      header <- "Household transmission model {.val {self$name}}"

      # append link to template documentation if any
      if (!is.na(get_template_dir(private$generator_$raw_name))) {
        header <- paste(header,
          sprintf("({.run [template](housepi::browse_templates(\"%s\"))})",
                  self$name)
        )
      }
      header <- c(i = paste0(header, "."))

      # UPDATE: to be rationalize in the future if more options are included
      if (private$generator_$is_discrete_time) {
        header <- c(header, "Model is set to {.emph discrete time} mode.")
      }

      # each section is a self-contained cli_bullets block (headline + subs)
      blocks <- list()
      # data
      if (length(private$get_all_required_data_wide_()) > 0) blocks <- c(blocks,
        list(fmt_assert_cli("Data", private$assert_data_wide_(FALSE)))
      )
      # parameters
      if (length(private$get_all_required_param_()) > 0) blocks <- c(blocks,
        list(fmt_assert_cli("Parameters", private$assert_param_(FALSE)))
      )
      # latent
      if (length(private$get_all_required_latent_()) > 0) blocks <- c(blocks,
        list(fmt_assert_cli("Latent", private$assert_latent_(FALSE)))
      )

      # inform about optional risks if dataset is suitable
      cols_extra <- setdiff(private$generator_$cols_data_wide,
                            data_reserved_cols)
      if (length(cols_extra) > 0) {
        # group risks (optional: default when definitions are absent)
        risks <- private$get_all_required_group_risks_()
        if (length(risks) > 0) {
          is_default <- length(template_missing_find(
            risks, names(private$group_risks_$defs)
          )) > 0
          blocks <- c(blocks, list(fmt_assert_cli(
            "Group risks",
            if (is_default) c(
              "Using default group risks.",
              paste("Available:", character_to_string_list(risks))
            ) else NULL,
            minor = TRUE
          )))
        }
        # contact risks (optional)
        if (private$generator_$has_contact_risk) {
          blocks <- c(blocks, list(fmt_assert_cli(
            "Contact risks",
            if (!private$contact_risks_$assert()) "Using default contact risks."
            else NULL,
            minor = TRUE
          )))
        }
      }

      # infectivity profile (required by the model, or defined anyway)
      if (private$generator_$has_infec_profile ||
            length(private$infec_profile_$defs) > 0) blocks <- c(blocks,
        list(fmt_assert_cli("Infectivity profile",
                            private$assert_infec_profile_(FALSE)))
      )

      # order blocks by gravity (error, minor, satisfied), stable within ties
      rank <- c(x = 1L, `!` = 2L, v = 3L)
      order_key <- vapply(blocks, function(b) rank[[names(b)[1]]], integer(1))
      blocks <- blocks[order(order_key)]
      order_key <- sort(order_key)

      # a single hint gathers all default (minor) sections after their block
      sms <- header
      minor <- which(order_key == 2L)
      last_minor <- if (length(minor)) max(minor) else 0L
      for (k in seq_along(blocks)) {
        sms <- c(sms, blocks[[k]])
        # aligned like a sub-line (invisible bullet), no clash with leading `i`
        if (k == last_minor) {
          sms <- c(sms, " " = "A custom choice can be provided.")
        }
      }

      print_bullets(sms)
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
    #' Dataset with associated latent variables for last state.
    #' @return `data.frame`
    export_augmented_data = function() {
      private$generator_$refactor_data_wide(
        .engine$export_augmented_data_wide(private$generator_$mcmc_ptr)
      )
    },

    #' @description
    #' All variables associated to given household for last state.
    #' Entries are named after the homonymous active fields, followed by the
    #' latent, observation and transmission log-likelihood contributions.
    #' An empty definition is dropped unless the model requires it.
    #' @param house_id `r rd_house_id_`
    #' @return `list`
    export_house_snapshot = function(house_id) {
      house_id <- as.character(house_id)
      # export household snapshot
      out <- .engine$export_house_snapshot(
        private$generator_$mcmc_ptr, house_id
      )
      out$data_wide <- private$generator_$refactor_data_wide(out$data_wide)
      ## data long
      if (all(dim(out$data_long) == 0)) out$data_long <- NULL
      else out$data_long <- private$generator_$refactor_data_long(out$data_long)
      ## parameters
      if (length(private$generator_$param_required) == 0 &&
            all(dim(out$parameters) == 0)) out$parameters <- NULL
      ## group risks
      if (length(private$generator_$group_risks_required) == 0 &&
            all(dim(out$group_risks) == 0)) out$group_risks <- NULL
      ## conctact risks
      if (!private$generator_$has_contact_risk &&
            all(dim(out$contact_risks) == 0)) out$contact_risks <- NULL
      ## infectivity profiles
      if (!private$generator_$has_infec_profile &&
            length(out$infectivity_profiles) == 0) {
        out$infectivity_profiles <- NULL
      } else {
        lst <- out$infectivity_profiles
        df <- data.frame(id = lst$id, discr = lst$discr, func = lst$func,
                         min = lst$min, max = lst$max)
        df$func_args <- lst$func_args
        out$infectivity_profiles <- df
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
  pvt$build_group_risks_()
  pvt$build_contact_risks_()

  gen <- pvt$generator_
  info_wide <- gen$info_data_wide
  utils::write.csv(pvt$augm_$defs,
                   file.path(outpath, "latent.csv"),
                   quote = FALSE, row.names = FALSE)
  utils::write.csv(pvt$contact_risks_$expose_to_cpp(info_wide),
                   file.path(outpath, "contact_risks.csv"),
                   quote = FALSE, row.names = FALSE)
  utils::write.csv(.engine$export_data_wide(gen$mcmc_ptr),
                   file.path(outpath, "data_wide.csv"),
                   quote = FALSE, row.names = FALSE)
  utils::write.csv(.engine$export_data_long(gen$mcmc_ptr),
                   file.path(outpath, "data_long.csv"),
                   quote = FALSE, row.names = FALSE)
  utils::write.csv(pvt$group_risks_$expose_to_cpp(info_wide),
                   file.path(outpath, "group_risks.csv"),
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
