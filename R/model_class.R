## housepi - Epidemiological Analysis of Household Outbreaks.
## Copyright 2024-2026 Institut Pasteur (Paris), France.
## SPDX-License-Identifier: GPL-3.0-or-later
##

########  HousepiModel  ########

# auxiliary environment
.housepi_env <- new.env(parent = emptyenv())
.housepi_env$hollow <- FALSE

#' @title
#' R6 class for compiling and running a household transmission model
#' @description
#' Compile the household transmission model, manage its data and definitions,
#' run inference and simulations, from a single class.
#' @examples
#' # created and modified through new_housepi_model() and the add_*() verbs
#' hh_model <- flu_1918_frederick |>
#'   new_housepi_model("final-size") |>
#'   add_parameter(import = "template")
#' hh_model
#' hh_model$parameters
#' @docType class
HousepiModel <- R6::R6Class(
  classname = "HousepiModel",
  cloneable = FALSE,
  lock_class = TRUE,

  ####  Private  ####

  private = list(
    gen_ = NULL,
    param_ = NULL,
    augm_ = NULL,
    group_risks_ = NULL,
    contact_risks_ = NULL,
    status_ = NULL,
    infec_profile_ = NULL,
    reporting_process_ = NULL,
    collapse_seen_ = character(),

    ####  Requirements  ####

    # Names of one source the latent declaration reads,
    # across latent, observation and infectivity profile definitions.
    # * read from the resolved tables -> a reference counts under its kind;
    # * each name -> labelled with the family of the definition reading it,
    #   member- and household-level observations alike.
    # @param type `character(1)` Source letter.
    required_of_declaration_ = function(type) {
      defs <- private$resolved_latent_()
      family <- vapply(sub("^house$", "member", defs[["kind"]]), kind2family,
                       character(1), USE.NAMES = FALSE)
      by_kind <- lapply(unique(family), function(f) {
        find_required_definitions(defs[family == f, , drop = FALSE], type,
                                  c("const", "func", "model", "min", "max"), f)
      })
      c(unlist(by_kind),
        find_required_definitions(private$resolved_profile_(), type,
                                  c("func", "min", "max"),
                                  private$infec_profile_$family))
    },

    # Latent variables required by the model and by other definitions.
    get_required_latent_ = function() {
      c(private$gen_$latent_required,      # by the model
        private$required_of_declaration_("A"))   # by the latent declaration
    },

    # Data columns required by the model and by other definitions.
    get_required_data_wide_ = function() {
      c(private$gen_$data_required,        # by the model
        private$group_risks_$data_required,      # by group risks
        private$contact_risks_$data_required,    # by contact risks
        private$required_of_declaration_("D"))   # by the latent declaration
    },

    # Group risks required by the model and by other definitions.
    get_required_group_risks_ = function() {
      c(private$gen_$group_risks_required, # by the model
        private$required_of_declaration_("G"))   # by the latent declaration
    },

    # Parameters required by the model and by other definitions.
    get_required_parameters_ = function() {
      c(private$gen_$param_required,       # by the model
        find_required_definitions(private$resolved_param_(), "P",
                                  c("const", "func", "min", "max"),
                                  private$param_$family), # by param functions
        private$group_risks_$param_required,     # by group risks
        private$contact_risks_$param_required,   # by contact risks
        private$required_of_declaration_("P"))   # by the latent declaration
    },

    ####  Assertions  ####

    # Whether household data is attached.
    has_data_wide_ = function() {
      length(private$gen_$cols_data_wide) > 0
    },

    # Assert the required latents are defined.
    # * `wall` -> error, else the messages are returned.
    assert_latent_ = function(wall = TRUE, warn_unused = TRUE) {
      private$augm_$assert(
        private$get_required_latent_(),
        wall = wall,
        warn_unused = warn_unused && !private$gen_$is_test_model()
      )
    },

    # Assert data_wide holds the columns required by the model and definitions.
    # * `wall` -> error, else the messages are returned;
    # * `simulating` -> generated observations are not required.
    assert_data_wide_ = function(wall = TRUE, simulating = FALSE) {
      # no data attached -> error at the wall
      if (isTRUE(wall) && !private$has_data_wide_()) stop_h(
        "No household data attached.",
        "Attach it with {.fn add_data_wide}."
      )
      required <- private$get_required_data_wide_()

      # simulating -> observations are generated, not required
      if (isTRUE(simulating)) {
        required <- setdiff(required, private$augm_$observation_names)
      }

      # a deterministic observation produces a surrogate data column
      cols <- private$gen_$cols_data_wide
      surrogate <- private$augm_$surrogate_names
      cols <- c(cols, intersect(surrogate, required))

      sms <- template_missing_stop(required, cols, "data",
                                   noun = "column", style = "field")
      wall_or_return_(sms, wall)
    },

    # Assert the model reads every group risk defined.
    assert_group_risks_ = function() {
      private$group_risks_$assert(private$get_required_group_risks_())
    },

    # Assert the model supports a feature.
    # @param label `character(1)` Feature name, plural: "Contact risks".
    assert_feature_supported_ = function(supported, label) {
      if (supported) return(invisible(self))
      stop_h(
        "{label} are not supported for model {.val {self$name}}.",
        "Choose a model that supports them."
      )
    },

    # Pin baselines to the reference value, keyed by parameter name.
    # * absolute risk -> baseline held at the reference value, declared or not;
    # * baseline defined as a function -> not pinned, user's own scale.
    # @return `list` of `character()` Risks pinning each baseline.
    baseline_pins_ = function() {
      anchors <- private$gen_$risk_anchors
      defs <- private$param_$defs
      pins <- list()
      for (name in names(anchors)) {
        # a baseline defined as a function carries the user's own scale
        row <- match(name, defs[["name"]])
        if (!is.na(row) && !is.na(defs[["func"]][row])) next
        risks <- anchors[[name]]$risks
        absolute <- risks[vapply(risks, private$group_risks_$is_absolute,
                                 logical(1))]
        if (isTRUE(anchors[[name]]$contacts) &&
              private$contact_risks_$is_absolute()) {
          absolute <- c(absolute, "contact")
        }
        if (length(absolute)) pins[[name]] <- absolute
      }
      pins
    },

    # Declared parameters, with the baselines absolute risks anchor
    # held at the reference value.
    # @return `HousepiParameter`
    effective_param_ = function() {
      pins <- names(private$baseline_pins_())
      if (!length(pins)) return(private$param_)
      param <- private$param_$clone(deep = TRUE)
      for (name in pins) param$add_quo(name, const = 1, .warn_over = FALSE)
      # a pin reads nothing, so leading keeps it ahead of whatever reads it
      is_pin <- param$defs[["name"]] %in% pins
      param$defs <- param$defs[order(!is_pin), , drop = FALSE]
      rownames(param$defs) <- NULL
      param
    },

    # Fingerprint of the definition row a baseline was declared with.
    # * not declared, or already the reference value -> `NA`.
    # @param name `character(1)` Baseline parameter name.
    declared_scale_ = function(name) {
      defs <- private$param_$defs
      row <- match(name, defs[["name"]])
      if (is.na(row)) return(NA_character_)
      is_reference <- is.na(defs[["prior"]][row]) &&
        is.na(defs[["func"]][row]) &&
        isTRUE(showcaser_extract_numeric(defs[["const"]][row]) == 1)
      if (is_reference) NA_character_
      else paste(unlist(defs[row, ]), collapse = "|")
    },

    # Warn once per declaration for each baseline an absolute risk absorbed,
    # discarding what its declaration asked for.
    # * `collapse_seen_` -> what was last reported;
    # * baseline released by its risks -> forgotten.
    warn_baseline_collapse_ = function() {
      pins <- private$baseline_pins_()
      defs <- private$param_$defs
      declared <- vapply(names(pins), private$declared_scale_, character(1))
      for (name in names(pins)) {
        if (is.na(declared[[name]]) ||
              identical(unname(private$collapse_seen_[name]),
                        declared[[name]])) next
        risks <- pins[[name]]
        value <- showcaser_extract_numeric(
          defs[["const"]][match(name, defs[["name"]])]
        )
        warn_h(
          paste("Absolute {cli::qty(risks)} risk{?s} {.val {risks}}",
                "{?has/have} no reference."),
          `!` = "Baseline {.val {name}} was turned into constant reference."
        )
      }
      private$collapse_seen_ <- declared
      invisible(self)
    },

    # Read a table of parameter values a run supplies.
    # * value for a baseline an absolute risk holds at the reference -> error.
    # @param table `data.frame` | path Values keyed by parameter name.
    # @param nrep `integer(1)` Rows kept at most, evenly spaced, all if `NA`.
    # @param thin_ess `logical(1)` Whether to thin by effective sample size.
    # @param arg `character(1)` Argument the table came from.
    read_param_table_ = function(table,
                                 nrep = NA_integer_,
                                 thin_ess = TRUE,
                                 arg = "draws") {
      table <- helper_read_param_table_(table, thin_ess = thin_ess)
      private$assert_no_baseline_value_(table, arg)
      helper_thin_even_(table, nrep)
    },

    # Assert no value a run supplies targets a baseline an absolute risk holds
    # at the reference.
    # * `NA` cell -> passes, it supplies nothing.
    # @param given `list` | `data.frame` Values keyed by parameter name.
    # @param arg `character(1)` Argument the values came from.
    assert_no_baseline_value_ = function(given, arg) {
      named <- intersect(names(given), names(private$baseline_pins_()))
      held <- named[!vapply(given[named], \(v) all(v[!is.na(v)] == 1),
                            logical(1))]
      if (!length(held)) return(invisible(self))
      stop_h(
        "{.arg {arg}} sets {.val {held}} to a value other than 1.",
        paste("An absolute risk leaves no scale for",
              "{cli::qty(held)}{?it/them}: fold the factor into the risk",
              "values, or name a {.code reference} group.")
      )
    },

    # Assert the required parameters are defined.
    # * `wall` -> error, else the messages are returned;
    # * `warn_unused` -> warning on unused ones.
    assert_param_ = function(wall = TRUE, warn_unused = TRUE) {
      private$effective_param_()$assert(
        private$get_required_parameters_(),
        wall = wall,
        warn_unused = warn_unused && !private$gen_$is_test_model()
      )
    },

    # Assert every reference names one declared definition, scanning the
    # resolved tables for terminals the dictionary left unsettled:
    # * name two holders claim, or nothing declares -> error;
    # * `wall` -> error, else the messages are returned.
    assert_sources_ = function(wall = TRUE) {
      # the full source-bearing column set, not the narrower requirement one:
      # a distribution or proposal argument carries references too
      cols <- c("const", "func", "prior", "model", "min", "max", "proposal")
      scan <- \(defs) find_required_definitions(defs, "?", cols, NA_character_)
      unresolved <- unique(c(scan(private$resolved_latent_()),
                             scan(private$resolved_profile_()),
                             scan(private$resolved_param_())))
      if (!length(unresolved)) return(character())

      sets <- private$source_sets_()
      declared <- unique(unlist(sets, use.names = FALSE))
      labels <- unname(source_labels_[c("A", "D", "P", "G")])
      nothing <- inline_cli_("No {.or {labels}} is declared under this name.")

      sms <- character()
      for (nm in unresolved) {
        claims <- names(sets)[vapply(sets, \(x) nm %in% x, logical(1))]
        sms <- c(sms,
                 inline_cli_("Could not resolve the source of {.val {nm}}."))
        if (length(claims) > 1) {
          sms <- c(sms,
                   inline_cli_("Declared both as {source_labels_[claims]}."),
                   "Rename one of them.")
        } else {
          sms <- c(sms, nothing,
                   find_closest_match(nm, declared, list_options = FALSE))
        }
      }
      wall_or_return_(sms, wall)
    },

    # Assert an infectivity profile is defined when the model needs one.
    # * `wall` -> error, else the messages are returned.
    assert_infec_profile_ = function(wall = TRUE) {
      if (!private$gen_$has_infec_profile ||
            length(private$infec_profile_$defs)) return(character())
      sms <- c("Missing infectivity profile definition.",
               "Define it with {.fn add_infectivity_profile}.")
      wall_or_return_(sms, wall)
    },

    ####  Build  ####

    # Assert no risk is keyed on an observation column, at either level.
    # A risk feeds the transmission model the observation is generated from,
    # so its strata would depend on the order of declaration.
    # @param module `R6` Group or contact risk holder.
    assert_no_obs_risk_collision_ = function(module) {
      obs <- union(private$augm_$observation_names,
                   private$augm_$surrogate_names)
      clash <- intersect(module$covariate_names, obs)
      if (length(clash)) stop_h(
        "A {module$family} is keyed on observation{?s} {.val {clash}}.",
        paste("An observation is generated by the model, which the risk feeds:",
              "it cannot stratify the risk it depends on."),
        "Key the {module$family} on a data column."
      )
      invisible(self)
    },

    # Assert and rebuild group or contact risks.
    # @param ... Stratification axes only the group risks know about.
    helper_build_risks_ = function(module, ...) {
      # proceed only if definitions exist
      if (!length(module$defs)) return(invisible(self))

      private$assert_no_obs_risk_collision_(module)

      module$implement(self$data_wide, ...)

      # append parameter if missing
      new_parameters <- module$parameters
      # find missing parameters
      mask <- setdiff(names(new_parameters), private$param_$defs$name)
      new_parameters <- new_parameters[mask]
      if (length(new_parameters)) {
        xquos <- convert_list_to_quosures(new_parameters)
        private$param_$add_quo_list(xquos, append = TRUE)
      }

      invisible(self)
    },

    # Assert and pre-build group risks.
    build_group_risks_ = function(wall = FALSE) {
      status_codes <- private$status_codes_()
      private$assert_group_risks_()
      private$helper_build_risks_(
        private$group_risks_, status_codes, wall,
        stratified = private$gen_$risk_status_stratified
      )
    },

    # Pre-build contact risks.
    build_contact_risks_ = function() {
      private$helper_build_risks_(private$contact_risks_)
    },

    # Factorize and push wide-format data to the engine.
    # The data is validated, and the caller rebuilds definitions via rebuild_().
    set_data_wide_ = function(data) {
      data_factorized <- private$gen_$factorize_data_wide(data)
      # soft check: data may be intentionally incomplete
      # hard check at the wall
      private$assert_data_wide_(wall = FALSE)
      .engine$import_data_wide(private$gen_$mcmc_ptr, data_factorized)
      invisible(self)
    },

    # Drop the risk-derived parameters, re-derived after a data change.
    # Reads the pre-change logic tables, so it runs before the next rebuild_().
    drop_stale_risk_params_ = function() {
      stale <- c(names(private$group_risks_$parameters),
                 names(private$contact_risks_$parameters))
      for (p in stale) private$param_$remove(p)
      invisible(self)
    },

    # Reset every definition holder to an empty object.
    reset_definitions_ = function() {
      private$param_ <- HousepiParameter$new()
      private$augm_ <- HousepiAugmented$new()
      private$group_risks_ <- HousepiGroupRisk$new()
      private$contact_risks_ <- HousepiContactRisk$new()
      private$status_ <- HousepiStatusProcess$new()
      private$infec_profile_ <- HousepiInfecProfile$new()
      private$reporting_process_ <- list(wide = NULL, long = NULL)
      private$collapse_seen_ <- character()
      invisible(self)
    },

    # Drop parameter definitions.
    drop_parameters_ = function() {
      private$param_ <- HousepiParameter$new()
      private$rebuild_()
      invisible(self)
    },

    # Drop group risks definitions, removing any associated parameters.
    drop_group_risks_ = function() {
      if (!is.null(private$group_risks_)) {
        for (p in names(private$group_risks_$parameters)) {
          private$param_$remove(p)
        }
        if (!is.null(private$param_$defs) && !nrow(private$param_$defs)) {
          private$drop_parameters_()
        }
      }
      private$group_risks_ <- HousepiGroupRisk$new()
      private$rebuild_()
      invisible(self)
    },

    # Ordered infected-status dictionary of the model, in order of precedence:
    # * data labels, the data being already coded;
    # * status process order;
    # * character labels declared on latent definitions.
    status_dictionary_ = function() {
      data_lab <- private$gen_$status_data_labels
      if (length(data_lab)) return(data_lab)
      if (length(private$status_$labels)) return(private$status_$labels)
      private$augm_$status_labels
    },

    # Engine infected-status codes, named after the labels a risk may key on.
    # * model knows labels -> labels win;
    # * integer-coded statuses -> their own dictionary, the codes the data
    #   carries and the status process draws, each standing for itself.
    status_codes_ = function() {
      dict <- private$status_dictionary_()
      if (length(dict)) return(stats::setNames(seq_along(dict), dict))
      codes <- sort(union(as.integer(private$gen_$infec_status_set),
                          private$status_$codes))
      stats::setNames(codes, as.character(codes))
    },

    # Declared names under each source, keyed by its source letter.
    source_sets_ = function() {
      gen <- private$gen_
      live <- !is.null(gen) && !is.null(gen$mcmc_ptr)
      # what the transmission model declares it needs, whether or not the
      # matching column or parameter has landed yet
      needs <- \(field) if (live) unname(gen[[field]]) else character()

      defs <- private$augm_$defs
      latent <- character()
      obs <- character()
      if (!is.null(defs)) {
        kind <- defs[["kind"]]
        latent <- defs[["name"]][kind == "latent"]
        obs <- defs[["name"]][kind %in% c("member", "house")]
      }
      cols <- character()
      if (live) cols <- setdiff(gen$cols_data_wide, data_reserved_cols)

      sets <- list(
        # an observation stands for a column, at either level:
        # datum in inference, generated value in simulation
        D = c(cols, needs("data_required"), obs),
        G = c(names(private$group_risks_$defs), needs("group_risks_required")),
        P = c(private$param_$defs[["name"]], needs("param_required")),
        A = c(latent, needs("latent_required"))
      )
      # latent shadowing the like-named column -> no ambiguity,
      # as a simulated output fed back as data
      sets$D <- setdiff(sets$D, sets$A)
      lapply(sets, function(nms) {
        nms <- unique(nms)
        nms[!is.na(nms) & nzchar(nms)]
      })
    },

    # Source letter a reference resolves to, for every name the model declares.
    # * name two holders claim -> `""`, no reference settles it.
    source_dictionary_ = function() {
      sets <- private$source_sets_()
      dict <- character()
      for (src in names(sets)) {
        nms <- sets[[src]]
        dict[intersect(nms, names(dict))] <- ""
        dict[setdiff(nms, names(dict))] <- src
      }
      dict
    },

    # Parameter definitions with sources and `discr` resolved on demand.
    resolved_param_ = function() {
      resolve_discr_flags_(
        resolve_source_types_(private$effective_param_()$defs,
                              private$source_dictionary_()),
        private$gen_$param_discr_required,
        is_prior = TRUE
      )
    },

    # Latent definitions with sources and `discr` resolved on demand.
    # Sources come first:
    # * bare read of a data column -> inherits its discreteness from `D`;
    # * member-level observation -> household values read once sources known.
    resolved_latent_ = function() {
      resolve_discr_flags_(
        helper_cast_house_level(
          resolve_source_types_(private$augm_$defs,
                                private$source_dictionary_())
        ),
        private$gen_$latent_discr_required,
        discr_col_types = private$gen_$discr_data_wide,
        discr_req_cols = private$gen_$data_wide_discr_required
      )
    },

    # Infectivity profile definitions with sources resolved on demand.
    resolved_profile_ = function() {
      resolve_source_types_(private$infec_profile_$defs,
                            private$source_dictionary_())
    },

    # Whether the latent declaration can be built as it stands.
    # * producer not landed yet -> waits for its own decorator;
    # * status block without a definition -> waits for the next decorator.
    # What the model requires stays out, supplied by its own decorator too.
    # @param defs `data.frame()` Declarations with status codes resolved.
    latent_declaration_complete_ = function(defs) {
      # the infectivity profile folds into the latent build
      need_latent <- private$required_of_declaration_("A")
      need_param <- private$required_of_declaration_("P")
      all(need_latent %in% private$augm_$defs[["name"]]) &&
        all(need_param %in% private$effective_param_()$defs[["name"]]) &&
        private$augm_$covers_statuses(defs, private$carried_statuses_())
    },

    # Re-push the R-side definitions to the engine, without hard validation,
    # then report the baselines the risks absorbed.
    # Hard validation runs in build_full_().
    # @param wall `logical(1)` Whether to raise an error.
    rebuild_ = function(wall = FALSE) {
      private$push_definitions_(wall)
      private$warn_baseline_collapse_()
      invisible(self)
    },

    # Derive the risks from the data,
    # then push every definition the engine can already take.
    # @param wall `logical(1)` Whether to raise an error.
    # @param derive `logical(1)` Whether to re-derive the risk tables.
    push_definitions_ = function(wall = FALSE, derive = TRUE) {
      # no data attached -> no risk to derive
      if (!length(private$gen_$cols_data_wide)) return(invisible(self))

      # settle the status dictionary first, read by the risks stratifying on it
      # and by the latent status codes below
      status_dict <- private$status_dictionary_()
      if (private$gen_$escaped %in% status_dict) stop_h(
        paste("Infected status {.val {private$gen_$escaped}} collides",
              "with the escaped-status label",
              "({.val {private$gen_$escaped}})."),
        "Rename the status or define a different {.code escaped =} label."
      )

      # risks derive from their own covariates (they append to param_$defs)
      if (derive && length(private$group_risks_$defs))
        private$build_group_risks_(wall)

      if (derive && length(private$contact_risks_$defs))
        private$build_contact_risks_()

      # defer the engine push while required data is incomplete
      if (length(private$assert_data_wide_(FALSE, simulating = TRUE)) > 0)
        return(invisible(self))

      # forward reference -> waits for the decorator that declares it,
      # never pushing an unsettled source
      if (!wall && length(private$assert_sources_(FALSE)) > 0)
        return(invisible(self))

      mcmc_ptr <- private$gen_$mcmc_ptr
      info_wide <- private$gen_$info_data_wide

      # parameters
      .engine$import_param(mcmc_ptr, private$resolved_param_())

      # group risks
      .engine$import_group_risk(
        mcmc_ptr,
        private$group_risks_$expose_to_cpp(info_wide, private$status_codes_())
      )

      # contact risks
      .engine$import_contact_risk(
        mcmc_ptr,
        private$contact_risks_$expose_to_cpp(info_wide)
      )

      # infectivity profile
      .engine$import_infec_profile(mcmc_ptr, private$resolved_profile_())

      # status process
      private$push_status_process_(derive)

      # latent
      latent_defs <- resolve_status_codes_(private$resolved_latent_(),
                                           status_dict)
      # unresolved definitions wait for completeness or for the wall
      if (wall || private$latent_declaration_complete_(latent_defs))
        .engine$import_latent(mcmc_ptr, latent_defs)

      # status dictionary drives the factor refactoring of the status column
      private$gen_$set_status_levels(status_dict)

      invisible(self)
    },

    # Check the status-process table against the data, then push it coded.
    # * labelled statuses -> read the model dictionary;
    # * codes -> stand for themselves.
    # @param derive `logical(1)` Whether to check the table against the data.
    push_status_process_ = function(derive = TRUE) {
      if (!private$status_$is_defined) return(invisible(self))
      levels <- private$status_draw_levels_()
      if (derive) private$status_$assert_fit(self$data_wide, levels)
      .engine$import_status_process(
        private$gen_$mcmc_ptr,
        private$status_$expose_to_cpp(private$gen_$info_data_wide,
                                      private$status_codes_()),
        length(levels)
      )
      invisible(self)
    },

    # Infected levels a member may draw:
    # * labelled status process -> the model dictionary;
    # * otherwise -> its integer statuses.
    status_draw_levels_ = function() {
      private$status_$draw_levels(private$status_codes_())
    },

    # Infected status codes a member can carry:
    # * status process -> every level it draws, whatever the data carries;
    # * otherwise -> the ones observed in the data.
    carried_statuses_ = function() {
      if (!private$status_$is_defined) return(private$gen_$infec_status_set)
      seq_along(private$status_draw_levels_())
    },

    # Build the model: refresh the run setup, enforce hard validation,
    # then push the final state to the engine.
    build_full_ = function(steps = 0L, adapt = 0L, burn = 0L,
                           thin = 1L, chains = 1L, threads = 1L,
                           id_draws = NULL, simulating = FALSE,
                           seed = draw_seed_(), ...) {
      # only run_inference resolves a selection, other builds store nothing
      if (!is.list(id_draws))
        id_draws <- list(id = character(), hh_id = character())

      .engine$import_setup(
        private$gen_$mcmc_ptr,
        char_setup_input_(c(list(
          steps = steps,
          adapt = adapt,
          burn = burn,
          thin = thin,
          chains = chains,
          threads = threads,
          id_draws = id_draws$id,
          house_draws = id_draws$hh_id,
          seed = seed
        ), list(...)))
      )
      # hard validation before MCMC
      private$assert_data_wide_(wall = TRUE, simulating = simulating)
      # a column not yet attached is the likeliest cause of an unsettled source
      private$assert_sources_(wall = TRUE)
      private$assert_param_(wall = TRUE)
      private$assert_group_risks_()
      private$assert_infec_profile_(wall = TRUE)
      # engine build reports unresolved read first
      private$rebuild_(wall = TRUE)
      # latent assertion reports plainly missing definitions
      private$assert_latent_(wall = TRUE)
      invisible(self)
    },

    # Build the model for simulation.
    # @todo build differently given a param table
    helper_simulation_build_ = function(draws = NULL,
                                        nrep = NA_integer_,
                                        thin_ess = TRUE,
                                        ...) {
      nrep <- assert_count_(nrep, "Number of replicates", allow_na = TRUE)

      # read param table
      param_table <- NULL
      if (is.null(draws)) {
        if (is.na(nrep)) nrep <- 1L
      } else {
        param_table <- private$read_param_table_(draws, nrep, thin_ess)
        nrep <- nrow(param_table)
        # an undeclared model builds on the table's first row: the engine
        # keeps it for the run, the caller's definitions stay as declared
        if (!nrow(self$parameters)) {
          declared <- private$param_$clone(deep = TRUE)
          on.exit(private$param_ <- declared, add = TRUE)
          param_dummy_list <- as.list(utils::head(param_table, 1))
          private$param_$add_quo_list(
            convert_list_to_quosures(param_dummy_list), append = TRUE
          )
          private$rebuild_()
        }
      }

      # transmission sets the latents it generates in simulation:
      # an undeclared one gets a bare row, whose bounds simulation never reads
      generated <- setdiff(private$gen_$latent_transmission,
                           private$augm_$defs[["name"]])
      if (length(generated)) {
        declared_augm <- private$augm_$clone(deep = TRUE)
        on.exit(private$augm_ <- declared_augm, add = TRUE)
        bare <- rep(list(list(min = 0, max = 1)), length(generated))
        private$augm_$add_quo_list(
          convert_list_to_quosures(stats::setNames(bare, generated)),
          append = TRUE, .kind = "latent"
        )
      }

      # internal build function
      private$build_full_(steps = nrep, simulating = TRUE, ...)

      # there may be no parameters at all
      if (!nrow(self$parameters)) return(nrep)

      # no column: a parameter is drawn from its prior
      if (is.null(param_table))
        param_table <- data.frame(row.names = seq_len(nrep))
      # attach table with parameters
      .engine$import_param_table(private$gen_$mcmc_ptr, param_table)

      nrep
    },

    # Set reporting process `slot`, "wide" or "long", to `func`.
    # * `func` NULL -> removes it.
    set_reporting_process_ = function(slot, func) {
      if (!is.null(func) && !is.function(func)) stop_h(
        "Reporting process {.arg {slot}} must be a function or {.code NULL}.",
        format_type_error(func)
      )
      private$reporting_process_[slot] <- list(func)
      invisible(self)
    },

    # Build the function relabelling each simulated `data_wide`,
    # then applying the wide reporting process.
    # * not `labelled`, no reporting process -> frames stay engine-coded,
    #   for consumers reading `hh_id` and `status` only.
    # @todo drop the relabelling of variant data
    helper_reporting_wide_ = function(labelled = TRUE) {
      gen <- private$gen_
      wide <- private$reporting_process_$wide
      if (is.function(wide)) \(x) wide(gen$refactor_data_wide(x))
      else if (labelled) \(x) gen$refactor_data_wide(x)
    },

    # Run a measure function on posterior simulations.
    # @todo loop in C++
    run_measure_ = function(draws,
                            apply_wide,
                            apply_long = NULL,
                            nrep = NA_integer_,
                            thin_ess = TRUE,
                            threads = 1L,
                            progress = TRUE,
                            labelled = TRUE) {
      param_table <- private$read_param_table_(draws, nrep, thin_ess)
      nsim <- nrow(param_table)

      # internal build function
      private$build_full_(steps = nsim, threads = threads)

      # attach param_table
      .engine$import_param_table(private$gen_$mcmc_ptr, param_table)

      long <- private$reporting_process_$long
      apply_long_func <-
        if (is.null(apply_long)) NULL else \(x) apply_long(long(x))
      with_progress_(progress, "simulate", 1, function(bar) {
        .engine$run_simulation(private$gen_$mcmc_ptr,
                               private$helper_reporting_wide_(labelled),
                               apply_wide, apply_long_func, bar)
      })
    }
  ),

  ####  Active  ####

  active = list(
    #' @field name `character(1)` Model name.
    name = function(value) {
      assert_read_only_(value, "name")
      private$gen_$name
    },

    #' @field escaped `character(1)`
    #' Label for status code `0` in a character or factor `status` column.
    escaped = function(value) {
      assert_read_only_(value, "escaped")
      private$gen_$escaped
    },

    #' @field setup `list` Last setup used by the build process.
    setup = function(value) {
      assert_read_only_(value, "setup")
      .engine$export_setup(private$gen_$mcmc_ptr)
    },

    #' @field data_wide `data.frame`
    #' Individual-level data, one row per individual, indexed by ID.
    data_wide = function(value) {
      assert_read_only_(value, "data_wide")
      # no data attached yet (deferred until add_data_wide() or the wall)
      if (!length(private$gen_$cols_data_wide)) return(NULL)
      .engine$export_data_wide(private$gen_$mcmc_ptr) |>
        private$gen_$refactor_data_wide()
    },

    #' @field data_long `data.frame`
    #' Additional data rows indexed by individual (long format).
    data_long = function(value) {
      assert_read_only_(value, "data_long")
      .engine$export_data_long(private$gen_$mcmc_ptr) |>
        private$gen_$refactor_data_long()
    },

    #' @field parameters `data.frame` Parameter definitions.
    parameters = function(value) {
      assert_read_only_(value, "parameters")
      private$resolved_param_() |>
        housepi_defin_showcaser() |>
        sort_definitions_()
    },

    #' @field group_risks `data.frame` Group risk definitions.
    group_risks = function(value) {
      assert_read_only_(value, "group_risks")
      private$group_risks_ |>
        helper_show_risk_table(setdiff(risk_reserved_cols, risk_pin_cols),
                               private$has_data_wide_()) |>
        sort_definitions_("group_risk", by_status = FALSE)
    },

    #' @field contact_risks `data.frame` Contact risk definitions.
    contact_risks = function(value) {
      assert_read_only_(value, "contact_risks")
      private$contact_risks_ |>
        helper_show_risk_table(c("param_name", "n", "rule_group"),
                               private$has_data_wide_()) |>
        sort_definitions_("param_name", by_status = FALSE)
    },

    #' @field latent `data.frame` Latent definitions.
    latent = function(value) {
      assert_read_only_(value, "latent")
      private$resolved_latent_() |>
        augm_defs_of_kind("latent") |>
        housepi_defin_showcaser() |>
        sort_definitions_()
    },

    #' @field observations `data.frame`
    #' Observation definitions.
    #' Column 'level' tells member ("member") from household ("house") entries.
    observations = function(value) {
      assert_read_only_(value, "observations")
      resolved <- private$resolved_latent_()
      do.call(rbind, lapply(c("member", "house"), function(level) {
        out <- augm_defs_of_kind(resolved, level)
        if (is.null(out)) return(NULL)
        # pinned to its datum, the like-named column
        out[["const"]] <- NULL
        cbind(out, level = level)
      })) |>
        housepi_defin_showcaser() |>
        sort_definitions_()
    },

    #' @field status_process `data.frame`
    #' Status process, one row per rule: its status, weight,
    #' the members it applies to, its rule group and covariates.
    status_process = function(value) {
      assert_read_only_(value, "status_process")
      data_wide <- self$data_wide
      levels <- if (!is.null(data_wide)) private$status_draw_levels_()
      private$status_$tidy_table(data_wide, levels)
    },

    #' @field infectivity_profiles `data.frame` Infectivity profile definitions.
    infectivity_profiles = function(value) {
      assert_read_only_(value, "infectivity_profiles")
      private$resolved_profile_() |>
        private$infec_profile_$showcase() |>
        sort_definitions_()
    },

    #' @field reporting_process `list`
    #' Reporting process: functions `wide` and `long`, `NULL` when absent.
    reporting_process = function(value) {
      assert_read_only_(value, "reporting_process")
      private$reporting_process_
    }
  ),

  ####  Public  ####

  public = list(
    #' @description
    #' Create a `HousepiModel`, the object inference and simulation run on.
    #' @param data `r rd_add_wide_`
    #' If omitted, attach it later with [add_data_wide()].
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
      # construction option for copying
      if (isTRUE(.housepi_env$hollow)) return(invisible(self))

      # check the data, NULL allowed and attached later with add_data_wide()
      data <- assert_household_data(data, model)

      # initialize internal fields
      private$reset_definitions_()
      # generate model wrapper
      private$gen_ <- HousepiGenerator$new(
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
      for (field in c("param_", "augm_", "group_risks_", "contact_risks_",
                      "status_", "infec_profile_", "reporting_process_",
                      "collapse_seen_")) {
        obj_priv[[field]] <- helper_clone(private[[field]])
      }

      # clone the generator (private fields except mcmc_ptr_)
      # then, re-init pointers
      new_gen <- HousepiGenerator$new()
      new_gen_priv <- new_gen$.__enclos_env__$private
      src_gen_priv <- private$gen_$.__enclos_env__$private
      for (field in names(new_gen_priv)) {
        val <- src_gen_priv[[field]]
        if (field != "mcmc_ptr_" && !is.function(val))
          new_gen_priv[[field]] <- helper_clone(val)
      }
      new_gen$init_pointers()
      obj_priv$gen_ <- new_gen

      # setup is not copied, build_full_() resupplies it at every run
      # share the data store, copy-on-write
      .engine$share_data_store(new_gen$mcmc_ptr, private$gen_$mcmc_ptr)

      # fresh pointer needs definitions pushed to the engine
      # cloned risk tables derive from the shared data, not re-derived
      obj_priv$push_definitions_(derive = FALSE)

      obj
    },

    #' @description
    #' Print the model and the state of its definitions.
    print = function() {
      header <- "Household transmission model {.val {self$name}}"

      # append link to template documentation if any
      if (!is.na(get_template_dir(private$gen_$raw_name)))
        header <- paste(header,
          "({.run [template](housepi::browse_templates(\"{self$name}\"))})"
        )
      header <- c(i = paste0(header, "."))

      # @todo rationalize if more options are added
      if (private$gen_$is_discrete_time)
        header <- c(header, "Model is set to {.emph discrete time} mode.")

      # each section is a self-contained cli_bullets block (headline + subs)
      blocks <- list()
      # data
      if (length(private$get_required_data_wide_()) > 0) blocks <- c(blocks,
        list(fmt_assert_cli("Data", private$assert_data_wide_(FALSE)))
      )
      # parameters
      if (length(private$get_required_parameters_()) > 0) blocks <- c(blocks,
        list(fmt_assert_cli("Parameters", private$assert_param_(FALSE)))
      )
      # latent
      if (length(private$get_required_latent_()) > 0) blocks <- c(blocks,
        list(fmt_assert_cli("Latent", private$assert_latent_(FALSE)))
      )
      # sources (reported only when one fails to resolve)
      unresolved <- private$assert_sources_(FALSE)
      if (length(unresolved)) blocks <- c(blocks,
        list(fmt_assert_cli("Sources", unresolved))
      )

      # inform about optional risks if dataset is suitable
      cols_extra <- setdiff(private$gen_$cols_data_wide,
                            data_reserved_cols)
      if (length(cols_extra)) {
        # group risks (optional: default when definitions are absent)
        risks <- private$get_required_group_risks_()
        if (length(risks)) {
          risks <- intersect(risks, template_missing_find(
            risks, names(private$group_risks_$defs)
          ))
          blocks <- c(blocks, list(fmt_assert_cli(
            "Group risks",
            if (length(risks)) c(
              "Using default group risks.",
              inline_cli_("Available: {.val {risks}}.")
            ) else NULL,
            minor = TRUE
          )))
        }
        # contact risks (optional)
        if (private$gen_$has_contact_risk) {
          blocks <- c(blocks, list(fmt_assert_cli(
            "Contact risks",
            if (!private$contact_risks_$assert()) "Using default contact risks."
            else NULL,
            minor = TRUE
          )))
        }
      }

      # infectivity profile (required by the model, or defined anyway)
      if (private$gen_$has_infec_profile ||
            length(private$infec_profile_$defs)) blocks <- c(blocks,
        list(fmt_assert_cli("Infectivity profile",
                            private$assert_infec_profile_(FALSE)))
      )

      # reporting process (only when attached)
      if (!all(vapply(private$reporting_process_, is.null, logical(1)))) {
        blocks <- c(blocks, list(fmt_assert_cli("Reporting process", NULL)))
      }

      # order blocks by gravity (error, minor, satisfied), stable within ties
      rank <- c(x = 1L, `!` = 2L, v = 3L)
      order_key <- vapply(blocks, \(b) rank[[names(b)[1]]], integer(1))
      blocks <- blocks[order(order_key)]
      order_key <- sort(order_key)

      # a single hint gathers all default (minor) sections after their block
      sms <- header
      minor <- which(order_key == 2L)
      last_minor <- if (length(minor)) max(minor) else 0L
      for (k in seq_along(blocks)) {
        sms <- c(sms, blocks[[k]])
        # aligned like a sub-line (invisible bullet), no clash with leading `i`
        if (k == last_minor)
          sms <- c(sms, " " = "A custom choice can be provided.")
      }

      print_bullets(sms)
      invisible(self)
    },

    #' @description
    #' Open the help page of the template model.
    help = function() {
      dir <- get_template_dir(private$gen_$raw_name)
      if (is.na(dir)) stop_h(
        "Model {.val {self$name}} is not a template and has no help page."
      )
      open_model_helper_(basename(dir))
      invisible(self)
    }
  )
)

# nocov start
# @todo create a standalone folder with the correct links to the headers
#' Export the definitions to the standalone folder, to debug the engine in C++.
#' @noRd
export_to_cpp_standalone <- function(hh_model, outpath = "inst/standalone/") {
  stopifnothousepi(hh_model)
  pvt <- hh_model$.__enclos_env__$private
  if (pvt$status_$is_defined) stop_h(
    "The standalone build reads no status process.",
    "Export a model without {.fn add_status_process}."
  )
  status_codes <- pvt$status_codes_()
  pvt$build_group_risks_()
  pvt$build_contact_risks_()

  gen <- pvt$gen_
  info_wide <- gen$info_data_wide
  # the definitions as the engine imports them, with their sources resolved
  latent_defs <- resolve_status_codes_(pvt$resolved_latent_(),
                                       pvt$status_dictionary_())
  utils::write.csv(latent_defs,
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
  utils::write.csv(pvt$group_risks_$expose_to_cpp(info_wide, status_codes),
                   file.path(outpath, "group_risks.csv"),
                   quote = FALSE, row.names = FALSE)
  utils::write.csv(pvt$resolved_profile_(),
                   file.path(outpath, "infec.csv"),
                   quote = FALSE, row.names = FALSE)
  utils::write.csv(pvt$resolved_param_(),
                   file.path(outpath, "param.csv"),
                   quote = FALSE, row.names = FALSE)
  utils::write.csv(char_setup_input_(hh_model$setup),
                   file.path(outpath, "setup.csv"),
                   quote = FALSE, row.names = FALSE)
}
# nocov end
