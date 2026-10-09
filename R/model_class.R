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

  ############################################################################
  ################################## PRIVATE #################################
  ############################################################################

  private = list(
    gen_ = NULL,
    param_ = NULL,
    augm_ = NULL,
    group_risks_ = NULL,
    contact_risks_ = NULL,
    infec_profile_ = NULL,
    reporting_process_ = NULL,
    collapse_seen_ = character(),

    ############################################################################
    ############################################################################
    ############################################################################

    # @description
    # Names of one source the latent declaration reads: latent, observation
    # and infectivity profile definitions. Read from the resolved tables, so a
    # reference counts under the kind it resolves to. Each
    # name is labelled with the family of the definition reading it, member-
    # and household-level observations alike.
    # @param type `character(1)`. Source letter.
    required_of_declaration_ = function(type) {
      defs <- private$resolved_latent_()
      family <- vapply(sub("^house$", "member", defs[["kind"]]), kind2family,
                       character(1), USE.NAMES = FALSE)
      by_kind <- lapply(unique(family), function(f) {
        find_required_definitions(defs[family == f, , drop = FALSE], type,
                                  c("value", "func", "model", "min", "max"), f)
      })
      c(unlist(by_kind),
        find_required_definitions(private$resolved_profile_(), type,
                                  c("func", "min", "max"),
                                  private$infec_profile_$family))
    },

    # @description
    # Latent variables required by implementation file and other definitions.
    get_required_latent_ = function() {
      c(private$gen_$latent_required,      # by the model
        private$required_of_declaration_("A"))   # by the latent declaration
    },

    # @description
    # Dataset headers required by implementation file and other definitions.
    get_required_data_wide_ = function() {
      c(private$gen_$data_required,        # by the model
        private$group_risks_$data_required,      # by group risks
        private$contact_risks_$data_required,    # by contact risks
        private$required_of_declaration_("D"))   # by the latent declaration
    },

    # @description
    # Group risks required by implementation file and other definitions.
    get_required_group_risks_ = function() {
      c(private$gen_$group_risks_required, # by the model
        private$required_of_declaration_("G"))   # by the latent declaration
    },

    # @description
    # Parameters required by implementation file and other definitions.
    get_required_parameters_ = function() {
      c(private$gen_$param_required,       # by the model
        find_required_definitions(private$resolved_param_(), "P",
                                  c("value", "func", "min", "max"),
                                  private$param_$family), # by param functions
        private$group_risks_$param_required,     # by group risks
        private$contact_risks_$param_required,   # by contact risks
        private$required_of_declaration_("P"))   # by the latent declaration
    },

    ############################################################################
    ############################################################################
    ############################################################################

    # @description
    # Whether household data is attached.
    has_data_wide_ = function() {
      length(private$gen_$cols_data_wide) > 0
    },

    # @description
    # Assert if latent definitions match the required.
    # @param wall `logical(1)` to raise error.
    assert_latent_ = function(wall = TRUE, warn_unused = TRUE) {
      private$augm_$assert(
        private$get_required_latent_(),
        wall = wall,
        warn_unused = warn_unused && !private$gen_$is_test_model()
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
      required <- private$get_required_data_wide_()

      # observations are generated during simulation: not required
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

    # @description
    # Assert group risk definitions from the model implementation script.
    assert_group_risks_ = function() {
      private$group_risks_$assert(private$get_required_group_risks_())
    },

    # @description
    # Assert the model supports a feature; definitions of it would be dropped.
    # @param supported `logical(1)` whether the model supports it.
    # @param label `character(1)` plural feature name.
    assert_feature_supported_ = function(supported, label) {
      if (supported) return(invisible(self))
      stop_h(
        "{label} are not supported for model {.val {self$name}}.",
        `!` = "Definitions will be dropped."
      )
    },

    # @description
    # Baselines an absolute risk anchors, keyed by parameter name: an absolute
    # risk leaves no scale for its baseline, which is held at the reference
    # value whether or not the parameter was declared. Each entry names the
    # risks that absorbed it.
    # @return `list` of `character()`
    baseline_pins_ = function() {
      if (is.null(private$gen_$mcmc_ptr)) return(list())
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

    # @description
    # Declared parameters, with every baseline an absolute risk anchors held
    # at the reference value.
    # @return `HousepiParameter`
    effective_param_ = function() {
      pins <- names(private$baseline_pins_())
      if (!length(pins)) return(private$param_)
      param <- private$param_$clone(deep = TRUE)
      for (name in pins) param$add_quo(
        name, value = 1, const = TRUE, .warn_over = FALSE
      )
      # a pin reads nothing, so leading keeps it ahead of whatever reads it
      is_pin <- param$defs[["name"]] %in% pins
      param$defs <- param$defs[order(!is_pin), , drop = FALSE]
      rownames(param$defs) <- NULL
      param
    },

    # @description
    # What a baseline was declared as, `NA` when it was not declared or
    # already is the reference value.
    # @param name `character(1)`. Baseline parameter name.
    # @return `character(1)` fingerprint of the whole definition row.
    declared_scale_ = function(name) {
      defs <- private$param_$defs
      row <- match(name, defs[["name"]])
      if (is.na(row)) return(NA_character_)
      is_reference <- is.na(defs[["prior"]][row]) &&
        is.na(defs[["func"]][row]) &&
        isTRUE(showcaser_extract_numeric(defs[["value"]][row]) == 1)
      if (is_reference) NA_character_
      else paste(unlist(defs[row, ]), collapse = "|")
    },

    # @description
    # Warn once per declaration for each baseline an absolute risk absorbed,
    # discarding what its declaration asked for. `collapse_seen_` holds what
    # was last reported: a baseline released by its risks is forgotten.
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
          defs[["value"]][match(name, defs[["name"]])]
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

    # @description
    # Read a table of parameter values a run supplies, refusing a value for a
    # baseline an absolute risk holds at the reference.
    # @param table `data.frame` or path. Values keyed by parameter name.
    # @param nrep `integer(1)`. Rows kept at most, evenly spaced, all if `NA`.
    # @param thin_ess `logical(1)`. Thin by effective sample size.
    # @param arg `character(1)`. Argument the table came from.
    read_param_table_ = function(table,
                                 nrep = NA_integer_,
                                 thin_ess = TRUE,
                                 arg = "param_table") {
      table <- helper_read_param_table_(table, thin_ess = thin_ess)
      private$assert_no_baseline_value_(table, arg)
      helper_thin_even_(table, nrep)
    },

    # @description
    # Refuse a value a run supplies for a baseline an absolute risk holds at
    # the reference. `NA` cells pass: they supply nothing.
    # @param given `list` or `data.frame`. Values keyed by parameter name.
    # @param arg `character(1)`. Argument the values came from.
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

    # @description
    # Assert if parameter definitions match required.
    # @param wall `logical(1)` to raise error.
    # @param warn_unused `logical(1)` warn if unused.
    assert_param_ = function(wall = TRUE, warn_unused = TRUE) {
      private$effective_param_()$assert(
        private$get_required_parameters_(),
        wall = wall,
        warn_unused = warn_unused && !private$gen_$is_test_model()
      )
    },

    # @description
    # Assert every reference names one declared definition. Scans the resolved
    # tables for terminals the dictionary left unsettled: a name two holders
    # claim, or a name nothing declares.
    # @param wall `logical(1)` to raise error.
    assert_sources_ = function(wall = TRUE) {
      # the full source-bearing column set, not the narrower requirement one:
      # a distribution or proposal argument carries references too
      cols <- c("value", "func", "prior", "model", "min", "max", "proposal")
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

    ############################################################################
    ############################################################################
    ############################################################################

    # @description
    # Assert an infectivity profile is defined when the model needs one.
    # @param wall `logical(1)` to raise error.
    assert_infec_profile_ = function(wall = TRUE) {
      if (!private$gen_$has_infec_profile ||
            length(private$infec_profile_$defs)) return(character())
      sms <- c("Missing infectivity profile definition.",
               "Define it with {.fn add_infectivity_profile}.")
      wall_or_return_(sms, wall)
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
      clash <- intersect(module$covariate_names, obs)
      if (length(clash)) stop_h(
        "A {module$family} is keyed on observation{?s} {.val {clash}}.",
        paste("An observation is generated by the model, which the risk feeds:",
              "it cannot stratify the risk it depends on."),
        "Key the {module$family} on a data column."
      )
      invisible(self)
    },

    # @description
    # Helper to assert and re-build group and contact risks.
    # `...` carries the stratification axes only the group risks know about.
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

    # @description
    # Assert and pre-build group risks.
    build_group_risks_ = function(wall = FALSE) {
      status_codes <- private$status_codes_()
      private$assert_group_risks_()
      private$helper_build_risks_(
        private$group_risks_, status_codes, wall,
        stratified = private$gen_$risk_status_stratified
      )
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
      data_factorized <- private$gen_$factorize_data_wide(data)
      # soft check: data may be intentionally incomplete
      # hard check at the wall
      private$assert_data_wide_(wall = FALSE)
      .engine$import_data_wide(private$gen_$mcmc_ptr, data_factorized)
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
      private$reporting_process_ <- list(wide = NULL, long = NULL)
      private$collapse_seen_ <- character()
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
        if (!is.null(private$param_$defs) && !nrow(private$param_$defs)) {
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
      data_lab <- private$gen_$status_data_labels
      if (length(data_lab)) return(data_lab)
      private$augm_$status_dictionary
    },

    # @description
    # Engine infected-status codes, named after the labels a risk may key on.
    # Labels win when the model knows any; a purely integer-coded status column
    # is its own dictionary, and the codes stand for themselves.
    status_codes_ = function() {
      dict <- private$status_dictionary_()
      if (length(dict)) return(stats::setNames(seq_along(dict), dict))
      codes <- as.integer(private$gen_$infec_status_set)
      stats::setNames(codes, as.character(codes))
    },

    # @description
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
        latent <- defs[["name"]][kind %in% c("latent", "status")]
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
      # a latent shadows the like-named column, as a simulated output fed
      # back as data: the latent is not an ambiguity
      sets$D <- setdiff(sets$D, sets$A)
      lapply(sets, function(nms) {
        nms <- unique(nms)
        nms[!is.na(nms) & nzchar(nms)]
      })
    },

    # @description
    # Source letter a reference resolves to, for every name the model
    # declares. A name two holders claim maps to `""`: no reference settles it.
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

    # @description
    # Parameter definitions with sources and `discr` resolved on demand.
    resolved_param_ = function() {
      resolve_discr_flags_(
        resolve_source_types_(private$effective_param_()$defs,
                              private$source_dictionary_()),
        private$gen_$param_discr_required,
        is_prior = TRUE
      )
    },

    # @description
    # Latent definitions with sources and `discr` resolved on demand.
    # Sources come first: a bare read of a data column inherits its
    # discreteness only once it carries the `D` letter, and a member-level
    # observation reads household values only once its sources are known.
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

    # @description
    # Infectivity profile definitions with sources resolved on demand.
    resolved_profile_ = function() {
      resolve_source_types_(private$infec_profile_$defs,
                            private$source_dictionary_())
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
      need_latent <- private$required_of_declaration_("A")
      need_param <- private$required_of_declaration_("P")
      all(need_latent %in% private$augm_$defs[["name"]]) &&
        all(need_param %in% private$effective_param_()$defs[["name"]]) &&
        private$augm_$covers_statuses(defs, private$gen_$infec_status_set)
    },

    # @description
    # Soft re-push of R-side definitions to the C++ engine, then report the
    # baselines the risks absorbed.
    # NOTE: Hard validation lives in `build_full_`.
    # @param wall `logical(1)` Whether to raise error.
    rebuild_ = function(wall = FALSE) {
      private$push_definitions_(wall)
      private$warn_baseline_collapse_()
      invisible(self)
    },

    # @description
    # Derive the risks from the data and push every definition the engine can
    # already take.
    # @param wall `logical(1)` Whether to raise error.
    # @param derive `logical(1)` Whether to re-derive the risk tables.
    push_definitions_ = function(wall = FALSE, derive = TRUE) {
      # guard: pre-init (before generator is created)
      if (is.null(private$gen_) ||
            is.null(private$gen_$mcmc_ptr)) {
        return(invisible(self))
      }
      # risks stratify on data: nothing to derive while none is attached
      if (!length(private$gen_$cols_data_wide)) {
        return(invisible(self))
      }

      # the status dictionary keys the risks that stratify on it, and the
      # latent status codes below: it must be settled before either reads it
      status_dict <- private$status_dictionary_()
      if (private$gen_$escaped %in% status_dict) stop_h(
        paste("Infected status {.val {private$gen_$escaped}} collides",
              "with the escaped-status label",
              "({.val {private$gen_$escaped}})."),
        "Rename the status or define a different {.code escaped =} label."
      )

      # risks derive from their own covariates (they append to param_$defs)
      if (derive && length(private$group_risks_$defs)) {
        private$build_group_risks_(wall)
      }
      if (derive && length(private$contact_risks_$defs)) {
        private$build_contact_risks_()
      }

      # defer the engine push while required data is incomplete
      if (length(private$assert_data_wide_(FALSE, simulating = TRUE)) > 0) {
        return(invisible(self))
      }
      # a forward reference waits for the decorator that declares it, rather
      # than pushing an unsettled source into the engine
      if (!wall && length(private$assert_sources_(FALSE)) > 0) {
        return(invisible(self))
      }

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

      # latent
      latent_defs <- resolve_status_codes_(private$resolved_latent_(),
                                           status_dict)
      # NOTE: unresolved definitions wait for completeness or for the wall
      if (wall || private$latent_declaration_complete_(latent_defs)) {
        .engine$import_latent(mcmc_ptr, latent_defs)
      }

      # status dictionary drives the factor refactoring of the status column
      private$gen_$set_status_levels(status_dict)

      invisible(self)
    },

    # @description
    # Full model build. Refresh runtime setup, enforce hard validation,
    # then push the final state to the engine.
    build_full_ = function(steps = 0L, adapt = 0L, burn = 0L,
                           thin = 1L, chains = 1L, threads = 1L,
                           id_draws = NULL, simulating = FALSE,
                           seed = draw_seed_(), ...) {
      # only run_inference resolves a selection; other builds store nothing
      if (!is.list(id_draws)) {
        id_draws <- list(id = character(), hh_id = character())
      }
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

    # Helper with simulation build.
    # FIXME: If there's a param table, build should change
    helper_simulation_build_ = function(param_table = NULL,
                                        nrep = NA_integer_,
                                        thin_ess = TRUE,
                                        ...) {
      nrep <- assert_count_(nrep, "Number of replicates", allow_na = TRUE)

      # read param table
      if (is.null(param_table)) {
        if (is.na(nrep)) nrep <- 1L
      } else {
        param_table <- private$read_param_table_(param_table, nrep, thin_ess)
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

      # internal build function
      private$build_full_(
        steps = nrep, simulating = TRUE, ...
      )

      # NOTE: there might be not parameters at all
      if (!nrow(self$parameters)) return(nrep)

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
      .engine$import_param_table(private$gen_$mcmc_ptr, param_table)

      nrep
    },

    # Set reporting process `slot` ("wide" or "long") to `func`, `NULL`
    # removing it.
    set_reporting_process_ = function(slot, func) {
      if (!is.null(func) && !is.function(func)) stop_h(
        "Reporting process {.arg {slot}} must be a function or {.code NULL}.",
        format_type_error(func)
      )
      private$reporting_process_[slot] <- list(func)
      invisible(self)
    },

    # Function relabelling each simulated `data_wide` and applying the wide
    # reporting process.
    # UPDATE: I shall drop necessity for refactoring with variant data.
    # `labelled = FALSE` hands the simulated frames over engine-coded when no
    # reporting process reads them, which only suits consumers reading `hh_id`
    # and `status`.
    helper_reporting_wide_ = function(labelled = TRUE) {
      gen <- private$gen_
      wide <- private$reporting_process_$wide
      if (is.function(wide)) \(x) wide(gen$refactor_data_wide(x))
      else if (labelled) \(x) gen$refactor_data_wide(x)
    },

    # @description
    # Run measure function on posterior simulations.
    # UPDATE: this has to be improved with loops in C++
    run_measure_ = function(param_table,
                            apply_wide,
                            apply_long = NULL,
                            nrep = NA_integer_,
                            thin_ess = TRUE,
                            threads = 1L,
                            progress = TRUE,
                            labelled = TRUE) {
      param_table <- private$read_param_table_(param_table, nrep, thin_ess)
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

  ############################################################################
  ################################## ACTIVE ##################################
  ############################################################################

  active = list(
    #' @field name `character(1)`
    #'  Model name.
    name = function(value) {
      if (!missing(value)) stop_read_only("name")
      private$gen_$name
    },

    #' @field escaped `character(1)`
    #' Label for status code `0` in a character or factor `status` column.
    escaped = function(value) {
      if (!missing(value)) stop_read_only("escaped")
      private$gen_$escaped
    },

    #' @field setup `list`
    #' Last setup used by the build process.
    setup = function(value) {
      if (!missing(value)) stop_read_only("setup")
      .engine$export_setup(private$gen_$mcmc_ptr)
    },

    #' @field data_wide `data.frame`
    #' Individual level data.
    #' One row per individual index by corresponding ID.
    data_wide = function(value) {
      if (!missing(value)) stop_read_only("data_wide")
      # no data attached yet (deferred until add_data_wide() or the wall)
      if (!length(private$gen_$cols_data_wide)) return(NULL)
      .engine$export_data_wide(private$gen_$mcmc_ptr) |>
        private$gen_$refactor_data_wide()
    },

    #' @field data_long `data.frame`
    #' Additional data rows indexed by individual (long format).
    data_long = function(value) {
      if (!missing(value)) stop_read_only("data_long")
      .engine$export_data_long(private$gen_$mcmc_ptr) |>
        private$gen_$refactor_data_long()
    },

    #' @field parameters `data.frame`
    #' Parameter definitions.
    parameters = function(value) {
      if (!missing(value)) stop_read_only("parameters")
      private$resolved_param_() |>
        housepi_defin_showcaser() |>
        sort_definitions_()
    },

    #' @field group_risks `data.frame`
    #' Group risk definitions.
    group_risks = function(value) {
      if (!missing(value)) stop_read_only("group_risks")
      private$group_risks_ |>
        helper_show_risk_table(setdiff(risk_reserved_cols, risk_pin_cols),
                               private$has_data_wide_()) |>
        sort_definitions_("group_risk", by_status = FALSE)
    },

    #' @field contact_risks `data.frame`
    #' Contact risk definitions.
    contact_risks = function(value) {
      if (!missing(value)) stop_read_only("contact_risks")
      private$contact_risks_ |>
        helper_show_risk_table(c("param_name", "n", "rule_group"),
                               private$has_data_wide_()) |>
        sort_definitions_("param_name", by_status = FALSE)
    },

    #' @field latent `data.frame`
    #' Latent definitions.
    latent = function(value) {
      if (!missing(value)) stop_read_only("latent")
      private$resolved_latent_() |>
        augm_defs_of_kind("latent") |>
        housepi_defin_showcaser() |>
        sort_definitions_()
    },

    #' @field observations `data.frame`
    #' Observation definitions.
    #' Column 'level' tells member ("member") from household ("house") entries.
    observations = function(value) {
      if (!missing(value)) stop_read_only("observations")
      resolved <- private$resolved_latent_()
      do.call(rbind, lapply(c("member", "house"), function(level) {
        out <- augm_defs_of_kind(resolved, level)
        if (is.null(out)) return(NULL)
        cbind(out, level = level)
      })) |>
        housepi_defin_showcaser() |>
        sort_definitions_()
    },

    #' @field status_process `data.frame`
    #' Status process: one row per infected status, with its code, label and
    #' probability.
    status_process = function(value) {
      if (!missing(value)) stop_read_only("status_process")
      defs <- augm_defs_of_kind(private$resolved_latent_(), "status")
      weights <- if (is.null(defs)) numeric() else
        unname(showcase_extract_args(defs[["model"]])[[1]])
      labels <- private$augm_$status_levels
      if (!length(labels)) labels <- rep(NA_character_, length(weights))
      data.frame(status = seq_along(weights),
                 label = labels,
                 weight = weights)
    },

    #' @field infectivity_profiles `data.frame`
    #' Infectivity profile definitions.
    infectivity_profiles = function(value) {
      if (!missing(value)) stop_read_only("infectivity_profiles")
      private$resolved_profile_() |>
        private$infec_profile_$showcase() |>
        sort_definitions_()
    },

    #' @field reporting_process `list`
    #' Reporting process: functions `wide` and `long`, `NULL` when absent.
    reporting_process = function(value) {
      if (!missing(value)) stop_read_only("reporting_process")
      private$reporting_process_
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
                      "infec_profile_", "reporting_process_",
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
        if (field != "mcmc_ptr_" && !is.function(val)) {
          new_gen_priv[[field]] <- helper_clone(val)
        }
      }
      new_gen$init_pointers()
      obj_priv$gen_ <- new_gen

      # @note setup is run-scoped, fully resupplied by build_full_ at every
      # run/evaluate, so it is intentionally not carried across a copy.
      # share data store (just copy-on-write)
      .engine$share_data_store(new_gen$mcmc_ptr, private$gen_$mcmc_ptr)

      # fresh pointer needs definitions pushed to the engine
      # @note cloned risk tables were derived on shared data so not re-derived
      obj_priv$push_definitions_(derive = FALSE)

      obj
    },

    #' @description
    #' Override the default printing behavior.
    print = function() {
      header <- "Household transmission model {.val {self$name}}"

      # append link to template documentation if any
      if (!is.na(get_template_dir(private$gen_$raw_name))) {
        header <- paste(header,
          "({.run [template](housepi::browse_templates(\"{self$name}\"))})"
        )
      }
      header <- c(i = paste0(header, "."))

      # UPDATE: to be rationalize in the future if more options are included
      if (private$gen_$is_discrete_time) {
        header <- c(header, "Model is set to {.emph discrete time} mode.")
      }

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
          is_default <- length(template_missing_find(
            risks, names(private$group_risks_$defs)
          )) > 0
          blocks <- c(blocks, list(fmt_assert_cli(
            "Group risks",
            if (is_default) c(
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
        if (k == last_minor) {
          sms <- c(sms, " " = "A custom choice can be provided.")
        }
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
#' @title Export to Cpp standalone (developer)
#' @description
#' Export definitions in standalone folder in order to run debuggin in Cpp.
#' Used for development.
#' @param hh_model `r rd_model_in_`
#' @details
#' This could be designed to be more helpful, e. g.
#' by creating a standalone folder somewhere with correct links to headers
#' @keywords internal
#' @noRd
export_to_cpp_standalone <- function(hh_model, outpath = "inst/standalone/") {
  stopifnothousepi(hh_model)
  pvt <- hh_model$.__enclos_env__$private
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
