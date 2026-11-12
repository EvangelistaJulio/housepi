## housepi - Epidemiological Analysis of Household Outbreaks.
## Copyright 2024-2026 Institut Pasteur (Paris), France.
## SPDX-License-Identifier: GPL-3.0-or-later
##

########  HousepiModel  ########

# auxiliary environment
.housepi_env <- new.env(parent = emptyenv())
.housepi_env$hollow <- FALSE

#' Read a parameter table from draws or a data.frame, thinned to the ESS.
#' * named vector or list -> one row;
#' * value without a name, or parameter without exactly one value -> error.
#' @param thin_ess `logical(1)` Whether to thin to the ESS.
#' @param arg `character(1)` Argument the table came from.
#' @noRd
helper_read_param_table_ <- function(ptab, thin_ess = TRUE, arg = "draws") {
  # extract parameter draws from table
  if (inherits(ptab, "HousepiPosteriorDraws")) {
    ptab <- draws_block_(ptab, "param")
  }

  if (posterior::is_draws(ptab)) {
    ptab <- posterior::as_draws_matrix(ptab)
    # reduce sample to effective sample size (ESS)
    if (isTRUE(thin_ess)) {
      # an ESS above the draw count keeps every draw, capped or not
      # degenerate draws give NA, asserted below
      ess_bulk <- suppressWarnings(posterior::summarise_draws(
        ptab, ess_bulk = posterior::ess_bulk
      )$ess_bulk)
      # check no NA values in ESS
      if (anyNA(ess_bulk)) stop_h(
        "ESS values return NA.",
        "Check your draws for degenerate parameters."
      )
      # keep at most ESS draws, evenly spaced over the rows
      ptab <- helper_thin_even_(as.data.frame(ptab), floor(min(ess_bulk)))
    }
  }

  # named values: a single draw
  if (is_named_values_(ptab)) {
    if (!rlang::is_named(ptab)) stop_h(
      "{.arg {arg}} holds a value without a parameter name."
    )
    several <- names(ptab)[lengths(ptab) != 1L]
    if (length(several)) stop_h(
      # argument name pasted, as a second glue value breaks the plural
      paste0("Parameter{?s} {.val {several}} {?does/do} not hold ",
             "exactly one value in {.arg ", arg, "}."),
      "Pass several draws as a {.cls data.frame}, one row per draw."
    )
    return(list2DF(lapply(as.list(ptab), unname)))
  }

  tryCatch({
    return(as.data.frame(ptab))
  }, error = \(e) stop_h("Type of {.arg {arg}} not supported.")
  )
}

#' Whether each declared status is an integer code, not a label.
#' @param st `character()` Declared statuses, without `NA`.
#' @noRd
is_status_code_ <- function(st) {
  as_int <- suppressWarnings(as.integer(st))
  !is.na(as_int) & as.character(as_int) == st
}

#' Split declared statuses into integer codes and labels.
#' @param st `character()` Declared statuses, without `NA`.
#' @return `list` Elements `codes` and `labels`, `character()` each.
#' @noRd
split_statuses_ <- function(st) {
  st <- as.character(st)
  is_code <- is_status_code_(st)
  list(codes = st[is_code], labels = st[!is_code])
}

#' Resolve the character `status =` labels of latent definitions
#' to codes 1..K against the model status dictionary.
#' * integer code -> passes through;
#' * label absent from the dictionary -> error.
#' @param dictionary `character()` Ordered infected-status labels.
#' @return `data.frame` | `NULL` Definitions, `status` coded as integers.
#' @noRd
resolve_status_codes_ <- function(defs, dictionary) {
  if (is.null(defs)) return(defs)
  st <- as.character(unlist(defs[["status"]], use.names = FALSE))
  is_lab <- !is.na(st) & !is_status_code_(st)
  if (any(is_lab)) {
    code <- match(st[is_lab], dictionary)
    unknown <- unique(st[is_lab][is.na(code)])
    if (length(unknown)) stop_h(
      "Unknown status label{?s} {.val {unknown}}.",
      find_closest_match(
        unknown, dictionary,
        otherwise = "Use a status label of the data or of the status process."
      )
    )
    st[is_lab] <- as.character(code)
    defs[["status"]] <- st
  }
  defs
}

#' Resolve unsettled sources against the model's own declarations.
#' The distribution and proposal columns carry encoded arguments too,
#' so they are rewritten alongside the value-bearing ones.
#' @param dict `character()` Source letter, named by declared name.
#' @noRd
resolve_source_types_ <- function(defs, dict) {
  if (!length(defs)) return(defs)
  cols <- c("const", "func", "prior", "model", "min", "max", "proposal")
  for (col in intersect(names(defs), cols)) {
    if (is.character(defs[[col]]))
      defs[[col]] <- housepi_resolve_definition(defs[[col]], dict)
  }
  defs
}

#' Resolve unspecified `discr` declarations against model and data.
#' * definition neither carrying a model/prior nor sampled -> `discr` FALSE.
#' @param discr_req_defs `character()` Definitions the model requires discrete.
#' @param discr_col_types `logical()` Storage type per data column, by name.
#' @param discr_req_cols `character()` Data columns the model requires discrete.
#' @param is_prior `logical(1)` Whether prior, otherwise model.
#' @return `data.frame` Definitions with a logical `discr` column.
#' @noRd
resolve_discr_flags_ <- function(defs, discr_req_defs,
                                 discr_col_types = NULL,
                                 discr_req_cols = character(),
                                 is_prior = FALSE) {
  if (is.null(defs)) return(defs)
  discr <- as.logical(defs[["discr"]])
  unresolved <- is.na(discr)
  if (!any(unresolved)) return(defs)

  # whether the transmission model requires the definition's name discrete
  by_model <- defs[["name"]] %in% discr_req_defs

  # whether the datum column the definition reads is discrete
  by_data <- rep(FALSE, length(discr))
  if (!is.null(discr_col_types)) {
    # name of the data column a definition is pinned to
    col <- as.character(defs[["const"]])
    hit <- grepl("^[A-Za-z0-9_.]+\\|D$", col)
    # NA if it reads none
    col <- ifelse(hit, sub("\\|D$", "", col), NA_character_)
    # the column storage type, or the model requiring that column discrete
    by_data <- unname(discr_col_types[col]) | col %in% discr_req_cols
    by_data[is.na(by_data)] <- FALSE
  }

  distr_col <- if (is_prior) "prior" else "model"
  # sampled: a prior, a model, or bare (no constant, no function)
  has_distr <- !is.na(defs[[distr_col]]) |
    (is.na(defs[["const"]]) & is.na(defs[["func"]]))
  mask <- has_distr & (by_model | by_data)
  discr[unresolved] <- mask[unresolved]
  defs[["discr"]] <- discr
  defs
}

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
#' hh_model$name
#' @docType class
HousepiModel <- R6::R6Class(
  classname = "HousepiModel",
  cloneable = FALSE,
  lock_class = TRUE,

  ####  Private  ####

  private = list(
    # `list` The compiled transmission model, as compile_model_() returns it.
    compiled_ = NULL,

    # Rcpp pointer to the MCMC runner.
    mcmc_ptr_ = NULL,

    # `character(1)` Label of status code 0, outside declared infected statuses.
    escaped_ = "escaped",

    # `list` Engine coding of the attached data, `dataset` and `long_dataset`,
    # as cpp_data_factorize() returns it without the data.
    coding_ = list(
      dataset = list(cols = character(), status = integer(),
                     status_labels = character()),
      long_dataset = list(cols = character())
    ),

    # `character()` Hash of the factorized data attached, keyed as `coding_`.
    data_hash_ = c(dataset = NA_character_, long_dataset = NA_character_),

    param_ = NULL,
    augm_ = NULL,
    group_rates_ = NULL,
    contact_rates_ = NULL,
    status_ = NULL,
    infec_profile_ = NULL,
    reporting_process_ = NULL,
    collapse_seen_ = character(),
    # engine, and hashes of the data and definitions, at the last wall push
    fingerprint_ = NULL,

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

    # Names the transmission model requires, from its engine `key` getter.
    required_of_model_ = function(key) {
      set_source_name_(key(private$mcmc_ptr_), "model")
    },

    # Latent variables required by the model and by other definitions.
    get_required_latent_ = function() {
      c(
        # by the model
        private$required_of_model_(.engine$latent_keys),
        # by the latent declaration
        private$required_of_declaration_("A")
      )
    },

    # Data columns required by the model and by other definitions.
    get_required_dataset_ = function() {
      c(
        # by the model
        private$required_of_model_(.engine$data_keys),
        # by group rates
        private$group_rates_$data_required,
        # by contact rates
        private$contact_rates_$data_required,
        # by the latent declaration
        private$required_of_declaration_("D")
      )
    },

    # Group rates required by the model and by other definitions.
    get_required_group_rates_ = function() {
      c(
        # by the model
        private$required_of_model_(.engine$group_rate_keys),
        # by the latent declaration
        private$required_of_declaration_("G")
      )
    },

    # Parameters required by the model and by other definitions.
    get_required_parameters_ = function() {
      c(
        # by the model
        private$required_of_model_(.engine$param_keys),
        # by param functions
        find_required_definitions(private$resolved_param_(), "P",
                                  c("const", "func", "min", "max"),
                                  private$param_$family),
        # by group rates
        private$group_rates_$param_required,
        # by contact rates
        private$contact_rates_$param_required,
        # by the latent declaration
        private$required_of_declaration_("P")
      )
    },

    ####  Assertions  ####

    # Whether household data is attached.
    has_dataset_ = function() {
      length(private$coding_$dataset$cols) > 0
    },

    # Whether the model steps one whole day at once.
    is_discrete_time_ = function() {
      "_HOUSEPI_DISCRETE_TIME" %in% private$compiled_$flags
    },

    # Assert the required latents are defined.
    # * `wall` -> error, else the messages are returned.
    assert_latent_ = function(wall = TRUE, warn_unused = TRUE) {
      private$augm_$assert(
        private$get_required_latent_(),
        wall = wall,
        warn_unused = warn_unused && !is_test_model_(private$compiled_$name)
      )
    },

    # Assert dataset holds the columns required by the model and definitions.
    # * `wall` -> error, else the messages are returned;
    # * `simulating` -> generated observations are not required.
    assert_dataset_ = function(wall = TRUE, simulating = FALSE) {
      # no data attached -> error at the wall
      if (isTRUE(wall) && !private$has_dataset_()) stop_h(
        "No household data attached.",
        "Attach it with {.fn add_dataset}."
      )
      required <- private$get_required_dataset_()

      # simulating -> observations are generated, not required
      # subsetting keeps the source names the message reads
      if (isTRUE(simulating))
        required <- required[!required %in% private$augm_$observation_names]

      # a deterministic observation produces a data column
      cols <- private$coding_$dataset$cols
      det_obs <- private$augm_$deterministic_obs_names
      cols <- c(cols, intersect(det_obs, required))

      sms <- definition_missing_stop(required, cols, "data",
                                     noun = "column", style = "field")
      wall_or_return_(sms, wall)
    },

    # Assert the model reads every group rate defined.
    assert_group_rates_ = function() {
      private$group_rates_$assert(private$get_required_group_rates_())
    },

    # Assert the model supports a feature.
    # @param label `character(1)` Feature name, plural: "Contact rates".
    assert_feature_supported_ = function(supported, label) {
      if (supported) return(invisible(self))
      stop_h(
        "{label} are not supported for model {.val {self$name}}.",
        "Choose a model that supports them."
      )
    },

    # Pin baselines to the reference value, keyed by parameter name.
    # * absolute rate -> baseline held at the reference value, declared or not;
    # * baseline defined as a function -> not pinned, user's own scale.
    # @return `list` of `character()` Rates pinning each baseline.
    baseline_pins_ = function() {
      anchors <- .engine$rate_anchors(private$mcmc_ptr_)
      defs <- private$param_$defs
      pins <- list()
      for (name in names(anchors)) {
        # a baseline defined as a function carries the user's own scale
        row <- match(name, defs[["name"]])
        if (!is.na(row) && !is.na(defs[["func"]][row])) next
        rates <- anchors[[name]]$rates
        absolute <- rates[vapply(rates, private$group_rates_$is_absolute,
                                 logical(1))]
        if (isTRUE(anchors[[name]]$contacts) &&
              private$contact_rates_$is_absolute()) {
          absolute <- c(absolute, "contact")
        }
        if (length(absolute)) pins[[name]] <- absolute
      }
      pins
    },

    # Declared parameters, with the baselines absolute rates anchor
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

    # Warn once per declaration for each baseline an absolute rate absorbed,
    # discarding what its declaration asked for.
    # * `collapse_seen_` -> what was last reported;
    # * baseline released by its rates -> forgotten.
    warn_baseline_collapse_ = function() {
      pins <- private$baseline_pins_()
      defs <- private$param_$defs
      declared <- vapply(names(pins), private$declared_scale_, character(1))
      for (name in names(pins)) {
        if (is.na(declared[[name]]) ||
              identical(unname(private$collapse_seen_[name]),
                        declared[[name]])) next
        rates <- pins[[name]]
        value <- showcaser_extract_numeric(
          defs[["const"]][match(name, defs[["name"]])]
        )
        warn_h(
          paste("Absolute {cli::qty(rates)} rate{?s} {.val {rates}}",
                "{?has/have} no reference."),
          `!` = "Baseline {.val {name}} was turned into constant reference."
        )
      }
      private$collapse_seen_ <- declared
      invisible(self)
    },

    # Read a table of parameter values a run supplies.
    # * value for a baseline an absolute rate holds at the reference -> error.
    # @param table `data.frame` | path Values keyed by parameter name.
    # @param nrep `integer(1)` Rows kept at most, evenly spaced, all if `NA`.
    # @param thin_ess `logical(1)` Whether to thin by effective sample size.
    # @param arg `character(1)` Argument the table came from.
    read_param_table_ = function(table,
                                 nrep = NA_integer_,
                                 thin_ess = TRUE,
                                 arg = "draws") {
      table <- helper_read_param_table_(table, thin_ess = thin_ess,
                                        arg = arg)
      private$assert_no_baseline_value_(table, arg)
      helper_thin_even_(table, nrep)
    },

    # Assert no value a run supplies targets a baseline an absolute rate holds
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
        paste("An absolute rate leaves no scale for",
              "{cli::qty(held)}{?it/them}: fold the factor into the rate",
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
        warn_unused = warn_unused && !is_test_model_(private$compiled_$name)
      )
    },

    # Assert the status dictionary is well defined:
    # * integer codes in one source, labels in another, data included
    #   -> error, no label is tied to a code;
    # * labelled data missing a declared label -> error;
    # * escaped label among the infected statuses -> error.
    assert_status_ = function() {
      coding <- private$coding_$dataset
      data_lab <- coding$status_labels
      sources <- list(
        `the data` = list(codes = if (!length(data_lab)) coding$status,
                          labels = data_lab),
        `the status process` = list(codes = private$status_$codes,
                                    labels = private$status_$labels),
        `latent and observation definitions` = private$augm_$statuses,
        `group rates` = private$group_rates_$statuses
      )
      has <- \(part) vapply(sources, \(s) length(s[[part]]) > 0, logical(1))
      code_in <- has("codes")
      label_in <- has("labels")
      if (any(code_in) && any(label_in)) stop_h(
        "Statuses are given both as integer codes and as labels.",
        paste("Integer codes come from {names(which(code_in))},",
              "labels from {names(which(label_in))}.")
      )
      unknown <- setdiff(unlist(lapply(sources[-1], `[[`, "labels")), data_lab)
      if (length(data_lab) && length(unknown)) stop_h(
        "Status label{?s} {.val {unknown}} {?is/are} not in the data.",
        find_closest_match(unknown, data_lab),
        paste("To keep a status no member has, make data column",
              "{.field status} a factor with that status among its levels.")
      )
      if (private$escaped_ %in% private$status_dictionary_()) stop_h(
        paste("Infected status {.val {private$escaped_}} collides",
              "with the escaped-status label."),
        "Rename the status or define a different {.code escaped =} label."
      )
      invisible(self)
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
          sms <- c(sms, nothing, find_closest_match(nm, declared))
        }
      }
      wall_or_return_(sms, wall)
    },

    # Assert an infectivity profile is defined when the model needs one.
    # * `wall` -> error, else the messages are returned.
    assert_infec_profile_ = function(wall = TRUE) {
      if (!.engine$has_infec_profile(private$mcmc_ptr_) ||
            length(private$infec_profile_$defs)) return(character())
      sms <- c("Missing infectivity profile definition.",
               "Define it with {.fn add_infectivity_profile}.")
      wall_or_return_(sms, wall)
    },

    ####  Build  ####

    # Assert no rate is keyed on an observation column, at either level.
    # A rate feeds the transmission model the observation is generated from,
    # so its strata would depend on the order of declaration.
    # @param module `R6` Group or contact rate holder.
    assert_no_obs_rate_collision_ = function(module) {
      obs <- union(private$augm_$observation_names,
                   private$augm_$deterministic_obs_names)
      clash <- intersect(module$covariate_names, obs)
      if (length(clash)) stop_h(
        "A {module$family} is keyed on observation{?s} {.val {clash}}.",
        paste("An observation is generated by the model, which the rate feeds:",
              "it cannot stratify the rate it depends on."),
        "Key the {module$family} on a data column."
      )
      invisible(self)
    },

    # Assert and rebuild group or contact rates.
    # * parameter no longer generated -> dropped;
    # * new parameter -> appended at its default.
    # @param ... Stratification axes only the group rates know about.
    helper_build_rates_ = function(module, ...) {
      # proceed only if definitions exist
      if (!length(module$defs)) return(invisible(self))

      private$assert_no_obs_rate_collision_(module)

      before <- module$parameters
      module$implement(get_dataset(self), ...)
      after <- module$parameters
      private$drop_rate_params_(before[setdiff(names(before), names(after))],
                                module$family)
      added <- after[setdiff(names(after), private$param_$defs$name)]
      if (length(added)) private$param_$add_quo_list(
        convert_list_to_quosures(added), append = TRUE
      )
      invisible(self)
    },

    # Assert and pre-build group rates.
    build_group_rates_ = function(wall = FALSE) {
      status_codes <- private$status_codes_()
      private$assert_group_rates_()
      private$helper_build_rates_(
        private$group_rates_, status_codes, wall,
        stratified = .engine$rate_status_stratified(private$mcmc_ptr_)
      )
    },

    # Pre-build contact rates.
    build_contact_rates_ = function() {
      private$helper_build_rates_(private$contact_rates_)
    },

    # Factorize and push the dataset to the engine.
    # The data is validated, and the caller rebuilds definitions via rebuild_().
    set_dataset_ = function(data) {
      data <- private$embed_data_(data, "dataset")
      # soft check: data may be intentionally incomplete
      # hard check at the wall
      private$assert_dataset_(wall = FALSE)
      .engine$import_dataset(private$mcmc_ptr_, data$ids, data$cells)
      invisible(self)
    },

    # Validate, factorize and push the long dataset to the engine.
    set_long_dataset_ = function(data) {
      if (!is.data.frame(data)) stop_h(
        "Long dataset should be an object of class {.cls data.frame}.",
        format_type_error(data)
      )
      assert_valid_data_columns(colnames(data), "id")

      data <- private$embed_data_(data, "long_dataset")
      .engine$import_long_dataset(private$mcmc_ptr_, data$ids, data$cells)
      invisible(self)
    },

    # Factorize data for the engine, recording its coding and hash.
    # @param format `character(1)` `"dataset"` or `"long_dataset"`.
    # @return `list` Identifier and data columns, for the engine to import.
    embed_data_ = function(data, format) {
      embed <- cpp_data_factorize(data, private$escaped_)
      private$coding_[[format]] <- embed[names(embed) != "data"]
      private$data_hash_[[format]] <- rlang::hash(embed$data)
      embed$data
    },

    # Refactor a data.frame back to user classes and levels,
    # the status to the model dictionary.
    # @param format `character(1)` `"dataset"` or `"long_dataset"`.
    # @todo drop the refactoring of variant data
    refactor_data_ = function(x, format) {
      coding <- private$coding_[[format]]
      cpp_data_refactor(x, coding$info, coding$cols,
                        private$status_dictionary_(), private$escaped_,
                        produced = private$augm_$defs[["name"]])
    },

    # Drop rate parameters `derived`, no longer generated by the rates.
    # * set by add_parameter() -> warning;
    # * at its default -> silent.
    # @param derived `list` Default definition per parameter name.
    # @param family `character(1)` Rate family, named in the warning.
    drop_rate_params_ = function(derived, family) {
      if (!length(derived)) return(invisible(self))
      defaults <- HousepiParameter$new()$add_quo_list(
        convert_list_to_quosures(derived)
      )$defs
      # all-NA columns differ in type between holders
      row_of <- \(defs, p) vapply(defs[defs$name == p, ], as.character, "")
      defs <- private$param_$defs
      overridden <- Filter(
        \(p) !identical(row_of(defs, p), row_of(defaults, p)), names(derived)
      )
      for (p in names(derived)) private$param_$remove(p)
      if (length(overridden)) warn_h(
        "Definition{?s} of parameter{?s} {.val {overridden}} dropped.",
        `!` = paste("The attached data has no {family} stratum",
                    "for {cli::qty(overridden)}{?it/them}."),
        "List the remaining parameters with {.fn get_parameters}."
      )
      invisible(self)
    },

    # Drop the rate-derived parameters, re-derived after a rate redefinition.
    # Reads the pre-change logic tables, so it runs before the next rebuild_().
    drop_stale_rate_params_ = function() {
      stale <- c(names(private$group_rates_$parameters),
                 names(private$contact_rates_$parameters))
      for (p in stale) private$param_$remove(p)
      invisible(self)
    },

    # Fix rate tables across later data changes, counts recomputed.
    # Reference, strata and parameter names stay those of the attached data.
    freeze_rates_ = function() {
      private$group_rates_$frozen <- TRUE
      private$contact_rates_$frozen <- TRUE
      invisible(self)
    },

    # Reset every definition holder to an empty object.
    reset_definitions_ = function() {
      private$param_ <- HousepiParameter$new()
      private$augm_ <- HousepiAugmented$new()
      private$group_rates_ <- HousepiGroupRate$new()
      private$contact_rates_ <- HousepiContactRate$new()
      private$status_ <- HousepiStatusProcess$new()
      private$infec_profile_ <- HousepiInfecProfile$new()
      private$reporting_process_ <- list(dataset = NULL, long_dataset = NULL)
      private$collapse_seen_ <- character()
      invisible(self)
    },

    # Drop parameter definitions.
    drop_parameters_ = function() {
      private$param_ <- HousepiParameter$new()
      private$rebuild_()
      invisible(self)
    },

    # Drop group rates definitions, removing any associated parameters.
    drop_group_rates_ = function() {
      if (!is.null(private$group_rates_)) {
        for (p in names(private$group_rates_$parameters)) {
          private$param_$remove(p)
        }
        if (!is.null(private$param_$defs) && !nrow(private$param_$defs)) {
          private$drop_parameters_()
        }
      }
      private$group_rates_ <- HousepiGroupRate$new()
      private$rebuild_()
      invisible(self)
    },

    # Ordered infected-status dictionary of the model, in order of precedence:
    # * data labels, the data being already coded;
    # * status process order;
    # * character labels declared on latent definitions.
    # The last two serve data carrying no infected status, see assert_status_().
    status_dictionary_ = function() {
      data_lab <- private$coding_$dataset$status_labels
      if (length(data_lab)) return(data_lab)
      if (length(private$status_$labels)) return(private$status_$labels)
      private$augm_$statuses$labels
    },

    # Engine infected-status codes, named after the labels a rate may key on.
    # * model knows labels -> labels win;
    # * integer-coded statuses -> their own dictionary, the codes the data
    #   carries and the status process draws, each standing for itself.
    status_codes_ = function() {
      dict <- private$status_dictionary_()
      if (length(dict)) return(stats::setNames(seq_along(dict), dict))
      codes <- sort(union(as.integer(private$coding_$dataset$status),
                          private$status_$codes))
      stats::setNames(codes, as.character(codes))
    },

    # Declared names under each source, keyed by its source letter.
    source_sets_ = function() {
      ptr <- private$mcmc_ptr_
      live <- !is.null(ptr)
      # what the transmission model declares it needs, whether or not the
      # matching column or parameter has landed yet
      needs <- \(key) if (live) key(ptr) else character()

      defs <- private$augm_$defs
      latent <- character()
      obs <- character()
      if (!is.null(defs)) {
        kind <- defs[["kind"]]
        latent <- defs[["name"]][kind == "latent"]
        obs <- defs[["name"]][kind %in% c("member", "house")]
      }
      cols <- character()
      if (live)
        cols <- setdiff(private$coding_$dataset$cols, data_reserved_cols)

      sets <- list(
        # an observation stands for a column, at either level:
        # datum in inference, generated value in simulation
        D = c(cols, needs(.engine$data_keys), obs),
        G = c(names(private$group_rates_$defs), needs(.engine$group_rate_keys)),
        P = c(private$param_$defs[["name"]], needs(.engine$param_keys)),
        A = c(latent, needs(.engine$latent_keys))
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
        .engine$discrete_param_keys(private$mcmc_ptr_),
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
        .engine$discrete_latent_keys(private$mcmc_ptr_),
        discr_col_types = private$coding_$dataset$discr,
        discr_req_cols = .engine$discrete_data_keys(private$mcmc_ptr_)
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
    # then report the baselines the rates absorbed.
    # Hard validation runs in build_full_().
    # @param wall `logical(1)` Whether to raise an error.
    rebuild_ = function(wall = FALSE) {
      private$push_definitions_(wall)
      private$warn_baseline_collapse_()
      invisible(self)
    },

    # Derive the rates from the data,
    # then push every definition the engine can already take.
    # @param wall `logical(1)` Whether to raise an error.
    # @param derive `logical(1)` Whether to re-derive the rate tables.
    push_definitions_ = function(wall = FALSE, derive = TRUE) {
      # no data attached -> no rate to derive
      if (!private$has_dataset_()) return(invisible(self))

      # settle the status dictionary first
      private$assert_status_()

      # rates derive from their own covariates (they append to param_$defs)
      if (derive) {
        private$build_group_rates_(wall)
        private$build_contact_rates_()
      }

      # defer the engine push while required data is incomplete
      if (length(private$assert_dataset_(FALSE, simulating = TRUE)) > 0)
        return(invisible(self))

      # forward reference -> waits for the decorator that declares it,
      # never pushing an unsettled source
      if (!wall && length(private$assert_sources_(FALSE)) > 0)
        return(invisible(self))

      mcmc_ptr <- private$mcmc_ptr_
      payload <- private$payload_()
      .engine$import_param(mcmc_ptr, payload$parameters)
      .engine$import_group_rate(mcmc_ptr, payload$group_rates)
      .engine$import_contact_rate(mcmc_ptr, payload$contact_rates)
      .engine$import_infec_profile(mcmc_ptr, payload$infectivity_profile)
      private$push_status_process_(payload$status_process, derive)
      # unresolved definitions wait for completeness or for the wall
      if (wall || private$latent_declaration_complete_(payload$latent))
        .engine$import_latent(mcmc_ptr, payload$latent)
      if (wall) private$fingerprint_ <- c(
        list(model = private$compiled_$module), as.list(private$data_hash_),
        hash_payload_(payload)
      )

      invisible(self)
    },

    # Encode every definition for the engine, one entry per family.
    # Latent status labels are coded against the model dictionary.
    # Status process:
    # * labelled statuses -> read the model dictionary;
    # * codes -> stand for themselves;
    # * undefined -> `NULL`.
    # @return `list` Tables the engine imports.
    payload_ = function() {
      info_dataset <- private$coding_$dataset$info
      codes <- private$status_codes_()
      status <- private$status_
      list(
        parameters = private$resolved_param_(),
        group_rates = private$group_rates_$expose_to_cpp(info_dataset, codes),
        contact_rates = private$contact_rates_$expose_to_cpp(info_dataset),
        infectivity_profile = private$resolved_profile_(),
        status_process = if (status$is_defined) list(
          table = status$expose_to_cpp(info_dataset, codes),
          levels = private$status_draw_levels_()
        ),
        latent = resolve_status_codes_(private$resolved_latent_(),
                                       private$status_dictionary_())
      )
    },

    # Check the status-process table against the data, then push it.
    # * `NULL` -> nothing pushed.
    # @param status `list` Status-process entry of payload_().
    # @param derive `logical(1)` Whether to check the table against the data.
    push_status_process_ = function(status, derive = TRUE) {
      if (is.null(status)) return(invisible(self))
      if (derive) private$status_$assert_fit(get_dataset(self), status$levels)
      .engine$import_status_process(private$mcmc_ptr_, status$table,
                                    length(status$levels))
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
      if (!private$status_$is_defined) return(private$coding_$dataset$status)
      seq_along(private$status_draw_levels_())
    },

    # Build the model: refresh the run settings, enforce hard validation,
    # then push the final state to the engine.
    build_full_ = function(steps = 0L, adapt = 0L,
                           thin = 1L, chains = 1L, threads = 1L,
                           id_draws = NULL, simulating = FALSE, nrep = NULL,
                           seed = draw_seed_(), ...) {
      # only run_inference resolves a selection, other builds store nothing
      if (!is.list(id_draws))
        id_draws <- list(id = character(), hh_id = character())

      .engine$import_settings(
        private$mcmc_ptr_,
        char_settings_input_(c(list(
          # a simulation runs one engine step per replicate
          steps = if (is.null(nrep)) steps else nrep,
          adapt = adapt,
          thin = thin,
          chains = chains,
          threads = threads,
          id_draws = id_draws$id,
          house_draws = id_draws$hh_id,
          seed = seed
        ), list(...)))
      )
      # hard validation before MCMC
      private$assert_dataset_(wall = TRUE, simulating = simulating)
      # a column not yet attached is the likeliest cause of an unsettled source
      private$assert_sources_(wall = TRUE)
      private$assert_param_(wall = TRUE)
      private$assert_group_rates_()
      private$assert_infec_profile_(wall = TRUE)
      # engine build reports unresolved read first
      private$rebuild_(wall = TRUE)
      # latent assertion reports plainly missing definitions
      private$assert_latent_(wall = TRUE)
      invisible(self)
    },

    # Stop on a model data key that simulation reads as a calendar date,
    # while transmission seeds its latents at day 0.
    # * key without a definition -> class of its data column;
    # * key defined by a definition -> data it reads,
    #   transmission latents excepted, simulation draws them.
    # @todo drop once an introduction process sets day 0 on the data calendar
    assert_simulation_days_ = function() {
      transm <- .engine$transmission_latent_keys(private$mcmc_ptr_)
      if (!length(transm)) return(invisible())
      defs <- private$resolved_latent_()
      defs <- defs[!defs[["name"]] %in% transm, , drop = FALSE]
      info <- private$coding_$dataset$info
      cols <- c("const", "func", "model", "min", "max")
      date_read_by <- function(name) {
        row <- defs[defs[["name"]] == name, , drop = FALSE]
        if (!nrow(row))
          return(if (identical(info[[name]]$class, "Date")) name)
        # a model reads its own draw
        reads <- setdiff(c(find_required_definitions(row, "D", cols, NA),
                           find_required_definitions(row, "A", cols, NA)),
                         name)
        for (x in reads) if (length(col <- date_read_by(x))) return(col)
      }
      for (key in .engine$data_keys(private$mcmc_ptr_)) {
        col <- date_read_by(key)
        if (!length(col)) next
        stop_h(
          paste0("Simulation does not support a calendar date as ",
                 "{.field {key}}",
                 if (col != key) ", read from data column {.field {col}}", "."),
          paste0("Give {.field {key}} in days since the first infection ",
                 "in each household.")
        )
      }
      invisible()
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
        if (!nrow(get_parameters(self))) {
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
      generated <- setdiff(
        .engine$transmission_latent_keys(private$mcmc_ptr_),
        private$augm_$defs[["name"]]
      )
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
      private$build_full_(nrep = nrep, simulating = TRUE, ...)
      private$assert_simulation_days_()

      # there may be no parameters at all
      if (!nrow(get_parameters(self))) return(nrep)

      # no column: a parameter is drawn from its prior
      if (is.null(param_table))
        param_table <- data.frame(row.names = seq_len(nrep))
      # attach table with parameters
      .engine$import_param_table(private$mcmc_ptr_, param_table)

      nrep
    },

    # Set reporting process `slot`, "dataset" or "long_dataset", to `func`.
    # * `func` NULL -> removes it.
    set_reporting_process_ = function(slot, func) {
      if (!is.null(func) && !is.function(func)) stop_h(
        "Reporting process {.arg {slot}} must be a function or {.code NULL}.",
        format_type_error(func)
      )
      private$reporting_process_[slot] <- list(func)
      invisible(self)
    },

    # Build the function relabelling each simulated `dataset`,
    # then applying the dataset reporting process.
    # * not `labelled`, no reporting process -> frames stay engine-coded,
    #   for consumers reading `hh_id` and `status` only.
    # @todo drop the relabelling of variant data
    helper_reporting_dataset_ = function(labelled = TRUE) {
      report <- private$reporting_process_$dataset
      if (is.function(report)) \(x) report(private$refactor_data_(x, "dataset"))
      else if (labelled) \(x) private$refactor_data_(x, "dataset")
    },

    # Simulate replicates from `draws`, each through the reporting process.
    # * `draws` NULL -> parameters drawn from their prior;
    # * `apply_dataset` NULL -> datasets returned;
    # * `apply_long_dataset` NULL, or no long dataset reporting process
    #   -> no long dataset output.
    # @param ... Run settings, passed to the build.
    # @return `list` Entries `dataset` and `long_dataset`, one per replicate,
    # and `parameters`, one row per replicate.
    # @todo loop in C++
    simulate_ = function(draws = NULL,
                         nrep = NA_integer_,
                         thin_ess = TRUE,
                         apply_dataset = NULL,
                         apply_long_dataset = NULL,
                         labelled = TRUE,
                         progress = FALSE,
                         ...) {
      private$helper_simulation_build_(
        draws = draws, nrep = nrep, thin_ess = thin_ess, ...
      )
      # load a measure on the long dataset
      report_long <- private$reporting_process_$long_dataset
      apply_long_dataset_func <- NULL
      if (is.function(report_long) && is.function(apply_long_dataset))
        apply_long_dataset_func <- \(x) apply_long_dataset(report_long(x))

      # run simulations
      output <- with_progress_(progress, "simulate", 1, function(bar) {
        .engine$run_simulation(private$mcmc_ptr_,
                               private$helper_reporting_dataset_(labelled),
                               apply_dataset, apply_long_dataset_func, bar)
      })
      if (!length(output$long_dataset)) output$long_dataset <- NULL
      output
    }
  ),

  ####  Active  ####

  active = list(
    #' @field name `character(1)` Model name.
    name = function(value) {
      assert_read_only_(value, "name")
      spaced_choice(private$compiled_$name)
    },

    #' @field escaped `character(1)`
    #' Label for status code `0` in a character or factor `status` column.
    escaped = function(value) {
      assert_read_only_(value, "escaped")
      private$escaped_
    },

    #' @field settings `list` Last settings used by the build process.
    settings = function(value) {
      assert_read_only_(value, "settings")
      .engine$export_settings(private$mcmc_ptr_)
    }
  ),

  ####  Public  ####

  public = list(
    #' @description
    #' Create a `HousepiModel`, the object inference and simulation run on.
    #' @param data `r rd_dataset_`
    #' If omitted, attach it later with [add_dataset()].
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

      # check the data, NULL allowed and attached later with add_dataset()
      data <- assert_household_data(data, model)
      if (length(escaped) != 1 || !is.character(escaped)) stop_h(
        "{.arg escaped} must be a non-empty {.cls character}.",
        format_type_error(escaped)
      )
      private$escaped_ <- escaped

      # initialize internal fields
      private$reset_definitions_()
      # compile the transmission model, then create its runner
      private$compiled_ <- compile_model_(
        model = model, as_template = as_template,
        force = force, quiet = quiet, ...
      )
      private$mcmc_ptr_ <- create_mcmc_runner_(private$compiled_)
      # import data, if provided (no user definitions to rebuild yet)
      if (!is.null(data)) private$set_dataset_(data)

      invisible(self)
    },

    #' @description Deep copy this object.
    #' @return A new independent object.
    copy = function() {
      .housepi_env$hollow <- TRUE
      on.exit(.housepi_env$hollow <- FALSE)
      obj <- HousepiModel$new()
      obj_priv <- obj$.__enclos_env__$private

      # clone every field but the runner, created afresh
      for (field in c("compiled_", "param_", "augm_", "group_rates_",
                      "contact_rates_", "status_", "infec_profile_",
                      "reporting_process_", "collapse_seen_", "escaped_",
                      "coding_", "data_hash_")) {
        obj_priv[[field]] <- helper_clone(private[[field]])
      }
      obj_priv$mcmc_ptr_ <- create_mcmc_runner_(obj_priv$compiled_)

      # settings are not copied, build_full_() resupplies them at every run
      # share the data store, copy-on-write
      .engine$share_data_store(obj_priv$mcmc_ptr_, private$mcmc_ptr_)

      # fresh pointer needs definitions pushed to the engine
      # cloned rate tables derive from the shared data, not re-derived
      obj_priv$push_definitions_(derive = FALSE)

      obj
    },

    #' @description
    #' Print the model and the state of its definitions.
    print = function() {
      header <- "Household transmission model {.val {self$name}}"

      # append link to template documentation if any
      if (!is.na(get_template_dir(private$compiled_$name)))
        header <- paste(header,
          "({.run [template](housepi::browse_templates(\"{self$name}\"))})"
        )
      header <- c(i = paste0(header, "."))

      # @todo rationalize if more options are added
      if (private$is_discrete_time_())
        header <- c(header, "Model is set to {.emph discrete time} mode.")

      # each section is a self-contained cli_bullets block (headline + subs)
      blocks <- list()
      # data
      if (length(private$get_required_dataset_()) > 0) blocks <- c(blocks,
        list(fmt_assert_cli("Data", private$assert_dataset_(FALSE)))
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

      # inform about optional rates if dataset is suitable
      cols_extra <- setdiff(private$coding_$dataset$cols,
                            data_reserved_cols)
      if (length(cols_extra)) {
        # group rates (optional: default when definitions are absent)
        rates <- private$get_required_group_rates_()
        if (length(rates)) {
          rates <- intersect(rates, definition_missing_find(
            rates, names(private$group_rates_$defs)
          ))
          blocks <- c(blocks, list(fmt_assert_cli(
            "Group rates",
            if (length(rates)) c(
              "Using default group rates.",
              inline_cli_("Available: {.val {rates}}.")
            ) else NULL,
            minor = TRUE
          )))
        }
        # contact rates (optional)
        if (.engine$has_contact_rate(private$mcmc_ptr_)) {
          blocks <- c(blocks, list(fmt_assert_cli(
            "Contact rates",
            if (!private$contact_rates_$assert()) "Using default contact rates."
            else NULL,
            minor = TRUE
          )))
        }
      }

      # infectivity profile (required by the model, or defined anyway)
      if (.engine$has_infec_profile(private$mcmc_ptr_) ||
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
      dir <- get_template_dir(private$compiled_$name)
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
  pvt$build_group_rates_()
  pvt$build_contact_rates_()

  # the definitions as the engine imports them, with their sources resolved
  payload <- pvt$payload_()
  utils::write.csv(payload$latent,
                   file.path(outpath, "latent.csv"),
                   quote = FALSE, row.names = FALSE)
  utils::write.csv(payload$contact_rates,
                   file.path(outpath, "contact_rates.csv"),
                   quote = FALSE, row.names = FALSE)
  utils::write.csv(.engine$export_dataset(pvt$mcmc_ptr_),
                   file.path(outpath, "dataset.csv"),
                   quote = FALSE, row.names = FALSE)
  utils::write.csv(.engine$export_long_dataset(pvt$mcmc_ptr_),
                   file.path(outpath, "long_dataset.csv"),
                   quote = FALSE, row.names = FALSE)
  utils::write.csv(payload$group_rates,
                   file.path(outpath, "group_rates.csv"),
                   quote = FALSE, row.names = FALSE)
  utils::write.csv(payload$infectivity_profile,
                   file.path(outpath, "infec.csv"),
                   quote = FALSE, row.names = FALSE)
  utils::write.csv(payload$parameters,
                   file.path(outpath, "param.csv"),
                   quote = FALSE, row.names = FALSE)
  utils::write.csv(char_settings_input_(hh_model$settings),
                   file.path(outpath, "settings.csv"),
                   quote = FALSE, row.names = FALSE)
}
# nocov end
