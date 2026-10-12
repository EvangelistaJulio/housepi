## housepi - Epidemiological Analysis of Household Outbreaks.
## Copyright 2024-2026 Institut Pasteur (Paris), France.
## SPDX-License-Identifier: GPL-3.0-or-later
##

# Covariates read by the `from_` columns among `cols`.
pair_covs_ <- function(cols) {
  sub("^from_", "", grep("^from_", cols, value = TRUE))
}

# The `from_` then the `to_` columns of covariates `covs`.
pair_cols_ <- function(covs) {
  c(paste0("from_", covs), paste0("to_", covs))
}

#' Build the logic table of contact risks, from/to paired columns.
#' @param covariates `list` Covariate data column names.
#' @return `data.frame` Columns `from_`, `to_` and the pair count `n`.
#' @noRd
contact_risk_logic_table <- function(data_wide, covariates) {
  # no hh_id -> households of size 1, no pairs
  if (!"hh_id" %in% colnames(data_wide)) return(NULL)

  # subset to needed columns, covariates coerced to character
  covariates <- as.character(unlist(covariates))
  data_wide <- data_wide[, c("hh_id", covariates), drop = FALSE]
  data_wide[] <- lapply(data_wide, as.character)

  # every 'from' meets every 'to' in a household:
  # * ordered pairs of two keys -> product of their per-household counts;
  # * self-pairs -> taken back off the diagonal.
  per_house <- table(data_wide$hh_id, row_keys_(data_wide, covariates))
  cnt <- crossprod(per_house)
  diag(cnt) <- diag(cnt) - colSums(per_house)

  # build the full observed grid
  observed <- unique(data_wide[, covariates, drop = FALSE])
  from_grid <- stats::setNames(observed, paste0("from_", covariates))
  to_grid <- stats::setNames(observed, paste0("to_", covariates))
  y <- merge(from_grid, to_grid, by = NULL)
  rownames(y) <- NULL

  # look up counts by the pair of covariate keys
  y$n <- as.integer(cnt[cbind(row_keys_(y, names(from_grid)),
                              row_keys_(y, names(to_grid)))])
  y$n[is.na(y$n)] <- 0L
  y[y$n > 0, , drop = FALSE]
}

#' Build the ordered member pairs a contact risk is checked against.
#' @param data_wide Data the pairs come from.
#' @param covs `character()` Covariates the rules key on.
#' @return `data.frame` Columns `from_` and `to_` of `covs`,
#' `.weight` the pairs each row holds, `.combo` its index.
#' @noRd
pair_cells_ <- function(data_wide, covs) {
  pairs <- contact_risk_logic_table(data_wide, covs)
  # no household structure -> no pair
  if (is.null(pairs)) pairs <- data.frame(
    lapply(stats::setNames(nm = pair_cols_(covs)), \(col) character()),
    n = integer()
  )
  pairs$.weight <- pairs$n
  pairs$.combo <- seq_len(nrow(pairs))
  pairs
}

#' Name the parameters of a contact risk logic table.
#' @param this_table `data.frame` Logic table with `from_` and `to_` columns.
#' @param prefix `character(1)` Prefix of parameter names.
#' @param reference Named vector with `from_`/`to_` prefixed names and values:
#' * `".auto"` -> the pair with the largest sample size;
#' * empty -> no reference.
#' @return `data.frame` Input with a `param_name` column.
#' @noRd
associate_contact_param_names <- function(this_table, prefix, reference,
                                          symmetric = TRUE) {
  tp <- "^to_"
  fp <- "^from_"
  colmask <- assess_covariate_values(this_table)
  from_cols <- grep(fp, colmask, value = TRUE)
  to_cols <- grep(tp, colmask, value = TRUE)
  covariates <- sub(fp, "", from_cols)

  # find reference index
  idx_refer <- NULL
  if (identical(reference, ".auto")) {
    idx_refer <- which.max(this_table$n)

  } else if (length(reference)) {
    idx_refer <- reference_rows_(this_table, colmask, reference)
    if (length(idx_refer) != 1) stop_h(
      "Provided reference contact risk is not valid.",
      `!` = paste(
        "The pair",
        "{reference_to_phrase(reference[grep(fp, names(reference))], fp)}",
        "{if (symmetric) '<->' else '-->'}",
        "{reference_to_phrase(reference[grep(tp, names(reference))], tp)}",
        "is never observed."
      )
    )
  }

  prefix <- format_risk_prefix_(prefix)

  # associate param names with one "from" and one "to" token
  this_table$param_name <- apply(
    this_table[colmask], 1, function(row) {
      ends <- c(paste(covariates, row[from_cols], sep = "_", collapse = "_"),
                paste(covariates, row[to_cols], sep = "_", collapse = "_"))
      # smaller key always "from"
      if (isTRUE(symmetric)) ends <- sort(ends)
      paste0(prefix, ends[1], "_to_", ends[2])
    }
  )
  if (length(idx_refer) == 1) {
    if (isTRUE(symmetric)) {
      ref_name <- this_table$param_name[idx_refer]
      this_table$param_name[this_table$param_name == ref_name] <- ".reference"
    } else {
      this_table$param_name[idx_refer] <- ".reference"
    }
  }

  this_table
}

#' Assert every `from_` column of a contact logic table has its `to_` partner,
#' and the reverse.
#' @return `data.frame` Input, invisibly.
#' @noRd
assert_paired_contact_cols <- function(strata) {
  cols <- risk_cov_cols_(strata)
  from <- pair_covs_(cols)
  to <- sub("^to_", "", grep("^to_", cols, value = TRUE))
  unpaired <- c(paste0("from_", setdiff(from, to), recycle0 = TRUE),
                paste0("to_", setdiff(to, from), recycle0 = TRUE))
  if (length(unpaired)) stop_h(
    "Contact strata column{?s} {.field {unpaired}} {?has/have} no partner.",
    paste("Add the {.code from_}/{.code to_} counterpart,",
          "set to {.code NA} to match any value.")
  )
  invisible(strata)
}

#' Add the reverse from/to pairs a custom contact logic table lacks,
#' with the same `param_name`.
#' @noRd
fill_symmetric_gaps <- function(tbl) {
  cov_cols <- risk_cov_cols_(tbl)
  nm <- colnames(tbl)
  flipped <- stats::setNames(tbl, ifelse(startsWith(nm, "from_"),
                                         sub("^from_", "to_", nm),
                                         sub("^to_", "from_", nm)))[nm]
  # a catch-all cell keys apart from every value
  key <- function(x) {
    cells <- lapply(x[cov_cols], \(v) ifelse(is.na(v), "\x1F", paste0("=", v)))
    do.call(paste, c(cells, sep = key_sep_))
  }
  gaps <- flipped[!key(flipped) %in% key(tbl), , drop = FALSE]
  if (!nrow(gaps)) return(tbl)
  tbl <- rbind(tbl, gaps)
  rownames(tbl) <- NULL
  tbl
}

#' Bind contact risk sub-tables into one `data.frame`, columns unified.
#' @noRd
tidyup_contact_risk_tables <- function(tables) {
  if (!length(tables)) return(data.frame())

  # covariates lead, the rule outcome (parameter, count, pins) closes the row
  present <- unique(unlist(lapply(tables, names)))
  tail_cols <- c("param_name", "n", "rule_group",
                 intersect(risk_pin_cols, present))
  all_cols <- c(sort(setdiff(present, risk_reserved_cols)),
                setdiff(tail_cols, "rule_group"))

  combined <- do.call(rbind, lapply(tables, pad_risk_cols_, cols = all_cols))
  rownames(combined) <- NULL
  combined <- risk_rule_groups_(combined)
  combined <- reorder_risks_logic_table(combined)
  rownames(combined) <- NULL
  as.data.frame(combined[c(setdiff(names(combined), tail_cols), tail_cols)])
}

# Validate and shape the automatic reference of contact risks.
# * reference naming one column of a cell -> names every column of that cell;
# * cell: the pair of one covariate, or of all covariates when `joint`.
auto_contact_reference_ <- function(by, reference, joint = FALSE) {
  if (identical(reference, ".auto")) return(".auto")

  valid_names <- pair_cols_(sort(by))
  # unnamed reference -> no covariate to match, checked before the convenience
  # branch, where names(NULL) %in% by is vacuously TRUE
  has_names <- !is.null(names(reference)) && all(nzchar(names(reference)))
  # a single covariate needs no name to disambiguate which one it labels
  if (!has_names && length(by) == 1 && length(reference) == 1) {
    names(reference) <- by
    has_names <- TRUE
  }
  if (!has_names) {
    if (length(reference) == length(by)) eg_values <- as.character(reference)
    else eg_values <- rep("value", length(by))
    # nolint next: object_usage_linter. read by the cli glue below
    pairs <- paste0(by, ' = "', eg_values, '"', collapse = ", ")
    stop_h(
      "Provided reference contact risk is not valid.",
      paste("{.arg reference} must be named using",
            "{cli::qty(by)}covariate{?s}:",
            "{.code reference = c({pairs})}.")
    )
  }
  # un-prefixed name -> sets both ends of the pair
  bare <- names(reference) %in% by
  if (any(bare))
    reference <- c(reference[!bare],
                   stats::setNames(rep(reference[bare], 2L),
                                   pair_cols_(names(reference)[bare])))

  wrong_entries <- setdiff(names(reference), valid_names)
  if (length(wrong_entries)) stop_h(
    "Unrecognized choice{?s} for reference: {.val {wrong_entries}}.",
    paste("Reference should match {.arg by}",
          "and use {.code from_}\\{.code to_} prefixes.")
  )
  cells <- if (isTRUE(joint)) list(by) else as.list(by)
  for (cell in cells) {
    cols <- pair_cols_(cell)
    named <- intersect(cols, names(reference))
    unnamed <- setdiff(cols, names(reference))
    if (length(named) && length(unnamed)) stop_h(
      "Reference names {.field {named}} but leaves {.field {unnamed}} unnamed.",
      "Name every column of the reference pair, or none of them."
    )
  }
  as.list(reference)
}

########  HousepiContactRisk  ########

#' Hold the contact-risk definition of a model.
#' @noRd
HousepiContactRisk <- R6::R6Class(
  classname = "HousepiContactRisk",
  cloneable = TRUE,

  ####  Private  ####

  private = list(
    deep_clone = function(name, value) helper_clone(value),

    # Add a definition in automatic mode.
    helper_auto_add_ = function(by, reference, prefix, joint, symmetric) {
      # validate covariates
      assert_risk_covariates_(by)

      # load reference and validate reference names
      if (length(reference))
        reference <- auto_contact_reference_(by, reference, joint)

      # store definition
      self$defs <- list(by = by, reference = reference, prefix = prefix,
                        joint = joint, symmetric = symmetric, mode = "auto")
      self$logic_tables <- list()
      self$param_overrides <- list()

      invisible(self)
    },

    # Add a definition in custom mode.
    helper_custom_add_ = function(strata, params, symmetric, reference) {
      assert_paired_contact_cols(strata)
      assert_valid_strata_(strata, params, reference, symmetric = symmetric)
      # cast covariate columns to character (consistent with auto mode)
      cov_cols <- risk_cov_cols_(strata)
      strata[cov_cols] <- lapply(strata[cov_cols], as.character)
      # store sentinel definition
      self$defs <- list(by = NA_character_, prefix = NA_character_,
                        reference = NULL, joint = TRUE,
                        symmetric = symmetric, mode = "custom")
      self$logic_tables <- list(strata)
      self$param_overrides <- if (is.null(params)) list() else params

      invisible(self)
    },

    # Build the logic tables of an automatic definition from the data.
    helper_implement_auto_ = function(data_wide, symmetric) {
      by <- self$defs[["by"]]
      assert_data_covariates(by, colnames(data_wide))

      # cast covariates, `joint` wraps them all together
      covariates <- as.list(by)
      if (isTRUE(self$defs[["joint"]])) covariates <- list(covariates)

      # default prefix
      prefix <- as.character(self$defs[["prefix"]])
      if (is.na(prefix)) prefix <- ""

      new_tables <- list()
      for (this_cov in covariates) {
        this_table <- contact_risk_logic_table(data_wide, this_cov)
        # missing covariate cell -> no level, its pairs match no rule
        cols <- pair_cols_(this_cov)
        if (!is.null(this_table))
          this_table <- this_table[stats::complete.cases(this_table[cols]), ,
                                   drop = FALSE]
        # extract reference entries for this sub-table, propagating .auto
        relevant_ref <- self$defs$reference
        if (!identical(relevant_ref, ".auto"))
          relevant_ref <- relevant_ref[names(relevant_ref) %in% cols]

        this_table <- associate_contact_param_names(
          this_table, prefix, relevant_ref, symmetric
        )
        new_tables[[length(new_tables) + 1]] <- this_table
      }
      self$logic_tables <- new_tables
      invisible(self)
    },

    # Assert each observed ordered pair resolves to one rule per group.
    assert_rules_fit_ = function(data_wide) {
      covs <- unique(pair_covs_(unlist(lapply(self$logic_tables, names))))
      assert_risk_fit_(self$logic_tables, pair_cells_(data_wide, covs),
                       "Contact risk", "ordered pair", self$defs)
      invisible(self)
    },

    # Validate custom mode logic table against data, then add counts.
    helper_implement_custom_ = function(data_wide, symmetric) {
      for (i in seq_along(self$logic_tables)) {
        # optionally fill symmetric gaps
        this_table <- self$logic_tables[[i]]
        # retrieve covariates
        this_cov <- pair_covs_(names(this_table))
        # validate that covariate columns exist in data
        assert_data_covariates(this_cov, colnames(data_wide))

        # value the data never carries -> error, before symmetrization,
        # as it would encode to a catch-all
        unknown <- unknown_strata_rows(
          this_table,
          known = lapply(stats::setNames(nm = this_cov),
                         \(cov) unique(data_wide[[cov]])),
          # both sides of a pair read the same data column
          cols = stats::setNames(rep(this_cov, 2L), pair_cols_(this_cov))
        )
        assert_no_dead_rules(this_table, unknown, risk_cov_cols_(this_table),
                             "param_name", "Contact risk", "ordered pair")

        if (symmetric) this_table <- fill_symmetric_gaps(this_table)
        this_table$n <- count_risk_rules_(this_table,
                                          pair_cells_(data_wide, this_cov))
        self$logic_tables[[i]] <- this_table
      }
      invisible(self)
    }
  ),

  ####  Active  ####

  active = list(
    # `list` Parameter definitions associated to contact risks.
    parameters = function(value) {
      risk_params_(self$logic_tables, self$param_overrides, self$raw_params,
                   self$is_absolute())
    },

    # `character()` Data columns required by contact risks.
    data_required = function(value) {
      required <- self$defs[["by"]]
      if (is.null(required) || all(is.na(required))) return(character(0))
      required <- required[!is.na(required)]
      set_source_name_(required, self$family)
    },

    # `character()` Data columns the risks stratify on,
    # declared through `by` or carried by a custom logic table.
    covariate_names = function(value) {
      custom <- lapply(lapply(self$logic_tables, names), pair_covs_)
      covs <- as.character(unlist(c(self$defs[["by"]], custom)))
      unique(covs[!is.na(covs)])
    },

    # `character()` Parameters required by contact risks.
    param_required = function(value) {
      required <- setdiff(unlist(lapply(self$logic_tables, `[[`, "param_name")),
                          ".reference")
      set_source_name_(as.character(required), self$family)
    },

    # `data.frame()` Contact risk definitions table.
    tidy_table = function(value) {
      tidyup_contact_risk_tables(self$logic_tables)
    }
  ),

  ####  Public  ####

  public = list(
    # `character(1)` Holder family.
    family = "contact risk",

    # `list` Current contact risk definition.
    defs = list(),

    # `list` | `NULL` Raw parameter definitions, `NULL` for the default.
    raw_params = NULL,

    # `list` Logic tables of the contact risk pairs.
    logic_tables = list(),

    # `list` Per-parameter MCMC overrides.
    param_overrides = list(),

    # Whether some rule group has no reference, making the risk absolute.
    is_absolute = function() {
      risk_is_absolute_(self$defs, self$logic_tables)
    },

    # Add the definition, arguments as in add_contact_risk().
    add_quo = function(by = NA_character_, ...,
                       reference = NULL,
                       joint = FALSE, symmetric = TRUE,
                       strata = NULL, params = NULL,
                       prior = NA_character_,
                       min = NA_real_, max = NA_real_,
                       proposal = NA_character_,
                       prefix = NA_character_,
                       .warn_over = TRUE) {
      # automatic definition
      by <- housepi_smart_unwrap(by, "by")
      prefix <- housepi_smart_unwrap(prefix, "prefix")
      reference <- housepi_smart_unwrap(reference, "reference", null_ok = TRUE)
      joint <- isTRUE(housepi_smart_unwrap(joint, "joint"))
      symmetric <- isTRUE(housepi_smart_unwrap(symmetric, "symmetric"))

      # custom definition
      strata <- housepi_smart_unwrap(strata, "strata", null_ok = TRUE)
      params <- housepi_smart_unwrap(params, "params", null_ok = TRUE)

      setup <- unwrap_risk_setup_(prior, proposal, min, max, list(...),
                                  self$family)
      mode <- risk_mode_(by, strata, params)
      if (.warn_over && length(self$defs)) warn_h(
        "Overwriting previous contact risk definition."
      )
      if (mode == "auto") {
        private$helper_auto_add_(by, reference, prefix, joint, symmetric)
      } else {
        private$helper_custom_add_(strata, params, symmetric, reference)
      }
      self$raw_params <- setup

      invisible(self)
    },

    ####  Build  ####

    # Build the logic tables from `data_wide`, then assert every pair fits.
    implement = function(data_wide) {
      symmetric <- isTRUE(self$defs[["symmetric"]])

      if (self$defs[["mode"]] == "custom") {
        private$helper_implement_custom_(data_wide, symmetric)
      } else {
        private$helper_implement_auto_(data_wide, symmetric)
      }
      private$assert_rules_fit_(data_wide)
      invisible(self)
    },

    # Whether a usable contact risk definition is held.
    assert = function() {
      length(self$defs) &&
        (!all(is.na(self$defs[["by"]])) || self$defs[["mode"]] == "custom")
    },

    # Encode the definition for the engine.
    # @param info `list` of `list`s Factor levels per data column.
    expose_to_cpp = function(info) {
      defs <- self$tidy_table
      defs[risk_pin_cols] <- NULL
      # both ends of a pair encode against the levels of one data column
      encode_strata_levels_(defs, info,
                            key = \(col) sub("^(from|to)_", "", col),
                            cols = risk_cov_cols_(defs))
    }
  )
)
