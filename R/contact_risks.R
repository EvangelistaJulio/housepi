## housepi - Epidemiological Analysis of Household Outbreaks.
##
## Copyright © 2024-2026 Institut Pasteur (Paris), CNRS and INSERM.
## Distributed under the terms of the GNU General Public License (GPLv3).
## See the LICENSE file for details.
##

# Covariates read by the `from_` columns among `cols`.
pair_covs_ <- function(cols) {
  sub("^from_", "", grep("^from_", cols, value = TRUE))
}

# The `from_` then the `to_` columns of covariates `covs`.
pair_cols_ <- function(covs) {
  c(paste0("from_", covs), paste0("to_", covs))
}

#' Return logic table of contact risks (from/to paired columns).
#' @param data_wide Data from which to extract contact risks.
#' @param covariates Covariates as list of data column names.
#' @return data.frame with from_X, to_X columns and count n.
#' @keywords internal
#' @noRd
contact_risk_logic_table <- function(data_wide, covariates) {
  # convention: no hh_id implies households of size 1: no pairs possible
  if (!"hh_id" %in% colnames(data_wide)) return(NULL)

  # subset to needed column, and coerce covariates to character
  covariates <- as.character(unlist(covariates))
  data_wide <- data_wide[, c("hh_id", covariates), drop = FALSE]
  data_wide[] <- lapply(data_wide, as.character)

  # within a household every 'from' meets every 'to', so the ordered pairs of
  # two covariate keys are the product of their per-household counts. Only the
  # members meeting themselves have to come back off the diagonal
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

#' Count the observed pairs every contact rule matches.
#' @param tbl `data.frame` logic table with from_X and to_X columns.
#' @param data_wide Data the pairs are counted over.
#' @param covariates `character()` Covariates the rules key on.
#' @return `integer()` One count per row of `tbl`.
#' @keywords internal
#' @noRd
count_contact_rule_pairs <- function(tbl, data_wide, covariates) {
  pairs <- contact_risk_logic_table(data_wide, covariates)
  # no household structure, no pair to count
  if (!NROW(pairs)) return(rep(0L, nrow(tbl)))
  count_resolved_rules(tbl, pairs[pair_cols_(covariates)], pairs$n)
}

#' Associate parameter names to contact risk logic table.
#' @param this_table `data.frame` with from_X and to_X columns.
#' @param prefix Prefix for parameter names.
#' @param reference Named vector with from_/to_ prefixed names and values.
#' `".auto"` picks the pair with the largest sample size.
#' Empty leaves the risk without a reference.
#' @return `data.frame` with param_name column added.
#' @keywords internal
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
      # order convention: smaller key is always "from"
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

#' Assert that every `from_` column of a contact logic table has its `to_`
#' partner, and the reverse.
#' @param strata `data.frame` logic table with from_X and to_X columns.
#' @return `data.frame` the input, invisibly.
#' @keywords internal
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

#' Fill symmetric gaps in a custom contact risk logic table.
#' For each row, if the reverse from/to pair is missing,
#' add it with the same param_name.
#' @param tbl `data.frame` with from_X, to_Y, and param_name columns.
#' @return `data.frame` with missing reverse pairs added.
#' @keywords internal
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

#' Combine list of contact risk sub-tables into a single data.frame.
#' @param tables `list` of `data.frame`s.
#' @return `data.frame` with columns unified across sub-tables.
#' @keywords internal
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
  combined <- compute_rule_groups(combined)
  combined <- reorder_risks_logic_table(combined)
  rownames(combined) <- NULL
  as.data.frame(combined[c(setdiff(names(combined), tail_cols), tail_cols)])
}

# Validate and shape the automatic reference of contact risks.
# A reference naming a column of a cell (the pair of one covariate, or of all
# of them when joint) names every column of that cell.
auto_contact_reference_ <- function(by, reference, joint = FALSE) {
  if (identical(reference, ".auto")) return(".auto")

  valid_names <- pair_cols_(sort(by))
  # an unnamed reference has no way to match a covariate: names(NULL) %in% by
  # is vacuously TRUE, so this must be checked before the convenience branch
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
  # convenience: an un-prefixed name sets both ends of the pair
  bare <- names(reference) %in% by
  if (any(bare)) {
    reference <- c(reference[!bare],
                   stats::setNames(rep(reference[bare], 2L),
                                   pair_cols_(names(reference)[bare])))
  }
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

############################################################################
############################################################################
########################### Housepi Contact Risk ###########################
############################################################################
############################################################################

#' @title Internal R6 holder for housepi contact-risk definitions
#' @description
#' Helper class managing contact-risk definitions of the model object.
#' @keywords internal
#' @noRd
HousepiContactRisk <- R6::R6Class(
  classname = "HousepiContactRisk",
  cloneable = TRUE,

  ############################################################################
  ################################## PRIVATE #################################
  ############################################################################

  private = list(
    deep_clone = function(name, value) helper_clone(value),

    # Add risk definition in automatic mode.
    helper_auto_add_ = function(by, reference, prefix, joint, symmetric) {
      # validate covariates
      assert_valid_covariates(by)

      # load reference and validate reference names
      if (length(reference)) {
        reference <- auto_contact_reference_(by, reference, joint)
      }

      # store definition
      self$defs <- list(by = by, reference = reference, prefix = prefix,
                        joint = joint, symmetric = symmetric, mode = "auto")
      self$logic_tables <- list()
      self$param_overrides <- list()

      invisible(self)
    },

    # Add risk definition in custom mode.
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

      # cast covariates: joint wraps all together
      covariates <- as.list(by)
      if (isTRUE(self$defs[["joint"]])) covariates <- list(covariates)

      # default prefix
      prefix <- as.character(self$defs[["prefix"]])
      if (is.na(prefix)) prefix <- ""

      new_tables <- list()
      for (this_cov in covariates) {
        this_table <- contact_risk_logic_table(data_wide, this_cov)
        # a missing covariate cell is no level: its pairs match no rule
        cols <- pair_cols_(this_cov)
        if (!is.null(this_table))
          this_table <- this_table[stats::complete.cases(this_table[cols]), ,
                                   drop = FALSE]
        # extract reference entries for this sub-table, propagating .auto
        relevant_ref <- self$defs$reference
        if (!identical(relevant_ref, ".auto")) {
          relevant_ref <- relevant_ref[names(relevant_ref) %in% cols]
        }
        this_table <- associate_contact_param_names(
          this_table, prefix, relevant_ref, symmetric
        )
        new_tables[[length(new_tables) + 1]] <- this_table
      }
      self$logic_tables <- new_tables
      invisible(self)
    },

    # Stop if a rule group leaves ordered member pairs unmatched.
    assert_covered_ = function(data_wide) {
      covs <- pair_covs_(unlist(lapply(self$logic_tables, names)))
      pairs <- contact_risk_logic_table(data_wide, unique(covs))
      missed <- Reduce(`|`, lapply(self$logic_tables, function(tbl) {
        unmatched_cells_(tbl, pairs[risk_cov_cols_(tbl)])
      }))
      stop_unmatched_(sum(pairs$n[missed]), "Contact risk", "ordered pair",
                      self$defs)
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

        # a value the data never carries would encode to a catch-all: it is
        # refused before symmetrization, on the rules as written
        unknown <- unknown_strata_rows(
          this_table,
          known = lapply(stats::setNames(nm = this_cov),
                         \(cov) unique(data_wide[[cov]])),
          # both sides of a pair read the same data column
          cols = stats::setNames(rep(this_cov, 2L), pair_cols_(this_cov))
        )
        assert_no_dead_rules(this_table, unknown, "Contact risk",
                             "ordered pair")

        if (symmetric) this_table <- fill_symmetric_gaps(this_table)
        this_table$n <- count_contact_rule_pairs(this_table, data_wide,
                                                 this_cov)
        assert_no_dead_rules(this_table, this_table$n == 0, "Contact risk",
                             "ordered pair")
        self$logic_tables[[i]] <- this_table
      }
      invisible(self)
    }
  ),

  ############################################################################
  ################################## ACTIVE ##################################
  ############################################################################

  active = list(
    #' @field parameters `list`.
    #' Parameter definitions associated to contact risks.
    parameters = function(value) {
      risk_params_(self$logic_tables, self$param_overrides, self$raw_params,
                   self$is_absolute())
    },

    #' @field data_required `character()`.
    #' Data required by contact risks.
    data_required = function(value) {
      required <- self$defs[["by"]]
      if (is.null(required) || all(is.na(required))) return(character(0))
      required <- required[!is.na(required)]
      set_source_name_(required, self$family)
    },

    #' @field covariate_names `character()`. Data columns the risks stratify on,
    #' declared through `by` or carried by a custom logic table.
    covariate_names = function(value) {
      custom <- lapply(lapply(self$logic_tables, names), pair_covs_)
      covs <- as.character(unlist(c(self$defs[["by"]], custom)))
      unique(covs[!is.na(covs)])
    },

    #' @field param_required `character()`.
    #' Parameters required by contact risks.
    param_required = function(value) {
      required <- setdiff(unlist(lapply(self$logic_tables, `[[`, "param_name")),
                          ".reference")
      set_source_name_(as.character(required), self$family)
    },

    #' @field tidy_table `data.frame()`.
    #' Contact risk definitions table.
    tidy_table = function(value) {
      tidyup_contact_risk_tables(self$logic_tables)
    }
  ),

  ############################################################################
  ################################## PUBLIC ##################################
  ############################################################################

  public = list(
    #' @field family `character(1)`. Holder family.
    family = "contact risk",

    #' @field defs `list`. Current contact risk definition.
    defs = list(),

    #' @field raw_params `list`, `NULL` for the default.
    #' Raw parameter definitions.
    raw_params = NULL,

    #' @field logic_tables `list`.
    #' List of tables for contact risk pair definitions.
    logic_tables = list(),

    #' @field param_overrides `list`.
    #' Per-parameter MCMC overrides.
    param_overrides = list(),

    #' @description
    #' Whether the risk carries absolute values, some rule group having no
    #' reference.
    #' @return `logical(1)`
    is_absolute = function() {
      risk_is_absolute_(self$defs, self$logic_tables)
    },

    #' @inherit add_contact_risk
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

    ############################################################################
    ############################################################################
    ############################################################################

    #' @description
    #' Implement contact risk definition according to data.
    #' @param data_wide A data.frame with the dataset.
    implement = function(data_wide) {
      symmetric <- isTRUE(self$defs[["symmetric"]])

      if (self$defs[["mode"]] == "custom") {
        private$helper_implement_custom_(data_wide, symmetric)
      } else {
        private$helper_implement_auto_(data_wide, symmetric)
      }
      private$assert_covered_(data_wide)
      invisible(self)
    },

    #' @description
    #' Assert that contact risk definitions are compatible.
    assert = function() {
      length(self$defs) &&
        (!all(is.na(self$defs[["by"]])) || self$defs[["mode"]] == "custom")
    },

    #' @description
    #' Tidy table of contact risk definitions for Cpp.
    #' @param info list of lists with factor levels.
    expose_to_cpp = function(info) {
      defs <- self$tidy_table
      defs[risk_pin_cols] <- NULL
      # both ends of a pair encode against the levels of one data column
      encode_risk_levels_(defs, info,
                          key = \(col) sub("^(from|to)_", "", col),
                          cols = risk_cov_cols_(defs))
    }
  )
)
