## housepi - Epidemiological Analysis of Household Outbreaks.
## Copyright 2024-2026 Institut Pasteur (Paris), France.
## SPDX-License-Identifier: GPL-3.0-or-later
##

# Columns of a status-process table no rule keys on.
status_reserved_cols <- c("weight", "n", "rule_group")

# Data covariates a status-process rule keys on, `status` aside.
status_covs_ <- function(strata) {
  strata_cov_cols_(strata, c(status_reserved_cols, "status"))
}

# Number of infected levels a status column draws:
# * labels -> their count;
# * codes -> the highest.
status_n_levels_ <- function(status) {
  if (is.character(status)) return(length(unique(stats::na.omit(status))))
  max(0L, status, na.rm = TRUE)
}

#' Validate a status-process table and number its rule groups.
#' @param strata `data.frame` Columns `status`, `weight` and covariates.
#' @return `data.frame` The table, `status` as labels or integer codes,
#' `rule_group` appended.
#' @noRd
assert_status_strata_ <- function(strata) {
  assert_strata_shape_(
    strata, c("status", "weight"), "The status process table",
    "Give each rule its status in {.field status} and its {.field weight}."
  )
  # a table read back from get_status_process() has its counts recomputed
  strata[["n"]] <- NULL
  weight <- strata[["weight"]]
  is_weight <- is.numeric(weight) && !anyNA(weight) &&
    all(is.finite(weight) & weight >= 0)
  if (!is_weight) stop_h(
    "Status weights must be non-negative numbers.",
    "Give each weight as a finite number, 0 or more."
  )
  # labels, or infected-status codes counted from 1
  status <- strata[["status"]]
  if (is.factor(status)) status <- as.character(status)
  is_code <- is.numeric(status) &&
    all(is.na(status) | (status >= 1 & status == round(status)))
  if (!is.character(status) && !is_code) stop_h(
    "Column {.field status} must hold status labels or whole numbers from 1.",
    format_type_error(status)
  )
  strata[["status"]] <- if (is_code) as.integer(status) else status
  n_levels <- status_n_levels_(strata[["status"]])
  if (n_levels > max_status_level) stop_h(
    paste("Status process implies {n_levels} statuses, exceeding",
          "the maximum of {max_status_level}.")
  )

  rule_cols <- strata_cov_cols_(strata, status_reserved_cols)
  assert_valid_covariates(setdiff(rule_cols, "status"),
                          c(data_reserved_cols, risk_reserved_cols))
  dup <- which(duplicated(strata[rule_cols]))
  n_dup <- length(dup)
  if (n_dup) stop_h(
    "The status process repeats {cli::qty(n_dup)}rule{?s} at row{?s} {dup}.",
    "Keep one weight per combination of covariates and status."
  )
  strata <- compute_rule_groups(strata, cov_cols = rule_cols,
                                value_col = "weight")
  # a rule group blind to status scales every status alike
  unkeyed <- vapply(split(strata[["status"]], strata[["rule_group"]]),
                    \(s) all(is.na(s)), logical(1))
  if (any(unkeyed)) stop_h(
    "A rule group of the status process does not key on {.field status}.",
    "Its weights scale every status alike: give its rules a status."
  )
  strata
}

#' Assert a status-process table leaves no member's status undrawable,
#' at every infected level a member may draw:
#' * each rule group has exactly one best rule for it;
#' * every rule resolves some member;
#' * some status carries a positive weight.
#' @param strata `data.frame` Table validated by `assert_status_strata_()`.
#' @param levels `character()` Infected levels a member may draw.
#' @return `data.frame` Table, invisibly.
#' @noRd
assert_status_fit_ <- function(strata, data_wide, levels) {
  covs <- status_covs_(strata)
  assert_data_covariates(covs, colnames(data_wide))
  cells <- member_cells_(data_wide, covs, levels)
  win <- assert_strata_fit_(list(strata), cells, status_reserved_cols,
                            "weight", "The status process", "member")[[1]]

  # weight of every cell, then total per covariate combination
  weight <- apply(win, 1, \(w) prod(strata[["weight"]][w]))
  total <- tapply(weight, cells$.combo, sum)
  zero <- as.integer(names(total)[total == 0])
  n <- sum(cells$.weight[match(zero, cells$.combo)])
  if (n) stop_h(
    "The status process gives every status a zero weight for {n} member{?s}.",
    "Give some status a positive weight for them."
  )
  invisible(strata)
}

# Status table of probabilities applying to every member: labelled statuses,
# or codes in declaration order.
status_strata_from_weights_ <- function(weights) {
  if (!length(weights)) stop_h(
    "No status probability provided.",
    "Provide at least one infected-status probability."
  )
  labels <- rlang::names2(weights)
  if (any(nzchar(labels)) && !all(nzchar(labels))) stop_h(
    "Status labels were given to some probabilities only.",
    "Name every status probability, or none."
  )
  status <- if (all(nzchar(labels))) labels else seq_along(weights)
  data.frame(status = status, weight = status_weight_values_(weights))
}

# Evaluate status probabilities to numbers, naming the one that is not.
status_weight_values_ <- function(weights) {
  labels <- rlang::names2(weights)
  vapply(seq_along(weights), function(k) {
    w <- tryCatch(rlang::eval_tidy(weights[[k]]), error = function(e) NULL)
    # nolint next: object_usage_linter. read by the cli glue below
    label <- if (nzchar(labels[k])) labels[k] else k
    if (!is.numeric(w) || length(w) != 1L) stop_h(
      "Status probability {.val {label}} is not a number.",
      "Give each status probability as a number."
    )
    as.numeric(w)
  }, numeric(1))
}

########  HousepiStatusProcess  ########

#' Hold the rule table drawing the infected status of members in simulation.
#' @noRd
HousepiStatusProcess <- R6::R6Class(
  classname = "HousepiStatusProcess",
  cloneable = TRUE,

  ####  Private  ####

  private = list(
    deep_clone = function(name, value) helper_clone(value)
  ),

  ####  Active  ####

  active = list(
    # `logical(1)` Whether a status process is declared.
    is_defined = function(value) !is.null(self$strata),

    # `character()` Status labels, in order of first appearance.
    # * integer statuses, or no status process -> empty.
    labels = function(value) {
      status <- self$strata[["status"]]
      if (!is.character(status)) return(character())
      unique(status[!is.na(status)])
    },

    # `integer()` Integer statuses drawn, from 1 to the highest code.
    # * labelled statuses, or no status process -> empty.
    codes = function(value) {
      status <- self$strata[["status"]]
      if (!is.numeric(status)) return(integer())
      seq_len(status_n_levels_(status))
    }
  ),

  ####  Public  ####

  public = list(
    # `data.frame` Rule table: covariates, `status`, `weight` and `rule_group`.
    strata = NULL,

    # Record the rule table, labels ordering the statuses by first appearance.
    # @param strata `data.frame` Columns `status`, `weight` and covariates.
    set = function(strata) {
      self$strata <- assert_status_strata_(strata)
      invisible(self)
    },

    # Infected levels a member may draw:
    # * labelled statuses -> every label of the model dictionary;
    # * otherwise -> the integer statuses up to the highest.
    # @param status_codes `integer()` Engine status codes, named by label.
    draw_levels = function(status_codes) {
      if (is.character(self$strata[["status"]])) return(names(status_codes))
      as.character(self$codes)
    },

    # Assert no member's status is undrawable, see assert_status_fit_().
    # @param levels `character()` Infected levels a member may draw.
    assert_fit = function(data_wide, levels) {
      assert_status_fit_(self$strata, data_wide, levels)
      invisible(self)
    },

    # Read the rule table, one row per rule: its status and weight,
    # the `n` members it applies to, its rule group and its covariates.
    # * `data_wide` NULL -> `n` unknown.
    # @param levels `character()` Infected levels a member may draw.
    tidy_table = function(data_wide, levels) {
      if (!self$is_defined) return(data.frame(
        status = character(), weight = numeric(), n = integer(),
        rule_group = integer()
      ))
      covs <- status_covs_(self$strata)
      # a member counts once per rule, at whichever status the rule wins
      n <- if (is.null(data_wide)) NA_integer_
      else rule_counts_(self$strata, member_cells_(data_wide, covs, levels),
                        status_reserved_cols)
      cbind(self$strata[c("status", "weight")], n = n,
            self$strata[c("rule_group", covs)])
    },

    # Encode the rule table for the engine:
    # * covariates -> over the data levels;
    # * labelled statuses -> over the engine codes.
    # @param info `list` Class and levels of every data column.
    # @param status_codes `integer()` Engine status codes, named by label.
    expose_to_cpp = function(info, status_codes) {
      encode_strata_levels_(
        code_status_col_(self$strata, status_codes), info,
        cols = status_covs_(self$strata)
      )
    }
  )
)
