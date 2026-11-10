## housepi - Epidemiological Analysis of Household Outbreaks.
## Copyright 2024-2026 Institut Pasteur (Paris), France.
## SPDX-License-Identifier: GPL-3.0-or-later
##

#' Stop when a value is assigned to a read-only active field.
#' @param value The `value` argument of the active binding, missing on read.
#' @param active `character(1)` Name of the active method.
#' @noRd
assert_read_only_ <- function(value, active) {
  if (!missing(value)) stop_h("Active field {.field {active}} is read only.")
  invisible()
}

#' Stop on a non-empty message when `wall` is set, else return it.
#' @noRd
wall_or_return_ <- function(sms, wall) {
  if (isTRUE(wall) && length(sms)) stop_h(sms) else sms
}

#' Seeds of the engine random number generators, drawn from the R stream.
#' @param n `integer(1)` Number of seeds.
#' @noRd
draw_seed_ <- function(n = 1L) {
  sample.int(.Machine$integer.max, n)
}

#' A count argument as a whole positive number.
#' @param what `character(1)` What the count is, opening the error.
#' @param allow_na `logical(1)` Let a single `NA` through, as `NA_integer_`.
#' @noRd
assert_count_ <- function(x, what, allow_na = FALSE) {
  if (allow_na && length(x) == 1L && is.na(x)) return(NA_integer_)
  if (!rlang::is_scalar_integerish(x) || is.na(x) || x < 1) stop_h(
    "{what} should be a positive integer."
  )
  as.integer(x)
}

#' Split the `id_draws` selection into individual and household identifiers.
#' * `NULL` -> nothing stored;
#' * `NA` or empty -> the default household is followed.
#' Status filtering and the default household pick are left to the engine.
#' @param id_draws Selection given to [run_inference()].
#' @noRd
resolve_id_draws_ <- function(id_draws) {
  out <- list(id = character(), hh_id = character())
  # explicit opt-out: forward the sentinel so the engine stores nothing
  if (is.null(id_draws)) {
    out$id <- ".none"
    return(out)
  }
  # unspecified: let the engine follow one random household
  if (!length(id_draws) ||
        (!is.list(id_draws) && length(id_draws) == 1L && is.na(id_draws))) {
    return(out)
  }
  # named list: split individual and household identifiers
  if (is.list(id_draws)) {
    extra <- setdiff(names(id_draws), c("id", "hh_id"))
    if (length(extra)) stop_h(
      "Unknown {.arg id_draws} name{?s}: {.val {extra}}."
    )
    if (!is.null(id_draws$id)) out$id <- as.character(id_draws$id)
    if (!is.null(id_draws$hh_id))
      out$hh_id <- as.character(id_draws$hh_id)
    return(out)
  }
  # bare vector: individual identifiers
  out$id <- as.character(id_draws)
  out
}

# Lower-case name with `-`, `_`, `.` and runs of spaces folded to one space.
spaced_choice <- function(x) {
  trimws(tolower(gsub("\\s+", " ", gsub("[-_. ]", " ", x))))
}

#' Assert undocumented advanced settings overrides from `...`.
#' @param allowed `character()` Accepted argument names.
#' @param ... Advanced settings overrides.
#' @return `list` Non-null overrides keyed by engine settings key.
#' @noRd
assert_advanced_settings_ <- function(allowed, ...) {
  dots <- list(...)
  for (nm in names(dots)) {
    if (!nm %in% allowed) stop_h(
      "Unexpected argument {.arg {nm}} passed through {.code ...}.",
      find_closest_match(nm, setdiff(allowed, names(dots)))
    )
  }
  invisible()
}

#' Match `value` to one of `choices`, exactly.
#' * `value` identical to `choices` -> the first, as for an unset formal.
#' @param arg `character(1)` Argument name, as printed.
#' @param help `character(1)` Help topic documenting the choices.
#' @noRd
match_choice_ <- function(value, choices, arg, help) {
  if (identical(value, choices)) return(choices[[1]])
  if (rlang::is_string(value) && value %in% choices) return(value)
  pointer <- paste0("See {.run [?{help}](utils::help(\"{help}\", ",
                    "\"housepi\"))} for the accepted values.")
  stop_h(
    "Unknown {.arg {arg}} {.val {value}}.",
    if (rlang::is_string(value))
      find_closest_match(value, choices, otherwise = pointer)
    else pointer
  )
}

#' Format settings for the import_settings Rcpp method.
#' @noRd
char_settings_input_ <- function(defs) {
  # change scientific notation locally
  withr::local_options(list(scipen = 500))
  # correct id / house trace vectors
  for (key in intersect(c("id_draws", "house_draws"), names(defs))) {
    y <- unlist(defs[[key]])
    if (!length(y)) defs[[key]] <- ""
    else defs[[key]] <- paste0(unique(y), collapse = ":")
  }
  data.frame(lapply(defs, as.character))
}

#' Assert, then resolve the sampler settings of an inference run:
#' * `steps` -> extended to store the last thinned draw;
#' * sample size -> per chain and in total.
#' Counts are whole numbers, `thin` up to `steps`.
#' @noRd
resolve_settings_ <- function(steps, thin, adapt, chains) {
  # a whole number past the integer range stays a double
  assert_whole <- function(x, what, min = 1) {
    if (rlang::is_scalar_integerish(x) && !is.na(x) && x >= min) return()
    # nolint next: object_usage_linter. read by the cli glue below
    kind <- if (min > 0) "positive" else "non-negative"
    stop_h("Wrong settings: {what} should be a {kind} integer.")
  }
  assert_whole(steps, "number of steps")
  assert_whole(thin, "thinning interval")
  assert_whole(adapt, "number of adaptive steps", min = 0)
  assert_whole(chains, "number of chains")
  if (thin > steps) stop_h(
    "Wrong settings: {.arg thin} cannot be greater than {.arg steps}."
  )
  # round the sample size up, then extend steps to fill it
  chain <- steps %/% thin + (steps %% thin > 0)
  list(steps = chain * thin, chain = chain, total = as.double(chain) * chains)
}
