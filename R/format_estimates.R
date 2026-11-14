## housepi - Epidemiological Analysis of Household Outbreaks.
## Copyright 2024-2026 Institut Pasteur (Paris), France.
## SPDX-License-Identifier: GPL-3.0-or-later
##

# Format one estimate and its credible-interval bounds into `fmt`.
# Only the first `slots` values reach `fmt`.
format_est_num_ <- function(pe, low, up, fmt, as_integer, slots) {
  if (is.na(pe)) return(NA_character_)
  fill <- function(x) do.call(sprintf, c(fmt, as.list(x[seq_len(slots)])))
  if (isTRUE(pe == low && pe == up)) return(fill(rep(pe, 3)))
  if (isTRUE(pe < low) || isTRUE(pe > up)) warn_h(
    "Found an inconsistent credible interval: {pe}, CrI = [{low}, {up}]."
  )

  # decimals for 2 significant digits of the scale
  # * interval covering zero -> scale is the width;
  # * else -> smaller of estimate and 10 times the width, NA bound skipped;
  # * zero scale, from an inconsistent interval -> integers
  scale <- up - low
  digits <- 0L
  if (!isTRUE(low <= 0 && up >= 0))
    scale <- min(abs(pe), 10 * scale, na.rm = TRUE)
  if (isFALSE(as_integer) && scale > 0)
    digits <- as.integer(max(0, 1 - floor(log10(scale))))
  fill(sprintf("%.*f", digits, c(pe, low, up)))
}

# Format estimates and bounds element-wise, after checking lengths and `fmt`.
format_est_vector_ <- function(estimate, cri_lower, cri_upper,
                               fmt, as_integer, slots) {
  if (length(estimate) != length(cri_lower) ||
        length(estimate) != length(cri_upper)) stop_h(
    "Point estimate and CrI bounds must be numeric of the same length."
  )

  if (length(gregexpr("%s", fmt, fixed = TRUE)[[1]]) != slots) stop_h(
    "Invalid choice for character {.arg fmt}.",
    "It must contain exactly {slots} {.code %s} placeholder{?s}."
  )

  mapply(format_est_num_, estimate, cri_lower, cri_upper,
         MoreArgs = list(fmt = fmt, as_integer = as_integer, slots = slots),
         USE.NAMES = FALSE)
}

# Format tabular estimates and bounds column-wise, after checking dimensions.
format_est_matrix_ <- function(estimate, cri_lower, cri_upper,
                               fmt, as_integer, slots) {
  if (!all(dim(estimate) == dim(cri_lower),
           dim(estimate) == dim(cri_upper)))
    stop_h("Tabular inputs must have the same dimensions.")

  res <- as.data.frame(Map(
    format_est_vector_, as.data.frame(estimate),
    as.data.frame(cri_lower), as.data.frame(cri_upper),
    MoreArgs = list(fmt = fmt, as_integer = as_integer, slots = slots)
  ))
  names(res) <- names(estimate)
  res
}

# Append column `formatted` to a table of estimates, every column kept.
# * column `cri_lower` or `cri_upper` -> housepi table, point estimate from
#   column `predicted` or `estimate`;
# * else -> `summarise_draws()` table, point estimate and bounds from
#   the global summary.
format_est_table_ <- function(table, fmt, as_integer) {
  hint <- NULL
  if (any(c("cri_lower", "cri_upper") %in% names(table))) {
    pe_col <- intersect(c("predicted", "estimate"), names(table))
    if (length(pe_col) != 1L) stop_h(
      "The table must hold one point estimate column.",
      "Name it {.field predicted} or {.field estimate}, not both."
    )
    read <- c(pe_col, "cri_lower", "cri_upper")
  } else {
    # same naming as posterior::quantile2()
    probs <- unname(cri_level_probs(get_global_cri_level()))
    read <- c(get_global_point_estimate(), paste0("q", probs * 100))
    hint <-
      "Summarise the draws at the current {.fn set_housepi_summary} setting."
  }
  absent <- setdiff(read, names(table))
  if (length(absent))
    stop_h("The table has no column {.field {absent}}.", hint)
  table[["formatted"]] <- format_est_vector_(
    table[[read[1]]], table[[read[2]]], table[[read[3]]],
    fmt, as_integer, slots = 3
  )
  table
}

#' Format estimates and credible intervals into printable strings
#'
#' @description
#' Format point estimates, alone or with their credible interval, into strings.
#' A table of estimates, given without bounds, comes back whole
#' with the strings in a new column `formatted`:
#' * a summary from [summarise_draws()] is read at the point estimate
#'   and credible level set by [set_housepi_summary()];
#' * a table with columns `cri_lower` and `cri_upper`, such as the
#'   formatted measures of [check_predictive()] or the estimates of
#'   [check_performance()], is read at those bounds and at its point
#'   estimate column, `predicted` or `estimate`.
#'
#' Estimate and bounds are rounded to the same number of decimals,
#' enough to keep two significant digits of the estimate
#' and to tell the bounds apart
#' \insertCite{gelman_bayesian_2026}{housepi}.
#' When the interval spans zero, the decimals keep two significant digits
#' of the interval width instead.
#' @param estimate `numeric()` | `data.frame`
#' Point estimates, or a table of estimates.
#' @param cri_lower,cri_upper `numeric()` Credible interval bounds.
#' If `NULL`, estimates are formatted alone.
#' @param fmt `character(1)`
#' Sprintf-style format with one `%s` placeholder per value, filled with
#' the point estimate, then the lower and upper bounds.
#' If `NULL`, `"%s [%s, %s]"` is used with bounds and `"%s"` without.
#' @param as_integer `logical(1)` Whether values are formatted as integers.
#' @return `character()` of the same length as the inputs.
#' A matrix or data frame of estimates gives a data frame of strings,
#' and a table of estimates the same table with column `formatted` added.
#' @references
#' \insertAllCited{}
#' @examples
#' pe <- c(9.9, 7.69e-02, 1.13)
#' cri_low <- c(2.6, 6.77e-02, NA)
#' cri_up <- c(23.7, 0.087, 1.53)
#' format_estimates(pe, cri_low, cri_up)
#' # custom format
#' format_estimates(pe, cri_low, cri_up, fmt = "%s (95%% CrI, %s-%s)")
#' # point estimates alone
#' format_estimates(pe)
#'
#' hh_fit <- flu_1918_frederick |>
#'   new_housepi_model("final-size") |>
#'   add_parameter(import = "template") |>
#'   run_inference()
#' # posterior summary, as set by set_housepi_summary()
#' hh_fit |> summarise_draws() |> format_estimates()
#' @export
format_estimates <- function(estimate, cri_lower = NULL, cri_upper = NULL,
                             fmt = NULL, as_integer = FALSE) {
  if (xor(is.null(cri_lower), is.null(cri_upper))) stop_h(
    "Credible interval bounds must come in pairs.",
    "Give both {.arg cri_lower} and {.arg cri_upper}, or neither."
  )
  table <- is.null(cri_lower) && is.data.frame(estimate)
  alone <- is.null(cri_lower) && !table
  if (is.null(fmt)) fmt <- if (alone) "%s" else "%s [%s, %s]"
  if (table) return(format_est_table_(estimate, fmt, as_integer))

  # estimate alone -> missing bounds, the estimate sets the digits
  if (alone) cri_lower <- cri_upper <- replace(estimate, TRUE, NA_real_)
  tabular <- is.data.frame(estimate) || is.matrix(estimate)
  format_est <- if (tabular) format_est_matrix_ else format_est_vector_
  format_est(estimate, cri_lower, cri_upper, fmt, as_integer,
             slots = if (alone) 1 else 3)
}
