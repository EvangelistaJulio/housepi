## housepi - Epidemiological Analysis of Household Outbreaks.
## Copyright 2024-2026 Institut Pasteur (Paris), France.
## SPDX-License-Identifier: GPL-3.0-or-later
##

# Format one estimate and its credible-interval bounds into `fmt`.
format_cri_num_ <- function(pe, low, up, fmt, as_integer) {
  if (is.na(pe)) return(NA_character_)
  if (isTRUE(pe == low && pe == up)) return(sprintf(fmt, pe, pe, pe))
  if (isTRUE(pe < low) || isTRUE(pe > up)) warn_h(
    "Found an inconsistent credible interval: {pe}, CrI = [{low}, {up}]."
  )

  # decimals for 2 significant digits of the scale
  # * interval covering zero -> scale is the width;
  # * else -> smaller of estimate and 10 times the width, NA bound skipped;
  # * zero scale, from an inconsistent interval -> integers
  scale <- up - low
  digits <- 0
  if (!isTRUE(low <= 0 && up >= 0))
    scale <- min(abs(pe), 10 * scale, na.rm = TRUE)
  if (isFALSE(as_integer) && scale > 0)
    digits <- max(0, 1 - floor(log10(scale)))

  pe_f <- formatC(pe, digits = digits, format = "f")
  low_f <- NA_character_
  up_f <- NA_character_
  if (!is.na(low)) low_f <- formatC(low, digits = digits, format = "f")
  if (!is.na(up)) up_f <- formatC(up, digits = digits, format = "f")

  sprintf(fmt, pe_f, low_f, up_f)
}

# Format estimates and bounds element-wise, after checking lengths and `fmt`.
format_cri_vector_ <- function(point_estimate, cri_lower, cri_upper,
                               fmt, as_integer) {
  if (length(point_estimate) != length(cri_lower) ||
        length(point_estimate) != length(cri_upper)) stop_h(
    "Point estimate and CrI bounds must be numeric of the same length."
  )

  if (length(gregexpr("%s", fmt, fixed = TRUE)[[1]]) != 3) stop_h(
    "Invalid choice for character {.arg str_format}.",
    "It must contain exactly three {.code %s} placeholders."
  )

  mapply(format_cri_num_, point_estimate, cri_lower, cri_upper,
         MoreArgs = list(fmt = fmt, as_integer = as_integer),
         USE.NAMES = FALSE)
}

# Format tabular estimates and bounds column-wise, after checking dimensions.
format_cri_matrix_ <- function(point_estimate, cri_lower, cri_upper,
                               fmt, as_integer) {
  if (!all(dim(point_estimate) == dim(cri_lower),
           dim(point_estimate) == dim(cri_upper)))
    stop_h("Tabular inputs must have the same dimensions.")

  res <- as.data.frame(Map(
    format_cri_vector_, as.data.frame(point_estimate),
    as.data.frame(cri_lower), as.data.frame(cri_upper),
    MoreArgs = list(fmt = fmt, as_integer = as_integer)
  ))
  names(res) <- names(point_estimate)
  res
}

#' Format credible interval into printable string
#'
#' @description
#' Format element-wise estimate and credible-interval into a string.
#' The function is vectorised over its arguments
#' and returns a character vector of the same length as the inputs.
#'
#' Estimate and bounds are rounded to the same number of decimals,
#' enough to keep at least two significant digits of the estimate
#' and to tell the bounds apart
#' \insertCite{gelman_bayesian_2026}{housepi}.
#' When the interval spans zero, the decimals keep two significant digits
#' of the interval width instead.
#' @param point_estimate `numeric()` Point estimate values.
#' @param cri_lower `numeric()` Lower bounds of the credible interval.
#' @param cri_upper `numeric()` Upper bounds of the credible interval.
#' @param fmt `character(1)`
#' Sprintf-style format with 3 `%s` placeholders,
#' filled with the point estimate, then the lower and upper bounds.
#' @param as_integer `logical(1)` Whether values are formatted as integers.
#' @return `character()`
#' @references
#' \insertAllCited{}
#' @examples
#' pe <- c(9.9, 7.69e-02, 1.13)
#' cri_low <- c(2.6, 6.77e-02, NA)
#' cri_up <- c(23.7, 0.087, 1.53)
#' format_credible_interval(pe, cri_low, cri_up)
#' # custom format
#' format_credible_interval(pe, cri_low, cri_up,
#'                          fmt = "%s (95%% CrI, %s-%s)")[2]
#' @export
format_credible_interval <- function(point_estimate, cri_lower, cri_upper,
                                     fmt = "%s [%s, %s]",
                                     as_integer = FALSE) {
  tabular <- is.data.frame(point_estimate) || is.matrix(point_estimate)
  format_cri <- if (tabular) format_cri_matrix_ else format_cri_vector_
  format_cri(point_estimate, cri_lower, cri_upper, fmt, as_integer)
}
