## housepi - Epidemiological Analysis of Household Outbreaks.
## Copyright 2024-2026 Institut Pasteur (Paris), France.
## SPDX-License-Identifier: GPL-3.0-or-later
##

#' Format a dataset or a long dataset for the import Rcpp method.
#' Column names cross verbatim: they are the keys the engine registers,
#' and the definitions sourcing them are written against the user's spelling.
#' @noRd
char_data_input_ <- function(df) {
  if (is.null(df)) NULL
  else data.frame(lapply(df, as.character), check.names = FALSE)
}

#' Code a character or factor status column into engine status codes.
#' * escaped label -> `0`;
#' * every other label -> `1..K` in level or appearance order,
#'   so a simulated factor round-trips to the same codes;
#' * escaped label unmatched -> every individual counted as infected.
#' @param escaped `character(1)` Label of status code 0.
#' @return `list` Entries `codes`, `integer()`, and `labels`, `character()`.
#' @noRd
status_to_codes_ <- function(x, escaped = "escaped") {
  lv <- if (is.factor(x)) levels(x) else unique(as.character(x))
  lv <- lv[!is.na(lv)]
  # an all-missing column carries no outcome at all, so it stays silent
  if (length(lv) && !escaped %in% lv) {
    guess <- lv[tolower(lv) %in% c("susceptible", "negative", "control",
                                   "no", "false", "uninfected", "healthy", "0")]
    warn_h(
      "No {.field status} value matches the escaped label {.val {escaped}}.",
      `!` = "Every individual counts as infected. Values found: {.val {lv}}.",
      if (length(guess) == 1) "Did you mean {.code escaped = \"{guess}\"}?"
      else paste("Name the value carried by individuals who were not infected",
                 "with {.code escaped =}.")
    )
  }
  labels <- as.character(setdiff(lv, escaped))
  code <- stats::setNames(c(0L, seq_along(labels)), c(escaped, labels))
  list(codes = unname(code[as.character(x)]), labels = labels)
}

#' Factorize the dates and non-numeric columns of `df` for the engine.
#' @param escaped `character(1)` Label of status code 0,
#' in a character or factor `status` column.
#' @noRd
cpp_data_factorize <- function(df, escaped = "escaped") {
  skip <- c("id", "hh_id", "status")
  cols <- colnames(df)
  info <- list()

  # character/factor status: factorize to integer codes (NA -> missing)
  status_labels <- character()
  if ("status" %in% cols && (is.character(df[["status"]]) ||
                               is.factor(df[["status"]]))) {
    coded <- status_to_codes_(df[["status"]], escaped)
    df[["status"]] <- coded$codes
    status_labels <- coded$labels
  }

  # extract all infectious status codes (now integer-coded)
  status <- integer()
  if ("status" %in% cols) {
    s <- suppressWarnings(as.integer(df[["status"]]))
    status <- sort(unique(s[!is.na(s) & s > 0]))
  }

  discr <- stats::setNames(logical(length(cols)), cols)

  for (col in setdiff(cols, skip)) {
    x <- df[[col]]
    discr[[col]] <- inherits(x, "Date") || is.integer(x) || is.logical(x)
    if (inherits(x, "Date")) {
      info[[col]] <- list(class = "Date")
      df[[col]] <- as.numeric(x)
    } else if (is.logical(x)) {
      # FALSE reads 0 and TRUE reads 1
      info[[col]] <- list(class = "logical")
      df[[col]] <- as.integer(x)
    } else if (is.character(x) || is.factor(x)) {
      cls <- if (is.character(x)) "character" else "factor"
      # a missing cell is missing data, not a category: it carries no level
      # and codes to NA, which the engine reads as an empty cell
      y <- factor(x, levels = unique(x))
      info[[col]] <- list(class = cls, levels = levels(y))
      df[[col]] <- as.integer(y)
    }
  }

  list(data = char_data_input_(df), info = info, cols = cols, discr = discr,
       status = status, status_labels = status_labels)
}

#' Restore the factor columns of data returned by the engine.
#' @param info `list` of `list`s Class and levels of the data columns.
#' @param headers `character()` Header columns.
#' @param status_levels `character()` Status levels.
#' @param escaped `character(1)` Label of status code 0.
#' @noRd
cpp_data_refactor <- function(data, info, headers,
                              status_levels = character(),
                              escaped = "escaped") {
  # re-factor
  for (col in names(info)) {
    cls <- info[[col]]$class
    if (identical(cls, "Date")) {
      data[[col]] <- as.Date(data[[col]])
    } else if (identical(cls, "logical")) {
      data[[col]] <- as.logical(data[[col]])
    } else {
      lv <- info[[col]]$levels
      data[[col]] <- factor(data[[col]],
                            levels = seq_along(lv),
                            labels = lv)
      if (identical(cls, "character")) {
        data[[col]] <- as.character(data[[col]])
      }
    }
  }

  if ("status" %in% colnames(data)) {
    s <- as.integer(data[["status"]])
    if (length(status_levels)) {
      # code 0 -> escaped label, code k -> k-th infected label
      lv <- c(escaped, status_levels)
      data[["status"]] <- factor(s + 1L, levels = seq_along(lv), labels = lv)
    } else {
      data[["status"]] <- s
    }
  }
  # re-order data and append any extra columns
  data[c(headers, setdiff(colnames(data), headers))]
}
