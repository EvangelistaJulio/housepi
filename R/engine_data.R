## housepi - Epidemiological Analysis of Household Outbreaks.
## Copyright 2024-2026 Institut Pasteur (Paris), France.
## SPDX-License-Identifier: GPL-3.0-or-later
##

#' Order the levels of a categorical column, missing values dropped.
#' * factor -> its levels, unused ones included;
#' * otherwise -> sorted values, C-locale order on every machine.
#' @noRd
level_order_ <- function(x) {
  if (is.factor(x)) levels(x)
  else sort(unique(as.character(x[!is.na(x)])), method = "radix")
}

#' Code a character or factor status column into engine status codes.
#' * escaped label -> `0`;
#' * every other label -> `1..K` in level_order_(),
#'   whatever the row order;
#' * escaped label unmatched -> every individual counted as infected.
#' @param escaped `character(1)` Label of status code 0.
#' @return `list` Entries `codes`, `integer()`, and `labels`, `character()`.
#' @noRd
status_to_codes_ <- function(x, escaped = "escaped") {
  lv <- level_order_(x)
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
  labels <- setdiff(lv, escaped)
  code <- stats::setNames(c(0L, seq_along(labels)), c(escaped, labels))
  list(codes = unname(code[as.character(x)]), labels = labels)
}

#' Factorize the dates and non-numeric columns of `df` for the engine.
#' Column names cross verbatim: they are the keys the engine registers,
#' and the definitions sourcing them are written against the user's spelling.
#' @param escaped `character(1)` Label of status code 0,
#' in a character or factor `status` column.
#' @return `list` Entry `data`, `list(ids = , cells = )`: identifiers as text,
#' every other column as doubles, beside the coding cpp_data_refactor() reads.
#' @noRd
cpp_data_factorize <- function(df, escaped = "escaped") {
  cols <- colnames(df)
  ids <- intersect(c("id", "hh_id"), cols)
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
    s <- df[["status"]]
    frac <- unique(s[!is.na(s) & s != round(s)])
    if (length(frac)) stop_h(
      "Status{?es} {.val {frac}} {?is/are} not {?a /}whole number{?s}.",
      "Give {.field status} as whole numbers, or as labels."
    )
    s <- suppressWarnings(as.integer(s))
    status <- sort(unique(s[!is.na(s) & s > 0]))
  }

  discr <- stats::setNames(logical(length(cols)), cols)

  # identifiers cross as text, their class restored on export
  for (col in ids) {
    x <- df[[col]]
    if (is.factor(x)) info[[col]] <- list(class = "factor", levels = levels(x))
    else if (is.numeric(x)) info[[col]] <- list(class = typeof(x))
  }

  for (col in setdiff(cols, c(ids, "status"))) {
    x <- df[[col]]
    # every column but a plain double reads whole values
    discr[[col]] <- !is.double(x) || inherits(x, "Date")
    if (inherits(x, "Date")) info[[col]] <- list(class = "Date")
    else if (is.integer(x)) info[[col]] <- list(class = "integer")
    # FALSE reads 0 and TRUE reads 1
    else if (is.logical(x)) info[[col]] <- list(class = "logical")
    else if (is.character(x) || is.factor(x)) {
      cls <- if (is.character(x)) "character" else "factor"
      # a missing cell is missing data, not a category: it carries no level
      # and codes to NA, which the engine reads as an empty cell
      y <- factor(x, levels = level_order_(x))
      info[[col]] <- list(class = cls, levels = levels(y))
      df[[col]] <- as.integer(y)
    }
  }

  # every column but the identifiers crosses as doubles
  as_frame <- \(x, f) data.frame(lapply(x, f), check.names = FALSE)
  list(data = list(ids = as_frame(df[ids], as.character),
                   cells = as_frame(df[setdiff(cols, ids)], as.double)),
       info = info, cols = cols, discr = discr,
       status = status, status_labels = status_labels)
}

#' Restore the column classes of data returned by the engine.
#' * identifier -> its input class, read from its text;
#' * integer column a definition writes -> engine type kept;
#' * categorical column -> labels from its codes.
#' @param info `list` of `list`s Class and levels of the data columns.
#' @param headers `character()` Header columns.
#' @param status_levels `character()` Status levels.
#' @param escaped `character(1)` Label of status code 0.
#' @param produced `character()` Columns a definition writes.
#' @noRd
cpp_data_refactor <- function(data, info, headers,
                              status_levels = character(),
                              escaped = "escaped",
                              produced = character()) {
  for (col in intersect(names(info), colnames(data))) {
    x <- data[[col]]
    lv <- info[[col]]$levels
    # identifiers cross as labels, categorical columns as codes
    if (!is.null(lv) && !col %in% c("id", "hh_id")) x <- lv[as.integer(x)]
    data[[col]] <- switch(
      info[[col]]$class,
      Date = as.Date(x),
      logical = as.logical(x),
      double = as.numeric(x),
      integer = if (col %in% produced) x else as.integer(x),
      character = x,
      factor = factor(x, levels = lv)
    )
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
