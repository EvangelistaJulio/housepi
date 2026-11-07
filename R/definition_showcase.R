## housepi - Epidemiological Analysis of Household Outbreaks.
## Copyright 2024-2026 Institut Pasteur (Paris), France.
## SPDX-License-Identifier: GPL-3.0-or-later
##

# Map `fn(fields, element)` over encoded definitions,
# `fields` the '|'-separated parts of a character element.
# * NA element -> `na`, a value or a function of the element;
# * other non-character element -> passed through.
showcase_map_ <- function(x, fn, na = NA_character_) {
  sapply(x, function(element) {
    if (is.na(element)) return(if (is.function(na)) na(element) else na)
    if (!is.character(element)) return(element)
    fn(strsplit(element, "\\|")[[1]], element)
  }, USE.NAMES = FALSE)
}

# Helper to extract sourced value from a definition.
showcaser_extract <- function(x, src) {
  showcase_map_(x, function(z, element) {
    if (length(z) != 2 || z[2] != src) element
    else as.character(z[1])
  })
}

# Extract the numeric value of a definition.
showcaser_extract_numeric <- function(x) {
  showcase_map_(x, function(z, element) {
    if (length(z) != 2 || z[2] != "N") return(element)
    if (z[1] %in% .reserved_math_consts) return(get(z[1]))
    as.numeric(z[1])
  }, na = \(e) if (is.logical(e)) e else NA_real_)
}

# Extract the function name of a definition.
showcaser_extract_func <- function(x) {
  showcase_map_(x, \(z, element) if (length(z) <= 2) element else z[1],
                na = NA_real_)
}

# Split on the outermost semicolon separators.
helper_split_external_semicol <- function(x) {
  chars <- strsplit(x, "")[[1]]
  depth <- 0
  # start of string
  split_pos <- c(0)
  for (i in seq_along(chars)) {
    if (chars[i] == "[") depth <- depth + 1
    else if (chars[i] == "]") depth <- depth - 1
    else if (chars[i] == ":" && depth == 0) split_pos <- c(split_pos, i)
  }
  split_pos <- c(split_pos, nchar(x) + 1)

  # extract sub-strings between split positions
  result <- character(length(split_pos) - 1)
  for (i in seq_along(result)) {
    start <- split_pos[i] + 1
    end <- split_pos[i + 1] - 1
    result[i] <- substr(x, start, end)
  }

  # remove outer brackets at the edges
  gsub("^\\[(.*)\\]$", "\\1", result)
}

# Extract the arguments of a definition.
showcase_extract_args <- function(x) {
  lapply(x, function(element) {
    if (is.na(element)) return(NA_real_)
    if (!is.character(element)) return(c(element, NA_real_))
    z <- strsplit(element, "\\|")[[1]]
    if (length(z) <= 2) return(NA_real_)
    unlist(stats::setNames(
      showcaser_extract_numeric(
        helper_split_external_semicol(paste(z[3:length(z)], collapse = "|"))
      ),
      helper_split_external_semicol(z[2])
    ))
  })
}

#' Convert encoded definitions into a readable format.
#' @return `data.frame()` Readable definitions.
#' @noRd
housepi_defin_showcaser <- function(defs) {
  defs <- as.data.frame(defs)
  output <- as.data.frame(lapply(defs, showcaser_extract_numeric))
  for (col in colnames(defs)) {
    # find function definition
    if (isTRUE(any(housepi_bar_count(output[[col]]) > 1))) {
      pos <- match(col, names(output))
      # format column with arguments
      # @todo build the argument columns more cleanly
      col_args <- data.frame(
        args = base::I(showcase_extract_args(output[[col]]))
      )
      class(col_args$args) <- "list"
      if (!col %in% c("prior", "model"))
        col_args <- stats::setNames(col_args, paste0(col, "_args"))

      right_side <- cbind(
        stats::setNames(list(showcaser_extract_func(output[[col]])), col),
        col_args
      )
      if (pos < ncol(output))
        right_side <- cbind(right_side, output[(pos + 1):ncol(output)])
      if (pos > 1) output <- cbind(output[1:(pos - 1)], right_side)
      else output <- right_side
    }
  }
  output
}

#' Order readable definitions for display:
#' * generic definitions, without `status`, first;
#' * then by status, numeric when every status is a number, else alphabetical;
#' * then by name.
#' @param defs `data.frame()` | `NULL` Readable definitions.
#' @param name `character(1)` Column holding the definition name.
#' @param by_status `logical(1)` Whether to order by `status` first.
#' `FALSE` for tables where `status` is a stratum, not a definition key.
#' @return `data.frame()` `defs`, rows reordered.
#' @noRd
sort_definitions_ <- function(defs, name = "name", by_status = TRUE) {
  if (NROW(defs) < 2L) return(defs)
  sort_key <- function(x) {
    num <- suppressWarnings(as.numeric(x))
    if (anyNA(num[!is.na(x)])) tolower(as.character(x)) else num
  }
  keys <- lapply(intersect(c(if (by_status) "status", name), names(defs)),
                 \(col) sort_key(defs[[col]]))
  rows <- do.call(order, c(keys, list(na.last = FALSE, method = "radix")))
  out <- defs[rows, , drop = FALSE]
  rownames(out) <- NULL
  out
}
