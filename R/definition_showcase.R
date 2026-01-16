## housepi - Epidemiological Analysis of Household Outbreaks.
##
## Copyright © 2026 Institut Pasteur (Paris), CNRS and INSERM.
## Distributed under the terms of the GNU General Public License (GPLv3).
## See the LICENSE file for details.
##

# Helper to extract sourced value from a definition.
showcaser_extract <- function(x, src) {
  sapply(x, function(element) {
    if (is.na(element)) return(NA_character_)
    if (!is.character(element)) return(element)
    z <- strsplit(element, "\\|")[[1]]
    if (length(z) != 2 || z[2] != src) return(element)
    as.character(z[1])
  }, USE.NAMES = FALSE)
}

# helper to extract numeric values from a definition
showcaser_extract_numeric <- function(x) {
  sapply(x, function(element) {
    if (is.na(element)) return(NA_real_)
    if (!is.character(element)) return(element)
    z <- strsplit(element, "\\|")[[1]]
    if (length(z) != 2 || z[2] != "N") return(element)
    as.numeric(z[1])
  }, USE.NAMES = FALSE)
}

# helper to extract function name from a definition
showcaser_extract_func <- function(x) {
  sapply(x, function(element) {
    if (is.na(element)) return(NA_real_)
    if (!is.character(element)) return(element)
    z <- strsplit(element, "\\|")[[1]]
    if (length(z) <= 2) return(element)
    z[1]
  }, USE.NAMES = FALSE)
}

# helper to split most external semicolon separator
helper_split_external_semicol <- function(x) {
  chars <- strsplit(x, "")[[1]]
  depth <- 0
  split_pos <- c(0)  # start of string
  for (i in seq_along(chars)) {
    if (chars[i] == "{") depth <- depth + 1
    else if (chars[i] == "}") depth <- depth - 1
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

  # remove outer braces at the edges
  gsub("^\\{(.*)\\}$", "\\1", result)
}

# helper to extract arguments from a definition
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

#' Convert encoded definition strings into a more human-readable format.
#' @param defs `data.frame()`. Encoded definition strings.
#' @return `data.frame()`. Readable definitions.
#' @keywords internal
#' @noRd
housepi_defin_showcaser <- function(defs) {
  defs <- as.data.frame(defs)
  output <- as.data.frame(lapply(defs, showcaser_extract_numeric))
  for (col in colnames(defs)) {
    # find function definition
    if (isTRUE(any(stringr::str_count(output[[col]], "\\|") > 1))) {
      pos <- match(col, names(output))
      # format column with arguments
      # UPDATE: this is not clean at all.
      col_args <- data.frame(
        args = base::I(showcase_extract_args(output[[col]]))
      )
      class(col_args$args) <- "list"
      if (!col %in% c("model", "prior")) {
        col_args <- stats::setNames(col_args, paste0(col, "_args"))
      }

      right_side_table <- cbind(
        stats::setNames(list(showcaser_extract_func(output[[col]])), col),
        col_args,
        output[(pos + 1):ncol(output)]
      )
      if (pos > 1) output <- cbind(output[1:(pos - 1)], right_side_table)
      else output <- right_side_table
    }
  }
  output
}
