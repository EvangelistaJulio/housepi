## housepi - Epidemiological Analysis of Household Outbreaks.
##
## Copyright © 2026 Institut Pasteur (Paris), CNRS and INSERM.
## Distributed under the terms of the GNU General Public License (GPLv3).
## See the LICENSE file for details.
##

# Helper to extract balanced {something} blocks at the top level.
helper_extract_defin_blocks <- function(s) {
  blocks <- list()
  depth <- 0
  start <- NULL
  for (i in seq_along(strsplit(s, "")[[1]])) {
    ch <- substr(s, i, i)
    if (ch == "{") {
      if (depth == 0) start <- i
      depth <- depth + 1
    } else if (ch == "}") {
      depth <- depth - 1
      if (depth == 0 && !is.null(start)) {
        blocks <- c(blocks, substr(s, start, i))
        start <- NULL
      }
    }
  }
  blocks
}

# Recursive helper function to evaluate blocks.
helper_evaluate_defin_block <- function(block, type) {
  # Remove curly braces
  block <- sub("^\\{", "", block)
  block <- sub("\\}$", "", block)

  # split the block into individual conditions
  conditions <- unlist(strsplit(block, ":"))
  output <- c()

  # Evaluate each condition
  for (cond in conditions) {
    if (grepl("\\{", cond)) {
      # recurse on nested blocks
      output <- c(output, housepi_walk_definition(cond, type))
    } else {
      pts <- unlist(strsplit(cond, "\\|"))
      if (length(pts) == 2 && pts[2] == type) output <- c(output, pts[1])
    }
  }
  output
}

#' Helper with type extraction from definitions.
#' @param input `character(1)`. Encoded function string.
#' @param type `character(1)`. Type to extract.
#' @return `character()`
#' @keywords internal
#' @noRd
housepi_walk_definition <- function(input, type) {
  output <- c()
  if (!is.na(input) && !is.null(input) && input != "") {
    blocks <- helper_extract_defin_blocks(input)
    if (length(blocks) > 0) {
      for (b in blocks) {
        output <- c(output, helper_evaluate_defin_block(b, type))
      }
      if (length(output) > 0) output <- unique(output)
    } else {
      # terminal node
      parts <- unlist(strsplit(input, "\\|"))
      if (length(parts) == 2 && parts[2] == type) output <- parts[1]
    }
  }
  output
}

# Definitions required by definitions according to the "source"
find_required_definitions <- function(defs, src, col_with_sources, family) {
  output <- c()
  if (length(defs) > 0) {
    for (col in col_with_sources) {
      for (x in defs[[col]]) {
        output <- c(output, housepi_walk_definition(x, src))
      }
    }
  }
  # remove duplicates and add source name
  set_source_name_(unique(output), family)
}
