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

#' Extract a typed argument from an encoded definition string.
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

# Parse implicit function/value.
helper_variable_core <- function(opts, family, value, func, func_args) {
  col_model <- family2model_col(family)
  if (is.null(value) && is.null(func)) return(opts)
  if (is.numeric(showcaser_extract_numeric(value))) {
    # implicit numeric definition
    opts$value <- value
    if (is.na(opts$const)) {
      opts$const <- isTRUE(all(is.na(c(opts[[col_model]], opts$proposal))))
    }
  } else {
    # NOTE: ensured just one between the 2 is not NULL
    x <- if (is.null(value)) func else value
    # distinguish purely implicit function to plain name
    y <- housepi_combine_to_call(func_name = x,
                                 func_args = func_args,
                                 is_func = is.null(value))
    if (housepi_is_sourced_value(y)) opts$value <- y
    else opts$func <- y

    # note: a priori a function is not constant
    if (is.na(opts$const)) opts$const <- FALSE
  }
  opts
}

# process model arguments
helper_variable_model <- function(opts, family, args) {
  col_model <- family2model_col(family)

  if (is.null(args)) return(opts)
  if (is.na(opts[[col_model]])) {
    sms <- sprintf("Spare arguments for %s '%s'.", family, opts$name)
    if (!is.na(opts$func)) sms <- c(sms, "Did you mean 'func_args'?")
    stop_h(sms)
  }
  opts[[col_model]] <- housepi_combine_to_call(opts[[col_model]], args)
  opts
}

# Interpret user choice for Monte Carlo move.
helper_variable_move <- function(opts, family, proposal_scale) {
  col_model <- family2model_col(family)

  # infer if sampling proposal but no constant parameter provided
  if (is.na(opts$const)
      && !(is.na(opts$proposal) && is.na(opts[[col_model]]))) {
    opts$const <- FALSE
  }

  if (isTRUE(opts$const)) {
    if (is.na(opts$value)) stop_h(
      sprintf("Missing value for constant %s '%s'.", family, opts$name)
    )
    opts$proposal <- NA_character_
    # drop prior if constant
    if (family == "parameter") opts[[col_model]] <- NA_character_
    # skip proposal if constant
    return(opts)
  } else {
    # check for prior if not constant (skip for func-defined parameters)
    if (family == "parameter" && is.na(opts[[col_model]])
        && is.na(opts$func)) stop_h(
      sprintf("No %s provided for parameter '%s'.", col_model, opts$name)
    )
  }

  if (is.na(opts$proposal)) {
    if (is.na(opts[[col_model]])) {
      if (is.na(opts$func) &&
            is.numeric(showcaser_extract_numeric(opts$value))) {
        # raise error if const is FALSE (not if NA or already TRUE)
        if (isFALSE(opts$const)) stop_h(
          sprintf("No %s to sample from for %s '%s'.",
                  col_model, family, opts$name)
        )
        opts$const <- TRUE
      }
    } else {
      # otherwise default proposal
      opts$proposal <- "sample"
    }
  } else {
    if (is.na(opts[[col_model]])) stop_h(
      sprintf("No %s to sample from for %s '%s'.", col_model, family, opts$name)
    )
    if (opts$proposal != "sample")  {
      if (is.na(proposal_scale)) opts$proposal_scale <- 1
      else if (is.numeric(proposal_scale) &&
                 proposal_scale > 0) opts$proposal_scale <- proposal_scale
      else stop_h("Proposal scale must be a numeric greater than 0.")
    }
  }

  opts
}
