## housepi - Epidemiological Analysis of Household Outbreaks.
##
## Copyright © 2024-2026 Institut Pasteur (Paris), CNRS and INSERM.
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
  if (is.null(value) && is.null(func)) return(opts)
  if (is.numeric(showcaser_extract_numeric(value))) {
    # implicit numeric definition
    opts$value <- value
    if (is.na(opts$const)) {
      opts$const <- isTRUE(all(is.na(c(opts[["prior"]], opts$proposal))))
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
  if (is.null(args)) return(opts)
  if (is.na(opts[["prior"]])) {
    sms <- "Spare arguments for {family} {.val {opts$name}}."
    if (!is.na(opts$func)) sms <- c(sms, "Did you mean {.arg func_args}?")
    stop_h(sms)
  }
  opts[["prior"]] <- housepi_combine_to_call(opts[["prior"]], args)
  opts
}

# Extract arguments carried by the dots.
# Fold any advanced `proposal_*` option into the proposal call.
helper_definition_args <- function(proposal, dots) {
  args <- housepi_extractor(!!!dots)
  nms <- names(args)
  # skip if empty
  if (length(nms) == 0) return(list(args = args, proposal = proposal))

  # find proposal advanced definitions
  mask_prop <- startsWith(nms, "proposal_")

  #>>>>
  # UPDATE: temporary stop advanced proposal setup
  unknown <- setdiff(nms[mask_prop], "proposal_scale")
  if (length(unknown) > 0) stop_h(
    "Unsupported advanced option{?s} {.arg {unknown}}.",
    "Only {.arg proposal_scale} is currently available."
  )
  #<<<<

  prop_args <- args[mask_prop]
  names(prop_args) <- stringr::str_remove(names(prop_args), "^proposal_")

  # fold user-defined proposal arguments into the proposal call
  if (!is.na(proposal) && proposal != "sample") {
    proposal <- housepi_combine_to_call(proposal, prop_args)
  }

  list(args = args[!mask_prop], proposal = proposal)
}

# Interpret user choice for Monte Carlo move.
helper_variable_move <- function(opts, family) {
  # infer if sampling proposal but no constant parameter provided
  if (is.na(opts$const)
      && !(is.na(opts$proposal) && is.na(opts[["prior"]]))) {
    opts$const <- FALSE
  }

  if (isTRUE(opts$const)) {
    if (is.na(opts$value)) stop_h(
      "Missing value for constant {family} {.val {opts$name}}."
    )
    opts$proposal <- NA_character_
    # drop prior if constant
    if (family == "parameter") opts[["prior"]] <- NA_character_
    # skip proposal if constant
    return(opts)
  } else {
    # check for prior if not constant (skip for func-defined parameters)
    if (family == "parameter" && is.na(opts[["prior"]])
        && is.na(opts$func)) stop_h(
      "No {.arg prior} provided for {family} {.val {opts$name}}."
    )
  }

  if (is.na(opts$proposal)) {
    if (is.na(opts[["prior"]])) {
      if (is.na(opts$func) &&
            is.numeric(showcaser_extract_numeric(opts$value))) {
        # raise error if const is FALSE (not if NA or already TRUE)
        if (isFALSE(opts$const)) stop_h(
          "No {.arg prior} to sample from for {family} {.val {opts$name}}."
        )
        opts$const <- TRUE
      }
    } else {
      # otherwise default proposal
      opts$proposal <- "sample"
    }
  } else {
    if (is.na(opts[["prior"]])) stop_h(
      "No {.arg prior} to sample from for {family} {.val {opts$name}}."
    )
  }

  opts
}
