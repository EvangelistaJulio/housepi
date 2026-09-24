## housepi - Epidemiological Analysis of Household Outbreaks.
##
## Copyright © 2024-2026 Institut Pasteur (Paris), CNRS and INSERM.
## Distributed under the terms of the GNU General Public License (GPLv3).
## See the LICENSE file for details.
##

# @description
# Helper to extract balanced [something] blocks at the top level.
helper_extract_defin_blocks <- function(s) {
  blocks <- list()
  depth <- 0
  start <- NULL
  for (i in seq_along(strsplit(s, "")[[1]])) {
    ch <- substr(s, i, i)
    if (ch == "[") {
      if (depth == 0) start <- i
      depth <- depth + 1
    } else if (ch == "]") {
      depth <- depth - 1
      if (depth == 0 && !is.null(start)) {
        blocks <- c(blocks, substr(s, start, i))
        start <- NULL
      }
    }
  }
  blocks
}

# @description
# Recursive helper function to evaluate blocks.
helper_evaluate_defin_block <- function(block, type) {
  # Remove square brackets
  block <- sub("^\\[", "", block)
  block <- sub("\\]$", "", block)

  # split the block into individual conditions
  conditions <- unlist(strsplit(block, ":"))
  output <- c()

  # Evaluate each condition
  for (cond in conditions) {
    if (grepl("\\[", cond)) {
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

#' Resolve unsettled source sigils in encoded definition strings.
#' A `name|?` terminal takes the letter the dictionary holds for that name;
#' a name the dictionary does not settle keeps the sentinel, for the assertion
#' to report.
#' @param input `character()`. Encoded definition strings.
#' @param dict named `character()`. Declared name to source letter.
#' @return `character()`
#' @keywords internal
#' @noRd
housepi_resolve_definition <- function(input, dict) {
  hit <- !is.na(input) & grepl("|?", input, fixed = TRUE)
  if (!any(hit)) return(input)
  terms <- gregexpr("[A-Za-z0-9_.]+\\|\\?", input[hit])
  regmatches(input[hit], terms) <- lapply(
    regmatches(input[hit], terms),
    function(terminal) {
      nm <- sub("\\|\\?$", "", terminal)
      src <- unname(dict[nm])
      ifelse(is.na(src) | !nzchar(src), terminal, paste0(nm, "|", src))
    }
  )
  input
}

# @description
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

# @description
# Catch the sibling defining-distribution keyword landing in the dots.
assert_no_sibling_dist_arg <- function(args, family, name = NULL) {
  dcol <- family2model_col(family)
  spare <- setdiff(c("prior", "model"), dcol)
  if (!spare %in% names(args)) return(invisible())
  holder <- if (is.null(name)) "{family}" else "{family} {.val {name}}"
  stop_h(
    paste0("Spare argument {.arg {spare}} for ", holder, "."),
    "Did you mean {.arg {dcol}} instead of {.arg {spare}}?"
  )
}

# @description
# Parse implicit function/value.
helper_variable_core <- function(opts, family, value, func) {
  if (is.null(value) && is.null(func)) return(opts)

  # NOTE: ensured just one between the 2 is not NULL
  x <- if (is.null(value)) func else value
  if (is.numeric(showcaser_extract_numeric(x))) {
    # implicit numeric definition, whichever argument carried it
    opts$value <- x
    if (is.na(opts$const)) {
      dcol <- family2model_col(family)
      opts$const <- isTRUE(all(is.na(c(opts[[dcol]], opts$proposal))))
    }
    return(opts)
  }

  # distinguish purely implicit function from plain name
  y <- housepi_combine_to_call(func_name = x,
                               func_args = character(),
                               is_func = is.null(value))
  # a value= given alongside a model/prior is a simulation-only override and
  # always lands in opts$value, whatever its shape; otherwise (including the
  # bulk `name = <expr>` shorthand, which is also passed in as `value`) the
  # shape alone decides, as a bare name/leaf and a compound call are handled
  # the same regardless of which argument carried it
  has_law <- !is.na(opts[[family2model_col(family)]])
  if ((!is.null(value) && has_law) || housepi_is_sourced_value(y)) {
    opts$value <- y
  } else {
    opts$func <- y
  }

  # note: a priori a function is not constant
  if (is.na(opts$const)) opts$const <- FALSE
  opts
}

# @description
# Attach the spare arguments to the object that takes them: the distribution
# sampling the variable, or the function computing it.
helper_variable_args <- function(opts, family, args, func) {
  if (length(args) == 0) return(opts)
  dcol <- family2model_col(family)
  col <- if (!is.na(opts[[dcol]])) dcol else if (!is.na(opts$func)) "func"
  headline <- paste("{cli::qty(args)}There {?was/were} {length(args)} spare",
                    "argument{?s} for {family} {.val {opts$name}}:",
                    "{.arg {names(args)}}.")
  if (is.null(col) && is.null(func)) stop_h(
    headline, "Provide the {.arg {dcol}} or the {.arg func} they belong to."
  )
  # a func folded into a constant carries no call to attach them to
  if (is.null(col)) col <- "func"
  if (is.na(opts[[col]]) || housepi_is_call(opts[[col]])) stop_h(
    headline, "The {.arg {col}} expression already carries its arguments."
  )
  opts[[col]] <- housepi_combine_to_call(opts[[col]], args)
  opts
}

# Random-walk proposal kernels (the marked family). All other proposal names
# are independence samplers: `"forward"` draws from the variable's own defining
# distribution, a distribution name draws from that distribution.
.random_walk_proposals <- c("rw", "rw_log", "rw_unif")

# Flatten dots holding a single quosures list (from `enquos()`).
flatten_definition_dots <- function(dots) {
  if (length(dots) == 1 && rlang::is_quosures(dots[[1]])) dots[[1]] else dots
}

# @description
# Extract arguments carried by the dots.
# Fold any advanced `proposal_*` option into the proposal call, routing
# `proposal_scale` to random walks and distribution arguments to independence
# samplers.
helper_definition_args <- function(proposal, dots) {
  args <- housepi_extractor(!!!flatten_definition_dots(dots))
  nms <- names(args)
  # skip if empty
  if (length(nms) == 0) return(list(args = args, proposal = proposal))

  # find proposal advanced definitions
  mask_prop <- startsWith(nms, "proposal_")
  prop_args <- args[mask_prop]
  names(prop_args) <- sub("^proposal_", "", names(prop_args))

  # validate proposal arguments against the chosen family
  if (length(prop_args) > 0) {
    arg_nms <- names(prop_args)
    is_rw <- !is.na(proposal) && proposal %in% .random_walk_proposals
    if (is.na(proposal)) stop_h(
      "Advanced {.arg proposal_*} options require an explicit {.arg proposal}."
    )
    if (identical(proposal, "forward")) stop_h(
      paste("Independence sampling from {.val forward}",
            "takes no {.arg proposal_*} arguments.")
    )
    if (is_rw && any(arg_nms != "scale")) stop_h(
      paste("Random-walk proposal {.val {proposal}}",
            "accepts only {.arg proposal_scale}."),
      x = "Unexpected {.arg {paste0('proposal_', setdiff(arg_nms, 'scale'))}}."
    )
    if (!is_rw && "scale" %in% arg_nms) stop_h(
      "{.arg proposal_scale} tunes a random walk.",
      paste("Use {.val rw}/{.val rw_log}/{.val rw_unif}, or pass",
            "distribution arguments (e.g. {.arg proposal_mean}) for",
            "independence sampling from {.val {proposal}}.")
    )
  }

  # fold user-defined proposal arguments into the proposal call
  if (!is.na(proposal) && proposal != "forward") {
    proposal <- housepi_combine_to_call(proposal, prop_args)
  }

  list(args = args[!mask_prop], proposal = proposal)
}

# @description
# Interpret user choice for Monte Carlo move.
helper_variable_move <- function(opts, family) {
  # prior for parameters/risks, model for latent/observations
  dist_word <- family2model_col(family)
  is_param <- family == "parameter"

  # infer if sampling proposal but no constant provided
  if (is.na(opts$const)
      && !(is.na(opts$proposal) && is.na(opts[[dist_word]]))) {
    opts$const <- FALSE
  }

  if (isTRUE(opts$const)) {
    if (is.na(opts$value)) stop_h(
      "Missing value for constant {family} {.val {opts$name}}."
    )
    opts$proposal <- NA_character_
    # drop prior if constant
    if (is_param) opts[[dist_word]] <- NA_character_
    # skip proposal if constant
    return(opts)
  } else if (is_param && is.na(opts[[dist_word]]) && is.na(opts$func)) {
    # missing prior/func for not constant parameter
    stop_h("No {.arg {dist_word}} provided for {family} {.val {opts$name}}.")
  }

  if (is.na(opts$proposal)) {
    if (is.na(opts[[dist_word]])) {
      if (is.na(opts$func) &&
            is.numeric(showcaser_extract_numeric(opts$value))) {
        # raise error if const is FALSE (not if NA or already TRUE)
        if (isFALSE(opts$const)) stop_h(paste(
          "No {.arg {dist_word}} to sample from",
          "for {family} {.val {opts$name}}."
        ))
        opts$const <- TRUE
      }
    } else {
      # otherwise default proposal: adaptive random walk
      opts$proposal <- "rw"
    }
  } else {
    if (!is.na(opts$func)) stop_h(
      "Deterministic {family} {.val {opts$name}} cannot be sampled.",
      "Drop {.arg proposal} and {.arg {dist_word}}, or drop {.arg func}."
    )
    if (is.na(opts[[dist_word]])) stop_h(
      "No {.arg {dist_word}} to sample from for {family} {.val {opts$name}}."
    )
  }

  opts
}

# Reduction functions: folding a per-member column yields one household value.
.reduction_functions <- c("sum", "prod", "min", "max")

# @description
# Settle the level a deterministic observation is read at
helper_variable_level <- function(opts, defs) {
  if (!identical(opts$kind, "member") ||
        !is.na(opts[["model"]]) ||
        !all(is.na(opts$status)) ||
        is.na(opts$func)) {
    return(opts)
  }
  # a bare reference folds nothing at all
  if (!housepi_is_call(opts$func)) return(opts)
  folded <- showcaser_extract_func(opts$func)
  if (!folded %in% .reduction_functions) return(opts)

  # a fold over a single per-member operand returns that operand: what the
  # definition asks for is the household reduction over its column
  operands <- showcase_extract_args(opts$func)[[1]]
  if (length(operands) != 1) return(opts)
  operand <- strsplit(as.character(operands[[1]]), "\\|")[[1]]
  if (length(operand) != 2 || operand[2] != "D") return(opts)

  sms <- paste("Member-level observation {.val {opts$name}}",
               "is a no-effect reduction {.fn {folded}}.")
  # a sibling household observation is already one value
  house <- character()
  if (!is.null(defs)) house <- defs[["name"]][defs[["kind"]] == "house"]
  if (operand[1] %in% house) stop_h(
    sms, "Drop redundant observation {.val {opts$name}}."
  )

  warn_h(
    sms, `!` = "Interpreted as an household-level observation instead.",
    paste("Silence this warning passing {.code level = \"house\"}",
          "or dropping {.fn {folded}}.")
  )
  opts$kind <- "house"
  opts
}

# @description
# Cast to the household level a deterministic member-level observation that
# reads household values only: parameters and household observations give
# every member the same value, one per household. A constant stays put.
# @param defs `data.frame()`. Definitions with sources resolved.
# @return `data.frame()`
helper_cast_house_level <- function(defs) {
  if (is.null(defs) || !"kind" %in% names(defs)) return(defs)
  reads <- function(r, src) {
    unlist(lapply(c("value", "func", "min", "max"), function(col) {
      housepi_walk_definition(defs[[col]][[r]], src)
    }))
  }
  is_candidate <- function(r) {
    identical(defs[["kind"]][r], "member") && is.na(defs[["model"]][r]) &&
      all(is.na(unlist(defs[["status"]][r]))) &&
      sum(defs[["name"]] == defs[["name"]][r]) == 1 &&
      length(c(reads(r, "A"), reads(r, "G"), reads(r, "?"))) == 0 &&
      length(c(reads(r, "D"), reads(r, "P"))) > 0
  }
  candidates <- Filter(is_candidate, seq_len(nrow(defs)))

  # a cast row is a household value to the rows reading it
  repeat {
    house <- defs[["name"]][defs[["kind"]] == "house"]
    cast <- Filter(function(r) {
      defs[["kind"]][r] == "member" && all(reads(r, "D") %in% house)
    }, candidates)
    if (length(cast) == 0) return(defs)
    defs[["kind"]][cast] <- "house"
  }
}
