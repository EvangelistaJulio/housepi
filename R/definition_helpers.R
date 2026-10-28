## housepi - Epidemiological Analysis of Household Outbreaks.
## Copyright 2024-2026 Institut Pasteur (Paris), France.
## SPDX-License-Identifier: GPL-3.0-or-later
##

# Extract balanced [...] blocks at the top level.
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

# Evaluate blocks, recursing on nested ones.
helper_evaluate_defin_block <- function(block, type) {
  # remove square brackets
  block <- sub("^\\[", "", block)
  block <- sub("\\]$", "", block)

  # split the block into individual conditions
  conditions <- unlist(strsplit(block, ":"))
  output <- c()

  # evaluate each condition
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
#' @param type `character(1)` Type to extract.
#' @noRd
housepi_walk_definition <- function(input, type) {
  output <- c()
  if (length(input) && !is.na(input) && nzchar(input)) {
    blocks <- helper_extract_defin_blocks(input)
    if (length(blocks)) {
      for (b in blocks) {
        output <- c(output, helper_evaluate_defin_block(b, type))
      }
      if (length(output)) output <- unique(output)
    } else {
      # terminal node
      parts <- unlist(strsplit(input, "\\|"))
      if (length(parts) == 2 && parts[2] == type) output <- parts[1]
    }
  }
  output
}

#' Resolve unsettled sources in encoded definition strings.
#' * `name|?` terminal -> the letter the dictionary holds for that name;
#' * name the dictionary does not settle -> keeps the sentinel, reported later.
#' @param dict `character()` Source letter, named by declared name.
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

# Names read by the definitions from source `src`, labelled with it.
find_required_definitions <- function(defs, src, col_with_sources, family) {
  output <- c()
  if (length(defs)) {
    for (col in col_with_sources) {
      for (x in defs[[col]]) {
        output <- c(output, housepi_walk_definition(x, src))
      }
    }
  }
  # remove duplicates and add source name
  set_source_name_(unique(output), family)
}

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

# Encode quosure `func`, a formula's right side included:
# * bare name -> alias;
# * string -> function name;
# * otherwise -> number or call.
helper_load_func_ <- function(func) {
  expr <- if (rlang::is_quosure(func)) rlang::quo_get_expr(func)
  if (rlang::is_formula(expr)) expr <- rlang::f_rhs(expr)
  is_alias <- rlang::is_symbol(strip_useless_parens(expr))
  func <- housepi_extractor(!!func)
  if (!is.null(func)) housepi_combine_to_call(func, character(), !is_alias)
}

# Assert `const` is a quosure or value other than a logical.
assert_const_value_ <- function(const, family, name) {
  expr <- if (rlang::is_quosure(const)) rlang::quo_get_expr(const) else const
  if (is.logical(expr)) stop_h(
    "{.arg const} of {family} {.val {name}} must be the value it is pinned to.",
    format_type_error(expr)
  )
  invisible(const)
}

# Assert a definition holds one of a constant, a law or a function:
# * none, unless a bounded latent -> error;
# * function beside a law, a constant or a proposal -> error;
# * constant beside a prior or a proposal -> error;
# * bounded constant with no law -> error.
# @param law `character(1)` Prior or model, `NA` if none.
# @param bounded `logical(1)` Whether `min` or `max` is given.
assert_definition_core_ <- function(family, name, const, func, law, proposal,
                                    bounded) {
  dcol <- family2model_col(family)
  # nolint next: object_usage_linter. read by the cli glue below
  keys <- switch(family, observation = c(dcol, "func"),
                 `household observation` = "func", c("const", dcol, "func"))
  # a bounded latent with no model is bare
  is_bare <- bounded && family == "latent"
  if (all(is.null(c(const, func)), is.na(law), !is_bare)) stop_h(
    "Incomplete definition for {family} {.val {name}}.",
    "Provide {cli::qty(keys)}{?/at least one of} {.arg {keys}}."
  )
  wrong <- "Wrong definition for {family} {.val {name}}."
  if (!is.null(func) && !is.na(law)) stop_h(
    wrong, "Provide only one of {.arg {c('func', dcol)}}."
  )
  if (!is.null(const) && !is.null(func)) stop_h(
    wrong, "Provide only one of {.arg {c('const', 'func')}}."
  )
  if (!is.null(func) && !is.na(proposal)) stop_h(
    "Deterministic {family} {.val {name}} cannot be sampled.",
    "Drop {.arg proposal}, or drop {.arg func}."
  )
  if (is.null(const)) return(invisible())
  if (dcol == "prior" && !is.na(law)) stop_h(
    wrong, "Provide only one of {.arg {c('const', 'prior')}}."
  )
  if (!is.na(proposal)) stop_h(
    "Constant {family} {.val {name}} cannot be sampled.",
    "Drop {.arg proposal} to keep it constant."
  )
  if (bounded && is.na(law)) stop_h(
    "Constant {family} {.val {name}} has no distribution to bound.",
    "Drop {.arg min} and {.arg max} to keep it constant."
  )
  invisible()
}

# Attach the spare arguments to the object that takes them: the distribution
# sampling the variable, or the function computing it.
helper_variable_args <- function(opts, family, args) {
  if (!length(args)) return(opts)
  dcol <- family2model_col(family)
  col <- if (!is.na(opts[[dcol]])) dcol else "func"
  headline <- paste("{cli::qty(args)}There {?was/were} {length(args)} spare",
                    "argument{?s} for {family} {.val {opts$name}}:",
                    "{.arg {names(args)}}.")
  if (is.na(opts[[col]])) stop_h(
    headline, "Provide the {.arg {dcol}} or the {.arg func} they belong to."
  )
  # an alias or a call -> already complete
  if (housepi_bar_count(opts[[col]]) > 0) stop_h(
    headline, "The {.arg {col}} expression already carries its arguments."
  )
  opts[[col]] <- housepi_combine_to_call(opts[[col]], args)
  opts
}

# Random-walk proposal kernels. Any other proposal is an independence sampler:
# * "forward" -> draws from the variable's own defining distribution;
# * distribution name -> draws from that distribution.
.random_walk_proposals <- c("rw", "rw_log", "rw_unif")

# Flatten dots holding a single quosures list (from `enquos()`).
flatten_definition_dots <- function(dots) {
  if (length(dots) == 1 && rlang::is_quosures(dots[[1]])) dots[[1]]
  else dots
}

# Extract the arguments in the dots of definition `name` of `family`.
# Advanced `proposal_*` options fold into the proposal call:
# * `proposal_scale` -> random walks;
# * distribution arguments -> independence samplers.
helper_definition_args <- function(proposal, dots, family, name = NULL) {
  args <- housepi_extractor(!!!flatten_definition_dots(dots))
  nms <- names(args)
  assert_unique_arg_names(nms, family, name)
  # skip if empty
  if (!length(nms)) return(list(args = args, proposal = proposal))

  # find proposal advanced definitions
  mask_prop <- startsWith(nms, "proposal_")
  prop_args <- args[mask_prop]
  names(prop_args) <- sub("^proposal_", "", names(prop_args))

  # validate proposal arguments against the chosen family
  if (length(prop_args)) {
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

# Interpret user choice for Monte Carlo move.
# * constant or function -> never sampled;
# * latent with neither model nor function -> bare, sampled within its bounds;
# * otherwise -> sampled from its prior or model.
helper_variable_move <- function(opts, family) {
  if (!is.na(opts$const) || !is.na(opts$func)) return(opts)
  # prior for parameters/rates, model for latent/observations
  dist_word <- family2model_col(family)
  if (identical(opts$kind, "latent") && is.na(opts[[dist_word]])) {
    headline <- paste("No {.arg {dist_word}} to sample from",
                      "for {family} {.val {opts$name}}.")
    if (identical(opts$proposal, "forward")) stop_h(
      headline, "Proposal {.val forward} draws from the {.arg {dist_word}}."
    )
    if (anyNA(c(opts$min, opts$max))) stop_h(
      headline,
      "Provide a {.arg {dist_word}}, or both {.arg min} and {.arg max}."
    )
  }
  opts
}

# Reduction functions: folding a per-member column yields one household value.
.reduction_functions <- c("sum", "prod", "min", "max")

# Settle the level a deterministic observation is read at.
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

  # fold over a single per-member operand -> household reduction of its column
  operands <- showcase_extract_args(opts$func)[[1]]
  if (length(operands) != 1) return(opts)
  operand <- strsplit(as.character(operands[[1]]), "\\|")[[1]]
  if (length(operand) != 2 || !operand[2] %in% c("D", "?")) return(opts)

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

# Cast to the household level a deterministic member-level observation reading
# household values only, parameters and household observations.
# * constant -> stays put.
# @param defs `data.frame()` Definitions with sources resolved.
helper_cast_house_level <- function(defs) {
  if (is.null(defs)) return(defs)
  reads <- function(r, src) {
    unlist(lapply(c("const", "func", "min", "max"),
                  \(col) housepi_walk_definition(defs[[col]][[r]], src)))
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
    if (!length(cast)) return(defs)
    defs[["kind"]][cast] <- "house"
  }
}
