## housepi - Epidemiological Analysis of Household Outbreaks.
## Copyright 2024-2026 Institut Pasteur (Paris), France.
## SPDX-License-Identifier: GPL-3.0-or-later
##

# Strata tables shared by group rates, contact rates and the status process:
# * a rule keys on covariates;
# * an NA cell matches any value;
# * the most specific matching rule of each rule group applies.

# Columns of strata table `tbl` a rule keys on: every one but the `reserved`.
strata_cov_cols_ <- function(tbl, reserved) setdiff(colnames(tbl), reserved)

# Assert rule table `tbl` is a data.frame with rows and the `required` columns.
# @param label `character(1)` Opening of the message.
# @param hint `character(1)` Hint following a missing column.
assert_strata_shape_ <- function(tbl, required, label, hint) {
  if (!is.data.frame(tbl) || !nrow(tbl)) stop_h(
    "{label} must be a data.frame with at least one row.",
    format_type_error(tbl)
  )
  missing <- setdiff(required, colnames(tbl))
  if (length(missing)) stop_h(
    "{label} has no column{?s} {.field {missing}}.", hint
  )
}

#' Assert no covariate takes a reserved name.
#' @param reserved `character()` Names no covariate may take.
#' @noRd
assert_valid_covariates <- function(covariates, reserved) {
  forbidden <- reserved[reserved %in% covariates]
  if (length(forbidden)) stop_h(
    "Impossible to stratify on covariate{?s}: {.field {forbidden}}.",
    "These names are reserved.",
    paste("Change {cli::qty(forbidden)}covariate{?s} and/or",
          "corresponding data column{?s}.")
  )
}

#' Assert the covariates are data columns.
#' @param covariates `list` Covariate data column names.
#' @param headers `character()` Data column names.
#' @noRd
assert_data_covariates <- function(covariates, headers) {
  missing <- setdiff(covariates, headers)
  if (length(missing)) stop_h(
    "Unknown covariate{?s}: {.field {missing}}.",
    find_closest_match(unlist(missing), headers,
                       otherwise = "Check your data headers.")
  )
}

#' Find the rules naming a value the data never carries.
#' @param known `list` Values the data carries, named after its columns.
#' @param cols `character()` Data column each column of `tbl` reads,
#' named after the columns of `tbl`.
#' @return `logical()` One flag per row of `tbl`.
#' @noRd
unknown_strata_rows <- function(tbl, known, cols) {
  unknown <- rep(FALSE, nrow(tbl))
  for (col in intersect(names(cols), colnames(tbl))) {
    # catch-all cell -> keys on nothing, names no unknown value
    vals <- as.character(tbl[[col]])
    unknown <- unknown |
      (!is.na(vals) & !vals %in% as.character(known[[cols[[col]]]]))
  }
  unknown
}

#' Assign a `rule_group` to each rule of a logic table.
#' A group is identified by a maximal pattern of non-NA covariates.
#' * rules with the same maximal pattern -> one group;
#' * rule fitting several equally specific groups -> ambiguous, error.
#' @param cov_cols `character()` Columns a rule keys on.
#' @param value_col `character(1)` Column holding what a rule applies.
#' @return `data.frame` Input with an integer `rule_group` column.
#' @noRd
compute_rule_groups <- function(tbl, cov_cols, value_col) {
  if (!nrow(tbl)) {
    tbl$rule_group <- integer(0)
    return(tbl)
  }

  # pattern of each row: names of non-NA covariate columns
  patterns <- lapply(seq_len(nrow(tbl)), function(i) {
    sort(cov_cols[!is.na(unlist(tbl[i, cov_cols, drop = FALSE]))])
  })

  # unique patterns for maximality check
  uniq_patterns <- unique(patterns)

  # a pattern is maximal when no other pattern strictly contains it
  has_strict_superset <- function(p) {
    any(vapply(uniq_patterns, \(q) length(q) > length(p) && all(p %in% q),
               logical(1)))
  }
  max_patterns <- uniq_patterns[
    !vapply(uniq_patterns, has_strict_superset, logical(1))
  ]

  # deterministic ordering: by serialized column set
  ord <- order(vapply(
    max_patterns, \(p) paste(p, collapse = key_sep_), character(1)
  ))
  max_patterns <- max_patterns[ord]

  # assign each rule to exactly one maximal pattern
  gids <- integer(nrow(tbl))
  for (i in seq_along(patterns)) {
    p <- patterns[[i]]
    candidates <- which(vapply(max_patterns, \(m) all(p %in% m), logical(1)))
    if (length(candidates) == 1) {
      gids[i] <- candidates
    } else {
      # nolint next: object_usage_linter. read by the cli glue below
      groups <- vapply(max_patterns[candidates], function(m) {
        cli::format_inline("[{.field {m}}]")
      }, character(1))
      stop_h(
        paste0("Ambiguous rule {i} ({.val {tbl[[value_col]][i]}}) ",
               if (!length(p)) "names no covariates"
               else "names only {.field {p}}",
               ", fitting {length(candidates)} equally specific groups: ",
               "{groups}."),
        "Add a disambiguating column to the rule."
      )
    }
  }
  tbl$rule_group <- gids
  tbl
}

#' Build the member cells a strata table is checked against.
#' Members holding the same covariate values resolve alike,
#' so each distinct combination stands once, weighted by its members.
#' * table read on infected members only -> checked at every infected status
#'   a member may take, whatever the data records.
#' @param covs `character()` Data columns the rules read, status excluded.
#' @param status_levels `character()` Infected-status labels,
#' empty for a table that does not read the status.
#' @return `data.frame` Columns `covs`, `status` when stratified,
#' `.weight` and the `.combo` index of each distinct combination.
#' @noRd
member_cells_ <- function(data_wide, covs, status_levels = character()) {
  if (!length(covs)) {
    cells <- data.frame(.weight = nrow(data_wide))
  } else {
    key <- row_keys_(data_wide, covs)
    first <- !duplicated(key)
    cells <- data_wide[first, covs, drop = FALSE]
    cells$.weight <- as.integer(table(key)[key[first]])
  }
  cells$.combo <- seq_len(nrow(cells))
  if (!length(status_levels)) return(cells)
  do.call(rbind, lapply(status_levels, function(lv) {
    cells$status <- rep(lv, nrow(cells))
    cells
  }))
}

#' Match every observed cell against every rule of a logic table.
#' @param cells `data.frame` Observed cells, columns named as the rule columns.
#' @return `matrix` Logical, one row per cell and one column per rule.
#' @noRd
rule_hits_ <- function(tbl, cells) {
  # * catch-all rule cell -> constrains nothing;
  # * missing data cell -> matches a catch-all only.
  Reduce(`&`, lapply(names(cells), function(col) {
    outer(as.character(cells[[col]]), as.character(tbl[[col]]),
          \(d, r) is.na(r) | (!is.na(d) & d == r))
  }), matrix(TRUE, nrow(cells), nrow(tbl)))
}

#' Find the observed cells some rule group of a logic table leaves unmatched.
#' @param cells `data.frame` Observed cells, columns named as the rule columns.
#' @param gid `integer()` Rule group of every rule.
#' @return `logical()` One flag per row of `cells`.
#' @noRd
unmatched_cells_ <- function(tbl, cells, gid) {
  # a table without rules matches no cell
  if (!nrow(tbl)) return(rep(TRUE, nrow(cells)))
  hit <- rule_hits_(tbl, cells)
  Reduce(`|`, lapply(unique(gid), function(g) {
    rowSums(hit[, gid == g, drop = FALSE]) == 0
  }))
}

#' Resolve every cell against the rules of a logic table.
#' Within each rule group, a cell goes to the most specific rule it matches,
#' as in the engine.
#' * rules tied at that specificity -> all win it.
#' @param cells `data.frame` Observed cells, columns named as the rule columns.
#' @param gid `integer()` Rule group of every rule.
#' @return `matrix` Logical, one row per cell and one column per rule.
#' @noRd
rule_winners_ <- function(tbl, cells, gid) {
  hit <- rule_hits_(tbl, cells)
  # 0 scores no match, otherwise one plus the rule specificity
  spec <- rowSums(!is.na(tbl[names(cells)]))
  score <- hit * rep(spec + 1L, each = nrow(cells))
  win <- hit
  for (g in unique(gid)) {
    in_g <- gid == g
    sub <- score[, in_g, drop = FALSE]
    top <- sub[cbind(seq_len(nrow(sub)), max.col(sub, ties.method = "first"))]
    win[, in_g] <- sub == top & top > 0
  }
  win
}

#' Count the units every rule of a logic table is resolved to.
#' A unit counts once per rule,
#' at whichever levels of its covariate combination the rule wins.
#' * rules tied at the top specificity -> share the cell.
#' @param tbl `data.frame` Logic table with its `rule_group` column.
#' @param cells `data.frame` Observed cells,
#' with the `.weight` units each holds
#' and the `.combo` index of its covariate combination.
#' @param reserved `character()` Columns of `tbl` no rule keys on.
#' @return `integer()` One count per rule.
#' @noRd
rule_counts_ <- function(tbl, cells, reserved) {
  cols <- strata_cov_cols_(tbl, reserved)
  win <- rule_winners_(tbl, cells[cols], tbl[["rule_group"]])
  by_combo <- rowsum(win * 1L, cells$.combo) > 0
  held <- cells$.weight[match(as.integer(rownames(by_combo)), cells$.combo)]
  as.integer(colSums(by_combo * held))
}

#' Collapse a named reference into a rendered phrase of `name = value` pairs.
#' @param reference `list()` Covariate values, named after the covariates.
#' @param strip `character(1)` Prefix removed from the names before printing.
#' @return `character(1)` Plain text, to be interpolated as a value.
#' @noRd
reference_to_phrase <- function(reference, strip = "^") {
  cli::ansi_collapse(mapply(function(nm, val) {
    cli::format_inline("{.field {nm}} = {.val {val}}")
  }, sub(strip, "", names(reference)), unlist(reference)))
}

# Render rule `i` of logic table `tbl` as its keyed cells and its value.
rule_phrase_ <- function(i, tbl, cov_cols, value_col) {
  cells <- as.list(tbl[i, cov_cols, drop = FALSE])
  keyed <- cells[!is.na(unlist(cells))]
  cli::format_inline(
    if (!length(keyed)) "catch-all" else reference_to_phrase(keyed),
    " ({.val {tbl[[value_col]][i]}})"
  )
}

# Stop when a strata table leaves `n` units unmatched.
# @param hint `character(1)` Rendered advice replacing the covering-rule one.
stop_unmatched_ <- function(n, label, unit, hint = NULL) {
  if (n == 0) return(invisible())
  if (is.null(hint)) hint <- cli::format_inline(
    "Add a rule covering them: an {.code NA} rule cell matches any value."
  )
  stop_h(
    "{label} leaves {n} {unit}{cli::qty(n)}{?s} unmatched by its rules.",
    "{hint}"
  )
}

#' Stop on the rules of a logic table that no unit resolves to.
#' @param dead `logical()` One flag per row of `tbl`.
#' @param cov_cols `character()` Columns a rule keys on.
#' @param value_col `character(1)` Column holding what a rule applies.
#' @param label `character(1)` The process, as it opens the message.
#' @param unit `character(1)` What a rule resolves: member or ordered pair.
#' @return `data.frame` Input, invisibly.
#' @noRd
assert_no_dead_rules <- function(tbl, dead, cov_cols, value_col, label, unit) {
  if (!any(dead)) return(invisible(tbl))
  # nolint next: object_usage_linter. read by the cli glue below
  rules <- vapply(which(dead), rule_phrase_, character(1),
                  tbl = tbl, cov_cols = cov_cols, value_col = value_col)
  stop_h(
    paste("{label} has {length(rules)} rule{?s} no {unit}",
          "resolves to: {rules}."),
    "Remove the rule, or check its values against the data."
  )
}

#' Stop on logic tables that do not resolve every unit to one rule per group,
#' checked over `cells`:
#' * every cell matches a rule of every rule group;
#' * every rule is the best match of some cell;
#' * no two rules of one group tie at a cell.
#' @param tables `list` of logic tables, each with its `rule_group` column.
#' @param cells `data.frame` Observed cells,
#' with the `.weight` units each holds
#' and the `.combo` index of its covariate combination.
#' @param reserved `character()` Columns of the tables no rule keys on.
#' @param value_col `character(1)` Column holding what a rule applies.
#' @param label `character(1)` The process, as it opens the message.
#' @param unit `character(1)` What a rule resolves: member or ordered pair.
#' @param hint `character(1)` Rendered advice on unmatched units,
#' `NULL` for the one to add a covering rule.
#' @return `list` Winner matrices from rule_winners_(), invisibly.
#' @noRd
assert_strata_fit_ <- function(tables, cells, reserved, value_col,
                               label, unit, hint = NULL) {
  cov_cols <- lapply(tables, strata_cov_cols_, reserved = reserved)
  gids <- lapply(tables, `[[`, "rule_group")
  # a combination is missed when one of its levels is
  missed <- Reduce(`|`, Map(function(tbl, cols, gid) {
    unmatched_cells_(tbl, cells[cols], gid)
  }, tables, cov_cols, gids))
  missed <- unique(cells$.combo[missed])
  stop_unmatched_(sum(cells$.weight[match(missed, cells$.combo)]), label,
                  unit, hint)

  wins <- Map(function(tbl, cols, gid) {
    rule_winners_(tbl, cells[cols], gid)
  }, tables, cov_cols, gids)
  for (k in seq_along(tables)) {
    tbl <- tables[[k]]
    win <- wins[[k]]
    gid <- gids[[k]]
    assert_no_dead_rules(tbl, colSums(win) == 0, cov_cols[[k]], value_col,
                         label, unit)
    # two rules of one group winning the same cell
    for (g in unique(gid)) {
      tied <- which(rowSums(win[, gid == g, drop = FALSE]) > 1)
      if (!length(tied)) next
      # nolint next: object_usage_linter. read by the cli glue below
      rules <- vapply(which(win[tied[1], ] & gid == g), rule_phrase_,
                      character(1), tbl = tbl, cov_cols = cov_cols[[k]],
                      value_col = value_col)
      stop_h("{label} has tied rules: {rules}.",
             "Add a rule more specific than both.")
    }
  }
  invisible(wins)
}

# Code rule columns as integers over the factor levels held in `info`,
# `key` mapping a column to its `info` entry.
# * logical column -> FALSE to 0, TRUE to 1;
# * column without levels -> as is;
# * NA rule cell -> catch-all, coded NA: the engine skips that covariate,
#   and it counts for no specificity.
encode_strata_levels_ <- function(defs, info, cols, key = identity) {
  for (col in cols) {
    entry <- info[[key(col)]]
    if (identical(entry$class, "logical")) {
      defs[[col]] <- as.integer(as.logical(defs[[col]]))
    } else if (!is.null(entry$levels)) {
      defs[[col]] <- as.integer(factor(defs[[col]], levels = entry$levels))
    }
  }
  defs
}

# Code the labelled statuses of logic table `tbl` against the engine codes.
code_status_col_ <- function(tbl, status_codes) {
  status <- tbl[["status"]]
  if (is.null(status) || is.numeric(status)) return(tbl)
  tbl[["status"]] <- unname(status_codes[as.character(status)])
  tbl
}
