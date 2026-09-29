## housepi - Epidemiological Analysis of Household Outbreaks.
##
## Copyright © 2024-2026 Institut Pasteur (Paris), CNRS and INSERM.
## Distributed under the terms of the GNU General Public License (GPLv3).
## See the LICENSE file for details.
##

# Reserved column names (non-covariate) in a logic table.
data_reserved_cols <- c("id", "hhid", "status")
# `value` pins a stratum during simulation, `const` fixes it in inference too.
risk_pin_cols <- c("value", "const")
risk_reserved_cols <- c("group_risk", "param_name", "n", "rule_group",
                        risk_pin_cols)

# Covariate columns of logic table `tbl`.
risk_cov_cols_ <- function(tbl) setdiff(colnames(tbl), risk_reserved_cols)

# Empty and typed prototype of every reserved column.
risk_reserved_proto <- list(group_risk = character(), param_name = character(),
                            n = integer(), rule_group = integer(),
                            value = double(), const = logical())

# Empty logic table with the reserved columns `cols`, in order.
empty_risk_table <- function(cols) {
  as.data.frame(risk_reserved_proto[cols])
}

# Typed `NA` padding for each of `cols`. A covariate column pads as character.
risk_na_pad <- function(cols) {
  lapply(stats::setNames(cols, cols), function(col) {
    proto <- risk_reserved_proto[[col]]
    if (is.null(proto)) NA_character_ else proto[NA_integer_]
  })
}

# Parameter definitions of the logic tables of one risk. A named parameter
# takes its override when given, else the shared `setup` or, missing that,
# the default for an `absolute` risk; its table pins its values.
risk_params_ <- function(tables, overrides, setup, absolute) {
  if (is.null(setup)) setup <- default_risk_mcmc_setup(absolute)
  out <- list()
  for (tbl in tables) {
    pins <- collect_table_pins(tbl)
    for (par_id in setdiff(tbl$param_name, ".reference")) {
      def <- if (!is.null(overrides[[par_id]])) overrides[[par_id]] else setup
      out[[par_id]] <- helper_pin_param_def(def, pins[[par_id]])
    }
  }
  out
}

# Whether a risk carries absolute values, some rule group having no
# reference. An automatic definition not yet built from the data reads it
# off its reference.
risk_is_absolute_ <- function(def, tables) {
  if (length(def) == 0) return(FALSE)
  if (identical(def$mode, "auto") && length(tables) == 0)
    return(is.null(def$reference))
  lacks_reference_(tables)
}

# Code rule columns as integers over the factor levels held in `info`.
# `key` maps a column to its `info` entry; columns without one stay as they
# are. An NA rule cell is a catch-all: it codes to NA, so the engine skips
# that covariate and it counts for no specificity.
encode_risk_levels_ <- function(defs, info, key = identity,
                                cols = colnames(defs)) {
  for (col in cols) {
    lv <- info[[key(col)]]$levels
    if (!is.null(lv)) {
      defs[[col]] <- as.integer(factor(defs[[col]], levels = lv))
    }
  }
  defs
}

# Reorder `df` to `cols`, padding the missing ones with typed `NA`.
pad_risk_cols_ <- function(df, cols) {
  miss <- setdiff(cols, names(df))
  if (length(miss) > 0) df[miss] <- risk_na_pad(miss)
  df[cols]
}

# Parameter-name prefix: empty, or the given prefix closed by an underscore.
format_risk_prefix_ <- function(prefix) {
  if (identical(prefix, "") || !is.character(prefix)) ""
  else paste0(prefix, "_")
}

# @description
# Tidy risk table for the read surface. Definitions expand from data.
# If none is attached the table is an empty shell with columns `cols`.
helper_show_risk_table <- function(module, cols, is_data_attached) {
  out <- module$tidy_table
  if (!is.null(out) && nrow(out) > 0) return(out)
  if (length(module$defs) > 0 && !is_data_attached) warn_h(
    "Cannot expand {module$family} definitions without data.",
    "Attach household data using {.fn add_data_wide}."
  )
  empty_risk_table(cols)
}

# Mode of a risk definition: automatic through `by`, custom through `strata`.
risk_mode_ <- function(by, strata, params) {
  hint <- paste("Choose {.arg by} for automatic definitions,",
                "{.arg strata} for custom.")
  is_auto <- !all(is.na(by))
  if (is_auto && !is.null(strata)) stop_h(
    "Impossible to provide at the same time {.arg {c('by', 'strata')}}.", hint
  )
  if (is_auto) {
    warn_unused_auto_params(params)
    return("auto")
  }
  if (is.null(strata)) stop_h(
    "Invalid risk definition: provide one between {.arg {c('by', 'strata')}}.",
    hint
  )
  "custom"
}

# Validate a custom strata table, with the `params` overriding its strata.
assert_valid_strata_ <- function(strata, params, reference,
                                 allow = character(), symmetric = FALSE) {
  assert_valid_logic_table(strata, allow)
  if (length(reference) > 0) stop_h(
    paste("Argument {.arg reference} is not used",
          "with a custom {.arg strata} table."),
    "Tag one row's {.field param_name} as {.val .reference} instead."
  )
  assert_one_reference_per_group(strata, symmetric)
  if (length(params) > 0) {
    assert_known_params_names(strata, params)
    assert_no_pin_params_clash(strata, params)
  }
  invisible(strata)
}

# Stop when a risk of definition `def` leaves `n` units unmatched.
stop_unmatched_ <- function(n, label, unit, def) {
  if (n == 0) return(invisible())
  stop_h(
    "{label} leaves {n} {unit}{cli::qty(n)}{?s} unmatched by its rules.",
    if (identical(def$mode, "auto")) {
      paste("Some {.field {def$by}} values are missing: change them,",
            "or pass a {.arg strata} table with a catch-all rule.")
    } else {
      "Add a rule covering them: an {.code NA} rule cell matches any value."
    }
  )
}

#' Assert that covariates choice is valid.
#' @param covariates Covariates as list of data column names.
#' @param allow Reserved names the caller accepts as covariates.
#' @keywords internal
#' @noRd
assert_valid_covariates <- function(covariates, allow = character()) {
  if (all(is.na(covariates))) return(invisible())

  if (isFALSE(is.character(covariates))) stop_h(
    "Argument {.arg by} must be of type character.",
    format_type_error(covariates)
  )

  # avoid group risks associated to reserved data/definition names
  reserved_data <- setdiff(c(data_reserved_cols, risk_reserved_cols), allow)
  forbidden <- reserved_data[reserved_data %in% covariates]
  if (length(forbidden) > 0) stop_h(
    "Impossible to associate risk to covariate{?s}: {.field {forbidden}}.",
    "These names are reserved.",
    paste("Change {cli::qty(forbidden)}covariate{?s} and/or",
          "corresponding data column{?s}.")
  )
}

#' Assert that covariates are present in data.
#' @param covariates Covariates as list of data column names.
#' @param headers All column names of the data.
#' @keywords internal
#' @noRd
assert_data_covariates <- function(covariates, headers) {
  missing <- setdiff(covariates, headers)
  if (length(missing) > 0) stop_h(
    "Unknown covariate{?s}: {.field {missing}}.",
    "Check your data headers."
  )
}

#' Assert that covariate values in data are safe.
#' @param this_table Logic table whose covariate columns are checked.
#' @return `character()`
#' @keywords internal
#' @noRd
assess_covariate_values <- function(this_table) {
  # column names to consider
  colmask <- risk_cov_cols_(this_table)
  for (col in colmask) {
    sapply(this_table[[col]], function(x) assert_valid_covariate(x, col))
  }
  colmask
}

#' Find the rules naming a value the data never carries.
#' @param tbl `data.frame` logic table with a `param_name` column.
#' @param known `list` Values the data carries, named after its columns.
#' @param cols `character()` Data column every column of `tbl` reads, named
#' after the columns of `tbl`.
#' @return `logical()` One flag per row of `tbl`.
#' @keywords internal
#' @noRd
unknown_strata_rows <- function(tbl, known, cols) {
  unknown <- rep(FALSE, nrow(tbl))
  for (col in intersect(names(cols), colnames(tbl))) {
    # a catch-all cell keys on nothing: it names no value to be unknown
    vals <- as.character(tbl[[col]])
    unknown <- unknown |
      (!is.na(vals) & !vals %in% as.character(known[[cols[[col]]]]))
  }
  unknown
}

#' Stop on the rules of a logic table that no unit resolves to.
#' @param tbl `data.frame` logic table with a `param_name` column.
#' @param dead `logical()` One flag per row of `tbl`.
#' @param label `character(1)` The risk, as it opens the message.
#' @param unit `character(1)` What a rule resolves: member or ordered pair.
#' @return `data.frame` the input, invisibly.
#' @keywords internal
#' @noRd
assert_no_dead_rules <- function(tbl, dead, label, unit) {
  if (!any(dead)) return(invisible(tbl))
  cov_cols <- risk_cov_cols_(tbl)
  # nolint next: object_usage_linter. read by the cli glue below
  rules <- vapply(which(dead), function(i) {
    cells <- as.list(tbl[i, cov_cols, drop = FALSE])
    keyed <- cells[!is.na(unlist(cells))]
    cli::format_inline(
      if (length(keyed) == 0) "catch-all" else reference_to_phrase(keyed),
      " ({.val {tbl$param_name[i]}})"
    )
  }, character(1))
  stop_h(
    paste("{label} has {length(rules)} rule{?s} no {unit}",
          "resolves to: {rules}."),
    "Remove the rule, or check its values against the data."
  )
}

#' Assert that a user-provided logic table is valid.
#' @param logic_table `data.frame` with a `param_name` column.
#' @param allow Reserved names the caller accepts as covariates.
#' @keywords internal
#' @noRd
assert_valid_logic_table <- function(logic_table, allow = character()) {
  if (!is.data.frame(logic_table)) stop_h(
    "Argument {.arg logic_table} must be a data.frame."
  )
  if (!"param_name" %in% colnames(logic_table)) stop_h(
    "Logic table must contain a {.field param_name} column."
  )
  if (nrow(logic_table) == 0) stop_h(
    "Logic table must contain at least one row."
  )
  # every row names its parameter, or tags the reference
  pnames <- as.character(logic_table$param_name)
  unnamed <- which(is.na(pnames) | !nzchar(pnames))
  if (length(unnamed) > 0) stop_h(
    paste("Strata {cli::qty(length(unnamed))}row{?s} without",
          "{.field param_name}: {unnamed}."),
    "Name the parameter of each row, or tag it {.val .reference}."
  )
  for (pname in setdiff(unique(pnames), ".reference")) {
    assert_valid_definition_name(pname, "parameter")
  }
  # validate covariate column names
  covariate_cols <- risk_cov_cols_(logic_table)
  assert_valid_covariates(covariate_cols, allow)
  # no covariate column, or only catch-all cells, duplicates the baseline
  if (all(is.na(logic_table[covariate_cols]))) stop_h(
    paste("Strata table sets no covariate value:",
          "its rules cannot tell members apart."),
    paste("Set a covariate value in a rule, or drop the risk:",
          "the baseline already treats every member alike.")
  )
  assert_valid_pins(logic_table)
}

#' Assert that the pin columns of a logic table are valid.
#' @param tbl `data.frame` logic table with a `param_name` column.
#' @return `data.frame` the input, invisibly.
#' @keywords internal
#' @noRd
assert_valid_pins <- function(tbl) {
  if (!any(risk_pin_cols %in% colnames(tbl))) return(invisible(tbl))
  # an all-NA column pins nothing, whichever type it was written with
  if (!all(is.na(tbl[["value"]])) && !is.numeric(tbl[["value"]])) stop_h(
    "Column {.field value} of the strata table must be numeric.",
    format_type_error(tbl[["value"]])
  )
  if (!all(is.na(tbl[["const"]])) && !is.logical(tbl[["const"]])) stop_h(
    "Column {.field const} of the strata table must be logical.",
    format_type_error(tbl[["const"]])
  )
  pins <- pad_risk_cols_(tbl, risk_pin_cols)
  bad <- tbl$param_name == ".reference" &
    (!is.na(pins$value) | !is.na(pins$const))
  if (any(bad)) stop_h(
    paste("The {.val .reference} is pinned at",
          "{cli::qty(sum(bad))}row{?s} {which(bad)}."),
    "The reference stratum is 1 by definition: it carries no pin."
  )
  bad <- pins$const %in% TRUE & is.na(pins$value)
  if (any(bad)) stop_h(
    paste("No {.field value} for fixed",
          "strat{?um/a} {.val {unique(tbl$param_name[bad])}}."),
    "A {.code const = TRUE} stratum needs the value it is fixed at."
  )
  # a stratum spread over several rules carries a single pin, whichever of
  # its rules the data keeps
  rows <- Filter(function(i) nrow(unique(pins[i, ])) > 1,
                 split(seq_len(nrow(tbl)), tbl$param_name))
  conflicts <- vapply(names(rows), function(nm) {
    cli::format_inline("{.val {nm}} (rows {rows[[nm]]})")
  }, character(1))
  if (length(conflicts) > 0) stop_h(
    "Conflicting {.field value} for {cli::qty(rows)}strat{?um/a} {conflicts}.",
    "A stratum spread over several rules carries a single {.field value}."
  )
  invisible(tbl)
}

#' @title Assign rule_group to each rule of logic table.
#' @description Rules belong to specificity groups.
#' A group is identified by the unique `maximal` non-NA covariates.
#' Rules with the same maximal pattern share a group.
#' Rules that are subset of multiple groups are ambiguous.
#' @param tbl `data.frame` logic table.
#' @return `data.frame` with integer `rule_group` column appended.
#' @keywords internal
#' @noRd
compute_rule_groups <- function(tbl) {
  if (nrow(tbl) == 0) {
    tbl$rule_group <- integer(0)
    return(tbl)
  }
  cov_cols <- risk_cov_cols_(tbl)

  # pattern of each row: names of non-NA covariate columns
  patterns <- lapply(seq_len(nrow(tbl)), function(i) {
    sort(cov_cols[!is.na(unlist(tbl[i, cov_cols, drop = FALSE]))])
  })

  # unique patterns for maximality check
  uniq_patterns <- unique(patterns)

  # a pattern is maximal when no other pattern strictly contains it
  has_strict_superset <- function(p) {
    any(vapply(uniq_patterns, function(q) {
      length(q) > length(p) && all(p %in% q)
    }, logical(1)))
  }
  max_patterns <- uniq_patterns[
    !vapply(uniq_patterns, has_strict_superset, logical(1))
  ]

  # deterministic ordering: by serialized column set
  ord <- order(vapply(
    max_patterns, function(p) paste(p, collapse = key_sep_), character(1)
  ))
  max_patterns <- max_patterns[ord]

  # assign each rule to exactly one maximal pattern
  gids <- integer(nrow(tbl))
  for (i in seq_along(patterns)) {
    p <- patterns[[i]]
    candidates <- which(vapply(max_patterns,
                               function(m) all(p %in% m), logical(1)))
    if (length(candidates) == 1) {
      gids[i] <- candidates
    } else {
      # nolint next: object_usage_linter. read by the cli glue below
      groups <- vapply(max_patterns[candidates], function(m) {
        cli::format_inline("[{.field {m}}]")
      }, character(1))
      stop_h(
        paste0("Ambiguous rule {i} ({.val {tbl$param_name[i]}}) ",
               if (length(p) == 0) "names no covariates"
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

#' Match every observed cell against every rule of a logic table.
#' @param tbl `data.frame` logic table.
#' @param cells `data.frame` observed cells, columns named as the rule columns.
#' @return `logical` matrix, one row per cell and one column per rule.
#' @keywords internal
#' @noRd
rule_hits_ <- function(tbl, cells) {
  # a catch-all rule cell constrains nothing; a missing data cell matches
  # nothing but a catch-all
  Reduce(`&`, lapply(names(cells), function(col) {
    outer(as.character(cells[[col]]), as.character(tbl[[col]]),
          function(d, r) is.na(r) | (!is.na(d) & d == r))
  }), matrix(TRUE, nrow(cells), nrow(tbl)))
}

#' Build the covariate cells a group risk is checked against.
#' @description Members holding the same covariate values resolve alike, so
#' each distinct combination stands once, weighted by the members holding it.
#' A risk read on infected members only is checked at every infected status a
#' member may take, whatever the data records.
#' @param data_wide `data.frame` Data the members come from.
#' @param covs `character()` Data columns the rules read, status excluded.
#' @param status_levels `character()` Infected-status labels, empty for a
#' risk that does not read the status.
#' @return `data.frame` with columns `covs`, `status` when stratified, the
#' `.weight` and the `.combo` index of each distinct combination.
#' @keywords internal
#' @noRd
risk_check_cells_ <- function(data_wide, covs, status_levels = character()) {
  if (length(covs) == 0) {
    cells <- data.frame(.weight = nrow(data_wide))
  } else {
    key <- row_keys_(data_wide, covs)
    first <- !duplicated(key)
    cells <- data_wide[first, covs, drop = FALSE]
    cells$.weight <- as.integer(table(key)[key[first]])
  }
  cells$.combo <- seq_len(nrow(cells))
  if (length(status_levels) == 0) return(cells)
  do.call(rbind, lapply(status_levels, function(lv) {
    cells$status <- rep(lv, nrow(cells))
    cells
  }))
}

#' Find the observed cells some rule group of a logic table leaves unmatched.
#' @param tbl `data.frame` logic table.
#' @param cells `data.frame` observed cells, columns named as the rule columns.
#' @return `logical()` One flag per row of `cells`.
#' @keywords internal
#' @noRd
unmatched_cells_ <- function(tbl, cells) {
  # a table without rules matches no cell
  if (nrow(tbl) == 0) return(rep(TRUE, nrow(cells)))
  hit <- rule_hits_(tbl, cells)
  gid <- compute_rule_groups(tbl)$rule_group
  Reduce(`|`, lapply(unique(gid), function(g) {
    rowSums(hit[, gid == g, drop = FALSE]) == 0
  }))
}

#' Count the units every rule of a logic table is resolved to.
#' @description Within each rule group, a cell goes to the most specific rule
#' it matches, as in the engine. Rules tied at that specificity share the cell.
#' @param tbl `data.frame` logic table.
#' @param cells `data.frame` observed cells, columns named as the rule columns.
#' @param weight `integer()` units each cell holds.
#' @return `integer()` One count per row of `tbl`.
#' @keywords internal
#' @noRd
count_resolved_rules <- function(tbl, cells, weight = 1L) {
  hit <- rule_hits_(tbl, cells)
  # 0 scores no match, otherwise one plus the rule specificity
  spec <- rowSums(!is.na(tbl[names(cells)]))
  score <- hit * rep(spec + 1L, each = nrow(cells))
  gid <- compute_rule_groups(tbl)$rule_group
  win <- hit
  for (g in unique(gid)) {
    in_g <- gid == g
    sub <- score[, in_g, drop = FALSE]
    top <- sub[cbind(seq_len(nrow(sub)), max.col(sub, ties.method = "first"))]
    win[, in_g] <- sub == top & top > 0
  }
  as.integer(colSums(win * weight))
}

#' @title Assert that each rule group designates at most one reference cell.
#' @description A reference sets the unit scale of its rule group,
#' so two distinct cells of one group cannot both carry the tag.
#' Symmetric contact rules reach one cell through either direction,
#' hence both orderings of a pair count as the same cell.
#' @param strata `data.frame` logic table with a `param_name` column.
#' @param symmetric `logical(1)` Treat reversed from/to pairs as one cell.
#' @return `data.frame` the input, invisibly.
#' @keywords internal
#' @noRd
assert_one_reference_per_group <- function(strata, symmetric = FALSE) {
  idx_refer <- which(strata$param_name == ".reference")
  if (length(idx_refer) < 2) return(invisible(strata))

  tbl <- compute_rule_groups(strata)[idx_refer, , drop = FALSE]
  cov_cols <- risk_cov_cols_(tbl)
  covariates <- pair_covs_(cov_cols)
  pair_wise <- isTRUE(symmetric)

  # cell identity: the rule itself, or the unordered pair when symmetric
  cell_key <- function(i) {
    row <- vapply(tbl[i, cov_cols, drop = FALSE], as.character, character(1))
    if (!pair_wise) return(paste(cov_cols, row, sep = "_", collapse = key_sep_))
    ends <- lapply(c("from_", "to_"), function(side) {
      paste(covariates, row[paste0(side, covariates)],
            sep = "_", collapse = key_sep_)
    })
    paste(sort(unlist(ends)), collapse = "\x1E")
  }
  keys <- vapply(seq_along(idx_refer), cell_key, character(1))

  # within a rule group, every tagged row must denote the same cell
  clash <- unlist(lapply(
    split(seq_along(idx_refer), tbl$rule_group),
    function(k) if (length(unique(keys[k])) > 1) idx_refer[k]
  ))
  if (length(clash) > 0) stop_h(
    paste("{cli::qty(length(clash))}Row{?s} {sort(clash)} of one rule group",
          "are tagged {.val .reference}."),
    "A rule group has a single reference, setting its unit scale."
  )
  invisible(strata)
}

#' Pins carried by a logic table, keyed by `param_name`.
#' A stratum spread over several rules carries a single pin, read off its first.
#' @param tbl `data.frame` logic table with a `param_name` column.
#' @return `list` of `list(value = , const = )`, one per pinned stratum.
#' @keywords internal
#' @noRd
collect_table_pins <- function(tbl) {
  if (!any(risk_pin_cols %in% colnames(tbl))) return(list())
  pins <- pad_risk_cols_(tbl, risk_pin_cols)
  out <- list()
  for (par_id in setdiff(unique(tbl$param_name), ".reference")) {
    first <- match(par_id, tbl$param_name)
    pin <- list(value = pins$value[first], const = pins$const[first])
    if (!is.na(pin$value) || !is.na(pin$const)) out[[par_id]] <- pin
  }
  out
}

#' Fold a stratum pin into its parameter definition.
#' A fixed stratum is a constant parameter: it carries no prior.
#' @param def `list` Shared MCMC setup of the risk.
#' @param pin `list` or `NULL` Pin of the stratum.
#' @return `list`
#' @keywords internal
#' @noRd
helper_pin_param_def <- function(def, pin) {
  if (is.null(pin)) return(def)
  if (isTRUE(pin$const)) return(list(value = pin$value, const = TRUE))
  c(def, list(value = pin$value))
}

#' Whether a set of logic tables leaves the scale of its risk free.
#' Every unit matches a rule of every rule group, so a rule group without a
#' `.reference` row sets the scale on its own. A table without rules carries
#' no scale.
#' @param tables `list` of logic tables with a `param_name` column.
#' @return `logical(1)`
#' @keywords internal
#' @noRd
lacks_reference_ <- function(tables) {
  any(vapply(tables, function(tbl) {
    gid <- compute_rule_groups(tbl)$rule_group
    !all(gid %in% gid[tbl$param_name == ".reference"])
  }, logical(1)))
}

#' Assert every `params` entry names a parameter the strata table declares.
#' @param strata `data.frame` logic table with a `param_name` column.
#' @param params `list` Per-parameter MCMC-setup overrides.
#' @return `list` the input `params`, invisibly.
#' @keywords internal
#' @noRd
assert_known_params_names <- function(strata, params) {
  unknown <- setdiff(names(params), setdiff(strata$param_name, ".reference"))
  if (length(unknown) > 0) stop_h(
    "Unknown parameter{?s} in {.arg params}: {.val {unknown}}."
  )
  invisible(params)
}

#' Warn that `params` is dropped: an automatic definition derives its own
#' parameters, so the list has nothing to override.
#' @param params `list` Per-parameter MCMC-setup overrides.
#' @return `list` the input `params`, invisibly.
#' @keywords internal
#' @noRd
warn_unused_auto_params <- function(params) {
  if (length(params) == 0) return(invisible(params))
  warn_h(
    "{.arg params} is not used with {.arg by}.",
    paste("Give {.val {names(params)}} a value with {.fn add_parameter},",
          "or name the parameters in a {.arg strata} table.")
  )
  invisible(params)
}

#' Assert that no stratum is both pinned and overridden through `params`.
#' @param strata `data.frame` logic table with a `param_name` column.
#' @param params `list` Per-parameter MCMC-setup overrides.
#' @return `list` the input `params`, invisibly.
#' @keywords internal
#' @noRd
assert_no_pin_params_clash <- function(strata, params) {
  clash <- intersect(names(params), names(collect_table_pins(strata)))
  if (length(clash) > 0) stop_h(
    paste("Pinned {cli::qty(clash)}strat{?um/a} also set through",
          "{.arg params}: {.val {clash}}."),
    paste("A pinned stratum carries no MCMC setup:",
          "drop it from {.arg params} or from the strata table.")
  )
  invisible(params)
}

#' Set ordering for risk logic-tables: `.reference` first, decreasing 'n'.
#' @param tbl `data.frame` Logic table with `param_name`, `n`.
#' @return `data.frame`
#' @keywords internal
#' @noRd
reorder_risks_logic_table <- function(tbl) {
  if (nrow(tbl) == 0) return(tbl)
  is_ref <- tbl$param_name == ".reference"
  n_val <- suppressWarnings(as.numeric(tbl[["n"]]))
  ord <- order(tbl$rule_group, !is_ref, -n_val, method = "radix")
  tbl[ord, , drop = FALSE]
}

#' Read the data status column as infected-status labels.
#' Escaped and missing members carry `NA`: they match no status rule and read
#' the reference risk.
#' @param x Status column, integer codes or a factor of labels.
#' @param status_codes `integer()` Engine status codes, named after the labels.
#' @return `character()`
#' @keywords internal
#' @noRd
status_as_labels <- function(x, status_codes) {
  labels <- names(status_codes)
  if (is.factor(x) || is.character(x)) {
    x <- as.character(x)
    return(ifelse(x %in% labels, x, NA_character_))
  }
  labels[match(suppressWarnings(as.integer(x)), status_codes)]
}

#' Return logic table of group risks
#' @param data_wide Data from which to extract group risks
#' @param this_cov Covariate as list of data column names.
#' @param status_levels `character()` Declared infected-status labels.
#' @param status_observed `logical(1)` Whether the data carries any status.
#' @return data.frame with logic table
#' @keywords internal
#' @noRd
group_risk_logic_table <- function(data_wide, this_cov,
                                   status_levels = character(),
                                   status_observed = TRUE) {
  # extract unique covariates and sort by list length
  # status strata are declared by the model, not read off the data
  uni_covariates <- lapply(this_cov, function(f) {
    if (identical(f, "status")) status_levels else sort(unique(data_wide[[f]]))
  })
  sorted_indices <- order(sapply(uni_covariates, length), decreasing = TRUE)
  this_cov <- this_cov[sorted_indices]
  uni_covariates <- uni_covariates[sorted_indices]
  names(uni_covariates) <- this_cov

  # create a dataframe of all possible multi-covariates
  y <- expand.grid(uni_covariates)

  # count every combination in one pass. A cell the data does not carry
  # matches no stratum, so those members are left out of the count
  carried <- data_wide[stats::complete.cases(data_wide[this_cov]), ,
                       drop = FALSE]
  counted <- table(row_keys_(carried, this_cov))
  multicov_count <- as.integer(counted[row_keys_(y, this_cov)])
  multicov_count[is.na(multicov_count)] <- 0L

  # keep the combinations the data carries; a member may draw any infected
  # status, so every declared level stands on each of them
  keep <- multicov_count > 0
  if ("status" %in% this_cov) {
    others <- setdiff(this_cov, "status")
    keep <- rep(TRUE, nrow(y))
    if (length(others) > 0) {
      held <- data_wide[stats::complete.cases(data_wide[others]), ,
                        drop = FALSE]
      keep <- row_keys_(y, others) %in% row_keys_(held, others)
    }
  }
  y <- y[keep, , drop = FALSE]
  # no status is observed: sample sizes unknown
  y$n <- if ("status" %in% this_cov && isFALSE(status_observed))
    rep(NA_integer_, nrow(y)) else multicov_count[keep]

  y
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

#' Associate parameter name to group risk logic table.
#' @param this_table
#' Data.frame containing conditions for the group risk.
#' @param prefix
#' Prefix to parameter names associated to the group risk.
#' @param reference
#' Choose reference providing corresponding covariates.
#' `".auto"` picks the group with the largest sample size.
#' Empty leaves the risk without a reference.
#' @return
#' input data.frame with logic table with parameter names.
#' @keywords internal
#' @noRd
associate_group_param_names <- function(this_table, prefix, reference) {
  # column names to consider
  colmask <- assess_covariate_values(this_table)

  # find reference index
  idx_refer <- NULL
  if (identical(reference, ".auto")) {
    if (all(is.na(this_table$n))) stop_h(
      "Impossible to pick an automatic reference.",
      `!` = "Sample sizes are unknown while the data carries no status.",
      paste("Name the reference stratum, as",
            "{.code reference = list(status = \"<label>\")}.")
    )
    idx_refer <- which.max(this_table$n)

  } else if (length(reference) > 0) {
    idx_refer <- reference_rows_(this_table, colmask, reference)
    if (length(idx_refer) != 1) stop_h(
      "Provided reference group is not valid.",
      `!` = paste("Individuals with {reference_to_phrase(reference)}",
                  "are never observed.")
    )
  }

  prefix <- format_risk_prefix_(prefix)

  # associate parameter names
  this_table$param_name <- apply(
    this_table[colmask], 1, function(row) {
      paste0(prefix, paste(paste(colmask, row, sep = "_"), collapse = "_"))
    }
  )
  if (length(idx_refer) == 1) this_table$param_name[idx_refer] <- ".reference"

  this_table
}

# Rows of `tbl` whose columns `cols` hold the values of `reference`.
reference_rows_ <- function(tbl, cols, reference) {
  unname(which(apply(tbl[cols], 1, function(row) {
    all(row == unlist(reference[cols]))
  })))
}

#' Combine list with group risk definition into a single dataframe
#' @param all_tables `list`. Definitions for each group risk.
#' @return `data.frame`
#' @keywords internal
#' @noRd
tidyup_group_risk_tables <- function(all_tables) {
  df_list <- lapply(names(all_tables), function(risk) {
    this_table <- all_tables[[risk]]
    if (length(this_table) == 0) return(data.frame(
      group_risk = risk, param_name = ".reference", rule_group = 1L
    ))
    # concat this risk's sub-tables, then compute rule_groups within risk
    cols <- unique(c("param_name", "n",
                     sort(unlist(lapply(this_table, names)))))
    tbl <- compute_rule_groups(do.call(rbind, lapply(this_table, pad_risk_cols_,
                                                     cols = cols)))
    cbind(group_risk = risk, reorder_risks_logic_table(tbl))
  })

  # all unique column names reordered, pin columns only when a pin exists
  present <- unique(unlist(lapply(df_list, names)))
  all_cols <- unique(c(setdiff(risk_reserved_cols, risk_pin_cols),
                       intersect(risk_pin_cols, present), sort(present)))

  combined_df <- do.call(rbind, lapply(df_list, pad_risk_cols_,
                                       cols = all_cols))
  rownames(combined_df) <- NULL
  as.data.frame(combined_df)
}

#' Return the user MCMC setup of a risk, `NULL` if every choice is missing.
#' @param dots `list` Distribution and proposal arguments.
#' @param family `character(1)` Holder family.
#' @param name `character(1)` or `NULL`. Name of the risk, when it has one.
#' @return `list` or `NULL`
#' @keywords internal
#' @noRd
unwrap_risk_setup_ <- function(prior, proposal, min, max, dots, family,
                               name = NULL) {
  prior <- housepi_smart_unwrap(prior, "prior")
  proposal <- housepi_smart_unwrap(proposal, "proposal")
  min <- housepi_smart_unwrap(min, "min")
  max <- housepi_smart_unwrap(max, "max")
  def_args <- helper_definition_args(proposal, dots, family, name)
  args <- def_args$args
  proposal <- def_args$proposal
  assert_no_sibling_dist_arg(args, family, name)

  is_default <- sapply(c(prior, args, min, max, proposal),
                       function(x) is.na(x) || is.null(x))
  # all defaults: resolved once the data tells absolute from relative
  if (all(is_default)) return(NULL)
  # nolint next: object_usage_linter. read by the cli glue below
  label <- if (is.null(name)) family else name
  if (is_default[1]) stop_h(
    "Found spare prior arguments for risk {.val {label}}.",
    "Either leave all as defaults or provide a {.arg prior}."
  )
  c(list(prior = prior), args, list(min = min, max = max, proposal = proposal))
}

#' Return the default MCMC setup of a risk.
#' Absolute risks are rates on their own; relative ones are ratios,
#' whose defaults follow \insertRef{cauchemez_household_2009}{housepi}.
#' @param absolute `logical(1)` Whether the risk carries absolute values.
#' @return `list`
#' @keywords internal
#' @noRd
default_risk_mcmc_setup <- function(absolute) {
  if (isTRUE(absolute)) list(prior = "uniform", proposal = "rw_log")
  else list(prior = "lognormal", meanlog = 0, sdlog = 3,
            min = 0, max = 10, proposal = "rw_log")
}

# Validate and shape the automatic reference of group risks.
auto_group_reference_ <- function(by, reference) {
  if (identical(reference, ".auto")) return(".auto")

  if (length(by) != length(reference)) stop_h(
    "Arguments {.arg {c('by', 'reference')}} should have equal dimensions."
  )
  else if (length(setdiff(names(reference), by)) != 0) stop_h(
    "Unrecognized names in reference list."
  )
  reference <- as.list(reference)
  stats::setNames(reference, by)
}

############################################################################
############################################################################
############################# Housepi Group Risk ###########################
############################################################################
############################################################################

#' @title Internal R6 holder for housepi group-risk definitions
#' @description
#' Helper class managing group-risk definitions of the model object.
#' @keywords internal
#' @noRd
HousepiGroupRisk <- R6::R6Class(
  classname = "HousepiGroupRisk",
  cloneable = TRUE,

  ############################################################################
  ################################## PRIVATE #################################
  ############################################################################

  private = list(
    deep_clone = function(name, value) helper_clone(value),

    # Add risk definition in automatic mode.
    helper_auto_add_ = function(name, by,
                                prefix = NA_character_,
                                reference = NULL,
                                joint = FALSE,
                                .warn_over = TRUE) {
      assert_valid_covariates(by, allow = "status")
      # remove previous definition if it exists
      rmv_name <- self$remove(name)

      warn_or_store_overwrite_(self, rmv_name, .warn_over)

      # load reference and validate reference names
      if (length(reference) > 0) {
        reference <- auto_group_reference_(by, reference)
      }

      # add new definition
      self$defs[[name]] <- list(by = by, reference = reference,
                                prefix = prefix, joint = joint,
                                mode = "auto")
      invisible(self)
    },

    # Add risk definition in custom mode.
    helper_custom_add_ = function(name, strata,
                                  params = NULL,
                                  reference = NULL,
                                  .warn_over = TRUE) {
      assert_valid_strata_(strata, params, reference, allow = "status")

      # handle replace
      if (name %in% names(self$defs) &&
            self$defs[[name]][["mode"]] == "auto") {
        # remove previous definition
        self$logic_tables[[name]] <- NULL
        self$param_overrides[[name]] <- NULL
        rmv_name <- self$remove(name)
        warn_or_store_overwrite_(self, rmv_name, .warn_over)
      }
      self$logic_tables[[name]] <- list(strata)
      self$param_overrides[[name]] <- params

      # store a sentinel definition
      self$defs[[name]] <- list(by = risk_cov_cols_(strata),
                                prefix = NA_character_, reference = NULL,
                                joint = FALSE, mode = "custom")
      invisible(self)
    },

    # Stop if a rule group of one risk leaves members unmatched, or if a rule
    # is the best match of no member. A risk read on infected members only is
    # checked at every level in `status_levels`, whatever the data records.
    assert_rules_fit_ = function(risk_id, data_wide,
                                 status_levels = character()) {
      tables <- self$logic_tables[[risk_id]]
      rule_cols <- lapply(tables, risk_cov_cols_)
      covs <- setdiff(unique(unlist(rule_cols)), "status")
      cells <- risk_check_cells_(data_wide, covs, status_levels)

      # a combination is missed when one of its levels is
      missed <- Reduce(`|`, Map(function(tbl, cols) {
        unmatched_cells_(tbl, cells[cols])
      }, tables, rule_cols))
      missed <- unique(cells$.combo[missed])
      label <- cli::format_inline("Group risk {.val {risk_id}}")
      stop_unmatched_(sum(cells$.weight[match(missed, cells$.combo)]), label,
                      "member", self$defs[[risk_id]])

      for (k in seq_along(tables)) {
        wins <- count_resolved_rules(tables[[k]], cells[rule_cols[[k]]],
                                     cells$.weight)
        assert_no_dead_rules(tables[[k]], wins == 0, label, "member")
      }
      invisible(self)
    },

    # Whether one risk stratifies on the infected status.
    keys_on_status_ = function(risk_id) {
      "status" %in% self$defs[[risk_id]][["by"]]
    },

    # Implement through automatic mode.
    helper_implement_auto_ = function(risk_id, data_wide,
                                      status_levels, status_observed) {
      this_def <- self$defs[[risk_id]] # definition of the risk
      covariates <- as.list(this_def[["by"]]) # data covariates
      covariates <- sort(as.character(covariates[!is.na(covariates)]))
      # status is carried by the model, not required as a data column
      assert_data_covariates(setdiff(covariates, "status"), colnames(data_wide))

      # cast to list if joint
      if (isTRUE(this_def[["joint"]])) covariates <- list(covariates)

      # default prefix
      prefix <- as.character(this_def[["prefix"]])
      if (is.na(prefix)) prefix <- substring(tolower(risk_id), 1, 3)

      self$logic_tables[[risk_id]] <- lapply(covariates, function(this_cov) {
        this_table <- group_risk_logic_table(data_wide, this_cov,
                                             status_levels, status_observed)
        associate_group_param_names(this_table, prefix, this_def$reference)
      })
    },

    # Validate custom mode logic table against data, then add counts.
    helper_implement_custom_ = function(risk_id, data_wide, status_observed) {
      for (ti in seq_along(self$logic_tables[[risk_id]])) {
        tbl <- self$logic_tables[[risk_id]][[ti]]
        cov_cols <- risk_cov_cols_(tbl)
        # status is carried by the model, not required as a data column
        assert_data_covariates(setdiff(cov_cols, "status"), colnames(data_wide))
        tbl$n <- count_resolved_rules(tbl, data_wide[cov_cols])
        # no status is observed: every declared stratum stands, sizes unknown
        if ("status" %in% cov_cols && isFALSE(status_observed))
          tbl$n <- NA_integer_
        self$logic_tables[[risk_id]][[ti]] <- tbl
      }
      invisible(self)
    }
  ),

  ############################################################################
  ################################## ACTIVE ##################################
  ############################################################################

  active = list(
    #' @field parameters `list`.
    #' Parameter definitions associated to group risks.
    parameters = function(value) {
      new_params <- list()
      for (risk_id in names(self$logic_tables)) {
        par <- risk_params_(self$logic_tables[[risk_id]],
                            self$param_overrides[[risk_id]],
                            self$raw_params[[risk_id]],
                            self$is_absolute(risk_id))
        new_params[names(par)] <- par
      }
      new_params
    },

    #' @field data_required `character()`. Data required by group risks.
    data_required = function(value) {
      required <- unique(unlist(lapply(self$defs, `[[`, "by")))
      required <- setdiff(required[!is.na(required)], "status")
      set_source_name_(required, self$family)
    },

    #' @field covariate_names `character()`. Data columns the risks stratify on,
    #' declared through `by` or carried by a custom logic table.
    covariate_names = function(value) {
      setdiff(self$raw_covariate_names, "status")
    },

    #' @field raw_covariate_names `character()`. Stratification axes, `status`
    #' included. `covariate_names` drops it: it is an always-present model axis,
    #' neither a required data column nor a name an observation can collide on.
    raw_covariate_names = function(value) {
      custom <- lapply(unlist(self$logic_tables, recursive = FALSE),
                       risk_cov_cols_)
      covs <- as.character(unlist(c(lapply(self$defs, `[[`, "by"), custom)))
      unique(covs[!is.na(covs)])
    },

    #' @field param_required `character()`. Parameters required by group risks.
    param_required = function(value) {
      tables <- unlist(self$logic_tables, recursive = FALSE)
      required <- setdiff(unlist(lapply(tables, `[[`, "param_name")),
                          ".reference")
      set_source_name_(required, self$family)
    },

    #' @field tidy_table `data.frame()`.
    #' Combined list with group risk and parameter definitions.
    tidy_table = function(value) {
      tidyup_group_risk_tables(self$logic_tables)
    }
  ),

  ############################################################################
  ################################## PUBLIC ##################################
  ############################################################################

  public = list(
    #' @field family `character(1)`. Holder family.
    family = "group risk",

    #' @field store_warn_over `character()`. Cached overwritten names.
    store_warn_over = character(),

    #' @field defs `list`. Current group risk definitions.
    defs = list(),

    #' @field raw_params `list`. Raw parameter definitions.
    raw_params = list(),

    #' @field param_overrides `list`.
    #' Per-parameter MCMC overrides keyed by risk name.
    param_overrides = list(),

    #' @field logic_tables `list`.
    #' Table of parameters associated to each group risk.
    logic_tables = list(),

    #' @description
    #' New object of class [HousepiGroupRisk][housepi::HousepiGroupRisk].
    initialize = function() {
      self$defs <- list()
      self$logic_tables <- list()
      self$raw_params <- list()
      self$param_overrides <- list()
      invisible(self)
    },

    #' @description
    #' Remove a given group_risk definition.
    #' @param name `character(1)`. Name of group risk definition to remove.
    remove = function(name) {
      if (any(names(self$defs) == name)) {
        # drop previous definition
        self$defs[name] <- NULL
        self$logic_tables[name] <- NULL
        self$raw_params[name] <- NULL
        self$param_overrides[name] <- NULL
        return(name)
      }

      c()
    },

    #' @description
    #' Whether a risk carries absolute values, some rule group having no
    #' reference.
    #' @param name `character(1)`. Name of the group risk definition.
    #' @return `logical(1)`
    is_absolute = function(name) {
      risk_is_absolute_(self$defs[[name]], self$logic_tables[[name]])
    },

    #' @inherit add_group_risk
    add_quo = function(name, by = NA_character_, ...,
                       reference = NULL, joint = FALSE,
                       strata = NULL, params = NULL,
                       prior = NA_character_,
                       min = NA_real_, max = NA_real_,
                       proposal = NA_character_,
                       prefix = NA_character_, .warn_over = TRUE) {
      name <- housepi_smart_unwrap(name, "name")
      by <- housepi_smart_unwrap(by, "by")
      # control valid name # NOTE: `by` is the injected default variable
      assert_valid_definition_name(name, self$family)

      # automatic definition
      prefix <- housepi_smart_unwrap(prefix, "prefix")
      reference <- housepi_smart_unwrap(reference, "reference", null_ok = TRUE)
      joint <- isTRUE(housepi_smart_unwrap(joint, "joint"))

      # custom definition
      strata <- housepi_smart_unwrap(strata, "strata", null_ok = TRUE)
      params <- housepi_smart_unwrap(params, "params", null_ok = TRUE)

      setup <- unwrap_risk_setup_(prior, proposal, min, max, list(...),
                                  self$family, name)
      if (risk_mode_(by, strata, params) == "auto") {
        private$helper_auto_add_(
          name = name, by = by, reference = reference, prefix = prefix,
          joint = joint, .warn_over = .warn_over
        )
      } else {
        private$helper_custom_add_(
          name = name, strata = strata, params = params,
          reference = reference, .warn_over = .warn_over
        )
      }

      # store shared MCMC setup (overwriting any existing)
      self$raw_params[name] <- list(setup)

      # reorder
      y <- sort(names(self$defs))
      self$defs <- self$defs[y]
      self$raw_params <- self$raw_params[y]

      invisible(self)
    },

    #' @inherit add_group_risk
    add_quo_list = function(xquos, append = FALSE) {
      housepi_load_defin_quos(self, "by", xquos, append)
      invisible(self)
    },

    ############################################################################
    ############################################################################
    ############################################################################

    #' @description
    #' Implement group risk definition according to data.
    #' @param data_wide
    #' A `data.frame` with the dataset to be considered.
    #' @param status_codes
    #' `integer()` Engine status codes, named after the labels the risks may
    #' stratify on.
    #' @param wall `logical(1)`. Whether to raise error.
    #' @param stratified `character()` Risks the model reads on infected
    #' members only.
    implement = function(data_wide, status_codes = integer(),
                         wall = FALSE, stratified = character()) {
      # preserve custom-mode entries and drop auto-mode entries
      mask_custom <- names(which(vapply(
        self$defs, function(d) d[["mode"]] == "custom", logical(1)
      )))
      self$logic_tables <- self$logic_tables[mask_custom]

      # read the status through the dictionary: the exported column is a factor
      # of labels once the levels are known, and integer codes before that
      status_levels <- names(status_codes)
      if ("status" %in% colnames(data_wide))
        data_wide[["status"]] <- status_as_labels(data_wide[["status"]],
                                                  status_codes)
      status_observed <- any(!is.na(data_wide[["status"]]))

      for (risk_id in names(self$defs)) {
        # a status-keyed risk waits for the statuses to be declared, so
        # decorator order stays free; the wall reports one that never came
        if (length(status_codes) == 0 && private$keys_on_status_(risk_id)) {
          if (isTRUE(wall)) stop_h(
            paste("A {self$family} is keyed on {.field status},",
                  "for which the model declares no level."),
            "Declare the infected statuses, or key the risk on a data column."
          )
          self$logic_tables[[risk_id]] <- list()
          next
        }
        if (self$defs[[risk_id]][["mode"]] == "auto") {
          private$helper_implement_auto_(risk_id, data_wide,
                                         status_levels, status_observed)
        } else {
          private$helper_implement_custom_(risk_id, data_wide, status_observed)
        }
        # the engine refuses a status key on a risk it does not stratify
        if (risk_id %in% stratified) {
          private$assert_rules_fit_(risk_id, data_wide, status_levels)
        } else if (!private$keys_on_status_(risk_id)) {
          private$assert_rules_fit_(risk_id, data_wide)
        }
      }

      invisible(self)
    },

    #' @description
    #' Assert that group risk definitions match the requirements.
    #' @param required `character()`. Required group risks.
    assert = function(required) {
      # load definitions
      defined <- names(self$defs)

      # UPDATE
      # point out group risks that are not used
      notused <- template_not_used_find(required, defined)
      # drop unused group risks
      if (length(notused) > 0) {
        sms <- c("Dropped unused group risk{?s}: {.val {notused}}.",
                 find_closest_match(notused, required))
        for (which in notused) self$remove(which)
        stop_h(sms)
      }

      template_missing_stop(required, defined, self$family)
    },

    #' @description
    #' Tidy table of group risks definition to be passed to Cpp.
    #' @param info
    #' `list` of `list`s. Info about class and levels of data columns.
    #' @param status_codes
    #' `integer()` Engine status codes, named after the labels.
    expose_to_cpp = function(info, status_codes) {
      defs <- self$tidy_table
      defs[c("n", risk_pin_cols)] <- NULL
      # status carries no entry in `info`: it is coded against the engine
      # codes the individuals hold, not against a data column
      if ("status" %in% colnames(defs)) {
        defs[["status"]] <- unname(
          status_codes[as.character(defs[["status"]])]
        )
      }
      encode_risk_levels_(defs, info,
                          cols = risk_cov_cols_(defs))
    }
  )
)
