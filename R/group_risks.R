## housepi - Epidemiological Analysis of Household Outbreaks.
##
## Copyright © 2024-2026 Institut Pasteur (Paris), CNRS and INSERM.
## Distributed under the terms of the GNU General Public License (GPLv3).
## See the LICENSE file for details.
##

risk_stop_hint <- paste("Choose {.arg by} for automatic definitions,",
                        "{.arg strata} for custom.")

# Reserved column names (non-covariate) in a logic table.
data_reserved_cols <- c("id", "hhid", "status")
# `value` pins a stratum during simulation, `const` fixes it in inference too.
risk_pin_cols <- c("value", "const")
risk_reserved_cols <- c("group_risk", "param_name", "n", "rule_group",
                        risk_pin_cols)

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

# Parameter definitions of one logic table. A named parameter takes its
# override when given, else the shared `fallback`; the table pins its values.
risk_table_params_ <- function(tbl, overrides, fallback) {
  pins <- collect_table_pins(tbl)
  out <- list()
  for (par_id in setdiff(tbl$param_name, ".reference")) {
    def <- if (!is.null(overrides[[par_id]])) overrides[[par_id]] else fallback
    out[[par_id]] <- helper_pin_param_def(def, pins[[par_id]])
  }
  out
}

# Code rule columns as integers over the factor levels held in `info`.
# `key` maps a column to its `info` entry; columns without one stay as they
# are. An NA rule cell is the wildcard: it codes to NA, so the engine skips
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

# Throw error if empty definition.
stop_empty_definition <- function() {
  stop_h(
    "Invalid risk definition: provide one between {.arg {c('by', 'strata')}}.",
    risk_stop_hint
  )
}

# Throw error for ambiguos risk definition.
stop_ambiguos_definition <- function() {
  stop_h(
    "Impossible to provide at the same time {.arg {c('by', 'logic_table')}}.",
    risk_stop_hint
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
  colmask <- setdiff(names(this_table), risk_reserved_cols)
  for (col in colmask) {
    sapply(this_table[[col]], function(x) assert_valid_covariate(x, col))
  }
  colmask
}

#' Drop the rules naming a value the data never carries.
#' A cell encodes against the observed levels, so a value absent from them
#' would read as `NA`, the wildcard: the rule would widen instead of dropping.
#' @param tbl `data.frame` logic table with a `param_name` column.
#' @param known `list` Values the data carries, named after its columns.
#' @param cols `character()` Data column every column of `tbl` reads, named
#' after the columns of `tbl`. Defaults to the columns of `known` themselves.
#' @param label `character(1)` The risk, as it appears in the warning.
#' @return `data.frame` the input without the dropped rules.
#' @keywords internal
#' @noRd
drop_unknown_strata_rows <- function(tbl, known, label, cols = NULL) {
  if (is.null(cols)) cols <- stats::setNames(nm = names(known))
  keep <- rep(TRUE, nrow(tbl))
  # the two sides of a contact pair read one data column: report it once
  missing_vals <- list()
  for (col in intersect(names(cols), colnames(tbl))) {
    field <- cols[[col]]
    # a wildcard cell keys on nothing: it names no value to be unknown
    vals <- as.character(tbl[[col]])
    unknown <- !is.na(vals) & !vals %in% as.character(known[[field]])
    if (!any(unknown)) next
    missing_vals[[field]] <- union(missing_vals[[field]],
                                   unique(vals[unknown]))
    keep <- keep & !unknown
  }
  for (field in names(missing_vals)) warn_h(
    paste("Unknown {cli::qty(missing_vals[[field]])}definition{?s}",
          "for {label} in data {.field {field}}:",
          "{.val {missing_vals[[field]]}}."),
    `!` = "Associated entries will be skipped."
  )
  tbl[keep, , drop = FALSE]
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
  # validate covariate column names
  covariate_cols <- setdiff(colnames(logic_table), risk_reserved_cols)
  assert_valid_covariates(covariate_cols, allow)
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
  # value
  value <- tbl[["value"]]
  if (is.null(value)) value <- rep(NA_real_, nrow(tbl))
  const <- tbl[["const"]]
  if (is.null(const)) const <- rep(NA, nrow(tbl))

  bad <- tbl$param_name == ".reference" & (!is.na(value) | !is.na(const))
  if (any(bad)) stop_h(
    paste("The {.val .reference} is pinned at",
          "{cli::qty(sum(bad))}row{?s} {which(bad)}."),
    "The reference stratum is 1 by definition: it carries no pin."
  )
  bad <- const %in% TRUE & is.na(value)
  if (any(bad)) stop_h(
    paste("No {.field value} for fixed",
          "strat{?um/a} {.val {unique(tbl$param_name[bad])}}."),
    "A {.code const = TRUE} stratum needs the value it is fixed at."
  )
  # a stratum spread over several rules carries a single pin, whichever of
  # its rules the data keeps
  pins <- data.frame(value, const)
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
  cov_cols <- setdiff(colnames(tbl), risk_reserved_cols)

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
    } else if (length(candidates) == 0) {
      stop_h("Rule grouping failed: no maximal pattern contains row's columns.")
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
  # a wildcard rule cell constrains nothing; a missing data cell matches
  # nothing but a wildcard
  hit <- Reduce(`&`, lapply(names(cells), function(col) {
    outer(as.character(cells[[col]]), as.character(tbl[[col]]),
          function(d, r) is.na(r) | (!is.na(d) & d == r))
  }))
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
  cov_cols <- setdiff(colnames(tbl), risk_reserved_cols)
  covariates <- sub("^from_", "", grep("^from_", cov_cols, value = TRUE))
  pair_wise <- isTRUE(symmetric) && length(covariates) > 0

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
  cols <- intersect(risk_pin_cols, colnames(tbl))
  if (length(cols) == 0) return(list())
  pins <- list()
  for (par_id in setdiff(unique(tbl$param_name), ".reference")) {
    pin <- list(value = NA_real_, const = NA)
    first <- match(par_id, tbl$param_name)
    pin[cols] <- as.list(tbl[first, cols, drop = FALSE])
    if (!is.na(pin$value) || !is.na(pin$const)) pins[[par_id]] <- pin
  }
  pins
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
#' A table the data left empty carries no rule, hence no scale.
#' @param tables `list` of logic tables with a `param_name` column.
#' @return `logical(1)`
#' @keywords internal
#' @noRd
lacks_reference_ <- function(tables) {
  length(tables) > 0 &&
    !any(vapply(tables, function(tbl) ".reference" %in% tbl$param_name,
                logical(1)))
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

  # no status is observed: every declared stratum stands, sample sizes unknown
  if ("status" %in% this_cov && isFALSE(status_observed)) {
    y$n <- NA_integer_
    return(y)
  }

  # create the final truth table and reverse it
  y <- y[multicov_count > 0, , drop = FALSE]
  # add counts
  y$n <- multicov_count[multicov_count > 0]

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
    # find definition row
    idx_refer <- unname(which(
      apply(this_table[colmask], 1,
            function(row) all(row == unlist(reference[colmask])))
    ))
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

#' Combine list with group risk definition into a single dataframe
#' @param all_tables `list`. Definitions for each group risk.
#' @return `data.frame`
#' @keywords internal
#' @noRd
tidyup_group_risk_tables <- function(all_tables) {
  # tidy version of group risk tables
  cbind_risk <- function(x, risk) cbind(group_risk = risk, x)
  df_list <- c()
  for (risk in names(all_tables)) {
    this_table <- all_tables[[risk]]
    if (length(this_table) == 0) {
      df <- list(data.frame(group_risk = risk, param_name = ".reference",
                            rule_group = 1L))
    } else {
      # concat this risk's sub-tables, then compute rule_groups within risk
      sub_cols <- unique(c("param_name", "n",
                           sort(unlist(lapply(this_table, names)))))
      sub_combined <- do.call(rbind, lapply(this_table, pad_risk_cols_,
                                            cols = sub_cols))
      sub_combined <- compute_rule_groups(sub_combined)
      sub_combined <- reorder_risks_logic_table(sub_combined)
      df <- list(cbind_risk(sub_combined, risk))
    }
    df_list <- c(df_list, df)
  }

  # all unique column names reordered, pin columns only when a pin exists
  present <- unique(unlist(lapply(df_list, names)))
  all_cols <- unique(c(setdiff(risk_reserved_cols, risk_pin_cols),
                       intersect(risk_pin_cols, present), sort(present)))

  # Add missing columns with NA
  combined_df <- do.call(rbind, lapply(df_list, pad_risk_cols_,
                                       cols = all_cols))
  rownames(combined_df) <- NULL

  as.data.frame(combined_df)
}

#' Return list with default MCMC setup if user choice is missing.
#' Absolute risks are rates on their own; relative ones are ratios,
#' whose defaults follow \insertRef{cauchemez_household_2009}{housepi}.
#' @return `list`
#' @keywords internal
#' @noRd
set_group_risks_mcmc_setup <- function(name, prior, args, min, max, proposal,
                                       absolute = FALSE) {
  # check wheter default configuration
  is_default <- sapply(c(prior, args, min, max, proposal),
                       function(x) is.na(x) || is.null(x))
  # if all defaults, provide default choice
  if (all(is_default)) {
    if (isTRUE(absolute)) list(prior = "uniform", proposal = "rw_log")
    else list(prior = "lognormal", meanlog = 0, sdlog = 3,
              min = 0, max = 10, proposal = "rw_log")
  } else {
    # if prior is default -> spare arguments
    if (is_default[1]) stop_h(
      "Found spare prior arguments for risk {.val {name}}.",
      "Either leave all as defaults or provide a {.arg prior}."
    )
    # else go with user choice
    c(list(prior = prior), args,
      list(min = min, max = max, proposal = proposal))
  }
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
      assert_valid_logic_table(strata, allow = "status")
      if (length(reference) > 0) stop_h(
        paste("Argument {.arg reference} is not used",
              "with a custom {.arg strata} table."),
        "Tag one row's {.field param_name} as {.val .reference} instead."
      )
      assert_one_reference_per_group(strata)
      # validate params keys against logic table param_names
      if (length(params) > 0) {
        assert_known_params_names(strata, params)
        assert_no_pin_params_clash(strata, params)
      }

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
      by <- setdiff(unique(unlist(lapply(self$logic_tables[[name]], colnames))),
                    risk_reserved_cols)
      self$defs[[name]] <- list(by = by, prefix = NA_character_,
                                reference = NULL, joint = FALSE,
                                mode = "custom")
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

      if (all(is.na(covariates))) {
        self$logic_tables[[risk_id]] <- list()
        return(invisible(self))
      }

      covariates <- sort(as.character(covariates[!is.na(covariates)]))
      # status is carried by the model, not required as a data column
      assert_data_covariates(setdiff(covariates, "status"), colnames(data_wide))

      # cast to list if joint
      if (isTRUE(this_def[["joint"]])) covariates <- list(covariates)

      # default prefix
      prefix <- as.character(this_def[["prefix"]])
      if (is.na(prefix)) prefix <- substring(tolower(risk_id), 1, 3)

      new_logic_table <- list()
      for (this_cov in covariates) {
        # 1. create logic table with group risks
        this_table <- group_risk_logic_table(data_wide, this_cov,
                                             status_levels, status_observed)
        # associate names for parameters
        this_table <- associate_group_param_names(
          this_table, prefix, this_def$reference
        )
        new_logic_table[[length(new_logic_table) + 1]] <- this_table

        # store the new risk logic table
        self$logic_tables[[risk_id]] <- new_logic_table
      }
    },

    # Validate custom mode logic table against data, then add counts.
    helper_implement_custom_ = function(risk_id, data_wide,
                                        status_levels, status_observed) {
      for (ti in seq_along(self$logic_tables[[risk_id]])) {
        tbl <- self$logic_tables[[risk_id]][[ti]]
        cov_cols <- setdiff(colnames(tbl), risk_reserved_cols)
        # status is carried by the model, not required as a data column
        assert_data_covariates(setdiff(cov_cols, "status"), colnames(data_wide))
        # validate that covariate values exist in data
        tbl <- drop_unknown_strata_rows(
          tbl,
          known = lapply(stats::setNames(nm = cov_cols), function(col) {
            # status levels are declared by the model, not read off the data
            if (identical(col, "status")) status_levels
            else unique(data_wide[[col]])
          }),
          label = cli::format_inline("risk {.val {risk_id}}")
        )
        # add counts: an emptied table keeps the column, typed
        tbl$n <- if (nrow(tbl) == 0) integer(0) else
          count_resolved_rules(tbl, data_wide[cov_cols])
        # no status is observed: every declared stratum stands, sizes unknown
        if ("status" %in% cov_cols && isFALSE(status_observed)) {
          tbl$n <- NA_integer_
        } else {
          # skip entries with no matching data
          tbl <- tbl[tbl$n > 0, , drop = FALSE]
        }
        self$logic_tables[[risk_id]][[ti]] <- tbl
      }
      # a table whose every stratum was dropped stratifies nothing:
      # it leaves no definition behind, so every member reads the reference
      self$logic_tables[[risk_id]] <- Filter(
        function(tbl) nrow(tbl) > 0, self$logic_tables[[risk_id]]
      )
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
        for (this_table in self$logic_tables[[risk_id]]) {
          par <- risk_table_params_(this_table,
                                    self$param_overrides[[risk_id]],
                                    self$raw_params[[risk_id]])
          new_params[names(par)] <- par
        }
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
                       function(tbl) setdiff(colnames(tbl), risk_reserved_cols))
      covs <- as.character(unlist(c(lapply(self$defs, `[[`, "by"), custom)))
      unique(covs[!is.na(covs)])
    },

    #' @field param_required `character()`. Parameters required by group risks.
    param_required = function(value) {
      required <- c()
      for (risk_id in names(self$logic_tables)) {
        for (this_table in self$logic_tables[[risk_id]]) {
          required <- c(required, this_table$param_name)
        }
      }
      required <- setdiff(unique(required), c(".reference"))
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
        # check if old definition was just placeholder
        placeholder <- FALSE
        if (all(is.na(self$defs[[name]]$by))) placeholder <- TRUE

        # drop previous definition
        self$defs[name] <- NULL
        self$logic_tables[name] <- NULL
        self$raw_params[name] <- NULL
        self$param_overrides[name] <- NULL
        if (!placeholder) return(name)

      }

      c()
    },

    #' @description
    #' Whether a risk carries absolute values, having no reference.
    #' @param name `character(1)`. Name of the group risk definition.
    #' @return `logical(1)`
    is_absolute = function(name) {
      def <- self$defs[[name]]
      if (is.null(def)) return(FALSE)
      if (identical(def$mode, "auto")) return(is.null(def$reference))
      lacks_reference_(self$logic_tables[[name]])
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

      # raw parameter definition
      prior <- housepi_smart_unwrap(prior, "prior")
      proposal <- housepi_smart_unwrap(proposal, "proposal")
      min <- housepi_smart_unwrap(min, "min")
      max <- housepi_smart_unwrap(max, "max")
      # extract args and proposal from the dots
      def_args <- helper_definition_args(proposal, list(...))
      args <- def_args$args
      proposal <- def_args$proposal
      assert_no_sibling_dist_arg(args, self$family, name)

      # distinguish between automatic and custom definitions
      is_auto <- !all(is.na(by))
      is_custom <- !is.null(strata)

      if (is_auto) {
        if (is_custom) stop_ambiguos_definition()
        warn_unused_auto_params(params)
        private$helper_auto_add_(
          name = name, by = by, reference = reference, prefix = prefix,
          joint = joint, .warn_over = .warn_over
        )
      } else if (is_custom) {
        private$helper_custom_add_(
          name = name, strata = strata, params = params,
          reference = reference, .warn_over = .warn_over
        )
      } else {
        stop_empty_definition()
      }

      # store shared MCMC setup (overwriting any existing)
      self$raw_params[[name]] <- set_group_risks_mcmc_setup(
        name = name, prior = prior, args = args,
        min = min, max = max, proposal = proposal,
        absolute = self$is_absolute(name)
      )

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
    implement = function(data_wide, status_codes = integer(),
                         wall = FALSE) {
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
          private$helper_implement_custom_(risk_id, data_wide,
                                           status_levels, status_observed)
        }
      }

      invisible(self)
    },

    #' @description
    #' Assert that group risk definitions match the requirements.
    #' @param required `character()`. Required group risks.
    #' @param wall `logical(1)`. Whether to raise error.
    assert = function(required, wall) {
      # load definitions
      defined <- names(self$defs)

      # UPDATE
      # point out group risks that are not used
      notused <- template_not_used_find(required, defined)
      # drop unused group risks
      if (length(notused) > 0) {
        sms <- c("Dropped unused group risk{?s}: {.val {notused}}.",
                 find_closest_match(notused, required))
        if (isTRUE(wall)) {
          for (which in notused) self$remove(which)
          stop_h(sms)
        }
        warn_h(sms)
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
      encode_risk_levels_(defs, info)
    }
  )
)
