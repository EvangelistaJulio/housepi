## housepi - Epidemiological Analysis of Household Outbreaks.
## Copyright 2024-2026 Institut Pasteur (Paris), France.
## SPDX-License-Identifier: GPL-3.0-or-later
##

# `value` pins a stratum during simulation, `const` fixes it in inference too.
risk_pin_cols <- c("value", "const")
risk_reserved_cols <- c("group_risk", "param_name", "n", "rule_group",
                        risk_pin_cols)

# Covariate columns of logic table `tbl`.
risk_cov_cols_ <- function(tbl) strata_cov_cols_(tbl, risk_reserved_cols)

# Assert the covariates of a risk are character and take no reserved name,
# except those in `allow`.
assert_risk_covariates_ <- function(covariates, allow = character()) {
  if (all(is.na(covariates))) return(invisible())
  if (isFALSE(is.character(covariates))) stop_h(
    "Argument {.arg by} must be of type character.",
    format_type_error(covariates)
  )
  assert_valid_covariates(
    covariates, setdiff(c(data_reserved_cols, risk_reserved_cols), allow)
  )
}

# Rule groups of risk logic table `tbl`, keyed on its covariate columns.
risk_rule_groups_ <- function(tbl) {
  compute_rule_groups(tbl, risk_cov_cols_(tbl), "param_name")
}

# Assert risk logic tables resolve every unit of `cells` to one rule per group,
# see assert_strata_fit_().
# * automatic definition `def` -> the hint names its `by` covariates.
assert_risk_fit_ <- function(tables, cells, label, unit, def) {
  hint <- if (identical(def$mode, "auto")) cli::format_inline(
    "Some {.field {def$by}} values are missing: change them, ",
    "or pass a {.arg strata} table with a catch-all rule."
  )
  assert_strata_fit_(lapply(tables, risk_rule_groups_), cells,
                     risk_reserved_cols, "param_name", label, unit, hint)
}

# Count the units of `cells` resolved to each rule of logic table `tbl`,
# see rule_counts_().
count_risk_rules_ <- function(tbl, cells) {
  rule_counts_(risk_rule_groups_(tbl), cells, risk_reserved_cols)
}

# Empty and typed prototype of every reserved column.
risk_reserved_proto <- list(group_risk = character(), param_name = character(),
                            n = integer(), rule_group = integer(),
                            value = double(), const = logical())

# Empty logic table with the reserved columns `cols`, in order.
empty_risk_table <- function(cols) {
  as.data.frame(risk_reserved_proto[cols])
}

# Typed `NA` padding for each of `cols`.
# * covariate column -> character.
risk_na_pad <- function(cols) {
  lapply(stats::setNames(cols, cols), function(col) {
    proto <- risk_reserved_proto[[col]]
    if (is.null(proto)) NA_character_ else proto[NA_integer_]
  })
}

# Parameter definitions of the logic tables of one risk.
# A named parameter takes, in order of precedence:
# * its override, when given;
# * the shared `setup`;
# * the default for an `absolute` risk.
# Its table pins its values.
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

# Whether some rule group of a risk has no reference, making it absolute.
# * automatic definition not yet built from data -> read off its reference.
risk_is_absolute_ <- function(def, tables) {
  if (!length(def)) return(FALSE)
  if (identical(def$mode, "auto") && !length(tables))
    return(is.null(def$reference))
  lacks_reference_(tables)
}

# Reorder `df` to `cols`, padding the missing ones with typed `NA`.
pad_risk_cols_ <- function(df, cols) {
  miss <- setdiff(cols, names(df))
  if (length(miss)) df[miss] <- risk_na_pad(miss)
  df[cols]
}

# Parameter-name prefix: empty, or the given prefix closed by an underscore.
format_risk_prefix_ <- function(prefix) {
  if (identical(prefix, "") || !is.character(prefix)) ""
  else paste0(prefix, "_")
}

# Tidy risk table for the read surface, definitions expanded from data.
# * no data attached -> empty shell with columns `cols`.
helper_show_risk_table <- function(module, cols, is_data_attached) {
  out <- module$tidy_table
  if (NROW(out)) return(out)
  if (length(module$defs) && !is_data_attached) warn_h(
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
    "Invalid risk definition: provide one of {.arg {c('by', 'strata')}}.",
    hint
  )
  "custom"
}

# Validate a custom strata table, with the `params` overriding its strata.
assert_valid_strata_ <- function(strata, params, reference,
                                 allow = character(), symmetric = FALSE) {
  assert_valid_logic_table(strata, allow)
  if (length(reference)) stop_h(
    paste("Argument {.arg reference} is not used",
          "with a custom {.arg strata} table."),
    "Tag one row's {.field param_name} as {.val .reference} instead."
  )
  assert_one_reference_per_group(strata, symmetric)
  if (length(params)) {
    assert_known_params_names(strata, params)
    assert_no_pin_params_clash(strata, params)
  }
  invisible(strata)
}

#' Assert the covariate values of logic table `this_table` are safe.
#' @return `character()` Covariate columns checked.
#' @noRd
assess_covariate_values <- function(this_table) {
  # column names to consider
  colmask <- risk_cov_cols_(this_table)
  for (col in colmask) lapply(this_table[[col]], assert_valid_covariate, col)
  colmask
}

#' Assert a user-provided logic table is valid.
#' @param allow `character()` Reserved names the caller accepts as covariates.
#' @noRd
assert_valid_logic_table <- function(logic_table, allow = character()) {
  assert_strata_shape_(
    logic_table, "param_name", "The strata table",
    "Give each rule its parameter in {.field param_name}."
  )
  # every row names its parameter, or tags the reference
  pnames <- as.character(logic_table$param_name)
  unnamed <- which(is.na(pnames) | !nzchar(pnames))
  if (length(unnamed)) stop_h(
    paste("Strata {cli::qty(length(unnamed))}row{?s} without",
          "{.field param_name}: {unnamed}."),
    "Name the parameter of each row, or tag it {.val .reference}."
  )
  for (pname in setdiff(unique(pnames), ".reference")) {
    assert_valid_definition_name(pname, "parameter")
  }
  # validate covariate column names
  covariate_cols <- risk_cov_cols_(logic_table)
  assert_risk_covariates_(covariate_cols, allow)
  # no covariate column, or only catch-all cells -> duplicates the baseline
  if (all(is.na(logic_table[covariate_cols]))) stop_h(
    paste("Strata table sets no covariate value:",
          "its rules cannot tell members apart."),
    paste("Set a covariate value in a rule, or drop the risk:",
          "the baseline already treats every member alike.")
  )
  assert_valid_pins(logic_table)
}

#' Assert the pin columns of a logic table are valid.
#' @return `data.frame` Input, invisibly.
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
  if (length(conflicts)) stop_h(
    "Conflicting {.field value} for {cli::qty(rows)}strat{?um/a} {conflicts}.",
    "A stratum spread over several rules carries a single {.field value}."
  )
  invisible(tbl)
}

#' Assert each rule group designates at most one reference cell.
#' A reference sets the unit scale of its rule group.
#' * symmetric contact rules -> both orderings of a pair are one cell.
#' @param symmetric `logical(1)` Whether reversed from/to pairs are one cell.
#' @return `data.frame` Input, invisibly.
#' @noRd
assert_one_reference_per_group <- function(strata, symmetric = FALSE) {
  idx_refer <- which(strata$param_name == ".reference")
  if (length(idx_refer) < 2) return(invisible(strata))

  tbl <- risk_rule_groups_(strata)[idx_refer, , drop = FALSE]
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
    \(k) if (length(unique(keys[k])) > 1) idx_refer[k]
  ))
  if (length(clash)) stop_h(
    paste("{cli::qty(length(clash))}Row{?s} {sort(clash)} of one rule group",
          "are tagged {.val .reference}."),
    "A rule group has a single reference, setting its unit scale."
  )
  invisible(strata)
}

#' Read the pins of a logic table, keyed by `param_name`.
#' * stratum spread over several rules -> one pin, read off its first rule.
#' @return `list` of `list(value = , const = )`, one per pinned stratum.
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
#' * fixed stratum -> constant parameter, no prior.
#' @param def `list` Shared MCMC setup of the risk.
#' @param pin `list` | `NULL` Pin of the stratum.
#' @noRd
helper_pin_param_def <- function(def, pin) {
  if (is.null(pin)) return(def)
  if (isTRUE(pin$const)) return(list(value = pin$value, const = TRUE))
  c(def, list(value = pin$value))
}

#' Whether a set of logic tables leaves the scale of its risk free.
#' * rule group without a `.reference` row -> sets the scale on its own,
#'   as every unit matches a rule of every rule group;
#' * table without rules -> carries no scale.
#' @param tables `list` of logic tables with a `param_name` column.
#' @noRd
lacks_reference_ <- function(tables) {
  any(vapply(tables, function(tbl) {
    gid <- risk_rule_groups_(tbl)$rule_group
    !all(gid %in% gid[tbl$param_name == ".reference"])
  }, logical(1)))
}

#' Assert every `params` entry names a parameter the strata table declares.
#' @param params `list` Per-parameter MCMC-setup overrides.
#' @return `list` `params`, invisibly.
#' @noRd
assert_known_params_names <- function(strata, params) {
  unknown <- setdiff(names(params), setdiff(strata$param_name, ".reference"))
  if (length(unknown)) stop_h(
    "Unknown parameter{?s} in {.arg params}: {.val {unknown}}."
  )
  invisible(params)
}

#' Warn that automatic mode drops `params`, as it derives its own parameters.
#' @param params `list` Per-parameter MCMC-setup overrides.
#' @return `list` `params`, invisibly.
#' @noRd
warn_unused_auto_params <- function(params) {
  if (!length(params)) return(invisible(params))
  warn_h(
    "{.arg params} is not used with {.arg by}.",
    paste("Give {.val {names(params)}} a value with {.fn add_parameter},",
          "or name the parameters in a {.arg strata} table.")
  )
  invisible(params)
}

#' Assert no stratum is both pinned and overridden through `params`.
#' @param params `list` Per-parameter MCMC-setup overrides.
#' @return `list` `params`, invisibly.
#' @noRd
assert_no_pin_params_clash <- function(strata, params) {
  clash <- intersect(names(params), names(collect_table_pins(strata)))
  if (length(clash)) stop_h(
    paste("Pinned {cli::qty(clash)}strat{?um/a} also set through",
          "{.arg params}: {.val {clash}}."),
    paste("A pinned stratum carries no MCMC setup:",
          "drop it from {.arg params} or from the strata table.")
  )
  invisible(params)
}

#' Order a risk logic table: `.reference` first, then by decreasing `n`.
#' @noRd
reorder_risks_logic_table <- function(tbl) {
  if (!nrow(tbl)) return(tbl)
  is_ref <- tbl$param_name == ".reference"
  n_val <- suppressWarnings(as.numeric(tbl[["n"]]))
  ord <- order(tbl$rule_group, !is_ref, -n_val, method = "radix")
  tbl[ord, , drop = FALSE]
}

#' Read the data status column as infected-status labels.
#' * escaped or missing member -> `NA`, which matches no status rule
#'   and reads the reference risk.
#' @param x Status column, integer codes or a factor of labels.
#' @param status_codes `integer()` Engine status codes, named after the labels.
#' @noRd
status_as_labels <- function(x, status_codes) {
  labels <- names(status_codes)
  if (is.factor(x) || is.character(x)) {
    x <- as.character(x)
    return(ifelse(x %in% labels, x, NA_character_))
  }
  labels[match(suppressWarnings(as.integer(x)), status_codes)]
}

#' Build the logic table of one group risk from the data.
#' @param this_cov `list` Covariate data column names.
#' @param status_levels `character()` Declared infected-status labels.
#' @param status_observed `logical(1)` Whether the data carries any status.
#' @noRd
group_risk_logic_table <- function(data_wide, this_cov,
                                   status_levels = character(),
                                   status_observed = TRUE) {
  # extract unique covariates and sort by list length
  # status strata are declared by the model, not read off the data
  uni_covariates <- lapply(this_cov, function(f) {
    if (identical(f, "status")) status_levels
    else sort(unique(data_wide[[f]]))
  })
  sorted_indices <- order(lengths(uni_covariates), decreasing = TRUE)
  this_cov <- this_cov[sorted_indices]
  uni_covariates <- uni_covariates[sorted_indices]
  names(uni_covariates) <- this_cov

  # create a dataframe of all possible multi-covariates
  y <- expand.grid(uni_covariates)

  # count every combination in one pass
  # * cell the data does not carry -> matches no stratum, members left out
  carried <- data_wide[stats::complete.cases(data_wide[this_cov]), ,
                       drop = FALSE]
  counted <- table(row_keys_(carried, this_cov))
  multicov_count <- as.integer(counted[row_keys_(y, this_cov)])
  multicov_count[is.na(multicov_count)] <- 0L

  # keep the combinations the data carries, each at every declared status
  # level, as a member may draw any infected status
  keep <- multicov_count > 0
  if ("status" %in% this_cov) {
    others <- setdiff(this_cov, "status")
    keep <- rep(TRUE, nrow(y))
    if (length(others)) {
      held <- data_wide[stats::complete.cases(data_wide[others]), ,
                        drop = FALSE]
      keep <- row_keys_(y, others) %in% row_keys_(held, others)
    }
  }
  y <- y[keep, , drop = FALSE]
  # no status observed -> sample sizes unknown
  y$n <- if ("status" %in% this_cov && isFALSE(status_observed))
    rep(NA_integer_, nrow(y)) else multicov_count[keep]

  y
}

#' Name the parameters of a group risk logic table.
#' @param this_table `data.frame` Conditions of the group risk.
#' @param prefix `character(1)` Prefix of parameter names.
#' @param reference Covariates of the reference group:
#' * `".auto"` -> the group with the largest sample size;
#' * empty -> no reference.
#' @return `data.frame` Input with a `param_name` column.
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

  } else if (length(reference)) {
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

#' Bind the group risk definitions into one `data.frame`.
#' @noRd
tidyup_group_risk_tables <- function(all_tables) {
  df_list <- lapply(names(all_tables), function(risk) {
    this_table <- all_tables[[risk]]
    if (!length(this_table)) return(data.frame(
      group_risk = risk, param_name = ".reference", rule_group = 1L
    ))
    # bind this risk's sub-tables, then compute its rule_groups
    cols <- unique(c("param_name", "n",
                     sort(unlist(lapply(this_table, names)))))
    tbl <- risk_rule_groups_(do.call(rbind, lapply(this_table, pad_risk_cols_,
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

#' Validate the user MCMC setup of a risk, returned as written.
#' It is encoded once the risk's parameters are declared.
#' * every choice missing -> `NULL`.
#' @param dots `list` Distribution and proposal arguments.
#' @param family `character(1)` Holder family.
#' @param name `character(1)` | `NULL` Name of the risk, when it has one.
#' @noRd
unwrap_risk_setup_ <- function(prior, proposal, min, max, dots, family,
                               name = NULL) {
  prior <- housepi_smart_unwrap(prior, "prior")
  proposal <- housepi_smart_unwrap(proposal, "proposal")
  min <- housepi_smart_unwrap(min, "min")
  max <- housepi_smart_unwrap(max, "max")
  assert_valid_law_name(prior, "prior")
  assert_valid_law_name(proposal, "proposal")
  dots <- flatten_definition_dots(dots)
  args <- helper_definition_args(proposal, dots, family, name)$args
  assert_no_sibling_dist_arg(args, family, name)

  setup <- c(list(prior = prior), as.list(args),
             list(min = min, max = max, proposal = proposal))
  is_default <- vapply(setup, \(x) isTRUE(is.na(x)), logical(1))
  # all defaults -> resolved once the data tells absolute from relative
  if (all(is_default)) return(NULL)
  # nolint next: object_usage_linter. read by the cli glue below
  label <- if (is.null(name)) family else name
  if (is_default[["prior"]]) stop_h(
    "Found spare prior arguments for risk {.val {label}}.",
    "Either leave all as defaults or provide a {.arg prior}."
  )
  c(list(prior = prior), dots, list(min = min, max = max, proposal = proposal))
}

#' Return the default MCMC setup of a risk.
#' * absolute -> rates on their own;
#' * relative -> ratios, as in \insertRef{cauchemez_household_2009}{housepi}.
#' @param absolute `logical(1)` Whether the risk carries absolute values.
#' @noRd
default_risk_mcmc_setup <- function(absolute) {
  if (isTRUE(absolute)) list(prior = "uniform", proposal = "rw_log")
  else list(prior = "lognormal", meanlog = 0, sdlog = 3, proposal = "rw_log")
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

########  HousepiGroupRisk  ########

#' Hold the group-risk definitions of a model.
#' @noRd
HousepiGroupRisk <- R6::R6Class(
  classname = "HousepiGroupRisk",
  cloneable = TRUE,

  ####  Private  ####

  private = list(
    deep_clone = function(name, value) helper_clone(value),

    # Add a definition in automatic mode.
    helper_auto_add_ = function(name, by,
                                prefix = NA_character_,
                                reference = NULL,
                                joint = FALSE,
                                .warn_over = TRUE) {
      assert_risk_covariates_(by, allow = "status")
      # remove previous definition if it exists
      rmv_name <- self$remove(name)

      warn_or_store_overwrite_(self, rmv_name, .warn_over)

      # load reference and validate reference names
      if (length(reference))
        reference <- auto_group_reference_(by, reference)

      # add new definition
      self$defs[[name]] <- list(by = by, reference = reference,
                                prefix = prefix, joint = joint,
                                mode = "auto")
      invisible(self)
    },

    # Add a definition in custom mode.
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

    # Assert one risk resolves every member to one rule per group.
    # * risk read on infected members only -> checked at every level in
    #   `status_levels`, whatever the data records.
    assert_rules_fit_ = function(risk_id, data_wide,
                                 status_levels = character()) {
      tables <- self$logic_tables[[risk_id]]
      covs <- setdiff(unique(unlist(lapply(tables, risk_cov_cols_))), "status")
      assert_risk_fit_(
        tables, member_cells_(data_wide, covs, status_levels),
        cli::format_inline("Group risk {.val {risk_id}}"), "member",
        self$defs[[risk_id]]
      )
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
        tbl$n <- count_risk_rules_(tbl, member_cells_(data_wide, cov_cols))
        # no status observed -> every declared stratum stands, sizes unknown
        if ("status" %in% cov_cols && isFALSE(status_observed))
          tbl$n <- NA_integer_
        self$logic_tables[[risk_id]][[ti]] <- tbl
      }
      invisible(self)
    }
  ),

  ####  Active  ####

  active = list(
    # `list` Parameter definitions associated to group risks.
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

    # `character()` Data columns required by group risks.
    data_required = function(value) {
      required <- unique(unlist(lapply(self$defs, `[[`, "by")))
      required <- setdiff(required[!is.na(required)], "status")
      set_source_name_(required, self$family)
    },

    # `character()` Data columns the risks stratify on,
    # declared through `by` or carried by a custom logic table.
    covariate_names = function(value) {
      setdiff(self$raw_covariate_names, "status")
    },

    # `character()` Stratification axes, `status` included.
    # `covariate_names` drops `status`, an always-present model axis,
    # neither a required data column nor a name an observation can collide on.
    raw_covariate_names = function(value) {
      custom <- lapply(unlist(self$logic_tables, recursive = FALSE),
                       risk_cov_cols_)
      covs <- as.character(unlist(c(lapply(self$defs, `[[`, "by"), custom)))
      unique(covs[!is.na(covs)])
    },

    # `character()` Parameters required by group risks.
    param_required = function(value) {
      tables <- unlist(self$logic_tables, recursive = FALSE)
      required <- setdiff(unlist(lapply(tables, `[[`, "param_name")),
                          ".reference")
      set_source_name_(required, self$family)
    },

    # `data.frame()` Combined list with group risk and parameter definitions.
    tidy_table = function(value) {
      tidyup_group_risk_tables(self$logic_tables)
    }
  ),

  ####  Public  ####

  public = list(
    # `character(1)` Holder family.
    family = "group risk",

    # `character()` Cached overwritten names.
    store_warn_over = character(),

    # `list` Current group risk definitions.
    defs = list(),

    # `list` Raw parameter definitions.
    raw_params = list(),

    # `list` Per-parameter MCMC overrides keyed by risk name.
    param_overrides = list(),

    # `list` Parameter table of each group risk.
    logic_tables = list(),

    # Create an empty holder.
    initialize = function() {
      self$defs <- list()
      self$logic_tables <- list()
      self$raw_params <- list()
      self$param_overrides <- list()
      invisible(self)
    },

    # Remove the group risk definition `name`.
    # @return `character(1)` | `NULL` Name removed, `NULL` when none matched.
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

    # Whether some rule group of risk `name` has no reference, making it
    # absolute.
    is_absolute = function(name) {
      risk_is_absolute_(self$defs[[name]], self$logic_tables[[name]])
    },

    # Add a definition, arguments as in add_group_risk().
    add_quo = function(name, by = NA_character_, ...,
                       reference = NULL, joint = FALSE,
                       strata = NULL, params = NULL,
                       prior = NA_character_,
                       min = NA_real_, max = NA_real_,
                       proposal = NA_character_,
                       prefix = NA_character_, .warn_over = TRUE) {
      name <- housepi_smart_unwrap(name, "name")
      by <- housepi_smart_unwrap(by, "by")
      # check the name, `by` being the injected default variable
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

    add_quo_list = function(xquos, append = FALSE) {
      housepi_load_defin_quos(self, "by", xquos, append)
      invisible(self)
    },

    ####  Build  ####

    # Build the logic tables from `data_wide`, then assert every member fits.
    # * keyed on status, none declared -> error at the `wall`, else skipped.
    # @param status_codes `integer()` Engine status codes,
    # named after the labels the risks may stratify on.
    # @param stratified `character()` Risks read on infected members only.
    implement = function(data_wide, status_codes = integer(),
                         wall = FALSE, stratified = character()) {
      # preserve custom-mode entries and drop auto-mode entries
      mask_custom <- names(which(vapply(
        self$defs, \(d) d[["mode"]] == "custom", logical(1)
      )))
      self$logic_tables <- self$logic_tables[mask_custom]

      # read the status through the dictionary:
      # * levels known -> a factor of labels;
      # * before -> integer codes.
      status_levels <- names(status_codes)
      if ("status" %in% colnames(data_wide))
        data_wide[["status"]] <- status_as_labels(data_wide[["status"]],
                                                  status_codes)
      status_observed <- any(!is.na(data_wide[["status"]]))

      for (risk_id in names(self$defs)) {
        # status-keyed risk -> waits for the statuses to be declared,
        # so decorator order stays free, the wall reporting one never declared
        if (!length(status_codes) && private$keys_on_status_(risk_id)) {
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

    # Assert the model reads every group risk defined.
    # * required risk left undefined -> 1 for every member.
    assert = function(required) {
      notused <- template_not_used_find(required, names(self$defs))
      if (length(notused)) stop_h(
        "Unused group risk{?s}: {.val {notused}}.",
        "The model reads no group risk of that name.",
        find_closest_match(notused, required)
      )
      invisible(self)
    },

    # Encode the definitions for the engine.
    # @param info `list` of `list`s Class and levels of the data columns.
    # @param status_codes `integer()` Engine status codes, named by label.
    expose_to_cpp = function(info, status_codes) {
      defs <- self$tidy_table
      defs[c("n", risk_pin_cols)] <- NULL
      # status has no entry in `info`, coded against the engine codes
      # members hold, not against a data column
      encode_strata_levels_(code_status_col_(defs, status_codes), info,
                            cols = risk_cov_cols_(defs))
    }
  )
)
