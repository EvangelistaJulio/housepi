## housepi - Epidemiological Analysis of Household Outbreaks.
## Copyright 2024-2026 Institut Pasteur (Paris), France.
## SPDX-License-Identifier: GPL-3.0-or-later
##

# `const` pins a stratum, in inference and simulation.
rate_pin_cols <- "const"
rate_reserved_cols <- c("group_rate", "param_name", "n", "rule_group",
                        rate_pin_cols)

# Covariate columns of logic table `tbl`.
rate_cov_cols_ <- function(tbl) strata_cov_cols_(tbl, rate_reserved_cols)

# Assert the covariates of a rate are character and take no reserved name,
# except those in `allow`.
assert_rate_covariates_ <- function(covariates, allow = character()) {
  if (all(is.na(covariates))) return(invisible())
  if (isFALSE(is.character(covariates))) stop_h(
    "Argument {.arg by} must be of type character.",
    format_type_error(covariates)
  )
  assert_valid_covariates(
    covariates, setdiff(c(data_reserved_cols, rate_reserved_cols), allow)
  )
}

# Rule groups of rate logic table `tbl`, keyed on its covariate columns.
rate_rule_groups_ <- function(tbl) {
  compute_rule_groups(tbl, rate_cov_cols_(tbl), "param_name")
}

# Assert rate logic tables resolve every unit of `cells` to one rule per group,
# see assert_strata_fit_().
# * automatic definition `def` -> the hint names its `by` covariates.
assert_rate_fit_ <- function(tables, cells, label, unit, def) {
  hint <- if (identical(def$mode, "auto")) cli::format_inline(
    "Some {.field {def$by}} values are missing: change them, ",
    "or pass a {.arg strata} table with a catch-all rule."
  )
  assert_strata_fit_(lapply(tables, rate_rule_groups_), cells,
                     rate_reserved_cols, "param_name", label, unit, hint)
}

# Count the units of `cells` resolved to each rule of logic table `tbl`,
# see rule_counts_().
count_rate_rules_ <- function(tbl, cells) {
  rule_counts_(rate_rule_groups_(tbl), cells, rate_reserved_cols)
}

# Empty and typed prototype of every reserved column.
rate_reserved_proto <- list(group_rate = character(), param_name = character(),
                            n = integer(), rule_group = integer(),
                            const = double())

# Empty logic table with the reserved columns `cols`, in order.
empty_rate_table <- function(cols) {
  as.data.frame(rate_reserved_proto[cols])
}

# Typed `NA` padding for each of `cols`.
# * covariate column -> character.
rate_na_pad <- function(cols) {
  lapply(stats::setNames(cols, cols), function(col) {
    proto <- rate_reserved_proto[[col]]
    if (is.null(proto)) NA_character_ else proto[NA_integer_]
  })
}

# Parameter definitions of the logic tables of one rate.
# A named parameter takes, in order of precedence:
# * its override, when given;
# * the shared `config`;
# * the default for an `absolute` rate.
# Its table pins its values.
rate_params_ <- function(tables, overrides, config, absolute) {
  if (is.null(config)) config <- default_rate_mcmc_config(absolute)
  out <- list()
  for (tbl in tables) {
    pins <- collect_table_pins(tbl)
    for (par_id in setdiff(tbl$param_name, ".reference")) {
      def <- if (!is.null(overrides[[par_id]])) overrides[[par_id]] else config
      out[[par_id]] <- helper_pin_param_def(def, pins[[par_id]])
    }
  }
  out
}

# Whether some rule group of a rate has no reference, making it absolute.
# * automatic definition not yet built from data -> read off its reference.
rate_is_absolute_ <- function(def, tables) {
  if (!length(def)) return(FALSE)
  if (identical(def$mode, "auto") && !length(tables))
    return(is.null(def$reference))
  lacks_reference_(tables)
}

# Reorder `df` to `cols`, padding the missing ones with typed `NA`.
pad_rate_cols_ <- function(df, cols) {
  miss <- setdiff(cols, names(df))
  if (length(miss)) df[miss] <- rate_na_pad(miss)
  df[cols]
}

# Parameter-name prefix: empty, or the given prefix closed by an underscore.
format_rate_prefix_ <- function(prefix) {
  if (identical(prefix, "")) "" else paste0(prefix, "_")
}

# Tidy rate table for the read surface, definitions expanded from data.
# * no data attached -> empty shell with columns `cols`.
helper_show_rate_table <- function(module, cols, is_data_attached) {
  out <- module$tidy_table
  if (NROW(out)) return(out)
  if (length(module$defs) && !is_data_attached) warn_h(
    "Cannot expand {module$family} definitions without data.",
    "Attach household data using {.fn add_dataset}."
  )
  empty_rate_table(cols)
}

# Mode of a rate definition: automatic through `by`, custom through `strata`.
rate_mode_ <- function(by, strata, params) {
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
    "Invalid rate definition: provide one of {.arg {c('by', 'strata')}}.",
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
  colmask <- rate_cov_cols_(this_table)
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
  covariate_cols <- rate_cov_cols_(logic_table)
  assert_rate_covariates_(covariate_cols, allow)
  # no covariate column, or only catch-all cells -> duplicates the baseline
  if (all(is.na(logic_table[covariate_cols]))) stop_h(
    paste("Strata table sets no covariate value:",
          "its rules cannot tell members apart."),
    paste("Set a covariate value in a rule, or drop the rate:",
          "the baseline already treats every member alike.")
  )
  assert_valid_pins(logic_table)
}

#' Assert the pin columns of a logic table are valid.
#' @return `data.frame` Input, invisibly.
#' @noRd
assert_valid_pins <- function(tbl) {
  if (!any(rate_pin_cols %in% colnames(tbl))) return(invisible(tbl))
  pins <- tbl[["const"]]
  # an all-NA column pins nothing, whichever type it was written with
  if (!all(is.na(pins)) && !is.numeric(pins)) stop_h(
    "Column {.field const} of the strata table must be numeric.",
    format_type_error(pins)
  )
  bad <- tbl$param_name == ".reference" & !is.na(pins)
  if (any(bad)) stop_h(
    paste("The {.val .reference} is pinned at",
          "{cli::qty(sum(bad))}row{?s} {which(bad)}."),
    "The reference stratum is 1 by definition: it carries no pin."
  )
  # a stratum spread over several rules carries a single pin, whichever of
  # its rules the data keeps
  rows <- Filter(function(i) length(unique(pins[i])) > 1,
                 split(seq_len(nrow(tbl)), tbl$param_name))
  conflicts <- vapply(names(rows), function(nm) {
    cli::format_inline("{.val {nm}} (rows {rows[[nm]]})")
  }, character(1))
  if (length(conflicts)) stop_h(
    "Conflicting {.field const} for {cli::qty(rows)}strat{?um/a} {conflicts}.",
    "A stratum spread over several rules carries a single {.field const}."
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

  tbl <- rate_rule_groups_(strata)[idx_refer, , drop = FALSE]
  cov_cols <- rate_cov_cols_(tbl)
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
#' @return `list` of `numeric(1)`, one per pinned stratum.
#' @noRd
collect_table_pins <- function(tbl) {
  if (!any(rate_pin_cols %in% colnames(tbl))) return(list())
  out <- list()
  for (par_id in setdiff(unique(tbl$param_name), ".reference")) {
    pin <- tbl[["const"]][match(par_id, tbl$param_name)]
    if (!is.na(pin)) out[[par_id]] <- pin
  }
  out
}

#' Fold a stratum pin into its parameter definition.
#' * pinned stratum -> constant parameter, no prior.
#' @param def `list` Shared MCMC configuration of the rate.
#' @param pin `numeric(1)` | `NULL` Pin of the stratum.
#' @noRd
helper_pin_param_def <- function(def, pin) {
  if (is.null(pin)) def else list(const = pin)
}

#' Whether a set of logic tables leaves the scale of its rate free.
#' * rule group without a `.reference` row -> sets the scale on its own,
#'   as every unit matches a rule of every rule group;
#' * table without rules -> carries no scale.
#' @param tables `list` of logic tables with a `param_name` column.
#' @noRd
lacks_reference_ <- function(tables) {
  any(vapply(tables, function(tbl) {
    gid <- rate_rule_groups_(tbl)$rule_group
    !all(gid %in% gid[tbl$param_name == ".reference"])
  }, logical(1)))
}

#' Assert every `params` entry names a parameter the strata table declares.
#' @param params `list` Per-parameter MCMC-configuration overrides.
#' @return `list` `params`, invisibly.
#' @noRd
assert_known_params_names <- function(strata, params) {
  declared <- setdiff(strata$param_name, ".reference")
  unknown <- setdiff(names(params), declared)
  if (length(unknown)) stop_h(
    "Unknown parameter{?s} in {.arg params}: {.val {unknown}}.",
    find_closest_match(
      unknown, declared,
      otherwise = "Name a {.field param_name} of the {.arg strata} table."
    )
  )
  invisible(params)
}

#' Warn that automatic mode drops `params`, as it derives its own parameters.
#' @param params `list` Per-parameter MCMC-configuration overrides.
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
#' @param params `list` Per-parameter MCMC-configuration overrides.
#' @return `list` `params`, invisibly.
#' @noRd
assert_no_pin_params_clash <- function(strata, params) {
  clash <- intersect(names(params), names(collect_table_pins(strata)))
  if (length(clash)) stop_h(
    paste("Pinned {cli::qty(clash)}strat{?um/a} also set through",
          "{.arg params}: {.val {clash}}."),
    paste("A pinned stratum carries no MCMC configuration:",
          "drop it from {.arg params} or from the strata table.")
  )
  invisible(params)
}

#' Order a rate logic table: `.reference` first, then by decreasing `n`.
#' @noRd
reorder_rates_logic_table <- function(tbl) {
  if (!nrow(tbl)) return(tbl)
  is_ref <- tbl$param_name == ".reference"
  n_val <- suppressWarnings(as.numeric(tbl[["n"]]))
  ord <- order(tbl$rule_group, !is_ref, -n_val, method = "radix")
  tbl[ord, , drop = FALSE]
}

#' Read the data status column as infected-status labels.
#' * escaped or missing member -> `NA`, which matches no status rule
#'   and reads the reference rate.
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

#' Build the logic table of one group rate from the data.
#' @param this_cov `list` Covariate data column names.
#' @param status_levels `character()` Declared infected-status labels.
#' @param status_observed `logical(1)` Whether the data carries any status.
#' @noRd
group_rate_logic_table <- function(dataset, this_cov,
                                   status_levels = character(),
                                   status_observed = TRUE) {
  # extract unique covariates and sort by list length
  # status strata are declared by the model, not read off the data
  uni_covariates <- lapply(this_cov, function(f) {
    if (identical(f, "status")) status_levels
    else sort(unique(dataset[[f]]))
  })
  sorted_indices <- order(lengths(uni_covariates), decreasing = TRUE)
  this_cov <- this_cov[sorted_indices]
  uni_covariates <- uni_covariates[sorted_indices]
  names(uni_covariates) <- this_cov

  # create a dataframe of all possible multi-covariates
  y <- expand.grid(uni_covariates)

  # count every combination in one pass
  # * cell the data does not carry -> matches no stratum, members left out
  carried <- dataset[stats::complete.cases(dataset[this_cov]), ,
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
      held <- dataset[stats::complete.cases(dataset[others]), ,
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

#' Name the parameters of a group rate logic table.
#' @param this_table `data.frame` Conditions of the group rate.
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

  prefix <- format_rate_prefix_(prefix)

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

#' Bind the group rate definitions into one `data.frame`.
#' @noRd
tidyup_group_rate_tables <- function(all_tables) {
  df_list <- lapply(names(all_tables), function(rate) {
    this_table <- all_tables[[rate]]
    if (!length(this_table)) return(data.frame(
      group_rate = rate, param_name = ".reference", rule_group = 1L
    ))
    # bind this rate's sub-tables, then compute its rule_groups
    cols <- unique(c("param_name", "n",
                     sort(unlist(lapply(this_table, names)))))
    tbl <- rate_rule_groups_(do.call(rbind, lapply(this_table, pad_rate_cols_,
                                                   cols = cols)))
    cbind(group_rate = rate, reorder_rates_logic_table(tbl))
  })

  # all unique column names reordered, pin columns only when a pin exists
  present <- unique(unlist(lapply(df_list, names)))
  all_cols <- unique(c(setdiff(rate_reserved_cols, rate_pin_cols),
                       intersect(rate_pin_cols, present), sort(present)))

  combined_df <- do.call(rbind, lapply(df_list, pad_rate_cols_,
                                       cols = all_cols))
  rownames(combined_df) <- NULL
  as.data.frame(combined_df)
}

#' Validate the user MCMC configuration of a rate, returned as written.
#' It is encoded once the rate's parameters are declared.
#' * every choice missing -> `NULL`.
#' @param dots `list` Distribution and proposal arguments.
#' @param family `character(1)` Holder family.
#' @param name `character(1)` | `NULL` Name of the rate, when it has one.
#' @noRd
unwrap_rate_config_ <- function(prior, proposal, min, max, dots, family,
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

  config <- c(list(prior = prior), as.list(args),
              list(min = min, max = max, proposal = proposal))
  is_default <- vapply(config, \(x) isTRUE(is.na(x)), logical(1))
  # all defaults -> resolved once the data tells absolute from relative
  if (all(is_default)) return(NULL)
  # nolint next: object_usage_linter. read by the cli glue below
  label <- if (is.null(name)) family else name
  if (is_default[["prior"]]) stop_h(
    "Found spare prior arguments for rate {.val {label}}.",
    "Either leave all as defaults or provide a {.arg prior}."
  )
  c(list(prior = prior), dots, list(min = min, max = max, proposal = proposal))
}

#' Return the default MCMC configuration of a rate.
#' * absolute -> rates on their own;
#' * relative -> ratios, as in \insertRef{cauchemez_household_2009}{housepi}.
#' @param absolute `logical(1)` Whether the rate carries absolute values.
#' @noRd
default_rate_mcmc_config <- function(absolute) {
  if (isTRUE(absolute)) list(prior = "uniform")
  else list(prior = "lognormal", meanlog = 0, sdlog = 3)
}

# Validate and shape the automatic reference of group rates.
auto_group_reference_ <- function(by, reference) {
  if (identical(reference, ".auto")) return(".auto")

  if (length(by) != length(reference)) stop_h(
    "Arguments {.arg {c('by', 'reference')}} should have equal dimensions."
  )
  unknown <- setdiff(names(reference), by)
  if (length(unknown)) stop_h(
    "Unknown covariate{?s} in {.arg reference}: {.val {unknown}}.",
    find_closest_match(unknown, by,
                       otherwise = "Name each entry after a {.arg by} column.")
  )
  reference <- as.list(reference)
  stats::setNames(reference, by)
}

########  HousepiGroupRate  ########

#' Hold the group-rate definitions of a model.
#' @noRd
HousepiGroupRate <- R6::R6Class(
  classname = "HousepiGroupRate",
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
      assert_rate_covariates_(by, allow = "status")
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
      self$defs[[name]] <- list(by = rate_cov_cols_(strata),
                                prefix = NA_character_, reference = NULL,
                                joint = FALSE, mode = "custom")
      invisible(self)
    },

    # Assert one rate resolves every member to one rule per group.
    # * rate read on infected members only -> checked at every level in
    #   `status_levels`, whatever the data records.
    assert_rules_fit_ = function(rate_id, dataset,
                                 status_levels = character()) {
      tables <- self$logic_tables[[rate_id]]
      covs <- setdiff(unique(unlist(lapply(tables, rate_cov_cols_))), "status")
      assert_rate_fit_(
        tables, member_cells_(dataset, covs, status_levels),
        cli::format_inline("Group rate {.val {rate_id}}"), "member",
        self$defs[[rate_id]]
      )
      invisible(self)
    },

    # Whether one rate stratifies on the infected status.
    keys_on_status_ = function(rate_id) {
      "status" %in% self$defs[[rate_id]][["by"]]
    },

    # Implement through automatic mode.
    helper_implement_auto_ = function(rate_id, dataset,
                                      status_levels, status_observed) {
      # definition of the rate
      this_def <- self$defs[[rate_id]]
      # data covariates
      covariates <- as.list(this_def[["by"]])
      covariates <- sort(as.character(covariates[!is.na(covariates)]))
      # status is carried by the model, not required as a data column
      assert_data_covariates(setdiff(covariates, "status"), colnames(dataset))

      # cast to list if joint
      if (isTRUE(this_def[["joint"]])) covariates <- list(covariates)

      # default prefix
      prefix <- as.character(this_def[["prefix"]])
      if (is.na(prefix)) prefix <- tolower(rate_id)

      self$logic_tables[[rate_id]] <- lapply(covariates, function(this_cov) {
        this_table <- group_rate_logic_table(dataset, this_cov,
                                             status_levels, status_observed)
        associate_group_param_names(this_table, prefix, this_def$reference)
      })
    },

    # Validate custom mode logic table against data, then add counts.
    helper_implement_custom_ = function(rate_id, dataset, status_observed) {
      for (ti in seq_along(self$logic_tables[[rate_id]])) {
        tbl <- self$logic_tables[[rate_id]][[ti]]
        cov_cols <- rate_cov_cols_(tbl)
        # status is carried by the model, not required as a data column
        assert_data_covariates(setdiff(cov_cols, "status"), colnames(dataset))
        tbl$n <- count_rate_rules_(tbl, member_cells_(dataset, cov_cols))
        # no status observed -> every declared stratum stands, sizes unknown
        if ("status" %in% cov_cols && isFALSE(status_observed))
          tbl$n <- NA_integer_
        self$logic_tables[[rate_id]][[ti]] <- tbl
      }
      invisible(self)
    }
  ),

  ####  Active  ####

  active = list(
    # `list` Parameter definitions associated to group rates.
    parameters = function(value) {
      new_params <- list()
      for (rate_id in names(self$logic_tables)) {
        par <- rate_params_(self$logic_tables[[rate_id]],
                            self$param_overrides[[rate_id]],
                            self$raw_params[[rate_id]],
                            self$is_absolute(rate_id))
        new_params[names(par)] <- par
      }
      new_params
    },

    # `character()` Data columns required by group rates.
    data_required = function(value) {
      required <- unique(unlist(lapply(self$defs, `[[`, "by")))
      required <- setdiff(required[!is.na(required)], "status")
      set_source_name_(required, self$family)
    },

    # `character()` Data columns the rates stratify on,
    # declared through `by` or carried by a custom logic table.
    covariate_names = function(value) {
      setdiff(self$raw_covariate_names, "status")
    },

    # `character()` Stratification axes, `status` included.
    # `covariate_names` drops `status`, an always-present model axis,
    # neither a required data column nor a name an observation can collide on.
    raw_covariate_names = function(value) {
      custom <- lapply(unlist(self$logic_tables, recursive = FALSE),
                       rate_cov_cols_)
      covs <- as.character(unlist(c(lapply(self$defs, `[[`, "by"), custom)))
      unique(covs[!is.na(covs)])
    },

    # `list` Statuses the rates key on, in order of first appearance,
    # see split_statuses_():
    # * automatic mode -> `reference` status;
    # * custom mode -> `status` column of the logic table.
    statuses = function(value) {
      st <- unlist(lapply(names(self$defs), function(rate_id) {
        ref <- self$defs[[rate_id]][["reference"]]
        if (self$defs[[rate_id]][["mode"]] == "auto")
          as.character(unlist(ref[names(ref) == "status"]))
        else lapply(self$logic_tables[[rate_id]],
                    \(tbl) as.character(tbl[["status"]]))
      }), use.names = FALSE)
      split_statuses_(unique(st[!is.na(st)]))
    },

    # `character()` Parameters required by group rates.
    param_required = function(value) {
      tables <- unlist(self$logic_tables, recursive = FALSE)
      required <- setdiff(unlist(lapply(tables, `[[`, "param_name")),
                          ".reference")
      set_source_name_(required, self$family)
    },

    # `data.frame()` Combined list with group rate and parameter definitions.
    tidy_table = function(value) {
      tidyup_group_rate_tables(self$logic_tables)
    }
  ),

  ####  Public  ####

  public = list(
    # `character(1)` Holder family.
    family = "group rate",

    # `character()` Cached overwritten names.
    store_warn_over = character(),

    # `list` Current group rate definitions.
    defs = list(),

    # `list` Raw parameter definitions.
    raw_params = list(),

    # `list` Per-parameter MCMC overrides keyed by rate name.
    param_overrides = list(),

    # `list` Parameter table of each group rate.
    logic_tables = list(),

    # Create an empty holder.
    initialize = function() {
      self$defs <- list()
      self$logic_tables <- list()
      self$raw_params <- list()
      self$param_overrides <- list()
      invisible(self)
    },

    # Remove the group rate definition `name`.
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

    # Whether some rule group of rate `name` has no reference, making it
    # absolute.
    is_absolute = function(name) {
      rate_is_absolute_(self$defs[[name]], self$logic_tables[[name]])
    },

    # Add a definition, arguments as in add_group_rate().
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

      config <- unwrap_rate_config_(prior, proposal, min, max, list(...),
                                    self$family, name)
      if (rate_mode_(by, strata, params) == "auto") {
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

      # store shared MCMC configuration (overwriting any existing)
      self$raw_params[name] <- list(config)

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

    # Build the logic tables from `dataset`, then assert every member fits.
    # * keyed on status, none declared -> error at the `wall`, else skipped.
    # @param status_codes `integer()` Engine status codes,
    # named after the labels the rates may stratify on.
    # @param stratified `character()` Rates read on infected members only.
    implement = function(dataset, status_codes = integer(),
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
      if ("status" %in% colnames(dataset))
        dataset[["status"]] <- status_as_labels(dataset[["status"]],
                                                status_codes)
      status_observed <- any(!is.na(dataset[["status"]]))

      for (rate_id in names(self$defs)) {
        # status-keyed rate -> waits for the statuses to be declared,
        # so decorator order stays free, the wall reporting one never declared
        if (!length(status_codes) && private$keys_on_status_(rate_id)) {
          if (isTRUE(wall)) stop_h(
            paste("A {self$family} is keyed on {.field status},",
                  "for which the model declares no level."),
            "Declare the infected statuses, or key the rate on a data column."
          )
          self$logic_tables[[rate_id]] <- list()
          next
        }
        if (self$defs[[rate_id]][["mode"]] == "auto") {
          private$helper_implement_auto_(rate_id, dataset,
                                         status_levels, status_observed)
        } else {
          private$helper_implement_custom_(rate_id, dataset, status_observed)
        }
        # the engine refuses a status key on a rate it does not stratify
        if (rate_id %in% stratified) {
          private$assert_rules_fit_(rate_id, dataset, status_levels)
        } else if (!private$keys_on_status_(rate_id)) {
          private$assert_rules_fit_(rate_id, dataset)
        }
      }

      invisible(self)
    },

    # Assert the model reads every group rate defined.
    # * required rate left undefined -> 1 for every member.
    assert = function(required) {
      notused <- definition_not_used_find(required, names(self$defs))
      if (length(notused)) stop_h(
        "Unused group rate{?s}: {.val {notused}}.",
        "The model reads no group rate of that name.",
        find_closest_match(
          notused, required,
          otherwise = "Print the model to show the group rates it reads."
        )
      )
      invisible(self)
    },

    # Encode the definitions for the engine.
    # @param info `list` of `list`s Class and levels of the data columns.
    # @param status_codes `integer()` Engine status codes, named by label.
    expose_to_cpp = function(info, status_codes) {
      defs <- self$tidy_table
      defs[rate_pin_cols] <- NULL
      # status has no entry in `info`, coded against the engine codes
      # members hold, not against a data column
      encode_strata_levels_(code_status_col_(defs, status_codes), info,
                            cols = rate_cov_cols_(defs))
    }
  )
)
