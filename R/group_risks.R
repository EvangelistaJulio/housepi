## housepi - Epidemiological Analysis of Household Outbreaks.
##
## Copyright © 2024-2026 Institut Pasteur (Paris), CNRS and INSERM.
## Distributed under the terms of the GNU General Public License (GPLv3).
## See the LICENSE file for details.
##

risk_stop_hint <- paste("Choose {.arg by} for automatic definitions,",
                        "{.arg strata} for custom.")

# Reserved column names (non-feature) in a logic table.
data_reserved_cols <- c("id", "house_id", "status")
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
  as.data.frame(risk_reserved_proto[cols], stringsAsFactors = FALSE)
}

# Typed `NA` padding for each of `cols`. A feature column pads as character.
risk_na_pad <- function(cols) {
  lapply(stats::setNames(cols, cols), function(col) {
    proto <- risk_reserved_proto[[col]]
    if (is.null(proto)) NA_character_ else proto[NA_integer_]
  })
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

#' Assert that features choice is valid.
#' @param features Features as list of data column names.
#' @param allow Reserved names the caller accepts as features.
#' @keywords internal
#' @noRd
assert_valid_features <- function(features, allow = character()) {
  if (all(is.na(features))) return(invisible())

  if (isFALSE(is.character(features))) stop_h(
    "Argument {.arg by} must be of type character.",
    format_type_error(features)
  )

  # avoid group risks associated to reserved data/definition names
  reserved_data <- setdiff(c(data_reserved_cols, risk_reserved_cols), allow)
  forbidden <- reserved_data[reserved_data %in% features]
  if (length(forbidden) > 0) stop_h(
    sprintf("Impossible to associate risk to %s %s",
            as_plural("feature", length(forbidden), "s:"),
            character_to_string_list(forbidden)),
    `!` = "These names are reserved.",
    paste("Please, change {as_plural('feature', length(forbidden))} and/or",
          "corresponding data {as_plural('column', length(forbidden))}.")
  )
}

#' Assert that features are present in data.
#' @param features Features as list of data column names.
#' @param headers All column names of the data.
#' @keywords internal
#' @noRd
assert_data_features <- function(features, headers) {
  # UPDATE: 'features' is not a useful terminology here
  missing <- setdiff(features, headers)
  if (length(missing) > 0) stop_h(
    sprintf("Unknown %s %s",
            as_plural("feature", length(missing), "s:"),
            character_to_string_list(missing)),
    "Please check your data headers."
  )
}

#' Assert that feature values in data are safe.
#' @param table Features as list of data column names.
#' @return `character()`
#' @keywords internal
#' @noRd
assess_feature_values <- function(this_table) {
  # column names to consider
  colmask <- setdiff(names(this_table), risk_reserved_cols)
  for (col in colmask) {
    sapply(this_table[[col]], function(x) assert_valid_feature(x, col))
  }
  colmask
}

#' Assert that a user-provided logic table is valid.
#' @param logic_table `data.frame` with a `param_name` column.
#' @param allow Reserved names the caller accepts as features.
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
  # validate feature column names
  feature_cols <- setdiff(colnames(logic_table), risk_reserved_cols)
  assert_valid_features(feature_cols, allow)
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
  if (!all(is.na(tbl$value)) && !is.numeric(tbl$value)) stop_h(
    "Column {.field value} of the strata table must be numeric.",
    format_type_error(tbl$value)
  )
  if (!all(is.na(tbl$const)) && !is.logical(tbl$const)) stop_h(
    "Column {.field const} of the strata table must be logical.",
    format_type_error(tbl$const)
  )
  value <- if (is.null(tbl$value)) rep(NA_real_, nrow(tbl)) else tbl$value
  const <- if (is.null(tbl$const)) rep(NA, nrow(tbl)) else tbl$const

  bad <- tbl$param_name == ".reference" & (!is.na(value) | !is.na(const))
  if (any(bad)) stop_h(
    sprintf("The {.val .reference} is pinned at %s %s.",
            as_plural("row", sum(bad)),
            character_to_phrase(as.character(which(bad)))),
    "The reference stratum is 1 by definition: it carries no pin."
  )
  bad <- const %in% TRUE & is.na(value)
  if (any(bad)) stop_h(
    sprintf("No {.field value} for fixed %s %s",
            as_plural("strat", sum(bad), "a", "um"),
            character_to_string_list(unique(tbl$param_name[bad]))),
    "A {.code const = TRUE} stratum needs the value it is fixed at."
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
  feat_cols <- setdiff(colnames(tbl), risk_reserved_cols)

  # pattern of each row: names of non-NA feature columns
  patterns <- lapply(seq_len(nrow(tbl)), function(i) {
    sort(feat_cols[!is.na(unlist(tbl[i, feat_cols, drop = FALSE]))])
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
    max_patterns, function(p) paste(p, collapse = "\x1F"), character(1)
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
      fmt <- function(x) {
        ifelse(length(x) == 0, "{{}}",
               paste0("{{", paste(x, collapse = ", "), "}}"))
      }
      stop_h(
        sprintf(
          "Ambiguous rule %d: pattern %s is subset of %d maximal patterns: %s.",
          i, fmt(p), length(candidates),
          paste(vapply(max_patterns[candidates], fmt, character(1)),
                collapse = " and ")
        ),
        "Add a disambiguating column to the rule."
      )
    }
  }
  tbl$rule_group <- gids
  tbl
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
  feat_cols <- setdiff(colnames(tbl), risk_reserved_cols)
  features <- sub("^from_", "", grep("^from_", feat_cols, value = TRUE))
  pair_wise <- isTRUE(symmetric) && length(features) > 0

  # cell identity: the rule itself, or the unordered pair when symmetric
  cell_key <- function(i) {
    row <- vapply(tbl[i, feat_cols, drop = FALSE], as.character, character(1))
    if (!pair_wise) return(paste(feat_cols, row, sep = "_", collapse = "\x1F"))
    ends <- lapply(c("from_", "to_"), function(side) {
      paste(features, row[paste0(side, features)], sep = "_", collapse = "\x1F")
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
    sprintf("%s %s of one rule group are tagged {.val .reference}.",
            as_plural("Row", length(clash)),
            character_to_phrase(sort(clash))),
    "A rule group has a single reference, setting its unit scale."
  )
  invisible(strata)
}

#' Pins carried by a logic table, keyed by `param_name`.
#' A stratum spread over several rules resolves to a single pin.
#' @param tbl `data.frame` logic table with a `param_name` column.
#' @return `list` of `list(value = , const = )`, one per pinned stratum.
#' @keywords internal
#' @noRd
collect_table_pins <- function(tbl) {
  cols <- intersect(risk_pin_cols, colnames(tbl))
  if (length(cols) == 0) return(list())
  pins <- list()
  for (par_id in setdiff(unique(tbl$param_name), ".reference")) {
    rows <- unique(tbl[tbl$param_name %in% par_id, cols, drop = FALSE])
    if (nrow(rows) > 1) stop_h(
      "Conflicting {.field {cols}} entries for stratum {.val {par_id}}.",
      "A stratum spread over several rules carries a single pin."
    )
    pin <- list(value = NA_real_, const = NA)
    pin[cols] <- as.list(rows)
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

#' Assert that no stratum is both pinned and overridden through `params`.
#' @param strata `data.frame` logic table with a `param_name` column.
#' @param params `list` Per-parameter MCMC-setup overrides.
#' @return `list` the input `params`, invisibly.
#' @keywords internal
#' @noRd
assert_no_pin_params_clash <- function(strata, params) {
  clash <- intersect(names(params), names(collect_table_pins(strata)))
  if (length(clash) > 0) stop_h(
    sprintf("Pinned %s also set through {.arg params}: %s",
            as_plural("strat", length(clash), "a", "um"),
            character_to_string_list(clash)),
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
#' @param this_feat Feature as list of data column names.
#' @param status_levels `character()` Declared infected-status labels.
#' @param status_observed `logical(1)` Whether the data carries any status.
#' @return data.frame with logic table
#' @keywords internal
#' @noRd
group_risk_logic_table <- function(data_wide, this_feat,
                                   status_levels = character(),
                                   status_observed = TRUE) {
  # extract unique features and sort by list length
  # status strata are declared by the model, not read off the data
  uni_features <- lapply(this_feat, function(f) {
    if (identical(f, "status")) status_levels else sort(unique(data_wide[[f]]))
  })
  sorted_indices <- order(sapply(uni_features, length), decreasing = TRUE)
  this_feat <- this_feat[sorted_indices]
  uni_features <- uni_features[sorted_indices]
  names(uni_features) <- this_feat

  # create a dataframe of all possible multi-features
  y <- expand.grid(uni_features)

  # count every combination in one pass. A cell the data does not carry
  # matches no stratum, so those members are left out of the tally
  strata_key <- function(x) do.call(paste, c(lapply(x[this_feat], as.character),
                                             list(sep = "\x1F")))
  carried <- data_wide[stats::complete.cases(data_wide[this_feat]), ,
                       drop = FALSE]
  tally <- table(strata_key(carried))
  multifeat_count <- as.integer(tally[strata_key(y)])
  multifeat_count[is.na(multifeat_count)] <- 0L

  # no status is observed: every declared stratum stands, sample sizes unknown
  if ("status" %in% this_feat && isFALSE(status_observed)) {
    y$n <- NA_integer_
    return(y)
  }

  # create the final truth table and reverse it
  y <- y[multifeat_count > 0, , drop = FALSE]
  # add counts
  y$n <- multifeat_count[multifeat_count > 0]

  y
}

#' Collapse a named reference into a phrase of `name = value` pairs.
#' @param reference `list()` Feature values, named after the features.
#' @param strip `character(1)` Prefix removed from the names before printing.
#' @return `character(1)`
#' @noRd
reference_to_phrase <- function(reference, strip = "^") {
  character_to_phrase(
    paste0("{.field ", sub(strip, "", names(reference)),
           "} = {.val ", unlist(reference), "}")
  )
}

#' Associate parameter name to group risk logic table.
#' @param this_table
#' Data.frame containing conditions for the group risk.
#' @param prefix
#' Prefix to parameter names associated to the group risk.
#' @param reference
#' Choose reference providing corresponding features.
#' `".auto"` picks the group with the largest sample size.
#' Empty leaves the risk without a reference.
#' @return
#' input data.frame with logic table with parameter names.
#' @keywords internal
#' @noRd
associate_group_param_names <- function(this_table, prefix, reference) {
  # column names to consider
  colmask <- assess_feature_values(this_table)

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
      `!` = sprintf("An individual with %s is never observed.",
                    reference_to_phrase(reference))
    )
  }

  # format prefix
  prefix <- ifelse(prefix == "" || !is.character(prefix),
                   "", paste0(prefix, "_"))

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
      sub_combined <- do.call(rbind, lapply(this_table, function(x) {
        miss <- setdiff(sub_cols, names(x))
        if (length(miss) > 0) x[miss] <- risk_na_pad(miss)
        x[sub_cols]
      }))
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
  combined_df <- do.call(rbind, lapply(df_list, function(df) {
    missing_cols <- setdiff(all_cols, names(df))
    if (length(missing_cols) > 0) df[missing_cols] <- risk_na_pad(missing_cols)
    df[all_cols]
  }))
  rownames(combined_df) <- NULL

  combined_df
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

    # Help defining the automatic reference.
    helper_auto_reference_ = function(by, reference) {
      if (identical(reference, ".auto")) return(".auto")

      if (length(by) != length(reference)) stop_h(
        "Arguments {.arg {c('by', 'reference')}} should have equal dimensions."
      )
      else if (length(setdiff(names(reference), by)) != 0) stop_h(
        "Unrecognized names in reference list."
      )
      reference <- as.list(reference)
      stats::setNames(reference, by)
    },

    # Add risk definition in automatic mode.
    helper_auto_add_ = function(name, by,
                                prefix = NA_character_,
                                reference = NULL,
                                joint = FALSE,
                                .warn_over = TRUE) {
      assert_valid_features(by, allow = "status")
      # remove previous definition if it exists
      rmv_name <- self$remove(name)

      if (.warn_over) template_overwrite_warn(rmv_name, self$family)
      else self$store_warn_over <- c(self$store_warn_over, rmv_name)

      # load reference and validate reference names
      if (length(reference) > 0) {
        reference <- private$helper_auto_reference_(by, reference)
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
        valid_names <- setdiff(strata$param_name, ".reference")
        unknown <- setdiff(names(params), valid_names)
        if (length(unknown) > 0) stop_h(
          sprintf("Unknown %s in {.arg params}: %s",
                  as_plural("parameter", length(unknown)),
                  character_to_string_list(unknown))
        )
        assert_no_pin_params_clash(strata, params)
      }

      # handle replace
      if (name %in% names(self$defs) &&
            self$defs[[name]][["mode"]] == "auto") {
        # remove previous definition
        self$logic_tables[[name]] <- NULL
        self$param_overrides[[name]] <- NULL
        rmv_name <- self$remove(name)
        if (.warn_over) template_overwrite_warn(rmv_name, self$family)
        else self$store_warn_over <- c(self$store_warn_over, rmv_name)
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
      features <- as.list(this_def[["by"]]) # data features

      if (all(is.na(features))) {
        self$logic_tables[[risk_id]] <- list()
        return(invisible(self))
      }

      features <- sort(as.character(features[!is.na(features)]))
      # status is carried by the model, not required as a data column
      assert_data_features(setdiff(features, "status"), colnames(data_wide))

      # cast to list if joint
      if (isTRUE(this_def[["joint"]])) features <- list(features)

      # default prefix
      prefix <- as.character(this_def[["prefix"]])
      if (is.na(prefix)) prefix <- substring(tolower(risk_id), 1, 3)

      new_logic_table <- list()
      for (this_feat in features) {
        # 1. create logic table with group risks
        this_table <- group_risk_logic_table(data_wide, this_feat,
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
        feat_cols <- setdiff(colnames(tbl), risk_reserved_cols)
        # status is carried by the model, not required as a data column
        assert_data_features(setdiff(feat_cols, "status"), colnames(data_wide))
        # validate that feature values exist in data
        for (col in feat_cols) {
          vals <- tbl[[col]][!is.na(tbl[[col]])]
          data_vals <- unique(data_wide[[col]])
          # status levels are declared by the model, not read off the data
          if (identical(col, "status")) data_vals <- status_levels
          missing_vals <- setdiff(as.character(vals), as.character(data_vals))
          if (length(missing_vals) > 0) warn_h(
            sprintf("Unknown %s for risk {.val %s} in data {.field %s}: %s",
                    as_plural("definition", length(missing_vals), "s"),
                    risk_id, col, character_to_string_list(missing_vals)),
            `!` = "Associated entries will be skipped."
          )
        }
        # add counts
        tbl$n <- apply(tbl[feat_cols], 1, function(row) {
          mask <- rep(TRUE, nrow(data_wide))
          for (col in feat_cols) {
            if (is.na(row[[col]])) next
            # `%in%`: an NA data cell matches no rule, it does not poison
            # the count of every stratum
            mask <- mask & (as.character(data_wide[[col]]) %in%
                              as.character(row[[col]]))
          }
          sum(mask)
        })
        # no status is observed: every declared stratum stands, sizes unknown
        if ("status" %in% feat_cols && isFALSE(status_observed)) {
          tbl$n <- NA_integer_
        } else {
          # skip entries with no matching data
          tbl <- tbl[tbl$n > 0, , drop = FALSE]
        }
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
        overrides <- self$param_overrides[[risk_id]]
        for (this_table in self$logic_tables[[risk_id]]) {
          pins <- collect_table_pins(this_table)
          for (par_id in setdiff(this_table$param_name, c(".reference"))) {
            def <- if (!is.null(overrides) && !is.null(overrides[[par_id]])) {
              overrides[[par_id]]
            } else {
              self$raw_params[[risk_id]]
            }
            new_params[[par_id]] <- helper_pin_param_def(def, pins[[par_id]])
          }
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

    #' @field feature_names `character()`. Data columns the risks stratify on,
    #' declared through `by` or carried by a custom logic table.
    feature_names = function(value) {
      setdiff(self$raw_feature_names, "status")
    },

    #' @field raw_feature_names `character()`. Stratification axes, `status`
    #' included. `feature_names` drops it: it is an always-present model axis,
    #' neither a required data column nor a name an observation can collide on.
    raw_feature_names = function(value) {
      custom <- lapply(unlist(self$logic_tables, recursive = FALSE),
                       function(tbl) setdiff(colnames(tbl), risk_reserved_cols))
      feats <- as.character(unlist(c(lapply(self$defs, `[[`, "by"), custom)))
      unique(feats[!is.na(feats)])
    },

    #' @field stratifies_on_status `logical(1)`.
    #' Whether a risk stratifies on the infected status.
    stratifies_on_status = function(value) {
      "status" %in% self$raw_feature_names
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
      !any(vapply(self$logic_tables[[name]],
                  function(tbl) ".reference" %in% tbl$param_name, logical(1)))
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
      reference <- housepi_smart_unwrap(reference, "reference")
      joint <- isTRUE(housepi_smart_unwrap(joint, "joint"))

      # custom definition
      strata <- housepi_smart_unwrap(strata, "strata")
      params <- housepi_smart_unwrap(params, "params")

      # raw parameter definition
      prior <- housepi_smart_unwrap(prior, "prior")
      proposal <- housepi_smart_unwrap(proposal, "proposal")
      min <- housepi_smart_unwrap(min, "min")
      max <- housepi_smart_unwrap(max, "max")
      dots <- list(...)
      if (length(dots) == 1 && rlang::is_quosures(dots[[1]])) dots <- dots[[1]]
      # extract args and proposal
      def_args <- helper_definition_args(proposal, dots)
      args <- def_args$args
      proposal <- def_args$proposal
      assert_no_sibling_dist_arg(args, self$family, name)

      # distinguish between automatic and custom definitions
      is_auto <- !all(is.na(by))
      is_custom <- !is.null(strata)

      if (is_auto && is_custom) stop_ambiguos_definition()
      else if (is_auto) private$helper_auto_add_(
        name = name, by = by, reference = reference, prefix = prefix,
        joint = joint, .warn_over = .warn_over
      )
      else if (is_custom) private$helper_custom_add_(
        name = name, strata = strata, params = params,
        reference = reference, .warn_over = .warn_over
      )
      else stop_empty_definition()

      # store shared MCMC setup (eventually overwrite)
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
        sms <- c(
          sprintf("Dropped unused group %s %s",
                  as_plural("risk", length(notused), "s:"),
                  character_to_string_list(notused)),
          find_closest_match(notused, required)
        )
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
      for (k in intersect(colnames(defs), names(info))) {
        # `exclude = NULL`: an NA cell is a level the engine codes like any
        # other, and dropping it here would shift every code by one
        defs[[k]] <- as.integer(
          factor(defs[[k]], levels = info[[k]]$levels, exclude = NULL)
        )
      }
      defs
    }
  )
)
