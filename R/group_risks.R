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
risk_reserved_cols <- c("group_risk", "param_name", "n", "rule_group")

# Empty logic table with the reserved columns `cols`, in order.
empty_risk_table <- function(cols) {
  proto <- list(group_risk = character(), param_name = character(),
                n = integer(), rule_group = integer())
  as.data.frame(proto[cols], stringsAsFactors = FALSE)
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
#' @keywords internal
#' @noRd
assert_valid_features <- function(features) {
  if (all(is.na(features))) return(invisible())

  if (isFALSE(is.character(features))) stop_h(
    "Argument {.arg by} must be of type character.",
    format_type_error(features)
  )

  # avoid group risks associated to reserved data/definition names
  reserved_data <- c(data_reserved_cols, risk_reserved_cols)
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
  # assert that features are present
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
#' @keywords internal
#' @noRd
assert_valid_logic_table <- function(logic_table) {
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
  assert_valid_features(feature_cols)
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

#' Return logic table of group risks
#' @param data_wide Data from which to extract group risks
#' @param this_feat Feature as list of data column names.
#' @return data.frame with logic table
#' @keywords internal
#' @noRd
group_risk_logic_table <- function(data_wide, this_feat) {
  # extract unique features and sort by list length
  uni_features <- lapply(this_feat, function(f) sort(unique(data_wide[[f]])))
  sorted_indices <- order(sapply(uni_features, length), decreasing = TRUE)
  this_feat <- this_feat[sorted_indices]
  uni_features <- uni_features[sorted_indices]
  names(uni_features) <- this_feat

  # create a dataframe of all possible multi-features
  y <- expand.grid(uni_features)

  # check that the multi-feature exists
  multifeat_count <- apply(y, 1, function(row) {
    sum(apply(data_wide[this_feat], 1, function(d_row) all(d_row == row)))
  })

  # create the final truth table and reverse it
  y <- y[multifeat_count > 0, , drop = FALSE]
  # add counts
  y$n <- multifeat_count[multifeat_count > 0]

  y
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
    idx_refer <- which.max(this_table$n)

  } else if (length(reference) > 0) {
    # find definition row
    idx_refer <- unname(which(
      apply(this_table[colmask], 1,
            function(row) all(row == unlist(reference[colmask])))
    ))
    if (length(idx_refer) != 1) stop_h("Provided reference group is not valid.")
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
        x[miss] <- NA_character_
        x[sub_cols]
      }))
      sub_combined <- compute_rule_groups(sub_combined)
      sub_combined <- reorder_risks_logic_table(sub_combined)
      df <- list(cbind_risk(sub_combined, risk))
    }
    df_list <- c(df_list, df)
  }

  # all unique column names reordered
  all_cols <- unique(c(risk_reserved_cols,
                       sort(unlist(lapply(df_list, names)))))

  # Add missing columns with NA
  combined_df <- do.call(rbind, lapply(df_list, function(df) {
    missing_cols <- setdiff(all_cols, names(df))
    df[missing_cols] <- NA_character_
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
      assert_valid_features(by)
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
      assert_valid_logic_table(strata)
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
                    "param_name")
      self$defs[[name]] <- list(by = by, prefix = NA_character_,
                                reference = NULL, joint = FALSE,
                                mode = "custom")
      invisible(self)
    },

    # Implement through automatic mode.
    helper_implement_auto_ = function(risk_id, data_wide, factor_info) {
      this_def <- self$defs[[risk_id]] # definition of the risk
      features <- as.list(this_def[["by"]]) # data features

      if (all(is.na(features))) {
        self$logic_tables[[risk_id]] <- list()
        return(invisible(self))
      }

      features <- sort(as.character(features[!is.na(features)]))
      assert_data_features(features, colnames(data_wide))

      # cast to list if joint
      if (isTRUE(this_def[["joint"]])) features <- list(features)

      # default prefix
      prefix <- as.character(this_def[["prefix"]])
      if (is.na(prefix)) prefix <- substring(tolower(risk_id), 1, 3)

      new_logic_table <- list()
      for (this_feat in features) {
        # 1. create logic table with group risks
        this_table <- group_risk_logic_table(data_wide, this_feat)
        # associate names for parameters
        this_table <- associate_group_param_names(
          this_table, prefix, this_def$reference
        )
        new_logic_table[[length(new_logic_table) + 1]] <- this_table

        # 2. create relative cpp script chunk
        # factorize according to external factor information
        for (col in intersect(names(factor_info), colnames(this_table))) {
          lv <- factor_info[[col]]$levels
          if (is.null(lv)) stop_h(
            "Wrong group risk definition: {.field {col}} has no levels."
          )
          this_table[[col]] <- factor(this_table[[col]],
                                      levels = lv,
                                      labels = seq_along(lv))
        }
        # concatenate new lines to the script
        # mask reference case
        mask <- this_table$param_name == ".reference"

        # store the new risk logic table
        self$logic_tables[[risk_id]] <- new_logic_table
      }
    },

    # Validate custom mode logic table against data, then add counts.
    helper_implement_custom_ = function(risk_id, data_wide) {
      for (ti in seq_along(self$logic_tables[[risk_id]])) {
        tbl <- self$logic_tables[[risk_id]][[ti]]
        feat_cols <- setdiff(colnames(tbl), risk_reserved_cols)
        # validate that feature columns exist in data
        assert_data_features(feat_cols, colnames(data_wide))
        # validate that feature values exist in data
        for (col in feat_cols) {
          vals <- tbl[[col]][!is.na(tbl[[col]])]
          data_vals <- unique(data_wide[[col]])
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
            mask <- mask & (as.character(data_wide[[col]]) ==
                              as.character(row[[col]]))
          }
          sum(mask)
        })
        # skip entries with no matching data
        tbl <- tbl[tbl$n > 0, , drop = FALSE]
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
          for (par_id in setdiff(this_table$param_name, c(".reference"))) {
            if (!is.null(overrides) && !is.null(overrides[[par_id]])) {
              new_params[[par_id]] <- overrides[[par_id]]
            } else {
              new_params[[par_id]] <- self$raw_params[[risk_id]]
            }
          }
        }
      }
      new_params
    },

    #' @field data_required `character()`. Data required by group risks.
    data_required = function(value) {
      required <- unique(unlist(lapply(self$defs, `[[`, "by")))
      required <- required[!is.na(required)]
      set_source_name_(required, self$family)
    },

    #' @field feature_names `character()`. Data columns the risks stratify on,
    #' declared through `by` or carried by a custom logic table.
    feature_names = function(value) {
      custom <- lapply(unlist(self$logic_tables, recursive = FALSE),
                       function(tbl) setdiff(colnames(tbl), risk_reserved_cols))
      feats <- as.character(unlist(c(lapply(self$defs, `[[`, "by"), custom)))
      unique(feats[!is.na(feats)])
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
    #' @param factor_info
    #' A `list` with the levels for each factorized feature.
    implement = function(data_wide, factor_info) {
      # preserve custom-mode entries and drop auto-mode entries
      mask_custom <- names(which(vapply(
        self$defs, function(d) d[["mode"]] == "custom", logical(1)
      )))
      self$logic_tables <- self$logic_tables[mask_custom]

      for (risk_id in names(self$defs)) {
        if (self$defs[[risk_id]][["mode"]] == "auto") {
          private$helper_implement_auto_(risk_id, data_wide, factor_info)
        } else {
          private$helper_implement_custom_(risk_id, data_wide)
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
    expose_to_cpp = function(info) {
      defs <- self$tidy_table
      defs$n <- NULL
      for (k in intersect(colnames(defs), names(info))) {
        defs[[k]] <- as.integer(factor(defs[[k]], levels = info[[k]]$levels))
      }
      defs
    }
  )
)
