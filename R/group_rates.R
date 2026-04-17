## housepi - Epidemiological Analysis of Household Outbreaks.
##
## Copyright © 2026 Institut Pasteur (Paris), CNRS and INSERM.
## Distributed under the terms of the GNU General Public License (GPLv3).
## See the LICENSE file for details.
##

rate_stop_hint <- "Choose 'by' for automatic definitions, 'strata' for custom."

# Reserved column names (non-feature) in a logic table.
rate_reserved_cols <- c("group_rate", "param_name", "n", "rule_group")

# Throw error if empty definition.
stop_empty_definition <- function() {
  stop_h("Invalid rate definition: provide one between 'by' and 'strata'.",
         rate_stop_hint)
}

# Throw error for ambiguos rate definition.
stop_ambiguos_definition <- function() {
  stop_h("Impossible to provide at the same time 'by' and 'logic_table'.",
         rate_stop_hint)
}

#' Assert that features choice is valid.
#' @param features Features as list of data column names.
#' @keywords internal
#' @noRd
assert_valid_features <- function(features) {
  if (all(is.na(features))) return(invisible())

  if (isFALSE(is.character(features))) stop_h(
    "Argument 'by' must be of type character.",
    format_type_error(features)
  )

  # avoid group rates associated to reserved data/definition names
  reserved_data <- c("id", "house_id", "status", rate_reserved_cols)
  forbidden <- reserved_data[reserved_data %in% features]
  if (length(forbidden) > 0) stop_h(
    sprintf("Impossible to associate rate to feature%s %s",
            as_plural(forbidden, "s:"),
            character_to_string_list(forbidden)),
    "These names are reserved.",
    sprintf("Please, change feature%s and/or corresponding data column%s.",
            as_plural(forbidden), as_plural(forbidden))
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
    sprintf("Unknown feature%s %s",
            as_plural(missing, "s:"),
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
  colmask <- setdiff(names(this_table), rate_reserved_cols)
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
    "Argument 'logic_table' must be a data.frame."
  )
  if (!"param_name" %in% colnames(logic_table)) stop_h(
    "Logic table must contain a 'param_name' column."
  )
  if (nrow(logic_table) == 0) stop_h(
    "Logic table must contain at least one row."
  )
  # validate feature column names
  feature_cols <- setdiff(colnames(logic_table), rate_reserved_cols)
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
  feat_cols <- setdiff(colnames(tbl), rate_reserved_cols)

  # pattern of each row: names of non-NA feature columns
  patterns <- lapply(seq_len(nrow(tbl)), function(i) {
    sort(feat_cols[!is.na(unlist(tbl[i, feat_cols, drop = FALSE]))])
  })

  # unique patterns for maximality check
  uniq_patterns <- unique(patterns)
  is_max <- vapply(uniq_patterns, function(p) {
    !any(vapply(uniq_patterns, function(q)
      length(q) > length(p) && all(p %in% q), logical(1)))
  }, logical(1))
  max_patterns <- uniq_patterns[is_max]

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
        ifelse(length(x) == 0, "{}",
               paste0("{", paste(x, collapse = ", "), "}"))
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

#' Return logic table of group rates
#' @param data_wide Data from which to extract group rates
#' @param this_feat Feature as list of data column names.
#' @return data.frame with logic table
#' @keywords internal
#' @noRd
group_rate_logic_table <- function(data_wide, this_feat) {
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

#' Associate parameter name to group rate logic table.
#' @param this_table
#' Data.frame containing conditions for the group rate.
#' @param prefix
#' Prefix to parameter names associated to the group rate.
#' @param reference
#' Choose reference providing corresponding features.
#' Default is the group with the largest sample size.
#' @return
#' input data.frame with logic table with parameter names.
#' @keywords internal
#' @noRd
associate_param_names <- function(this_table, prefix, reference) {
  # column names to consider
  colmask <- assess_feature_values(this_table)

  # extract reference index
  if (length(reference) == 0) {
    idx_refer <- which.max(this_table$n)

  } else {
    # find row
    idx_refer <- unname(which(
      apply(this_table[colmask], 1,
            function(row) all(row == unlist(reference[colmask])))
    ))

    if (length(idx_refer) != 1) stop_h("Provided reference group is not valid.")
  }

  # format prefix and associate parameter names
  prefix <- ifelse(
    (prefix == "") || (!is.character(prefix)),
    "", paste0(prefix, "_")
  )

  # associate parameter names
  this_table$param_name <- apply(
    this_table[colmask], 1, function(row) {
      paste0(prefix, paste(paste(colmask, row, sep = "_"), collapse = "_"))
    }
  )
  this_table$param_name[idx_refer] <- ".reference"

  this_table
}

#' Combine list with group rate definition into a single dataframe
#' @param all_tables `list`. Definitions for each group rate.
#' @return `data.frame`
#' @keywords internal
#' @noRd
tidyup_group_rate_tables <- function(all_tables) {
  # tidy version of group rate tables
  cbind_rate <- function(x, rate) cbind(group_rate = rate, x)
  df_list <- c()
  for (rate in names(all_tables)) {
    this_table <- all_tables[[rate]]
    if (length(this_table) == 0) {
      df <- list(data.frame(group_rate = rate, param_name = ".reference",
                            rule_group = 1L))
    } else {
      # concat this rate's sub-tables, then compute rule_groups within rate
      sub_cols <- unique(c("param_name", "n",
                           sort(unlist(lapply(this_table, names)))))
      sub_combined <- do.call(rbind, lapply(this_table, function(x) {
        miss <- setdiff(sub_cols, names(x))
        x[miss] <- NA_character_
        x[sub_cols]
      }))
      sub_combined <- compute_rule_groups(sub_combined)
      df <- list(cbind_rate(sub_combined, rate))
    }
    df_list <- c(df_list, df)
  }

  # all unique column names reordered
  all_cols <- unique(c(rate_reserved_cols,
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
#' @references
#' Cauchemez, et al. (2009)
#' Household Transmission of 2009 Pandemic Influenza A (H1N1)
#' Virus in the United States.
#' New England Journal of Medicine. 361: 2619-2627.
#' \doi{10.1056/NEJMoa0905498}
#' @return `list`
#' @keywords internal
#' @noRd
set_group_rates_mcmc_setup <- function(prior, args, min, max,
                                       proposal, proposal_scale) {
  # check wheter default configuration
  is_default <- sapply(c(prior, args, min, max, proposal, proposal_scale),
                       function(x) is.na(x) || is.null(x))
  # if all defaults, provide default choice
  if (all(is_default)) list(
    prior = "lognormal", meanlog = 0, sdlog = 3,
    min = 0, max = 10, proposal = "lognormal"
  )
  # else go with user choice
  else c(list(prior = prior), args,
         list(min = min, max = max,
              proposal = proposal, proposal_scale = proposal_scale))
}

############################################################################
############################################################################
############################# Housepi Group Rate ###########################
############################################################################
############################################################################

#' @title Housepi Group Rate
#' @description Help with group rate definitions of model Rcpp object.
#' @keywords internal
#' @noRd
HousepiGroupRate <- R6::R6Class(
  classname = "HousepiGroupRate",
  cloneable = TRUE,

  ############################################################################
  ################################## PRIVATE #################################
  ############################################################################

  private = list(
    deep_clone = function(name, value) helper_clone(value),

    # Add rate definition in automatic mode.
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

      # load reference
      if (!is.null(reference)) {
        if (length(by) != length(reference)) stop_h(
          "Arguments 'by' and 'reference' should have the same dimension."
        )
        else if (length(setdiff(names(reference), by)) != 0) stop_h(
          "Unrecognized names in reference list."
        )
        reference <- as.list(reference)
        reference <- stats::setNames(reference, by)
      }

      # add new definition
      self$defs[[name]] <- list(by = by, reference = reference,
                                prefix = prefix, joint = joint,
                                mode = "auto")
      invisible(self)
    },

    # Add rate definition in custom mode.
    helper_custom_add_ = function(name, strata,
                                  params = NULL,
                                  .warn_over = TRUE) {
      assert_valid_logic_table(strata)
      # validate params keys against logic table param_names
      if (length(params) > 0) {
        valid_names <- setdiff(strata$param_name, ".reference")
        unknown <- setdiff(names(params), valid_names)
        if (length(unknown) > 0) stop_h(
          sprintf("Unknown parameter%s %s in 'params'.",
                  as_plural(unknown),
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
    helper_implement_auto_ = function(rate_id, data_wide, factor_info) {
      this_def <- self$defs[[rate_id]] # definition of the rate
      features <- as.list(this_def[["by"]]) # data features

      if (all(is.na(features))) {
        self$logic_tables[[rate_id]] <- list()
        return(invisible(self))
      }

      features <- sort(as.character(features[!is.na(features)]))
      assert_data_features(features, colnames(data_wide))

      # cast to list if joint
      if (isTRUE(this_def[["joint"]])) features <- list(features)

      # default prefix
      prefix <- as.character(this_def[["prefix"]])
      if (is.na(prefix)) prefix <- substring(tolower(rate_id), 1, 3)

      new_logic_table <- list()
      for (this_feat in features) {
        # 1. create logic table with group rates
        this_table <- group_rate_logic_table(data_wide, this_feat)
        # associate names for parameters
        this_table <- associate_param_names(
          this_table, prefix, this_def$reference
        )
        new_logic_table[[length(new_logic_table) + 1]] <- this_table

        # 2. create relative cpp script chunk
        # factorize according to external factor information
        for (col in names(factor_info)) {
          lv <- factor_info[[col]]$levels
          if (is.null(lv)) stop_h(
            "Wrong group rate definition: '", col, "' has no levels."
          )

          if (col %in% colnames(this_table)) {
            this_table[[col]] <- factor(this_table[[col]],
                                        levels = lv,
                                        labels = seq_along(lv))
          }
        }
        # concatenate new lines to the script
        # mask reference case
        mask <- this_table$param_name == ".reference"

        # store the new rate logic table
        self$logic_tables[[rate_id]] <- new_logic_table
      }
    },

    # Validate custom mode logic table against data, then add counts.
    helper_implement_custom_ = function(rate_id, data_wide) {
      for (ti in seq_along(self$logic_tables[[rate_id]])) {
        tbl <- self$logic_tables[[rate_id]][[ti]]
        feat_cols <- setdiff(colnames(tbl), rate_reserved_cols)
        # validate that feature columns exist in data
        assert_data_features(feat_cols, colnames(data_wide))
        # validate that feature values exist in data
        for (col in feat_cols) {
          vals <- tbl[[col]][!is.na(tbl[[col]])]
          data_vals <- unique(data_wide[[col]])
          missing_vals <- setdiff(as.character(vals), as.character(data_vals))
          if (length(missing_vals) > 0) warning_h(
            sprintf("Rate '%s' definition value%s not found in data '%s': %s",
                    rate_id, as_plural(missing_vals, "s"), col,
                    character_to_string_list(missing_vals)),
            "Corresponding logic table entry will be skipped."
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
        self$logic_tables[[rate_id]][[ti]] <- tbl
      }
      invisible(self)
    }
  ),

  ############################################################################
  ################################## ACTIVE ##################################
  ############################################################################

  active = list(
    #' @field parameters `list`.
    #' Parameter definitions associated to group rates.
    parameters = function(value) {
      if (!missing(value)) stop_read_only("parameters")
      new_params <- list()
      for (rate_id in names(self$logic_tables)) {
        overrides <- self$param_overrides[[rate_id]]
        for (this_table in self$logic_tables[[rate_id]]) {
          for (par_id in setdiff(this_table$param_name, c(".reference"))) {
            if (!is.null(overrides) && !is.null(overrides[[par_id]])) {
              new_params[[par_id]] <- overrides[[par_id]]
            } else {
              new_params[[par_id]] <- self$raw_params[[rate_id]]
            }
          }
        }
      }
      new_params
    },

    #' @field data_required `character()`. Data required by group rates.
    data_required = function(value) {
      if (!missing(value)) stop_read_only("data_required")
      required <- unique(unlist(lapply(self$defs, `[[`, "by")))
      required <- required[!is.na(required)]
      set_source_name_(required, self$family)
    },

    #' @field param_required `character()`. Parameters required by group rates.
    param_required = function(value) {
      if (!missing(value)) stop_read_only("param_required")
      required <- c()
      for (rate_id in names(self$logic_tables)) {
        for (this_table in self$logic_tables[[rate_id]]) {
          required <- c(required, this_table$param_name)
        }
      }
      required <- setdiff(unique(required), c(".reference"))
      set_source_name_(required, self$family)
    },

    #' @field tidy_table `data.frame()`.
    #' Combined list with group rate and parameter definitions.
    tidy_table = function(value) {
      if (!missing(value)) stop_read_only("tidy_table")
      tidyup_group_rate_tables(self$logic_tables)
    }
  ),

  ############################################################################
  ################################## PUBLIC ##################################
  ############################################################################

  public = list(
    #' @field family `character(1)`. Module family.
    family = "group rate",

    #' @field store_warn_over `character()`. Cached overwritten names.
    store_warn_over = character(),

    #' @field defs `list`. Current group rate definitions.
    defs = list(),

    #' @field raw_params `list`. Raw parameter definitions.
    raw_params = list(),

    #' @field param_overrides `list`.
    #' Per-parameter MCMC overrides keyed by rate name.
    param_overrides = list(),

    #' @field logic_tables `list`.
    #' Table of parameters associated to each group rate.
    logic_tables = list(),

    #' @description
    #' New object of class [HousepiGroupRate][housepi::HousepiGroupRate].
    initialize = function() {
      self$defs <- list()
      self$logic_tables <- list()
      self$raw_params <- list()
      self$param_overrides <- list()
      invisible(self)
    },

    #' @description
    #' Find group rates required by the model.
    #' @param refer_file Path to reference file.
    require = function(refer_file) {
      output <- template_extract_keys(refer_file, "GroupRateKey")
      set_source_name_(output, "model")
    },

    #' @description
    #' Remove a given group_rate definition.
    #' @param name `character(1)`. Name of group rate definition to remove.
    #' @param warn_removal `logical(1)`. Show warnings for failed removal.
    remove = function(name, warn_removal = FALSE) {
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

      } else if (isTRUE(warn_removal)) {
        warning_h(sprintf("Attempt to remove unknown group_rate '%s'.", name))
      }

      c()
    },

    #' @inherit add_group_rate
    add_quo = function(name,
                       by = NA_character_, prefix = NA_character_,
                       reference = NULL, joint = FALSE,
                       strata = NULL, params = NULL,
                       prior = NA_character_, ...,
                       min = NA_real_, max = NA_real_,
                       proposal = NA_character_, proposal_scale = NA_real_,
                       .warn_over = TRUE) {
      name <- housepi_smart_unwrap(name, "name")
      by <- housepi_smart_unwrap(by, "by")
      # control valid name # NOTE: `by` is the injected default variable
      assert_valid_definition_name(name, self$family, !is.na(by))

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
      proposal_scale <- housepi_smart_unwrap(proposal_scale, "proposal_scale")
      min <- housepi_smart_unwrap(min, "min")
      max <- housepi_smart_unwrap(max, "max")
      dots <- list(...)
      if (length(dots) == 1 && rlang::is_quosures(dots[[1]])) dots <- dots[[1]]
      args <- housepi_extractor(!!!dots)

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
        .warn_over = .warn_over
      )
      else stop_empty_definition()

      # store shared MCMC setup (eventually overwrite)
      self$raw_params[[name]] <- set_group_rates_mcmc_setup(
        prior = prior, args = args, min = min, max = max,
        proposal = proposal, proposal_scale = proposal_scale
      )

      # reorder
      y <- sort(names(self$defs))
      self$defs <- self$defs[y]
      self$raw_params <- self$raw_params[y]

      invisible(self)
    },

    #' @inherit load_group_rate_list
    load_quosures = function(xquos, append = FALSE) {
      if (isFALSE(append)) self$initialize()
      housepi_load_defin_quos(self, "by", xquos)
      invisible(self)
    },

    ############################################################################
    ############################################################################
    ############################################################################

    #' @description
    #' Implement group rate definition according to data.
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

      for (rate_id in names(self$defs)) {
        if (self$defs[[rate_id]][["mode"]] == "auto") {
          private$helper_implement_auto_(rate_id, data_wide, factor_info)
        } else {
          private$helper_implement_custom_(rate_id, data_wide)
        }
      }

      invisible(self)
    },

    #' @description
    #' Assert that group rate definitions match the requirements.
    #' @param required `character()`. Required group rates.
    #' @param capital `logical(1)`. Whether to raise error.
    assert = function(required, capital) {
      # load definitions
      defined <- names(self$defs)

      # UPDATE
      # point out group rates that are not used
      notused <- template_not_used_find(required, defined)
      # drop unused group rates
      if (length(notused) > 0) {
        sms <- c(
          sprintf("Dropped unused group rate%s %s",
                  as_plural(notused, "s:"),
                  character_to_string_list(notused)),
          find_closest_match(notused, required, min_dist = 3)
        )
        if (isTRUE(capital)) {
          for (which in notused) self$remove(which)
          stop_h(sms)
        }
        warning_h(sms)
      }

      template_missing_stop(required, defined, self$family)
    },

    #' @description
    #' Tidy table of group rates definition to be passed to Cpp.
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
