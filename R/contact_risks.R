## housepi - Epidemiological Analysis of Household Outbreaks.
##
## Copyright © 2024-2026 Institut Pasteur (Paris), CNRS and INSERM.
## Distributed under the terms of the GNU General Public License (GPLv3).
## See the LICENSE file for details.
##

#' Return logic table of contact risks (from/to paired columns).
#' @param data_wide Data from which to extract contact risks.
#' @param features Features as list of data column names.
#' @return data.frame with from_X, to_X columns and count n.
#' @keywords internal
#' @noRd
contact_risk_logic_table <- function(data_wide, features) {
  # convention: no house_id implies households of size 1: no pairs possible
  if (!"house_id" %in% colnames(data_wide)) return(NULL)
  # use existing "id" column or assign a row id
  if (!"id" %in% colnames(data_wide)) data_wide$id <- seq_len(nrow(data_wide))

  # subset to needed column, and coerce features to character
  features <- as.character(unlist(features))
  data_wide <- data_wide[, c("house_id", "id", features), drop = FALSE]
  data_wide[] <- lapply(data_wide, as.character)

  # vectorised self-join to retrieve every within-household ordered pair
  pairs <- merge(data_wide, data_wide, by = "house_id",
                 suffixes = c(".from", ".to"))
  pairs <- pairs[pairs$id.from != pairs$id.to, , drop = FALSE]

  # hash using unprintable separators '\x1E' and '\x1F'
  make_key <- function(x, cols) do.call(paste, c(x[cols], list(sep = "\x1F")))

  # hash-based counting of pair types
  pair_key <- paste(make_key(pairs, paste0(features, ".from")),
                    make_key(pairs, paste0(features, ".to")),
                    sep = "\x1E")
  cnt <- table(pair_key)

  # build the full observed grid
  observed <- unique(data_wide[, features, drop = FALSE])
  from_grid <- stats::setNames(observed, paste0("from_", features))
  to_grid <- stats::setNames(observed, paste0("to_", features))
  y <- merge(from_grid, to_grid, by = NULL)
  rownames(y) <- NULL

  # look up counts by the same key
  grid_key <- paste(make_key(y, names(from_grid)),
                    make_key(y, names(to_grid)),
                    sep = "\x1E")
  y$n <- as.integer(cnt[grid_key])
  y$n[is.na(y$n)] <- 0L
  y[y$n > 0, ]
}

#' Associate parameter names to contact risk logic table.
#' @param this_table `data.frame` with from_X and to_X columns.
#' @param prefix Prefix for parameter names.
#' @param reference Named vector with from_/to_ prefixed names and values.
#' @return `data.frame` with param_name column added.
#' @keywords internal
#' @noRd
associate_contact_param_names <- function(this_table, prefix, reference,
                                          symmetric = TRUE) {
  colmask <- assess_feature_values(this_table)
  from_cols <- grep("^from_", colmask, value = TRUE)
  to_cols <- grep("^to_", colmask, value = TRUE)
  features <- sub("^from_", "", from_cols)

  # find reference index
  idx_refer <- NULL
  if (length(reference) == 0) {
    idx_refer <- which.max(this_table$n)

  } else if (!all(reference == ".none")) {
    # find definition row
    refer_mask <- apply(this_table[colmask], 1, function(row) {
      all(row == unlist(reference[colmask]))
    })
    idx_refer <- unname(which(refer_mask))
    if (length(idx_refer) != 1) stop_h(
      "Provided reference contact risk is not valid.",
      `!` = sprintf(
        "The pair %s %s %s is never observed.",
        character_to_phrase(reference$from, quote = TRUE, and = FALSE),
        ifelse(symmetric, "<->", "-->"),
        character_to_phrase(reference$to, quote = TRUE, and = FALSE)
      )
    )
  }

  # format prefix
  prefix <- ifelse(prefix == "" || !is.character(prefix),
                   "", paste0(prefix, "_"))

  # associate param names with one "from" and one "to" token
  this_table$param_name <- apply(
    this_table[colmask], 1, function(row) {
      from_part <- paste(features, row[from_cols], sep = "_", collapse = "_")
      to_part <- paste(features, row[to_cols], sep = "_", collapse = "_")
      if (isTRUE(symmetric) && from_part > to_part) {
        # order convention: smaller key is always "from"
        tmp <- from_part
        from_part <- to_part
        to_part <- tmp
      }
      paste0(prefix, from_part, "_to_", to_part)
    }
  )
  if (length(idx_refer) == 1) {
    if (isTRUE(symmetric)) {
      ref_name <- this_table$param_name[idx_refer]
      this_table$param_name[this_table$param_name == ref_name] <- ".reference"
    } else {
      this_table$param_name[idx_refer] <- ".reference"
    }
  }

  this_table
}

#' Fill symmetric gaps in a custom contact risk logic table.
#' For each row, if the reverse from/to pair is missing,
#' add it with the same param_name.
#' @param tbl `data.frame` with from_X, to_Y, and param_name columns.
#' @return `data.frame` with missing reverse pairs added.
#' @keywords internal
#' @noRd
fill_symmetric_gaps <- function(tbl) {
  feat_cols <- setdiff(colnames(tbl), risk_reserved_cols)
  from_cols <- grep("^from_", feat_cols, value = TRUE)
  features <- sub("^from_", "", from_cols)

  new_rows <- list()
  for (i in seq_len(nrow(tbl))) {
    rev_row <- tbl[i, , drop = FALSE]
    for (f in features) {
      rev_row[[paste0("from_", f)]] <- tbl[i, paste0("to_", f)]
      rev_row[[paste0("to_", f)]] <- tbl[i, paste0("from_", f)]
    }

    def_exists <- any(apply(tbl[feat_cols], 1, function(row) {
      rv <- unlist(rev_row[feat_cols])
      all(mapply(function(a, b) (is.na(a) && is.na(b)) || isTRUE(a == b),
                 row, rv))
    }))
    if (!def_exists) new_rows[[length(new_rows) + 1]] <- rev_row
  }

  if (length(new_rows) > 0) {
    tbl <- rbind(tbl, do.call(rbind, new_rows))
    rownames(tbl) <- NULL
  }
  tbl
}

############################################################################
############################################################################
########################### Housepi Contact Risk ###########################
############################################################################
############################################################################

#' Combine list of contact risk sub-tables into a single data.frame.
#' @param tables `list` of `data.frame`s.
#' @return `data.frame` with columns unified across sub-tables.
#' @keywords internal
#' @noRd
tidyup_contact_risk_tables <- function(tables) {
  if (length(tables) == 0) return(data.frame())

  all_cols <- unique(c("param_name", "n",
                       sort(unlist(lapply(tables, names)))))

  combined <- do.call(rbind, lapply(seq_along(tables), function(i) {
    df <- tables[[i]]
    missing_cols <- setdiff(all_cols, names(df))
    df[missing_cols] <- NA_character_
    df[all_cols]
  }))
  rownames(combined) <- NULL
  combined <- compute_rule_groups(combined)
  combined <- reorder_risks_logic_table(combined)
  rownames(combined) <- NULL
  combined
}

#' @title Internal R6 holder for housepi contact-risk definitions
#' @description
#' Helper class managing contact-risk definitions of the model object.
#' @keywords internal
#' @noRd
HousepiContactRisk <- R6::R6Class(
  classname = "HousepiContactRisk",
  cloneable = TRUE,

  ############################################################################
  ################################## PRIVATE #################################
  ############################################################################

  private = list(
    deep_clone = function(name, value) helper_clone(value),

    # Help defining the automatic reference.
    helper_auto_reference_ = function(by, reference) {
      if (isTRUE(any(reference == ".none"))) {
        if (isTRUE(all(reference == ".none"))) return(".none")
        else stop_h("Ambigous entry {.val .none} in reference vector.")
      }

      valid_names <- c(paste0("from_", sort(by)), paste0("to_", sort(by)))
      # convenience: allow un-prefixed names
      if (length(reference) == length(by) && all(names(reference) %in% by)) {
        reference <- c(stats::setNames(reference[by], paste0("from_", by)),
                       stats::setNames(reference[by], paste0("to_", by)))
      }
      wrong_entries <- setdiff(names(reference), valid_names)
      if (length(wrong_entries) != 0) stop_h(
        sprintf("Unrecognized %s for reference: %s",
                as_plural("choice", length(wrong_entries)),
                character_to_string_list(wrong_entries)),
        paste("Reference should match {.arg by}",
              "and use {.code from_}\\{.code to_} prefixes.")
      )
      as.list(reference)
    },

    # Add risk definition in automatic mode.
    helper_auto_add_ = function(by, reference, prefix, joint, symmetric) {
      # validate features
      assert_valid_features(by)

      # load reference and validate reference names
      if (length(reference) > 0) {
        reference <- private$helper_auto_reference_(by, reference)
      }

      # store definition
      self$defs <- list(by = by, reference = reference, prefix = prefix,
                        joint = joint, symmetric = symmetric, mode = "auto")
      self$logic_tables <- list()
      self$param_overrides <- list()

      invisible(self)
    },

    # Add risk definition in custom mode.
    helper_custom_add_ = function(strata, params, symmetric, reference) {
      assert_valid_logic_table(strata)
      none_ref <- length(reference) > 0 && isTRUE(all(reference == ".none"))
      if (!none_ref && !".reference" %in% strata$param_name) stop_h(
        "Custom logic table must contain a {.val .reference} row.",
        paste("Designate one row's param_name as {.val .reference},",
              "or pass {.code reference = \".none\"} to opt out.")
      )
      # cast feature columns to character (consistent with auto mode)
      feat_cols <- setdiff(colnames(strata), risk_reserved_cols)
      strata[feat_cols] <- lapply(strata[feat_cols], as.character)
      # store sentinel definition
      self$defs <- list(by = NA_character_, prefix = NA_character_,
                        reference = NULL, joint = TRUE,
                        symmetric = symmetric, mode = "custom")
      self$logic_tables <- list(strata)
      self$param_overrides <- if (is.null(params)) list() else params

      invisible(self)
    },

    # Validate custom mode logic table against data, then add counts.
    helper_implement_custom_ = function(data_wide, symmetric) {
      for (i in seq_along(self$logic_tables)) {
        # optionally fill symmetric gaps
        this_table <- self$logic_tables[[i]]
        # retrieve covariates
        this_feat <- sub("^from_", "",
                         grep("^from_", names(this_table), value = TRUE))
        # validate that feature columns exist in data
        assert_data_features(this_feat, colnames(data_wide))

        if (symmetric) this_table <- fill_symmetric_gaps(this_table)
        # FIXME: this is a manual fix, I can't find the source of extra 'n'
        this_table$n <- NULL
        # FIXME: this is not counting the wildcars...
        this_table <- merge(
          this_table,
          contact_risk_logic_table(data_wide, this_feat),
          by = setdiff(names(this_table), risk_reserved_cols),
          all.x = TRUE
        )
        this_table$n[is.na(this_table$n)] <- 0
        self$logic_tables[[i]] <- this_table
      }
      invisible(self)
    }
  ),

  ############################################################################
  ################################## ACTIVE ##################################
  ############################################################################

  active = list(
    #' @field parameters `list`.
    #' Parameter definitions associated to contact risks.
    parameters = function(value) {
      new_params <- list()
      for (ptbl in self$logic_tables) {
        for (par_id in setdiff(ptbl$param_name, ".reference")) {
          if (!is.null(self$param_overrides[[par_id]])) {
            new_params[[par_id]] <- self$param_overrides[[par_id]]
          } else {
            new_params[[par_id]] <- self$raw_params
          }
        }
      }
      new_params
    },

    #' @field data_required `character()`.
    #' Data required by contact risks.
    data_required = function(value) {
      required <- self$defs[["by"]]
      if (is.null(required) || all(is.na(required))) return(character(0))
      required <- required[!is.na(required)]
      set_source_name_(required, self$family)
    },

    #' @field feature_names `character()`. Data columns the risks stratify on,
    #' declared through `by` or carried by a custom logic table.
    feature_names = function(value) {
      custom <- lapply(self$logic_tables, function(tbl) {
        sub("^from_", "", grep("^from_", names(tbl), value = TRUE))
      })
      feats <- as.character(unlist(c(self$defs[["by"]], custom)))
      unique(feats[!is.na(feats)])
    },

    #' @field param_required `character()`.
    #' Parameters required by contact risks.
    param_required = function(value) {
      required <- c()
      for (ptbl in self$logic_tables) required <- c(required, ptbl$param_name)
      required <- setdiff(unique(required), ".reference")
      if (length(required) == 0) return(character(0))
      set_source_name_(required, self$family)
    },

    #' @field tidy_table `data.frame()`.
    #' Contact risk definitions table.
    tidy_table = function(value) {
      tidyup_contact_risk_tables(self$logic_tables)
    }
  ),

  ############################################################################
  ################################## PUBLIC ##################################
  ############################################################################

  public = list(
    #' @field family `character(1)`. Holder family.
    family = "contact risk",

    #' @field defs `list`. Current contact risk definition.
    defs = list(),

    #' @field raw_params `list`. Raw parameter definitions.
    raw_params = list(),

    #' @field logic_tables `list`.
    #' List of tables for contact risk pair definitions.
    logic_tables = list(),

    #' @field param_overrides `list`.
    #' Per-parameter MCMC overrides.
    param_overrides = list(),

    #' @description
    #' New object of class [HousepiContactRisk][housepi::HousepiContactRisk].
    initialize = function() {
      self$defs <- list()
      self$raw_params <- list()
      self$logic_tables <- list()
      self$param_overrides <- list()
      invisible(self)
    },

    #' @inherit add_contact_risk
    add_quo = function(by = NA_character_, ...,
                       reference = NULL,
                       joint = FALSE, symmetric = TRUE,
                       strata = NULL, params = NULL,
                       prior = NA_character_,
                       min = NA_real_, max = NA_real_,
                       proposal = NA_character_,
                       prefix = NA_character_,
                       .warn_over = TRUE) {
      # automatic definition
      by <- housepi_smart_unwrap(by, "by")
      prefix <- housepi_smart_unwrap(prefix, "prefix")
      reference <- housepi_smart_unwrap(reference, "reference")
      joint <- isTRUE(housepi_smart_unwrap(joint, "joint"))
      symmetric <- isTRUE(housepi_smart_unwrap(symmetric, "symmetric"))

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

      # distinguish between automatic and custom definitions
      is_auto <- !all(is.na(by))
      is_custom <- !is.null(strata)
      if (is_auto && is_custom) stop_ambiguos_definition()

      # warn overwrite
      if (.warn_over && length(self$defs) > 0) warn_h(
        "Overwriting previous contact risk definition."
      )
      if (is_auto) {
        private$helper_auto_add_(by, reference, prefix, joint, symmetric)
      } else if (is_custom) {
        private$helper_custom_add_(strata, params, symmetric, reference)
      } else {
        stop_empty_definition()
      }

      # store MCMC parameter setup
      self$raw_params <- set_group_risks_mcmc_setup(
        prior = prior, args = args,
        min = min, max = max, proposal = proposal
      )

      invisible(self)
    },

    ############################################################################
    ############################################################################
    ############################################################################

    #' @description
    #' Implement contact risk definition according to data.
    #' @param data_wide A data.frame with the dataset.
    #' @param factor_info A list with factor levels.
    implement = function(data_wide, factor_info) {
      symmetric <- isTRUE(self$defs[["symmetric"]])

      # custom mode
      if (self$defs[["mode"]] == "custom") {
        private$helper_implement_custom_(data_wide, symmetric)
        return(invisible(self))
      }

      self$logic_tables <- list()

      by <- self$defs[["by"]]
      if (is.null(by) || all(is.na(by))) return(invisible(self))
      assert_data_features(by, colnames(data_wide))

      # cast features: joint wraps all together
      features <- as.list(by)
      if (isTRUE(self$defs[["joint"]])) features <- list(features)

      # default prefix
      prefix <- as.character(self$defs[["prefix"]])
      if (is.na(prefix)) prefix <- ""

      new_tables <- list()
      for (this_feat in features) {
        this_table <- contact_risk_logic_table(data_wide, this_feat)
        # extract reference entries for this sub-table, propagating .none
        ref <- self$defs$reference
        if (length(ref) > 0 && isTRUE(all(ref == ".none"))) {
          relevant_ref <- ".none"
        } else {
          ref_cols <- c(paste0("from_", this_feat), paste0("to_", this_feat))
          relevant_ref <- self$defs$reference[
            names(self$defs$reference) %in% ref_cols
          ]
        }
        this_table <- associate_contact_param_names(
          this_table, prefix, relevant_ref, symmetric
        )
        new_tables[[length(new_tables) + 1]] <- this_table
      }
      self$logic_tables <- new_tables
      invisible(self)
    },

    #' @description
    #' Assert that contact risk definitions are compatible.
    assert = function() {
      length(self$defs) > 0 &&
        (!all(is.na(self$defs[["by"]])) || self$defs[["mode"]] == "custom")
    },

    #' @description
    #' Tidy table of contact risk definitions for Cpp.
    #' @param info list of lists with factor levels.
    expose_to_cpp = function(info) {
      defs <- self$tidy_table
      if (nrow(defs) == 0) return(defs)

      # encode from_X and to_X columns using factor info
      features <- self$defs[["by"]]
      # for custom mode (by is NA), infer features from column names
      if (is.null(features) || all(is.na(features))) {
        from_cols <- grep("^from_", colnames(defs), value = TRUE)
        features <- sub("^from_", "", from_cols)
      }
      if (is.null(features) || length(features) == 0) return(defs)

      for (this_feat in features[!is.na(features)]) {
        if (is.null(info[[this_feat]])) next
        lv <- info[[this_feat]]$levels
        for (prefix in c("from_", "to_")) {
          col <- paste0(prefix, this_feat)
          if (col %in% colnames(defs)) {
            defs[[col]] <- as.integer(factor(defs[[col]], levels = lv))
          }
        }
      }
      defs
    }
  )
)
