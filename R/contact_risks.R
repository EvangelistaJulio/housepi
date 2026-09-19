## housepi - Epidemiological Analysis of Household Outbreaks.
##
## Copyright © 2024-2026 Institut Pasteur (Paris), CNRS and INSERM.
## Distributed under the terms of the GNU General Public License (GPLv3).
## See the LICENSE file for details.
##

#' Return logic table of contact risks (from/to paired columns).
#' @param data_wide Data from which to extract contact risks.
#' @param covariates Covariates as list of data column names.
#' @return data.frame with from_X, to_X columns and count n.
#' @keywords internal
#' @noRd
contact_risk_logic_table <- function(data_wide, covariates) {
  # convention: no house_id implies households of size 1: no pairs possible
  if (!"house_id" %in% colnames(data_wide)) return(NULL)

  # subset to needed column, and coerce covariates to character
  covariates <- as.character(unlist(covariates))
  data_wide <- data_wide[, c("house_id", covariates), drop = FALSE]
  data_wide[] <- lapply(data_wide, as.character)

  # within a household every 'from' meets every 'to', so the ordered pairs of
  # two covariate keys are the product of their per-household counts. Only the
  # members meeting themselves have to come back off the diagonal
  per_house <- table(data_wide$house_id, row_keys_(data_wide, covariates))
  cnt <- crossprod(per_house)
  diag(cnt) <- diag(cnt) - colSums(per_house)

  # build the full observed grid
  observed <- unique(data_wide[, covariates, drop = FALSE])
  from_grid <- stats::setNames(observed, paste0("from_", covariates))
  to_grid <- stats::setNames(observed, paste0("to_", covariates))
  y <- merge(from_grid, to_grid, by = NULL)
  rownames(y) <- NULL

  # look up counts by the pair of covariate keys
  y$n <- as.integer(cnt[cbind(row_keys_(y, names(from_grid)),
                              row_keys_(y, names(to_grid)))])
  y$n[is.na(y$n)] <- 0L
  y[y$n > 0, ]
}

#' Count the observed pairs every contact rule matches.
#' @param tbl `data.frame` logic table with from_X and to_X columns.
#' @param data_wide Data the pairs are counted over.
#' @param covariates `character()` Covariates the rules key on.
#' @return `integer()` One count per row of `tbl`.
#' @keywords internal
#' @noRd
count_contact_rule_pairs <- function(tbl, data_wide, covariates) {
  pairs <- contact_risk_logic_table(data_wide, covariates)
  # no household structure, no pair to count
  if (is.null(pairs) || nrow(pairs) == 0) return(rep(0L, nrow(tbl)))

  cols <- c(paste0("from_", covariates), paste0("to_", covariates))
  count_resolved_rules(tbl, pairs[cols], pairs$n)
}

#' Associate parameter names to contact risk logic table.
#' @param this_table `data.frame` with from_X and to_X columns.
#' @param prefix Prefix for parameter names.
#' @param reference Named vector with from_/to_ prefixed names and values.
#' `".auto"` picks the pair with the largest sample size.
#' Empty leaves the risk without a reference.
#' @return `data.frame` with param_name column added.
#' @keywords internal
#' @noRd
associate_contact_param_names <- function(this_table, prefix, reference,
                                          symmetric = TRUE) {
  tp <- "^to_"
  fp <- "^from_"
  colmask <- assess_covariate_values(this_table)
  from_cols <- grep(fp, colmask, value = TRUE)
  to_cols <- grep(tp, colmask, value = TRUE)
  covariates <- sub(fp, "", from_cols)

  # find reference index
  idx_refer <- NULL
  if (identical(reference, ".auto")) {
    idx_refer <- which.max(this_table$n)

  } else if (length(reference) > 0) {
    # find definition row
    refer_mask <- apply(this_table[colmask], 1, function(row) {
      all(row == unlist(reference[colmask]))
    })
    idx_refer <- unname(which(refer_mask))
    if (length(idx_refer) != 1) stop_h(
      "Provided reference contact risk is not valid.",
      `!` = sprintf(
        "The pair %s %s %s is never observed.",
        reference_to_phrase(reference[grep(fp, names(reference))], fp),
        ifelse(symmetric, "<->", "-->"),
        reference_to_phrase(reference[grep(tp, names(reference))], tp)
      )
    )
  }

  prefix <- format_risk_prefix_(prefix)

  # associate param names with one "from" and one "to" token
  this_table$param_name <- apply(
    this_table[colmask], 1, function(row) {
      from_part <- paste(covariates, row[from_cols], sep = "_", collapse = "_")
      to_part <- paste(covariates, row[to_cols], sep = "_", collapse = "_")
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
  cov_cols <- setdiff(colnames(tbl), risk_reserved_cols)
  from_cols <- grep("^from_", cov_cols, value = TRUE)
  covariates <- sub("^from_", "", from_cols)

  new_rows <- list()
  for (i in seq_len(nrow(tbl))) {
    rev_row <- tbl[i, , drop = FALSE]
    for (f in covariates) {
      rev_row[[paste0("from_", f)]] <- tbl[i, paste0("to_", f)]
      rev_row[[paste0("to_", f)]] <- tbl[i, paste0("from_", f)]
    }

    def_exists <- any(apply(tbl[cov_cols], 1, function(row) {
      rv <- unlist(rev_row[cov_cols])
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

#' Combine list of contact risk sub-tables into a single data.frame.
#' @param tables `list` of `data.frame`s.
#' @return `data.frame` with columns unified across sub-tables.
#' @keywords internal
#' @noRd
tidyup_contact_risk_tables <- function(tables) {
  if (length(tables) == 0) return(data.frame())

  present <- unique(unlist(lapply(tables, names)))
  all_cols <- unique(c("param_name", "n",
                       intersect(risk_pin_cols, present), sort(present)))

  combined <- do.call(rbind, lapply(tables, pad_risk_cols_, cols = all_cols))
  rownames(combined) <- NULL
  combined <- compute_rule_groups(combined)
  combined <- reorder_risks_logic_table(combined)
  rownames(combined) <- NULL
  combined
}

# Validate and shape the automatic reference of contact risks.
auto_contact_reference_ <- function(by, reference) {
  if (identical(reference, ".auto")) return(".auto")

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
}

############################################################################
############################################################################
########################### Housepi Contact Risk ###########################
############################################################################
############################################################################

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

    # Add risk definition in automatic mode.
    helper_auto_add_ = function(by, reference, prefix, joint, symmetric) {
      # validate covariates
      assert_valid_covariates(by)

      # load reference and validate reference names
      if (length(reference) > 0) {
        reference <- auto_contact_reference_(by, reference)
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
      if (length(reference) > 0) stop_h(
        paste("Argument {.arg reference} is not used",
              "with a custom {.arg strata} table."),
        "Tag one row's {.field param_name} as {.val .reference} instead."
      )
      assert_one_reference_per_group(strata, symmetric)
      if (length(params) > 0) assert_no_pin_params_clash(strata, params)
      # cast covariate columns to character (consistent with auto mode)
      cov_cols <- setdiff(colnames(strata), risk_reserved_cols)
      strata[cov_cols] <- lapply(strata[cov_cols], as.character)
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
        this_cov <- sub("^from_", "",
                        grep("^from_", names(this_table), value = TRUE))
        # validate that covariate columns exist in data
        assert_data_covariates(this_cov, colnames(data_wide))

        # rules naming unobserved levels are dropped before symmetrization
        # so a dead rule leaves no counterpart behind
        this_table <- drop_unknown_strata_rows(
          this_table,
          known = lapply(stats::setNames(nm = this_cov),
                         function(cov) unique(data_wide[[cov]])),
          label = self$family,
          # both sides of a pair read the same data column
          cols = stats::setNames(rep(this_cov, 2L),
                                 c(paste0("from_", this_cov),
                                   paste0("to_", this_cov)))
        )

        if (symmetric) this_table <- fill_symmetric_gaps(this_table)
        this_table$n <- count_contact_rule_pairs(this_table, data_wide,
                                                 this_cov)
        self$logic_tables[[i]] <- this_table
      }
      # a table whose every rule was dropped stratifies nothing:
      # it leaves definition behind, so every pair reads the reference
      self$logic_tables <- Filter(function(tbl) nrow(tbl) > 0,
                                  self$logic_tables)
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
        par <- risk_table_params_(ptbl, self$param_overrides, self$raw_params)
        new_params[names(par)] <- par
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

    #' @field covariate_names `character()`. Data columns the risks stratify on,
    #' declared through `by` or carried by a custom logic table.
    covariate_names = function(value) {
      custom <- lapply(self$logic_tables, function(tbl) {
        sub("^from_", "", grep("^from_", names(tbl), value = TRUE))
      })
      covs <- as.character(unlist(c(self$defs[["by"]], custom)))
      unique(covs[!is.na(covs)])
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

    #' @description
    #' Whether the risk carries absolute values, having no reference.
    #' @return `logical(1)`
    is_absolute = function() {
      if (length(self$defs) == 0) return(FALSE)
      if (identical(self$defs$mode, "auto"))
        return(is.null(self$defs$reference))
      !any(vapply(self$logic_tables,
                  function(tbl) ".reference" %in% tbl$param_name, logical(1)))
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
      reference <- housepi_smart_unwrap(reference, "reference", null_ok = TRUE)
      joint <- isTRUE(housepi_smart_unwrap(joint, "joint"))
      symmetric <- isTRUE(housepi_smart_unwrap(symmetric, "symmetric"))

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
      assert_no_sibling_dist_arg(args, self$family)

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
        name = self$family, prior = prior, args = args,
        min = min, max = max, proposal = proposal,
        absolute = self$is_absolute()
      )

      invisible(self)
    },

    ############################################################################
    ############################################################################
    ############################################################################

    #' @description
    #' Implement contact risk definition according to data.
    #' @param data_wide A data.frame with the dataset.
    implement = function(data_wide) {
      symmetric <- isTRUE(self$defs[["symmetric"]])

      # custom mode
      if (self$defs[["mode"]] == "custom") {
        private$helper_implement_custom_(data_wide, symmetric)
        return(invisible(self))
      }

      self$logic_tables <- list()

      by <- self$defs[["by"]]
      if (is.null(by) || all(is.na(by))) return(invisible(self))
      assert_data_covariates(by, colnames(data_wide))

      # cast covariates: joint wraps all together
      covariates <- as.list(by)
      if (isTRUE(self$defs[["joint"]])) covariates <- list(covariates)

      # default prefix
      prefix <- as.character(self$defs[["prefix"]])
      if (is.na(prefix)) prefix <- ""

      new_tables <- list()
      for (this_cov in covariates) {
        this_table <- contact_risk_logic_table(data_wide, this_cov)
        # extract reference entries for this sub-table, propagating .auto
        relevant_ref <- self$defs$reference
        if (!identical(relevant_ref, ".auto")) {
          ref_cols <- c(paste0("from_", this_cov), paste0("to_", this_cov))
          relevant_ref <- relevant_ref[names(relevant_ref) %in% ref_cols]
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
      defs[risk_pin_cols] <- NULL

      # encode from_X and to_X columns using factor info
      covariates <- self$defs[["by"]]
      # for custom mode (by is NA), infer covariates from column names
      if (is.null(covariates) || all(is.na(covariates))) {
        from_cols <- grep("^from_", colnames(defs), value = TRUE)
        covariates <- sub("^from_", "", from_cols)
      }
      if (is.null(covariates) || length(covariates) == 0) return(defs)

      cols <- paste0(rep(c("from_", "to_"), each = length(covariates)),
                     covariates)
      encode_risk_levels_(defs, info,
                          key = function(col) sub("^(from|to)_", "", col),
                          cols = intersect(cols, colnames(defs)))
    }
  )
)
