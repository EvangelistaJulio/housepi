## housepi - Epidemiological Analysis of Household Outbreaks.
##
## Copyright © 2026 Institut Pasteur (Paris), CNRS and INSERM.
## Distributed under the terms of the GNU General Public License (GPLv3).
## See the LICENSE file for details.
##

#' Return logic table of contact rates (from/to paired columns).
#' @param data_wide Data from which to extract contact rates.
#' @param features Features as list of data column names.
#' @return data.frame with from_X, to_X columns and count n.
#' @keywords internal
#' @noRd
contact_rate_logic_table <- function(data_wide, features) {
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

  y
}

#' Associate parameter names to contact rate logic table.
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
  if (length(reference) == 0) {
    idx_refer <- which.max(this_table$n)
  } else {
    refer_mask <- apply(this_table[colmask], 1, function(row) {
      all(row == unlist(reference[colmask]))
    })
    idx_refer <- unname(which(refer_mask))
    if (length(idx_refer) != 1) stop_h(
      "Provided reference group is not valid.",
      "Reference should be a vector with names starting with 'from_X', 'to_Y'."
    )
  }

  # format prefix
  prefix <- ifelse(prefix == "" || !is.character(prefix),
                   "", paste0(prefix, "_"))

  # build param names with one "from" and one "to" token
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
      paste0(prefix, "from_", from_part, "_to_", to_part)
    }
  )
  if (length(idx_refer) == 1) {
    this_table$param_name[idx_refer] <- ".reference"
  }

  this_table
}

#' Fill symmetric gaps in a custom contact rate logic table.
#' For each row, if the reverse from/to pair is missing,
#' add it with the same param_name.
#' @param tbl `data.frame` with from_X, to_Y, and param_name columns.
#' @return `data.frame` with missing reverse pairs added.
#' @keywords internal
#' @noRd
fill_symmetric_gaps <- function(tbl) {
  feat_cols <- setdiff(colnames(tbl), c("param_name", "n"))
  from_cols <- grep("^from_", feat_cols, value = TRUE)
  features <- sub("^from_", "", from_cols)

  new_rows <- list()
  for (i in seq_len(nrow(tbl))) {
    rev_row <- tbl[i, , drop = FALSE]
    for (f in features) {
      rev_row[[paste0("from_", f)]] <- tbl[i, paste0("to_", f)]
      rev_row[[paste0("to_", f)]] <- tbl[i, paste0("from_", f)]
    }
    exists <- any(apply(tbl[feat_cols], 1, function(row) {
      all(row == rev_row[feat_cols])
    }))
    if (!exists) new_rows[[length(new_rows) + 1]] <- rev_row
  }

  if (length(new_rows) > 0) {
    tbl <- rbind(tbl, do.call(rbind, new_rows))
    rownames(tbl) <- NULL
  }
  tbl
}

############################################################################
############################################################################
########################### Housepi Contact Rate ###########################
############################################################################
############################################################################

#' Combine list of contact rate sub-tables into a single data.frame.
#' @param tables `list` of `data.frame`s.
#' @return `data.frame` with columns unified across sub-tables.
#' @keywords internal
#' @noRd
tidyup_contact_rate_tables <- function(tables) {
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
  combined
}

#' @title Housepi Contact Rate
#' @description Help with contact rate definitions of model Rcpp object.
#' @keywords internal
#' @noRd
HousepiContactRate <- R6::R6Class(
  classname = "HousepiContactRate",
  cloneable = TRUE,

  ############################################################################
  ################################## PRIVATE #################################
  ############################################################################

  private = list(
    deep_clone = function(name, value) helper_clone(value),

    # Add rate definition in automatic mode.
    helper_auto_add_ = function(by, reference, prefix, joint, symmetric) {
      # validate features
      assert_valid_features(by)

      # validate reference names
      if (!is.null(reference)) {
        valid_names <- c(paste0("from_", sort(by)), paste0("to_", sort(by)))
        # convenience: allow un-prefixed names
        if (length(reference) == length(by) && all(names(reference) %in% by)) {
          reference <- c(stats::setNames(reference[by], paste0("from_", by)),
                         stats::setNames(reference[by], paste0("to_", by)))
        }
        if (length(setdiff(names(reference), valid_names)) != 0) stop_h(
          "Reference does not use 'from_' and 'to_' prefixes matching 'by'."
        )
        reference <- as.list(reference)
      }

      # store definition
      self$defs <- list(by = by, reference = reference, prefix = prefix,
                        joint = joint, symmetric = symmetric, mode = "auto")
      self$logic_tables <- list()
      self$param_overrides <- list()
      invisible(self)
    },

    # Add rate definition in custom mode.
    helper_custom_add_ = function(logic_table, params, symmetric) {
      assert_valid_logic_table(logic_table)
      if (!".reference" %in% logic_table$param_name) stop_h(
        "Custom logic table must contain a '.reference' row.",
        "Designate one row's param_name as '.reference'."
      )
      # store sentinel definition
      self$defs <- list(by = NA_character_, prefix = NA_character_,
                        reference = NULL, joint = TRUE,
                        symmetric = symmetric, mode = "custom")
      self$logic_tables <- list(logic_table)
      self$param_overrides <- if (is.null(params)) list() else params
      invisible(self)
    }
  ),

  ############################################################################
  ################################## ACTIVE ##################################
  ############################################################################

  active = list(
    #' @field parameters `list`.
    #' Parameter definitions associated to contact rates.
    parameters = function(value) {
      if (!missing(value)) stop_read_only("parameters")
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
    #' Data required by contact rates.
    data_required = function(value) {
      if (!missing(value)) stop_read_only("data_required")
      required <- self$defs[["by"]]
      if (is.null(required) || all(is.na(required))) return(character(0))
      required <- required[!is.na(required)]
      set_source_name_(required, self$family)
    },

    #' @field param_required `character()`.
    #' Parameters required by contact rates.
    param_required = function(value) {
      if (!missing(value)) stop_read_only("param_required")
      required <- c()
      for (ptbl in self$logic_tables) required <- c(required, ptbl$param_name)
      required <- setdiff(unique(required), ".reference")
      if (length(required) == 0) return(character(0))
      set_source_name_(required, self$family)
    },

    #' @field tidy_table `data.frame()`.
    #' Contact rate definitions table.
    tidy_table = function(value) {
      if (!missing(value)) stop_read_only("tidy_table")
      tidyup_contact_rate_tables(self$logic_tables)
    }
  ),

  ############################################################################
  ################################## PUBLIC ##################################
  ############################################################################

  public = list(
    #' @field family `character(1)`. Module family.
    family = "contact rate",

    #' @field defs `list`. Current contact rate definition.
    defs = list(),

    #' @field raw_params `list`. Raw parameter definitions.
    raw_params = list(),

    #' @field logic_tables `list`.
    #' List of tables for contact rate pair definitions.
    logic_tables = list(),

    #' @field param_overrides `list`.
    #' Per-parameter MCMC overrides.
    param_overrides = list(),

    #' @description
    #' New object of class [HousepiContactRate][housepi::HousepiContactRate].
    initialize = function() {
      self$defs <- list()
      self$raw_params <- list()
      self$logic_tables <- list()
      self$param_overrides <- list()
      invisible(self)
    },

    #' @inherit add_contact_rate
    add_quo = function(by = NA_character_,
                       reference = NULL, prefix = "ctr",
                       joint = TRUE, symmetric = TRUE,
                       logic_table = NULL, params = NULL,
                       prior = NA_character_, ...,
                       min = NA_real_, max = NA_real_,
                       proposal = NA_character_,
                       proposal_scale = NA_real_,
                       .warn_over = TRUE) {
      # automatic definition
      by <- housepi_smart_unwrap(by, "by")
      prefix <- housepi_smart_unwrap(prefix, "prefix")
      reference <- housepi_smart_unwrap(reference, "reference")
      joint <- isTRUE(housepi_smart_unwrap(joint, "joint"))
      symmetric <- isTRUE(housepi_smart_unwrap(symmetric, "symmetric"))

      # custom definition
      logic_table <- housepi_smart_unwrap(logic_table, "logic_table")
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
      is_custom <- !is.null(logic_table)
      if (is_auto && is_custom) stop_ambiguos_definition()

      # warn overwrite
      if (.warn_over && length(self$defs) > 0) warning_h(
        "Overwriting previous contact rate definition."
      )
      if (is_auto) {
        private$helper_auto_add_(by, reference, prefix, joint, symmetric)
      } else if (is_custom) {
        private$helper_custom_add_(logic_table, params, symmetric)
      } else {
        return(invisible(self))
      }

      # store MCMC parameter setup
      self$raw_params <- set_group_rates_mcmc_setup(
        prior = prior, args = args, min = min, max = max,
        proposal = proposal, proposal_scale = proposal_scale
      )

      invisible(self)
    },

    ############################################################################
    ############################################################################
    ############################################################################

    #' @description
    #' Implement contact rate definition according to data.
    #' @param data_wide A data.frame with the dataset.
    #' @param factor_info A list with factor levels.
    implement = function(data_wide, factor_info) {
      symmetric <- isTRUE(self$defs[["symmetric"]])

      # custom mode
      if (self$defs[["mode"]] == "custom") {
        if (!symmetric) return(invisible(self))
        # optionally fill symmetric gaps
        for (i in seq_along(self$logic_tables)) {
          self$logic_tables[[i]] <- fill_symmetric_gaps(self$logic_tables[[i]])
        }
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
      if (is.na(prefix)) prefix <- "ctr"

      new_tables <- list()
      for (this_feat in features) {
        this_table <- contact_rate_logic_table(data_wide, this_feat)
        # extract reference entries for this sub-table
        ref_cols <- c(paste0("from_", this_feat), paste0("to_", this_feat))
        relevant_ref <- self$defs$reference[
          names(self$defs$reference) %in% ref_cols
        ]
        this_table <- associate_contact_param_names(
          this_table, prefix, relevant_ref, symmetric
        )
        new_tables[[length(new_tables) + 1]] <- this_table
      }
      self$logic_tables <- new_tables
      invisible(self)
    },

    #' @description
    #' Assert that contact rate definitions are compatible.
    #' @param mcmc_ptr Rcpp pointer to MCMC object.
    assert = function(mcmc_ptr) {
      has_defs <- length(self$defs) > 0 &&
        (!all(is.na(self$defs[["by"]])) || self$defs[["mode"]] == "custom")
      has_model <- .engine$has_contact_rates(mcmc_ptr)
      if (has_defs && !has_model) stop_h(
        "Contact rates are defined but the loaded model does not use them."
      )
    },

    #' @description
    #' Tidy table of contact rate definitions for Cpp.
    #' @param info list of lists with factor levels.
    expose_to_cpp = function(info) {
      defs <- self$tidy_table
      if (nrow(defs) == 0) return(defs)

      # encode from_X and to_X columns using factor info
      features <- self$defs[["by"]]
      if (is.null(features)) return(defs)

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
