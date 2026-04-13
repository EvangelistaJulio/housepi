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
  # unique values per feature
  uni_values <- lapply(features, function(f) sort(unique(data_wide[[f]])))
  names(uni_values) <- features

  # build from and to grids
  from_grid <- do.call(expand.grid, stats::setNames(
    uni_values, paste0("from_", features)
  ))
  to_grid <- do.call(expand.grid, stats::setNames(
    uni_values, paste0("to_", features)
  ))

  # full cross of from x to
  y <- merge(from_grid, to_grid, by = NULL)

  # convert factors to character
  for (col in colnames(y)) y[[col]] <- as.character(y[[col]])

  y
}

#' Associate parameter names to contact rate logic table.
#' @param this_table `data.frame` with from_X and to_X columns.
#' @param prefix Prefix for parameter names.
#' @param reference Named vector with from_/to_ prefixed names and values.
#' @return `data.frame` with param_name column added.
#' @keywords internal
#' @noRd
associate_contact_param_names <- function(this_table, prefix, reference) {
  colmask <- assess_feature_values(this_table)

  # find reference index
  if (length(reference) == 0) {
    idx_refer <- integer(0)
  } else {
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
  if (length(idx_refer) == 1) {
    this_table$param_name[idx_refer] <- ".reference"
  }

  this_table
}

############################################################################
############################################################################
########################### Housepi Contact Rate ###########################
############################################################################
############################################################################

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
    helper_auto_add_ = function(by, reference, prefix) {
      # validate features
      assert_valid_features(by)

      # validate reference names
      if (!is.null(reference)) {
        valid_names <- c(paste0("from_", sort(by)), paste0("to_", sort(by)))
        # convenience: allow un-prefixed names (apply to both from_ and to_)
        if (length(reference) == length(by) &&
              all(names(reference) %in% by)) {
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
                        mode = "auto")
      self$logic_table <- NULL
      self$param_overrides <- list()
      invisible(self)
    },

    # Add rate definition in custom mode.
    helper_custom_add_ = function(logic_table, params) {
      assert_valid_logic_table(logic_table)
      # store sentinel definition
      self$defs <- list(by = NA_character_, prefix = NA_character_,
                        reference = NULL, mode = "custom")
      self$logic_table <- logic_table
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
      if (!is.null(self$logic_table)) {
        for (par_id in setdiff(self$logic_table$param_name, ".reference")) {
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
      if (is.null(self$logic_table)) return(character(0))
      required <- setdiff(unique(self$logic_table$param_name), ".reference")
      set_source_name_(required, self$family)
    },

    #' @field tidy_table `data.frame()`.
    #' Contact rate definitions table.
    tidy_table = function(value) {
      if (!missing(value)) stop_read_only("tidy_table")
      if (is.null(self$logic_table)) return(data.frame())
      self$logic_table
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

    #' @field logic_table `data.frame`.
    #' Table of parameters associated to contact rate pairs.
    logic_table = NULL,

    #' @field param_overrides `list`.
    #' Per-parameter MCMC overrides.
    param_overrides = list(),

    #' @description
    #' New object of class [HousepiContactRate][housepi::HousepiContactRate].
    initialize = function() {
      self$defs <- list()
      self$raw_params <- list()
      self$logic_table <- NULL
      self$param_overrides <- list()
      invisible(self)
    },

    #' @inherit add_contact_rate
    add_quo = function(by = NA_character_,
                       reference = NULL, prefix = "cr",
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
      if (is_auto) private$helper_auto_add_(by, reference, prefix)
      else if (is_custom) private$helper_custom_add_(logic_table, params)
      else return(invisible(self))

      # store MCMC parameter setup
      self$raw_params <- set_group_rates_mcmc_setup(
        prior = prior, args = args, min = min, max = max,
        proposal = proposal, proposal_scale = proposal_scale
      )

      invisible(self)
    },

    #' @description
    #' Implement contact rate definition according to data.
    #' @param data_wide A data.frame with the dataset.
    #' @param factor_info A list with levels for each factorized feature.
    implement = function(data_wide, factor_info) {
      # skip if logic_table mode -- already stored directly
      if (self$defs[["mode"]] == "custom") return(invisible(self))

      self$logic_table <- NULL

      by <- self$defs[["by"]]
      if (is.null(by) || all(is.na(by))) return(invisible(self))
      assert_data_features(by, colnames(data_wide))

      # build the paired logic table
      this_table <- contact_rate_logic_table(data_wide, by)

      # default prefix
      prefix <- as.character(self$defs[["prefix"]])
      if (is.na(prefix)) prefix <- "cr"

      # associate parameter names
      this_table <- associate_contact_param_names(
        this_table, prefix, self$defs$reference
      )

      self$logic_table <- this_table
      invisible(self)
    },

    #' @description
    #' Assert that contact rate definitions are compatible with model.
    #' @param mcmc_ptr Rcpp pointer to MCMC object.
    assert = function(mcmc_ptr) {
      has_defs <- length(self$defs) > 0 &&
        (!all(is.na(self$defs[["by"]])) || self$defs[["mode"]] == "custom")
      has_model <- .engine$has_contact_rates(mcmc_ptr)
      if (has_defs && !has_model) stop_h(
        "Contact rates are defined but the loaded model does not use them.",
        "Check 'has_contact_rates()' in your model implementation."
      )
    },

    #' @description
    #' Tidy table of contact rate definitions to be passed to Cpp.
    #' @param info list of lists with factor levels for data columns.
    expose_to_cpp = function(info) {
      defs <- self$tidy_table
      if (nrow(defs) == 0) return(defs)

      # integer-encode feature columns using factor info
      for (k in intersect(colnames(defs), names(info))) {
        defs[[k]] <- as.integer(factor(defs[[k]], levels = info[[k]]$levels))
      }
      # also encode from_X and to_X using the base column info
      for (by_col in self$defs[["by"]]) {
        if (!is.null(info[[by_col]])) {
          lv <- info[[by_col]]$levels
          for (prefix in c("from_", "to_")) {
            col <- paste0(prefix, by_col)
            if (col %in% colnames(defs)) {
              defs[[col]] <- as.integer(factor(defs[[col]], levels = lv))
            }
          }
        }
      }
      defs
    }
  )
)
