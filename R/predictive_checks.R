## housepi - Epidemiological Analysis of Household Outbreaks.
##
## Copyright © 2024-2026 Institut Pasteur (Paris), CNRS and INSERM.
## Distributed under the terms of the GNU General Public License (GPLv3).
## See the LICENSE file for details.
##

rd_thin_ess_ <- "`logical(1)` Whether to thin `param_table` to its
bulk effective sample size (ESS).
Ignored if `param_table` has no posterior draws."
rd_max_nrep_ <- "`integer(1)` Larger number of replicates.
The parameter table (optionally ESS-thinned) is thinned to `max_nrep` rows."
rd_min_hhsize_ <- "`integer(1)` Minimum household size retained in the outputs.
If `NA`, no lower bound is applied."
rd_max_hhsize_ <- "`integer(1)` Maximum household size retained in the outputs.
If `NA`, no upper bound is applied."
rd_measure_dots_ <- "Named arguments. `function`s to compute measures from data.
Each function should return one between `numeric()`, `list` or `data.frame`."
rd_point_estimate_ <- "`character(1)` Summary statistic used as the predicted
point estimate across simulations.
One between `\"mean\"` and `\"median\"`."

############################################################################
############################################################################
############################################################################

# UPDATE: is this still necessary?
#' @importFrom rlang .data
NULL

helper_process_data <- function(data_wide) {
  df <- as.data.frame(data_wide)
  df$status <- as.integer(df$status)
  df$house_id <- as.character(df$house_id)
  df
}

#' @title Compute observed attack rate by household size
#' @description
#' Retrieve attack rates from household data.
#' @param data_wide `r rd_add_wide_`
#' @return `data.frame`
#' @details
#' The attack rate is computed for all households
#' with size \eqn{n \ge 1} according to
#' \eqn{\frac{\mathbb{E}_{n}[n^{+}]}{n}}
#' @examples
#' measure_ar_by_hhsize(flu_1918_frederick)
#' @export
measure_ar_by_hhsize <- function(data_wide) {
  .engine$measure_ar_by_hhsize(helper_process_data(data_wide))
}

#' @title Compute observed secondary attack rate by household size
#' @description
#' Retrieve secondary attack rates (SAR) from household data.
#' @param data_wide `r rd_add_wide_`
#' @return `data.frame`
#' @details
#' The secondary attack rate is computed for all households
#' with size \eqn{n > 1} according to
#' \eqn{\frac{\mathbb{E}_{n}[n^{+}] - 1}{n - 1}}
#' @examples
#' measure_sar_by_hhsize(flu_1918_frederick)
#' @export
measure_sar_by_hhsize <- function(data_wide) {
  .engine$measure_sar_by_hhsize(helper_process_data(data_wide))
}

#' @title Tabulate observed household cases by household size
#' @description
#' Retrieve cases by household size from household data.
#' @param data_wide `r rd_add_wide_`
#' @return `data.frame`
#' @examples
#' measure_cases_by_hhsize(flu_1918_frederick)
#' @export
measure_cases_by_hhsize <- function(data_wide) {
  .engine$measure_cases_by_hhsize(helper_process_data(data_wide))
}

############################################################################
############################################################################
######################### Housepi Posterior Prediction #####################
############################################################################
############################################################################

# Validate summary statistic choice.
ppc_set_point_estimate <- function(value) {
  match.arg(value, c("mean", "median"))
}

#' Apply function `func` to corresponding elements
#' across a list of data.frame that share the same dimensions.
#' @param frame_list `list` of `data.frame`s. Homogeneous dataframes.
#' @param func `function`. To be applied to corresponding data elements.
#' @return A `data.frame` of the same dimensions as the input dataframes.
#' Columns and rownames of the first data.frame in the list.
#' @param name `character(1)` Measure name.
#' @keywords internal
#' @noRd
ppc_frame_apply <- function(frame_list, func, name = NA_character_) {
  # ensure all elements share the same dimensions
  elem_n1 <- frame_list[[1]]
  dims_match <- vapply(frame_list[-1],
                       function(f) identical(dim(f), dim(elem_n1)),
                       logical(1))
  if (!all(dims_match)) stop_h(
    paste0("Not all computed measures share identical dimensions",
           ifelse(is.na(name), ".", " in measure {.field {name}}.")),
    "Check output of your measure function."
  )

  if (is.null(dim(elem_n1))) {
    # scalar or vector
    output <- sapply(seq_along(elem_n1),
                     function(x) func(sapply(frame_list, `[[`, x)))
    names(output) <- names(elem_n1)
    if (is.list(elem_n1)) output <- as.list(output)

  } else {
    # matrix or data.frame
    # NOTE: I need operations to be column-wise for eventual mixed types
    nr <- nrow(elem_n1)

    output <- as.data.frame(stats::setNames(
      lapply(seq_len(ncol(elem_n1)), function(x) {
        apply(matrix(unlist(lapply(frame_list, `[[`, x)), nrow = nr), 1L, func)
      }),
      colnames(elem_n1)
    ), check.names = FALSE)
    rownames(output) <- rownames(elem_n1)
    if (!is.data.frame(elem_n1)) output <- as.matrix(output)
  }

  output
}

#' Apply a numeric statistic across corresponding elements of dataframes.
#' Falls back to NA when all values are NA, and to the unique value when
#' the inputs are non-numeric but all identical.
#' @param frame_list `list` of `data.frame`s. Homogeneous dataframes.
#' @param stat_fn `function`. Numeric statistic applied to each element.
#' @param label `character(1)`. Statistic name used in error messages.
#' @param name `character(1)` Measure name.
#' @return `data.frame`
#' @keywords internal
#' @noRd
ppc_frame_apply_stat <- function(frame_list, stat_fn, label,
                                 name = NA_character_) {
  ppc_frame_apply(frame_list, function(x) {
    if (all(is.na(x))) NA_real_
    else if (all(is.numeric(x))) stat_fn(x)
    else if (length(unique(x)) == 1) unique(x)
    else stop_h(
      paste0("Undefined {.fn {label}} for non-numeric, non-unique values",
             ifelse(is.na(name), ".", " in measure {.field {name}}.")),
      "Check output of your measure function."
    )
  }, name = name)
}

#' Compute a summary statistic across different dataframes.
#' @param frame_list `list` of `data.frame`s.
#' @param point_estimate `character(1)`. Either `"mean"` or `"median"`.
#' @param name `character(1)` Measure name.
#' @return `data.frame`
#' @keywords internal
#' @noRd
ppc_frame_apply_summary <- function(frame_list,
                                    point_estimate = "median",
                                    name = NA_character_) {
  point_estimate <- ppc_set_point_estimate(point_estimate)
  fn <- switch(point_estimate,
               mean = function(x) mean(x, na.rm = TRUE),
               median = function(x) stats::median(x, na.rm = TRUE))
  ppc_frame_apply_stat(frame_list, fn, label = point_estimate, name = name)
}

#' Compute the credible interval across different dataframes.
#' @param frame_list `list` of `data.frame`s. Homogeneous dataframes.
#' @param cri_level `numeric(1)` Level for credible interval.
#' It must be a value in (0.5, 1.0).
#' @param upper `logical(1)` Compute the upper credible interval.
#' If `FALSE`, compute the lower.
#' @param name `character(1)` Measure name.
#' @return `data.frame`
#' @keywords internal
#' @noRd
ppc_frame_apply_cri <- function(frame_list, cri_level = 0.95,
                                upper = TRUE, name = NA_character_) {
  p <- cri_level_probs(cri_level)[[ifelse(isTRUE(upper), "upper", "lower")]]
  ppc_frame_apply_stat(frame_list,
                       function(x) stats::quantile(x, probs = p, na.rm = TRUE),
                       label = "credible interval",
                       name = name)
}

# Element-wise mask of observed values that fall outside [cri_lower, cri_upper].
# Non-numeric or NA elements are reported as FALSE (not a mismatch).
ppc_observed_outside_cri <- function(x) {
  .vec_outside <- function(obs, low, up) {
    if (!is.numeric(obs) || !is.numeric(low) || !is.numeric(up)) {
      return(rep(FALSE, length(obs)))
    }
    !is.na(obs) & !is.na(low) & !is.na(up) & (obs < low | obs > up)
  }
  if (!is.data.frame(x$observed)) {
    return(.vec_outside(x$observed, x$cri_lower, x$cri_upper))
  }
  as.data.frame(Map(.vec_outside, x$observed, x$cri_lower, x$cri_upper))
}

# extract predicted quantities from cases-by-house-size table list
ppc_extract_predicted_results <- function(x, point_estimate = "median", cri,
                                          name = NULL) {
  list(predicted = ppc_frame_apply_summary(x, point_estimate, name = name),
       cri_lower = ppc_frame_apply_cri(x, cri, upper = FALSE, name = name),
       cri_upper = ppc_frame_apply_cri(x, cri, upper = TRUE, name = name))
}

#' @title Posterior predictive check container for user-defined measures
#' @description
#' Class to extract posterior predictive checks over measures on simulations.
#' @docType class
HousepiCustomPPC <- R6::R6Class(
  classname = "HousepiCustomPPC",
  cloneable = TRUE,
  lock_objects = FALSE,

  ############################################################################
  ################################## PRIVATE #################################
  ############################################################################

  private = list(
    deep_clone = function(name, value) helper_clone(value),
    cri_level_ = 0.95,
    point_estimate_ = "mean",
    nrep_ = NA_integer_,
    thin_ess_ = TRUE,
    is_posterior_ = FALSE,
    model_name_ = NA_character_,
    measure_names_ = character(0),
    measures_ = list(),

    # Compare observed values against the credible interval of each measure.
    # Collapses all mismatches into a single warning: one finding line per
    # measure that has observed entries outside its [cri_lower, cri_upper]
    # bounds, followed by one closing hint.
    check_observed_in_cri_ = function() {
      cri_pct <- round(100 * private$cri_level_)
      findings <- c()
      for (m in private$measure_names_) {
        mismatch <- ppc_observed_outside_cri(self[[m]])
        len_bad <- sum(unlist(mismatch))
        if (len_bad > 0L) findings <- c(findings, sprintf(
          "%s observed %s outside of %d%% CrI for measure {.field %s}.",
          format_there_were_was(len_bad), as_plural("value", len_bad),
          cri_pct, m
        ))
      }
      if (length(findings) > 0L) {
        names(findings) <- rep("!", length(findings))
        warn_h(c(
          "Posterior predictive checks flagged the following:",
          findings,
          "Observed data may be poorly captured by the model."
        ))
      }
      invisible(self)
    }
  ),

  ############################################################################
  ################################## ACTIVE ##################################
  ############################################################################

  active = list(
    #' @field point_estimate `character(1)`. Summary statistic name.
    point_estimate = function(value) {
      if (!missing(value)) stop_read_only("point_estimate")
      private$point_estimate_
    },
    #' @field cri_level `numeric(1)`. Credible interval level.
    cri_level = function(value) {
      if (!missing(value)) stop_read_only("cri_level")
      private$cri_level_
    },
    #' @field nrep `integer(1)`. Number of replicates used for the predictions.
    nrep = function(value) {
      if (!missing(value)) stop_read_only("nrep")
      private$nrep_
    },
    #' @field thin_ess `logical(1)`. Whether thinning to ESS was applied.
    thin_ess = function(value) {
      if (!missing(value)) stop_read_only("thin_ess")
      private$thin_ess_
    },
    #' @field is_posterior `logical(1)`.
    #' Whether the predictions are posterior predictions.
    is_posterior = function(value) {
      if (!missing(value)) stop_read_only("is_posterior")
      private$is_posterior_
    },
    #' @field model_name `character(1)`. Name of the underlying model.
    model_name = function(value) {
      if (!missing(value)) stop_read_only("model_name")
      private$model_name_
    }
  ),

  ############################################################################
  ################################### PUBLIC #################################
  ############################################################################

  public = list(
    #' @description
    #' Given a model and a posterior sample of parameters,
    #' performs simulation and compute predicted values of a set of measures.
    #' @param hhmodel `r rd_model_in_`
    #' @param param_table `r rd_param_table_`
    #' @param obs_process `r rd_simul_observ_`
    #' @param obs_process_long `r rd_simul_obs_long_`
    #' @param ... `r rd_measure_dots_`
    #' @param thin_ess `r rd_thin_ess_`
    #' @param max_nrep `r rd_max_nrep_`
    #' @param measures_long `function` | `list` of `function`s.
    #' Applied to `data_long` (derived via `obs_process_long`).
    #' @param threads `r rd_thread_`
    #' @param progress `r rd_progress_`
    #' @param point_estimate `r rd_point_estimate_`
    initialize = function(hhmodel,
                          param_table,
                          obs_process = NULL,
                          obs_process_long = NULL,
                          ...,
                          measures_long = NULL,
                          thin_ess = TRUE,
                          max_nrep = 200L,
                          threads = 1L,
                          progress = TRUE,
                          point_estimate = "median") {
      if (!is.na(max_nrep) && (!is.numeric(max_nrep) || max_nrep < 1)) stop_h(
        "{.arg max_nrep} must be a positive integer or NA."
      )
      # associate summary statistic and credible interval
      private$point_estimate_ <- ppc_set_point_estimate(point_estimate)
      private$cri_level_ <- get_global_cri_level()

      # load and classify measure functions (wide vs long)
      m_wide <- list(...)
      has_wide <- length(m_wide) > 0
      # FIXME this breaks if it's already list
      m_long <- list()
      if (!is.null(measures_long)) m_long <- list(measures_long)
      has_long <- length(m_long) > 0

      if (!has_wide && !has_long) stop_h("No measure was provided.")
      # FIXME: temporary forcing obs process long
      if (has_long && is.null(obs_process_long)) stop_h(
        "{.arg obs_process_long} is required when long measures are provided."
      )

      # store (and test) results from observed data wide and/or long
      observations <- list()

      # wide-measure application
      apply_wide <- NULL
      if (has_wide) {
        observations <- c(
          observations, lapply(m_wide, function(func) func(hhmodel$data_wide))
        )
        apply_wide <- function(x) lapply(m_wide, function(func) func(x))
      }
      # long-measure application
      apply_long <- NULL
      if (has_long) {
        observations <- c(
          observations, lapply(m_long, function(func) func(hhmodel$data_long))
        )
        apply_long <- function(x) lapply(m_long, function(func) func(x))
      }

      # check measures
      if (any(lapply(observations, NROW) == 0)) stop_h(
        "One or more measures returned empty observation on data."
      )

      # run measures on posterior simulations
      output <- hhmodel$.__enclos_env__$private$run_measure_(
        param_table,
        apply_wide = apply_wide,
        apply_long = apply_long,
        obs_proc_wide = obs_process,
        obs_proc_long = obs_process_long,
        thin_ess = thin_ess,
        max_nrep = max_nrep,
        threads = threads,
        progress = progress
      )

      # parse results
      all_measure_names <- c(names(m_wide), names(m_long))

      parsed <- c()
      if (has_wide) parsed <- c(parsed,
        lapply(names(m_wide), function(x) lapply(output$wide, `[[`, x))
      )
      if (has_long) parsed <- c(parsed,
        lapply(names(m_long), function(x) lapply(output$long, `[[`, x))
      )
      names(parsed) <- all_measure_names

      predictions <- Map(
        ppc_extract_predicted_results,
        x = parsed, name = names(parsed),
        MoreArgs = list(point_estimate = private$point_estimate_,
                        cri = private$cri_level_)
      )

      # expose each measure as a read-only active binding
      for (mn in all_measure_names) {
        private$measures_[[mn]] <- list(
          observed = observations[[mn]],
          predicted = predictions[[mn]]$predicted,
          cri_lower = predictions[[mn]]$cri_lower,
          cri_upper = predictions[[mn]]$cri_upper
        )
        # make read-only active bindings
        local({
          m_name <- mn
          setter <- function(value) {
            if (!missing(value)) stop_read_only(m_name)
            private$measures_[[m_name]]
          }
          makeActiveBinding(m_name, setter, self)
        })
      }
      # store decorations for print and read-only access
      private$model_name_ <- hhmodel$name
      private$is_posterior_ <-
        inherits(param_table, "HousepiPosteriorDraws") ||
        posterior::is_draws(param_table)
      private$measure_names_ <- all_measure_names
      private$nrep_ <- length(parsed[[1]])
      private$thin_ess_ <- isTRUE(thin_ess)

      private$check_observed_in_cri_()

      invisible(self)
    },

    #' @description
    #' Override the default printing behavior.
    print = function() {
      kind <- if (private$is_posterior_) "Posterior predictions" else
        "Predictions"
      sms <- c(i = "{kind} for model {.val {private$model_name_}}:")
      for (m in private$measure_names_) {
        sms <- c(sms, paste0(cli::symbol$bullet, " ", m))
      }
      print_bullets(c(
        sms,
        paste0("N. replicates: {private$nrep_}",
               ifelse(private$thin_ess_, " (ESS)", "")),
        "Point estimate: {private$point_estimate_}",
        "Credible interval: {round(100 * private$cri_level_)}%"
      ))
      invisible(self)
    }
  )
)

############################################################################
############################################################################
############################ Predicted Cases By Size #######################
############################################################################
############################################################################

# @description
# Total households and households with at least one case, per household size.
#
# @return `data.frame`
helper_house_counts <- function(data_wide) {
  .engine$count_houses_by_hhsize(helper_process_data(data_wide))
}

# @description
# Return a list with all filtered case tables.
#
# @return `list`
helper_filter_tables <- function(cases_count,
                                 min_hhsize = NA,
                                 max_hhsize = NA) {
  # load observed table to use as a type shape template
  obsdf <- cases_count$observed

  # remove rows without counts
  mask_row <- rowSums(obsdf[colnames(obsdf) != "hhsize"], na.rm = TRUE) != 0

  # remove column with zero cases if none was observed
  mask_col <- colnames(obsdf) != "0" |
    colSums(obsdf == 0, na.rm = TRUE) != nrow(obsdf)

  # implement max household size mask (also caps the case columns, as a
  # household cannot hold more cases than its size)
  if (!is.na(max_hhsize)) {
    col_max <- which(colnames(obsdf) == max_hhsize)
    mask_col <- mask_col & (seq_len(ncol(obsdf)) <= col_max)
    mask_row <- mask_row & (obsdf$hhsize <= max_hhsize)
  }

  # implement min household size mask (rows only: any case column stays valid)
  if (!is.na(min_hhsize)) {
    mask_row <- mask_row & (obsdf$hhsize >= min_hhsize)
  }

  list(
    observed = obsdf[mask_row, mask_col],
    predicted = cases_count$predicted[mask_row, mask_col],
    cri_lower = cases_count$cri_lower[mask_row, mask_col],
    cri_upper = cases_count$cri_upper[mask_row, mask_col]
  )
}

# @description
# Function to combine cases table.
#
# @param obsdf A `data.frame` with observed cases.
# @param preddf A `data.frame` with predicted cases.
#
# @return `data.frame`
helper_combine_obs_pred <- function(obsdf, preddf) {
  # 1. curate observed dataframe
  obsdf$distr <- "observed"
  rownames(obsdf) <- seq(1, 2 * nrow(obsdf), 2)

  # 2. curate predicted dataframe
  preddf$hhsize <- obsdf$hhsize
  preddf$distr <- "predicted"
  rownames(preddf) <- seq(2, 2 * nrow(preddf), 2)

  # 3. combine by row bind
  combined <- rbind(obsdf, preddf)[as.character(1:(2 * nrow(obsdf))), ]
  # reset row names
  rownames(combined) <- NULL
  # reorder columns
  y <- colnames(combined)
  aux <- c("hhsize", "distr")
  combined <- combined[, c(aux, y[!(y %in% aux)])]

  # format columns to extract labels with labelizor
  colnames(combined)[3:ncol(combined)] <- sapply(
    colnames(combined)[3:ncol(combined)],
    function(x) paste("cases", x, sep = "_")
  )

  combined
}

#' @title Standard epidemiological posterior predictive check
#' @description
#' Posterior predictive check on the standard epidemiological measures by
#' household size: cases, attack rate and secondary attack rate, comparing
#' observed and predicted values.
#' @docType class
HousepiPPC <- R6::R6Class(
  classname = "HousepiPPC",
  inherit = HousepiCustomPPC,
  cloneable = TRUE,
  lock_objects = FALSE,

  ############################################################################
  ################################## PRIVATE #################################
  ############################################################################

  private = list(
    deep_clone = function(name, value) helper_clone(value),
    min_hhsize_ = NA_integer_,
    max_hhsize_ = NA_integer_,
    house_counts_ = NULL,

    # Restrict a by-hhsize data.frame to the set [min_hhsize_, max_hhsize_].
    filter_hhsize_ = function(df) {
      if (!is.na(private$min_hhsize_)) {
        df <- df[df$hhsize >= private$min_hhsize_, ]
      }
      if (!is.na(private$max_hhsize_)) {
        df <- df[df$hhsize <= private$max_hhsize_, ]
      }
      df
    }
  ),

  ############################################################################
  ################################## ACTIVE ##################################
  ############################################################################

  active = list(
    #' @field formatted_attack_rate `data.frame`.
    #' Observed and predicted attack rates by household size.
    #' Sizes outside the set household size bounds are not returned.
    formatted_attack_rate = function(value) {
      if (!missing(value)) stop_read_only("formatted_attack_rate")
      # convert internal list to data.frame
      output <- data.frame(
        hhsize = self$attack_rate$observed$hhsize,
        observed = self$attack_rate$observed$ar,
        predicted = self$attack_rate$predicted$ar,
        cri_lower = self$attack_rate$cri_lower$ar,
        cri_upper = self$attack_rate$cri_upper$ar,
        counts = private$house_counts_$total[
          match(self$attack_rate$observed$hhsize,
                private$house_counts_$hhsize)
        ]
      )
      # filter to household size bounds if applicable
      private$filter_hhsize_(output)
    },

    #' @field formatted_sec_attack_rate `data.frame`.
    #' Observed and predicted secondary attack rates by household size.
    #' Sizes outside the set household size bounds are not returned.
    formatted_sec_attack_rate = function(value) {
      if (!missing(value)) stop_read_only("formatted_sec_attack_rate")
      # convert internal list to data.frame
      output <- data.frame(
        hhsize = self$sec_attack_rate$observed$hhsize,
        observed = self$sec_attack_rate$observed$sar,
        predicted = self$sec_attack_rate$predicted$sar,
        cri_lower = self$sec_attack_rate$cri_lower$sar,
        cri_upper = self$sec_attack_rate$cri_upper$sar,
        counts = private$house_counts_$at_least_1_case[
          match(self$sec_attack_rate$observed$hhsize,
                private$house_counts_$hhsize)
        ]
      )
      # filter to household size bounds if applicable
      private$filter_hhsize_(output)
    },

    #' @field formatted_cases_count `data.frame`.
    #' Formatted observed and predicted cases by household size.
    #' Sizes outside the set household size bounds are not returned.
    formatted_cases_count = function(value) {
      if (!missing(value)) stop_read_only("formatted_cases_count")
      # filter according to household size bounds
      cases <- helper_filter_tables(self$cases_count,
                                    private$min_hhsize_, private$max_hhsize_)

      # empty frame as a placeholder for Credible Interval
      na_df <- cases$observed
      na_df[] <- lapply(na_df, function(col) ifelse(is.numeric(col), NA, col))
      na_df$hhsize <- cases$observed$hhsize

      list(
        predicted = helper_combine_obs_pred(cases$observed, cases$predicted),
        cri_lower = helper_combine_obs_pred(na_df, cases$cri_lower),
        cri_upper = helper_combine_obs_pred(na_df, cases$cri_upper)
      )
    }
  ),

  ############################################################################
  ################################## PUBLIC ##################################
  ############################################################################

  public = list(
    #' @param hhmodel `r rd_model_in_`
    #' @param param_table `r rd_param_table_`
    #' @param obs_process `r rd_simul_observ_`
    #' @param thin_ess `r rd_thin_ess_`
    #' @param max_nrep `r rd_max_nrep_`
    #' @param min_hhsize `r rd_min_hhsize_`
    #' @param max_hhsize `r rd_max_hhsize_`
    #' @param point_estimate `r rd_point_estimate_`
    #' @param threads `r rd_thread_`
    #' @param progress `r rd_progress_`
    initialize = function(hhmodel,
                          param_table,
                          obs_process = NULL,
                          thin_ess = TRUE,
                          max_nrep = 200L,
                          min_hhsize = NA_integer_,
                          max_hhsize = NA_integer_,
                          point_estimate = "median",
                          threads = 1L,
                          progress = TRUE) {
      # ininitialize parent class HousepiCustomPPC
      super$initialize(
        hhmodel = hhmodel, param_table = param_table,
        attack_rate = measure_ar_by_hhsize,
        sec_attack_rate = measure_sar_by_hhsize,
        cases_count = measure_cases_by_hhsize,
        obs_process = obs_process,
        thin_ess = thin_ess, max_nrep = max_nrep,
        point_estimate = point_estimate,
        threads = threads, progress = progress
      )
      private$house_counts_ <- helper_house_counts(hhmodel$data_wide)
      self$set_min_hhsize(min_hhsize)
      self$set_max_hhsize(max_hhsize)

      invisible(self)
    },

    #' @description Set minimum household size to consider in result outputs.
    #' @param n `integer(1)`. Minimum household size to consider.
    #' If `NA`, restore observed minimum.
    set_min_hhsize = function(n) {
      private$min_hhsize_ <- n
      invisible(self)
    },

    #' @description Restore minimum household size to observed minimum.
    drop_min_hhsize = function() {
      self$set_min_hhsize(NA_integer_)
      invisible(self)
    },

    #' @description Set maximum household size to consider in result outputs.
    #' @param n `integer(1)`. Maximum household size to consider.
    #' If `NA`, restore observed maximum.
    set_max_hhsize = function(n) {
      private$max_hhsize_ <- n
      invisible(self)
    },

    #' @description Restore maximum household size to observed maximum.
    drop_max_hhsize = function() {
      self$set_max_hhsize(NA_integer_)
      invisible(self)
    },

    #' @description Perform Chi-squared test and return list with results.
    #' @param lower_tail
    #' `logical(1)`. Consider lower tail probabilities `P[X<=x]` (default).
    #' If `FALSE`, consider upper tail `P[X>x]`.
    #' @return `list`
    chisq_test = function(lower_tail = FALSE) {
      cases <- helper_filter_tables(self$cases_count,
                                    private$min_hhsize_, private$max_hhsize_)
      res <- (cases$observed - cases$predicted)^2 / cases$predicted
      x_val <- sum(res, na.rm = TRUE)
      ddof <- prod(dim(res)) - nrow(res) - sum(is.na(res))
      list(chisq_value = x_val, ddof = ddof,
           p_value = stats::pchisq(x_val, df = ddof, lower.tail = lower_tail))
    }
  )
)

#' @title Custom posterior predictive check on a housepi model
#' @description
#' Posterior predictive check on user-supplied measures, comparing observed
#' and predicted values according to a given table of parameters.
#' @param hhmodel `r rd_model_in_`
#' @param param_table `r rd_param_table_`
#' @param ... `r rd_measure_dots_`
#' @param obs_process `r rd_simul_observ_`
#' @param obs_process_long `r rd_simul_obs_long_`
#' @param thin_ess `r rd_thin_ess_`
#' @param max_nrep `r rd_max_nrep_`
#' @param point_estimate `r rd_point_estimate_`
#' @param threads `r rd_thread_`
#' @param progress `r rd_progress_`
#' @return An object of class [HousepiCustomPPC].
#' @export
run_custom_predictive_check <- function(hhmodel,
                                        param_table,
                                        ...,
                                        obs_process = NULL,
                                        obs_process_long = NULL,
                                        thin_ess = TRUE,
                                        max_nrep = 200L,
                                        point_estimate = "median",
                                        threads = 1L,
                                        progress = TRUE) {
  stopifnothousepi(hhmodel)
  HousepiCustomPPC$new(
    hhmodel$copy(), param_table, ...,
    thin_ess = thin_ess, max_nrep = max_nrep,
    obs_process = obs_process, obs_process_long = obs_process_long,
    point_estimate = point_estimate, threads = threads, progress = progress
  )
}

#' @title Standard epidemiological posterior predictive check
#' @description
#' Posterior predictive check on standard epidemiological measures
#' by household size. It compares observed and predicted values of
#' \describe{
#' \item{`"attack_rate"`}{see [measure_ar_by_hhsize];}
#' \item{`"sec_attack_rate"`}{see [measure_sar_by_hhsize];}
#' \item{`"cases_count"`}{see [measure_cases_by_hhsize].}
#' }
#' @param hhmodel `r rd_model_in_`
#' @param param_table `r rd_param_table_`
#' @param obs_process `r rd_simul_observ_`
#' @param thin_ess `r rd_thin_ess_`
#' @param max_nrep `r rd_max_nrep_`
#' @param min_hhsize `r rd_min_hhsize_`
#' @param max_hhsize `r rd_max_hhsize_`
#' @param point_estimate `r rd_point_estimate_`
#' @param threads `r rd_thread_`
#' @param progress `r rd_progress_`
#' @return
#' An object of class [HousepiPPC]
#' containing methods to compare observed and predicted cases.
#' @export
run_epi_predictive_check <- function(hhmodel,
                                     param_table,
                                     obs_process = NULL,
                                     thin_ess = TRUE,
                                     max_nrep = 200L,
                                     min_hhsize = NA_integer_,
                                     max_hhsize = NA_integer_,
                                     point_estimate = "median",
                                     threads = 1L,
                                     progress = TRUE) {
  stopifnothousepi(hhmodel)
  HousepiPPC$new(
    hhmodel$copy(), param_table,
    obs_process = obs_process, thin_ess = thin_ess,
    max_nrep = max_nrep, min_hhsize = min_hhsize, max_hhsize = max_hhsize,
    point_estimate = point_estimate,
    threads = threads, progress = progress
  )
}

#' Format credible interval into printable string
#'
#' @description
#' Format element-wise estimate and credible-interval into a string.
#' The function is vectorised over its arguments
#' and returns a character vector of the same length as the inputs.
#' @param point_estimate `numeric()`. Point estimate values.
#' @param cri_lower `numeric()`. Lower bounds of the credible interval.
#' @param cri_upper `numeric()`. Upper bounds of the credible interval.
#' @param fmt `character(1)`.
#' Sprintf-style format with 3 '%s' placeholders.
#' Default is `"point_estimate [cri_lower, cri_upper]"`.
#' @param as_integer `logical(1)`. Whether values are formatted as integers.
#' Default is `FALSE`.
#' @return `character()`
#' @examples
#' pe <- c(9.9, 7.69e-02, 1.13)
#' cri_low <- c(2.6, 6.77e-02, NA)
#' cri_up <- c(23.7, 0.087, 1.53)
#' format_credible_interval(pe, cri_low, cri_up)
#' # custom format
#' format_credible_interval(pe, cri_low, cri_up,
#'                          fmt = "%s (95%% CrI, %s-%s)")[2]
#' @export
format_credible_interval <- function(point_estimate, cri_lower, cri_upper,
                                     fmt = "%s [%s, %s]",
                                     as_integer = FALSE) {
  .fmt_vector <- function(point_estimate, cri_lower, cri_upper) {
    if (length(point_estimate) != length(cri_lower) ||
          length(point_estimate) != length(cri_upper)) stop_h(
      "Point estimate and CrI bounds must be numeric of the same length."
    )

    if (length(gregexpr("%s", fmt, fixed = TRUE)[[1]]) != 3) stop_h(
      "Invalid choice for character {.arg str_format}.",
      "It must contain exactly three {.code %s} placeholders."
    )

    .fmt_num <- function(pe, low, up) {
      if (is.na(pe)) return(NA_character_)
      if (pe == low && pe == up) return(sprintf(fmt, pe, pe, pe))
      if (isTRUE(pe < low) || isTRUE(pe > up)) warn_h(sprintf(
        "Found an inconsistent credible interval: %g, CrI = [%g, %g].",
        pe, low, up
      ))

      if (isTRUE(as_integer)) digits <- 0
      else if (pe == 0) digits <- max(0, ceiling(-log10(abs(low))) + 1)
      else digits <- max(0, ceiling(-log10(abs(pe))) + 1)

      pe_f <- formatC(pe, digits = digits, format = "f")
      low_f <- ifelse(is.na(low), NA_character_,
                      formatC(low, digits = digits, format = "f"))
      up_f <- ifelse(is.na(up), NA_character_,
                     formatC(up, digits = digits, format = "f"))

      sprintf(fmt, pe_f, low_f, up_f)
    }

    mapply(.fmt_num, point_estimate, cri_lower, cri_upper, USE.NAMES = FALSE)
  }

  if (is.data.frame(point_estimate) || is.matrix(point_estimate)) {
    # check shapes
    if (!all(dim(point_estimate) == dim(cri_lower),
             dim(point_estimate) == dim(cri_upper))) stop_h(
      "Tabular inputs must have the same dimensions."
    )

    res <- as.data.frame(Map(.fmt_vector, as.data.frame(point_estimate),
                             as.data.frame(cri_lower),
                             as.data.frame(cri_upper)))
    names(res) <- names(point_estimate)
    return(res)
  }
  .fmt_vector(point_estimate, cri_lower, cri_upper)
}
