## housepi - Epidemiological Analysis of Household Outbreaks.
## Copyright 2024-2026 Institut Pasteur (Paris), France.
## SPDX-License-Identifier: GPL-3.0-or-later
##

rd_min_hh_size_ <- "`integer(1)` Minimum household size retained in the outputs.
If `NA`, no lower bound is applied."
rd_max_hh_size_ <- "`integer(1)` Maximum household size retained in the outputs.
If `NA`, no upper bound is applied."
rd_measure_dots_ <- "Named `function`s computing measures from data.
Each returns a `numeric()`, a `list` or a `data.frame`.
Names must be unique and cannot match a member of the returned object."

########  Measures  ########

# @todo check whether .data is still needed
#' @importFrom rlang .data
NULL

helper_process_data <- function(data_wide, escaped = "escaped") {
  df <- as.data.frame(data_wide)
  # an empty call returns the empty shape of the measure
  if (!ncol(df)) return(df)
  for (col in c("hh_id", "status")) if (!col %in% colnames(df)) stop_h(
    "Column {.field {col}} is missing from the data.",
    paste("Measures by household size need one row per individual, with a",
          "household identifier and a follow-up outcome.")
  )
  if (is.character(df$status) || is.factor(df$status))
    df$status <- status_to_codes_(df$status, escaped)$codes

  df$status <- as.integer(df$status)
  df$hh_id <- match(df$hh_id, unique(df$hh_id))
  df
}

#' @title Compute observed attack rate by household size
#' @description
#' Retrieve attack rates from household data.
#' @param data_wide `r rd_add_wide_`
#' @param escaped `r rd_mdl_escaped_`
#' @return `data.frame`
#' @details
#' The attack rate is computed for all households
#' with size \eqn{n \ge 1} according to
#' \eqn{\frac{\mathbb{E}_{n}[n^{+}]}{n}}
#' @examples
#' measure_ar_by_hh_size(flu_1918_frederick)
#' @export
measure_ar_by_hh_size <- function(data_wide, escaped = "escaped") {
  # normalise eagerly: a lazy argument would raise inside the engine wrapper,
  # which restyles the condition as if the engine had authored it
  df <- helper_process_data(data_wide, escaped)
  .engine$measure_ar_by_hh_size(df)
}

#' @title Compute observed secondary attack rate by household size
#' @description
#' Retrieve secondary attack rates (SAR) from household data.
#' @param data_wide `r rd_add_wide_`
#' @param escaped `r rd_mdl_escaped_`
#' @return `data.frame`
#' @details
#' The secondary attack rate is computed for all households
#' with size \eqn{n > 1} according to
#' \eqn{\frac{\mathbb{E}_{n}[n^{+}] - 1}{n - 1}}
#' @examples
#' measure_sar_by_hh_size(flu_1918_frederick)
#' @export
measure_sar_by_hh_size <- function(data_wide, escaped = "escaped") {
  df <- helper_process_data(data_wide, escaped)
  .engine$measure_sar_by_hh_size(df)
}

#' @title Tabulate observed cases distribution by household size
#' @description
#' Retrieve the distribution of household cases ("final size")
#' by household size from household data.
#' @param data_wide `r rd_add_wide_`
#' @param escaped `r rd_mdl_escaped_`
#' @return `data.frame` One row per household size `hh_size`,
#' then one column per number of cases, holding the number of households
#' of that size with that number of cases.
#' A number of cases above the household size is `NA`.
#' @examples
#' measure_cases_by_hh_size(flu_1918_frederick)
#' @export
measure_cases_by_hh_size <- function(data_wide, escaped = "escaped") {
  df <- helper_process_data(data_wide, escaped)
  .engine$measure_cases_by_hh_size(df)
}

#' @title Compute the case-weighted geometric mean of household size
#' @description
#' Retrieve the geometric mean of household size from household data,
#' weighted by the number of secondary cases in each household.
#' It is a data-driven choice for the reference household size `hh_size_ref`
#' of the built-in models.
#' @param data_wide `r rd_add_wide_`
#' @param escaped `r rd_mdl_escaped_`
#' @return `numeric(1)` The measure, `NA` if no household has two cases or more.
#' @details
#' For households \eqn{h} of size \eqn{n_h} with \eqn{n_h^{+}} cases,
#' the measure is
#' \eqn{\exp\left(\frac{\sum_h w_h \log n_h}{\sum_h w_h}\right)}
#' with weights \eqn{w_h = \max(n_h^{+} - 1, 0)}.
#' @examples
#' measure_geom_mean_hh_size(flu_1918_frederick)
#' @export
measure_geom_mean_hh_size <- function(data_wide, escaped = "escaped") {
  df <- helper_process_data(data_wide, escaped)
  .engine$measure_geom_mean_hh_size(df)
}

########  HousepiPredictiveCheck  ########

#' Apply `func` to corresponding elements of same-shaped data.frames.
#' @param frame_list `list` of `data.frame`s Same-shaped data.frames.
#' @param func `function` Applied to corresponding elements.
#' @param name `character(1)` Measure name.
#' @return `data.frame` Same shape, with the columns and row names of the first.
#' @noRd
ppc_frame_apply <- function(frame_list, func, name = NA_character_) {
  # ensure all elements share the same dimensions
  elem_n1 <- frame_list[[1]]
  dims_match <- vapply(frame_list[-1], \(f) identical(dim(f), dim(elem_n1)),
                       logical(1))
  if (!all(dims_match)) stop_h(
    paste0("Not all computed measures share identical dimensions",
           if (is.na(name)) "." else " in measure {.field {name}}."),
    "Check output of your measure function."
  )

  if (is.null(dim(elem_n1))) {
    # scalar or vector
    output <- sapply(seq_along(elem_n1), \(x) func(sapply(frame_list, `[[`, x)))
    names(output) <- names(elem_n1)
    if (is.list(elem_n1)) output <- as.list(output)

  } else {
    # matrix or data.frame
    # column-wise, for possibly mixed types
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

#' Apply a numeric statistic across corresponding elements of data.frames.
#' * every value NA -> NA;
#' * non-numeric, all identical -> that value.
#' @param frame_list `list` of `data.frame`s Same-shaped data.frames.
#' @param stat_fn `function` Numeric statistic applied to each element.
#' @param label `character(1)` Statistic name used in error messages.
#' @param name `character(1)` Measure name.
#' @noRd
ppc_frame_apply_stat <- function(frame_list, stat_fn, label,
                                 name = NA_character_) {
  ppc_frame_apply(frame_list, function(x) {
    if (all(is.na(x))) NA_real_
    else if (all(is.numeric(x))) stat_fn(x)
    else if (length(unique(x)) == 1) unique(x)
    else stop_h(
      paste0("Undefined {.fn {label}} for non-numeric, non-unique values",
             if (is.na(name)) "." else " in measure {.field {name}}."),
      "Check output of your measure function."
    )
  }, name = name)
}

#' Point estimate across same-shaped data.frames.
#' @param point_estimate `character(1)` Either `"mean"` or `"median"`.
#' @param name `character(1)` Measure name.
#' @noRd
ppc_frame_apply_summary <- function(frame_list,
                                    point_estimate = "median",
                                    name = NA_character_) {
  point_estimate <- validate_point_estimate(point_estimate)
  fn <- switch(point_estimate,
               mean = \(x) mean(x, na.rm = TRUE),
               median = \(x) stats::median(x, na.rm = TRUE))
  ppc_frame_apply_stat(frame_list, fn, label = point_estimate, name = name)
}

#' Credible interval bound across same-shaped data.frames.
#' @param cri_level `numeric(1)` Credible level, between 0.5 and 1.
#' @param upper `logical(1)` Whether to compute the upper bound, else the lower.
#' @param name `character(1)` Measure name.
#' @noRd
ppc_frame_apply_cri <- function(frame_list, cri_level = 0.9,
                                upper = TRUE, name = NA_character_) {
  p <- cri_level_probs(cri_level)[[if (isTRUE(upper)) "upper" else "lower"]]
  ppc_frame_apply_stat(frame_list,
                       \(x) stats::quantile(x, probs = p, na.rm = TRUE),
                       label = "credible interval",
                       name = name)
}

# Element-wise mask of observed values that fall outside [cri_lower, cri_upper].
# * non-numeric or NA element -> FALSE, not a mismatch.
ppc_observed_outside_cri <- function(x) {
  .vec_outside <- function(obs, low, up) {
    if (!is.numeric(obs) || !is.numeric(low) || !is.numeric(up))
      return(rep(FALSE, length(obs)))
    !is.na(obs) & !is.na(low) & !is.na(up) & (obs < low | obs > up)
  }
  if (!is.data.frame(x$observed))
    return(.vec_outside(x$observed, x$cri_lower, x$cri_upper))
  as.data.frame(Map(.vec_outside, x$observed, x$cri_lower, x$cri_upper))
}

# Extract predicted quantities from a list of cases-by-size tables.
ppc_extract_predicted_results <- function(x, point_estimate = "median", cri,
                                          name = NULL) {
  list(predicted = ppc_frame_apply_summary(x, point_estimate, name = name),
       cri_lower = ppc_frame_apply_cri(x, cri, upper = FALSE, name = name),
       cri_upper = ppc_frame_apply_cri(x, cri, upper = TRUE, name = name))
}

# Assert measure names are given, unique and free of the members of `self`,
# as each measure becomes a read-only field of the check object.
assert_measure_names_ <- function(nms, taken = character()) {
  if (is.null(nms) || anyNA(nms) || !all(nzchar(nms))) stop_h(
    "Every measure in {.arg ...} must be named."
  )
  duplicates <- unique(nms[duplicated(nms)])
  if (length(duplicates)) stop_h(
    "{cli::qty(duplicates)}Repeated measure name{?s} {.field {duplicates}}."
  )
  clash <- intersect(nms, taken)
  if (length(clash)) stop_h(
    "{cli::qty(clash)}Measure name{?s} {.field {clash}} {?is/are} reserved.",
    "Rename the measure."
  )
  invisible(nms)
}

#' @title Posterior predictive check container for user-defined measures
#' @description
#' Posterior predictive checks of measures on simulated surveys.
#' @examples
#' hh_model <- flu_1918_frederick |>
#'   new_housepi_model("final-size") |>
#'   add_parameter(import = "template")
#' hh_fit <- run_inference(hh_model, chains = 2, progress = FALSE)
#'
#' # a user-defined measure: the overall attack rate
#' overall_ar <- function(data_wide) mean(data_wide$status > 0)
#' ppc <- check_predictive(hh_model, hh_fit,
#'                         sar = measure_sar_by_hh_size,
#'                         overall = overall_ar,
#'                         epi = FALSE, nrep = 50, progress = FALSE)
#' ppc$overall
#' @docType class
HousepiPredictiveCheck <- R6::R6Class(
  classname = "HousepiPredictiveCheck",
  cloneable = TRUE,
  lock_objects = FALSE,

  ####  Private  ####

  private = list(
    deep_clone = function(name, value) helper_clone(value),
    cri_level_ = 0.9,
    point_estimate_ = "median",
    nrep_ = NA_integer_,
    thin_ess_ = TRUE,
    is_posterior_ = FALSE,
    model_name_ = NA_character_,
    secs_ = NA_real_,
    measure_names_ = character(0),
    # help topic of each built-in measure, by measure name
    measure_topics_ = character(0),
    measures_ = list(),
    # whether the simulated frames reach the measures with their labels
    # restored, as user-supplied measures always do
    labelled_ = TRUE,

    # Compare observed values against the credible interval of each measure,
    # in a single warning: one finding line per measure with observed entries
    # outside its [cri_lower, cri_upper] bounds, then one closing hint.
    check_observed_in_cri_ = function() {
      cri_pct <- round(100 * private$cri_level_)
      findings <- c()
      for (m in private$measure_names_) {
        mismatch <- ppc_observed_outside_cri(self[[m]])
        len_bad <- sum(unlist(mismatch))
        if (len_bad > 0L) findings <- c(findings, inline_cli_(
          "{cli::qty(len_bad)}There {?was/were} {len_bad} observed value{?s}",
          " outside of {cri_pct}% CrI for measure {.field {m}}."
        ))
      }
      if (length(findings)) {
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

  ####  Active  ####

  active = list(
    #' @field point_estimate `character(1)` Summary statistic name.
    point_estimate = function(value) {
      assert_read_only_(value, "point_estimate")
      private$point_estimate_
    },
    #' @field cri_level `numeric(1)` Credible interval level.
    cri_level = function(value) {
      assert_read_only_(value, "cri_level")
      private$cri_level_
    },
    #' @field nrep `integer(1)` Number of replicates used for the predictions.
    nrep = function(value) {
      assert_read_only_(value, "nrep")
      private$nrep_
    },
    #' @field thin_ess `logical(1)` Whether thinning to ESS was applied.
    thin_ess = function(value) {
      assert_read_only_(value, "thin_ess")
      private$thin_ess_
    },
    #' @field is_posterior `logical(1)`
    #' Whether the predictions are posterior predictions.
    is_posterior = function(value) {
      assert_read_only_(value, "is_posterior")
      private$is_posterior_
    },
    #' @field model_name `character(1)` Name of the underlying model.
    model_name = function(value) {
      assert_read_only_(value, "model_name")
      private$model_name_
    }
  ),

  ####  Public  ####

  public = list(
    #' @description
    #' Simulate from a model and parameter draws,
    #' then compute the predicted values of a set of measures.
    #' @param hh_model `r rd_model_in_`
    #' @param draws `r rd_draws_`
    #' @param ... `r rd_measure_dots_`
    #' @param nrep `r rd_nrep_`
    #' @param thin_ess `r rd_thin_ess_`
    #' @param measures_long `function` | `list` of `function`s
    #' Applied to `data_long`, derived by the `long` reporting process.
    #' @param threads `r rd_thread_`
    #' @param progress `r rd_progress_`
    initialize = function(hh_model,
                          draws,
                          ...,
                          measures_long = NULL,
                          nrep = NA_integer_,
                          thin_ess = TRUE,
                          threads = 1L,
                          progress = TRUE) {
      start <- proc.time()[["elapsed"]]
      nrep <- assert_count_(nrep, "Number of replicates", allow_na = TRUE)
      # associate summary statistic and credible interval
      private$point_estimate_ <- get_global_point_estimate()
      private$cri_level_ <- get_global_cri_level()

      # load and classify measure functions (wide vs long)
      m_wide <- list(...)
      if (length(m_wide))
        assert_measure_names_(names(m_wide), ls(self, all.names = TRUE))
      # @todo a list of functions breaks here
      m_long <- list()
      if (!is.null(measures_long)) m_long <- list(measures_long)

      reporting <- get_reporting_process(hh_model)
      if (length(m_long) && is.null(reporting$long)) stop_h(
        "Measures on {.code data_long} need a long reporting process.",
        "Add one with {.code add_reporting_process(long = )}."
      )

      # store (and test) results from observed data wide and/or long
      observations <- list()

      # wide-measure application
      apply_wide <- NULL
      if (length(m_wide)) {
        observations <- c(
          observations, lapply(m_wide, \(func) func(get_data_wide(hh_model)))
        )
        apply_wide <- \(x) lapply(m_wide, \(func) func(x))
      }
      # long-measure application
      apply_long <- NULL
      if (length(m_long)) {
        observations <- c(
          observations, lapply(m_long, \(func) func(get_data_long(hh_model)))
        )
        apply_long <- \(x) lapply(m_long, \(func) func(x))
      }

      # check measures
      if (any(lapply(observations, NROW) == 0)) stop_h(
        "One or more measures returned empty observation on data."
      )

      # run measures on posterior simulations
      output <- hh_model$.__enclos_env__$private$run_measure_(
        draws,
        apply_wide = apply_wide,
        apply_long = apply_long,
        nrep = nrep,
        thin_ess = thin_ess,
        threads = threads,
        progress = progress,
        labelled = private$labelled_
      )

      # parse results
      all_measure_names <- c(names(m_wide), names(m_long))

      parsed <- c()
      if (length(m_wide)) parsed <- c(parsed,
        lapply(names(m_wide), \(x) lapply(output$wide, `[[`, x))
      )
      if (length(m_long)) parsed <- c(parsed,
        lapply(names(m_long), \(x) lapply(output$long, `[[`, x))
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
            assert_read_only_(value, m_name)
            private$measures_[[m_name]]
          }
          makeActiveBinding(m_name, setter, self)
        })
      }
      # store decorations for print and read-only access
      private$model_name_ <- hh_model$name
      private$is_posterior_ <-
        inherits(draws, "HousepiPosteriorDraws") || posterior::is_draws(draws)
      private$measure_names_ <- all_measure_names
      private$nrep_ <- length(parsed[[1]])
      private$thin_ess_ <- isTRUE(thin_ess)

      private$check_observed_in_cri_()
      private$secs_ <- proc.time()[["elapsed"]] - start

      invisible(self)
    },

    #' @description
    #' Print the measures, observed against predicted.
    print = function() {
      kind <- if (private$is_posterior_) "Posterior predictions" else
        "Predictions"
      # built-in measures link to their help page
      measures <- vapply(private$measure_names_, function(m) {
        topic <- private$measure_topics_[m]
        if (is.na(topic)) m
        else cli::format_inline(
          "{.run [{m}](utils::help(\"{topic}\", \"housepi\"))}"
        )
      }, character(1), USE.NAMES = FALSE)
      print_rows_("{kind} for model {.val {private$model_name_}}.", c(
        Measures = "{measures}",
        Replicates = paste0("{private$nrep_}",
                            if (private$thin_ess_) " (ESS)" else ""),
        `Point estimate` = "{private$point_estimate_}",
        `Credible interval` = "{round(100 * private$cri_level_)}%",
        `Run time` = "{format_secs_(private$secs_)}"
      ))
      invisible(self)
    }
  )
)

########  HousepiEpiPredictiveCheck  ########

# Total households and households with at least one case, per household size.
helper_house_counts <- function(data_wide, escaped = "escaped") {
  df <- helper_process_data(data_wide, escaped)
  .engine$count_houses_by_hh_size(df)
}

# List every case table, filtered.
helper_filter_tables <- function(cases_distr,
                                 min_hh_size = NA,
                                 max_hh_size = NA) {
  # load observed table to use as a type shape template
  obsdf <- cases_distr$observed

  # remove rows without counts
  mask_row <- rowSums(obsdf[colnames(obsdf) != "hh_size"], na.rm = TRUE) != 0

  # remove column with zero cases if none was observed
  mask_col <- colnames(obsdf) != "0" |
    colSums(obsdf == 0, na.rm = TRUE) != nrow(obsdf)

  # implement max household size mask (also caps the case columns, as a
  # household cannot hold more cases than its size)
  if (!is.na(max_hh_size)) {
    col_max <- which(colnames(obsdf) == max_hh_size)
    mask_col <- mask_col & (seq_len(ncol(obsdf)) <= col_max)
    mask_row <- mask_row & (obsdf$hh_size <= max_hh_size)
  }

  # implement min household size mask (rows only: any case column stays valid)
  if (!is.na(min_hh_size)) mask_row <- mask_row & (obsdf$hh_size >= min_hh_size)

  list(
    observed = obsdf[mask_row, mask_col],
    predicted = cases_distr$predicted[mask_row, mask_col],
    cri_lower = cases_distr$cri_lower[mask_row, mask_col],
    cri_upper = cases_distr$cri_upper[mask_row, mask_col]
  )
}

# Bind the observed and predicted case tables.
# @param obsdf `data.frame` Observed cases.
# @param preddf `data.frame` Predicted cases.
helper_combine_obs_pred <- function(obsdf, preddf) {
  # 1. curate observed dataframe
  obsdf$distr <- "observed"
  rownames(obsdf) <- seq(1, 2 * nrow(obsdf), 2)

  # 2. curate predicted dataframe
  preddf$hh_size <- obsdf$hh_size
  preddf$distr <- "predicted"
  rownames(preddf) <- seq(2, 2 * nrow(preddf), 2)

  # 3. combine by row bind
  combined <- rbind(obsdf, preddf)[as.character(1:(2 * nrow(obsdf))), ]
  # reset row names
  rownames(combined) <- NULL
  # reorder columns
  y <- colnames(combined)
  aux <- c("hh_size", "distr")
  combined <- combined[, c(aux, y[!y %in% aux])]

  # format columns to extract labels with labelizor
  colnames(combined)[3:ncol(combined)] <- sapply(
    colnames(combined)[3:ncol(combined)], \(x) paste("cases", x, sep = "_")
  )

  combined
}

#' @title Standard epidemiological posterior predictive check
#' @description
#' Posterior predictive check of standard epidemiological measures
#' by household size, observed against predicted:
#' cases distribution, attack rate and secondary attack rate.
#' @examples
#' hh_model <- flu_1918_frederick |>
#'   new_housepi_model("final-size") |>
#'   add_parameter(import = "template")
#' hh_fit <- run_inference(hh_model, chains = 2, progress = FALSE)
#'
#' ppc <- check_predictive(hh_model, hh_fit, nrep = 50, progress = FALSE)
#' ppc$format_attack_rate()
#' @docType class
HousepiEpiPredictiveCheck <- R6::R6Class(
  classname = "HousepiEpiPredictiveCheck",
  inherit = HousepiPredictiveCheck,
  cloneable = TRUE,
  lock_objects = FALSE,

  ####  Private  ####

  private = list(
    deep_clone = function(name, value) helper_clone(value),
    min_hh_size_ = NA_integer_,
    max_hh_size_ = NA_integer_,
    house_counts_ = NULL,
    measure_topics_ = c(attack_rate = "measure_ar_by_hh_size",
                        sec_attack_rate = "measure_sar_by_hh_size",
                        cases_distr = "measure_cases_by_hh_size"),

    # Restrict a by-hh_size data.frame to the set [min_hh_size_, max_hh_size_].
    filter_hh_size_ = function(df) {
      if (!is.na(private$min_hh_size_))
        df <- df[df$hh_size >= private$min_hh_size_, ]
      if (!is.na(private$max_hh_size_))
        df <- df[df$hh_size <= private$max_hh_size_, ]
      df
    }
  ),

  ####  Public  ####

  public = list(
    #' @param hh_model `r rd_model_in_`
    #' @param draws `r rd_draws_`
    #' @param ... `r rd_measure_dots_`
    #' They are checked alongside standard epidemiological measures.
    #' @param nrep `r rd_nrep_`
    #' @param thin_ess `r rd_thin_ess_`
    #' @param min_hh_size `r rd_min_hh_size_`
    #' @param max_hh_size `r rd_max_hh_size_`
    #' @param threads `r rd_thread_`
    #' @param progress `r rd_progress_`
    initialize = function(hh_model,
                          draws,
                          ...,
                          nrep = NA_integer_,
                          thin_ess = TRUE,
                          min_hh_size = NA_integer_,
                          max_hh_size = NA_integer_,
                          threads = 1L,
                          progress = TRUE) {
      clash <- intersect(names(list(...)), names(private$measure_topics_))
      if (length(clash)) stop_h(
        paste("{cli::qty(clash)}Measure name{?s} {.field {clash}} already used",
              "by a standard epidemiological measure."),
        "Rename the measure or set {.code epi = FALSE}."
      )
      # built-in measures read only hh_id and status:
      # * no wide reporting process, no custom measure -> frames stay
      #   engine-coded, the relabelling skipped once per replicate.
      esc <- hh_model$escaped
      private$labelled_ <- !is.null(get_reporting_process(hh_model)$wide) ||
        length(list(...))

      # initialize parent class HousepiPredictiveCheck
      super$initialize(
        hh_model = hh_model, draws = draws,
        attack_rate = \(x) measure_ar_by_hh_size(x, esc),
        sec_attack_rate = \(x) measure_sar_by_hh_size(x, esc),
        cases_distr = \(x) measure_cases_by_hh_size(x, esc),
        ...,
        nrep = nrep, thin_ess = thin_ess,
        threads = threads, progress = progress
      )
      private$house_counts_ <- helper_house_counts(get_data_wide(hh_model), esc)
      self$set_min_hh_size(min_hh_size)
      self$set_max_hh_size(max_hh_size)

      invisible(self)
    },

    #' @description Set minimum household size to consider in result outputs.
    #' @param n `integer(1)` Minimum household size to consider.
    #' If `NA`, restore observed minimum.
    set_min_hh_size = function(n) {
      private$min_hh_size_ <- n
      invisible(self)
    },

    #' @description Restore minimum household size to observed minimum.
    drop_min_hh_size = function() {
      self$set_min_hh_size(NA_integer_)
      invisible(self)
    },

    #' @description Set maximum household size to consider in result outputs.
    #' @param n `integer(1)` Maximum household size to consider.
    #' If `NA`, restore observed maximum.
    set_max_hh_size = function(n) {
      private$max_hh_size_ <- n
      invisible(self)
    },

    #' @description Restore maximum household size to observed maximum.
    drop_max_hh_size = function() {
      self$set_max_hh_size(NA_integer_)
      invisible(self)
    },

    #' @description Format observed and predicted attack rates
    #' by household size.
    #' Sizes outside the set household size bounds are not returned.
    #' @return `data.frame`
    format_attack_rate = function() {
      # convert internal list to data.frame
      output <- data.frame(
        hh_size = self$attack_rate$observed$hh_size,
        observed = self$attack_rate$observed$ar,
        predicted = self$attack_rate$predicted$ar,
        cri_lower = self$attack_rate$cri_lower$ar,
        cri_upper = self$attack_rate$cri_upper$ar,
        counts = private$house_counts_$total[
          match(self$attack_rate$observed$hh_size,
                private$house_counts_$hh_size)
        ]
      )
      # filter to household size bounds if applicable
      private$filter_hh_size_(output)
    },

    #' @description Format observed and predicted secondary attack rates
    #' by household size.
    #' Sizes outside the set household size bounds are not returned.
    #' @return `data.frame`
    format_sec_attack_rate = function() {
      # convert internal list to data.frame
      output <- data.frame(
        hh_size = self$sec_attack_rate$observed$hh_size,
        observed = self$sec_attack_rate$observed$sar,
        predicted = self$sec_attack_rate$predicted$sar,
        cri_lower = self$sec_attack_rate$cri_lower$sar,
        cri_upper = self$sec_attack_rate$cri_upper$sar,
        counts = private$house_counts_$at_least_1_case[
          match(self$sec_attack_rate$observed$hh_size,
                private$house_counts_$hh_size)
        ]
      )
      # filter to household size bounds if applicable
      private$filter_hh_size_(output)
    },

    #' @description Format observed and predicted cases distributions
    #' by household size.
    #' Sizes outside the set household size bounds are not returned.
    #' @return `list` of `data.frame`: `predicted`, `cri_lower`, `cri_upper`.
    format_cases_distr = function() {
      # filter according to household size bounds
      cases <- helper_filter_tables(self$cases_distr,
                                    private$min_hh_size_, private$max_hh_size_)

      # empty frame as a placeholder for Credible Interval
      na_df <- cases$observed
      na_df[] <- lapply(na_df, \(col) if (is.numeric(col)) NA else col)
      na_df$hh_size <- cases$observed$hh_size

      list(
        predicted = helper_combine_obs_pred(cases$observed, cases$predicted),
        cri_lower = helper_combine_obs_pred(na_df, cases$cri_lower),
        cri_upper = helper_combine_obs_pred(na_df, cases$cri_upper)
      )
    }
  )
)

#' @title Posterior predictive check of a housepi model
#' @description
#' Posterior predictive check comparing observed values of household measures
#' against their values on surveys simulated from parameter draws
#' \insertCite{gelman_posterior_1996}{housepi}.
#' It checks standard epidemiological measures by household size:
#' \describe{
#' \item{`"attack_rate"`}{see [measure_ar_by_hh_size];}
#' \item{`"sec_attack_rate"`}{see [measure_sar_by_hh_size];}
#' \item{`"cases_distr"`}{see [measure_cases_by_hh_size].}
#' }
#' Named measure functions in `...` are checked alongside them,
#' under names that differ from these.
#' Each replicate goes through the reporting process of the model,
#' see [add_reporting_process()].
#' @param hh_model `r rd_model_in_`
#' @param draws `r rd_draws_`
#' @param ... `r rd_measure_dots_`
#' @param nrep `r rd_nrep_`
#' @param thin_ess `r rd_thin_ess_`
#' @param epi `logical(1)`
#' Whether to check standard epidemiological measures.
#' If `FALSE`, only the measures in `...` are checked,
#' and the formatted outputs and household size bounds of
#' [HousepiEpiPredictiveCheck] are unavailable.
#' @param min_hh_size `r rd_min_hh_size_` Only with `epi = TRUE`.
#' @param max_hh_size `r rd_max_hh_size_` Only with `epi = TRUE`.
#' @param threads `r rd_thread_`
#' @param progress `r rd_progress_`
#' @return
#' An object of class [HousepiEpiPredictiveCheck] if `epi = TRUE`,
#' of class [HousepiPredictiveCheck] otherwise.
#' @details
#' `r rd_global_summary_`
#' @references
#' \insertAllCited{}
#' @examples
#' hh_model <- flu_1918_frederick |>
#'   new_housepi_model("final-size") |>
#'   add_parameter(import = "template")
#' hh_fit <- run_inference(hh_model, chains = 2, progress = FALSE)
#'
#' ppc <- check_predictive(hh_model, hh_fit, nrep = 50, progress = FALSE)
#' ppc
#'
#' # a user-defined measure: the overall attack rate
#' overall_ar <- function(data_wide) mean(data_wide$status > 0)
#' ppc <- check_predictive(hh_model, hh_fit, overall = overall_ar,
#'                         nrep = 50, progress = FALSE)
#' ppc
#' ppc$overall$observed
#'
#' # only the user-defined measures
#' ppc <- check_predictive(hh_model, hh_fit, overall = overall_ar,
#'                         epi = FALSE, nrep = 50, progress = FALSE)
#' ppc
#' @export
check_predictive <- function(hh_model,
                             draws,
                             ...,
                             epi = TRUE,
                             nrep = NA_integer_,
                             thin_ess = TRUE,
                             min_hh_size = NA_integer_,
                             max_hh_size = NA_integer_,
                             threads = 1L,
                             progress = TRUE) {
  stopifnothousepi(hh_model)
  if (!rlang::is_bool(epi)) stop_h(
    "{.arg epi} should be {.code TRUE} or {.code FALSE}."
  )
  if (isTRUE(epi)) return(HousepiEpiPredictiveCheck$new(
    hh_model$copy(), draws, ..., nrep = nrep,
    thin_ess = thin_ess, min_hh_size = min_hh_size, max_hh_size = max_hh_size,
    threads = threads, progress = progress
  ))
  if (!length(list(...))) stop_h(
    "No measure to check.",
    "Pass measure functions in {.arg ...} or keep {.code epi = TRUE}."
  )
  if (!all(is.na(c(min_hh_size, max_hh_size)))) stop_h(
    paste("Arguments {.arg min_hh_size} and {.arg max_hh_size} bound only",
          "standard epidemiological measures."),
    "Set {.code epi = TRUE} or restrict household sizes inside the measures."
  )
  HousepiPredictiveCheck$new(
    hh_model$copy(), draws, ...,
    nrep = nrep, thin_ess = thin_ess,
    threads = threads, progress = progress
  )
}

#' Format credible interval into printable string
#'
#' @description
#' Format element-wise estimate and credible-interval into a string.
#' The function is vectorised over its arguments
#' and returns a character vector of the same length as the inputs.
#' @param point_estimate `numeric()` Point estimate values.
#' @param cri_lower `numeric()` Lower bounds of the credible interval.
#' @param cri_upper `numeric()` Upper bounds of the credible interval.
#' @param fmt `character(1)`
#' Sprintf-style format with 3 `%s` placeholders,
#' filled with the point estimate, then the lower and upper bounds.
#' @param as_integer `logical(1)` Whether values are formatted as integers.
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
      if (isTRUE(pe < low) || isTRUE(pe > up)) warn_h(
        "Found an inconsistent credible interval: {pe}, CrI = [{low}, {up}]."
      )

      if (isTRUE(as_integer)) digits <- 0
      else if (pe == 0) digits <- max(0, ceiling(-log10(abs(low))) + 1)
      else digits <- max(0, ceiling(-log10(abs(pe))) + 1)

      pe_f <- formatC(pe, digits = digits, format = "f")
      low_f <- if (is.na(low)) NA_character_
      else formatC(low, digits = digits, format = "f")
      up_f <- if (is.na(up)) NA_character_
      else formatC(up, digits = digits, format = "f")

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
