## housepi - Epidemiological Analysis of Household Outbreaks.
## Copyright 2024-2026 Institut Pasteur (Paris), France.
## SPDX-License-Identifier: GPL-3.0-or-later
##

rd_min_hh_size_ <- "`integer(1)` Minimum household size retained in the outputs.
If `NA` or `NULL`, no lower bound is applied."
rd_max_hh_size_ <- "`integer(1)` Maximum household size retained in the outputs.
If `NA` or `NULL`, no upper bound is applied."
rd_measure_dots_ <- "Named `function`s computing measures from data.
Each returns a `numeric()`, a `list` or a `data.frame`.
Names must be unique and cannot match a member of the returned object.
A measure with an `escaped` argument receives the model's escaped label."
rd_long_measures_ <- "`list`
Named `function`s computing measures from the long dataset.
Each replicate's long dataset comes from the `long_dataset` reporting process.
Arguments, returned values and names follow the rules of `...`."

########  Measures  ########

# @todo check whether .data is still needed
#' @importFrom rlang .data
NULL

helper_process_data <- function(dataset, escaped = "escaped") {
  df <- as.data.frame(dataset)
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

#' @title Measure attack rate by household size
#' @description
#' Retrieve attack rates from household data.
#' @param data `r rd_dataset_`
#' @param escaped `r rd_mdl_escaped_`
#' @return `data.frame`
#' @details
#' The attack rate is computed for all households
#' with size \eqn{n \ge 1} according to
#' \eqn{\frac{\mathbb{E}_{n}[n^{+}]}{n}}
#' @examples
#' measure_ar_by_hh_size(flu_1918_frederick)
#' @export
measure_ar_by_hh_size <- function(data, escaped = "escaped") {
  # normalise eagerly: a lazy argument would raise inside the engine wrapper,
  # which restyles the condition as if the engine had authored it
  df <- helper_process_data(data, escaped)
  .engine$measure_ar_by_hh_size(df)
}

#' @title Measure secondary attack rate by household size
#' @description
#' Retrieve secondary attack rates (SAR) from household data.
#' @param data `r rd_dataset_`
#' @param escaped `r rd_mdl_escaped_`
#' @return `data.frame`
#' @details
#' The secondary attack rate is computed for all households
#' with size \eqn{n > 1} according to
#' \eqn{\frac{\mathbb{E}_{n}[n^{+}] - 1}{n - 1}}
#' @examples
#' measure_sar_by_hh_size(flu_1918_frederick)
#' @export
measure_sar_by_hh_size <- function(data, escaped = "escaped") {
  df <- helper_process_data(data, escaped)
  .engine$measure_sar_by_hh_size(df)
}

#' @title Measure cases distribution by household size
#' @description
#' Retrieve the distribution of household cases ("final size")
#' by household size from household data.
#' @param data `r rd_dataset_`
#' @param escaped `r rd_mdl_escaped_`
#' @return `data.frame` One row per household size `hh_size`,
#' then one column per number of cases, holding the number of households
#' of that size with that number of cases.
#' A number of cases above the household size is `NA`.
#' @examples
#' measure_cases_by_hh_size(flu_1918_frederick)
#' @export
measure_cases_by_hh_size <- function(data, escaped = "escaped") {
  df <- helper_process_data(data, escaped)
  .engine$measure_cases_by_hh_size(df)
}

########  HousepiPredictiveCheck  ########

#' Apply `func` to corresponding elements of same-shaped data.frames.
#' @param frame_list `list` of `data.frame`s Same-shaped data.frames.
#' @param func `function` Applied to corresponding elements.
#' @param name `character(1)` Measure name.
#' @return `data.frame` Same shape, with the columns and row names of the first.
#' @noRd
ppc_frame_apply <- function(frame_list, func, name = NA_character_) {
  # ensure all elements share the same dimensions, vectors their length
  elem_n1 <- frame_list[[1]]
  same_shape <- function(f) {
    identical(dim(f), dim(elem_n1)) && length(f) == length(elem_n1)
  }
  dims_match <- vapply(frame_list[-1], same_shape, logical(1))
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

    output <- as.data.frame(lapply(seq_len(ncol(elem_n1)), function(x) {
      apply(matrix(unlist(lapply(frame_list, \(f) f[, x])), nrow = nr),
            1L, func)
    }), col.names = seq_len(ncol(elem_n1)))
    if (!is.data.frame(elem_n1)) output <- as.matrix(output)
    dimnames(output) <- dimnames(elem_n1)
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
# @param measures `list` Measures, keyed by the argument passing them.
assert_measure_names_ <- function(measures, taken = character()) {
  for (arg in names(measures)) {
    nms <- rlang::names2(measures[[arg]])
    if (anyNA(nms) || !all(nzchar(nms))) stop_h(
      "Every measure in {.arg {arg}} must be named."
    )
  }
  nms <- unlist(lapply(unname(measures), names))
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

# Split the `...` of a predictive check into measures
# and run settings, as `run_simulation()` takes them.
measure_dots_ <- function(...) {
  dots <- list(...)
  setting <- rlang::names2(dots) %in% "time_step"
  list(measures = dots[!setting], settings = dots[setting])
}

#' Whether `x` holds posterior draws, from a run or as `posterior::draws`.
#' @noRd
is_draws_ <- function(x) {
  inherits(x, "HousepiPosteriorDraws") || posterior::is_draws(x)
}

#' Source of the parameter values a check simulates at,
#' `hh_model` built at the wall.
#' * `NULL` -> `"prior"`;
#' * run fitted to `hh_model` and its data, or `posterior::draws`
#'   -> `"posterior"`;
#' * run fitted to another model or data -> `"table"`;
#' * table or named values -> `"table"`.
#' @noRd
draws_source_ <- function(draws, hh_model) {
  if (is.null(draws)) return("prior")
  if (!is_draws_(draws)) return("table")
  if (!inherits(draws, "HousepiPosteriorDraws")) return("posterior")
  differ <- differing_entries_(hh_model$.__enclos_env__$private$fingerprint_,
                               draws$.__enclos_env__$private$fingerprint_)
  if (length(differ)) "table" else "posterior"
}

#' @title Prior or posterior predictive check of user-defined measures
#' @description
#' Prior or posterior predictive checks of measures
#' on simulated household studies.
#' @examples
#' hh_model <- flu_1918_frederick |>
#'   new_housepi_model("final-size") |>
#'   add_parameter(import = "template")
#' hh_fit <- run_inference(hh_model)
#'
#' # a user-defined measure: the overall attack rate
#' overall_ar <- function(data) mean(data$status > 0)
#' ppc <- check_predictive(hh_model, hh_fit,
#'                         sar = measure_sar_by_hh_size,
#'                         overall = overall_ar,
#'                         epi = FALSE, nrep = 50)
#' ppc$overall
#' @docType class
HousepiPredictiveCheck <- R6::R6Class(
  classname = "HousepiPredictiveCheck",
  cloneable = TRUE,
  lock_objects = FALSE,

  ####  Private  ####

  private = list(
    cri_level_ = 0.9,
    point_estimate_ = "median",
    nrep_ = NA_integer_,
    thin_ess_ = TRUE,
    source_ = NA_character_,
    model_name_ = NA_character_,
    secs_ = NA_real_,
    warnings_ = NULL,
    measure_names_ = character(0),
    # help topic of each built-in measure, by measure name
    measure_topics_ = character(0),
    measures_ = list(),
    # whether the simulated frames reach the measures with their labels
    # restored, as user-supplied measures always do
    labelled_ = TRUE,

    # Name the check after the source of its parameter values.
    kind_ = function() {
      switch(private$source_,
        prior = "Prior predictive check",
        posterior = "Posterior predictive check",
        table = "Predictive check"
      )
    },

    # Compare observed values against the credible interval of each measure,
    # in a single warning: one finding line per measure with observed entries
    # outside its [cri_lower, cri_upper] bounds, then one closing hint.
    # Each finding is recorded as a row of `warnings_`, raised once.
    check_observed_in_cri_ = function() {
      cri_pct <- round(100 * private$cri_level_)
      findings <- character()
      for (m in private$measure_names_) {
        mismatch <- ppc_observed_outside_cri(self[[m]])
        len_bad <- sum(unlist(mismatch))
        if (len_bad > 0L) findings <- c(findings, inline_cli_(
          "{cli::qty(len_bad)}There {?was/were} {len_bad} observed value{?s}",
          " outside of {cri_pct}% CrI for measure {.field {m}}."
        ))
      }
      private$warnings_ <- data.frame(
        message = cli::ansi_strip(findings), count = rep(1L, length(findings)),
        stringsAsFactors = FALSE
      )
      if (length(findings)) {
        names(findings) <- rep("!", length(findings))
        kind <- private$kind_()
        warn_h(c(
          "{kind} flagged the following:",
          findings,
          "Observed data may be poorly captured by the model."
        ), .class = "housepi_warning_ppc")
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
    #' @field thin_ess `logical(1)`
    #' Whether posterior draws were thinned to their ESS.
    thin_ess = function(value) {
      assert_read_only_(value, "thin_ess")
      private$thin_ess_
    },
    #' @field source `r rd_check_source_()`
    source = function(value) {
      assert_read_only_(value, "source")
      private$source_
    },
    #' @field model_name `character(1)` Name of the underlying model.
    model_name = function(value) {
      assert_read_only_(value, "model_name")
      private$model_name_
    },
    #' @field warnings `r rd_warnings_`
    #' They hold the measures with observed values outside their credible
    #' interval, compared in the field of each measure.
    warnings = function(value) {
      assert_read_only_(value, "warnings")
      private$warnings_
    }
  ),

  ####  Public  ####

  public = list(
    #' @description
    #' Simulate from a model and parameter draws, or from its prior,
    #' then compute the predicted values of a set of measures.
    #' @param hh_model `r rd_model_in_`
    #' @param draws `r rd_draws_`
    #' @param ... `r rd_measure_dots_`
    #' @param nrep `r rd_nrep_`
    #' `r rd_nrep_prior_` `NA` is then an error.
    #' @param thin_ess `r rd_thin_ess_`
    #' @param long_measures `r rd_long_measures_`
    #' @param threads `r rd_thread_`
    #' @param progress `r rd_progress_`
    initialize = function(hh_model,
                          draws = NULL,
                          ...,
                          long_measures = NULL,
                          nrep = NA_integer_,
                          thin_ess = TRUE,
                          threads = 1L,
                          progress = TRUE) {
      start <- proc.time()[["elapsed"]]
      nrep <- assert_count_(nrep, "Number of replicates", allow_na = TRUE)
      if (is.null(draws) && is.na(nrep)) stop_h(
        "A prior predictive check needs {.arg nrep}.",
        "Set the number of replicates, or pass posterior {.arg draws}."
      )
      # associate summary statistic and credible interval
      private$point_estimate_ <- get_global_point_estimate()
      private$cri_level_ <- get_global_cri_level()

      # measures on the dataset and on the long dataset
      dots <- measure_dots_(...)
      measures <- list(dataset = dots$measures, long_dataset = long_measures)
      assert_measure_names_(
        stats::setNames(measures, c("...", "long_measures")),
        ls(self, all.names = TRUE)
      )
      reporting <- get_reporting_process(hh_model)
      if (length(measures$long_dataset) && is.null(reporting$long_dataset))
        stop_h(
          "Measures on the long dataset need a long dataset reporting process.",
          "Add one with {.code add_reporting_process(long_dataset = )}."
        )

      # measure with an `escaped` formal -> called with the model's label
      esc <- hh_model$escaped
      bind_escaped <- function(f) {
        if (!is.function(f) || !"escaped" %in% names(formals(f))) return(f)
        \(x) f(x, escaped = esc)
      }
      measures <- lapply(measures, \(fs) lapply(fs, bind_escaped))

      # each side applies its measures to a frame
      # * no long measure -> long dataset reporting process skipped
      appliers <- lapply(measures, \(fs) function(x) lapply(fs, \(f) f(x)))
      if (!length(measures$long_dataset)) appliers$long_dataset <- NULL
      observations <- c(
        appliers$dataset(get_dataset(hh_model)),
        if (length(measures$long_dataset))
          appliers$long_dataset(get_long_dataset(hh_model))
      )
      if (any(lapply(observations, NROW) == 0)) stop_h(
        "One or more measures returned empty observation on data."
      )

      # run measures on simulations
      output <- do.call(hh_model$.__enclos_env__$private$simulate_, c(list(
        draws = draws,
        nrep = nrep,
        thin_ess = thin_ess,
        apply_dataset = appliers$dataset,
        apply_long_dataset = appliers$long_dataset,
        labelled = private$labelled_,
        progress = progress,
        threads = threads
      ), dots$settings))

      # values of each measure, one per replicate
      replicates <- output$dataset
      if (!is.null(output$long_dataset))
        replicates <- Map(c, replicates, output$long_dataset)
      all_measure_names <- names(observations)
      parsed <- lapply(stats::setNames(nm = all_measure_names),
                       \(m) lapply(replicates, `[[`, m))

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
      # simulation built the model at the wall, stamping its fingerprint
      private$source_ <- draws_source_(draws, hh_model)
      private$measure_names_ <- all_measure_names
      private$nrep_ <- length(replicates)
      private$thin_ess_ <- isTRUE(thin_ess) && is_draws_(draws)

      private$check_observed_in_cri_()
      private$secs_ <- proc.time()[["elapsed"]] - start

      invisible(self)
    },

    #' @description
    #' Print the measures, observed against predicted.
    print = function() {
      # built-in measures link to their help page
      measures <- vapply(private$measure_names_, function(m) {
        topic <- private$measure_topics_[m]
        if (is.na(topic)) m
        else cli::format_inline(
          "{.run [{m}](utils::help(\"{topic}\", \"housepi\"))}"
        )
      }, character(1), USE.NAMES = FALSE)
      kind <- private$kind_()
      print_rows_("{kind} for model {.val {private$model_name_}}.", c(
        Measures = "{measures}",
        Replicates = "{private$nrep_}",
        `Point estimate` = "{private$point_estimate_}",
        `Credible interval` = "{round(100 * private$cri_level_)}%",
        `Run time` = "{format_secs_(private$secs_)}"
      ), warnings = private$warnings_)
      invisible(self)
    }
  )
)

########  HousepiEpiPredictiveCheck  ########

# Total households and households with at least one case, per household size.
helper_house_counts <- function(dataset, escaped = "escaped") {
  df <- helper_process_data(dataset, escaped)
  .engine$count_houses_by_hh_size(df)
}

# Restrict a by-hh_size data.frame to sizes in [min_hh_size, max_hh_size],
# rows renumbered.
# * `NA` bound -> no bound on that side.
ppc_filter_hh_size <- function(df, min_hh_size = NA, max_hh_size = NA) {
  keep <- (is.na(min_hh_size) | df$hh_size >= min_hh_size) &
    (is.na(max_hh_size) | df$hh_size <= max_hh_size)
  `rownames<-`(df[keep, ], NULL)
}

# Format a rate by household size, observed against predicted.
# @param measure `list` Measure: `observed`, `predicted`, `cri_lower`,
# `cri_upper`.
# @param rate `character(1)` Column holding the rate in each frame.
# @param counts `data.frame` Household size, then number of households.
ppc_format_rate <- function(measure, rate, counts,
                            min_hh_size = NA, max_hh_size = NA) {
  hh_size <- measure$observed$hh_size
  df <- data.frame(
    hh_size = hh_size,
    observed = measure$observed[[rate]],
    predicted = measure$predicted[[rate]],
    cri_lower = measure$cri_lower[[rate]],
    cri_upper = measure$cri_upper[[rate]]
  )
  df$outside <- ppc_observed_outside_cri(df)
  df$counts <- counts[[2]][match(hh_size, counts$hh_size)]
  ppc_filter_hh_size(df, min_hh_size, max_hh_size)
}

# Format the cases distribution by household size, observed against predicted,
# one row per household size and number of cases.
# * more cases than members -> row dropped;
# * household size with no household -> rows dropped;
# * no household without cases -> rows with zero cases dropped.
ppc_format_cases_distr <- function(measure, min_hh_size = NA,
                                   max_hh_size = NA) {
  cases <- setdiff(names(measure$observed), "hh_size")
  .stack <- \(x) unlist(x[cases], use.names = FALSE)
  df <- data.frame(
    hh_size = rep(measure$observed$hh_size, length(cases)),
    cases = rep(as.integer(cases), each = nrow(measure$observed)),
    observed = .stack(measure$observed),
    predicted = .stack(measure$predicted),
    cri_lower = .stack(measure$cri_lower),
    cri_upper = .stack(measure$cri_upper)
  )
  df <- df[!is.na(df$observed), ]
  houses <- tapply(df$observed, df$hh_size, sum)
  df <- df[df$hh_size %in% names(houses)[houses > 0], ]
  if (all(df$observed[df$cases == 0] == 0)) df <- df[df$cases != 0, ]
  df$outside <- ppc_observed_outside_cri(df)
  ppc_filter_hh_size(df[order(df$hh_size, df$cases), ],
                     min_hh_size, max_hh_size)
}

#' @title Standard epidemiological prior or posterior predictive check
#' @description
#' Prior or posterior predictive check of standard epidemiological measures
#' by household size, observed against predicted:
#' cases distribution, attack rate and secondary attack rate.
#' On a case-ascertained model, such as [template_case_ascertained],
#' households enter the study through a case,
#' and the attack rate is not checked.
#'
#' Each checked measure comes with a method formatting it as a `data.frame`
#' by household size, restricted to the household size bounds.
#' Columns `observed`, `predicted`, `cri_lower` and `cri_upper` hold
#' the observed value, the predicted point estimate and credible interval.
#' Column `outside` flags observed values beyond that interval.
#' \describe{
#' \item{`format_attack_rate()`}{attack rates, one row per household size,
#' with `counts` the number of households of each size;}
#' \item{`format_sec_attack_rate()`}{secondary attack rates,
#' one row per household size,
#' with `counts` the number of households with a case;}
#' \item{`format_cases_distr()`}{number of households with each number of
#' `cases`, one row per household size and number of cases.}
#' }
#' @examples
#' hh_model <- flu_1918_frederick |>
#'   new_housepi_model("final-size") |>
#'   add_parameter(import = "template")
#' hh_fit <- run_inference(hh_model)
#'
#' ppc <- check_predictive(hh_model, hh_fit, nrep = 50)
#' ppc$format_attack_rate()
#' @docType class
HousepiEpiPredictiveCheck <- R6::R6Class(
  classname = "HousepiEpiPredictiveCheck",
  inherit = HousepiPredictiveCheck,
  cloneable = TRUE,
  lock_objects = FALSE,

  ####  Private  ####

  private = list(
    min_hh_size_ = NA_integer_,
    max_hh_size_ = NA_integer_,
    house_counts_ = NULL,
    measure_topics_ = c(attack_rate = "measure_ar_by_hh_size",
                        sec_attack_rate = "measure_sar_by_hh_size",
                        cases_distr = "measure_cases_by_hh_size"),
    # body of method `format_<measure>`, evaluated in the object
    formats_ = list(
      attack_rate = quote(ppc_format_rate(
        self$attack_rate, "ar", private$house_counts_[c("hh_size", "total")],
        private$min_hh_size_, private$max_hh_size_
      )),
      sec_attack_rate = quote(ppc_format_rate(
        self$sec_attack_rate, "sar",
        private$house_counts_[c("hh_size", "at_least_1_case")],
        private$min_hh_size_, private$max_hh_size_
      )),
      cases_distr = quote(ppc_format_cases_distr(
        self$cases_distr, private$min_hh_size_, private$max_hh_size_
      ))
    )
  ),

  ####  Public  ####

  public = list(
    #' @param hh_model `r rd_model_in_`
    #' @param draws `r rd_draws_`
    #' @param ... `r rd_measure_dots_`
    #' They are checked alongside standard epidemiological measures.
    #' @param long_measures `r rd_long_measures_`
    #' @param nrep `r rd_nrep_`
    #' `r rd_nrep_prior_` `NA` is then an error.
    #' @param thin_ess `r rd_thin_ess_`
    #' @param min_hh_size `r rd_min_hh_size_`
    #' @param max_hh_size `r rd_max_hh_size_`
    #' @param threads `r rd_thread_`
    #' @param progress `r rd_progress_`
    initialize = function(hh_model,
                          draws = NULL,
                          ...,
                          long_measures = NULL,
                          nrep = NA_integer_,
                          thin_ess = TRUE,
                          min_hh_size = NA_integer_,
                          max_hh_size = NA_integer_,
                          threads = 1L,
                          progress = TRUE) {
      clash <- intersect(c(names(list(...)), names(long_measures)),
                         names(private$measure_topics_))
      if (length(clash)) stop_h(
        paste("{cli::qty(clash)}Measure name{?s} {.field {clash}} already used",
              "by a standard epidemiological measure."),
        "Rename the measure or set {.code epi = FALSE}."
      )
      # built-in measures read only hh_id and status:
      # * no dataset reporting process, no custom measure -> frames stay
      #   engine-coded, the relabelling skipped once per replicate.
      esc <- hh_model$escaped
      private$labelled_ <- !is.null(get_reporting_process(hh_model)$dataset) ||
        length(measure_dots_(...)$measures) || length(long_measures)

      # case-ascertained model -> attack rate not checked
      standard <- list(attack_rate = measure_ar_by_hh_size,
                       sec_attack_rate = measure_sar_by_hh_size,
                       cases_distr = measure_cases_by_hh_size)
      mcmc_ptr <- hh_model$.__enclos_env__$private$mcmc_ptr_
      if (.engine$has_case_ascertainment(mcmc_ptr)) standard$attack_rate <- NULL

      # formatter of each checked measure, enclosed in the object:
      # a clone rebinds it to its own fields
      for (m in names(standard)) self[[paste0("format_", m)]] <-
        rlang::new_function(list(), private$formats_[[m]],
                            self$.__enclos_env__)

      # initialize parent class HousepiPredictiveCheck
      rlang::inject(super$initialize(
        hh_model = hh_model, draws = draws, !!!standard, ...,
        long_measures = long_measures, nrep = nrep, thin_ess = thin_ess,
        threads = threads, progress = progress
      ))
      private$house_counts_ <- helper_house_counts(get_dataset(hh_model), esc)
      self$set_min_hh_size(min_hh_size)
      self$set_max_hh_size(max_hh_size)

      invisible(self)
    },

    #' @description Set the minimum household size of the formatted outputs.
    #' @param n `r rd_min_hh_size_`
    set_min_hh_size = function(n) {
      private$min_hh_size_ <- if (is.null(n)) NA_integer_ else n
      invisible(self)
    },

    #' @description Set the maximum household size of the formatted outputs.
    #' @param n `r rd_max_hh_size_`
    set_max_hh_size = function(n) {
      private$max_hh_size_ <- if (is.null(n)) NA_integer_ else n
      invisible(self)
    }
  )
)

#' @title Prior or posterior predictive check of a housepi model
#' @description
#' Predictive check comparing observed values of household measures
#' against their values on household studies simulated from parameter values
#' \insertCite{rubin_bayesianly_1984}{housepi}.
#' Without `draws`, the parameters are drawn from their prior,
#' as in [run_simulation()], for a prior predictive check.
#' With posterior draws, it is a posterior predictive check.
#' Draws of [run_inference()] on another model or other data
#' make a predictive check at given values, not a posterior one.
#' It checks standard epidemiological measures by household size:
#' \describe{
#' \item{`"attack_rate"`}{see [measure_ar_by_hh_size],
#' not checked on a case-ascertained model;}
#' \item{`"sec_attack_rate"`}{see [measure_sar_by_hh_size];}
#' \item{`"cases_distr"`}{see [measure_cases_by_hh_size].}
#' }
#' Named measure functions in `...` are checked
#' alongside them, under names that differ from these.
#' @param hh_model `r rd_model_in_`
#' @param draws `r rd_draws_`
#' @param ... `r rd_measure_dots_`
#' @param nrep `r rd_nrep_`
#' `r rd_nrep_prior_` `NA` is then an error.
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
#' hh_fit <- run_inference(hh_model)
#'
#' ppc <- check_predictive(hh_model, hh_fit, nrep = 50)
#' ppc
#'
#' # prior predictive check, before any inference
#' ppc <- check_predictive(hh_model, nrep = 50)
#' ppc
#'
#' # a user-defined measure: the overall attack rate
#' overall_ar <- function(data) mean(data$status > 0)
#' ppc <- check_predictive(hh_model, hh_fit, overall = overall_ar,
#'                         nrep = 50)
#' ppc
#' ppc$overall$observed
#'
#' # only the user-defined measures
#' ppc <- check_predictive(hh_model, hh_fit, overall = overall_ar,
#'                         epi = FALSE, nrep = 50)
#' ppc
#' @export
check_predictive <- function(hh_model,
                             draws = NULL,
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
  if (!length(measure_dots_(...)$measures)) stop_h(
    "No measure to check.",
    "Pass measure functions in {.arg ...}, or keep {.code epi = TRUE}."
  )
  if (!all(is.na(c(min_hh_size, max_hh_size)))) stop_h(
    paste("Arguments {.arg min_hh_size} and {.arg max_hh_size} bound only",
          "standard epidemiological measures."),
    "Set {.code epi = TRUE} or restrict household sizes inside the measures."
  )
  HousepiPredictiveCheck$new(
    hh_model$copy(), draws, ..., nrep = nrep, thin_ess = thin_ess,
    threads = threads, progress = progress
  )
}
