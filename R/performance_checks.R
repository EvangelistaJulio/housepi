## housepi - Epidemiological Analysis of Household Outbreaks.
## Copyright 2024-2026 Institut Pasteur (Paris), France.
## SPDX-License-Identifier: GPL-3.0-or-later
##

# Roxygen text shared by several blocks, read at the use site as `r <name>`.

rd_perf_houses_ <- "`integer()` One or more household counts to sweep.
If `NA`, the number of households carried by the attached data."
rd_perf_truth_ <- "Parameter values every study is simulated at,
one per estimated parameter, keyed by parameter name,
as in `c(alpha = 0.1, beta = 0.2)` or a one-row `data.frame`.
Posterior draws, from [run_inference()] or as `posterior::draws`,
are collapsed to their point estimate."

########  Performance  ########

#' Simulation values: `truth` per declared parameter.
#' * posterior draws -> collapsed to their point estimate;
#' * other table holding several values for a parameter -> error;
#' * name not declared -> warning, value dropped;
#' * `NA` value -> dropped;
#' * estimated parameter without a value -> error,
#'   as a performance check simulates every replicate at the same known values.
#' @param estimated `character()` Parameters the check is about.
#' @noRd
performance_values_ <- function(hh_model, truth, point_estimate, estimated) {
  draws <- inherits(truth, "HousepiPosteriorDraws") ||
    posterior::is_draws(truth)
  several <- if (!draws) names(truth)[lengths(truth) > 1L]
  if (length(several)) stop_h(
    "{.arg truth} holds several values for parameter{?s} {.val {several}}.",
    "A performance check simulates every replicate at the same values.",
    "To check inference across posterior draws, use {.fn check_calibration}."
  )
  values <- if (!is.null(truth)) {
    hh_model$.__enclos_env__$private$read_param_table_(
      truth, thin_ess = FALSE, arg = "truth"
    ) |> performance_collapse_truth_(point_estimate)
  }
  declared <- get_parameters(hh_model)[["name"]]
  unknown <- setdiff(names(values), declared)
  if (length(unknown)) warn_h(
    "No parameter named {.val {unknown}} in the model.",
    `!` = "Its {.arg truth} value is ignored.",
    find_closest_match(
      unknown, declared,
      otherwise = "Show the model parameters with {.fn get_parameters}."
    )
  )
  values <- values[intersect(names(values), declared)]
  values <- values[!is.na(values)]
  unpinned <- setdiff(estimated, names(values))
  if (length(unpinned)) stop_h(
    "No simulation value for parameter{?s} {.val {unpinned}}.",
    "A performance check simulates at known values: pass them as {.arg truth}."
  )
  values
}

#' Collapse a table of draws to one named vector of point estimates.
#' @noRd
performance_collapse_truth_ <- function(table, point_estimate) {
  stat <- if (identical(point_estimate, "mean")) mean else stats::median
  vapply(table, function(x) stat(as.numeric(x)), numeric(1L))
}

#' Performance measures per household count and parameter,
#' each with its Monte Carlo standard error (Morris et al. 2019, table 6).
#' * one replicate -> every Monte Carlo standard error `NA`;
#' * zero truth -> relative bias `NA`.
#' @param estimates `data.frame` One row per replicate and parameter.
#' @noRd
performance_summary_ <- function(estimates) {
  sweep_by_(estimates, c("houses", "variable"), function(d) {
    n <- nrow(d)
    truth <- d[["truth"]][[1L]]
    scale <- if (truth == 0) NA_real_ else truth
    error <- d[["estimate"]] - truth
    sq_error <- error^2
    width <- d[["cri_upper"]] - d[["cri_lower"]]
    coverage <- mean(d[["covers"]])
    # standard error of a mean over the replicates
    mcse <- \(x) stats::sd(x) / sqrt(n)
    data.frame(
      d[1L, c("houses", "variable")], truth = truth,
      estimate = mean(d[["estimate"]]),
      bias = mean(error), bias_mcse = mcse(error),
      rel_bias = mean(error) / scale, rel_bias_mcse = mcse(error) / abs(scale),
      mse = mean(sq_error), mse_mcse = mcse(sq_error),
      coverage = coverage,
      coverage_mcse =
        if (n > 1L) sqrt(coverage * (1 - coverage) / n) else NA_real_,
      cri_width = mean(width), cri_width_mcse = mcse(width),
      ess_bulk = mean(d[["ess_bulk"]]),
      nrep = n, n_flagged = sum(d[["flagged"]])
    )
  })
}

#' @title Resample the household structure of a study
#' @description
#' Draw whole households with replacement, so the household-size distribution
#' and the within-household composition of the covariates are preserved
#' while the number of households is scaled.
#' Results then describe a population of households
#' shaped exactly like those of `data`.
#' @param data `r rd_dataset_`
#' @param houses `integer(1)` Number of households to draw.
#' If `NA`, the number of households carried by `data`.
#' @return `data.frame` The resampled households,
#' `hh_id` renumbered sequentially,
#' so a household drawn twice stands as two independent households,
#' and `id` regenerated as a row number when the input carries one.
#' @examples
#' structure <- sample_house_structure(flu_1918_frederick, houses = 20)
#' new_housepi_model(structure, "final_size") |>
#'   summarise_house_cases()
#' @export
sample_house_structure <- function(data, houses = NA_integer_) {
  data <- as.data.frame(data)
  if (!"hh_id" %in% names(data)) stop_h(
    "Column {.field hh_id} is missing from the data.",
    "Provide an household identifier in order to allow sampling."
  )
  blocks <- split(data, data[["hh_id"]], drop = TRUE)
  houses <- sweep_house_counts_(houses, length(blocks), one = TRUE)
  drawn <- blocks[sample.int(length(blocks), houses, replace = TRUE)]
  out <- do.call(rbind, drawn)
  rownames(out) <- NULL
  out[["hh_id"]] <- rep(seq_along(drawn), vapply(drawn, nrow, integer(1L)))
  if ("id" %in% names(out)) out[["id"]] <- seq_len(nrow(out))
  out
}

#' @title Performance check of a housepi model
#' @description
#' Container for the sweep run by [check_performance()]:
#' household studies simulated at known parameter values and fitted back,
#' replicated at one or more household counts.
#' @examples
#' performance <- flu_1918_frederick |>
#'   new_housepi_model("final-size") |>
#'   add_parameter(import = "template") |>
#'   check_performance(c(alpha = 0.15, beta = 0.2, gamma = 1),
#'                     nrep = 2, houses = 30)
#' performance
#' performance$summarise()
#' @docType class
HousepiPerformanceCheck <- R6::R6Class(
  classname = "HousepiPerformanceCheck",

  ####  Private  ####

  private = list(
    name_ = NA_character_,
    variables_ = NULL,
    truth_ = NULL,
    point_estimate_ = NA_character_,
    cri_level_ = NA_real_,
    # `HousepiSweep` Replicates fitted.
    sweep_ = NULL,

    # Reduce one fit to a row per estimated parameter: point estimate,
    # posterior sd, credible interval, whether it covers `truth`, diagnostics.
    reduce_ = function(fit, truth) {
      variables <- private$variables_
      pdraws <- get_parameter_draws(fit)
      summ <- sweep_point_summary_(
        lapply(stats::setNames(nm = variables), \(v) as.numeric(pdraws[[v]])),
        private$point_estimate_, private$cri_level_
      )
      truth <- unname(truth[variables])
      covers <- summ[["cri_lower"]] <= truth & truth <= summ[["cri_upper"]]
      cbind(summ["variable"], truth = truth,
            summ[c("estimate", "sd", "cri_lower", "cri_upper")],
            covers = covers,
            sweep_diagnostics_(fit, variables),
            row.names = NULL)
    }
  ),

  ####  Active  ####

  active = list(
    #' @field estimates `data.frame`
    #' One row per household count, replicate and parameter
    #' \describe{
    #'   \item{truth}{Value the study was simulated at;}
    #'   \item{estimate}{Posterior point estimate;}
    #'   \item{sd}{Posterior standard deviation;}
    #'   \item{cri_lower, cri_upper}{Credible interval bounds;}
    #'   \item{covers}{Whether the credible interval covers the truth;}
    #'   \item{rhat, ess_bulk}{Convergence diagnostics of the fit;}
    #'   \item{flagged}{`TRUE` when the diagnostics flagged the parameter:
    #'     `r paste(diagnostic_flags_, collapse = "; ")`.}
    #' }
    estimates = function(value) {
      assert_read_only_(value, "estimates")
      private$sweep_$estimates
    },

    #' @field replicates `data.frame`
    #' One row per household count and replicate
    #' \describe{
    #'   \item{seed}{Seed of the fit's random number generators;}
    #'   \item{inference_secs}{Wall-clock of the sampling loop;}
    #'   \item{min_ess_bulk}{Smallest bulk effective sample size
    #'     over the parameters of the fit;}
    #'   \item{flagged}{Whether any parameter of the fit was flagged.}
    #' }
    replicates = function(value) {
      assert_read_only_(value, "replicates")
      private$sweep_$replicates
    },

    #' @field overview `data.frame`
    #' Sweep cost per household count, averaged over the replicates
    #' \describe{
    #'   \item{nrep, n_flagged}{Number of replicates run, and of those flagged;}
    #'   \item{mean_inference_secs}{Mean wall-clock of a sampling loop;}
    #'   \item{secs_per_ess}{Seconds of sampling per effective draw
    #'     of the worst-mixing parameter, pooled over the replicates:
    #'     a cost of this run, machine and thread count, not of the model.}
    #' }
    overview = function(value) {
      assert_read_only_(value, "overview")
      private$sweep_$overview
    },

    #' @field truth named `numeric()`
    #' Values the studies were simulated at.
    truth = function(value) {
      assert_read_only_(value, "truth")
      private$truth_
    },

    #' @field houses `integer()` Household counts swept.
    houses = function(value) {
      assert_read_only_(value, "houses")
      private$sweep_$houses
    },

    #' @field nrep `integer()`
    #' Number of replicates at each household count, in the order of `houses`.
    nrep = function(value) {
      assert_read_only_(value, "nrep")
      private$sweep_$nrep
    },

    #' @field point_estimate `character(1)`
    #' Summary statistic taken as the posterior point estimate.
    point_estimate = function(value) {
      assert_read_only_(value, "point_estimate")
      private$point_estimate_
    },

    #' @field cri_level `numeric(1)`
    #' Credible level the intervals and the coverage are computed at.
    cri_level = function(value) {
      assert_read_only_(value, "cri_level")
      private$cri_level_
    },

    #' @field warnings `r rd_warnings_`
    #' They pool the replicates of the sweep.
    warnings = function(value) {
      assert_read_only_(value, "warnings")
      private$sweep_$warnings
    }
  ),

  ####  Public  ####

  public = list(
    #' @description
    #' Create the performance check object.
    #' Use [check_performance()] instead of calling `$new()`.
    #' @param hh_model `r rd_model_in_`
    #' @param truth `r rd_perf_truth_`
    #' @param nrep `r rd_sweep_nrep_`
    #' @param houses `r rd_perf_houses_`
    #' @param structure_sampler `r rd_sweep_sampler_`
    #' @param steps `r rd_smpl_steps_`
    #' @param thin `r rd_smpl_thin_`
    #' @param adapt `r rd_smpl_adapt_`
    #' @param chains `r rd_smpl_chain_`
    #' @param threads `r rd_sweep_threads_`
    #' @param progress `r rd_progress_`
    #' @param ... `r rd_dots_config_`
    initialize = function(hh_model, truth = NULL,
                          nrep = 50L, houses = NA_integer_,
                          structure_sampler = sample_house_structure,
                          steps = 2000L, thin = 10L, adapt = 1000L,
                          chains = 4L, threads = 1L,
                          progress = TRUE, ...) {
      nrep <- assert_count_(nrep, "Number of replicates")
      private$point_estimate_ <- get_global_point_estimate()
      private$cri_level_ <- get_global_cri_level()
      private$name_ <- hh_model$.__enclos_env__$private$compiled_$name

      private$variables_ <- sweep_estimated_(hh_model)
      values <- performance_values_(hh_model, truth, private$point_estimate_,
                                    private$variables_)
      private$truth_ <- values[private$variables_]
      private$sweep_ <- HousepiSweep$new(
        hh_model, sampler = structure_sampler, houses = houses, nrep = nrep,
        values = sweep_param_table_(values, nrep),
        fit_args = list(steps = steps, thin = thin, adapt = adapt,
                        chains = chains),
        threads = threads, dots = list(...), reducer = private$reduce_,
        progress = progress, label = "Performance replicates",
        what = "Performance check"
      )
      invisible(self)
    },

    #' @description Compute the performance measures
    #' per household count and parameter.
    #' @return `data.frame`
    #' \describe{
    #'   \item{bias, rel_bias}{Mean error, absolute and relative to the truth;}
    #'   \item{mse}{Mean squared error;}
    #'   \item{coverage}{Share of replicates whose interval covers the truth;}
    #'   \item{cri_width}{Mean width of the credible interval;}
    #'   \item{bias_mcse, rel_bias_mcse, mse_mcse, coverage_mcse,
    #'     cri_width_mcse}{Monte Carlo standard error of each measure,
    #'     `NA` for a single replicate;}
    #'   \item{nrep, n_flagged}{Number of replicates aggregated,
    #'     and of those flagged.}
    #' }
    summarise = function() performance_summary_(private$sweep_$estimates),

    #' @description Same as `$summarise()`.
    #' @return See `$summarise()`.
    summarize = function() self$summarise(),

    #' @description Print a summary of the sweep.
    print = function() {
      model <- spaced_choice(private$name_)
      private$sweep_$print_rows(
        "Performance check of model {.val {model}}.",
        c(`Point estimate` = "{private$point_estimate_}",
          `Credible interval` = "{round(100 * private$cri_level_)}%")
      )
      invisible(self)
    }
  )
)

#' @title Check the performance of inference at known parameter values
#' @description
#' Simulate outbreaks in the households of the attached data
#' at known parameter values, fit each resulting household study back
#' with the same model, and measure the performance of the posterior
#' against the values that generated it
#' \insertCite{morris_using_2019}{housepi}.
#' `r rd_sweep_structure_`
#' The sweep is repeated at one or more household counts,
#' which shows whether a bias shrinks with the size of the study.
#'
#' Every study is simulated at `truth`,
#' and the priors of the same model drive inference,
#' so no separate simulation model is needed.
#' Repeated sampling at one value measures frequentist properties
#' of the posterior.
#' To check inference at values drawn anew for each replicate,
#' see [check_calibration()].
#' @param hh_model `r rd_model_in_`
#' @param truth `r rd_perf_truth_`
#' @param nrep `r rd_sweep_nrep_`
#' @param houses `r rd_perf_houses_`
#' @param structure_sampler `r rd_sweep_sampler_`
#' @param steps `r rd_smpl_steps_`
#' @param thin `r rd_smpl_thin_`
#' @param adapt `r rd_smpl_adapt_`
#' @param chains `r rd_smpl_chain_`
#' @param threads `r rd_sweep_threads_`
#' @param progress `r rd_progress_`
#' @param ... `r rd_dots_config_`
#' @return An object of class [HousepiPerformanceCheck].
#' @details
#' The first replicate that fails stops the check,
#' with an error naming its household count and replicate.
#' `r rd_sweep_reporting_`
#' `r rd_global_summary_` `r rd_sweep_seed_`
#' @examples
#' performance <- flu_1918_frederick |>
#'   new_housepi_model("final-size") |>
#'   add_parameter(import = "template") |>
#'   check_performance(c(alpha = 0.15, beta = 0.2, gamma = 1),
#'                     nrep = 2, houses = 30)
#' performance$summarise()
#' performance$overview
#' @references
#' \insertAllCited{}
#' @export
check_performance <- function(hh_model,
                              truth = NULL,
                              nrep = 50L,
                              houses = NA_integer_,
                              structure_sampler = sample_house_structure,
                              steps = 2000L, thin = 10L, adapt = 1000L,
                              chains = 4L, threads = 1L,
                              progress = TRUE,
                              ...) {
  stopifnothousepi(hh_model)
  assert_advanced_settings_("time_step", ...)
  HousepiPerformanceCheck$new(
    hh_model$copy(), truth = truth, nrep = nrep, houses = houses,
    structure_sampler = structure_sampler,
    steps = steps, thin = thin, adapt = adapt,
    chains = chains, threads = threads, progress = progress, ...
  )
}
