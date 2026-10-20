## housepi - Epidemiological Analysis of Household Outbreaks.
## Copyright 2024-2026 Institut Pasteur (Paris), France.
## SPDX-License-Identifier: GPL-3.0-or-later
##

# Roxygen text shared by several blocks, read at the use site as `r <name>`.

rd_recov_houses_ <- "`integer()` One or more household counts to sweep.
If `NA`, the number of households carried by the attached data."
rd_recov_draws_ <- "Values the surveys are simulated at,
one for every estimated parameter.
Several draws are collapsed to their point estimate."

########  Recovery  ########

#' Simulation values: the point estimate of `draws`, per declared parameter.
#' @noRd
recovery_pins_ <- function(hh_model, draws, point_estimate) {
  if (is.null(draws)) return(numeric())
  given <- recovery_collapse_truth_(draws, point_estimate)
  hh_model$.__enclos_env__$private$assert_no_baseline_value_(given, "draws")
  declared <- get_parameters(hh_model)[["name"]]
  unknown <- setdiff(names(given), declared)
  if (length(unknown)) warn_h(
    "No parameter named {.val {unknown}} in the model.",
    `!` = "Its {.arg draws} value is ignored."
  )
  given <- given[intersect(names(given), declared)]
  given[!is.na(given)]
}

#' Collapse draws to one named vector of point estimates.
#' @noRd
recovery_collapse_truth_ <- function(draws, point_estimate) {
  tab <- helper_read_param_table_(draws, thin_ess = FALSE)
  stat <- if (identical(point_estimate, "mean")) mean else stats::median
  vapply(tab, function(x) stat(as.numeric(x)), numeric(1L))
}

#' Simulation values of the estimated parameters.
#' * estimated parameter without a value -> error,
#' as a recovery check simulates every replicate at the same known values.
#' @noRd
recovery_truth_ <- function(hh_model, pins) {
  estimated <- sweep_estimated_(hh_model)
  unpinned <- setdiff(estimated, names(pins))
  if (length(unpinned)) stop_h(
    "No simulation value for parameter{?s} {.val {unpinned}}.",
    "A recovery check simulates at known values: pass them as {.arg draws}."
  )
  pins[estimated]
}

#' Recovery statistics per household count and parameter.
#' @param estimates `data.frame` One row per replicate and parameter.
#' @noRd
recovery_summary_ <- function(estimates) {
  sweep_by_(estimates, c("houses", "variable"), function(d) {
    truth <- d[["truth"]]
    error <- d[["estimate"]] - truth
    data.frame(
      d[1L, c("houses", "variable")], truth = truth[[1L]],
      estimate = mean(d[["estimate"]]), bias = mean(error),
      rel_bias = if (any(truth == 0)) NA_real_ else mean(error / truth),
      rmse = sqrt(mean(error^2)),
      coverage = mean(d[["covers"]]),
      cri_width = mean(d[["cri_upper"]] - d[["cri_lower"]]),
      ess_bulk = mean(d[["ess_bulk"]]),
      nrep = nrow(d), n_flagged = sum(d[["flagged"]])
    )
  })
}

#' @title Resample the household structure of a survey
#' @description
#' Draw whole households with replacement, so the household-size distribution
#' and the within-household composition of the covariates are preserved
#' while the number of households is scaled.
#' @param data_wide `r rd_add_wide_`
#' @param houses `integer(1)` Number of households to draw.
#' If `NA`, the number of households carried by `data_wide`.
#' @return `data.frame` The survey, `hh_id` renumbered sequentially,
#' so a household drawn twice stands as two independent households,
#' and `id` regenerated as a row number when the input carries one.
#' @examples
#' structure <- sample_house_structure(flu_1918_frederick, houses = 20)
#' new_housepi_model(structure, "final_size") |>
#'   summarise_house_cases()
#' @export
sample_house_structure <- function(data_wide, houses = NA_integer_) {
  data <- as.data.frame(data_wide)
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

#' @title Parameter recovery check of a housepi model
#' @description
#' Container for the sweep run by [check_recovery()]:
#' surveys simulated at known parameter values and fitted back,
#' replicated at one or more household counts.
#' @examples
#' recovery <- flu_1918_frederick |>
#'   new_housepi_model("final-size") |>
#'   add_parameter(import = "template") |>
#'   check_recovery(data.frame(alpha = 0.15, beta = 0.2, gamma = 1),
#'                  nrep = 2, houses = 30, progress = FALSE)
#' recovery
#' recovery$summary
#' @docType class
HousepiRecoveryCheck <- R6::R6Class(
  classname = "HousepiRecoveryCheck",

  ####  Private  ####

  private = list(
    name_ = NA_character_,
    nrep_ = NA_integer_,
    houses_ = NULL,
    truth_ = NULL,
    point_estimate_ = NA_character_,
    cri_level_ = NA_real_,
    estimates_ = NULL,
    replicates_ = NULL,
    warnings_ = NULL,
    secs_ = NA_real_,

    # Summarise the draws per parameter: point estimate, posterior sd,
    # credible interval, whether it covers `truth`, convergence diagnostics.
    summarise_draws_ = function(draws, truth) {
      stat <- switch(private$point_estimate_, mean = mean, stats::median)
      probs <- cri_level_probs(private$cri_level_)
      pdraws <- get_parameter_draws(draws)
      summ <- as.data.frame(posterior::summarise_draws(
        pdraws,
        estimate = stat,
        sd = stats::sd,
        cri_lower = \(x) unname(stats::quantile(x, probs[["lower"]])),
        cri_upper = \(x) unname(stats::quantile(x, probs[["upper"]]))
      ))
      summ <- summ[summ[["variable"]] %in% names(private$truth_), ]
      truth <- unname(truth[summ[["variable"]]])
      covers <- summ[["cri_lower"]] <= truth & truth <= summ[["cri_upper"]]
      cbind(summ["variable"], truth = truth,
            summ[c("estimate", "sd", "cri_lower", "cri_upper")],
            covers = covers,
            sweep_diagnostics_(draws, summ[["variable"]]),
            row.names = NULL)
    }
  ),

  ####  Active  ####

  active = list(
    #' @field estimates `data.frame`
    #' One row per household count, replicate and parameter
    #' \describe{
    #'   \item{truth}{Value the survey was simulated at;}
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
      private$estimates_
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
      private$replicates_
    },

    #' @field summary `data.frame`
    #' Recovery statistics per household count and parameter
    #' \describe{
    #'   \item{bias, rel_bias}{Mean error, absolute and relative to the truth;}
    #'   \item{rmse}{Root mean squared error;}
    #'   \item{coverage}{Share of replicates whose interval covers the truth;}
    #'   \item{cri_width}{Mean width of the credible interval;}
    #'   \item{nrep, n_flagged}{Replicates aggregated, and those flagged.}
    #' }
    summary = function(value) {
      assert_read_only_(value, "summary")
      recovery_summary_(private$estimates_)
    },

    #' @field overview `data.frame`
    #' Sweep cost per household count, averaged over the replicates
    #' \describe{
    #'   \item{nrep, n_flagged}{Replicates run, and those flagged;}
    #'   \item{mean_inference_secs}{Mean wall-clock of a sampling loop;}
    #'   \item{secs_per_ess}{Seconds of sampling per effective draw
    #'     of the worst-mixing parameter, pooled over the replicates:
    #'     a cost of this run, machine and thread count, not of the model.}
    #' }
    overview = function(value) {
      assert_read_only_(value, "overview")
      sweep_overview_(private$replicates_)
    },

    #' @field truth named `numeric()`
    #' Values the surveys were simulated at.
    truth = function(value) {
      assert_read_only_(value, "truth")
      private$truth_
    },

    #' @field houses `integer()` Household counts swept.
    houses = function(value) {
      assert_read_only_(value, "houses")
      private$houses_
    },

    #' @field nrep `integer(1)` Replicates run at each household count.
    nrep = function(value) {
      assert_read_only_(value, "nrep")
      private$nrep_
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

    #' @field warnings `data.frame`
    #' Distinct warnings raised over the sweep, with their counts.
    warnings = function(value) {
      assert_read_only_(value, "warnings")
      private$warnings_
    }
  ),

  ####  Public  ####

  public = list(
    #' @description
    #' Create the recovery check object.
    #' Use [check_recovery()] instead of calling `$new()`.
    #' @param hh_model `r rd_model_in_`
    #' @param draws `r rd_draws_`
    #' `r rd_recov_draws_`
    #' @param nrep `r rd_sweep_nrep_`
    #' @param houses `r rd_recov_houses_`
    #' @param structure_sampler `r rd_sweep_sampler_`
    #' @param steps `r rd_smpl_steps_`
    #' @param burn `r rd_smpl_burn_`
    #' @param thin `r rd_smpl_thin_`
    #' @param adapt `r rd_smpl_adapt_`
    #' @param chains `r rd_smpl_chain_`
    #' @param threads `r rd_sweep_threads_`
    #' @param progress `r rd_progress_`
    #' @param ... `r rd_dots_config_`
    initialize = function(hh_model, draws = NULL,
                          nrep = 50L, houses = NA_integer_,
                          structure_sampler = sample_house_structure,
                          steps = 1000L, burn = 0L,
                          thin = 5L, adapt = 500L,
                          chains = 4L, threads = 1L,
                          progress = TRUE, ...) {
      private$nrep_ <- assert_count_(nrep, "Number of replicates")
      private$point_estimate_ <- get_global_point_estimate()
      private$cri_level_ <- get_global_cri_level()
      private$name_ <- hh_model$.__enclos_env__$private$gen_$raw_name

      pins <- recovery_pins_(hh_model, draws, private$point_estimate_)
      private$truth_ <- recovery_truth_(hh_model, pins)
      structure <- sweep_base_structure_(hh_model)
      private$houses_ <- sweep_houses_(structure, houses,
                                       is.null(structure_sampler))

      fit_args <- list(steps = steps, burn = burn, thin = thin, adapt = adapt,
                       chains = chains)
      swept <- sweep_run_(
        hh_model, structure, structure_sampler, private$houses_, private$nrep_,
        pins, fit_args, threads, list(...), private$summarise_draws_,
        progress, "Recovery replicates"
      )
      private$estimates_ <- swept$estimates
      private$replicates_ <- swept$replicates
      private$warnings_ <- swept$warnings
      private$secs_ <- swept$secs
      sweep_report_warnings_(private$warnings_, "Parameter recovery")
      invisible(self)
    },

    #' @description Print a summary of the sweep.
    print = function() {
      model <- spaced_choice(private$name_)
      flagged <- sum(private$replicates_[["flagged"]])
      fit_secs <- mean(private$replicates_[["inference_secs"]])
      print_rows_("Parameter recovery check of model {.val {model}}.", c(
        Households = "{private$houses_}",
        `Replicates each` = "{private$nrep_}",
        `Point estimate` = "{private$point_estimate_}",
        `Credible interval` = "{round(100 * private$cri_level_)}%",
        `Flagged by diagnostics` = "{flagged}",
        `Run time` = "{format_secs_(private$secs_)}",
        `Inference time per fit` = "{format_secs_(fit_secs)}"
      ))
      invisible(self)
    }
  )
)

#' @title Check that inference recovers known parameters
#' @description
#' Simulate surveys at known parameter values, fit each one back with the same
#' model, and report how well the posterior recovers the values that generated
#' it.
#' The sweep is repeated at one or more household counts,
#' which shows whether a bias shrinks with the size of the study.
#'
#' The surveys are simulated at `draws`,
#' and the priors of the same model drive inference,
#' so no separate simulation model is needed.
#' @param hh_model `r rd_model_in_`
#' @param draws `r rd_draws_`
#' `r rd_recov_draws_`
#' @param nrep `r rd_sweep_nrep_`
#' @param houses `r rd_recov_houses_`
#' @param structure_sampler `r rd_sweep_sampler_`
#' @param steps `r rd_smpl_steps_`
#' @param burn `r rd_smpl_burn_`
#' @param thin `r rd_smpl_thin_`
#' @param adapt `r rd_smpl_adapt_`
#' @param chains `r rd_smpl_chain_`
#' @param threads `r rd_sweep_threads_`
#' @param progress `r rd_progress_`
#' @param ... `r rd_dots_config_`
#' @return An object of class [HousepiRecoveryCheck].
#' @details
#' The first replicate that fails stops the check,
#' with an error naming its household count and replicate.
#' `r rd_sweep_reporting_`
#' `r rd_global_summary_` `r rd_sweep_seed_`
#' @examples
#' recovery <- flu_1918_frederick |>
#'   new_housepi_model("final-size") |>
#'   add_parameter(import = "template") |>
#'   check_recovery(data.frame(alpha = 0.15, beta = 0.2, gamma = 1),
#'                  nrep = 2, houses = 30, progress = FALSE)
#' recovery$summary
#' recovery$overview
#' @export
check_recovery <- function(hh_model,
                           draws = NULL,
                           nrep = 50L,
                           houses = NA_integer_,
                           structure_sampler = sample_house_structure,
                           steps = 1000L, burn = 0L,
                           thin = 5L, adapt = 500L,
                           chains = 4L, threads = 1L,
                           progress = TRUE,
                           ...) {
  stopifnothousepi(hh_model)
  assert_advanced_setup_(c("adapt_window", "time_step"), ...)
  HousepiRecoveryCheck$new(
    hh_model$copy(), draws = draws, nrep = nrep, houses = houses,
    structure_sampler = structure_sampler,
    steps = steps, burn = burn, thin = thin, adapt = adapt,
    chains = chains, threads = threads, progress = progress, ...
  )
}
