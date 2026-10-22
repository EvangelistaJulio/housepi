## housepi - Epidemiological Analysis of Household Outbreaks.
## Copyright 2024-2026 Institut Pasteur (Paris), France.
## SPDX-License-Identifier: GPL-3.0-or-later
##

# Roxygen text shared by several blocks, read at the use site as `r <name>`.

rd_calib_houses_ <- "`integer(1)` Household count of the surveys.
If `NA`, the number of households carried by the attached data."
rd_calib_rank_draws_ <- "`integer(1)`
Number of posterior draws each simulating value is ranked among,
thinned evenly from the stored draws of all chains."

########  Calibration  ########

#' Level of the uniformity test, simultaneous over the whole ECDF.
#' @noRd
calibration_alpha_ <- 0.05

#' Rank of `truth` among `size` draws thinned evenly from `x`,
#' ties broken uniformly at random.
#' A calibrated sampler gives ranks uniform on `0..size`.
#' @noRd
calibration_rank_ <- function(x, truth, size) {
  x <- helper_thin_even_(x, size)
  sum(x < truth) + sample.int(sum(x == truth) + 1L, 1L) - 1L
}

#' Ranks at or below each of `0..size - 1`,
#' `n` times the ECDF of the fractional ranks `(rank + 1) / (size + 1)`
#' at `z = k / (size + 1)`.
#' @noRd
calibration_counts_ <- function(ranks, size) {
  cumsum(tabulate(ranks + 1L, nbins = size + 1L))[seq_len(size)]
}

#' Logarithm of the statistic of Säilynoja et al. (2022),
#' twice the smallest binomial tail probability of the ECDF counts
#' over the evaluation points.
#' @noRd
calibration_log_gamma_ <- function(ranks, size) {
  n <- length(ranks)
  z <- seq_len(size) / (size + 1)
  counts <- calibration_counts_(ranks, size)
  log(2) + min(stats::pbinom(counts, n, z, log.p = TRUE),
               stats::pbinom(counts - 1, n, z, lower.tail = FALSE,
                             log.p = TRUE))
}

#' Bounds on the ECDF counts of `n` uniform ranks,
#' each holding with probability `1 - gamma` on its own.
#' @noRd
calibration_band_ <- function(n, size, gamma) {
  z <- seq_len(size) / (size + 1)
  list(lower = stats::qbinom(gamma / 2, n, z),
       upper = stats::qbinom(1 - gamma / 2, n, z))
}

#' Probability that the ECDF counts of `n` uniform ranks stay inside the band
#' at every evaluation point, Säilynoja et al. (2022), Sec. 2.3.
#' The counts form a Markov chain: of the ranks above `z[k - 1]`,
#' each falls at or below `z[k]` with probability `1 / (size + 2 - k)`.
#' @noRd
calibration_inside_ <- function(n, size, gamma) {
  band <- calibration_band_(n, size, gamma)
  states <- 0L
  prob <- 1
  for (k in seq_len(size)) {
    to <- seq(band$lower[[k]], band$upper[[k]])
    step <- outer(states, to, function(from, at) {
      stats::dbinom(at - from, n - from, 1 / (size + 2 - k))
    })
    prob <- drop(prob %*% step)
    states <- to
  }
  sum(prob)
}

#' Logarithm of the `alpha` quantile of the statistic under uniform ranks,
#' Modrák et al. (2025): the largest pointwise level whose band holds
#' simultaneously with probability at least `1 - alpha`.
#' Uniformity is rejected when the statistic falls below it,
#' so with probability at most `alpha` for uniform ranks.
#' The band probability steps down as the level rises:
#' bisection keeps the lower end on the side where the band holds,
#' starting from the union bound.
#' @noRd
calibration_log_threshold_ <- function(n, size, alpha) {
  holds <- \(lgamma) calibration_inside_(n, size, exp(lgamma)) >= 1 - alpha
  lower <- log(alpha / size)
  upper <- 0
  while (upper - lower > 1e-9) {
    mid <- (lower + upper) / 2
    if (holds(mid)) lower <- mid else upper <- mid
  }
  lower
}

#' Uniformity verdict per household count and parameter.
#' @param threshold `function` Log-threshold of a number of replicates.
#' @noRd
calibration_summary_ <- function(estimates, size, threshold) {
  sweep_by_(estimates, c("houses", "variable"), function(d) {
    log_gamma <- calibration_log_gamma_(d[["rank"]], size)
    log_threshold <- threshold(nrow(d))
    data.frame(
      d[1L, c("houses", "variable")], nrep = nrow(d),
      n_flagged = sum(d[["flagged"]]),
      log_gamma = log_gamma, log_gamma_threshold = log_threshold,
      uniform = log_gamma >= log_threshold
    )
  })
}

#' ECDF of the fractional ranks and its simultaneous band,
#' per household count and parameter.
#' @param threshold `function` Log-threshold of a number of replicates.
#' @noRd
calibration_ecdf_ <- function(estimates, size, threshold) {
  sweep_by_(estimates, c("houses", "variable"), function(d) {
    n <- nrow(d)
    band <- calibration_band_(n, size, exp(threshold(n)))
    data.frame(
      houses = d[["houses"]][[1L]], variable = d[["variable"]][[1L]],
      z = seq_len(size) / (size + 1),
      ecdf = calibration_counts_(d[["rank"]], size) / n,
      lower = band$lower / n, upper = band$upper / n
    )
  })
}

#' @title Simulation calibration check of a housepi model
#' @description
#' Container for the sweep run by [check_calibration()]:
#' surveys simulated at parameter values drawn from the prior and fitted back,
#' each value ranked among the posterior draws of its fit.
#' @examples
#' calibration <- flu_1918_frederick |>
#'   new_housepi_model("final-size") |>
#'   add_parameter(import = "template") |>
#'   check_calibration(nrep = 2, houses = 30, progress = FALSE)
#' calibration
#' calibration$summary
#' @docType class
HousepiCalibrationCheck <- R6::R6Class(
  classname = "HousepiCalibrationCheck",

  ####  Private  ####

  private = list(
    name_ = NA_character_,
    nrep_ = NA_integer_,
    houses_ = NULL,
    rank_draws_ = NA_integer_,
    variables_ = NULL,
    estimates_ = NULL,
    replicates_ = NULL,
    warnings_ = NULL,
    secs_ = NA_real_,
    thresholds_ = list(),

    # Rank the simulating value of every estimated parameter,
    # with the convergence diagnostics of the fit.
    # * draws too correlated for `rank_draws` independent ones -> flagged.
    summarise_ranks_ = function(draws, truth) {
      pdraws <- get_parameter_draws(draws)
      size <- private$rank_draws_
      variables <- private$variables_
      rank <- vapply(variables, function(v) {
        calibration_rank_(as.numeric(pdraws[[v]]), truth[[v]], size)
      }, integer(1L))
      diag <- sweep_diagnostics_(draws, variables)
      diag[["flagged"]] <- diag[["flagged"]] |
        !((diag[["ess_bulk"]] >= size) %in% TRUE)
      data.frame(variable = variables, truth = unname(truth[variables]),
                 rank = unname(rank), diag, row.names = NULL)
    },

    # Log-threshold of the uniformity test at `n` replicates, computed once.
    threshold_ = function(n) {
      key <- as.character(n)
      if (is.null(private$thresholds_[[key]]))
        private$thresholds_[[key]] <- calibration_log_threshold_(
          n, private$rank_draws_, calibration_alpha_
        )

      private$thresholds_[[key]]
    }
  ),

  ####  Active  ####

  active = list(
    #' @field estimates `data.frame`
    #' One row per replicate and parameter
    #' \describe{
    #'   \item{truth}{Value drawn from the prior the survey was simulated at;}
    #'   \item{rank}{Number of ranked posterior draws below the truth,
    #'     ties split at random;}
    #'   \item{rhat, ess_bulk}{Convergence diagnostics of the fit;}
    #'   \item{flagged}{`TRUE` when the diagnostics flagged the parameter:
    #'     `r paste(diagnostic_flags_, collapse = "; ")`,
    #'     or when its bulk effective sample size is below `rank_draws`.}
    #' }
    estimates = function(value) {
      assert_read_only_(value, "estimates")
      private$estimates_
    },

    #' @field replicates `data.frame`
    #' One row per replicate
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

    #' @field ecdf `data.frame`
    #' Empirical CDF of the fractional ranks `(rank + 1) / (rank_draws + 1)`
    #' per parameter, with its simultaneous band
    #' \describe{
    #'   \item{z}{Point the ECDF is evaluated at;}
    #'   \item{ecdf}{Share of the ranks at or below `z`;}
    #'   \item{lower, upper}{Band the ECDF of uniform ranks stays inside,
    #'     at every `z` at once, with probability at least 0.95.}
    #' }
    ecdf = function(value) {
      assert_read_only_(value, "ecdf")
      calibration_ecdf_(private$estimates_, private$rank_draws_,
                        private$threshold_)
    },

    #' @field summary `data.frame`
    #' Uniformity of the ranks per parameter
    #' \describe{
    #'   \item{nrep, n_flagged}{Replicates ranked, and those flagged;}
    #'   \item{log_gamma}{Logarithm of the uniformity statistic,
    #'     twice the smallest binomial tail probability of the ECDF;}
    #'   \item{log_gamma_threshold}{Value of `log_gamma` below which
    #'     uniformity is rejected, at level 0.05;}
    #'   \item{uniform}{Whether the ranks pass the test.}
    #' }
    summary = function(value) {
      assert_read_only_(value, "summary")
      calibration_summary_(private$estimates_, private$rank_draws_,
                           private$threshold_)
    },

    #' @field overview `data.frame`
    #' Sweep cost, averaged over the replicates
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

    #' @field houses `integer(1)` Household count of the surveys.
    houses = function(value) {
      assert_read_only_(value, "houses")
      private$houses_
    },

    #' @field nrep `integer(1)` Replicates run.
    nrep = function(value) {
      assert_read_only_(value, "nrep")
      private$nrep_
    },

    #' @field rank_draws `integer(1)`
    #' Posterior draws each simulating value is ranked among.
    rank_draws = function(value) {
      assert_read_only_(value, "rank_draws")
      private$rank_draws_
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
    #' Create the calibration check object.
    #' Use [check_calibration()] instead of calling `$new()`.
    #' @param hh_model `r rd_model_in_`
    #' @param nrep `r rd_sweep_nrep_`
    #' @param houses `r rd_calib_houses_`
    #' @param structure_sampler `r rd_sweep_sampler_`
    #' @param rank_draws `r rd_calib_rank_draws_`
    #' @param steps `r rd_smpl_steps_`
    #' @param burn `r rd_smpl_burn_`
    #' @param thin `r rd_smpl_thin_`
    #' @param adapt `r rd_smpl_adapt_`
    #' @param chains `r rd_smpl_chain_`
    #' @param threads `r rd_sweep_threads_`
    #' @param progress `r rd_progress_`
    #' @param ... `r rd_dots_config_`
    initialize = function(hh_model, nrep = 200L, houses = NA_integer_,
                          structure_sampler = sample_house_structure,
                          rank_draws = 99L, steps = 1000L, burn = 0L,
                          thin = 5L, adapt = 500L,
                          chains = 4L, threads = 1L,
                          progress = TRUE, ...) {
      private$nrep_ <- assert_count_(nrep, "Number of replicates")
      private$rank_draws_ <- assert_count_(rank_draws, "Number of ranked draws")
      stored <- resolve_setup_(steps, burn, thin, adapt, chains)$total
      if (stored < rank_draws) stop_h(
        "The fit stores {stored} draw{?s}, fewer than {.arg rank_draws}.",
        paste("Raise {.arg steps} or {.arg chains}, or lower {.arg thin} or",
              "{.arg rank_draws}.")
      )
      private$name_ <- hh_model$.__enclos_env__$private$compiled_$name
      private$variables_ <- sweep_estimated_(hh_model)
      structure <- sweep_base_structure_(hh_model)
      private$houses_ <- sweep_houses_(structure, houses,
                                       is.null(structure_sampler), one = TRUE)

      # no column: simulation draws every parameter anew from its prior
      pins <- numeric()
      fit_args <- list(steps = steps, burn = burn, thin = thin, adapt = adapt,
                       chains = chains)
      swept <- sweep_run_(
        hh_model, structure, structure_sampler, private$houses_, private$nrep_,
        pins, fit_args, threads, list(...), private$summarise_ranks_,
        progress, "Calibration replicates"
      )
      private$estimates_ <- swept$estimates
      private$replicates_ <- swept$replicates
      private$warnings_ <- swept$warnings
      private$secs_ <- swept$secs
      sweep_report_warnings_(private$warnings_, "Simulation calibration")
      invisible(self)
    },

    #' @description Print a summary of the sweep.
    print = function() {
      model <- spaced_choice(private$name_)
      flagged <- sum(private$replicates_[["flagged"]])
      fit_secs <- mean(private$replicates_[["inference_secs"]])
      print_rows_("Simulation calibration check of model {.val {model}}.", c(
        Households = "{private$houses_}",
        Replicates = "{private$nrep_}",
        `Posterior draws ranked among` = "{private$rank_draws_}",
        `Flagged by diagnostics` = "{flagged}",
        `Run time` = "{format_secs_(private$secs_)}",
        `Inference time per fit` = "{format_secs_(fit_secs)}"
      ))
      invisible(self)
    }
  )
)

#' @title Check that inference is calibrated against the prior
#' @description
#' Simulation-based calibration: draw parameter values from their prior,
#' simulate a survey at each draw, fit it back with the same model,
#' and rank the drawn value among the posterior draws of its fit.
#' When inference is correct, the ranks are uniform over the replicates.
#' Ranks piling at both ends point to a posterior too narrow,
#' in the middle to one too wide, on one side to a biased one.
#'
#' Every parameter carrying a prior and free to move is drawn anew
#' for each replicate, and its `value` is ignored.
#' Any other parameter is held at its value.
#' @param hh_model `r rd_model_in_`
#' @param nrep `r rd_sweep_nrep_`
#' @param houses `r rd_calib_houses_`
#' @param structure_sampler `r rd_sweep_sampler_`
#' @param rank_draws `r rd_calib_rank_draws_`
#' @param steps `r rd_smpl_steps_`
#' @param burn `r rd_smpl_burn_`
#' @param thin `r rd_smpl_thin_`
#' @param adapt `r rd_smpl_adapt_`
#' @param chains `r rd_smpl_chain_`
#' @param threads `r rd_sweep_threads_`
#' @param progress `r rd_progress_`
#' @param ... `r rd_dots_config_`
#' @return An object of class [HousepiCalibrationCheck].
#' @details
#' Ranks follow \insertCite{modrak_simulation_2025;textual}{housepi}:
#' the drawn value is ranked among `rank_draws` posterior draws thinned evenly
#' from the fit, and ties are split at random.
#' Correlated draws make the ranks pile at the ends,
#' even when inference is correct,
#' so a parameter whose bulk effective sample size
#' \insertCite{vehtari_rank_normalization_2021}{housepi} falls below
#' `rank_draws` is flagged.
#'
#' Uniformity is tested on the empirical CDF of the ranks, with the
#' simultaneous band and the statistic of
#' \insertCite{sailynoja_graphical_2022;textual}{housepi}.
#' Uniformity is rejected when the statistic falls below its 5th percentile
#' under uniform ranks \insertCite{modrak_simulation_2025}{housepi},
#' so a calibrated sampler is rejected with probability at most 0.05.
#' The check itself goes back to
#' \insertCite{talts_validating_2020;textual}{housepi}.
#'
#' `r rd_sweep_reporting_`
#'
#' `r rd_sweep_seed_`
#' @references
#' \insertAllCited{}
#' @examples
#' calibration <- flu_1918_frederick |>
#'   new_housepi_model("final-size") |>
#'   add_parameter(import = "template") |>
#'   check_calibration(nrep = 2, houses = 30, progress = FALSE)
#' calibration$summary
#' calibration$overview
#' @export
check_calibration <- function(hh_model,
                              nrep = 200L,
                              houses = NA_integer_,
                              structure_sampler = sample_house_structure,
                              rank_draws = 99L,
                              steps = 1000L, burn = 0L,
                              thin = 5L, adapt = 500L,
                              chains = 4L, threads = 1L,
                              progress = TRUE,
                              ...) {
  stopifnothousepi(hh_model)
  assert_advanced_setup_(c("adapt_window", "time_step"), ...)
  HousepiCalibrationCheck$new(
    hh_model$copy(), nrep = nrep, houses = houses,
    structure_sampler = structure_sampler, rank_draws = rank_draws,
    steps = steps, burn = burn, thin = thin, adapt = adapt,
    chains = chains, threads = threads, progress = progress, ...
  )
}
