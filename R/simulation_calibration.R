## housepi - Epidemiological Analysis of Household Outbreaks.
## Copyright 2024-2026 Institut Pasteur (Paris), France.
## SPDX-License-Identifier: GPL-3.0-or-later
##

# Roxygen text shared by several blocks, read at the use site as `r <name>`.

rd_calib_draws_ <- "`HousepiPosteriorDraws`
Posterior draws of the model on its attached data, from [run_inference()].
If `NULL`, the parameters are drawn from their prior."
rd_calib_houses_ <- "`integer(1)`
Number of households in each simulated study.
If `NA`, the number of households carried by the attached data.
With posterior `draws`, each study is fitted alongside the attached households."
rd_calib_rank_levels_ <- "`integer(1)`
Number of values a rank takes, from 0 to `rank_levels - 1`.
Each simulating value is ranked among `rank_levels - 1` posterior draws,
thinned evenly from the stored draws pooled over the chains."

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

#' Simulation values of a posterior check, `nrep` independent draws,
#' and the engine seeds of the runs they come from.
#' * `draws` not from run_inference() -> error;
#' * `draws` fitted to another model or data -> error;
#' * fewer independent draws than `nrep` -> error.
#' @return `list` Entries `values` and `seed`.
#' @noRd
calibration_posterior_ <- function(hh_model, draws, nrep) {
  if (!inherits(draws, "HousepiPosteriorDraws")) stop_h(
    "Argument {.arg draws} should come from {.fn run_inference}.",
    "Posterior calibration needs the draws of the model on its attached data."
  )
  pvt <- hh_model$.__enclos_env__$private
  fitted <- draws$.__enclos_env__$private
  # the wall build stamps the fingerprint the draws are compared with
  pvt$build_full_()
  differ <- differing_entries_(pvt$fingerprint_, fitted$fingerprint_)
  if (length(differ)) stop_h(
    "The draws were not fitted to this model and its data.",
    "They differ in {.field {differ}}.",
    "Pass the draws of {.fn run_inference} on this model."
  )
  values <- pvt$read_param_table_(draws, nrep, thin_ess = TRUE)
  # nolint next: object_usage_linter. read by the cli glue below
  n_draws <- nrow(values)
  if (n_draws < nrep) stop_h(
    "The draws hold {n_draws} independent draw{?s}, fewer than {.arg nrep}.",
    paste("Raise {.arg steps} or {.arg chains} of {.fn run_inference},",
          "or lower {.arg nrep}.")
  )
  rownames(values) <- NULL
  list(values = values, seed = fitted$seed_)
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
#' household studies simulated at parameter values drawn from the prior,
#' or from the posterior of the attached data, and fitted back,
#' each value ranked among the posterior draws of its fit.
#' @examples
#' calibration <- flu_1918_frederick |>
#'   new_housepi_model("final-size") |>
#'   add_parameter(import = "template") |>
#'   check_calibration(nrep = 2, houses = 30)
#' calibration
#' calibration$summarise()
#' @docType class
HousepiCalibrationCheck <- R6::R6Class(
  classname = "HousepiCalibrationCheck",

  ####  Private  ####

  private = list(
    name_ = NA_character_,
    source_ = NA_character_,
    # engine seed of each run the posterior draws come from
    draws_seed_ = NULL,
    # household ids of the attached data each study is fitted alongside
    observed_ids_ = character(),
    rank_levels_ = NA_integer_,
    # estimated parameters, then `log_lik`
    variables_ = NULL,
    # `HousepiSweep` Replicates fitted.
    sweep_ = NULL,

    # Reduce one fit to a row per ranked quantity: rank of `truth`
    # among the posterior draws, diagnostics.
    # * draws too correlated for `size` independent ones -> flagged.
    reduce_ = function(fit, truth) {
      size <- private$ranked_draws_()
      variables <- private$variables_
      params <- setdiff(variables, "log_lik")
      pdraws <- get_parameter_draws(fit)
      log_lik <- private$log_lik_draws_(fit)
      draws <- c(lapply(stats::setNames(nm = params),
                        \(v) as.numeric(pdraws[[v]])),
                 list(log_lik = log_lik$log_lik))
      sm <- diagnose_block_(log_lik, posterior::nchains(log_lik))
      diag <- rbind(sweep_diagnostics_(fit, params),
                    data.frame(sm[c("rhat", "ess_bulk")],
                               flagged = sm$bad_rhat | sm$bad_ess))
      rank <- vapply(variables, function(v) {
        calibration_rank_(draws[[v]], truth[[v]], size)
      }, integer(1L))
      diag[["flagged"]] <- diag[["flagged"]] |
        is.na(diag[["ess_bulk"]]) | diag[["ess_bulk"]] < size
      data.frame(variable = variables, truth = unname(truth[variables]),
                 rank = unname(rank), diag, row.names = NULL)
    },

    # Sum the pointwise log-likelihood draws over the simulated households,
    # chains kept, as variable `log_lik`.
    log_lik_draws_ = function(fit) {
      draws <- posterior::as_draws_matrix(get_log_lik_draws(fit))
      simulated <- setdiff(colnames(draws), private$observed_ids_)
      posterior::draws_df(log_lik = rowSums(draws[, simulated, drop = FALSE]),
                          .nchains = posterior::nchains(draws))
    },

    # Name the check after the source of its parameter values.
    kind_ = function() {
      switch(private$source_,
        prior = "Prior calibration check",
        posterior = "Posterior calibration check"
      )
    },

    # Posterior draws a truth is ranked among, one fewer than the rank levels.
    ranked_draws_ = function() private$rank_levels_ - 1L,

    # Log-threshold of the uniformity test, one per number of replicates.
    # @return `function` Log-threshold of `n` replicates.
    threshold_ = function() {
      n <- unique(self$nrep)
      log_threshold <- vapply(n, calibration_log_threshold_, numeric(1L),
                              size = private$ranked_draws_(),
                              alpha = calibration_alpha_)
      \(m) log_threshold[[match(m, n)]]
    }
  ),

  ####  Active  ####

  active = list(
    #' @field estimates `data.frame`
    #' One row per replicate and ranked quantity,
    #' each estimated parameter and the log-likelihood `log_lik`
    #' \describe{
    #'   \item{truth}{Value the study was simulated at,
    #'     drawn from the prior or a posterior draw.
    #'     For `log_lik`, the log-likelihood of the simulated households
    #'     at that value;}
    #'   \item{rank}{Number of ranked posterior draws below the truth,
    #'     from 0 to `rank_levels - 1`, ties split at random;}
    #'   \item{rhat, ess_bulk}{Convergence diagnostics of the fit;}
    #'   \item{flagged}{`TRUE` when the diagnostics flagged the quantity:
    #'     `r paste(diagnostic_flags_, collapse = "; ")`,
    #'     or when its bulk effective sample size is below `rank_levels - 1`.
    #'     This adds to the diagnostics only when `rank_levels - 1`
    #'     exceeds 100 draws per chain.}
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
    #'     over the ranked quantities of the fit;}
    #'   \item{flagged}{Whether any ranked quantity of the fit was flagged.}
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

    #' @field source `r rd_check_source_(table = FALSE)`
    source = function(value) {
      assert_read_only_(value, "source")
      private$source_
    },

    #' @field rank_levels `integer(1)` Number of values a rank takes.
    rank_levels = function(value) {
      assert_read_only_(value, "rank_levels")
      private$rank_levels_
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
    #' Create the calibration check object.
    #' Use [check_calibration()] instead of calling `$new()`.
    #' @param hh_model `r rd_model_in_`
    #' @param draws `r rd_calib_draws_`
    #' @param nrep `r rd_sweep_nrep_`
    #' @param houses `r rd_calib_houses_`
    #' @param structure_sampler `r rd_sweep_sampler_`
    #' @param rank_levels `r rd_calib_rank_levels_`
    #' @param steps `r rd_smpl_steps_`
    #' @param thin `r rd_smpl_thin_`
    #' @param adapt `r rd_smpl_adapt_`
    #' @param chains `r rd_smpl_chain_`
    #' @param threads `r rd_sweep_threads_`
    #' @param progress `r rd_progress_`
    #' @param ... `r rd_dots_config_`
    initialize = function(hh_model, draws = NULL,
                          nrep = 200L, houses = NA_integer_,
                          structure_sampler = sample_house_structure,
                          rank_levels = 100L,
                          steps = 2000L, thin = 10L, adapt = 1000L,
                          chains = 4L, threads = 1L,
                          progress = TRUE, ...) {
      nrep <- assert_count_(nrep, "Number of replicates")
      rank_levels <- assert_count_(rank_levels, "Number of rank levels")
      if (rank_levels < 2L) stop_h(
        "Number of rank levels should be at least 2.",
        "A rank is taken among {.code rank_levels - 1} posterior draws."
      )
      private$rank_levels_ <- rank_levels
      stored <- resolve_settings_(steps, thin, adapt, chains)$total
      # nolint next: object_usage_linter. read by the cli glue below
      ranked <- private$ranked_draws_()
      if (stored < ranked) stop_h(
        "The fit stores {stored} draw{?s}, too few to rank among {ranked}.",
        "A rank is taken among {.code rank_levels - 1} posterior draws.",
        paste("Raise {.arg steps} or {.arg chains}, or lower {.arg thin} or",
              "{.arg rank_levels}.")
      )
      private$name_ <- hh_model$.__enclos_env__$private$compiled_$name
      private$variables_ <- c(sweep_estimated_(hh_model), "log_lik")

      if (is.null(draws)) {
        private$source_ <- "prior"
        # no column: simulation draws every parameter anew from its prior
        values <- sweep_param_table_(numeric(), nrep)
        conditioned <- NULL
      } else {
        private$source_ <- "posterior"
        # read before the sweep attaches the studies to the model
        posterior <- calibration_posterior_(hh_model, draws, nrep)
        values <- posterior$values
        private$draws_seed_ <- posterior$seed
        # a latent is no data: its like-named column leaves, as in the studies
        conditioned <- get_dataset(hh_model)
        latent <- get_latent(hh_model)[["name"]]
        conditioned <- conditioned[setdiff(names(conditioned), latent)]
        private$observed_ids_ <- sweep_observed_ids_(conditioned)
      }
      private$sweep_ <- HousepiSweep$new(
        hh_model, sampler = structure_sampler, houses = houses, nrep = nrep,
        values = values,
        fit_args = list(steps = steps, thin = thin, adapt = adapt,
                        chains = chains),
        threads = threads, dots = list(...), reducer = private$reduce_,
        progress = progress, label = "Calibration replicates",
        what = private$kind_(), conditioned = conditioned, one = TRUE,
        log_lik = TRUE
      )
      invisible(self)
    },

    #' @description Test the uniformity of the ranks per parameter.
    #' @return `data.frame`
    #' \describe{
    #'   \item{nrep, n_flagged}{Number of replicates ranked,
    #'     and of those flagged;}
    #'   \item{log_gamma}{Logarithm of the uniformity statistic,
    #'     twice the smallest binomial tail probability of the ECDF;}
    #'   \item{log_gamma_threshold}{Value of `log_gamma` below which
    #'     uniformity is rejected, at level 0.05;}
    #'   \item{uniform}{Whether the ranks pass the test.}
    #' }
    summarise = function() {
      calibration_summary_(private$sweep_$estimates, private$ranked_draws_(),
                           private$threshold_())
    },

    #' @description Same as `$summarise()`.
    #' @return See `$summarise()`.
    summarize = function() self$summarise(),

    #' @description Evaluate the empirical CDF of the fractional ranks
    #' `(rank + 1) / rank_levels` per parameter,
    #' with its simultaneous band.
    #' @return `data.frame`
    #' \describe{
    #'   \item{z}{Point the ECDF is evaluated at;}
    #'   \item{ecdf}{Share of the ranks at or below `z`;}
    #'   \item{lower, upper}{Band the ECDF of uniform ranks stays inside,
    #'     at every `z` at once, with probability at least 0.95.}
    #' }
    evaluate_ecdf = function() {
      calibration_ecdf_(private$sweep_$estimates, private$ranked_draws_(),
                        private$threshold_())
    },

    #' @description Print a summary of the sweep.
    print = function() {
      model <- spaced_choice(private$name_)
      kind <- private$kind_()
      conditioned <- if (private$source_ == "posterior") c(
        `Conditioned on` = "{length(private$observed_ids_)} observed households"
      )
      private$sweep_$print_rows(
        "{kind} of model {.val {model}}.",
        c(conditioned,
          `Rank levels` = "{private$rank_levels_}")
      )
      invisible(self)
    }
  )
)

#' @title Check that inference is calibrated against the prior or the data
#' @description
#' Simulation-based calibration: draw parameter values from their prior,
#' simulate outbreaks in the households of the attached data at each draw,
#' fit the resulting household study back with the same model,
#' and rank the drawn value among the posterior draws of its fit.
#' When inference is correct, the ranks are uniform over the replicates.
#' Ranks piling at both ends point to a posterior too narrow,
#' in the middle to one too wide, on one side to a biased one.
#'
#' With posterior `draws`, the check is conditional on the attached data
#' \insertCite{sailynoja_posterior_2026}{housepi}.
#' Each replicate simulates a study at one posterior draw,
#' fits it alongside the attached households,
#' and ranks the draw among the posterior draws of that fit.
#' It checks inference where the posterior of the data puts its mass,
#' rather than over the whole prior.
#' The draws must be those of [run_inference()] on this model and its data.
#'
#' `r rd_sweep_structure_`
#' Every parameter carrying a prior and free to move is drawn anew
#' for each replicate, from its prior or from the posterior draws,
#' and its `value` is ignored.
#' Any other parameter is held at its value.
#' @param hh_model `r rd_model_in_`
#' @param draws `r rd_calib_draws_`
#' @param nrep `r rd_sweep_nrep_`
#' @param houses `r rd_calib_houses_`
#' @param structure_sampler `r rd_sweep_sampler_`
#' @param rank_levels `r rd_calib_rank_levels_`
#' @param steps `r rd_smpl_steps_`
#' @param thin `r rd_smpl_thin_`
#' @param adapt `r rd_smpl_adapt_`
#' @param chains `r rd_smpl_chain_`
#' @param threads `r rd_sweep_threads_`
#' @param progress `r rd_progress_`
#' @param ... `r rd_dots_config_`
#' @return An object of class [HousepiCalibrationCheck].
#' @details
#' Ranks follow \insertCite{modrak_simulation_2025;textual}{housepi}:
#' the drawn value is ranked among `rank_levels - 1` posterior draws
#' thinned evenly from the fit, pooled over the chains,
#' and ties are split at random, so a rank takes `rank_levels` values.
#' Correlated draws make the ranks pile at the ends,
#' even when inference is correct,
#' so a parameter whose bulk effective sample size
#' \insertCite{vehtari_rank_normalization_2021}{housepi} falls below
#' `rank_levels - 1` is flagged.
#' The diagnostics of the fit already flag a bulk effective sample size
#' below 100 per chain,
#' so this flag matters only when `rank_levels - 1` exceeds 100 times
#' the number of chains.
#'
#' Next to the parameters, the log-likelihood of the simulated households
#' is ranked as `log_lik` \insertCite{modrak_simulation_2025}{housepi}.
#' Its truth is evaluated at the values the study was simulated at,
#' and each of its draws sums the pointwise log-likelihood of
#' [get_log_lik_draws()] over the same households.
#' It reads the data and every parameter at once,
#' so it catches a fit that ignores the data,
#' which the ranks of the parameters cannot see.
#' The log-likelihood is evaluated by the model itself:
#' a likelihood that leaves out part of the data leaves it out of the truth
#' as well, and goes unnoticed.
#' On a model with latent variables, `log_lik` adds their log-density,
#' as in [get_log_lik_draws()],
#' and its truth reads the latent values the study was simulated with.
#'
#' Uniformity is tested on the empirical CDF of the ranks, with the
#' simultaneous band and the statistic of
#' \insertCite{sailynoja_graphical_2022;textual}{housepi}.
#' Uniformity is rejected when the statistic falls below its 5th percentile
#' under uniform ranks \insertCite{modrak_simulation_2025}{housepi},
#' so a calibrated sampler is rejected with probability at most 0.05.
#'
#' Posterior draws are thinned to their bulk effective sample size,
#' then `nrep` of them are picked evenly spaced, one per replicate,
#' so the run must hold at least `nrep` independent draws.
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
#'   check_calibration(nrep = 2, houses = 30)
#' calibration$summarise()
#' calibration$overview
#'
#' # conditional on the attached data
#' hh_model <- flu_1918_frederick |>
#'   new_housepi_model("final-size") |>
#'   add_parameter(import = "template")
#' hh_fit <- run_inference(hh_model, chains = 2)
#' calibration <- check_calibration(hh_model, hh_fit, nrep = 2, houses = 30)
#' calibration
#' @export
check_calibration <- function(hh_model,
                              draws = NULL,
                              nrep = 200L,
                              houses = NA_integer_,
                              structure_sampler = sample_house_structure,
                              rank_levels = 100L,
                              steps = 2000L, thin = 10L, adapt = 1000L,
                              chains = 4L, threads = 1L,
                              progress = TRUE,
                              ...) {
  stopifnothousepi(hh_model)
  assert_advanced_settings_("time_step", ...)
  HousepiCalibrationCheck$new(
    hh_model$copy(), draws = draws, nrep = nrep, houses = houses,
    structure_sampler = structure_sampler, rank_levels = rank_levels,
    steps = steps, thin = thin, adapt = adapt,
    chains = chains, threads = threads, progress = progress, ...
  )
}
