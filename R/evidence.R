## housepi - Epidemiological Analysis of Household Outbreaks.
##
## Copyright © 2024-2026 Institut Pasteur (Paris), CNRS and INSERM.
## Distributed under the terms of the GNU General Public License (GPLv3).
## See the LICENSE file for details.
##

## Roxygen text shared by several blocks, read at the use site as `r <name>`.

rd_pwr_frac_ <- "\\eqn{\\frac{r^5}{{n_\\mathrm{rungs}}^5}}{r^5 / n_rungs^5}"
rd_evi_scheme_ <- sprintf("`character(1)` Scheme used for evidence estimation.
One between
\\describe{
  \\item{`\"adaptive\"`}{uses adaptive thermodynamic ladder (default).}
  \\item{`\"power-fractions\"`}{uses inverse temperatures %s.}
}", rd_pwr_frac_)
rd_evi_rungs_ <- "`integer(1)` Number of rungs in thermodynamic ladder."
rd_evi_method_ <- "`character(1)` Estimator computed from the sampled ladder.
One between
\\describe{
  \\item{`\"stepping-stone\"`}{stepping-stone sampling (default).}
  \\item{`\"thermodynamic\"`}{thermodynamic integration.}
}"

############################################################################
############################################################################
############################################################################

#' Warn once per session that log-evidence estimation is experimental.
#' @noRd
warn_logevidence_experimental_ <- function() {
  rlang::warn(
    "Log-evidence estimation is experimental and may change.",
    .frequency = "once",
    .frequency_id = "housepi-logevidence-experimental"
  )
}

#' @title Power-posterior ladder of a housepi model
#' @description
#' Container for the power-posterior ladder sampled by [run_power_posteriors()].
#' Estimate the log-evidence from it with [evaluate_log_evidence()].
#' @references
#' \insertRef{friel_improving_2014}{housepi}
#' @examples
#' ladder <- flu_1918_frederick |>
#'   new_housepi_model("final_size") |>
#'   add_parameter(import = "template") |>
#'   run_power_posteriors(rungs = 10, steps = 200, chains = 2, progress = FALSE)
#' ladder
#' @docType class
HousepiPowerPosteriors <- R6::R6Class(
  classname = "HousepiPowerPosteriors",

  ############################################################################
  ################################## PRIVATE #################################
  ############################################################################

  private = list(
    name_ = NA_character_,
    scheme_ = NA_character_,
    inv_temp_ = NULL,
    mean_log_lik_ = NULL,
    var_log_lik_ = NULL,
    draws_ = NULL,

    # Extract log-likelihood trace at rung k.
    ll_at_ = function(k) {
      posterior::extract_variable(private$draws_, sprintf("log_lik[%d]", k))
    },

    # Stepping-stone estimator over the ladder.
    # Tempered increments max-stabilised via log-mean-exp.
    # @references
    # \insertRef{xie_improving_2011}{housepi}
    stepping_stone_ = function() {
      db <- diff(private$inv_temp_, lag = 1L)
      rss_k <- vapply(seq_along(db), function(k) {
        ll <- private$ll_at_(k)
        ll_peak <- max(ll)
        w <- exp(db[k] * (ll - ll_peak))
        db[k] * ll_peak + log(mean(w))
      }, numeric(1L))
      sum(rss_k)
    },

    # Variance of the stepping-stone estimator.
    # Delta-method approximation with max-stabilisation.
    # @references
    # \insertRef{xie_improving_2011}{housepi}
    stepping_stone_var_ = function() {
      db <- diff(private$inv_temp_, lag = 1L)
      var_rss_k <- vapply(seq_along(db), function(k) {
        ll <- private$ll_at_(k)
        ll_peak <- max(ll)
        w <- exp(db[k] * (ll - ll_peak))
        v <- var(w)
        v / (posterior::ess_mean(w) * mean(w)^2)
      }, numeric(1L))
      sum(var_rss_k)
    },

    # Thermodynamic integration over the ladder.
    # Using trapezoidal rule with variance-based correction.
    # @references
    # \insertRef{friel_improving_2014}{housepi}
    thermodynamic_ = function() {
      if (any(!is.finite(posterior::as_draws_matrix(private$draws_)))) {
        stop_h(
          paste("Thermodynamic integration is not applicable:",
                "ladder holds non-finite log-likelihood draws."),
          "Try to use {.code method = \"stepping-stone\"} instead."
        )
      }
      avg <- private$mean_log_lik_
      var <- private$var_log_lik_
      db <- diff(private$inv_temp_, lag = 1L)
      k <- seq_along(db)
      sum((avg[k + 1L] + avg[k]) * db / 2 +
            (var[k] - var[k + 1L]) * db^2 / 12)
    }
  ),

  ############################################################################
  ################################## ACTIVE ##################################
  ############################################################################

  active = list(
    #' @field scheme `character(1)`
    #' Power-posterior scheme used to sample the ladder.
    scheme = function(value) {
      if (!missing(value)) stop_read_only("scheme")
      private$scheme_
    },

    #' @field draws `posterior::draws_df`
    #'Draws log-likelihood draws under corresponding power-posterior.
    #' `r rd_natural_log_`
    draws = function(value) {
      if (!missing(value)) stop_read_only("draws")
      private$draws_
    },

    #' @field detail `data.frame`
    #' Detail of the thermodynamic ladder
    #' \describe{
    #'   \item{inv_temp}{Inverse temperature}
    #'   \item{mean_log_lik}{Mean of the complete log-likelihood}
    #'   \item{var_log_lik}{Variance of the complete log-likelihood}
    #' }
    #' Mean log-likelihood can legitimately be `-Inf` at inverse temperature 0.
    #' `r rd_natural_log_`
    detail = function(value) {
      if (!missing(value)) stop_read_only("detail")
      data.frame(inv_temp = private$inv_temp_,
                 mean_log_lik = private$mean_log_lik_,
                 var_log_lik = private$var_log_lik_,
                 row.names = NULL)
    }
  ),

  ############################################################################
  ################################## PUBLIC ##################################
  ############################################################################

  public = list(
    #' @description
    #' Internal: [run_power_posteriors()] builds this object, don't call
    #' `$new()` directly.
    #' @param generator Internal model compiler.
    #' @param scheme `r rd_evi_scheme_`
    #' @param rungs `r rd_evi_rungs_`
    initialize = function(generator, scheme = "adaptive", rungs = 30L) {
      ptr <- generator$mcmc_ptr
      output <- .engine$run_power_posteriors(ptr, as.character(scheme),
                                             as.integer(rungs))
      private$name_ <- as.character(generator$raw_name)
      # NOTE: scheme may fall back from "adaptive" to "power-fractions"
      private$scheme_ <- as.character(output$scheme)
      private$inv_temp_ <- output$inv_temp
      private$mean_log_lik_ <- output$mean_log_lik
      private$var_log_lik_ <- output$var_log_lik
      # one tempered log-likelihood column `log_lik[r]` per rung
      private$draws_ <- tag_rcpp_draws_(output$log_lik_draws,
                                        draws_index_(.engine$export_setup(ptr)))
      invisible(self)
    },

    #' @description Estimate the log-evidence from the stored ladder.
    #' @param method `r rd_evi_method_`
    #' @return named `numeric(2)` (`estimate`, `sd`)
    #' `sd` is `NA` if method provides no standard-error estimate.
    #' `r rd_natural_log_`
    #' @references
    #' \insertRef{friel_improving_2014}{housepi}
    #' \insertRef{xie_improving_2011}{housepi}
    evaluate_log_evidence = function(
      method = c("stepping-stone", "thermodynamic")
    ) {
      switch(match.arg(method),
        "stepping-stone" = c(
          estimate = private$stepping_stone_(),
          sd = sqrt(private$stepping_stone_var_())
        ),
        "thermodynamic" = c(
          estimate = private$thermodynamic_(),
          sd = NA_real_
        )
      )
    },

    #' @description Override the default printing behavior.
    print = function() {
      model <- spaced_template_choice(private$name_)
      rungs <- length(private$inv_temp_) - 1L
      print_bullets(c(
        i = "Power posteriors draws of model {.val {model}}.",
        "Scheme: {.val {private$scheme_}}. Rungs: {rungs}."
      ))
      invisible(self)
    }
  )
)

#' @title Sample the power-posterior ladder of a housepi model
#' @description
#' Sample the power-posterior ladder used to estimate the log-evidence.
#' See [evaluate_log_evidence()] to estimate the log-evidence.
#' @param hhmodel `r rd_model_in_`
#' @param scheme `r rd_evi_scheme_`
#' @param rungs `r rd_evi_rungs_`
#' @param steps `r rd_smpl_steps_`
#' @param burn `r rd_smpl_burn_`
#' @param thin `r rd_smpl_thin_`
#' @param adapt `r rd_smpl_adapt_`
#' @param chains `r rd_smpl_chain_`
#' @param threads `r rd_thread_`
#' @param progress `r rd_progress_`
#' @param ... `r rd_dots_config_`
#' @return An object of class `HousepiPowerPosteriors`.
#' @references
#' \insertRef{friel_improving_2014}{housepi}
#' @examples
#' ladder <- flu_1918_frederick |>
#'   new_housepi_model("final_size") |>
#'   add_parameter(import = "template") |>
#'   run_power_posteriors(rungs = 10, steps = 200, chains = 2, progress = FALSE)
#' evaluate_log_evidence(ladder)
#' @export
run_power_posteriors <- function(hhmodel,
                                 scheme = c("adaptive", "power-fractions"),
                                 rungs = 30L,
                                 steps = 1000L, burn = 0L,
                                 thin = 5L, adapt = 500L,
                                 chains = 4L, threads = 1L,
                                 progress = TRUE,
                                 ...) {
  stopifnothousepi(hhmodel)
  warn_logevidence_experimental_()
  assert_advanced_setup_("adapt_window", ...)
  scheme <- match.arg(scheme)
  pvt <- hhmodel$.__enclos_env__$private
  # internal build function
  pvt$build_full_(steps = steps, burn = burn,
                  thin = thin, adapt = adapt,
                  chains = chains, threads = threads,
                  progress = progress, ...)
  # sample the ladder from the internal Rcpp object
  HousepiPowerPosteriors$new(pvt$generator_, scheme = scheme, rungs = rungs)
}

#' @title Estimate the log-evidence
#' @description
#' Estimate the logarithm of the evidence (aka marginal likelihood) from a
#' power-posteriors ladder sampled with [run_power_posteriors()].
#' @param hhdraws A `HousepiPowerPosteriors` object.
#' @param method `r rd_evi_method_`
#' @return named `numeric(2)` (`estimate`, `sd`)
#' `sd` is `NA` if method provides no standard-error estimate.
#' `r rd_natural_log_`
#' @references
#' \insertRef{friel_improving_2014}{housepi}
#' \insertRef{xie_improving_2011}{housepi}
#' @examples
#' ladder <- flu_1918_frederick |>
#'   new_housepi_model("final_size") |>
#'   add_parameter(import = "template") |>
#'   run_power_posteriors(rungs = 10, steps = 200, chains = 2, progress = FALSE)
#' evaluate_log_evidence(ladder)
#' evaluate_log_evidence(ladder, method = "thermodynamic")
#' @export
evaluate_log_evidence <- function(
  hhdraws, method = c("stepping-stone", "thermodynamic")
) {
  if (inherits(hhdraws, "HousepiModel")) {
    stop_h(
      "Impossible to evaluate log-evidence on {.cls HousepiModel} directly.",
      paste("Sample a ladder first with {.fn run_power_posteriors},",
            "then pass it to {.fn evaluate_log_evidence}.")
    )
  }
  if (!inherits(hhdraws, "HousepiPowerPosteriors")) stop_h(
    "{.arg hhdraws} must be of class {.cls HousepiPowerPosteriors}.",
    "Create one with {.fn run_power_posteriors}."
  )
  hhdraws$evaluate_log_evidence(method = method)
}
