## housepi - Epidemiological Analysis of Household Outbreaks.
##
## Copyright © 2024-2026 Institut Pasteur (Paris), CNRS and INSERM.
## Distributed under the terms of the GNU General Public License (GPLv3).
## See the LICENSE file for details.
##

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

    # Thermodynamic integration over the ladder.
    # Using trapezoidal rule with variance-based correction.
    # @references
    # \insertRef{friel_improving_2014}{housepi}
    thermodynamic_ = function() {
      avg <- private$mean_log_lik_
      var <- private$var_log_lik_
      dt <- diff(private$inv_temp_, lag = 1)
      k <- seq_along(dt)
      sum((avg[k + 1L] + avg[k]) * dt / 2 +
            (var[k] - var[k + 1L]) * dt^2 / 12)
    },

    # Stepping-stone estimator over the ladder.
    # Tempered increments max-stabilised via log-mean-exp.
    # @references
    # \insertRef{xie_improving_2011}{housepi}
    stepping_stone_ = function() {
      dt <- diff(private$inv_temp_, lag = 1)
      contrib <- vapply(seq_along(dt), function(k) {
        log_lik <- posterior::extract_variable(private$draws_,
                                               sprintf("log_lik[%d]", k))
        peak <- max(log_lik)
        dt[k] * peak + log(mean(exp(dt[k] * (log_lik - peak))))
      }, numeric(1))
      sum(contrib)
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
    #' @description Sample the power-posterior ladder and store the draws.
    #' @param generator `HousepiGenerator` Helper class for model generation.
    #' @param scheme `r rd_evi_scheme_`
    #' @param rungs `r rd_evi_rungs_`
    #' @references
    #' \insertRef{friel_improving_2014}{housepi}
    initialize = function(generator, scheme = "power-fractions", rungs = 40L) {
      ptr <- generator$mcmc_ptr
      output <- .engine$run_power_posteriors(ptr, as.character(scheme),
                                             as.integer(rungs))
      private$name_ <- as.character(generator$raw_name)
      private$scheme_ <- as.character(scheme)
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
    #' @return `numeric(1)` Log-evidence estimate.
    #' @references
    #' \insertRef{friel_improving_2014}{housepi}
    #' \insertRef{xie_improving_2011}{housepi}
    evaluate_log_evidence = function(
      method = c("thermodynamic", "stepping-stone")
    ) {
      method <- match.arg(method)
      switch(method,
        "thermodynamic" = private$thermodynamic_(),
        "stepping-stone" = private$stepping_stone_()
      )
    },

    #' @description Override the default printing behavior.
    print = function() {
      fmt_cli_info(c(
        sprintf("Power posteriors draws of model '%s'",
                spaced_template_choice(private$name_)),
        sprintf("Scheme: '%s'. Rungs: %d.",
                private$scheme_, length(private$inv_temp_) - 1L)
      ))
      invisible(self)
    }
  )
)
