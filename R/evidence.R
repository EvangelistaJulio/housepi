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

#' @title Log-evidence estimate of a housepi model
#' @description
#' Container for log-evidence estimate returned by [evaluate_log_evidence()].
#' @docType class
HousepiEvidence <- R6::R6Class(
  classname = "HousepiEvidence",

  ############################################################################
  ################################## PRIVATE #################################
  ############################################################################

  private = list(
    name_ = NA_character_,
    estimate_ = NA_real_,
    scheme_ = NA_character_,
    detail_ = NULL
  ),

  ############################################################################
  ################################## ACTIVE ##################################
  ############################################################################

  active = list(
    #' @field estimate `numeric(1)`
    #' Log-evidence (marginal log-likelihood) estimate.
    estimate = function(value) {
      if (!missing(value)) stop_read_only("estimate")
      private$estimate_
    },

    #' @field scheme `character(1)`
    #' Power-posterior scheme used for the estimate.
    scheme = function(value) {
      if (!missing(value)) stop_read_only("scheme")
      private$scheme_
    },

    #' @field detail `data.frame`
    #' Detail of the thermodynamic ladder
    #' \describe{
    #'   \item{inv_temper}{Inverse temperature}
    #'   \item{mean_log_lik}{Mean of the complete log-likelihood}
    #'   \item{var_log_lik}{Variance of the complete log-likelihood}
    #' }
    detail = function(value) {
      if (!missing(value)) stop_read_only("detail")
      private$detail_
    }
  ),

  ############################################################################
  ################################## PUBLIC ##################################
  ############################################################################

  public = list(
    #' @description Estimate the log-evidence and store the result.
    #' @param generator `HousepiGenerator` Helper class for model generation.
    #' @param scheme `r rd_evi_scheme_`
    #' @param rungs `r rd_evi_rungs_`
    initialize = function(generator, scheme = "power-fractions", rungs = 40) {
      output <- .engine$evaluate_log_evidence(generator$mcmc_ptr,
                                              as.character(scheme),
                                              as.integer(rungs))
      private$name_ <- as.character(generator$raw_name)
      private$estimate_ <- output$estimate
      private$scheme_ <- as.character(output$scheme)
      mask <- c("inv_temper", "mean_log_lik", "var_log_lik")
      private$detail_ <- as.data.frame(output[mask])
      invisible(self)
    },

    #' @description Override the default printing behavior.
    print = function() {
      fmt_cli_info(c(
        sprintf("Log-evidence of model '%s': %.4g.",
                spaced_template_choice(private$name_), private$estimate_),
        sprintf("Power posteriors, scheme '%s': %d rungs.",
                private$scheme_, nrow(private$detail_) - 1L),
        "Thermodynamic ladder in `$detail`."
      ))
      invisible(self)
    }
  )
)
