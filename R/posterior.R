#
#        __
#       / /_  ____  __  __ ____ __  ____  __
#      / __ \/ __ \/ / / / ___/ _ \/ __ \/ /
#     / / / / /_/ / /_/ (__  )  __/ /_/ / /
#    /_/ /_/\____/\____/____/\___/ ____/_/
#                               /_/
#
# housepi - Epidemiological Analysis of Household Outbreaks.
#
# Author: Francesco Camaglia
# Copyright © 2025 Institut Pasteur (Paris), CNRS and INSERM.
#
# Distributed under the terms of the GNU General Public License (GPLv3).
# See the LICENSE file for details.
#
# posterior.R
#

#' @importFrom posterior summarize_draws
#' @importFrom posterior subset_draws
NULL

############################################################################
############################################################################
#########################  Housepi Posterior Draws  ########################
############################################################################
############################################################################

#' @title Housepi Posterior Draws
#' @description
#' Posterior draws as a results of housepi model inference.
#'
#' @docType class
HousepiPosteriorDraws <- R6Class(
  classname = "HousepiPosteriorDraws",

  ############################################################################
  ################################## PRIVATE #################################
  ############################################################################

  private = list(
    # @field
    .name = NA_character_,
    .param_draws = NULL,
    .prior_draws = NULL,
    .augm_draws = NULL,
    .total_log_lik_draws = NULL,
    .point_log_lik_draws = NULL,
    .chain_tags = NULL,
    .iter_names = NULL,
    .all_iter_names = NULL,

    # Store auxiliary tags from inference setup.
    .process_setup = function(setup) {
      n_steps <- setup[["steps"]]
      n_burnin <- setup[["burnin"]]
      n_thin <- setup[["thin"]]
      n_chains <- setup[["chains"]]
      c_sample_size <- (n_steps - n_burnin) %/% n_thin

      private$.chain_tags <- rep(seq_len(n_chains), each = c_sample_size)
      private$.iter_names <- n_burnin + seq_len(c_sample_size) * n_thin
      private$.all_iter_names <- rep(private$.iter_names, times = n_chains)
      invisible(self)
    },

    # @description
    # Convert draws export from `Rcpp_ModelHousepi` to `posterior.`
    #
    # @param x
    # `matrix` or `data.frame` with draws as exported from Rcpp core.
    #
    # @return
    # `posterior::draws_df`
    .cpp_draws_to_posterior = function(x) {
      if (!is.matrix(x)) stop_h(
        "Something went wrong while exporting posterior draws."
      )
      if (length(x) == 0) return(NULL)

      # split into a list of matrix
      list_trs <- lapply(split(seq_len(nrow(x)), private$.chain_tags),
                         function(idx) x[idx, drop = FALSE])

      # create empty array
      n_chain <- length(list_trs)
      par_names <- colnames(x)
      n_iter <- length(private$.all_iter_names)
      array_trs <- array(
        NA_real_, dim = c(n_iter, n_chain, length(par_names)),
        dimnames = list(private$.all_iter_names, seq_len(n_chain), par_names)
      )

      # associate draws from Cpp
      for (i in seq_along(list_trs)) array_trs[, i, ] <- list_trs[[i]]
      output <- posterior::as_draws_df(array_trs)
      output$.iteration <- private$.iter_names
      output
    }
  ),

  ############################################################################
  ################################## ACTIVE ##################################
  ############################################################################

  active = list(
    #' @field param_draws `posterior::draws_df`
    #' Posterior draws data.frame of parameters.
    param_draws = function(value) {
      if (!missing(value)) stop_h("'param_draws' is read only.")
      private$.param_draws
    },

    #' @field prior_draws `posterior::draws_df`
    #' Posterior draws data.frame of values of parameter priors.
    prior_draws = function(value) {
      if (!missing(value)) stop_h("'prior_draws' is read only.")
      posterior::mutate_variables(
        private$.prior_draws,
        lprior = rowSums(posterior::as_draws_matrix(private$.prior_draws))
      )
    },

    #' @field augm_draws `posterior::draws_df`
    #' Posterior draws data.frame of selected individuals augmented variables.
    augm_draws = function(value) {
      if (!missing(value)) stop_h("'augm_draws' is read only.")
      private$.augm_draws
    },

    #' @field log_lik_draws `posterior::draws_df`
    #' Posterior draws data.frame of households pointwise log-likelihood.
    log_lik_draws = function(value) {
      if (!missing(value)) stop_h("'log_lik_draws' is read only.")
      private$.point_log_lik_draws
    },

    #' @field total_log_lik_draws `posterior::draws_df`
    #' Posterior draws data.frame of the total log-likelihood.
    total_log_lik_draws = function(value) {
      if (!missing(value)) stop_h("'total_log_lik_draws' is read only.")
      private$.total_log_lik_draws
    },

    #' @field complete_draws `posterior::draws_df`
    #' Parameters draws complete of pointwise log-likelihood and priors.
    complete_draws = function(value) {
      if (!missing(value)) stop_h("'complete_draws' is read only.")
      # rename pointwise log-likelihood
      aux <- private$.point_log_lik_draws
      n_vars <- posterior::nvariables(aux)
      names(aux)[1:n_vars] <- paste0("log_lik[", 1:n_vars, "]")
      # bind all together
      posterior::bind_draws(private$.param_draws, aux, self$prior_draws)
    }
  ),

  ############################################################################
  ################################## PUBLIC ##################################
  ############################################################################

  public = list(
    #' store posterior::draws_df internally
    #' @param name `character(1)` Model name.
    #' @param rcpp_object `Rcpp_ModelHousepi`
    initialize = function(name, rcpp_object) {
      # inherit model name
      private$.name <- as.character(name)
      # process inference setup
      private$.process_setup(rcpp_object$export_setup())

      # store cpp draws concverted to posterior
      private$.total_log_lik_draws <- private$.cpp_draws_to_posterior(
        rcpp_object$flush_total_log_lik_draws()
      )
      private$.param_draws <- private$.cpp_draws_to_posterior(
        rcpp_object$flush_param_draws()
      )
      private$.prior_draws <- private$.cpp_draws_to_posterior(
        rcpp_object$flush_param_lprior_draws()
      )
      private$.augm_draws <- private$.cpp_draws_to_posterior(
        rcpp_object$flush_augm_draws()
      )
      private$.point_log_lik_draws <- private$.cpp_draws_to_posterior(
        rcpp_object$flush_house_log_lik_draws()
      )

      # free auxiliary memory
      private$.chain_tags <- NULL
      private$.iter_names <- NULL

      invisible()
    },

    #' @description
    #' Override the default printing behavior.
    print = function() {
      cli::cli_alert_info(sprintf(
        "Posterior draws of model {.pkg %s}.\n%s",
        spaced_template_choice(private$.name),
        sprintf(
          "%s iterations, %s chains%s%s.",
          posterior::niterations(private$.total_log_lik_draws),
          posterior::nchains(private$.total_log_lik_draws),
          ifelse(is.null(private$.param_draws), "",
                 sprintf(", %s parameters",
                         posterior::nvariables(private$.param_draws))),
          ifelse(is.null(private$.augm_draws), "",
                 sprintf(", %s augmented",
                         posterior::nvariables(private$.augm_draws)))
        )
      ))
      invisible()
    },

    #' @description
    #' Summarize all posterior draws.
    #' @return `posterior::draws_df`
    summarize_all_draws = function() {
      x1 <- x2 <- x3 <- NULL
      if (posterior::is_draws(private$.param_draws)) {
        x1 <- summarize_draws(private$.param_draws)
      }
      if (posterior::is_draws(private$.augm_draws)) {
        x2 <- summarize_draws(private$.augm_draws)
      }
      if (posterior::is_draws(private$.total_log_lik_draws)) {
        x3 <- summarize_draws(private$.total_log_lik_draws)
      }
      rbind(x1, x2, x3)
    }
  )
)

############################################################################
############################################################################
############################################################################

#' @title Get posterior draws
#'
#' @description
#' Posterior draws dataframe.
#'
#' @param hhdraws `HousepiPosteriorDraws` object
#' as generated by `HousepiModel` inference.
#'
#' @return `posterior::draws_df`
#'
#' @rdname get_draws
#' @export
get_param_draws <- function(hhdraws) {
  stopifnot(inherits(hhdraws, "HousepiPosteriorDraws"))
  hhdraws$param_draws
}

#' @rdname get_draws
#' @export
get_augm_draws <- function(hhdraws) {
  stopifnot(inherits(hhdraws, "HousepiPosteriorDraws"))
  hhdraws$augm_draws
}

#' @rdname get_draws
#' @export
get_log_lik_draws <- function(hhdraws) {
  stopifnot(inherits(hhdraws, "HousepiPosteriorDraws"))
  hhdraws$log_lik_draws
}

#' @rdname get_draws
#' @export
get_total_log_lik_draws <- function(hhdraws) {
  stopifnot(inherits(hhdraws, "HousepiPosteriorDraws"))
  hhdraws$total_log_lik_draws
}

#' @rdname get_draws
#' @export
get_prior_draws <- function(hhdraws) {
  stopifnot(inherits(hhdraws, "HousepiPosteriorDraws"))
  hhdraws$prior_draws
}

#' @rdname get_draws
#' @export
get_complete_draws <- function(hhdraws) {
  stopifnot(inherits(hhdraws, "HousepiPosteriorDraws"))
  hhdraws$complete_draws
}
