## housepi - Epidemiological Analysis of Household Outbreaks.
##
## Copyright © 2024-2026 Institut Pasteur (Paris), CNRS and INSERM.
## Distributed under the terms of the GNU General Public License (GPLv3).
## See the LICENSE file for details.
##

#' @importFrom posterior summarise_draws
#' @export
posterior::summarise_draws

#' @importFrom posterior summarize_draws
#' @export
posterior::summarize_draws

# Thresholds for the post-inference diagnostic warnings.
diagnostic_ <- list(
  ss_min = 400,    # min sample size # Vehtari 2021
  rhat_max = 1.05,
  neff_min = 0.1   # ESS ratio (see bayesplot::mcmc_neff)
)

# Validate a credible interval level.
validate_cri_level <- function(value) {
  if (length(value) != 1 || !is.numeric(value) || is.na(value)) stop_h(
    "Credible interval level must be a single non-NA numeric value."
  )
  if ((value <= 0.5) || (value >= 1)) stop_h(
    "Credible interval level must be a value between 0.5 and 1."
  )
  value
}

#' @title Set global credible interval level
#' @description
#' Set credible interval level used by housepi objects.
#' @param value `numeric(1)` Must lie in `(0.5, 1)`. Default is `0.95`.
#' @return Invisibly returns the previous value.
#' @examples
#' old <- set_global_cri_level(0.90)
#' set_global_cri_level(old)
#' @export
set_global_cri_level <- function(value) {
  value <- validate_cri_level(value)
  old <- getOption("housepi.cri_level", 0.95)
  options(housepi.cri_level = value)
  invisible(old)
}

# Return the global credible interval level (option "housepi.cri_level").
# Not exported: housepi internals must always pull through this helper.
get_global_cri_level <- function() {
  validate_cri_level(getOption("housepi.cri_level", 0.95))
}

# Lower/upper quantile probabilities corresponding to a credible interval level.
cri_level_probs <- function(cri_level) {
  half_tail <- 0.5 * (1 - cri_level)
  c(lower = half_tail, upper = 1. - half_tail)
}

# Prepend "housepi_draws" class tag for summarise_draws S3 method dispatch.
tag_as_housepi_draws <- function(x) {
  if (is.null(x)) return(x)
  class(x) <- unique(c("housepi_draws", class(x)))
  x
}

#' @title Summarise housepi posterior draws using the global credible interval
#' @description
#' S3 method on [posterior::summarise_draws()] for objects tagged
#' `"housepi_draws"` (returned by housepi getters and active fields on
#' [`HousepiPosteriorDraws`]). When no summary functions are passed via
#' `...`, the housepi defaults are used and the quantile pair tracks
#' [set_global_cri_level()]. Otherwise dispatch falls through to standard
#' [posterior::summarise_draws()] behaviour.
#' @param .x A housepi-tagged draws object, or a [`HousepiPosteriorDraws`].
#' @param ... Summary functions, as in [posterior::summarise_draws()].
#' @param .args,.num_args,.cores Forwarded to [posterior::summarise_draws()].
#' @return A `draws_summary` `data.frame`.
#' @name summarise_draws.housepi_draws
#' @exportS3Method posterior::summarise_draws
summarise_draws.housepi_draws <- function(
    .x, ..., .args = list(),
    .num_args = getOption("posterior.num_args", list()),
    .cores = 1) {
  # strip the tag so the forwarded call dispatches to the next class
  class(.x) <- setdiff(class(.x), "housepi_draws")
  if (length(list(...)) > 0) return(
    posterior::summarise_draws(.x, ..., .args = .args,
                               .num_args = .num_args, .cores = .cores)
  )
  probs <- unname(cri_level_probs(get_global_cri_level()))
  quantile_fn <- function(x) posterior::quantile2(x, probs = probs)
  posterior::summarise_draws(
    .x,
    "mean", "median", "sd", "mad",
    quantile_fn,
    rhat = posterior::rhat,
    ess_bulk = posterior::ess_bulk,
    ess_tail = posterior::ess_tail,
    .args = .args, .num_args = .num_args, .cores = .cores
  )
}

#' @rdname summarise_draws.housepi_draws
#' @exportS3Method posterior::summarise_draws
summarise_draws.HousepiPosteriorDraws <- function(
    .x, ..., .args = list(),
    .num_args = getOption("posterior.num_args", list()),
    .cores = 1) {
  params <- .x$parameter_draws
  augm <- .x$augmented_draws
  draws <- if (posterior::nvariables(augm) == 0) params
  else if (posterior::nvariables(params) == 0) augm
  else tag_as_housepi_draws(posterior::bind_draws(params, augm))
  posterior::summarise_draws(draws, ..., .args = .args,
                             .num_args = .num_args, .cores = .cores)
}

#' @rdname summarise_draws.housepi_draws
#' @exportS3Method posterior::summarize_draws
summarize_draws.housepi_draws <- summarise_draws.housepi_draws

#' @rdname summarise_draws.housepi_draws
#' @exportS3Method posterior::summarize_draws
summarize_draws.HousepiPosteriorDraws <- summarise_draws.HousepiPosteriorDraws

############################################################################
############################################################################
#########################  Housepi Posterior Draws  ########################
############################################################################
############################################################################

#' @title Posterior draws from housepi model inference
#' @description
#' Posterior draws as a result of housepi inference on a household model.
#' This class stores posterior draws of model parameters,
#' priors, total and household point-wise log-likelihood.
#' If requested so during inference, it also stores augmented variables.
#' @docType class
HousepiPosteriorDraws <- R6::R6Class(
  classname = "HousepiPosteriorDraws",

  ############################################################################
  ################################## PRIVATE #################################
  ############################################################################

  private = list(
    # @field
    name_ = NA_character_,
    param_draws_ = NULL,
    prior_draws_ = NULL,
    augm_draws_ = NULL,
    total_log_lik_draws_ = NULL,
    point_log_lik_draws_ = NULL,
    chain_tags_ = NULL,
    iter_names_ = NULL,

    ############################################################################
    ############################################################################
    ############################################################################

    # Store auxiliary tags from inference setup.
    pre_process_setup_ = function(setup) {
      steps <- setup[["steps"]]
      burn <- setup[["burn"]]
      thin <- setup[["thin"]]
      chains <- setup[["chains"]]
      chain_size <- setup_sample_size_(steps = steps, burn = burn,
                                       thin = thin, chains = chains)$chain
      private$chain_tags_ <- rep(seq_len(chains), each = chain_size)
      private$iter_names_ <- burn + seq_len(chain_size) * thin
      invisible(self)
    },

    # @description Convert draws from housepi Rcpp to `posterior.`
    # @param df `data.frame`. Draws exported from Rcpp core.
    # @return `posterior::draws_df`
    to_posterior_draws_ = function(df) {
      if (length(df) == 0) return(NULL)
      n_chains <- length(unique(private$chain_tags_))
      df$.chain <- private$chain_tags_
      df$.iteration <- rep(private$iter_names_, times = n_chains)
      df$.draw <- seq_len(nrow(df))
      tag_as_housepi_draws(posterior::as_draws_df(df))
    },

    # @description Return empty draws.
    # @return `posterior::draws_df`
    empty_draws_ = function() {
      tag_as_housepi_draws(posterior::as_draws_df(data.frame()))
    },

    # Stan-style post-inference diagnostics.
    # See https://mc-stan.org/learn-stan/diagnostics-warnings.html
    check_diagnostics_ = function() {
      niter <- posterior::niterations(private$total_log_lik_draws_)
      nchains <- posterior::nchains(private$total_log_lik_draws_)
      ss <- niter * nchains
      if (ss < diagnostic_[["ss_min"]]) {
        warn_h(
          sprintf("Total sample size was below %d.",
                  diagnostic_[["ss_min"]]),
          "Convergence diagnostics may be unreliable and were skipped."
        )
        return(invisible())
      }

      private$check_diagnostics_one_(private$param_draws_, "parameter")
      private$check_diagnostics_one_(private$augm_draws_, "augmented variable")
      invisible()
    },

    # Emits at most one warning per family combination.
    check_diagnostics_one_ = function(d, family) {
      if (is.null(d) || posterior::nvariables(d) == 0L) return(invisible())

      rhat_max <- diagnostic_[["rhat_max"]]
      neff_min <- diagnostic_[["neff_min"]]
      ss <- posterior::niterations(d) * posterior::nchains(d)
      sm <- posterior::summarise_draws(d)

      bad_rhat <- which(!is.na(sm$rhat) & sm$rhat > rhat_max)
      bad_bulk <- which(!is.na(sm$ess_bulk) & sm$ess_bulk / ss < neff_min)
      both <- intersect(bad_rhat, bad_bulk)
      only_rhat <- setdiff(bad_rhat, both)
      only_bulk <- setdiff(bad_bulk, both)

      sms <- c()
      if (length(both) > 0L) sms <- c(sms,
        sprintf("%s %s with R-hat > %.2f and Bulk ESS < %.0f%%:%s",
                format_there_were_was(length(both)),
                as_plural(family, length(both)),
                rhat_max, 100. * neff_min,
                character_to_string_list(sm$variable[both])),
        "Chains may have not mixed; mean and median may be unreliable."
      )
      if (length(only_rhat) > 0L) sms <- c(sms,
        sprintf("%s %s with R-hat > %.2f:%s",
                format_there_were_was(length(only_rhat)),
                as_plural(family, length(only_rhat)),
                rhat_max,
                character_to_string_list(sm$variable[only_rhat])),
        "Chains may have not mixed."
      )
      if (length(only_bulk) > 0L) sms <- c(sms,
        sprintf("%s %s with Bulk Effective Sample Size (ESS) < %.0f%%:%s",
                format_there_were_was(length(only_bulk)),
                as_plural(family, length(only_bulk)),
                100. * neff_min,
                character_to_string_list(sm$variable[only_bulk])),
        "Mean and median may be unreliable."
      )
      if (length(sms) > 0L) warn_h(sms)
      invisible()
    }
  ),

  ############################################################################
  ################################## ACTIVE ##################################
  ############################################################################

  active = list(
    #' @field parameter_draws `posterior::draws_df`.
    #' Draws of parameters.
    parameter_draws = function(value) {
      if (!missing(value)) stop_read_only("param_draws")
      if (is.null(private$param_draws_)) private$empty_draws_()
      else private$param_draws_
    },

    #' @field prior_draws `posterior::draws_df`.
    #' Draws of parameter priors.
    prior_draws = function(value) {
      if (!missing(value)) stop_read_only("prior_draws")
      if (is.null(private$prior_draws_)) private$empty_draws_()
      else posterior::mutate_variables(
        private$prior_draws_,
        lprior = rowSums(posterior::as_draws_matrix(private$prior_draws_))
      )
    },

    #' @field augmented_draws `posterior::draws_df`.
    #' Draws selected individuals augmented variables.
    augmented_draws = function(value) {
      if (!missing(value)) stop_read_only("augm_draws")
      if (is.null(private$augm_draws_)) private$empty_draws_()
      else private$augm_draws_
    },

    #' @field log_lik_draws `posterior::draws_df`.
    #' Draws of households pointwise log-likelihood.
    log_lik_draws = function(value) {
      if (!missing(value)) stop_read_only("log_lik_draws")
      private$point_log_lik_draws_
    },

    #' @field total_log_lik_draws `posterior::draws_df`.
    #' Draws of the total log-likelihood.
    total_log_lik_draws = function(value) {
      if (!missing(value)) stop_read_only("total_log_lik_draws")
      private$total_log_lik_draws_
    },

    #' @field complete_draws `posterior::draws_df`.
    #' Draws of parameters, priors and pointwise log-likelihood values.
    complete_draws = function(value) {
      if (!missing(value)) stop_read_only("complete_draws")
      # rename pointwise log-likelihood
      aux <- private$point_log_lik_draws_
      n_vars <- posterior::nvariables(aux)
      names(aux)[1:n_vars] <- paste0("log_lik[", 1:n_vars, "]")
      # bind all together
      # NOTE: bind_draws strips the tag, so re-stamp
      tag_as_housepi_draws(
        posterior::bind_draws(private$param_draws_, aux, self$prior_draws)
      )
    }
  ),

  ############################################################################
  ################################## PUBLIC ##################################
  ############################################################################

  public = list(
    #' @description
    #' Store posterior::draws_df internally.
    #' @param name `character(1)`. Model name.
    #' @param ptr `Rcpp::XPtr<MCMC>`
    initialize = function(name, ptr) {
      # inherit model name
      private$name_ <- as.character(name)
      # process inference setup
      private$pre_process_setup_(.engine$export_setup(ptr))

      # store cpp draws converted to posterior
      private$total_log_lik_draws_ <- private$to_posterior_draws_(
        .engine$flush_total_log_lik_draws(ptr)
      )
      private$param_draws_ <- private$to_posterior_draws_(
        .engine$flush_param_draws(ptr)
      )
      private$prior_draws_ <- private$to_posterior_draws_(
        .engine$flush_param_log_prior_draws(ptr)
      )
      private$augm_draws_ <- private$to_posterior_draws_(
        .engine$flush_augm_draws(ptr)
      )
      private$point_log_lik_draws_ <- private$to_posterior_draws_(
        .engine$flush_house_log_lik_draws(ptr)
      )

      # free auxiliary memory
      private$chain_tags_ <- NULL
      private$iter_names_ <- NULL

      # Stan-style automatic diagnostics
      private$check_diagnostics_()

      invisible()
    },

    #' @description Override the default printing behavior.
    print = function() {
      niter <- posterior::niterations(private$total_log_lik_draws_)
      nchains <- posterior::nchains(private$total_log_lik_draws_)

      if (is.null(private$param_draws_)) nparams <- 0
      else nparams <- posterior::nvariables(private$param_draws_)

      if (is.null(private$augm_draws_)) naugm <- 0
      else naugm <- posterior::nvariables(private$augm_draws_)

      fmt_cli_info(c(
        sprintf("Posterior draws of model '%s'.",
                spaced_template_choice(private$name_)),
        sprintf(
          "%s iteration%s, %s chain%s%s%s.",
          niter, ifelse(niter > 1, "s", ""),
          nchains, ifelse(nchains > 1, "s", ""),
          ifelse(nparams == 0, "",
                 sprintf(", %s parameter%s",
                         nparams, ifelse(nparams > 1, "s", ""))),
          ifelse(naugm == 0, "", sprintf(", %s augmented", naugm))
        )
      ))
      invisible()
    }
  )
)

############################################################################
############################################################################
############################################################################

#' @title Get posterior draws of respective variables
#' @description Return a posterior draws dataframe with the MCMC samples.
#' @param hhdraws `HousepiPosteriorDraws` object.
#' Generated by `HousepiModel` inference.
#' @return `posterior::draws_df`
#' @rdname get_draws
#' @export
get_parameter_draws <- function(hhdraws) {
  stopifnot(inherits(hhdraws, "HousepiPosteriorDraws"))
  hhdraws$parameter_draws
}

#' @rdname get_draws
#' @export
get_augmented_draws <- function(hhdraws) {
  stopifnot(inherits(hhdraws, "HousepiPosteriorDraws"))
  hhdraws$augmented_draws
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
