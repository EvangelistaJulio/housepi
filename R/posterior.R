## housepi - Epidemiological Analysis of Household Outbreaks.
##
## Copyright © 2026 Institut Pasteur (Paris), CNRS and INSERM.
## Distributed under the terms of the GNU General Public License (GPLv3).
## See the LICENSE file for details.
##

#' @importFrom posterior subset_draws
NULL

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

#' @title Set the global credible interval level
#' @description
#' Set credible interval level used by housepi functions.
#' It also drives
#' the quantile bounds returned by [summarise_draws()] when no summary
#' function is supplied.
#' @param value `numeric(1)`. Must lie in `(0.5, 1)`. Default is `0.95`.
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

#' @title Summarise posterior draws using the housepi credible interval
#' @description
#' Drop-in replacement for [posterior::summarise_draws()] that, when called
#' without summary functions, returns `median`, the lower/upper credible
#' interval bounds derived from the global credible interval level
#' (see [set_global_cri_level()]), `ess_bulk` and `rhat`.
#' When summary functions are supplied through `...`, the call is forwarded
#' to [posterior::summarise_draws()] unchanged.
#' @param .x A `draws` object (or anything coercible to one).
#' @param ... Summary functions, as in [posterior::summarise_draws()].
#' If empty, housepi defaults are used.
#' @param .args,.num_args,.cores Forwarded to [posterior::summarise_draws()].
#' @return A `draws_summary` `data.frame`.
#' @examples
#' if (requireNamespace("posterior", quietly = TRUE)) {
#'   d <- posterior::example_draws()
#'   summarise_draws(d)
#' }
#' @export
summarise_draws <- function(.x, ..., .args = list(),
                            .num_args = getOption("posterior.num_args", list()),
                            .cores = 1) {
  funs <- list(...)
  if (length(funs) == 0) {
    probs <- unname(cri_level_probs(get_global_cri_level()))
    return(posterior::summarise_draws(
      .x,
      mean,
      median,
      sd,
      mad,
      ~ posterior::quantile2(.x, probs = probs),
      rhat = posterior::rhat,
      ess_bulk = posterior::ess_bulk,
      ess_tail = posterior::ess_tail,
      .args = .args,
      .num_args = .num_args,
      .cores = .cores
    ))
  }
  posterior::summarise_draws(.x, ..., .args = .args,
                             .num_args = .num_args, .cores = .cores)
}

#' @rdname summarise_draws
#' @export
summarize_draws <- summarise_draws

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
      n_steps <- setup[["steps"]]
      n_burn <- setup[["burn"]]
      n_thin <- setup[["thin"]]
      n_chains <- setup[["chains"]]
      c_sample_size <- (n_steps - n_burn) %/% n_thin

      private$chain_tags_ <- rep(seq_len(n_chains), each = c_sample_size)
      private$iter_names_ <- n_burn + seq_len(c_sample_size) * n_thin
      invisible(self)
    },

    # @description Convert draws from housepi Rcpp to `posterior.`
    # @param x `matrix` or `data.frame`. Draws exported from Rcpp core.
    # @return `posterior::draws_df`
    to_posterior_draws_ = function(x) {
      if (length(x) == 0) return(NULL)

      # split into a list of matrix
      list_trs <- lapply(split(seq_len(nrow(x)), private$chain_tags_),
                         function(idx) x[idx, , drop = FALSE])

      # create empty array
      n_chains <- length(list_trs)
      par_names <- colnames(x)
      n_iter <- length(private$iter_names_)
      array_trs <- array(
        NA_real_, dim = c(n_iter, n_chains, length(par_names)),
        dimnames = list(private$iter_names_, seq_len(n_chains), par_names)
      )

      # associate draws from Cpp
      for (i in seq_along(list_trs)) array_trs[, i, ] <- list_trs[[i]]
      output <- posterior::as_draws_df(array_trs)
      output$.iteration <- rep(private$iter_names_, times = n_chains)
      output
    },

    # @description Return empty draws.
    # @return `posterior::draws_df`
    empty_draws_ = function() {
      posterior::as_draws_df(data.frame())
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
      posterior::bind_draws(private$param_draws_, aux, self$prior_draws)
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
