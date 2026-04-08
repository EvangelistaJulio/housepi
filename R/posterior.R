## housepi - Epidemiological Analysis of Household Outbreaks.
##
## Copyright © 2026 Institut Pasteur (Paris), CNRS and INSERM.
## Distributed under the terms of the GNU General Public License (GPLv3).
## See the LICENSE file for details.
##

#' @importFrom posterior summarize_draws
#' @importFrom posterior subset_draws
NULL

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
      if (!is.matrix(x)) stop_h(
        "Something went wrong while exporting posterior draws."
      )
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
    #' @field param_draws `posterior::draws_df`.
    #' Draws of parameters.
    param_draws = function(value) {
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

    #' @field augm_draws `posterior::draws_df`.
    #' Draws selected individuals augmented variables.
    augm_draws = function(value) {
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

      cli::cli_alert_info(sprintf(
        "Posterior draws of model {.pkg %s}.\n%s",
        spaced_template_choice(private$name_),
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
    },

    #' @description Summarize all posterior draws.
    #' @return `posterior::draws_df`
    summarize_all_draws = function() {
      x1 <- x2 <- x3 <- NULL
      if (posterior::is_draws(private$param_draws_)) {
        x1 <- summarize_draws(private$param_draws_)
      }
      if (posterior::is_draws(private$augm_draws_)) {
        x2 <- summarize_draws(private$augm_draws_)
      }
      if (posterior::is_draws(private$total_log_lik_draws_)) {
        x3 <- summarize_draws(private$total_log_lik_draws_)
      }
      rbind(x1, x2, x3)
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
