#      __
#     / /_  ____  __  __ ____ __  ____  __
#    / __ \/ __ \/ / / / ___/ _ \/ __ \/ /
#   / / / / /_/ / /_/ (__  )  __/ /_/ / /
#  /_/ /_/\____/\____/____/\___/ ____/_/
#                             /_/
#
# Copyright (C) 2025 CAMAGLIA Francesco - Institut Pasteur, Paris (France).
#
# posterior.R
#

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
HousepiPosteriorDraws <- R6::R6Class(
  classname = "HousepiPosteriorDraws",

  ############################################################################
  ################################## PRIVATE #################################
  ############################################################################

  private = list(
    # @field
    .name = NA_character_,
    .param_df = NULL,
    .prior_df = NULL,
    .augm_df = NULL,
    .full_log_lik_df = NULL,
    .point_log_lik_df = NULL,

    # split into a list of matrix
    .split_draws_by_chain = function(x) {
      if (length(x) == 0) stop_h("Traces are empty.")
      lapply(
        split(seq_len(nrow(x)), x[, ".chain"]),
        function(idx) x[idx, colnames(x) != ".chain", drop = FALSE]
      )
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
      # nothing to be done
      if (posterior::is_draws(x)) return(x)

      if (is.matrix(x) || (is.data.frame(x))) {
        list_trs <- private$.split_draws_by_chain(x)

        # create empty array
        n_chain <- length(list_trs)
        par_names <- colnames(list_trs[[1]])
        n_iter <- nrow(list_trs[[1]])
        iter_names <- rownames(list_trs[[1]])
        array_trs <- array(
          NA_real_, dim = c(n_iter, n_chain, length(par_names)),
          dimnames = list(iter_names, seq_len(n_chain), par_names)
        )

        # associate draws from Cpp
        for (i in seq_along(list_trs)) array_trs[, i, ] <- list_trs[[i]]
        output <- posterior::as_draws_df(array_trs)
        output$.iteration <- as.integer(rownames(x))
        return(output)
      }

      stop_h(
        "Unrecognized type for 'Rcpp_ModelHousepi' draws."
      )
    },

    # @description
    # Convert draws export from `Rcpp_ModelHousepi` to `coda.`
    #
    # @param trs
    # Traces as exported from Rcpp `housepi` model.
    #
    # @return `coda::mcmc` | `coda::mcmc.list`
    .cpp_draws_to_coda = function(x) {
      require_suggested_library("coda")
      if ((coda::is.mcmc(x) || coda::is.mcmc.list(x))) {
        return(x)

      } else if (is.matrix(x) || (is.data.frame(x))) {
        list_trs <- split_draws_by_chain_(x)
        sweeps <- as.integer(rownames(list_trs[[1]]))
        list_coda_mcmc <- lapply(list_trs, function(x) {
          coda::mcmc(
            data = x,
            start = sweeps[1],
            end = sweeps[length(sweeps)],
            thin = sweeps[2] - sweeps[1]
          )
        })
        if (length(list_coda_mcmc) == 1) return(list_coda_mcmc[[1]])
        return(coda::as.mcmc.list(list_coda_mcmc))
      }

      stop_h(
        "Unrecognized type for 'Rcpp_ModelHousepi' draws export."
      )
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
      private$.param_df
    },

    #' @field prior_draws `posterior::draws_df`
    #' Posterior draws data.frame of values of parameter priors.
    prior_draws = function(value) {
      if (!missing(value)) stop_h("'prior_draws' is read only.")
      private$.prior_df
    },

    #' @field augm_draws `posterior::draws_df`
    #' Posterior draws data.frame of selected individuals augmented variables.
    augm_draws = function(value) {
      if (!missing(value)) stop_h("'augm_draws' is read only.")
      private$.augm_df
    },

    #' @field log_lik_draws `posterior::draws_df`
    #' Posterior draws data.frame of households pointwise loglikelihood.
    log_lik_draws = function(value) {
      if (!missing(value)) stop_h("'log_lik_draws' is read only.")
      private$.point_log_lik_df
    },

    #' @field total_log_lik_draws `posterior::draws_df`
    #' Posterior draws data.frame of the total loglikelihood.
    total_log_lik_draws = function(value) {
      if (!missing(value)) stop_h("'total_log_lik_draws' is read only.")
      private$.full_log_lik_df
    },

    #' @field complete_draws `posterior::draws_df`
    #' Parameters draws complete of pointwise log-likelihood and priors.
    complete_draws <- function(value) {
      if (!missing(value)) stop_h("'complete_draws' is read only.")
      # rename pointwise loglikelihood
      aux <- private$.point_log_lik_df
      n_vars <- posterior::nvariables(aux)
      names(aux)[1:n_vars] <- paste0("log_lik[", 1:n_vars, "]")
      # bind all together
      posterior::mutate_variables(
        posterior::bind_draws(private$.param_draws, aux, private$.prior_draws),
        lprior = rowSums(posterior::as_draws_matrix(private$.prior_draws))
      )
    }
  ),

  ############################################################################
  ################################## PUBLIC ##################################
  ############################################################################

  public = list(
    #' store posterior::draws_df internally
    #' @param name `character(1)` Model name.
    #' @param rcpp_object `Rcpp_ModelHousepi`
    #' @noRd
    initialize = function(name, rcpp_object) {
      private$.name <- as.character(name)
      all_traces <- private$.cpp_draws_to_posterior(rcpp_object$export_draws())
      vars <- posterior::variables(all_traces)

      # augmented draws
      id_draws <- rcpp_object$export_setup()$id_draws

      if (length(id_draws) > 0) {
        mask <- apply(sapply(id_draws,
                             function(x) stringr::str_detect(vars, x)),
                      1, base::any)
        augm_vars <- vars[mask]
        private$.augm_df <- posterior::subset_draws(all_traces,
                                                    variable = augm_vars)
        private$.augm_df$.iteration <- all_traces$.iteration
        vars <- setdiff(vars, augm_vars)
      }

      # loglikelihood draws
      private$.full_log_lik_df <- posterior::subset_draws(
        all_traces, variable = c("log_lik")
      )
      private$.full_log_lik_df$.iteration <- all_traces$.iteration
      vars <- setdiff(vars, c("log_lik"))

      # prior parameters draws
      prior_vars <- vars[stringr::str_detect(vars, "lprior_")]
      private$.prior_df <- posterior::subset_draws(all_traces,
                                                   variable = prior_vars)
      private$.prior_df$.iteration <- all_traces$.iteration
      vars <- setdiff(vars, prior_vars)

      # parameters draws
      private$.param_df <- posterior::subset_draws(all_traces, variable = vars)
      private$.param_df$.iteration <- all_traces$.iteration

      # household pointwise loglikelihood
      private$.point_log_lik_df <- private$.cpp_draws_to_posterior(
        rcpp_object$export_house_loglike_draws()
      )

      # release memory in original object
      rcpp_object$clear_memory()
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
          posterior::niterations(self$total_log_lik_draws),
          posterior::nchains(self$total_log_lik_draws),
          ifelse(is.null(self$param_draws), "",
                 sprintf(", %s parameters",
                         posterior::nvariables(self$param_draws))),
          ifelse(is.null(self$augm_draws), "",
                 sprintf(", %s augmented",
                         posterior::nvariables(self$augm_draws)))
        )
      ))
      invisible()
    },

    #' @description
    #' Summarize all posterior draws.
    #' @return `posterior::draws_df`
    summarize_all_draws = function() {
      x1 <- x2 <- x3 <- NULL
      if (posterior::is_draws(self$param_draws)) {
        x1 <- posterior::summarize_draws(self$param_draws)
      }
      if (posterior::is_draws(self$augm_draws)) {
        x2 <- posterior::summarize_draws(self$augm_draws)
      }
      if (posterior::is_draws(self$total_log_lik_draws)) {
        x3 <- posterior::summarize_draws(self$total_log_lik_draws)
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
