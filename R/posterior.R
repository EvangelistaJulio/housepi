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
#' @keywords internal
HousepiPosteriorDraws <- R6::R6Class(
  classname = "HousepiPosteriorDraws",

  ############################################################################
  ################################## PRIVATE #################################
  ############################################################################

  private = list(
    # @field
    .param_df = NULL,
    .augm_df = NULL,
    .loglike_df = NULL,
    .house_df = NULL,

    # split into a list of matrix
    .split_draws_by_chain = function(x) {
      if (length(x) == 0) stop_h("Traces are empty.")
      lapply(
        split(seq_len(nrow(x)), x[, "..chain.."]),
        function(idx) x[idx, colnames(x) != "..chain..", drop = FALSE]
      )
    },

    # @description
    # Convert draws export from `Rcpp_ModelHousepi` to `posterior.`
    #
    # @param trs `matrix` | `data.frame`\cr
    # Traces as exported from Rcpp core.
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
      if (missing(value)) private$.param_df
      else stop_h("'param_draws' is read only.")
    },

    #' @field augm_draws `posterior::draws_df`
    #' Posterior draws data.frame of selected individuals augmented variables.
    augm_draws = function(value) {
      if (missing(value)) private$.augm_df
      else stop_h("'augm_draws' is read only.")
    },

    #' @field loglike_draws `posterior::draws_df`
    #' Posterior draws data.frame of full loglikelihood.
    loglike_draws = function(value) {
      if (missing(value)) private$.loglike_df
      else stop_h("'loglike_draws' is read only.")
    },

    #' @field house_draws `posterior::draws_df`
    #' Posterior draws data.frame of households pointwise loglikelihood.
    house_draws = function(value) {
      if (missing(value)) private$.house_df
      else stop_h("'house_draws' is read only.")
    }
  ),

  ############################################################################
  ################################## PUBLIC ##################################
  ############################################################################

  public = list(
    #' store posterior::draws_df internally
    #' @param rcpp_object `Rcpp_ModelHousepi`
    #' @noRd
    initialize = function(rcpp_object) {
      all_traces <- private$.cpp_draws_to_posterior(rcpp_object$export_draws())
      vars <- posterior::variables(all_traces)

      # augmented draws
      id_draws <- rcpp_object$export_setup()$id_draws

      if (length(id_draws) > 0) {
        augm_vars <- vars[apply(
          sapply(id_draws, function(x) stringr::str_detect(vars, x)),
          1, any
        )]
        private$.augm_df <- posterior::subset_draws(all_traces, 
                                                    variable = augm_vars)
        vars <- setdiff(vars, augm_vars)
      }

      # loglikelihood draws
      private$.loglike_df <- posterior::subset_draws(all_traces, 
                                                     variable = c("loglikelihood"))
      vars <- setdiff(vars, c("loglikelihood"))

      # parameters draws
      private$.param_df <- posterior::subset_draws(all_traces, variable = vars)

      # household pointwise loglikelihood
      private$.house_df <- private$.cpp_draws_to_posterior(
        rcpp_object$export_house_loglike_draws()
      )

      # release memory in original object
      rcpp_object$clear_memory()
      invisible()
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
#' @param object `HousepiPosteriorDraws` object
#' as generated by `HousepiModel` inference.
#'
#' @return `posterior::draws_df`
#'
#' @rdname get_draws
#' @export
get_param_draws <- function(object) {
  stopifnot(inherits(object, "HousepiPosteriorDraws"))
  object$param_draws
}

#' @rdname get_draws
#' @export
get_augm_draws <- function(object) {
  stopifnot(inherits(object, "HousepiPosteriorDraws"))
  object$augm_draws
}

#' @rdname get_draws
#' @export
get_loglike_draws <- function(object) {
  stopifnot(inherits(object, "HousepiPosteriorDraws"))
  object$loglike_draws
}

#' @rdname get_draws
#' @export
get_house_loglike_draws <- function(object) {
  stopifnot(inherits(object, "HousepiPosteriorDraws"))
  object$house_draws
}