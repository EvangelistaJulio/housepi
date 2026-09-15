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
# \insertRef{vehtari_rank_normalization_2021}{housepi}
diagnostic_ <- list(
  # min sample size
  ss_min = 400,
  # max Rhat
  rhat_max = 1.05,
  # ESS ratio
  # NOTE: 10% bayesplot/HMC default seems unrealistic for Metropolis
  neff_min = 0.01
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
  class(x) <- unique(c("housepi_draws", class(x)))
  x
}

# @description Return empty draws.
# @return `posterior::draws_df`
empty_draws_ <- function() {
  tag_as_housepi_draws(posterior::as_draws_df(data.frame()))
}

# Per-draw chain and iteration index implied by an MCMC `setup` list.
draws_index_ <- function(setup) {
  chain_size <- setup_sample_size_(
    steps = setup[["steps"]], burn = setup[["burn"]],
    thin = setup[["thin"]], chains = setup[["chains"]]
  )$chain
  list(chain = rep(seq_len(setup[["chains"]]), each = chain_size),
       iteration = setup[["burn"]] + seq_len(chain_size) * setup[["thin"]])
}

# Tag an Rcpp draws data.frame (rows = draws, cols = variables) as a housepi
# `posterior::draws_df`, attaching the `.chain`/`.iteration`/`.draw` bookkeeping
# from a `draws_index_()` result.
tag_rcpp_draws_ <- function(df, index) {
  if (length(df) == 0) return(NULL)
  df$.chain <- index$chain
  df$.iteration <- rep(index$iteration, times = length(unique(index$chain)))
  df$.draw <- seq_len(nrow(df))
  tag_as_housepi_draws(posterior::as_draws_df(df))
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
  augm <- .x$latent_draws
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
#' This class stores posterior draws of model parameters, the joint log-prior
#' (`lprior`) and the pointwise per-household log-likelihood (`log_lik`).
#' If requested during inference, it also stores latent variables.
#'
#' For models with latent variables the chain targets the joint
#' (latent, parameter) posterior. Write \eqn{x} for a single household's
#' data, \eqn{\theta} for the parameter vector and \eqn{z} for that household's
#' latent variables. The unnormalised joint log-density splits into the two
#' bundles this class exposes,
#' \deqn{\log p(\theta, z, x) =
#'   \underbrace{\log p(\theta)}_{\texttt{lprior}}
#'   + \underbrace{\log p^\mathrm{lat}(z \mid \theta)
#'     + \log p^\mathrm{transm}(x \mid z, \theta)
#'     + \log p^\mathrm{obs}(x \mid z, \theta)}_{\texttt{log\_lik}},}
#' where \eqn{p(\theta)} is the parameter prior, \eqn{p^\mathrm{lat}} the
#' latent-variable model, \eqn{p^\mathrm{transm}} the transmission likelihood
#' and \eqn{p^\mathrm{obs}} the observation likelihood. The transmission term is
#' the contract of the per-household `log_lik`: with no latent variable or
#' observation node the remaining terms vanish and `log_lik` contracts to
#' transmission alone.
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
    param_prior_draws_ = NULL,
    latent_draws_ = NULL,
    point_log_lik_draws_ = NULL,
    draws_index_ = NULL,

    ############################################################################
    ############################################################################
    ############################################################################

    # Store the per-draw chain/iteration index from the inference setup.
    pre_process_setup_ = function(setup) {
      private$draws_index_ <- draws_index_(setup)
      invisible(self)
    },

    # @description Convert draws from housepi Rcpp to `posterior.`
    # @param df `data.frame`. Draws exported from Rcpp core.
    # @return `posterior::draws_df`
    to_posterior_draws_ = function(df) {
      tag_rcpp_draws_(df, private$draws_index_)
    },

    # Stan-style post-inference diagnostics.
    # See https://mc-stan.org/learn-stan/diagnostics-warnings.html
    check_diagnostics_ = function() {
      niter <- posterior::niterations(private$point_log_lik_draws_)
      nchains <- posterior::nchains(private$point_log_lik_draws_)
      ss <- niter * nchains
      if (ss < diagnostic_[["ss_min"]]) {
        warn_h(
          "Total sample size was below {diagnostic_[['ss_min']]}.",
          `!` = "Convergence diagnostics were skipped.",
          "Estimates may be unreliable."
        )
        return(invisible())
      }

      private$check_diagnostics_one_(private$param_draws_, "parameter")
      private$check_diagnostics_one_(private$latent_draws_, "latent variable")
      invisible()
    },

    # Emits at most one warning per draws_kind combination.
    check_diagnostics_one_ = function(d, draws_kind) {
      if (is.null(d) || posterior::nvariables(d) == 0L) return(invisible())

      rhat_max <- diagnostic_[["rhat_max"]]
      neff_min <- diagnostic_[["neff_min"]]
      ss <- posterior::niterations(d) * posterior::nchains(d)
      sm <- posterior::summarise_draws(d)

      # Round R-hat to 2 dp so the check matches the 3-sig-fig value shown by
      # summarise_draws() (e.g. 1.0512 prints as 1.05 and must not be flagged).
      bad_rhat <- which(!is.na(sm$rhat) & round(sm$rhat, 2) > rhat_max)
      bad_bulk <- which(!is.na(sm$ess_bulk) & sm$ess_bulk / ss < neff_min)
      both <- intersect(bad_rhat, bad_bulk)
      only_rhat <- setdiff(bad_rhat, both)
      only_bulk <- setdiff(bad_bulk, both)

      findings <- c()
      if (length(both) > 0L) findings <- c(findings,
        sprintf("%s %s with R-hat > %.2f and Bulk ESS < %.0f%%: %s",
                format_there_were_was(length(both)),
                as_plural(draws_kind, length(both)),
                rhat_max, 100. * neff_min,
                character_to_string_list(sm$variable[both]))
      )
      if (length(only_rhat) > 0L) findings <- c(findings,
        sprintf("%s %s with R-hat > %.2f: %s",
                format_there_were_was(length(only_rhat)),
                as_plural(draws_kind, length(only_rhat)),
                rhat_max,
                character_to_string_list(sm$variable[only_rhat]))
      )
      if (length(only_bulk) > 0L) findings <- c(findings,
        sprintf("%s %s with Bulk Effective Sample Size (ESS) < %.0f%%: %s",
                format_there_were_was(length(only_bulk)),
                as_plural(draws_kind, length(only_bulk)),
                100. * neff_min,
                character_to_string_list(sm$variable[only_bulk]))
      )
      if (length(findings) == 0L) return(invisible())
      names(findings) <- rep("!", length(findings))

      # One hint assembled from the failure modes present: a high R-hat means
      # chains may not have mixed, a low ESS means the sampler mixes slowly.
      hint <- paste(c(
        if (length(bad_rhat) > 0L) "chains may have not mixed",
        if (length(bad_bulk) > 0L) "the sampler explores the posterior slowly"
      ), collapse = "; ")
      substr(hint, 1L, 1L) <- toupper(substr(hint, 1L, 1L))
      warn_h(c(
        "Convergence diagnostics flagged the following:",
        findings, paste0(hint, ".")
      ))
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
      if (is.null(private$param_draws_)) return(empty_draws_())
      private$param_draws_
    },

    #' @field latent_draws `posterior::draws_df`.
    #' Draws of the selected individuals' latent variables.
    latent_draws = function(value) {
      if (!missing(value)) stop_read_only("latent_draws")
      if (is.null(private$latent_draws_)) return(empty_draws_())
      private$latent_draws_
    },

    #' @field log_lik_draws `posterior::draws_df`.
    #' Pointwise per-household log-likelihood
    #' \eqn{
    #'   \log p^\mathrm{obs}(y \mid z, \theta) +
    #'   \log p^\mathrm{transm}(y \mid z, \theta) +
    #'   \log p^\mathrm{lat}(z \mid \theta),
    #' }
    #' the observation, transmission, and latent model contributions.
    #' `r rd_natural_log_`
    log_lik_draws = function(value) {
      if (!missing(value)) stop_read_only("log_lik_draws")
      private$point_log_lik_draws_
    },

    #' @field lprior_draws `posterior::draws_df`.
    #' Single column `lprior`: the parameter log-prior \eqn{\log p(\theta)}. The
    #' latent-variable log-density enters the `log_lik` bundle, not `lprior`.
    #' `r rd_natural_log_`
    lprior_draws = function(value) {
      if (!missing(value)) stop_read_only("lprior_draws")
      has_param <- !is.null(private$param_prior_draws_) &&
        posterior::nvariables(private$param_prior_draws_) > 0
      # with no parameter prior the parameter log-prior is an empty product,
      # log p(theta) = 0; keep a proper zero column over the pointwise draws
      template <- if (has_param) private$param_prior_draws_
      else private$point_log_lik_draws_
      if (is.null(template) || posterior::ndraws(template) == 0) {
        return(empty_draws_())
      }
      lp <- if (has_param) {
        rowSums(posterior::as_draws_matrix(private$param_prior_draws_))
      } else {
        rep(0, posterior::ndraws(template))
      }
      out <- posterior::mutate_variables(template, lprior = lp)
      tag_as_housepi_draws(posterior::subset_draws(out, variable = "lprior"))
    },

    #' @field complete_draws `posterior::draws_df`.
    #' Parameters, pointwise `log_lik` and the joint `lprior`.
    #' Note `log_lik`/`lprior` are evaluated over the joint (latent,
    #' parameter) posterior (see `log_lik_draws` and `lprior_draws`).
    #' `r rd_natural_log_`
    complete_draws = function(value) {
      if (!missing(value)) stop_read_only("complete_draws")
      # rename pointwise log-likelihood standard convention
      aux <- private$point_log_lik_draws_
      n_vars <- posterior::nvariables(aux)
      names(aux)[1:n_vars] <- paste0("log_lik[", 1:n_vars, "]")
      # bind parameters, pointwise log_lik and the parameter-only lprior
      # NOTE: bind_draws strips the tag, so re-stamp
      tag_as_housepi_draws(
        posterior::bind_draws(private$param_draws_, aux, self$lprior_draws)
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
      private$param_draws_ <- private$to_posterior_draws_(
        .engine$flush_param_draws(ptr)
      )
      # per-parameter prior values kept internal: feed `lprior`, not a field
      private$param_prior_draws_ <- private$to_posterior_draws_(
        .engine$flush_param_log_prior_draws(ptr)
      )
      private$latent_draws_ <- private$to_posterior_draws_(
        .engine$flush_latent_draws(ptr)
      )
      # pointwise per-household log-likelihood = transmission (+ observation)
      # (+ latent model, when latents are defined)
      private$point_log_lik_draws_ <- private$to_posterior_draws_(
        .engine$flush_house_log_lik_draws(ptr)
      )

      # free auxiliary memory
      private$draws_index_ <- NULL

      # Stan-style automatic diagnostics
      private$check_diagnostics_()

      invisible()
    },

    #' @description Override the default printing behavior.
    print = function() {
      niter <- posterior::niterations(private$point_log_lik_draws_)
      nchains <- posterior::nchains(private$point_log_lik_draws_)

      if (is.null(private$param_draws_)) nparams <- 0
      else nparams <- posterior::nvariables(private$param_draws_)

      if (is.null(private$latent_draws_)) nlatent <- 0
      else nlatent <- posterior::nvariables(private$latent_draws_)

      model <- spaced_template_choice(private$name_)
      print_bullets(c(
        i = "Posterior draws of model {.val {model}}.",
        paste0(
          "{niter} iteration{?s}, {nchains} chain{?s}",
          ifelse(nparams > 0, ", {nparams} parameter{?s}", ""),
          ifelse(nlatent > 0, ", {nlatent} latent", ""), "."
        )
      ))
      invisible()
    }
  )
)

############################################################################
############################################################################
############################################################################

# Keep the requested variables, all of them when none is named. A name selects
# itself and its indexed elements `name[i]`; `prefix` is a regex the column may
# start with, `(.+_)?` selecting the `<id>_<latent>` traces of a latent name.
subset_draws_ <- function(out, variable, prefix = "") {
  if (is.null(variable)) return(out)
  vars <- posterior::variables(out)
  escaped <- gsub("([][{}()+*^$|\\\\?.])", "\\\\\\1", variable)
  hits <- lapply(sprintf("^%s%s(\\[.+\\])?$", prefix, escaped),
                 grep, x = vars, value = TRUE)
  unmatched <- variable[lengths(hits) == 0L]
  if (length(unmatched) > 0L) stop_h(
    "No draws of variable{?s} {.val {unmatched}}."
  )
  tag_as_housepi_draws(
    posterior::subset_draws(out, variable = unique(unlist(hits)))
  )
}

#' @title Get posterior draws of respective variables
#' @description Return a posterior draws dataframe with the MCMC samples.
#' `get_parameter_draws()` and `get_latent_draws()` warn and return an empty
#' draws object when no such draws were stored.
#' @param hhdraws `HousepiPosteriorDraws` object.
#' Generated by `HousepiModel` inference.
#' @param variable `character()`. Variables to return, all of them when `NULL`.
#' A name also returns its indexed elements, as `log_lik` returns `log_lik[1]`,
#' `log_lik[2]`, ... Latent draws are named `<id>_<latent>`: a latent name
#' returns its draws for every stored individual, a full name only that
#' individual's. For finer selections use [posterior::subset_draws()].
#' @return `posterior::draws_df`
#' For `log_lik`/`lprior` columns: `r rd_natural_log_`
#' @rdname get_draws
#' @export
get_parameter_draws <- function(hhdraws, variable = NULL) {
  stopifnot(inherits(hhdraws, "HousepiPosteriorDraws"))
  out <- hhdraws$parameter_draws
  if (posterior::nvariables(out) == 0) warn_h(
    "Empty draws: no parameter was defined in the model."
  )
  subset_draws_(out, variable)
}

#' @rdname get_draws
#' @export
get_latent_draws <- function(hhdraws, variable = NULL) {
  stopifnot(inherits(hhdraws, "HousepiPosteriorDraws"))
  out <- hhdraws$latent_draws
  if (posterior::nvariables(out) == 0) warn_h(
    "Empty draws: latent variables were not defined or not stored.",
    "Set {.arg id_draws} in {.fn run_inference} to store latent draws."
  )
  subset_draws_(out, variable, prefix = "(.+_)?")
}

#' @rdname get_draws
#' @export
get_log_lik_draws <- function(hhdraws, variable = NULL) {
  stopifnot(inherits(hhdraws, "HousepiPosteriorDraws"))
  subset_draws_(hhdraws$log_lik_draws, variable)
}

#' @rdname get_draws
#' @export
get_complete_draws <- function(hhdraws, variable = NULL) {
  stopifnot(inherits(hhdraws, "HousepiPosteriorDraws"))
  subset_draws_(hhdraws$complete_draws, variable)
}
