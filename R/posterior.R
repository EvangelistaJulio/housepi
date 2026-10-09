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
  neff_min = 0.01,
  # min absolute ESS
  ess_min = 100
)

# What each diagnostic flags, worded once for the warnings and the docs.
diagnostic_flags_ <- list(
  short = sprintf("Total sample size was below %d", diagnostic_[["ss_min"]]),
  rhat = sprintf("R-hat > %.2f", diagnostic_[["rhat_max"]]),
  ess = sprintf(
    "Bulk Effective Sample Size (ESS) below %d or %.0f%% of the draws",
    diagnostic_[["ess_min"]], 100. * diagnostic_[["neff_min"]]
  )
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

# Return the global credible interval level (option "housepi.cri_level").
# Not exported: housepi internals must always pull through this helper.
get_global_cri_level <- function() {
  validate_cri_level(getOption("housepi.cri_level", 0.9))
}

# Lower/upper quantile probabilities corresponding to a credible interval level.
cri_level_probs <- function(cri_level) {
  0.5 * c(lower = 1. - cri_level, upper = 1. + cri_level)
}

# Validate a posterior point estimate choice.
validate_point_estimate <- function(value) {
  if (length(value) != 1 || !is.character(value) || is.na(value)) stop_h(
    "Point estimate must be a single non-NA character value."
  )
  if (!value %in% c("mean", "median")) stop_h(
    "Point estimate must be one between {.val mean} and {.val median}."
  )
  value
}

# Return the global point estimate (option "housepi.point_estimate").
# Not exported: housepi internals must always pull through this helper.
get_global_point_estimate <- function() {
  validate_point_estimate(getOption("housepi.point_estimate", "median"))
}

#' @title Set the global posterior summary
#' @description
#' Set the point estimate and the credible level housepi objects report.
#' An argument left out keeps the value in place.
#' @param point_estimate `character(1)` One between `"mean"` and `"median"`.
#' Default is `"median"`.
#' @param cri_level `numeric(1)` Must lie in `(0.5, 1)`. Default is `0.9`.
#' @return Invisibly returns the previous values, as a `list` the function
#' takes back to restore them.
#' @details
#' A rejected value leaves both settings as they were.
#' @examples
#' old <- set_housepi_summary(point_estimate = "mean")
#' do.call(set_housepi_summary, old)
#' @export
set_housepi_summary <- function(point_estimate = NULL, cri_level = NULL) {
  old <- list(point_estimate = get_global_point_estimate(),
              cri_level = get_global_cri_level())
  # both values are validated before either is written
  new <- c(
    if (!is.null(point_estimate))
      list(housepi.point_estimate = validate_point_estimate(point_estimate)),
    if (!is.null(cri_level))
      list(housepi.cri_level = validate_cri_level(cri_level))
  )
  if (length(new)) options(new)
  invisible(old)
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
  counts <- setup[c("steps", "burn", "thin", "adapt", "chains")]
  chain_size <- do.call(resolve_setup_, counts)$chain
  list(chain = rep(seq_len(setup[["chains"]]), each = chain_size),
       iteration = setup[["burn"]] + seq_len(chain_size) * setup[["thin"]])
}

# Tag an Rcpp draws data.frame (rows = draws, cols = variables) as a housepi
# `posterior::draws_df`, attaching the `.chain`/`.iteration`/`.draw` bookkeeping
# from a `draws_index_()` result.
tag_rcpp_draws_ <- function(df, index) {
  if (!length(df)) return(NULL)
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
#' [set_housepi_summary()]. Otherwise dispatch falls through to
#' [posterior::summarise_draws()] behaviour.
#' @param .x A housepi-tagged draws object, or a [`HousepiPosteriorDraws`].
#' @param ... Summary functions, as in [posterior::summarise_draws()].
#' @param .args,.num_args,.cores Forwarded to [posterior::summarise_draws()].
#' @return A `draws_summary` `data.frame`.
#' @examples
#' results <- flu_1918_frederick |>
#'   new_housepi_model("final-size") |>
#'   add_parameter(import = "template") |>
#'   run_inference(chains = 2, progress = FALSE)
#'
#' # all parameters and stored latent variables
#' summarise_draws(results)
#' # a subset of the draws
#' summarise_draws(get_parameter_draws(results, "alpha"))
#' @name summarise_draws.housepi_draws
#' @exportS3Method posterior::summarise_draws
summarise_draws.housepi_draws <- function(
  .x, ..., .args = list(),
  .num_args = getOption("posterior.num_args", list()),
  .cores = 1
) {
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
  .cores = 1
) {
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

# Per-variable R-hat and Bulk ESS of a draws block over a sample of size `ss`,
# and which of them fail: `NA` where a diagnostic is undefined.
diagnose_block_ <- function(d, ss) {
  if (is.null(d) || posterior::nvariables(d) == 0L) return(NULL)
  # a capped ESS on a short sample is reported as a low ESS finding
  sm <- as.data.frame(suppressWarnings(posterior::summarise_draws(
    d, rhat = posterior::rhat, ess_bulk = posterior::ess_bulk
  )))
  # Round R-hat to 2 dp so the check matches the 3-sig-fig value shown by
  # summarise_draws() (e.g. 1.0512 prints as 1.05 and must not be flagged).
  sm$bad_rhat <- round(sm$rhat, 2) > diagnostic_[["rhat_max"]]
  sm$bad_ess <- sm$ess_bulk <
    max(diagnostic_[["ess_min"]], diagnostic_[["neff_min"]] * ss)
  sm
}

# Findings of one diagnosed block, under a header naming the flagged draws,
# one line per failure combination, with a hint on the failure modes present;
# `NULL` when nothing is flagged.
check_diagnostics_one_ <- function(sm, draws_kind) {
  if (is.null(sm)) return(NULL)
  bad_rhat <- which(sm$bad_rhat %in% TRUE)
  bad_bulk <- which(sm$bad_ess %in% TRUE)
  both <- intersect(bad_rhat, bad_bulk)
  only_rhat <- setdiff(bad_rhat, both)
  only_bulk <- setdiff(bad_bulk, both)

  finding <- function(idx, test) {
    n <- length(idx)
    if (n > 0L) inline_cli_(
      "{cli::qty(n)}There {?was/were} {n} {draws_kind}{cli::qty(n)}{?s}",
      " with {test}: {.val {sm$variable[idx]}}."
    )
  }
  findings <- c(
    finding(both, paste(diagnostic_flags_$rhat, "and", diagnostic_flags_$ess)),
    finding(only_rhat, diagnostic_flags_$rhat),
    finding(only_bulk, diagnostic_flags_$ess)
  )
  if (!length(findings)) return(NULL)
  names(findings) <- rep("!", length(findings))

  # a high R-hat means chains may not have mixed,
  # a low ESS means the sampler mixes slowly
  hint <- paste(c(
    if (length(bad_rhat)) "chains may have not mixed",
    if (length(bad_bulk)) "the sampler explores the posterior slowly"
  ), collapse = "; ")
  substr(hint, 1L, 1L) <- toupper(substr(hint, 1L, 1L))
  list(header = paste0("Convergence diagnostics flagged ", draws_kind,
                       " draws."),
       findings = findings, hint = paste0(hint, "."))
}

# Stan-style post-inference diagnostics: the message blocks to report and the
# per-parameter R-hat, Bulk ESS and whether the parameter is flagged.
# A short sample reports its findings under one block of its own.
# See https://mc-stan.org/learn-stan/diagnostics-warnings.html
check_diagnostics_ <- function(draws) {
  niter <- posterior::niterations(draws$point_log_lik)
  nchains <- posterior::nchains(draws$point_log_lik)
  ss <- niter * nchains
  short <- ss < diagnostic_[["ss_min"]]

  sm <- lapply(draws[c("param", "latent")], diagnose_block_, ss = ss)
  found <- Filter(Negate(is.null), list(
    check_diagnostics_one_(sm$param, "parameter"),
    check_diagnostics_one_(sm$latent, "latent variable")
  ))
  # the header is interpolated here, so it is recorded as literal text
  blocks <- if (short) list(c(
    paste0(diagnostic_flags_$short, "."),
    unlist(lapply(found, `[[`, "findings")),
    "Estimates may be unreliable."
  )) else lapply(found, function(f) {
    c(f$header, f$findings, f$hint)
  })

  param <- sm$param
  if (!is.null(param)) param <- data.frame(
    param[c("variable", "rhat", "ess_bulk")],
    flagged = short | !param$bad_rhat %in% FALSE | !param$bad_ess %in% FALSE
  )
  list(blocks = blocks, param = param)
}

# Distinct headers of a list of message blocks, with their counts.
# The header states the defect; its bullets are detail.
count_diagnostics_ <- function(blocks) {
  headers <- vapply(blocks, `[[`, character(1L), 1L)
  counts <- table(headers)
  data.frame(message = as.character(names(counts)),
             count = as.integer(counts),
             row.names = NULL, stringsAsFactors = FALSE)
}

#' @title Posterior draws from housepi model inference
#' @description
#' Posterior draws returned by [run_inference()] on a household model.
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
#' @examples
#' results <- flu_1918_frederick |>
#'   new_housepi_model("final-size") |>
#'   add_parameter(import = "template") |>
#'   run_inference(chains = 2, progress = FALSE)
#' results
#' results$parameter_draws
#' @docType class
HousepiPosteriorDraws <- R6::R6Class(
  classname = "HousepiPosteriorDraws",

  ############################################################################
  ################################## PRIVATE #################################
  ############################################################################

  private = list(
    # @field
    name_ = NA_character_,
    draws_ = list(),
    secs_ = NA_real_,
    warnings_ = NULL,
    # per-parameter R-hat, Bulk ESS and convergence verdict
    diagnostics_ = NULL
  ),

  ############################################################################
  ################################## ACTIVE ##################################
  ############################################################################

  active = list(
    #' @field parameter_draws `posterior::draws_df`.
    #' Draws of parameters.
    parameter_draws = function(value) {
      if (!missing(value)) stop_read_only("param_draws")
      if (is.null(private$draws_$param)) return(empty_draws_())
      private$draws_$param
    },

    #' @field latent_draws `posterior::draws_df`.
    #' Draws of the selected individuals' latent variables.
    latent_draws = function(value) {
      if (!missing(value)) stop_read_only("latent_draws")
      if (is.null(private$draws_$latent)) return(empty_draws_())
      private$draws_$latent
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
      private$draws_$point_log_lik
    },

    #' @field lprior_draws `posterior::draws_df`.
    #' Single column `lprior`: the parameter log-prior \eqn{\log p(\theta)}. The
    #' latent-variable log-density enters the `log_lik` bundle, not `lprior`.
    #' `r rd_natural_log_`
    lprior_draws = function(value) {
      if (!missing(value)) stop_read_only("lprior_draws")
      if (is.null(private$draws_$lprior)) return(empty_draws_())
      private$draws_$lprior
    },

    #' @field complete_draws `posterior::draws_df`.
    #' Parameters, pointwise `log_lik` and the joint `lprior`.
    #' Note `log_lik`/`lprior` are evaluated over the joint (latent,
    #' parameter) posterior (see `log_lik_draws` and `lprior_draws`).
    #' `r rd_natural_log_`
    complete_draws = function(value) {
      if (!missing(value)) stop_read_only("complete_draws")
      # rename pointwise log-likelihood standard convention
      aux <- private$draws_$point_log_lik
      n_vars <- posterior::nvariables(aux)
      names(aux)[1:n_vars] <- paste0("log_lik[", 1:n_vars, "]")
      # bind parameters, pointwise log_lik and the parameter-only lprior
      # NOTE: bind_draws strips the tag, so re-stamp
      tag_as_housepi_draws(
        posterior::bind_draws(private$draws_$param, aux, private$draws_$lprior)
      )
    },

    #' @field warnings `data.frame`
    #' Convergence diagnostics the run flagged
    #' \describe{
    #'   \item{message}{Line stating what was flagged}
    #'   \item{count}{Times it was raised}
    #' }
    #' Recorded whether or not the run reported them. Per-variable R-hat and
    #' Bulk ESS behind a flag are read from [summarise_draws()].
    warnings = function(value) {
      if (!missing(value)) stop_read_only("warnings")
      private$warnings_
    },

    #' @field inference_secs `numeric(1)`
    #' Wall-clock of the sampling loop, leaving out the model build and the
    #' conversion of the draws.
    inference_secs = function(value) {
      if (!missing(value)) stop_read_only("inference_secs")
      private$secs_
    }
  ),

  ############################################################################
  ################################## PUBLIC ##################################
  ############################################################################

  public = list(
    #' @description
    #' Internal: [run_inference()] builds this object,
    #' don't call `$new()` directly.
    #' @param name Model name.
    #' @param ptr Internal engine handle.
    #' @param secs Wall-clock the engine timed the sampling loop at.
    initialize = function(name, ptr, secs = NA_real_) {
      # inherit model name
      private$name_ <- as.character(name)
      private$secs_ <- as.numeric(secs)
      # per-draw chain/iteration index from the inference setup
      idx <- draws_index_(.engine$export_setup(ptr))

      # store cpp draws converted to posterior
      tag <- \(df) tag_rcpp_draws_(df, idx)
      lp <- .engine$flush_param_log_prior_draws(ptr)
      private$draws_ <- list(
        param = tag(.engine$flush_param_draws(ptr)),
        lprior = tag(if (length(lp)) data.frame(lprior = rowSums(lp))),
        latent = tag(.engine$flush_latent_draws(ptr)),
        point_log_lik = tag(.engine$flush_house_log_lik_draws(ptr))
      )

      # Stan-style automatic diagnostics, recorded and reported alike:
      # a quiet run holds the warnings back, the record stands either way
      diag <- check_diagnostics_(private$draws_)
      private$diagnostics_ <- diag$param
      private$warnings_ <- count_diagnostics_(diag$blocks)
      for (block in diag$blocks) warn_h(block)

      invisible()
    },

    #' @description Override the default printing behavior.
    print = function() {
      niter <- posterior::niterations(private$draws_$point_log_lik)
      nchains <- posterior::nchains(private$draws_$point_log_lik)

      if (is.null(private$draws_$param)) nparams <- 0
      else nparams <- posterior::nvariables(private$draws_$param)

      if (is.null(private$draws_$latent)) nlatent <- 0
      else nlatent <- posterior::nvariables(private$draws_$latent)

      model <- spaced_choice(private$name_)
      print_bullets(c(
        i = "Posterior draws of model {.val {model}}.",
        paste0(
          "{niter} iteration{?s}, {nchains} chain{?s}",
          if (nparams > 0) ", {nparams} parameter{?s}" else "",
          if (nlatent > 0) ", {nlatent} latent" else "", "."
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
  if (length(unmatched)) stop_h(
    "No draws of variable{?s} {.val {unmatched}}."
  )
  tag_as_housepi_draws(
    posterior::subset_draws(out, variable = unique(unlist(hits)))
  )
}

# Ensure object of class HousepiPosteriorDraws.
stopifnotdraws <- function(hhdraws) {
  if (!inherits(hhdraws, "HousepiPosteriorDraws")) stop_h(
    "{.arg hhdraws} must be of class {.cls HousepiPosteriorDraws}.",
    "Create one with {.fn new_housepi_model}."
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
#' `log_lik[2]`, ... from `get_complete_draws()`. `get_log_lik_draws()` names
#' its columns by `hh_id`: a household id returns that household's draws.
#' Latent draws are named `<id>_<latent>`: a latent name
#' returns its draws for every stored individual, a full name only that
#' individual's. For finer selections use [posterior::subset_draws()].
#' @return `posterior::draws_df`
#' For `log_lik`/`lprior` columns: `r rd_natural_log_`
#' @examples
#' results <- flu_1918_frederick |>
#'   new_housepi_model("final-size") |>
#'   add_parameter(import = "template") |>
#'   run_inference(chains = 2, progress = FALSE)
#' get_parameter_draws(results, c("alpha", "beta"))
#'
#' # parameters, pointwise log_lik and lprior side by side
#' complete <- get_complete_draws(results)
#' posterior::variables(complete) |> head()
#' @rdname get_draws
#' @export
get_parameter_draws <- function(hhdraws, variable = NULL) {
  stopifnotdraws(hhdraws)
  out <- hhdraws$parameter_draws
  if (posterior::nvariables(out) == 0) warn_h(
    "Empty draws: no parameter was defined in the model."
  )
  subset_draws_(out, variable)
}

#' @rdname get_draws
#' @export
get_latent_draws <- function(hhdraws, variable = NULL) {
  stopifnotdraws(hhdraws)
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
  stopifnotdraws(hhdraws)
  subset_draws_(hhdraws$log_lik_draws, variable)
}

#' @rdname get_draws
#' @export
get_complete_draws <- function(hhdraws, variable = NULL) {
  stopifnotdraws(hhdraws)
  subset_draws_(hhdraws$complete_draws, variable)
}
