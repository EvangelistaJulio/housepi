## housepi - Epidemiological Analysis of Household Outbreaks.
## Copyright 2024-2026 Institut Pasteur (Paris), France.
## SPDX-License-Identifier: GPL-3.0-or-later
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
  # R-hat flagged at or above
  rhat_bound = 1.01,
  # min Bulk and Tail ESS per chain: 50 per split chain
  ess_chain_min = 100
)

# What each diagnostic flags, worded once for the warnings and the docs.
diagnostic_flags_ <- list(
  rhat = sprintf("R-hat of %.2f or more, or undefined",
                 diagnostic_[["rhat_bound"]]),
  ess = sprintf(
    "Bulk or Tail Effective Sample Size (ESS) below %d per chain, or undefined",
    diagnostic_[["ess_chain_min"]]
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

# Global credible interval level, option "housepi.cri_level".
# housepi internals read it through this helper only.
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
    "Point estimate must be one of {.val mean} and {.val median}."
  )
  value
}

# Global point estimate, option "housepi.point_estimate".
# housepi internals read it through this helper only.
get_global_point_estimate <- function() {
  validate_point_estimate(getOption("housepi.point_estimate", "median"))
}

#' @title Set the global posterior summary
#' @description
#' Set the point estimate and the credible level housepi objects report.
#' An argument left out keeps the value in place.
#' @param point_estimate `character(1)` One of `"mean"` and `"median"`.
#' The package default is `"median"`.
#' @param cri_level `numeric(1)` Between 0.5 and 1, both excluded.
#' The package default is `0.9`.
#' @return `list` The previous values, invisibly,
#' which the function takes back to restore them.
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

# Empty `posterior::draws_df`.
empty_draws_ <- function() {
  tag_as_housepi_draws(posterior::as_draws_df(data.frame()))
}

# Stored draws of one block of a `HousepiPosteriorDraws`, empty when absent.
# @param block `character(1)` `"param"`, `"latent"`, `"point_log_lik"`
# or `"lprior"`.
draws_block_ <- function(hh_fit, block) {
  out <- hh_fit$.__enclos_env__$private$draws_[[block]]
  if (is.null(out)) empty_draws_() else out
}

# Per-draw chain index implied by an MCMC `settings` list.
draws_chain_ <- function(settings) {
  chain_size <- settings[["steps"]] %/% settings[["thin"]]
  rep(seq_len(settings[["chains"]]), each = chain_size)
}

# Tag an Rcpp draws data.frame, draws by variables,
# as a housepi `posterior::draws_df`,
# with `.chain` from a draws_chain_() result,
# `.iteration` and `.draw` numbered by posterior.
tag_rcpp_draws_ <- function(df, chain) {
  if (!length(df)) return(NULL)
  df$.chain <- chain
  tag_as_housepi_draws(posterior::as_draws_df(df))
}

#' @title Summarise housepi posterior draws using the global credible interval
#' @description
#' S3 method on [posterior::summarise_draws()] for objects tagged
#' `"housepi_draws"`, returned by housepi getters and by the active fields
#' of [`HousepiPosteriorDraws`].
#' Without summary functions in `...`, the housepi defaults are used,
#' and the quantile pair tracks [set_housepi_summary()].
#' Otherwise the call falls through to [posterior::summarise_draws()].
#' @param .x A housepi-tagged draws object, or a [`HousepiPosteriorDraws`].
#' @param ... Summary functions, as in [posterior::summarise_draws()].
#' @param .args,.num_args,.cores Forwarded to [posterior::summarise_draws()].
#' @return A `draws_summary` `data.frame`.
#' @examples
#' hh_fit <- flu_1918_frederick |>
#'   new_housepi_model("final-size") |>
#'   add_parameter(import = "template") |>
#'   run_inference()
#'
#' # all parameters and stored latent variables
#' summarise_draws(hh_fit)
#' # a subset of the draws
#' summarise_draws(get_parameter_draws(hh_fit, "alpha"))
#' @name summarise_draws.housepi_draws
#' @exportS3Method posterior::summarise_draws
summarise_draws.housepi_draws <- function(
  .x, ..., .args = list(),
  .num_args = getOption("posterior.num_args", list()),
  .cores = 1
) {
  # strip the tag so the forwarded call dispatches to the next class
  class(.x) <- setdiff(class(.x), "housepi_draws")
  if (length(list(...)) > 0)
    return(posterior::summarise_draws(.x, ..., .args = .args,
                                      .num_args = .num_args, .cores = .cores))

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
  params <- draws_block_(.x, "param")
  augm <- draws_block_(.x, "latent")
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

########  HousepiPosteriorDraws  ########

# Per-variable R-hat, Bulk and Tail ESS of a draws block over `nchains` chains,
# and which of them fail.
# * undefined diagnostic -> failed.
diagnose_block_ <- function(d, nchains) {
  if (is.null(d) || posterior::nvariables(d) == 0L) return(NULL)
  # a capped ESS on a short sample is reported as a low ESS finding
  sm <- as.data.frame(suppressWarnings(posterior::summarise_draws(
    d, rhat = posterior::rhat, ess_bulk = posterior::ess_bulk,
    ess_tail = posterior::ess_tail
  )))
  sm$bad_rhat <- is.na(sm$rhat) | sm$rhat >= diagnostic_[["rhat_bound"]]
  ess_min <- diagnostic_[["ess_chain_min"]] * nchains
  sm$bad_ess <- is.na(sm$ess_bulk) | is.na(sm$ess_tail) |
    sm$ess_bulk < ess_min | sm$ess_tail < ess_min
  sm
}

# Findings of one diagnosed block, under a header naming the flagged draws,
# one line per failure combination, with a hint on the failure modes present.
# * nothing flagged -> `NULL`.
check_diagnostics_one_ <- function(sm, draws_kind) {
  if (is.null(sm)) return(NULL)
  bad_rhat <- which(sm$bad_rhat)
  bad_bulk <- which(sm$bad_ess)
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

# Stan-style post-inference diagnostics: the message blocks to report,
# and the per-parameter R-hat, Bulk and Tail ESS and flag.
# See https://mc-stan.org/learn-stan/diagnostics-warnings.html
check_diagnostics_ <- function(draws) {
  nchains <- posterior::nchains(draws$point_log_lik)
  sm <- lapply(draws[c("param", "latent")], diagnose_block_, nchains = nchains)
  found <- Filter(Negate(is.null), list(
    check_diagnostics_one_(sm$param, "parameter"),
    check_diagnostics_one_(sm$latent, "latent variable")
  ))
  # the header is interpolated here, so it is recorded as literal text
  blocks <- lapply(found, function(f) c(f$header, f$findings, f$hint))

  param <- sm$param
  if (!is.null(param)) param <- data.frame(
    param[c("variable", "rhat", "ess_bulk", "ess_tail")],
    flagged = param$bad_rhat | param$bad_ess
  )
  list(blocks = blocks, param = param)
}

# Distinct headers of a list of message blocks, with their counts.
# The header states the defect, its bullets the detail.
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
#' This class stores posterior draws of model parameters,
#' the joint log-prior (`lprior`)
#' and the pointwise per-household log-likelihood (`log_lik`).
#' If requested during inference, it also stores latent variables.
#' Read them with [get_parameter_draws()], [get_latent_draws()],
#' [get_log_lik_draws()] and [get_complete_draws()].
#'
#' For models with latent variables,
#' the chain targets the joint (latent, parameter) posterior.
#' Write \eqn{x} for a single household's data,
#' \eqn{\theta} for the parameter vector
#' and \eqn{z} for that household's latent variables.
#' The unnormalised joint log-density splits into the two bundles
#' this class exposes,
#' \deqn{\log p(\theta, z, x) =
#'   \underbrace{\log p(\theta)}_{\texttt{lprior}}
#'   + \underbrace{\log p^\mathrm{lat}(z \mid \theta)
#'     + \log p^\mathrm{transm}(x \mid z, \theta)
#'     + \log p^\mathrm{obs}(x \mid z, \theta)}_{\texttt{log\_lik}},}
#' where \eqn{p(\theta)} is the parameter prior,
#' \eqn{p^\mathrm{lat}} the latent-variable model,
#' \eqn{p^\mathrm{transm}} the transmission likelihood
#' and \eqn{p^\mathrm{obs}} the observation likelihood.
#' The transmission term is the contract of the per-household `log_lik`:
#' with no latent variable or observation node, the remaining terms vanish
#' and `log_lik` contracts to transmission alone.
#' @examples
#' hh_fit <- flu_1918_frederick |>
#'   new_housepi_model("final-size") |>
#'   add_parameter(import = "template") |>
#'   run_inference()
#' hh_fit
#' get_parameter_draws(hh_fit)
#' @docType class
HousepiPosteriorDraws <- R6::R6Class(
  classname = "HousepiPosteriorDraws",

  ####  Private  ####

  private = list(
    name_ = NA_character_,
    draws_ = list(),
    secs_ = NA_real_,
    # engine settings of the run, seed, threads and chains aside
    settings_ = NULL,
    # engine, and hashes of the data and definitions the run was fitted on
    fingerprint_ = NULL,
    # engine seed of each run the chains come from
    seed_ = NA_real_,
    warnings_ = NULL,
    # per-parameter R-hat, Bulk and Tail ESS and convergence verdict
    diagnostics_ = NULL,

    # Record Stan-style diagnostics of the draws, then report them.
    # * quiet run -> warnings held back, record kept.
    diagnose_ = function() {
      diag <- check_diagnostics_(private$draws_)
      private$diagnostics_ <- diag$param
      private$warnings_ <- count_diagnostics_(diag$blocks)
      for (block in diag$blocks) {
        warn_h(block, .class = "housepi_warning_diagnostics")
      }
    }
  ),

  ####  Active  ####

  active = list(
    #' @field warnings `r rd_warnings_`
    #' They hold the convergence diagnostics the run flagged.
    #' Per-variable R-hat, Bulk and Tail ESS behind a flag
    #' are read from [summarise_draws()].
    warnings = function(value) {
      assert_read_only_(value, "warnings")
      private$warnings_
    },

    #' @field inference_secs `numeric(1)`
    #' Wall-clock of the sampling loop,
    #' leaving out the model build and the conversion of the draws.
    inference_secs = function(value) {
      assert_read_only_(value, "inference_secs")
      private$secs_
    }
  ),

  ####  Public  ####

  public = list(
    #' @description
    #' Create the posterior draws object.
    #' Use [run_inference()] instead of calling `$new()`.
    #' @param name Model name.
    #' @param ptr Internal engine handle.
    #' @param secs Wall-clock the engine timed the sampling loop at.
    #' @param fingerprint Engine, and hashes of the data and definitions fitted.
    initialize = function(name, ptr, secs = NA_real_, fingerprint = NULL) {
      # inherit model name
      private$name_ <- as.character(name)
      private$secs_ <- as.numeric(secs)
      private$fingerprint_ <- fingerprint
      settings <- .engine$export_settings(ptr)
      private$seed_ <- settings$seed
      private$settings_ <-
        settings[setdiff(names(settings), c("seed", "threads", "chains"))]
      # per-draw chain index from the inference settings
      chain <- draws_chain_(settings)

      # store cpp draws converted to posterior
      tag <- \(df) tag_rcpp_draws_(df, chain)
      lp <- .engine$flush_param_log_prior_draws(ptr)
      # parameter columns in the order of get_parameters()
      param <- .engine$flush_param_draws(ptr)
      param <- param[sort_definitions_(data.frame(name = names(param)))$name]
      private$draws_ <- list(
        param = tag(param),
        lprior = tag(if (length(lp)) data.frame(lprior = rowSums(lp))),
        latent = tag(.engine$flush_latent_draws(ptr)),
        point_log_lik = tag(.engine$flush_house_log_lik_draws(ptr))
      )
      private$diagnose_()
      invisible()
    },

    #' @description Print a summary of the run.
    print = function() {
      niter <- posterior::niterations(private$draws_$point_log_lik)
      nchains <- posterior::nchains(private$draws_$point_log_lik)

      if (is.null(private$draws_$param)) nparams <- 0
      else nparams <- posterior::nvariables(private$draws_$param)

      if (is.null(private$draws_$latent)) nlatent <- 0
      else nlatent <- posterior::nvariables(private$draws_$latent)

      model <- spaced_choice(private$name_)
      print_rows_("Posterior draws of model {.val {model}}.", c(
        Iterations = "{niter}",
        Chains = "{nchains}",
        if (nparams > 0) c(Parameters = "{nparams}"),
        if (nlatent > 0) c(Latent = "{nlatent}"),
        `Inference time` = "{format_secs_(private$secs_)}"
      ), warnings = private$warnings_)
      invisible()
    }
  )
)

#' @title Combine posterior draws of one model
#' @description
#' Bind the chains of several runs of [run_inference()],
#' so a run grows without sampling its chains again.
#' @param ... Objects of class [HousepiPosteriorDraws].
#' @return An object of class [HousepiPosteriorDraws]
#' holding the chains of all runs, numbered in sequence.
#' @details
#' The runs must share the model, its definitions, the attached data,
#' and the sampler settings other than the number of chains.
#' They must also store the latent draws of the same members,
#' see `id_draws` in [run_inference()].
#' Convergence diagnostics are computed again over all chains,
#' and reported as [run_inference()] reports them.
#' Runs after the same `set.seed()` call draw the same chains,
#' and do not combine.
#' @examples
#' hh_model <- flu_1918_frederick |>
#'   new_housepi_model("final-size") |>
#'   add_parameter(import = "template")
#' hh_fit <- run_inference(hh_model)
#' hh_fit_2 <- run_inference(hh_model)
#' c(hh_fit, hh_fit_2)
#' @name combine_draws
#' @export
c.HousepiPosteriorDraws <- function(...) {
  runs <- list(...)
  pooled <- c("draws_", "secs_", "seed_", "warnings_", "diagnostics_")
  comb <- combine_r6_(runs, pooled, "runs")
  privs <- comb$privs
  seeds <- unlist(lapply(privs, `[[`, "seed_"))
  if (anyDuplicated(seeds)) stop_h(
    "The runs repeat a seed.",
    "Runs after the same {.fn set.seed} call draw the same chains."
  )
  # chains bound per block of draws
  bind <- function(block) {
    parts <- lapply(privs, \(pvt) pvt$draws_[[block]])
    if (is.null(parts[[1L]])) return(NULL)
    vars <- lapply(parts, posterior::variables)
    if (!all(vapply(vars, identical, logical(1L), vars[[1L]]))) stop_h(
      "The runs store draws of different variables.",
      "Pass the same {.arg id_draws} to {.fn run_inference}."
    )
    tag_as_housepi_draws(
      do.call(posterior::bind_draws, c(parts, along = "chain"))
    )
  }
  pvt <- comb$pvt
  pvt$draws_ <- lapply(stats::setNames(nm = names(pvt$draws_)), bind)
  pvt$seed_ <- seeds
  pvt$diagnose_()
  comb$out
}

########  Getters  ########

# Keep the requested variables, all of them when none is named.
# * name -> itself and its indexed elements `name[i]`;
# * `prefix` -> regex the column may start with, `(.+_)?` selecting
#   the `<id>_<latent>` traces of a latent name.
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

# Assert the object is a `HousepiPosteriorDraws`.
stopifnotdraws <- function(hh_fit) {
  if (!inherits(hh_fit, "HousepiPosteriorDraws")) stop_h(
    "{.arg hh_fit} must be of class {.cls HousepiPosteriorDraws}.",
    "Create one with {.fn new_housepi_model}."
  )
}

#' @title Get posterior draws of respective variables
#' @description Return a posterior draws data frame with the MCMC samples.
#' `get_parameter_draws()` and `get_latent_draws()`
#' warn and return an empty draws object when no such draws were stored.
#' @details
#' Write \eqn{x} for a single household's data,
#' \eqn{\theta} for the parameter vector
#' and \eqn{z} for that household's latent variables.
#' \describe{
#'   \item{`get_parameter_draws()`}{Draws of \eqn{\theta}.}
#'   \item{`get_latent_draws()`}{Draws of \eqn{z}
#'     for the individuals selected by `id_draws` in [run_inference()].}
#'   \item{`get_log_lik_draws()`}{Pointwise per-household log-likelihood
#'     \eqn{\log p^\mathrm{lat}(z \mid \theta)
#'       + \log p^\mathrm{transm}(x \mid z, \theta)
#'       + \log p^\mathrm{obs}(x \mid z, \theta)},
#'     the latent, transmission and observation model contributions,
#'     one column per household.}
#'   \item{`get_complete_draws()`}{Parameters, the pointwise log-likelihood
#'     renamed `log_lik[1]`, `log_lik[2]`, ... in household order,
#'     and the single column `lprior`,
#'     the parameter log-prior \eqn{\log p(\theta)}.
#'     The latent-variable log-density enters `log_lik`, not `lprior`.
#'     Both are evaluated along the joint (latent, parameter) chain.
#'     `lprior` alone is read with `get_complete_draws(hh_fit, "lprior")`.}
#' }
#' See [HousepiPosteriorDraws] for how the joint log-density
#' splits into `lprior` and `log_lik`.
#' @param hh_fit `HousepiPosteriorDraws` Object returned by [run_inference()].
#' @param variable `character()` Variables to return, all of them when `NULL`.
#' A name also returns its indexed elements,
#' e.g. `log_lik[1]`, `log_lik[2]`, ... for `log_lik` in `get_complete_draws()`.
#' `get_log_lik_draws()` names its columns by `hh_id`,
#' so a household id returns that household's draws.
#' Latent draws are named `<id>_<latent>`,
#' so a latent name returns its draws for every stored individual,
#' and a full name only that individual's.
#' For finer selections use [posterior::subset_draws()].
#' @return `posterior::draws_df`
#' For `log_lik`/`lprior` columns: `r rd_natural_log_`
#' @examples
#' hh_fit <- flu_1918_frederick |>
#'   new_housepi_model("final-size") |>
#'   add_parameter(import = "template") |>
#'   run_inference()
#' get_parameter_draws(hh_fit, c("alpha", "beta"))
#'
#' # parameters, pointwise log_lik and lprior side by side
#' complete <- get_complete_draws(hh_fit)
#' posterior::variables(complete) |> head()
#' @rdname get_draws
#' @export
get_parameter_draws <- function(hh_fit, variable = NULL) {
  stopifnotdraws(hh_fit)
  out <- draws_block_(hh_fit, "param")
  if (posterior::nvariables(out) == 0) warn_h(
    "Empty draws: no parameter was defined in the model."
  )
  subset_draws_(out, variable)
}

#' @rdname get_draws
#' @export
get_latent_draws <- function(hh_fit, variable = NULL) {
  stopifnotdraws(hh_fit)
  out <- draws_block_(hh_fit, "latent")
  if (posterior::nvariables(out) == 0) warn_h(
    "Empty draws: latent variables were not defined or not stored.",
    "Set {.arg id_draws} in {.fn run_inference} to store latent draws."
  )
  subset_draws_(out, variable, prefix = "(.+_)?")
}

#' @rdname get_draws
#' @export
get_log_lik_draws <- function(hh_fit, variable = NULL) {
  stopifnotdraws(hh_fit)
  subset_draws_(draws_block_(hh_fit, "point_log_lik"), variable)
}

#' @rdname get_draws
#' @export
get_complete_draws <- function(hh_fit, variable = NULL) {
  stopifnotdraws(hh_fit)
  draws <- hh_fit$.__enclos_env__$private$draws_
  # rename the pointwise log-likelihood after the standard convention
  aux <- draws$point_log_lik
  n_vars <- posterior::nvariables(aux)
  names(aux)[1:n_vars] <- paste0("log_lik[", 1:n_vars, "]")
  # bind parameters, pointwise log_lik and the parameter-only lprior
  # bind_draws strips the tag, re-stamped
  posterior::bind_draws(draws$param, aux, draws$lprior) |>
    tag_as_housepi_draws() |>
    subset_draws_(variable)
}
