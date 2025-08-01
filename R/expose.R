#      __
#     / /_  ____  __  __ ____ __  ____  __
#    / __ \/ __ \/ / / / ___/ _ \/ __ \/ /
#   / / / / /_/ / /_/ (__  )  __/ /_/ / /
#  /_/ /_/\____/\____/____/\___/ ____/_/
#                             /_/
#
# Copyright (C) 2025 CAMAGLIA Francesco - Institut Pasteur, Paris (France).
#
# expose.R
#

#' @title New Housepi Model
#' @description
#' Create a new instance of class 'HousepiModel' 
#' used to run inference and simulation routines.
#'
#' @param choice `character(1)`\cr
#' Choose a model providing one between:
#' - template name;
#' - path/to/folder with implementation;
#' - path/to/folder with previously compiled model;
#' - "blank": template name for blank model (default).
#' @param as_template `logical(1)`\cr
#' If `TRUE` store the dynamic library as a local template.
#' IF `FALSE` recompile the model at every new session.
#' @param force `logical(1)`\cr
#' If `TRUE` compile even if a cached model already exists.
#' @param name `chracter(1)`\cr
#' Custom name of the model (default is inferred from `choice`).
#' @param cri_level `double(1)`\cr
#' Level for the credible interval.
#' Value between 0 and 1 (default is 0.95).
#' @param quiet `logical(1)`\cr
#' If `TRUE` hide `pkgbuild` compilation commands (default).
#' 
#' @return `HousepiModel` object.
#' @export
new_housepi_model <- function(choice = NA, as_template = FALSE, force = FALSE,
                              name = NA, cri_level = 0.95, quiet = TRUE) {
  HousepiModel$new(choice = choice, as_template = as_template,
                   force = force, name = name, cri_level = cri_level,
                   quiet = quiet)
}

############################################################################
############################################################################
############################################################################

#' @title Add dataset
#' @description
#' Make an internal copy of the provided household dataset.
#' 
#' @param object `HousepiModel`\cr
#' A HousepiModel object.
#' @param ... Ignored.
#' #' @param dataset `data.frame`\cr
#' Input dataset for analysis. 
#' 
#' @return `HousepiModel` object with updated dataset.
#' 
#' @export
add_dataset <- function(object, ...) {
  UseMethod("add_dataset")
}
#' @exportS3Method add_dataset HousepiModel
add_dataset.HousepiModel <- function(object, dataset, ...) {
  object$clone(deep = TRUE)$add_dataset(dataset)
}

#' @title Add datalist
#' @description
#' Maje an internal copy of the provided datalist.
#' 
#' @param object `HousepiModel`\cr
#' A HousepiModel object.
#' @param ... Ignored.
#' #' @param datalist `data.frame`\cr
#' Input datalist for analysis. 
#' 
#' @return `HousepiModel` object with updated datalist.
#' 
#' @export
add_datalist <- function(object, ...) {
  UseMethod("add_datalist")
}
#' @exportS3Method add_datalist HousepiModel
add_datalist.HousepiModel <- function(object, datalist, ...) {
  object$clone(deep = TRUE)$add_datalist(datalist)
}

#' Re-export dplyr::summarise
#'
#' @importFrom dplyr summarise
#' @export
dplyr::summarise

#' @title Summarise Housepi Model Statuses
#' @description Summarise statuses in the dataset.
#' 
#' @param object `HousepiModel`\cr
#' A HousepiModel object with a dataset.
#' @param ... Ignored.
#' 
#' @return
#' `data.frame`
#' 
#' @importFrom dplyr summarise
#' @exportS3Method summarise HousepiModel
summarise.HousepiModel <- function(object, ...) {
  object$summarise()
}

############################################################################
############################################################################
############################################################################

#' @title Drop Parameters from Housepi Model
#' @description
#' Drop current parameter definitions.
#' 
#' @param object `HousepiModel`\cr
#' A HousepiModel object.
#' @param ... Ignored.
#'
#' @return `HousepiModel` object with cleared parameters.
#' 
#' @export
drop_parameters <- function(object, ...) {
  UseMethod("drop_parameters")
}
#' @exportS3Method drop_parameters HousepiModel
drop_parameters.HousepiModel <- function(object, ...) {
  object$clone(deep = TRUE)$drop_parameters(...)
}

#' @title Drop Augmented Variables from Housepi Model
#' @description
#' Drop current augmented variable definitions.
#' 
#' @param object `HousepiModel`\cr
#' A HousepiModel object.
#' @param ... Ignored.
#'
#' @return `HousepiModel` object with cleared augmented variables.
#' 
#' @export
drop_augmented <- function(object, ...) {
  UseMethod("drop_augmented")
}
#' @exportS3Method drop_augmented HousepiModel
drop_augmented.HousepiModel <- function(object, ...) {
  object$clone(deep = TRUE)$drop_augmented(...)
}

#' @title Drop Group Rates from Housepi Model
#' @description
#' Drop current group rates definitions.
#' 
#' @param object `HousepiModel`\cr
#' A HousepiModel object.
#' @param ... Ignored.
#'
#' @return `HousepiModel` object with cleared group rates.
#' 
#' @export
drop_group_rates <- function(object, ...) {
  UseMethod("drop_group_rates")
}
#' @exportS3Method drop_group_rates HousepiModel
drop_group_rates.HousepiModel <- function(object, ...) {
  object$clone(deep = TRUE)$drop_group_rates(...)
}

#' @title Drop All Definitions from Housepi Model
#' @description
#' Drop all current definitions.
#' 
#' @param object `HousepiModel`\cr
#' A HousepiModel object.
#' @param ... Ignored.
#'
#' @return `HousepiModel` object without definitions.
#' 
#' @export
drop_all <- function(object, ...) {
  UseMethod("drop_all")
}
#' @exportS3Method drop_all HousepiModel
drop_all.HousepiModel <- function(object, ...) {
  object$clone(deep = TRUE)$drop_all(...)
}

############################################################################
############################################################################
############################################################################

#' @title Get Parameter Definitions
#' @description 
#' Get parameter formatted definitions.
#' 
#' @param object `HousepiModel`\cr
#' A HousepiModel object.
#' @param ... Ignored.
#' 
#' @return `data.frame`
#' 
#' @export
get_parameters <- function(object, ...) {
  UseMethod("get_parameters")
}
#' @exportS3Method get_parameters HousepiModel
get_parameters.HousepiModel <- function(object, ...) {
  object$parameters
}

#' @title Add Parameter Definition
#' @description
#' Add a new parameter to houspepi model definitions.
#' If the parameter already exists, it is updated.
#' 
#' @param object `HousepiModel`\cr
#' A HousepiModel object.
#' @param ... Ignored.
#' @param name `character(1)`\cr
#' Name of the parameter.
#' @param value `double(1)`\cr
#' Initial value of the parameter.
#' @param const `logical(1)`\cr
#' If `FALSE` the parameter is treated as a constant.
#' @param discr `logical(1)`\cr
#' If `TRUE` the parameter is discrete.
#' @param prior `character(1)`\cr
#' Name of the prior distribution.
#' @param args named `list`\cr
#' Arguments for the prior distribution.
#' @param rule `character(1)`\cr
#' Name of the inference proposal rule.
#' @param scale `double(1)`\cr
#' Scale parameter for the proposal rule.
#' Changed automatically during inference if `n_adpat`>0.
#' @param min `double(1)`\cr
#' Value of the minimum.
#' @param max `double(1)`\cr
#' Value of the maximum.
#' 
#' @return `HousepiModel` object with updated definitions.
#' 
#' @export
add_parameter <- function(object, ...) {
  UseMethod("add_parameter")
}
#' @exportS3Method add_parameter HousepiModel
add_parameter.HousepiModel <- function(object,
                                       name = NA, value = NA,
                                       const = NA, discr = FALSE,
                                       prior = NA, args = list(),
                                       rule = NA, scale = NA,
                                       min = NA, max = NA,
                                       ...) {
  object$clone(deep = TRUE)$add_parameter(
    name = name, value = value,
    const = const, discr = discr, 
    prior = prior, args = args,
    rule = rule, scale = scale,
    min = min, max = max
  )
}

#' @title Load List of Parameter Definitions
#' @description
#' Load parameter definitions either as arguments or from predefined list.
#'
#' @param object `HousepiModel`\cr
#' A HousepiModel object.
#' @param ...
#' Named arguments with respective list of definitions.
#' @param from named `list`\cr
#' Named list of parameter definitions.
#' @param append `logical(1)`\cr
#' If `TRUE` appends to current definitions.
#' If `FALSE` drop all pre-existing definitions (default).
#' 
#' @return `HousepiModel` object with updated definitions.
#' 
#' @export
load_parameter_list <- function(object, ..., from = NA, append = FALSE) {
  object$clone(deep = TRUE)$load_parameter_list(
    ..., from = from, append = append
  )
}

#' @title Load Parameter Template Definitions
#' @description
#' Load parameter definitions from model template.
#'
#' @param object `HousepiModel`\cr
#' A HousepiModel object.
#' @param ... Ignored.
#' @param append `logical(1)`\cr
#' If `TRUE` append missing parameter definitions.
#' If `FALSE` drop all pre-existing definitions (default).
#' 
#' @return `HousepiModel` object with updated definitions.
#' 
#' @export
load_parameter_template <- function(object, append = FALSE, ...) {
  object$clone(deep = TRUE)$load_parameter_template(append = append)
}

############################################################################
############################################################################
############################################################################

#' @title Get augmented definitions
#' @description 
#' Get augmented variable/function formatted definitions.
#' 
#' @param object `HousepiModel`\cr
#' A HousepiModel object.
#' @param ... Ignored.
#' 
#' @return `data.frame`
#' 
#' @export
get_augmented <- function(object, ...) {
  UseMethod("get_augmented")
}
#' @exportS3Method get_augmented HousepiModel
get_augmented.HousepiModel <- function(object, ...) {
  object$augmented
}

#' @title Add Augmented Definition
#' @description
#' Add a new augmented variable/function to houspepi model definitions.
#' If the augmented variablw already exists, it is updated.
#' 
#' @param object `HousepiModel`\cr
#' A HousepiModel object.
#' @param ... Ignored.
#' @param name `character(1)`\cr
#' Name of the augmented variable.
#' @param status `integer(1)`\cr
#' The associated status of the definition.
#' Default NA is for all infectious.
#' @param value `double(1)` | `character(1)`\cr
#' One between:
#' - Scalar value;
#' - Parameter name (optionally followed by `@P`);
#' - Data name (followed by `@D`).
#' @param init_rand `character(1)`\cr
#' If `TRUE` draw initial vale from model distribution.
#' @param const `logical(1)`\cr
#' If `FALSE` the augmented variable is treated as a constant.
#' @param discr `logical(1)`\cr
#' IF `TRUE` the augmenter variable is discrete.
#' @param func `character()`\cr
#' Name of the augmented function.
#' @param func_args named `list`\cr
#' Arguments of the augmented function.
#' @param model `character(1)`\cr
#' Name of the model distribution.
#' @param args named `list`\cr
#' Arguments for the model distribution.
#' @param min `double(1)` | `character(1)`\cr
#' Minimum value. One between:
#' - Scalar value;
#' - Parameter name (optionally followed by `@P`);
#' - Data name (followed by `@D`).
#' @param max `double(1)` | `character(1)`\cr
#' Maxmimum value. One between:
#' - Scalar value;
#' - Parameter name (optionally followed by `@P`);
#' - Data name (followed by `@D`).
#' @param rule `character(1)`\cr
#' Name of the inference proposal rule.
#' @param scale `double(1)`\cr
#' Scale parameter for the proposal rule.
#' @param from `double(1)` | `character(1)`\cr
#' Advanced options for non random walk moves.
#' `FIXME: missing`
#' 
#' @return `HousepiModel` object with updated definitions.
#' 
#' @export
add_augmented <- function(object, ...) {
  UseMethod("add_augmented")
}
#' @exportS3Method add_augmented HousepiModel
add_augmented.HousepiModel <- function(object, 
                                       name = NA, status = NA,
                                       const = NA, discr = FALSE,
                                       func = NA, func_args = NULL,
                                       value = NA, init_rand = FALSE,
                                       model = NA, args = list(),
                                       min = NA, max = NA,
                                       rule = NA, scale = NA, from = NA,
                                       ...) {
  object$clone(deep = TRUE)$add_augmented(
    name = name, status = status,
    const = const, discr = discr,
    func = func, func_args = func_args,
    value = value, init_rand = init_rand,
    model = model, args = args,
    min = min, max = max,
    rule = rule, scale = scale, from = from
  )
}

#' @title Load List of Augmented Definitions
#' @description
#' Load augmented variable/function definitions 
#' either as arguments or from predefined list.
#'
#' @param object `HousepiModel`\cr
#' A HousepiModel object.
#' @param ... Ignored.
#' Named arguments with respective list of definition details.
#' @param from named `list`\cr
#' Named list of parameter definitions.
#' @param append `logical(1)`\cr
#' If `TRUE` appends to current definitions.
#' If `FALSE` drop all pre-existing definitions (default).
#' 
#' @return `HousepiModel` object with updated definitions.
#' 
#' @export
load_augmented_list <- function(object, ..., from = NA, append = FALSE) {
  object$clone(deep = TRUE)$load_augmented_list(
    ..., from = from, append = append
  )
}

#' @title Load Augmented Template Definitions
#' @description
#' Load augmented variable/function definitions from template.
#'
#' @param object `HousepiModel`\cr
#' A HousepiModel object.
#' @param ... Ignored.
#' @param append
#' If `TRUE` append missing augmented variable definitions.
#' If `FALSE` drop all pre-existing definitions (default).
#' 
#' @return `HousepiModel` object with updated definitions.
#' 
#' @export
load_augmented_template <- function(object, append = FALSE, ...) {
  object$clone(deep = TRUE)$load_augmented_template(append = append)
}

############################################################################
############################################################################
############################################################################

#' @title Get group rates definitions
#' @description 
#' Get group rates formatted definitions.
#' 
#' @param object `HousepiModel`\cr
#' A HousepiModel object.
#' @param ... Ignored.
#' 
#' @return `data.frame`
#' 
#' @export
get_group_rates <- function(object, ...) {
  UseMethod("get_group_rates")
}
#' @exportS3Method get_group_rates HousepiModel
get_group_rates.HousepiModel <- function(object, ...) {
  # internal build requires deep clone to avoid inplace
  object$clone(deep = TRUE)$group_rates
}

#' @title Add Group Rate Definition
#' @description
#' Add definitions of a group rate by given dataset feature(s).
#' Reference value has unitary rate and can be specified
#' (default is most represented group).
#' 
#' @param object `HousepiModel`\cr
#' A HousepiModel object.
#' @param ... Ignored.
#' @param name `character(1)`\cr
#' Name of the group rate.
#' @param by `list`\cr
#' Feature(s) associated to group rate as dataset headers.
#' Corresponding dataset column have to be categorical.
#' @param prefix `character(1)`\cr
#' Prefix to parameter names associated to the group rate.
#' @param reference named `list`\cr
#' Choose reference group providing corresponding features.
#' Default is the group with the largest sample size.
#' @param independent `logical(1)`\cr
#' If `TRUE` then features are considered independently.
#' @param prior `character(1)`\cr
#' Name of the prior distribution.
#' Default associated parameter definitions
#' have prior "lnorm" and `args = list(meanlog = 0, sdlog = 3)`.
#' @param args named `list`\cr
#' Arguments for the prior distribution.
#' @param min `double(1)`\cr
#' Value of the minimum.
#' @param max `double(1)`\cr
#' Value of the maximum.
#' @param rule `character(1)`\cr
#' Name of the inference proposal rule.
#' @param scale `double(1)`\cr
#' Scale parameter for the proposal rule.
#' Changed automatically during inference if `n_adpat` > 0.
#'
#' @references
#' Cauchemez, et al. (2009)
#' Household Transmission of 2009 Pandemic Influenza A (H1N1) Virus in the United States.
#' New England Journal of Medicine. 361: 2619-2627. doi:10.1056/NEJMoa0905498
#' 
#' @export 
add_group_rate <- function(object, ...) {
  UseMethod("add_group_rate")
}
#' @exportS3Method add_group_rate HousepiModel
add_group_rate.HousepiModel <- function(object,
                                        name, by, prefix = NA, 
                                        reference = NULL, independent = TRUE,
                                        prior = "lnorm", 
                                        args = list(meanlog = 0, sdlog = 3),
                                        min = NA, max = 10,
                                        rule = "lnorm", scale = 1,
                                        ...) {
  object$clone(deep = TRUE)$add_group_rate(
    name = name, by = by, prefix = prefix,
    reference = reference, independent = independent,
    prior = prior, args = args,
    min = min, max = max,
    rule = rule, scale = scale
  )
}

#' @title Load List of Group Rate Definitions
#' @description
#' Load group rate definitions either as direct arguments
#' or from predefined list.
#' 
#' @param object `HousepiModel`\cr
#' A HousepiModel object.
#' @param ...
#' Named arguments with respective list of definition details.
#' @param from named `list`\cr
#' Named list of parameter definitions.
#' @param append `logical(1)`\cr
#' If `TRUE` appends to current definitions.
#' If `FALSE` drop all pre-existing definitions (default).
#' 
#' @export
load_group_rate_list <- function(object, ..., from = NA, append = FALSE) {
  object$clone(deep = TRUE)$load_group_rate_list(
    ..., from = from, append = append
  )
}

############################################################################
############################################################################
############################################################################

#' @title Add Infectivity Profile Definition
#' @description
#' Add definition of infectivity profile.
#'
#' @param object `HousepiModel`\cr
#' A HousepiModel object.
#' @param ... Ignored.
#' @param name `character(1)`\cr
#' Name of the infectivity profile.
#' @param discr `logical(1)`\cr
#' If `TRUE` the profile is discrete.
#' @param args named `list`\cr
#' Arguments (eventually with types) for the profile.
#' @param max `double(1)` | `character(1)`\cr
#' One between:
#' - Scalar value;
#' - Parameter name (optionally followed by `@P`).
#' 
#' @export
add_infec_profile <- function(object, 
                              name, discr = FALSE,
                              args = list(), max = NA,
                              ...) {
  object$clone(deep=TRUE)$add_infec_profile(
    name = name, discr = discr,
    args = args, max = max
  )
}

############################################################################
############################################################################
############################################################################

#' @title Compute Loglikelihood of Housepi Model
#' 
#' @description
#' Compute the log-likelihood of a `HousepiModel` object.
#'
#' @param object `HousepiModel`\cr
#' A HousepiModel object.
#' @param dataset `data.frame`\cr
#' Input dataset for analysis. 
#' If `NULL`, uses the dataset stored in the model. 
#' @param seed `integer(1)` \cr
#' Seed for the random number generator.
#' If `NA` (default), uses R seed.
#' @param ... Ignored.
#'
#' @return `numeric(1)`
#' 
#' @export
loglikelihood <- function(object, ...) {
  UseMethod("loglikelihood")
}
#' @exportS3Method loglikelihood HousepiModel
loglikelihood.HousepiModel <- function(object, dataset = NULL, seed = NA, ...) {
  object$loglikelihood(dataset = dataset, seed = seed)
}

#' @title Run Inference on Housepi Model
#' 
#' @description
#' Run inference with Bayesian Markov Chain Monte Carlo
#' and store internally the chain traces.
#' 
#' @param object `HousepiModel`\cr
#' A HousepiModel object.
#' @param dataset `data.frame` \cr
#' Input dataset for analysis.
#' If `NULL`, uses the dataset stored in the model. 
#' @param sample_size `integer(1)` \cr
#' Number of samples to draw from the posterior distribution in each chain.
#' @param n_steps `integer(1)` \cr
#' Total number of MCMC steps to perform.
#' @param n_burnin `integer(1)` \cr
#' Number of initial steps to discard as burn-in period.
#' @param thin `integer(1)` \cr
#' Thinning interval: number of steps to skip between each saved sample.
#' @param n_adapt `integer(1)` \cr
#' Number of adaptation steps before starting the main sampling.
#' @param n_chains `integer(1)` \cr
#' Number of independent MCMC chains to run.
#' @param n_threads `integer(1)` \cr
#' Number of parallel threads to use for computations.
#' @param rand_init `logical(1)` \cr
#' Controls random initialization of parameters:
#' - If `TRUE` (default), initial values are ndomly drawn from the prior distribution;
#' - If `FALSE` initial value definitions are used.
#' @param id_draws `character(s)` \cr
#' Individual identifiers for augmented trace plots.
#' @param seed `integer(1)` \cr
#' Seed for the random number generator.
#' If `NA` (default), uses R seed.
#' @param progress `logical(1)` | `numeric(1)`\cr
#' Controls progress visualization:
#' - If `TRUE` (default), prints progress every second;
#' - If an integer `n`, prints progress every `n` seconds;
#' - If `FALSE` or `0`, suppresses progress display.
#' @param ... Ignored.
#' 
#' @return `HousepiModel` object with updated draws.
#' 
#' @export 
infer <- function(object, ...) {
  UseMethod("infer")
}
#' @exportS3Method infer HousepiModel
infer.HousepiModel <- function(object, dataset = NULL,
                               sample_size = NA, n_steps = 10000,
                               n_burnin = 0, thin = 1, n_adapt = 0,
                               n_chains = 1, n_threads = 1,
                               rand_init = TRUE, id_draws = NA,
                               seed = NA, progress = TRUE,
                               ...) {
  object$clone(deep = TRUE)$infer(
    dataset = dataset, sample_size = sample_size,
    n_steps = n_steps, n_burnin = n_burnin, thin = thin,
    n_adapt = n_adapt, n_chains = n_chains,
    n_threads = n_threads, rand_init = rand_init,
    id_draws = id_draws, seed = seed, progress = progress
  )
}

#' @title Estimate Evidence of Housepi Model
#' @description
#' Estimate the model evidence (or marginal likelihood).
#' It uses the method of power posteriors with different available schemes.
#'
#' @param dataset `data.frame` \cr
#' Input dataset for simulation.
#' @param scheme `character(1)` \cr
#' Scheme used for evidence estimation.
#' - If "power-fractions" (default) uses inverse temperatures (r / n_rungs)^5.
#' - If "adaptive" uses adaptive thermodynamic ladder.
#' @param n_rungs `integer(1)` \cr
#' Number of rungs in the thermodynamic ladder.
#' @param sample_size `integer(1)` \cr
#' Number of samples to draw from the posterior distribution in each chain.
#' @param n_steps `integer(1)` \cr
#' Total number of MCMC steps to perform.
#' @param n_burnin `integer(1)` \cr
#' Number of initial steps to discard as burn-in period.
#' @param thin `integer(1)` \cr
#' Thinning interval: number of steps to skip between each saved sample.
#' @param n_adapt `integer(1)` \cr
#' Number of adaptation steps before starting the main sampling.
#' @param n_chains `integer(1)` \cr
#' Number of independent MCMC chains to run.
#' @param n_threads `integer(1)` \cr
#' Number of parallel threads to use for computations.
#' @param seed `integer(1)` \cr
#' Seed for the random number generator.
#' - If a specific integer is provided, sets the random number generator seed.
#' - If `NA` (default), uses the system clock to generate a seed.
#' @param progress `logical(1)` | `numeric(1)`\cr
#' Controls progress visualization:
#' - If `TRUE` (default), prints progress every second;
#' - If an integer `n`, prints progress every `n` seconds;
#' - If `FALSE` or `0`, suppresses progress display.
#'
#' @return `double(1)`
#'
#' @references
#' Friel N, Hurn M, Wyse J.
#' Improving power posterior estimation of statistical evidence.
#' Stat Comput. 2014;24: 709-723. doi:10.1007/s11222-013-9397-1
#' 
#' @export
evidence <- function(object, ...) {
  UseMethod("evidence")
}
#' @exportS3Method evidence HousepiModel
evidence.HousepiModel <- function(object, dataset = NULL,
                                  scheme = "power-fractions", n_rungs = 40,
                                  sample_size = NA, n_steps = 10000,
                                  n_burnin = 0, thin = 1, n_adapt = 0,
                                  n_chains = 1, n_threads = 1,
                                  seed = NA, progress = TRUE,
                                  ...) {
  object$evidence(dataset = dataset, scheme = scheme, n_rungs = n_rungs,
                  sample_size = sample_size, n_steps = n_steps,
                  n_burnin = n_burnin, thin = thin, n_adapt = n_adapt,
                  n_chains = n_chains, n_threads = n_threads,
                  seed = seed, progress = progress)
}

############################################################################
############################################################################
############################################################################

#' Re-export posterior::as_draws_df
#'
#' @importFrom posterior as_draws_df
#' @export
posterior::as_draws_df

#' Re-export posterior::as_draws_matrix
#'
#' @importFrom posterior as_draws_matrix
#' @export
posterior::as_draws_matrix

#' @title As Posterior Draws DataFrame
#' 
#' @description
#' Posterior draws dataframe of last inference.
#' 
#' @param object `HousepiModel`\cr
#' A HousepiModel object.
#' @param ... Ignored.
#' 
#' @return `posterior::draws_df`
#' 
#' @importFrom posterior as_draws_df
#' @exportS3Method as_draws_df HousepiModel
as_draws_df.HousepiModel <- function(object, ...) {
  object$as_draws_df()
}

#' @title As Posterior Draws Matrix
#' 
#' @description
#' Posterior draws matrix of last inference.
#' 
#' @param object `HousepiModel`\cr
#' A HousepiModel object.
#' @param ... Ignored.
#' 
#' @return `posterior::draws_matrix`
#' 
#' @importFrom posterior as_draws_matrix
#' @exportS3Method as_draws_matrix HousepiModel
as_draws_matrix.HousepiModel <- function(object, ...) {
  object$as_draws_matrix()
}

#' @title DataFrame of Household Loglikelihood Draws
#' 
#' @description 
#' Posterior draws of pointwise loglikelihood for all households.
#' 
#' @param object `HousepiModel`\cr
#' A HousepiModel object.
#' 
#' @return `posterior::draws_df`
#' 
#' @export
as_house_loglike_draws_df <- function(object) {
  object$as_house_loglike_draws_df()
}

#' @title Matrix of Household Loglikelihood Draws
#' 
#' @description 
#' Posterior draws of pointwise loglikelihood for all households.
#' 
#' @param object `HousepiModel`\cr
#' A HousepiModel object.
#' 
#' @return `posterior::draws_matrix`
#' 
#' @export
as_house_loglike_draws_matrix <- function(object) {
  object$as_house_loglike_draws_matrix()
}

############################################################################
############################################################################
############################################################################

#' Re-export stats::simulate
#'
#' @importFrom stats simulate
#' @export
stats::simulate

#' @title Simulate Household Dataset
#' @description
#' Run simulation and return data.frame with simulated dataset.
#'
#' @param object `HousepiModel`\cr
#' A HousepiModel object.
#' @param dataset `data.frame` \cr
#' Input dataset for simulation.
#' @param seed `integer(1)` \cr
#' Seed for the random number generator.
#' - If a specific integer is provided, sets the random number generator seed.
#' - If `NA` (default), uses the system clock to generate a seed.
#' @param progress `logical(1)` | `numeric(1)`\cr
#' Controls progress visualization:
#' - If `TRUE` (default), prints progress every second;
#' - If an integer `n`, prints progress every `n` seconds;
#' - If `FALSE` or `0`, suppresses progress display.
#' @param ... Ignored.
#'
#' @return `data.frame`
#' 
#' @importFrom stats simulate
#' @exportS3Method simulate HousepiModel
simulate.HousepiModel <- function(object, dataset = NULL, 
                                  seed = NA, progress = FALSE,
                                  ...) {
  object$simulate(
    dataset = dataset, seed = seed, progress = progress
  )
}

