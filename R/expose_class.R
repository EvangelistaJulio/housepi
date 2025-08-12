#      __
#     / /_  ____  __  __ ____ __  ____  __
#    / __ \/ __ \/ / / / ___/ _ \/ __ \/ /
#   / / / / /_/ / /_/ (__  )  __/ /_/ / /
#  /_/ /_/\____/\____/____/\___/ ____/_/
#                             /_/
#
# Copyright (C) 2025 CAMAGLIA Francesco - Institut Pasteur, Paris (France).
#
# expose_class.R
#

#' @title New housepi model
#' @description
#' Create a new instance of class 'HousepiModel' 
#' used to run inference and simulation routines.
#'
#' @param choice `r rd_mdl_choice_`
#' @param as_template `r rd_mdl_astmpl_`
#' @param force `r rd_mdl_force_`
#' @param quiet `r rd_mdl_quiet_`
#' 
#' @return `r rd_model_in_`
#' 
#' @export
new_housepi_model <- function(choice = NA, as_template = FALSE, 
                              force = FALSE, quiet = TRUE) {
  HousepiModel$new(choice = choice, as_template = as_template,
                   force = force, quiet = quiet)
}

############################################################################
############################################################################
############################################################################

#' @title Add dataset
#' @description
#' Make an internal copy of the provided household dataset.
#' 
#' @param object `r rd_model_in_`
#' @param dataset `r rd_add_dataset_`
#' 
#' @return `r rd_model_out_`
#' 
#' @export
add_dataset  <- function(object, dataset) {
  stopifnot(inherits(object, "HousepiModel"))
  object$clone(deep = TRUE)$add_dataset(dataset)
}

#' @title Add datalist
#' @description
#' Maje an internal copy of the provided datalist.
#' 
#' @param object `r rd_model_in_`
#' @param datalist `r rd_add_datalist_`
#' 
#' @return `r rd_model_out_`
#' 
#' @export
add_datalist <- function(object, datalist) {
  stopifnot(inherits(object, "HousepiModel"))
  object$clone(deep = TRUE)$add_datalist(datalist)
}

#' Re-export dplyr::summarise
#'
#' @importFrom dplyr summarise
#' @export
dplyr::summarise
#' @title Summarise loaded household dataset 
#' @description Summarise statuses in the dataset.
#' 
#' @param object `r rd_model_in_`
#' @param ... Ignored.
#' 
#' @return `data.frame`
#' 
#' @exportS3Method summarise HousepiModel
summarise.HousepiModel <- function(object, ...) {
  object$summarise()
}

############################################################################
############################################################################
############################################################################

#' @title Drop parameters from housepi model
#' @description
#' Drop parameter definitions.
#' 
#' @param object `r rd_model_in_`
#' @param ... Ignored.
#'
#' @return `r rd_model_out_`
#' 
#' @export
drop_parameters <-  function(object, ...) {
  stopifnot(inherits(object, "HousepiModel"))
  object$clone(deep = TRUE)$drop_parameters(...)
}

#' @title Drop augmented variables from housepi model
#' @description
#' Drop augmented variable definitions.
#' 
#' @param object `r rd_model_in_`
#' @param ... Ignored.
#'
#' @return `r rd_model_out_`
#' 
#' @export
drop_augmented <-  function(object, ...) {
  stopifnot(inherits(object, "HousepiModel"))
  object$clone(deep = TRUE)$drop_augmented(...)
}

#' @title Drop group rates from housepi model
#' @description
#' Drop group rates definitions.
#' 
#' @param object `r rd_model_in_`
#' @param ... Ignored.
#'
#' @return `r rd_model_out_`
#' 
#' @export
drop_group_rates <-  function(object, ...) {
  stopifnot(inherits(object, "HousepiModel"))
  object$clone(deep = TRUE)$drop_group_rates(...)
}

#' @title Drop all definitions from housepi model
#' @description
#' Drop all definitions.
#' 
#' @param object `r rd_model_in_`
#'
#' @return `r rd_model_out_`
#' 
#' @export
drop_all <- function(object) {
  stopifnot(inherits(object, "HousepiModel"))
  object$clone(deep = TRUE)$drop_all()
}

############################################################################
############################################################################
############################################################################

#' @title Get parameter definitions
#' @description 
#' Get parameter formatted definitions.
#' 
#' @param object `r rd_model_in_`
#' @param ... Ignored.
#' 
#' @return `data.frame`
#' 
#' @export
get_parameters <- function(object, ...) {
  stopifnot(inherits(object, "HousepiModel"))
  object$parameters
}

#' @title Add parameter definition
#' @description
#' Add a new parameter to houspepi model definitions.
#' If the parameter already exists, it is overwritten.
#' 
#' @param object `r rd_model_in_`
#' @param name `r rd_var_name_`
#' @param value `r rd_var_value_`
#' @param const `r rd_var_const_`
#' @param discr `r rd_par_discr_`
#' @param prior `r rd_par_model_`
#' @param args `r rd_par_args_`
#' @param min `r rd_var_min_`
#' @param max `r rd_var_max_`
#' @param rule `r rd_par_rule_`
#' @param scale `r rd_var_scale_`
#' 
#' @return `r rd_model_out_`
#' 
#' @export
add_parameter <- function(object,
                          name = NA, value = NA,
                          const = NA, discr = FALSE,
                          prior = NA, args = list(),
                          min = NA, max = NA,
                          rule = NA, scale = 1) {
  stopifnot(inherits(object, "HousepiModel"))
  object$clone(deep = TRUE)$add_parameter(
    name = name, value = value,
    const = const, discr = discr, 
    prior = prior, args = args,
    min = min, max = max,
    rule = rule, scale = scale
  )
}

#' @title Load list of parameter definitions
#' @description
#' Load parameter definitions either as arguments or from predefined list.
#'
#' @param object `r rd_model_in_`
#' @param ... `r rd_list_args_`
#' @param from `r rd_list_from_`
#' @param append `r rd_list_append_`
#' 
#' @return `r rd_model_out_`
#' 
#' @export
load_parameter_list <- function(object, ..., from = NA, append = FALSE) {
  stopifnot(inherits(object, "HousepiModel"))
  object$clone(deep = TRUE)$load_parameter_list(
    ..., from = from, append = append
  )
}

#' @title Load parameter template definitions
#' @description
#' Load parameter definitions from model template.
#' 
#' @param object `r rd_model_in_`
#' @param append `r rd_tmpl_append_`
#' 
#' @return `r rd_model_out_`
#' 
#' @export
load_parameter_template <- function(object, append = FALSE) {
  stopifnot(inherits(object, "HousepiModel"))
  object$clone(deep = TRUE)$load_parameter_template(append = append)
}

############################################################################
############################################################################
############################################################################

#' @title Get augmented definitions
#' @description 
#' Get augmented variable/function formatted definitions.
#' 
#' @param object `r rd_model_in_`
#' @param ... Ignored.
#' 
#' @return `data.frame`
#' 
#' @export
get_augmented <- function(object, ...) {
  stopifnot(inherits(object, "HousepiModel"))
  object$augmented
}

#' @title Add augmented definition
#' @description
#' Add a new augmented variable/function to houspepi model definitions.
#' If the augmented variablw already exists, it is updated.
#' 
#' @param object `r rd_model_in_`
#' @param name `r rd_var_name_`
#' @param status `r rd_aug_status_`
#' @param value `r rd_aug_val_`
#' @param init_rand `r rd_aug_rndin_`
#' @param const `r rd_var_const_`
#' @param discr `r rd_aug_discr_`
#' @param func `r rd_aug_func_`
#' @param func_args `r rd_aug_funargs_`
#' @param model `r rd_aug_model_`
#' @param args `r rd_aug_args_`
#' @param min `r rd_aug_min_`
#' @param max `r rd_aug_max_`
#' @param rule `r rd_var_rule_`
#' @param scale `r rd_var_scale_`
#' @param from `r rd_aug_from_`
#' @param expr `r rd_aug_expr_`
#' 
#' @return `r rd_model_out_` 
#' 
#' @export
add_augmented <- function(object, 
                          name = NA, status = NA,
                          const = NA, discr = FALSE,
                          func = NA, func_args = NULL,
                          value = NA, init_rand = FALSE,
                          model = NA, args = list(),
                          min = NA, max = NA,
                          rule = NA, scale = 1, from = NA,
                          expr = NA) {
  stopifnot(inherits(object, "HousepiModel"))
  object$clone(deep = TRUE)$add_augmented(
    name = name, status = status,
    const = const, discr = discr,
    func = func, func_args = func_args,
    value = value, init_rand = init_rand,
    model = model, args = args,
    min = min, max = max,
    rule = rule, scale = scale, from = from,
    expr = expr
  )
}

#' @title Load list of augmented definitions
#' @description
#' Load augmented variable/function definitions 
#' either as arguments or from predefined list.
#'
#' @param object `r rd_model_in_`
#' @param ... `r rd_list_args_` 
#' @param from `r rd_list_from_`
#' @param append `r rd_list_append_`
#' 
#' @return `r rd_model_out_` 
#' 
#' @export
load_augmented_list <- function(object, ..., from = NA, append = FALSE) {
  stopifnot(inherits(object, "HousepiModel"))
  object$clone(deep = TRUE)$load_augmented_list(
    ..., from = from, append = append
  )
}

#' @title Load augmented template definitions
#' @description
#' Load augmented variable/function definitions from template.
#'
#' @param object `r rd_model_in_`
#' @param append `r rd_tmpl_append_`
#' 
#' @return `r rd_model_out_` 
#' 
#' @export
load_augmented_template <- function(object, append = FALSE) {
  stopifnot(inherits(object, "HousepiModel"))
  object$clone(deep = TRUE)$load_augmented_template(append = append)
}

############################################################################
############################################################################
############################################################################

#' @title Get group rates definitions
#' @description 
#' Get group rates formatted definitions.
#' 
#' @param object `r rd_model_in_`
#' @param ... Ignored.
#' 
#' @return `data.frame`
#' 
#' @export
get_group_rates <- function(object, ...) {
  stopifnot(inherits(object, "HousepiModel"))
  # internal build requires deep clone to avoid inplace
  object$clone(deep = TRUE)$group_rates
}

#' @title Add group rate definition
#' @description
#' Add definitions of a group rate by given dataset feature(s).
#' The reference value is equal to 1 and can be specified.
#' 
#' @param object `r rd_model_in_`
#' @param name `r rd_var_name_`
#' @param by `r rd_grr_by_`
#' @param prefix `r rd_grr_prefix_`
#' @param reference `r rd_grr_refer_`
#' @param independent `r rd_grr_indep_`
#' @param prior `r rd_par_model_`
#' Default is "lnorm". 
#' @param args `r rd_par_args_`
#' Default is `list(meanlog = 0, sdlog = 3)`.
#' @param min `r rd_var_min_`
#' Default is 0.
#' @param max `r rd_var_max_`
#' Default is 10.
#' @param rule `r rd_var_rule_`
#' @param scale `r rd_var_scale_`
#'
#' @references
#' Cauchemez, et al. (2009)
#' Household Transmission of 2009 Pandemic Influenza A (H1N1) Virus in the United States.
#' New England Journal of Medicine. 361: 2619-2627. doi:10.1056/NEJMoa0905498
#' 
#' @return `r rd_model_out_` 
#' 
#' @export 
add_group_rate <-  function(object,
                            name, by, prefix = NA, 
                            reference = NULL, independent = TRUE,
                            prior = "lnorm", 
                            args = list(meanlog = 0, sdlog = 3),
                            min = 0, max = 10,
                            rule = "lnorm", scale = 1) {
  stopifnot(inherits(object, "HousepiModel"))
  object$clone(deep = TRUE)$add_group_rate(
    name = name, by = by, prefix = prefix,
    reference = reference, independent = independent,
    prior = prior, args = args,
    min = min, max = max,
    rule = rule, scale = scale
  )
}

#' @title Load list of rroup rate definitions
#' @description
#' Load group rate definitions either as direct arguments
#' or from predefined list.
#' 
#' @param object `r rd_model_in_`
#' @param ... ` r rd_list_args_`
#' @param from `r rd_list_from_`
#' @param append `r rd_list_append_`
#' 
#' @return `r rd_model_out_` 
#' 
#' @export
load_group_rate_list <- function(object, ..., from = NA, append = FALSE) {
  stopifnot(inherits(object, "HousepiModel"))
  object$clone(deep = TRUE)$load_group_rate_list(
    ..., from = from, append = append
  )
}

############################################################################
############################################################################
############################################################################

#' @title Add infectivity profile definition
#' @description
#' Add definition of infectivity profile.
#'
#' @param object `r rd_model_in_`
#' @param name `r rd_var_name_`
#' @param discr `r rd_inf_discr_`
#' @param args `r rd_inf_args_`
#' @param max `r rd_aug_max_`
#' 
#' @return `r rd_model_out_` 
#' 
#' @export
add_infec_profile <- function(object, 
                              name, discr = FALSE,
                              args = list(), max = NA) {
  stopifnot(inherits(object, "HousepiModel"))
  object$clone(deep=TRUE)$add_infec_profile(
    name = name, discr = discr,
    args = args, max = max
  )
}

############################################################################
############################################################################
############################################################################

#' @title Compute loglikelihood of housepi model
#' @description
#' Compute the log-likelihood of a `HousepiModel` object.
#'
#' @param object `r rd_model_in_`
#' @param dataset `r rd_smpl_dataset_`
#' @param seed `r rd_seed_`
#'
#' @return `numeric(1)`
#' 
#' @export
loglikelihood <- function(object, dataset = NULL, seed = NA) {
  stopifnot(inherits(object, "HousepiModel"))
  object$loglikelihood(dataset = dataset, seed = seed)
}

#' @title Run inference on housepi model
#' @description
#' Run inference with Bayesian Markov Chain Monte Carlo.
#' 
#' @param object `r rd_model_in_`
#' @param dataset `r rd_smpl_dataset_`
#' @param sample_size `r rd_smpl_size_`
#' @param n_steps `r rd_smpl_steps_`
#' @param n_burnin `r rd_smpl_burn_`
#' @param thin `r rd_smpl_thin_`
#' @param n_adapt `r rd_smpl_adapt_` 
#' @param n_chains `r rd_smpl_chain_`
#' @param n_threads `r rd_smpl_thread_`
#' @param rand_init `r rd_smpl_rndin_`
#' @param id_draws `r rd_smpl_iddrw_`
#' @param seed `r rd_seed_`
#' @param progress `r rd_progress_`
#' 
#' @return `HousepiPosteriorDraws`\cr
#' Container with posterior draws.
#' 
#' @export 
infer <- function(object, dataset = NULL,
                  sample_size = NA, n_steps = 5000,
                  n_burnin = 0, thin = 1, n_adapt = 0,
                  n_chains = 1, n_threads = 1,
                  rand_init = TRUE, id_draws = NA,
                  seed = NA, progress = TRUE) {
  stopifnot(inherits(object, "HousepiModel"))
  object$infer(
    dataset = dataset, sample_size = sample_size,
    n_steps = n_steps, n_burnin = n_burnin, thin = thin,
    n_adapt = n_adapt, n_chains = n_chains,
    n_threads = n_threads, rand_init = rand_init,
    id_draws = id_draws, seed = seed, progress = progress
  )
}

#' @title Estimate evidence of housepi model
#' @description
#' Estimate the model evidence (or marginal likelihood).
#' It uses the method of power posteriors with different available schemes.
#'
#' @param dataset `r rd_smpl_dataset_`
#' @param scheme `r rd_evi_scheme_`
#' @param n_rungs `r rd_evi_rungs_`
#' @param sample_size `r rd_smpl_size_`
#' @param n_steps `r rd_smpl_steps_`
#' @param n_burnin `r rd_smpl_burn_`
#' @param thin `r rd_smpl_thin_`
#' @param n_adapt `r rd_smpl_adapt_`
#' @param n_chains `r rd_smpl_chain_`
#' @param n_threads `r rd_smpl_thread_`
#' @param seed `r rd_seed_`
#' @param progress `r rd_progress_`
#'
#' @references
#' Friel N, Hurn M, Wyse J.
#' Improving power posterior estimation of statistical evidence.
#' Stat Comput. 2014;24: 709-723. doi:10.1007/s11222-013-9397-1
#' 
#' @return `numeric(1)`
#' 
#' @export
evidence <- function(object, dataset = NULL,
                     scheme = "power-fractions", n_rungs = 40,
                     sample_size = NA, n_steps = 5000,
                     n_burnin = 0, thin = 1, n_adapt = 0,
                     n_chains = 1, n_threads = 1,
                     seed = NA, progress = TRUE) {
  stopifnot(inherits(object, "HousepiModel"))
  object$evidence(dataset = dataset, scheme = scheme, n_rungs = n_rungs,
                  sample_size = sample_size, n_steps = n_steps,
                  n_burnin = n_burnin, thin = thin, n_adapt = n_adapt,
                  n_chains = n_chains, n_threads = n_threads,
                  seed = seed, progress = progress)
}

############################################################################
############################################################################
############################################################################

#' Re-export stats::simulate
#'
#' @importFrom stats simulate
#' @export
stats::simulate
#' @title Simulate household dataset
#' @description
#' Generate a simulated dataset.
#'
#' @param object `r rd_model_in_`
#' @param dataset `r rd_smpl_dataset_`
#' @param seed `r rd_seed_`
#' @param progress `r rd_progress_`
#' @param ... Ignored.
#'
#' @return `data.frame`
#' 
#' @exportS3Method simulate HousepiModel
simulate.HousepiModel <- function(object, dataset = NULL, 
                                  seed = NA, progress = FALSE,
                                  ...) {
  stopifnot(inherits(object, "HousepiModel"))
  object$simulate(
    dataset = dataset, seed = seed, progress = progress
  )
}
  