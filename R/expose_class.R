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
#' Create a new object of class 'HousepiModel'
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

#' @title Add data wide
#' @description
#' Make an internal copy of the provided household data.
#'
#' @param object `r rd_model_in_`
#' @param data_wide `r rd_add_data_wide_`
#'
#' @return `r rd_model_out_`
#'
#' @export
add_data_wide  <- function(object, data_wide) {
  stopifnot(inherits(object, "HousepiModel"))
  object$clone(deep = TRUE)$add_data_wide(data_wide)
}

#' @title Add data long
#' @description
#' Male an internal copy of the provided data in long format.
#'
#' @param object `r rd_model_in_`
#' @param data_long `r rd_add_data_long_`
#'
#' @return `r rd_model_out_`
#'
#' @export
add_data_long <- function(object, data_long) {
  stopifnot(inherits(object, "HousepiModel"))
  object$clone(deep = TRUE)$add_data_long(data_long)
}

#' @title Summarise cases in household data
#' @description
#' Summarise cases as given from the statuses in the laoded household data.
#'
#' @param object `r rd_model_in_`
#' @param ... Ignored.
#'
#' @return `data.frame`
#'
#' @export
summarise_cases <- function(object, ...) {
  object$summarise_cases()
}

#' @title Summarise house sizes in household data
#' @description
#' Summarise number of members (size) in the laoded household data.
#'
#' @param object `r rd_model_in_`
#' @param ... Ignored.
#'
#' @return `data.frame`
#'
#' @export
summarise_house_sizes <- function(object, ...) {
  object$summarise_house_sizes()
}

############################################################################
############################################################################
############################################################################

#' @title Drop definitions
#' @description
#' Drop current formatted definitions.
#'
#' @param object `r rd_model_in_`
#' @param ... Ignored.
#'
#' @return `r rd_model_out_`
#'
#' @rdname drop_definitions
#' @export
drop_parameters <-  function(object, ...) {
  stopifnot(inherits(object, "HousepiModel"))
  object$clone(deep = TRUE)$drop_parameters(...)
}

#' @rdname drop_definitions
#' @export
drop_augmented <-  function(object, ...) {
  stopifnot(inherits(object, "HousepiModel"))
  object$clone(deep = TRUE)$drop_augmented(...)
}

#' @rdname drop_definitions
#' @export
drop_group_rates <-  function(object, ...) {
  stopifnot(inherits(object, "HousepiModel"))
  object$clone(deep = TRUE)$drop_group_rates(...)
}

#' @rdname drop_definitions
#' @export
drop_all <- function(object) {
  stopifnot(inherits(object, "HousepiModel"))
  object$clone(deep = TRUE)$drop_all()
}

############################################################################
############################################################################
############################################################################

#' @title Get definitions
#' @description
#' Get current formatted definitions.
#'
#' @param object `r rd_model_in_`
#' @param ... Ignored.
#'
#' @return `data.frame`
#'
#' @rdname get_definitions
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
                          name = NA, value = NULL,
                          const = NA, discr = FALSE,
                          prior = NA, args = NULL,
                          min = NULL, max = NULL,
                          rule = NA, scale = NA) {
  stopifnot(inherits(object, "HousepiModel"))
  object$clone(deep = TRUE)$add_parameter(
    name = name, value = value,
    const = const, discr = discr,
    prior = prior, args = args,
    min = min, max = max,
    rule = rule, scale = scale
  )
}

#' @title Load list of definitions
#' @description
#' Load definitions either as arguments or from predefined list.
#'
#' @param object `r rd_model_in_`
#' @param ... `r rd_list_args_`
#' @param from `r rd_list_from_`
#' @param append `r rd_list_append_`
#'
#' @return `r rd_model_out_`
#'
#' @rdname load_definitions_list
#' @export
load_parameter_list <- function(object, ..., from = NA, append = FALSE) {
  stopifnot(inherits(object, "HousepiModel"))
  object$clone(deep = TRUE)$load_parameter_list(
    ..., from = from, append = append
  )
}

#' @title Load template definitions
#' @description
#' Load default definitions from model template.
#'
#' @param object `r rd_model_in_`
#' @param append `r rd_tmpl_append_`
#'
#' @return `r rd_model_out_`
#'
#' @rdname load_template_definitions
#' @export
load_parameter_template <- function(object, append = FALSE) {
  stopifnot(inherits(object, "HousepiModel"))
  object$clone(deep = TRUE)$load_parameter_template(append = append)
}

############################################################################
############################################################################
############################################################################

#' @rdname get_definitions
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
#' @param expr `r rd_aug_expr_`
#'
#' @return `r rd_model_out_`
#'
#' @export
add_augmented <- function(object,
                          name = NA, status = NA,
                          value = NULL, const = NA, discr = FALSE,
                          func = NA, func_args = NULL,
                          model = NA, args = NULL,
                          min = NULL, max = NULL,
                          rule = NA, scale = NA,
                          expr = NULL) {
  stopifnot(inherits(object, "HousepiModel"))
  object$clone(deep = TRUE)$add_augmented(
    name = name, status = status,
    value = value, const = const, discr = discr,
    func = func, func_args = func_args,
    model = model, args = args,
    min = min, max = max,
    rule = rule, scale = scale,
    expr = expr
  )
}

#' @rdname load_definitions_list
#' @export
load_augmented_list <- function(object, ..., from = NA, append = FALSE) {
  stopifnot(inherits(object, "HousepiModel"))
  object$clone(deep = TRUE)$load_augmented_list(
    ..., from = from, append = append
  )
}

#' @rdname load_template_definitions
#' @export
load_augmented_template <- function(object, append = FALSE) {
  stopifnot(inherits(object, "HousepiModel"))
  object$clone(deep = TRUE)$load_augmented_template(append = append)
}

############################################################################
############################################################################
############################################################################

#' @rdname get_definitions
#' @export
get_group_rates <- function(object, ...) {
  stopifnot(inherits(object, "HousepiModel"))
  # internal build requires deep clone to avoid inplace
  object$clone(deep = TRUE)$group_rates
}

#' @title Add group rate definition
#' @description
#' Add definitions of a group rate by given data feature(s).
#' The reference value is equal to 1 and can be specified.
#'
#' @param object `r rd_model_in_`
#' @param name `r rd_var_name_`
#' @param by `r rd_grr_by_`
#' @param prefix `r rd_grr_prefix_`
#' @param reference `r rd_grr_refer_`
#' @param independent `r rd_grr_indep_`
#' @param prior `r rd_par_model_`
#' Default is "lognormal".
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
                            prior = "lognormal",
                            args = list(meanlog = 0, sdlog = 3),
                            min = 0, max = 10,
                            rule = "lognormal", scale = 1) {
  stopifnot(inherits(object, "HousepiModel"))
  object$clone(deep = TRUE)$add_group_rate(
    name = name, by = by, prefix = prefix,
    reference = reference, independent = independent,
    prior = prior, args = args,
    min = min, max = max,
    rule = rule, scale = scale
  )
}

#' @rdname load_definitions_list
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
                              args = NULL, max = NULL) {
  stopifnot(inherits(object, "HousepiModel"))
  object$clone(deep = TRUE)$add_infec_profile(
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
#' @param seed `r rd_seed_`
#'
#' @return `numeric(1)`
#'
#' @export
loglikelihood <- function(object, seed = NA) {
  stopifnot(inherits(object, "HousepiModel"))
  object$loglikelihood(seed = seed)
}

#' @title Run inference on housepi model
#' @description
#' Run inference with Bayesian Markov Chain Monte Carlo.
#'
#' @param object `r rd_model_in_`
#' @param steps `r rd_smpl_steps_`
#' @param burnin `r rd_smpl_burn_`
#' @param thin `r rd_smpl_thin_`
#' @param adapt `r rd_smpl_adapt_`
#' @param chains `r rd_smpl_chain_`
#' @param threads `r rd_smpl_thread_`
#' @param id_draws `r rd_smpl_iddrw_`
#' @param seed `r rd_seed_`
#' @param progress `r rd_progress_`
#'
#' @return `HousepiPosteriorDraws`\cr
#' Container with posterior draws.
#'
#' @export
infer <- function(object,
                  steps = 5000, burnin = 0,
                  thin = 1, adapt = 0,
                  chains = 1, threads = 1,
                  id_draws = NA, seed = NA, progress = TRUE) {
  stopifnot(inherits(object, "HousepiModel"))
  object$infer(
    steps = steps, burnin = burnin,
    thin = thin, adapt = adapt,
    chains = chains, threads = threads,
    id_draws = id_draws, seed = seed, progress = progress
  )
}

#' @title Estimate logevidence of housepi model
#' @description
#' Estimate the logarithm of the marginal likelihood (aka evidence).
#' It uses the method of power posteriors with different available schemes.
#'
#' @param object `r rd_model_in_`
#' @param scheme `r rd_evi_scheme_`
#' @param rungs `r rd_evi_rungs_`
#' @param steps `r rd_smpl_steps_`
#' @param burnin `r rd_smpl_burn_`
#' @param thin `r rd_smpl_thin_`
#' @param adapt `r rd_smpl_adapt_`
#' @param chains `r rd_smpl_chain_`
#' @param threads `r rd_smpl_thread_`
#' @param seed `r rd_seed_`
#' @param progress `r rd_progress_`
#'
#' @references
#' Friel N, Hurn M, Wyse J.
#' Improving power posterior estimation of statistical evidence.
#' Stat Comput. 2014;24: 709-723. doi:10.1007/s11222-013-9397-1
#'
#' @return `list()`
#'
#' @keywords internal
logevidence <- function(object,
                        scheme = "power-fractions", rungs = 40,
                        steps = 5000, burnin = 0,
                        thin = 1, adapt = 0,
                        chains = 1, threads = 1,
                        seed = NA, progress = TRUE) {
  stopifnot(inherits(object, "HousepiModel"))
  object$logevidence(scheme = scheme, rungs = rungs,
                     steps = steps, burnin = burnin,
                     thin = thin, adapt = adapt,
                     chains = chains, threads = threads,
                     seed = seed, progress = progress)
}

############################################################################
############################################################################
############################################################################

#' @title Simulate household data
#' @description
#' Generate a simulated outbreak.
#' It assumes the same community structure of the loaded data.
#'
#' @param object `r rd_model_in_`
#' @param nsim `r rd_simul_nsim_`
#' @param seed `r rd_seed_`
#' @param progress `r rd_progress_`
#' @param obs_process `r rd_simul_observ_`
#' @param ... Ignored.
#'
#' @return `data.frame`
#'
#' @importFrom stats simulate
#' @export
simulate.HousepiModel <- function(object,
                                  nsim = 1, seed = NA,
                                  progress = FALSE,
                                  obs_process = NULL,
                                  ...) {
  stopifnot(inherits(object, "HousepiModel"))
  object$simulate(
    nsim = nsim, seed = seed, progress = progress,
    obs_process = obs_process
  )
}
