#
#        __
#       / /_  ____  __  __ ____ __  ____  __
#      / __ \/ __ \/ / / / ___/ _ \/ __ \/ /
#     / / / / /_/ / /_/ (__  )  __/ /_/ / /
#    /_/ /_/\____/\____/____/\___/ ____/_/
#                               /_/
#
# housepi - Epidemiological Analysis of Household Outbreaks.
#
# Author: Francesco Camaglia
# Copyright © 2025 Institut Pasteur (Paris), CNRS and INSERM.
#
# Distributed under the terms of the GNU General Public License (GPLv3).
# See the LICENSE file for details.
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
#' @param ... `r rd_mdl_flags_`
#'
#' @return `r rd_model_in_`
#'
#' @export
new_housepi_model <- function(choice = NA_character_, as_template = FALSE,
                              force = FALSE, quiet = FALSE, ...) {
  HousepiModel$new(choice = choice, as_template = as_template,
                   force = force, quiet = quiet, ...)
}

############################################################################
############################################################################
############################################################################

#' @title Add data wide
#' @description
#' Make an internal copy of the provided household data.
#'
#' @param hhmodel `r rd_model_in_`
#' @param data_wide `r rd_add_data_wide_`
#'
#' @return `r rd_model_out_`
#'
#' @export
add_data_wide  <- function(hhmodel, data_wide) {
  stopifnot(inherits(hhmodel, "HousepiModel"))
  hhmodel$clone(deep = TRUE)$add_data_wide(data_wide)
}

#' @title Add data long
#' @description
#' Male an internal copy of the provided data in long format.
#'
#' @param hhmodel `r rd_model_in_`
#' @param data_long `r rd_add_data_long_`
#'
#' @return `r rd_model_out_`
#'
#' @export
add_data_long <- function(hhmodel, data_long) {
  stopifnot(inherits(hhmodel, "HousepiModel"))
  hhmodel$clone(deep = TRUE)$add_data_long(data_long)
}

#' @title Summarise cases in household data
#' @description
#' Summarise cases as given from the statuses in the laoded household data.
#'
#' @param hhmodel `r rd_model_in_`
#' @param ... Ignored.
#'
#' @return `data.frame`
#'
#' @export
summarise_cases <- function(hhmodel, ...) {
  hhmodel$summarise_cases()
}

#' @title Summarise house sizes in household data
#' @description
#' Summarise number of members (size) in the laoded household data.
#'
#' @param hhmodel `r rd_model_in_`
#' @param ... Ignored.
#'
#' @return `data.frame`
#'
#' @export
summarise_house_sizes <- function(hhmodel, ...) {
  hhmodel$summarise_house_sizes()
}

############################################################################
############################################################################
############################################################################

#' @title Drop definitions
#' @description
#' Drop current formatted definitions.
#'
#' @param hhmodel `r rd_model_in_`
#' @param ... Ignored.
#'
#' @return `r rd_model_out_`
#'
#' @rdname drop_definitions
#' @export
drop_parameters <-  function(hhmodel, ...) {
  stopifnot(inherits(hhmodel, "HousepiModel"))
  hhmodel$clone(deep = TRUE)$drop_parameters(...)
}

#' @rdname drop_definitions
#' @export
drop_augmented <-  function(hhmodel, ...) {
  stopifnot(inherits(hhmodel, "HousepiModel"))
  hhmodel$clone(deep = TRUE)$drop_augmented(...)
}

#' @rdname drop_definitions
#' @export
drop_group_rates <-  function(hhmodel, ...) {
  stopifnot(inherits(hhmodel, "HousepiModel"))
  hhmodel$clone(deep = TRUE)$drop_group_rates(...)
}

#' @rdname drop_definitions
#' @export
drop_all <- function(hhmodel) {
  stopifnot(inherits(hhmodel, "HousepiModel"))
  hhmodel$clone(deep = TRUE)$drop_all()
}

############################################################################
############################################################################
############################################################################

#' @title Get definitions
#' @description
#' Get current formatted definitions.
#'
#' @param hhmodel `r rd_model_in_`
#' @param ... Ignored.
#'
#' @return `data.frame`
#'
#' @rdname get_definitions
#' @export
get_parameters <- function(hhmodel, ...) {
  stopifnot(inherits(hhmodel, "HousepiModel"))
  hhmodel$parameters
}

#' @title Add parameter definition
#' @description
#' Add a new parameter to houspepi model definitions.
#' If the parameter already exists, it is overwritten.
#'
#' @param hhmodel `r rd_model_in_`
#' @param name `r rd_var_name_`
#' @param value `r rd_var_value_`
#' @param const `r rd_var_const_`
#' @param discr `r rd_par_discr_`
#' @param prior `r rd_par_model_`
#' @param ... `r rd_par_args_`
#' @param min `r rd_var_min_`
#' @param max `r rd_var_max_`
#' @param proposal `r rd_par_proposal_`
#' @param proposal_scale `r rd_var_proposal_scale_`
#'
#' @return `r rd_model_out_`
#'
#' @export
add_parameter <- function(hhmodel,
                          name = NA_character_, value = NULL,
                          const = NA, discr = FALSE,
                          prior = NA_character_, ...,
                          min = NULL, max = NULL,
                          proposal = NA_character_,
                          proposal_scale = NA_real_) {
  stopifnot(inherits(hhmodel, "HousepiModel"))
  new_model <- hhmodel$clone(deep = TRUE)
  new_model$.__enclos_env__$private$.param$add_quo(
    name = name, value = rlang::enquo(value),
    const = const, discr = discr,
    prior = prior, rlang::enquos(...),
    min = rlang::enquo(min), max = rlang::enquo(max),
    proposal = proposal, proposal_scale = proposal_scale
  )
  new_model
}

#' @title Load list of definitions
#' @description
#' Load definitions as named arguments.
#'
#' @param hhmodel `r rd_model_in_`
#' @param ... `r rd_list_args_`
#' @param append `r rd_list_append_`
#'
#' @return `r rd_model_out_`
#'
#' @rdname load_definitions_list
#' @export
load_parameter_list <- function(hhmodel, ..., append = FALSE) {
  stopifnot(inherits(hhmodel, "HousepiModel"))
  new_model <- hhmodel$clone(deep = TRUE)
  new_model$.__enclos_env__$private$.param$load_quosures(
    rlang::enquos(...), append = append
  )
  new_model
}

#' @title Load template definitions
#' @description
#' Load default definitions from model template.
#'
#' @param hhmodel `r rd_model_in_`
#' @param append `r rd_tmpl_append_`
#'
#' @return `r rd_model_out_`
#'
#' @rdname load_template_definitions
#' @export
load_parameter_template <- function(hhmodel, append = FALSE) {
  stopifnot(inherits(hhmodel, "HousepiModel"))
  hhmodel$clone(deep = TRUE)$load_parameter_template(append = append)
}

############################################################################
############################################################################
############################################################################

#' @rdname get_definitions
#' @export
get_augmented <- function(hhmodel, ...) {
  stopifnot(inherits(hhmodel, "HousepiModel"))
  hhmodel$augmented
}

#' @title Add augmented definition
#' @description
#' Add a new augmented variable/function to houspepi model definitions.
#' If the augmented variablw already exists, it is updated.
#'
#' @param hhmodel `r rd_model_in_`
#' @param name `r rd_var_name_`
#' @param value `r rd_aug_val_`
#' @param status `r rd_aug_status_`
#' @param const `r rd_var_const_`
#' @param discr `r rd_aug_discr_`
#' @param func `r rd_aug_func_`
#' @param func_args `r rd_aug_funargs_`
#' @param model `r rd_aug_model_`
#' @param ... `r rd_aug_args_`
#' @param min `r rd_aug_min_`
#' @param max `r rd_aug_max_`
#' @param proposal `r rd_var_proposal_`
#' @param proposal_scale `r rd_var_proposal_scale_`
#'
#' @return `r rd_model_out_`
#'
#' @export
add_augmented <- function(hhmodel,
                          name = NA_character_, value = NULL,
                          status = NA_integer_, const = NA, discr = FALSE,
                          func = NULL, func_args = NULL,
                          model = NA_character_, ...,
                          min = NULL, max = NULL,
                          proposal = NA_character_,
                          proposal_scale = NA_real_) {
  stopifnot(inherits(hhmodel, "HousepiModel"))
  new_model <- hhmodel$clone(deep = TRUE)
  new_model$.__enclos_env__$private$.augm$add_quo(
    name = name, value = rlang::enquo(value),
    status = status, const = const, discr = discr,
    func = rlang::enquo(func), func_args = rlang::enquos(func_args),
    model = model, rlang::enquos(...),
    min = rlang::enquo(min), max = rlang::enquo(max),
    proposal = proposal, proposal_scale = proposal_scale
  )
  new_model
}

#' @rdname load_definitions_list
#' @export
load_augmented_list <- function(hhmodel, ..., append = FALSE) {
  stopifnot(inherits(hhmodel, "HousepiModel"))
  new_model <- hhmodel$clone(deep = TRUE)
  new_model$.__enclos_env__$private$.augm$load_quosures(
    rlang::enquos(...), append = append
  )
  new_model
}

#' @rdname load_template_definitions
#' @export
load_augmented_template <- function(hhmodel, append = FALSE) {
  stopifnot(inherits(hhmodel, "HousepiModel"))
  hhmodel$clone(deep = TRUE)$load_augmented_template(append = append)
}

############################################################################
############################################################################
############################################################################

#' @rdname get_definitions
#' @export
get_group_rates <- function(hhmodel, ...) {
  stopifnot(inherits(hhmodel, "HousepiModel"))
  # internal build requires deep clone to avoid inplace
  hhmodel$clone(deep = TRUE)$group_rates
}

#' @title Add group rate definition
#' @description
#' Add definitions of a group rate by given data feature(s).
#' Group rates are relative to a reference class that can be specified.
#' If no MCMC options are provided, the following configuration is assumed
#' - `prior`: "lognormal"
#' - `meanlog`: 0
#' - `sdlog`: 3
#' - `min`: 0
#' - `max` 10
#'
#' @param hhmodel `r rd_model_in_`
#' @param name `r rd_var_name_`
#' @param by `r rd_grr_by_`
#' @param prefix `r rd_grr_prefix_`
#' @param reference `r rd_grr_refer_`
#' @param independent `r rd_grr_indep_`
#' @param prior `r rd_par_model_`
#' @param ... `r rd_par_args_`
#' @param min `r rd_var_min_`
#' @param max `r rd_var_max_`
#' @param proposal `r rd_var_proposal_`
#' @param proposal_scale `r rd_var_proposal_scale_`
#'
#' @references
#' Cauchemez, et al. (2009)
#' Household Transmission of 2009 Pandemic Influenza A (H1N1) Virus in the United States.
#' New England Journal of Medicine. 361: 2619-2627. doi:10.1056/NEJMoa0905498
#'
#' @return `r rd_model_out_`
#'
#' @export
add_group_rate <-  function(hhmodel,
                            name, by, prefix = NA_character_,
                            reference = NULL, independent = TRUE,
                            prior = NA_character_, ...,
                            min = NULL, max = NULL,
                            proposal = NA_character_,
                            proposal_scale = NA_real_) {
  stopifnot(inherits(hhmodel, "HousepiModel"))
  new_model <- hhmodel$clone(deep = TRUE)
  new_model$.__enclos_env__$private$.group_rates$add_quo(
    name = name, by = by, prefix = prefix,
    reference = reference, independent = independent,
    prior = prior, rlang::enquos(...),
    min = rlang::enquo(min), max = rlang::enquo(max),
    proposal = proposal, proposal_scale = proposal_scale
  )
  new_model
}

#' @rdname load_definitions_list
#' @export
load_group_rate_list <- function(hhmodel, ..., append = FALSE) {
  stopifnot(inherits(hhmodel, "HousepiModel"))
  new_model <- hhmodel$clone(deep = TRUE)
  if (isFALSE(append)) new_model$drop_group_rates()
  new_model <- new_model
  new_model$.__enclos_env__$private$.group_rates$load_quosures(
    rlang::enquos(...), append = append
  )
  new_model
}

############################################################################
############################################################################
############################################################################

#' @rdname get_definitions
#' @export
get_infec_profile <- function(hhmodel, ...) {
  stopifnot(inherits(hhmodel, "HousepiModel"))
  hhmodel$infec_profile
}

#' @title Add infectivity profile definition
#' @description
#' Add definition of infectivity profile.
#'
#' @param hhmodel `r rd_model_in_`
#' @param name `r rd_var_name_`
#' @param status `r rd_aug_status_`
#' @param discr `r rd_inf_discr_`
#' @param ... `r rd_inf_args_`
#' @param max `r rd_aug_max_`
#'
#' @return `r rd_model_out_`
#'
#' @export
add_infec_profile <- function(hhmodel,
                              name, status = NA_integer_,
                              discr = FALSE,
                              ..., max = NULL) {
  stopifnot(inherits(hhmodel, "HousepiModel"))
  new_model <- hhmodel$clone(deep = TRUE)
  new_model$.__enclos_env__$private$.infec_profile$add_quo(
    name = name, status = status, discr = discr,
    rlang::enquos(...), max = rlang::enquo(max)
  )
  new_model
}

############################################################################
############################################################################
############################################################################

#' @title Evaluate log-likelihood of housepi model
#' @description
#' Compute the log-likelihood of a `HousepiModel` object.
#'
#' @param hhmodel `r rd_model_in_`
#' @param house_id `r rd_house_id_`
#' If `NA` (default) it computes for the whole data.
#' @param seed `r rd_seed_`
#'
#' @return `numeric(1)`
#'
#' @export
evaluate_log_lik <- function(hhmodel,
                             house_id = NA_character_,
                             seed = NA_integer_) {
  stopifnot(inherits(hhmodel, "HousepiModel"))
  hhmodel$evaluate_log_lik(house_id = house_id, seed = seed)
}

#' @title Run inference for housepi model
#' @description
#' Run inference with Bayesian Markov Chain Monte Carlo.
#'
#' @param hhmodel `r rd_model_in_`
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
run_inference <- function(hhmodel,
                          steps = 200, burnin = 0,
                          thin = 1, adapt = 0,
                          chains = 1, threads = 1,
                          id_draws = NULL, seed = NA_integer_,
                          progress = TRUE) {
  stopifnot(inherits(hhmodel, "HousepiModel"))
  hhmodel$run_inference(steps = steps, burnin = burnin,
                        thin = thin, adapt = adapt,
                        chains = chains, threads = threads,
                        id_draws = id_draws, seed = seed,
                        progress = progress)
}

#' @title Run simulation for housepi model
#' @description
#' Run simulation to generate a household survey.
#' It assumes the same community structure of the loaded data.
#'
#' @param hhmodel `r rd_model_in_`
#' @param nsim `r rd_simul_nsim_`
#' @param seed `r rd_seed_`
#' @param progress `r rd_progress_`
#' @param obs_process `r rd_simul_observ_`
#'
#' @return `data.frame`
#'
#' @export
run_simulation <- function(hhmodel,
                           nsim = 1, seed = NA_integer_,
                           progress = FALSE, obs_process = NULL) {
  stopifnot(inherits(hhmodel, "HousepiModel"))
  hhmodel$run_simulation(nsim = nsim, seed = seed,
                         progress = progress,
                         obs_process = obs_process)
}

#' @title Estimate logevidence of housepi model
#' @description
#' Estimate the logarithm of the marginal likelihood (aka evidence).
#' It uses the method of power posteriors with different available schemes.
#'
#' @param hhmodel `r rd_model_in_`
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
#' @noRd
# logevidence <- function(hhmodel,
#                         scheme = "power-fractions", rungs = 40,
#                         steps = 5000, burnin = 0,
#                         thin = 1, adapt = 0,
#                         chains = 1, threads = 1,
#                         seed = NA_integer_, progress = TRUE) {
#   stopifnot(inherits(hhmodel, "HousepiModel"))
#   hhmodel$logevidence(scheme = scheme, rungs = rungs,
#                       steps = steps, burnin = burnin,
#                       thin = thin, adapt = adapt,
#                       chains = chains, threads = threads,
#                       seed = seed, progress = progress)
# }
