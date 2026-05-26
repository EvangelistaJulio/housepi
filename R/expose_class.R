## housepi - Epidemiological Analysis of Household Outbreaks.
##
## Copyright © 2024-2026 Institut Pasteur (Paris), CNRS and INSERM.
## Distributed under the terms of the GNU General Public License (GPLv3).
## See the LICENSE file for details.
##

#' @title Create a new household-transmission model instance
#' @description
#' Create a new object of class 'HousepiModel'
#' used to run inference and simulation routines.
#' @param data `r rd_add_wide_`
#' @param model `r rd_mdl_choice_`
#' @param as_template `r rd_mdl_astmpl_`
#' @param force `r rd_mdl_force_`
#' @param quiet `r rd_mdl_quiet_`
#' @param ... `r rd_mdl_flags_`
#' @return `r rd_model_in_`
#' @examples
#' # build a model from a template using a built-in dataset
#' hhmodel <- new_housepi_model(flu_1918_frederick, "chain_binomial")
#' print(hhmodel)
#'
#' # 'blank' template provides a minimal model for testing
#' new_housepi_model(model = "blank")
#' @export
new_housepi_model <- function(data = NULL, model = NA_character_,
                              as_template = FALSE,
                              force = FALSE, quiet = FALSE, ...) {
  HousepiModel$new(data = data, model = model,
                   as_template = as_template,
                   force = force, quiet = quiet, ...)
}

############################################################################
############################################################################
############################################################################

# Check hosuepi model.
stopifnothousepi <- function(hhmodel) {
  if (!inherits(hhmodel, "HousepiModel")) stop_h(
    "'hhmodel' must be of class 'HousepiModel'",
    "Create one with new_housepi_model()."
  )
}

#' @title Attach individual-level long-format data to a housepi model
#' @description
#' Make an internal copy of the provided data in long format.
#' @param hhmodel `r rd_model_in_`
#' @param data_long `r rd_add_data_long_`
#' @return `r rd_model_out_`
#' @export
add_data_long <- function(hhmodel, data_long) {
  stopifnothousepi(hhmodel)
  new_obj <- hhmodel$copy()
  new_obj$add_data_long(data_long)
  new_obj$.__enclos_env__$private$rebuild_()
  new_obj
}

#' @title Summarise observed cases in the loaded household data
#' @description
#' Summarise cases as given from the statuses in the loaded household data.
#' @param hhmodel `r rd_model_in_`
#' @return `data.frame`
#' @examples
#' hhmodel <- new_housepi_model(flu_1918_frederick, "chain_binomial")
#' summarise_house_cases(hhmodel)
#' @export
summarise_house_cases <- function(hhmodel) {
  hhmodel$summarise_house_cases()
}

#' @rdname summarise_house_cases
#' @export
summarize_house_cases <- summarise_house_cases

#' @title Summarise household sizes in the loaded household data
#' @description
#' Summarise number of members (size) in the loaded household data.
#' @param hhmodel `r rd_model_in_`
#' @return `data.frame`
#' @examples
#' hhmodel <- new_housepi_model(flu_1918_frederick, "chain_binomial")
#' summarise_house_sizes(hhmodel)
#' @export
summarise_house_sizes <- function(hhmodel) {
  hhmodel$summarise_house_sizes()
}

#' @rdname summarise_house_sizes
#' @export
summarize_house_sizes <- summarise_house_sizes

############################################################################
############################################################################
############################################################################

#' @title Drop currently loaded definitions from a housepi model
#' @description Remove current definitions from the housepi model.
#' Each variant clears a single category (parameters, augmented variables,
#' group rates, contact rates, infectivity profile) while `drop_all()` clears
#' them all.
#' @param hhmodel `r rd_model_in_`
#' @return `r rd_model_out_`
#' @examples
#' hhmodel <- new_housepi_model(flu_1918_frederick, "chain_binomial") |>
#'   add_parameter_template()
#' hhmodel |>
#'   drop_parameters() |>
#'   get_parameters()
#' @rdname drop_definitions
#' @export
drop_parameters <-  function(hhmodel) {
  stopifnothousepi(hhmodel)
  new_obj <- hhmodel$copy()
  new_obj$drop_parameters()
  new_obj
}

#' @rdname drop_definitions
#' @export
drop_augmented <-  function(hhmodel) {
  stopifnothousepi(hhmodel)
  new_obj <- hhmodel$copy()
  new_obj$drop_augmented()
  new_obj
}

#' @rdname drop_definitions
#' @export
drop_group_rates <-  function(hhmodel) {
  stopifnothousepi(hhmodel)
  new_obj <- hhmodel$copy()
  new_obj$drop_group_rates()
  new_obj
}

#' @rdname drop_definitions
#' @export
drop_contact_rates <-  function(hhmodel) {
  stopifnothousepi(hhmodel)
  new_obj <- hhmodel$copy()
  new_obj$drop_contact_rates()
  new_obj
}

#' @rdname drop_definitions
#' @export
drop_infectivity_profile <-  function(hhmodel) {
  stopifnothousepi(hhmodel)
  new_obj <- hhmodel$copy()
  new_obj$drop_infectivity_profile()
  new_obj
}

#' @rdname drop_definitions
#' @export
drop_all <- function(hhmodel) {
  stopifnothousepi(hhmodel)
  new_obj <- hhmodel$copy()
  new_obj$drop_all()
  new_obj
}

############################################################################
############################################################################
############################################################################

#' @title Retrieve currently loaded data and definitions from a housepi model
#' @description
#' Return the definitions currently stored in the housepi model.
#' Each variant returns one category (data, parameters, augmented variables,
#' group rates, contact rates, infectivity profile).
#' @param hhmodel `r rd_model_in_`
#' @return `data.frame`
#' @examples
#' hhmodel <- new_housepi_model(flu_1918_frederick, "chain_binomial") |>
#'   add_parameter_template()
#' get_data_wide(hhmodel) |>
#'   head()
#' get_parameters(hhmodel)
#' @rdname get_definitions
#' @export
get_data_wide <- function(hhmodel) {
  stopifnothousepi(hhmodel)
  hhmodel$data_wide
}

#' @rdname get_definitions
#' @export
get_data_long <- function(hhmodel) {
  stopifnothousepi(hhmodel)
  hhmodel$data_long
}

#' @rdname get_definitions
#' @export
get_parameters <- function(hhmodel) {
  stopifnothousepi(hhmodel)
  hhmodel$parameters
}

#' @rdname get_definitions
#' @export
get_augmented <- function(hhmodel) {
  stopifnothousepi(hhmodel)
  hhmodel$augmented
}

#' @rdname get_definitions
#' @export
get_group_rates <- function(hhmodel) {
  stopifnothousepi(hhmodel)
  # NOTE: internal build requires deep clone to avoid in-place
  hhmodel$copy()$group_rates
}

#' @rdname get_definitions
#' @export
get_contact_rates <- function(hhmodel) {
  stopifnothousepi(hhmodel)
  # NOTE: internal build requires deep clone to avoid in-place
  hhmodel$copy()$contact_rates
}

#' @rdname get_definitions
#' @export
get_infec_profile <- function(hhmodel) {
  stopifnothousepi(hhmodel)
  hhmodel$infec_profile
}

############################################################################
############################################################################
############################################################################

#' @title Add or overwrite a single parameter definition
#' @description
#' Add a new parameter to housepi model definitions.
#' If the parameter already exists, it is overwritten.
#' @param hhmodel `r rd_model_in_`
#' @param name `r rd_var_name_`
#' @param value `r rd_var_value_`
#' @param const `r rd_var_const_`
#' @param discr `r rd_par_discr_`
#' @param prior `r rd_par_model_`
#' @param ... `r rd_par_args_`
#' @param min `r rd_var_min_`
#' @param func `r rd_aug_func_`
#' @param func_args `r rd_aug_funargs_`
#' @param max `r rd_var_max_`
#' @param proposal `r rd_par_proposal_`
#' @param proposal_scale `r rd_proposal_scale_`
#' @return `r rd_model_out_`
#' @examples
#' hhmodel <- flu_1918_frederick |>
#'   new_housepi_model("chain_binomial") |>
#'   # alpha : community transmission rate, with uniform prior
#'   add_parameter("alpha", prior = "uniform", proposal = "lognormal") |>
#'   # beta : person-to-person transmission rate, fixed to a constant
#'   add_parameter("beta", value = 0.3, const = TRUE)
#' @export
add_parameter <- function(hhmodel,
                          name = NA_character_, value = NULL,
                          const = NA, discr = FALSE,
                          prior = NA_character_, ...,
                          func = NULL, func_args = NULL,
                          min = NULL, max = NULL,
                          proposal = NA_character_,
                          proposal_scale = NA_real_) {
  stopifnothousepi(hhmodel)
  new_obj <- hhmodel$copy()
  new_obj$.__enclos_env__$private$param_$add_quo(
    name = name, value = rlang::enquo(value),
    const = const, discr = discr,
    prior = prior, rlang::enquos(...),
    func = rlang::enquo(func), func_args = rlang::enquos(func_args),
    min = rlang::enquo(min), max = rlang::enquo(max),
    proposal = proposal, proposal_scale = proposal_scale,
    .assert = TRUE
  )
  new_obj$.__enclos_env__$private$rebuild_()
  new_obj
}

#' @title Load default definitions from the model template
#' @description
#' Load default parameter or augmented definitions from model template.
#' @param hhmodel `r rd_model_in_`
#' @param append `r rd_tmpl_append_`
#' @return `r rd_model_out_`
#' @examples
#' hhmodel <- new_housepi_model(flu_1918_frederick, "chain_binomial") |>
#'   add_parameter_template()
#' get_parameters(hhmodel)
#' @rdname add_template_definitions
#' @export
add_parameter_template <- function(hhmodel, append = FALSE) {
  stopifnothousepi(hhmodel)
  new_obj <- hhmodel$copy()
  new_obj$add_parameter_template(append = append)
  new_obj
}

############################################################################
############################################################################
############################################################################

#' @title Add or update an augmented variable definition
#' @description
#' Add a new augmented variable/function to housepi model definitions.
#' If the augmented variable already exists, it is updated.
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
#' @param proposal_scale `r rd_proposal_scale_`
#' @return `r rd_model_out_`
#' @examples
#' hhmodel <- blank_house |>
#'   new_housepi_model("augmented_example") |>
#'   # define an incubation period drawn from a gamma distribution
#'   add_augmented("incub_period",
#'                 model = "gamma", shape = 2, scale = 1, max = 14)
#' @export
add_augmented <- function(hhmodel,
                          name = NA_character_, value = NULL,
                          status = NA_integer_, const = NA, discr = FALSE,
                          func = NULL, func_args = NULL,
                          model = NA_character_, ...,
                          min = NULL, max = NULL,
                          proposal = NA_character_,
                          proposal_scale = NA_real_) {
  stopifnothousepi(hhmodel)
  new_obj <- hhmodel$copy()
  new_obj$.__enclos_env__$private$augm_$add_quo(
    name = name, value = rlang::enquo(value),
    status = status, const = const, discr = discr,
    func = rlang::enquo(func), func_args = rlang::enquos(func_args),
    model = model, rlang::enquos(...),
    min = rlang::enquo(min), max = rlang::enquo(max),
    proposal = proposal, proposal_scale = proposal_scale,
    .assert = TRUE
  )
  new_obj$.__enclos_env__$private$rebuild_()
  new_obj
}

#' @rdname add_template_definitions
#' @export
add_augmented_template <- function(hhmodel, append = FALSE) {
  stopifnothousepi(hhmodel)
  new_obj <- hhmodel$copy()
  new_obj$add_augmented_template(append = append)
  new_obj
}

############################################################################
############################################################################
############################################################################

#' @title Add a stratified group rate definition
#' @description `r rd_grr_desc_`
#' @param hhmodel `r rd_model_in_`
#' @param name `r rd_var_name_`
#' @param by `r rd_grr_by_`
#' @param reference `r rd_grr_refer_`
#' @param joint `r rd_grr_joint_`
#' @param strata `r rd_grr_strata_`
#' @param params `r rd_grr_params_`
#' @param prior `r rd_par_model_`
#' @param ... `r rd_par_args_`
#' @param min `r rd_var_min_`
#' @param max `r rd_var_max_`
#' @param proposal `r rd_var_proposal_`
#' @param proposal_scale `r rd_proposal_scale_`
#' @param prefix `r rd_grr_prefix_`
#' @references
#' \insertRef{cauchemez_household_2009}{housepi}
#' @return `r rd_model_out_`
#' @export
add_group_rate <-  function(hhmodel,
                            name, by = NA_character_,
                            reference = NULL, joint = FALSE,
                            strata = NULL, params = NULL,
                            prior = NA_character_, ...,
                            min = NA_real_, max = NA_real_,
                            proposal = NA_character_,
                            proposal_scale = NA_real_,
                            prefix = NA_character_) {
  stopifnothousepi(hhmodel)
  new_obj <- hhmodel$copy()
  new_obj$.__enclos_env__$private$group_rates_$add_quo(
    name = name, by = by,
    reference = reference, joint = joint,
    strata = strata, params = params,
    prior = prior, rlang::enquos(...),
    min = rlang::enquo(min), max = rlang::enquo(max),
    proposal = proposal, proposal_scale = proposal_scale,
    prefix = prefix
  )
  new_obj$.__enclos_env__$private$rebuild_()
  new_obj
}

############################################################################
############################################################################
############################################################################

#' @title Add a pairwise contact rate definition
#' @description `r rd_ctr_desc_`
#' @param hhmodel `r rd_model_in_`
#' @param by `r rd_ctr_by_`
#' @param reference `r rd_ctr_refer_`
#' @param joint `r rd_ctr_joint_`
#' @param symmetric `r rd_ctr_symmetric_`
#' @param strata `r rd_ctr_strata_`
#' @param params `r rd_ctr_params_`
#' @param prior `r rd_par_model_`
#' @param ... `r rd_par_args_`
#' @param min `r rd_var_min_`
#' @param max `r rd_var_max_`
#' @param proposal `r rd_var_proposal_`
#' @param proposal_scale `r rd_proposal_scale_`
#' @param prefix `r rd_ctr_prefix_`
#' @return `r rd_model_out_`
#' @export
add_contact_rate <- function(hhmodel, by = NA_character_,
                             reference = NULL,
                             joint = FALSE, symmetric = TRUE,
                             strata = NULL, params = NULL,
                             prior = NA_character_, ...,
                             min = NA_real_, max = NA_real_,
                             proposal = NA_character_,
                             proposal_scale = NA_real_,
                             prefix = "ctr") {
  stopifnothousepi(hhmodel)
  new_obj <- hhmodel$copy()
  new_obj$.__enclos_env__$private$contact_rates_$add_quo(
    by = by, reference = reference,
    joint = joint, symmetric = symmetric,
    strata = strata, params = params,
    prior = prior, rlang::enquos(...),
    min = rlang::enquo(min), max = rlang::enquo(max),
    proposal = proposal, proposal_scale = proposal_scale,
    prefix = prefix
  )
  new_obj$.__enclos_env__$private$rebuild_()
  new_obj
}

############################################################################
############################################################################
############################################################################

#' @title Add an infectivity profile definition
#' @description
#' Add definition of infectivity profile.
#' @param hhmodel `r rd_model_in_`
#' @param func `r rd_var_name_`
#' @param status `r rd_aug_status_`
#' @param discr `r rd_inf_discr_`
#' @param ... `r rd_inf_args_`
#' @param max `r rd_aug_max_`
#' @return `r rd_model_out_`
#' @export
add_infec_profile <- function(hhmodel,
                              func, status = NA_integer_,
                              discr = FALSE,
                              ..., max = NULL) {
  stopifnothousepi(hhmodel)
  new_obj <- hhmodel$copy()
  new_obj$.__enclos_env__$private$infec_profile_$add_quo(
    func = func, status = status, discr = discr,
    rlang::enquos(...), max = rlang::enquo(max),
    .assert = TRUE
  )
  new_obj$.__enclos_env__$private$rebuild_()
  new_obj
}

############################################################################
############################################################################
############################################################################

#' @title Load multiple definitions at once via named arguments
#' @description
#' Load definitions as named arguments.
#' Each variant targets one definition category (parameters, augmented
#' variables, group rates, infectivity profile).
#' @param hhmodel `r rd_model_in_`
#' @param ... `r rd_list_args_`
#' @param append `r rd_list_append_`
#' @return `r rd_model_out_`
#' @examples
#' # define some constant paramters
#' flu_1918_frederick |>
#'   new_housepi_model("chain_binomial") |>
#'   add_parameter_list(alpha = 0.05, beta = 0.3,
#'                      gamma = 1, size_ref = 4) |>
#'   get_parameters()
#' @rdname add_definitions_list
#' @export
add_parameter_list <- function(hhmodel, ..., append = FALSE) {
  stopifnothousepi(hhmodel)
  new_obj <- hhmodel$copy()
  new_obj$.__enclos_env__$private$param_$add_quo_list(
    rlang::enquos(...), append = append
  )
  new_obj$.__enclos_env__$private$rebuild_()
  new_obj
}

#' @rdname add_definitions_list
#' @export
add_augmented_list <- function(hhmodel, ..., append = FALSE) {
  stopifnothousepi(hhmodel)
  new_obj <- hhmodel$copy()
  new_obj$.__enclos_env__$private$augm_$add_quo_list(
    rlang::enquos(...), append = append
  )
  new_obj$.__enclos_env__$private$rebuild_()
  new_obj
}

#' @rdname add_definitions_list
#' @export
add_group_rate_list <- function(hhmodel, ..., append = FALSE) {
  stopifnothousepi(hhmodel)
  new_obj <- hhmodel$copy()
  # NOTE: dropping group rates might cause parameter drops
  if (isFALSE(append)) new_obj$drop_group_rates()
  new_obj$.__enclos_env__$private$group_rates_$add_quo_list(
    rlang::enquos(...), append = append
  )
  new_obj$.__enclos_env__$private$rebuild_()
  new_obj
}

#' @rdname add_definitions_list
#' @export
add_infec_profile_list <- function(hhmodel, ..., append = FALSE) {
  stopifnothousepi(hhmodel)
  new_obj <- hhmodel$copy()
  new_obj$.__enclos_env__$private$infec_profile_$add_quo_list(
    rlang::enquos(...), append = append
  )
  new_obj$.__enclos_env__$private$rebuild_()
  new_obj
}

############################################################################
############################################################################
############################################################################

#' @title Evaluate the log-likelihood of a housepi model at current values
#' @description
#' Compute the log-likelihood of a `HousepiModel` object.
#' @param hhmodel `r rd_model_in_`
#' @param house_id `r rd_house_id_`
#' If `NA` (default) it computes for the whole data.
#' @return `numeric(1)`
#' @examples
#' flu_1918_frederick |>
#'   new_housepi_model("chain_binomial") |>
#'   add_parameter_template() |>
#'   evaluate_log_lik()
#' @export
evaluate_log_lik <- function(hhmodel, house_id = NA_character_) {
  stopifnothousepi(hhmodel)
  hhmodel$evaluate_log_lik(house_id = house_id)
}

#' @title Run Bayesian MCMC inference on a housepi model
#' @description
#' Run inference with Bayesian Markov Chain Monte Carlo.
#' @param hhmodel `r rd_model_in_`
#' @param steps `r rd_smpl_steps_`
#' @param burn `r rd_smpl_burn_`
#' @param thin `r rd_smpl_thin_`
#' @param adapt `r rd_smpl_adapt_`
#' @param chains `r rd_smpl_chain_`
#' @param threads `r rd_thread_`
#' @param id_draws `r rd_smpl_iddrw_`
#' @param progress `r rd_progress_`
#' @return An object of class `HousepiPosteriorDraws`.
#' @examples
#' flu_1918_frederick |>
#'   new_housepi_model("chain_binomial") |>
#'   add_parameter_template() |>
#'   run_inference(steps = 50, burn = 10, thin = 2,
#'                 chains = 2, progress = FALSE) |>
#'   get_parameter_draws() |>
#'   posterior::summarize_draws()
#' @export
run_inference <- function(hhmodel,
                          steps = 1000, burn = 0,
                          thin = 5, adapt = 500,
                          chains = 4, threads = 1,
                          id_draws = NULL,
                          progress = TRUE) {
  stopifnothousepi(hhmodel)
  hhmodel$run_inference(steps = steps, burn = burn,
                        thin = thin, adapt = adapt,
                        chains = chains, threads = threads,
                        id_draws = id_draws,
                        progress = progress)
}

#' @title Simulate a household survey from a housepi model
#' @description
#' Run simulation to generate a household survey.
#' It assumes the same community structure of the loaded data.
#' @param hhmodel `r rd_model_in_`
#' @param obs_process `r rd_simul_observ_`
#' @param obs_process_long `r rd_simul_obs_long_`
#' @param param_table `r rd_param_table_`
#' @param nrep `r rd_simul_nrep_`
#' @return `data.frame` or `list` (when `obs_process_long` is provided)
#' @examples
#' flu_1918_frederick |>
#'   new_housepi_model("chain_binomial") |>
#'   add_parameter_template() |>
#'   run_simulation() |>
#'   head()
#' @export
run_simulation <- function(hhmodel,
                           obs_process = NULL,
                           obs_process_long = NULL,
                           param_table = NULL,
                           nrep = NA_integer_) {
  stopifnothousepi(hhmodel)
  hhmodel$run_simulation(obs_process = obs_process,
                         obs_process_long = obs_process_long,
                         param_table = param_table,
                         nrep = nrep)
}

#' @title Estimate the log-evidence (marginal log-likelihood) of a housepi model
#' @description
#' Estimate the logarithm of the marginal likelihood (aka evidence).
#' It uses the method of power posteriors with different available schemes.
#' @param hhmodel `r rd_model_in_`
#' @param scheme `r rd_evi_scheme_`
#' @param rungs `r rd_evi_rungs_`
#' @param steps `r rd_smpl_steps_`
#' @param burn `r rd_smpl_burn_`
#' @param thin `r rd_smpl_thin_`
#' @param adapt `r rd_smpl_adapt_`
#' @param chains `r rd_smpl_chain_`
#' @param threads `r rd_thread_`
#' @param progress `r rd_progress_`
#' @references
#' \insertRef{friel_improving_2014}{housepi}
#' @return `list()`
#' @noRd
# logevidence <- function(hhmodel,
#                         scheme = "power-fractions", rungs = 40,
#                         steps = 5000, burn = 0,
#                         thin = 1, adapt = 0,
#                         chains = 1, threads = 1,
#                         progress = TRUE) {
#   stopifnothousepi(hhmodel)
#   hhmodel$logevidence(scheme = scheme, rungs = rungs,
#                       steps = steps, burn = burn,
#                       thin = thin, adapt = adapt,
#                       chains = chains, threads = threads,
#                       progress = progress)
# }
