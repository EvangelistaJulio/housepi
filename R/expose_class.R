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
#' @param escaped `r rd_mdl_escaped_`
#' @param ... `r rd_mdl_flags_`
#' @return `r rd_model_in_`
#' @examples
#' # build a model from a template using a built-in dataset
#' hhmodel <- new_housepi_model(flu_1918_frederick, "final_size")
#' print(hhmodel)
#'
#' # 'test' template provides a minimal model for testing
#' new_housepi_model(model = "test")
#' @export
new_housepi_model <- function(data = NULL, model = NA_character_,
                              as_template = FALSE,
                              force = FALSE, quiet = FALSE,
                              escaped = "escaped", ...) {
  HousepiModel$new(data = data, model = model,
                   as_template = as_template,
                   force = force, quiet = quiet, escaped = escaped, ...)
}

############################################################################
############################################################################
############################################################################

# Check hosuepi model.
stopifnothousepi <- function(hhmodel) {
  if (!inherits(hhmodel, "HousepiModel")) stop_h(
    "{.arg hhmodel} must be of class {.cls HousepiModel}.",
    "Create one with {.fn new_housepi_model}."
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
  pvt <- new_obj$.__enclos_env__$private
  pvt$generator_$import_data_long(data_long)
  pvt$rebuild_()
  new_obj
}

#' @title Summarise observed cases in the loaded household data
#' @description
#' Summarise cases as given from the statuses in the loaded household data.
#' @param hhmodel `r rd_model_in_`
#' @return `r rd_summ_cases_`
#' @examples
#' hhmodel <- new_housepi_model(flu_1918_frederick, "final_size")
#' summarise_house_cases(hhmodel)
#' @export
summarise_house_cases <- function(hhmodel) {
  stopifnothousepi(hhmodel)
  pvt <- hhmodel$.__enclos_env__$private
  .engine$summarise_house_cases(pvt$generator_$mcmc_ptr)
}

#' @rdname summarise_house_cases
#' @export
summarize_house_cases <- summarise_house_cases

#' @title Summarise household sizes in the loaded household data
#' @description
#' Summarise number of members (size) in the loaded household data.
#' @param hhmodel `r rd_model_in_`
#' @return `r rd_summ_sizes_`
#' @examples
#' hhmodel <- new_housepi_model(flu_1918_frederick, "final_size")
#' summarise_house_sizes(hhmodel)
#' @export
summarise_house_sizes <- function(hhmodel) {
  stopifnothousepi(hhmodel)
  pvt <- hhmodel$.__enclos_env__$private
  .engine$summarise_house_sizes(pvt$generator_$mcmc_ptr)
}

#' @rdname summarise_house_sizes
#' @export
summarize_house_sizes <- summarise_house_sizes

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
#' hhmodel <- new_housepi_model(flu_1918_frederick, "final_size") |>
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

definition_helper_ <- c(
  i = "Supply a {.arg name} for a single definition.",
  i = "Provide multiple definitions as {.code name = list()} pairs.",
  i = "Otherwise import an external list using {.args import}."
)

# Guard the inline (`...`) bulk form.
assert_bulk_defs_ <- function(dots) {
  if (length(dots) == 0L) stop_h(definition_helper_)
  dots
}

# Build the quosure list for the explicit `import` bulk form.
# Reject `import` combined with single-mode or inline `...` pairs.
defs_to_quos_ <- function(import, key, dots, env) {
  if (!is.null(key) || length(dots) > 0L) stop_h(
    "Mixing imported definitions with single and multiple definitions",
    definition_helper_
  )
  convert_list_to_quosures(import, env)
}

# Load default templated definitions shipped with the package.
add_template_defs_ <- function(pvt, family, append) {
  holder <- if (family == "parameter") pvt$param_ else pvt$augm_
  if (pvt$generator_$is_test_model()) {
    warn_h("No template {family} definitions available for test model.")
  } else {
    defs <- get_template_def_list(pvt$generator_$source_dir, family)
    holder$add_quo_list(convert_list_to_quosures(defs), append = isTRUE(append))
  }
  pvt$rebuild_()
}

rd_add_modes_ <- "There three mutually exclusive definition modes:
\\describe{
  \\item{single}{
  Supply a `name` to adds a single definition using the explicit arguments;
  }
  \\item{multiple}{
  Pass pairs `name = list()` (or `name = value`) for multiple definitions;
  }
  \\item{import}{
  Pass a pre-built named list of definitions through `import`.
  }
}"
rd_add_details_ <- "Definition names colliding with reserved argument
are deprecated.
The can be supplied only through single mode."

#' @title Add or overwrite parameter definitions
#' @description
#' Add parameter(s) to a housepi model.
#' If a parameter with the same name exists, the former gets overwritten.
#' `r rd_add_modes_`
#' @details
#' `r rd_add_details_`
#' @param hhmodel `r rd_model_in_`
#' @param name `r rd_var_name_`
#' @param value `r rd_var_value_`
#' @param const `r rd_var_const_`
#' @param discr `r rd_par_discr_`
#' @param prior `r rd_var_prior_`
#' @param ... `r rd_par_args_`
#' In mode 'multiple', the `name = value` definitions.
#' @param min `r rd_var_min_`
#' @param func `r rd_aug_func_`
#' @param func_args `r rd_aug_funargs_`
#' @param max `r rd_var_max_`
#' @param proposal `r rd_par_proposal_`
#' @param proposal_scale `r rd_proposal_scale_`
#' @param import `r rd_import_`
#' @param append `r rd_list_append_`
#' @return `r rd_model_out_`
#' @examples
#' # the single, inline-bulk and list-bulk forms compose freely in a pipe
#' hhmodel <- flu_1918_frederick |>
#'   new_housepi_model("final_size") |>
#'   # single: alpha, with a uniform prior and a lognormal proposal
#'   add_parameter("alpha", prior = "uniform", proposal = "lognormal") |>
#'   # single: beta, fixed to a constant
#'   add_parameter("beta", value = 0.3, const = TRUE) |>
#'   # inline bulk: several definitions as `name = value` pairs
#'   add_parameter(gamma = 1, size_ref = 4)
#'
#' # list bulk: load a pre-built named list of definitions through `import`
#' my_import <- list(gamma = 1, size_ref = 4)
#' hhmodel <- flu_1918_frederick |>
#'   new_housepi_model("final_size") |>
#'   add_parameter(import = my_import)
#' @export
add_parameter <- function(hhmodel,
                          name = NULL, value = NULL, ...,
                          const = NA, discr = FALSE,
                          prior = NA_character_,
                          func = NULL, func_args = NULL,
                          min = NULL, max = NULL,
                          proposal = NA_character_,
                          proposal_scale = NA_real_,
                          import = NULL, append = TRUE) {
  stopifnothousepi(hhmodel)
  new_obj <- hhmodel$copy()
  pvt <- new_obj$.__enclos_env__$private
  if (!is.null(import)) pvt$param_$add_quo_list(
    defs_to_quos_(import, name, rlang::enquos(...), rlang::caller_env()),
    append = append
  )
  else if (!is.null(name)) pvt$param_$add_quo(
    name = name, value = rlang::enquo(value),
    const = const, discr = discr,
    prior = prior, rlang::enquos(...),
    func = rlang::enquo(func), func_args = rlang::enquos(func_args),
    min = rlang::enquo(min), max = rlang::enquo(max),
    proposal = proposal, proposal_scale = proposal_scale,
    .assert = TRUE
  )
  else pvt$param_$add_quo_list(
    assert_bulk_defs_(rlang::enquos(...)),
    append = append
  )
  pvt$rebuild_()
  new_obj
}

#' @title Load default definitions from the model template
#' @description
#' Load default parameter or augmented definitions from model template.
#' @param hhmodel `r rd_model_in_`
#' @param append `r rd_tmpl_append_`
#' @return `r rd_model_out_`
#' @examples
#' hhmodel <- new_housepi_model(flu_1918_frederick, "final_size") |>
#'   add_parameter_template()
#' get_parameters(hhmodel)
#' @rdname add_template_definitions
#' @export
add_parameter_template <- function(hhmodel, append = TRUE) {
  stopifnothousepi(hhmodel)
  new_obj <- hhmodel$copy()
  add_template_defs_(new_obj$.__enclos_env__$private, "parameter", append)
  new_obj
}

############################################################################
############################################################################
############################################################################

#' @title Add or update augmented variable definitions
#' @description
#' Add augmented variables/functions to housepi model.
#' If an augmented definitions with the same name and status exists,
#' the former gets overwritten.
#' `r rd_add_modes_`
#' @details
#' `r rd_add_details_`
#' @param hhmodel `r rd_model_in_`
#' @param name `r rd_var_name_`
#' @param value `r rd_aug_val_`
#' @param status `r rd_aug_status_`
#' @param const `r rd_var_const_`
#' @param discr `r rd_aug_discr_`
#' @param func `r rd_aug_func_`
#' @param func_args `r rd_aug_funargs_`
#' @param prior `r rd_var_prior_`
#' @param ... `r rd_aug_args_`
#' In mode 'multiple', the `name = value` definitions.
#' @param min `r rd_aug_min_`
#' @param max `r rd_aug_max_`
#' @param proposal `r rd_var_proposal_`
#' @param proposal_scale `r rd_proposal_scale_`
#' @param import `r rd_import_`
#' @param append `r rd_list_append_`
#' @return `r rd_model_out_`
#' @examples
#' hhmodel <- dummy_household |>
#'   new_housepi_model("test") |>
#'   # define an incubation period drawn from a gamma distribution
#'   add_augmented("incub_period",
#'                 prior = "gamma", shape = 2, scale = 1, max = 14)
#' @export
add_augmented <- function(hhmodel,
                          name = NULL, value = NULL, ...,
                          status = NA_integer_, const = NA, discr = FALSE,
                          func = NULL, func_args = NULL,
                          prior = NA_character_,
                          min = NULL, max = NULL,
                          proposal = NA_character_,
                          proposal_scale = NA_real_,
                          import = NULL, append = TRUE) {
  stopifnothousepi(hhmodel)
  new_obj <- hhmodel$copy()
  pvt <- new_obj$.__enclos_env__$private
  if (!is.null(import)) pvt$augm_$add_quo_list(
    defs_to_quos_(import, name, rlang::enquos(...), rlang::caller_env()),
    append = append, status = status
  )
  else if (!is.null(name)) pvt$augm_$add_quo(
    name = name, value = rlang::enquo(value),
    status = status, const = const, discr = discr,
    func = rlang::enquo(func), func_args = rlang::enquos(func_args),
    prior = prior, rlang::enquos(...),
    min = rlang::enquo(min), max = rlang::enquo(max),
    proposal = proposal, proposal_scale = proposal_scale,
    .assert = TRUE
  )
  else pvt$augm_$add_quo_list(
    assert_bulk_defs_(rlang::enquos(...)),
    append = append, status = status
  )
  pvt$rebuild_()
  new_obj
}

#' @rdname add_template_definitions
#' @export
add_augmented_template <- function(hhmodel, append = TRUE) {
  stopifnothousepi(hhmodel)
  new_obj <- hhmodel$copy()
  add_template_defs_(new_obj$.__enclos_env__$private, "augmented", append)
  new_obj
}

#' @title Add an individual observation
#' @description
#' Add observation for positive status individual associated to data \eqn{x}.
#' If model `obs` is provided, adds to household likelihood during evaluations.
#' Sample \eqn{x} from the model during simulations.
#' If observation is deterministic it is considered as a datum.
#' Note that indicator contribution to log-likelihood are expressily forbidden.
#' `r rd_add_modes_`
#' @details
#' `r rd_add_details_`
#' @param hhmodel `r rd_model_in_`
#' @param name `r rd_var_name_`
#' @param obs `r rd_obs_dist_`
#' @param ... `r rd_obs_args_`
#' In mode 'multiple', the `name = value` definitions.
#' @param status `r rd_aug_status_`
#' @param discr `r rd_aug_discr_`
#' @param min `r rd_aug_min_`
#' @param max `r rd_aug_max_`
#' @param import `r rd_import_`
#' @param append `r rd_list_append_`
#' @return `r rd_model_out_`
#' @export
add_observation <- function(hhmodel, name = NULL,
                            obs = NA_character_, ...,
                            status = NA_integer_, discr = FALSE,
                            min = NULL, max = NULL,
                            import = NULL, append = TRUE) {
  stopifnothousepi(hhmodel)
  new_obj <- hhmodel$copy()
  pvt <- new_obj$.__enclos_env__$private
  if (!is.null(import)) pvt$augm_$add_quo_list(
    defs_to_quos_(import, name, rlang::enquos(...), rlang::caller_env()),
    append = append, .kind = "obs"
  )
  else if (!is.null(name)) pvt$augm_$add_quo(
    name = name, value = NULL, status = status,
    const = NA, discr = discr, prior = NA_character_,
    obs = obs, rlang::enquos(...),
    min = rlang::enquo(min), max = rlang::enquo(max),
    .assert = TRUE, .kind = "obs"
  )
  else pvt$augm_$add_quo_list(
    assert_bulk_defs_(rlang::enquos(...)),
    append = append, .kind = "obs"
  )
  pvt$rebuild_()
  new_obj
}

#' @title Add a status process
#' @description
#' Declare how a member's infected status is generated during simulation, as a
#' discrete distribution over the infected statuses. The probabilities are given
#' in declaration order (the \eqn{k}-th probability is that of the \eqn{k}-th
#' infected status); when named, the names become the status labels.
#' During evaluation the status is fixed observed data, so it contributes a
#' constant and the probabilities must be data or numeric constants.
#' @param hhmodel `r rd_model_in_`
#' @param ... Probabilities of the infected statuses, in order. When named, the
#' names label the corresponding statuses.
#' @param dist `character(1)`. Discrete distribution generating the status.
#' @return `r rd_model_out_`
#' @export
add_status <- function(hhmodel, ..., dist = "categorical") {
  stopifnothousepi(hhmodel)
  new_obj <- hhmodel$copy()
  pvt <- new_obj$.__enclos_env__$private
  weights <- rlang::enquos(...)
  if (length(weights) == 0) stop_h(
    "No status probability provided.",
    "Provide at least one infected-status probability."
  )
  # a single reserved node carries the status process (there is one status)
  pvt$augm_$add_quo(
    name = ".status", value = NULL, status = NA_integer_,
    const = NA, discr = FALSE, prior = NA_character_,
    obs = dist, weights,
    .assert = TRUE, .kind = "status"
  )
  # record the ordered status labels (the dictionary lives here, not in C++)
  pvt$augm_$set_status_levels(names(weights))
  pvt$rebuild_()
  new_obj
}

#' @title Add household-level observations
#' @description
#' Add household-level observations for all members.
#' Household-level observations must be reducible to a constant value.
#' `r rd_add_modes_`
#' @details
#' `r rd_add_details_`
#' @param hhmodel `r rd_model_in_`
#' @param name `r rd_var_name_`
#' @param value `r rd_aug_val_`
#' @param ... `r rd_obs_args_`
#' In mode 'multiple', the `name = value` definitions.
#' @param status `r rd_aug_status_`
#' @param discr `r rd_aug_discr_`
#' @param min `r rd_aug_min_`
#' @param max `r rd_aug_max_`
#' @param import `r rd_import_`
#' @param append `r rd_list_append_`
#' @return `r rd_model_out_`
#' @export
add_house_observ <- function(hhmodel, name = NULL, value = NULL, ...,
                             status = NA_integer_, discr = FALSE,
                             min = NULL, max = NULL,
                             import = NULL, append = TRUE) {
  stopifnothousepi(hhmodel)
  new_obj <- hhmodel$copy()
  pvt <- new_obj$.__enclos_env__$private
  if (!is.null(import)) pvt$augm_$add_quo_list(
    defs_to_quos_(import, name, rlang::enquos(...), rlang::caller_env()),
    append = append, .kind = "house"
  )
  else if (!is.null(name)) pvt$augm_$add_quo(
    name = name, value = rlang::enquo(value),
    status = status, discr = discr,
    min = rlang::enquo(min), max = rlang::enquo(max),
    .assert = TRUE, .kind = "house"
  )
  else pvt$augm_$add_quo_list(
    assert_bulk_defs_(rlang::enquos(...)),
    append = append, .kind = "house"
  )
  pvt$rebuild_()
  new_obj
}

############################################################################
############################################################################
############################################################################

#' @title Add a stratified group rate definition
#' @description `r rd_grr_desc_`
#' `r rd_add_modes_`
#' @details
#' `r rd_add_details_`
#' @param hhmodel `r rd_model_in_`
#' @param name `r rd_var_name_`
#' @param by `r rd_grr_by_`
#' @param reference `r rd_grr_refer_`
#' @param joint `r rd_grr_joint_`
#' @param strata `r rd_grr_strata_`
#' @param params `r rd_grr_params_`
#' @param prior `r rd_var_prior_`
#' @param ... `r rd_par_args_`
#' @param min `r rd_var_min_`
#' @param max `r rd_var_max_`
#' @param proposal `r rd_var_proposal_`
#' @param proposal_scale `r rd_proposal_scale_`
#' @param prefix `r rd_grr_prefix_`
#' @param import `r rd_import_`
#' @param append `r rd_list_append_`
#' @references
#' \insertRef{cauchemez_household_2009}{housepi}
#' @return `r rd_model_out_`
#' @export
add_group_rate <-  function(hhmodel,
                            name = NULL, by = NA_character_, ...,
                            reference = NULL, joint = FALSE,
                            strata = NULL, params = NULL,
                            prior = NA_character_,
                            min = NA_real_, max = NA_real_,
                            proposal = NA_character_,
                            proposal_scale = NA_real_,
                            prefix = NA_character_,
                            import = NULL, append = TRUE) {
  stopifnothousepi(hhmodel)
  new_obj <- hhmodel$copy()
  pvt <- new_obj$.__enclos_env__$private
  if (is.null(name) && isFALSE(append)) {
    # NOTE: dropping group rates might cause parameter drops
    pvt$drop_group_rates_()
  }
  if (!is.null(import)) pvt$group_rates_$add_quo_list(
    defs_to_quos_(import, name, rlang::enquos(...), rlang::caller_env()),
    append = append
  )
  else if (!is.null(name)) pvt$group_rates_$add_quo(
    name = name, by = by,
    reference = reference, joint = joint,
    strata = strata, params = params,
    prior = prior, rlang::enquos(...),
    min = rlang::enquo(min), max = rlang::enquo(max),
    proposal = proposal, proposal_scale = proposal_scale,
    prefix = prefix
  )
  else pvt$group_rates_$add_quo_list(
    assert_bulk_defs_(rlang::enquos(...)),
    append = append
  )
  pvt$rebuild_()
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
#' @param prior `r rd_var_prior_`
#' @param ... `r rd_par_args_`
#' @param min `r rd_var_min_`
#' @param max `r rd_var_max_`
#' @param proposal `r rd_var_proposal_`
#' @param proposal_scale `r rd_proposal_scale_`
#' @param prefix `r rd_ctr_prefix_`
#' @return `r rd_model_out_`
#' @export
add_contact_rate <- function(hhmodel, by = NA_character_, ...,
                             reference = NULL,
                             joint = FALSE, symmetric = TRUE,
                             strata = NULL, params = NULL,
                             prior = NA_character_,
                             min = NA_real_, max = NA_real_,
                             proposal = NA_character_,
                             proposal_scale = NA_real_,
                             prefix = NA_character_) {
  # FIXME: missing an import for contact-rate
  stopifnothousepi(hhmodel)
  new_obj <- hhmodel$copy()
  pvt <- new_obj$.__enclos_env__$private
  pvt$contact_rates_$add_quo(
    by = by, reference = reference,
    joint = joint, symmetric = symmetric,
    strata = strata, params = params,
    prior = prior, rlang::enquos(...),
    min = rlang::enquo(min), max = rlang::enquo(max),
    proposal = proposal, proposal_scale = proposal_scale,
    prefix = prefix
  )
  pvt$rebuild_()
  new_obj
}

############################################################################
############################################################################
############################################################################

#' @title Add an infectivity profile definition
#' @description
#' Add definition of infectivity profile.
#' `r rd_add_modes_`
#' @details
#' `r rd_add_details_`
#' @param hhmodel `r rd_model_in_`
#' @param func `r rd_var_name_`
#' @param status `r rd_aug_status_`
#' @param discr `r rd_inf_discr_`
#' @param ... `r rd_inf_args_`
#' In mode 'multiple', the `func` definitions.
#' @param max `r rd_aug_max_`
#' @param import `r rd_import_`
#' @param append `r rd_list_append_`
#' @return `r rd_model_out_`
#' @export
add_infec_profile <- function(hhmodel,
                              func = NULL, ...,
                              status = NA_integer_, discr = FALSE,
                              max = NULL,
                              import = NULL, append = TRUE) {
  stopifnothousepi(hhmodel)
  new_obj <- hhmodel$copy()
  pvt <- new_obj$.__enclos_env__$private
  if (!is.null(import)) pvt$infec_profile_$add_quo_list(
    defs_to_quos_(import, func, rlang::enquos(...), rlang::caller_env()),
    append = append
  )
  else if (!is.null(func)) pvt$infec_profile_$add_quo(
    func = func, status = status, discr = discr,
    rlang::enquos(...), max = rlang::enquo(max),
    .assert = TRUE
  )
  else pvt$infec_profile_$add_quo_list(
    assert_bulk_defs_(rlang::enquos(...)),
    append = append
  )
  pvt$rebuild_()
  new_obj
}

############################################################################
############################################################################
############################################################################

#' @title Evaluate the log-likelihood of a housepi model at current values
#' @description
#' Compute the complete-data log-likelihood of a `HousepiModel`:
#' transmission plus observation.
#' For models with augmented variables
#' this is conditioned on current augmented values.
#' @param hhmodel `r rd_model_in_`
#' @param house_id `r rd_house_id_`
#' If `NA` (default) it computes for the whole data.
#' @return `numeric(1)`. The complete-data log-likelihood.
#' `r rd_natural_log_`
#' @examples
#' flu_1918_frederick |>
#'   new_housepi_model("final_size") |>
#'   add_parameter_template() |>
#'   evaluate_log_lik()
#' @export
evaluate_log_lik <- function(hhmodel, house_id = NA_character_) {
  stopifnothousepi(hhmodel)
  pvt <- hhmodel$.__enclos_env__$private
  pvt$build_full_()
  .engine$evaluate_log_lik(pvt$generator_$mcmc_ptr, house_id)
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
#' @param ... `r rd_dots_config_`
#' @return An object of class `HousepiPosteriorDraws`.
#' @examples
#' flu_1918_frederick |>
#'   new_housepi_model("final_size") |>
#'   add_parameter_template() |>
#'   run_inference(steps = 50, burn = 10, thin = 2,
#'                 chains = 2, progress = FALSE) |>
#'   get_parameter_draws() |>
#'   posterior::summarize_draws()
#' @export
run_inference <- function(hhmodel,
                          steps = 1000L, burn = 0L,
                          thin = 5L, adapt = 500L,
                          chains = 4L, threads = 1L,
                          id_draws = NULL,
                          progress = TRUE,
                          ...) {
  stopifnothousepi(hhmodel)
  assert_advanced_setup_(list(...), "inv_temp")
  pvt <- hhmodel$.__enclos_env__$private
  pvt$assert_posterior_memory_(
    steps = steps, burn = burn, thin = thin,
    chains = chains, id_draws = id_draws
  )
  pvt$build_full_(
    steps = steps, burn = burn, thin = thin, adapt = adapt,
    chains = chains, threads = threads, id_draws = id_draws,
    progress = progress, ...
  )
  .engine$run_inference(pvt$generator_$mcmc_ptr)
  HousepiPosteriorDraws$new(pvt$generator_$raw_name,
                            pvt$generator_$mcmc_ptr)
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
#' @param ... `r rd_dots_config_`
#' @return `data.frame` or `list` (when `obs_process_long` is provided)
#' @examples
#' flu_1918_frederick |>
#'   new_housepi_model("final_size") |>
#'   add_parameter_template() |>
#'   run_simulation() |>
#'   head()
#' @export
run_simulation <- function(hhmodel,
                           obs_process = NULL,
                           obs_process_long = NULL,
                           param_table = NULL,
                           nrep = NA_integer_,
                           ...) {
  stopifnothousepi(hhmodel)
  # simulation specific build
  assert_advanced_setup_(list(...), "time_step")
  pvt <- hhmodel$.__enclos_env__$private
  nrep <- pvt$helper_simulation_build_(
    param_table = param_table, nrep = nrep, ...
  )

  # combine observation process to measures and refactor
  obs_func <- pvt$helper_obs_process_(obs_process)
  # run simulations
  output <- .engine$run_simulation(
    pvt$generator_$mcmc_ptr,
    obs_func, NULL,
    obs_process_long
  )

  # drop empty wide/long
  if (length(output$wide) == 0) output$wide <- NULL
  if (length(output$long) == 0) output$long <- NULL
  if (length(output) == 1) {
    output <- output[[1]]
    if (nrep == 1) output <- output[[1]]
  } else if (nrep == 1) {
    output$wide <- output$wide[[1]]
    output$long <- output$long[[1]]
  }

  output
}
