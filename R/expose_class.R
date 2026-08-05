## housepi - Epidemiological Analysis of Household Outbreaks.
##
## Copyright © 2024-2026 Institut Pasteur (Paris), CNRS and INSERM.
## Distributed under the terms of the GNU General Public License (GPLv3).
## See the LICENSE file for details.
##

##########################################################################
############################# DOCUMENTATION ##############################
##########################################################################

## Model surface ##

rd_model_in_ <- "`HousepiModel` object."
rd_model_out_ <- paste(
  rd_model_in_,
  "Deep clone of input with new definitions."
)
rd_add_wide_ <- "`data.frame` Household data used for analysis (wide format).
Optional at construction; if omitted, attach it later with
[add_data_wide()][housepi::add_data_wide]."
rd_house_id_ <- "`character(1)` ID of the household of interest."
rd_mdl_choice_ <- "`character(1)` Choice for the model.
One between:
\\describe{
  \\item{`<template>`}{a built-in template name;}
  \\item{`path/to/source`}{folder with implementation file `\"model.cpp\"`;}
  \\item{`path/to/compiled`}{folder with a previously compiled model;}
  \\item{`\"test\"`}{template name for test model;}
  \\item{`NA`}{fall back to `\"test\"` (default).}
}"
rd_mdl_astmpl_ <- "`logical(1)` Store dynamic library as a local template.
If `FALSE` recompile the model at every new session (default)."
rd_mdl_force_ <- "`logical(1)` Compile even if cached model already exists."
rd_mdl_quiet_ <- "`logical(1)` If `TRUE` hide compilation commands."
rd_mdl_escaped_ <- "`character(1)` Label for status code `0` (individuals not
among the declared infected statuses) when the `status` column is
character/factor."
rd_mdl_flags_ <- "Named arguments. Additional compilation options:
\\describe{
  \\item{`discr`}{`logical(1)` If `TRUE` compile to discrete time.}
  \\item{`flags`}{`character()` Compilation flags.}
  \\item{`libs`}{`character()` Links to libraries.}
  \\item{`paths`}{`character()` Include paths.}
}"

## Definition vocabulary ##

rd_var_name_ <- "`character(1)` Name used by internal calls."
rd_var_const_ <- "`logical(1)` Treat definition as constant.
If `TRUE`, additional Monte Carlo details are ignored."
rd_var_discr_ <- function(what) {
  sprintf("`logical(1)` Treat definition as discrete.
If applicable, the %s is discretized.", what)
}
rd_var_prior_ <- "`character(1)` Prior distribution."
rd_var_args_ <- function(what, proposal = FALSE) {
  paste0(sprintf("named arguments for the %s distribution.", what),
         if (proposal) rd_proposal_opts_)
}
rd_var_min_ <- "`numeric(1)` Truncate to a minimum value.
Relative to the location argument `loc`, when one is given."
rd_var_max_ <- "`numeric(1)` Truncate to a maximum value.
Relative to the location argument `loc`, when one is given."
rd_var_func_ <- "One between:
- `character(1)` Name of augmented function;
- rlang expression."
rd_var_funargs_ <- "named `list` Arguments of the augmented function."
rd_var_proposal_ <- "`character(1)` Proposal used in the Monte Carlo.
Random walks perturb the current value:
\\describe{
  \\item{`\"rw\"`}{normal increments (default);}
  \\item{`\"rw_log\"`}{lognormal log-increments (positive values only);}
  \\item{`\"rw_unif\"`}{uniform (bounded) increments.}
}
Independence samplers draw a fresh value from a distribution:
\\describe{
  \\item{`\"forward\"`}{drawn from the variable's own defining distribution
    (its prior for a parameter, its model for an augmented variable);}
  \\item{a distribution name}{drawn from that distribution, e.g.
    `proposal = \"normal\"` with `proposal_mean`/`proposal_sd`.}
}
`NA` defaults to `\"rw\"`.
"
rd_proposal_opts_ <- "
Advanced Monte Carlo options may also be passed through `...`, each prefixed
`proposal_`. For a random walk, `proposal_scale` sets the step scale: if
omitted it is adapted during inference when `adapt > 0`; if a positive numeric
it is fixed and never adapted. For independence sampling from a named
distribution, the distribution parameters are supplied the same way (e.g.
`proposal_mean`, `proposal_sd`)."
# FIXME:
rd_types_ <- "
One between:
\\describe{
  \\item{augmented}{`@A`;}
  \\item{data}{`@D`;}
  \\item{group risk}{`@G`;}
  \\item{parameter}{`@P`;}
  \\item{scalar}{`@N`.}
}"
rd_aug_status_ <- "`integer(1)` Individual positive status.
Used to specify status-specific definitions.
If `NA`, associate definition to all infectious individuals (default)."
rd_aug_val_ <- paste("`numeric(1)` | `character(1)`. Initial value.",
                     rd_types_)
rd_aug_min_ <- paste(rd_var_min_, rd_types_)
rd_aug_max_ <- paste(rd_var_max_, rd_types_)

## Definition modes ##

rd_add_modes_ <- "There are three mutually exclusive definition modes:
\\describe{
  \\item{single}{
  Supply a `name` to add a single definition using the explicit arguments;
  }
  \\item{multiple}{
  Pass pairs `name = list()` (or `name = value`) for multiple definitions;
  }
  \\item{import}{
  Pass a pre-built named list of definitions through `import`.
  }
}"
rd_add_details_ <- "A definition name colliding with reserved argument
is deprecated.
It can be supplied only via single mode definition."
rd_import_ <- "Named `list` of pre-built definitions to import.
Cannot be combined with single and/or multiple definitions."
rd_tmpl_import_ <- "Named `list` of pre-built definitions to import,
or `\"template\"` for the definitions shipped with the model.
Cannot be combined with single and/or multiple definitions."
rd_list_append_ <- "`logical(1)` Append to current definitions.
If `FALSE` drop all pre-existing definitions."

## Risks ##

rd_risk_desc_ <- function(kind) {
  sprintf("Add definitions of a %s risk by given data covariate(s).
Risks are relative to a reference class that can be specified:
each parameter multiplies the baseline risk of its reference class.
If no MCMC options are provided, the prior is assumed to be
\\eqn{\\mathrm{lognormal}[\\mu=0, \\sigma=3]} with support \\eqn{[0, 10]}.
Bypass automatic generation by providing a custom pre-built logic table.",
          kind)
}
rd_risk_prefix_ <- "`character(1)` Prefix to names of associated parameters."
rd_grr_prefix_ <- paste(
  rd_risk_prefix_,
  "If `NA` uses the first 3 letters of risk name (default)."
)
rd_risk_joint_ <- "`logical(1)` Treat covariates as joint.
If `FALSE`, covariates are independent (default)."
rd_risk_strata_ <-
  "`data.frame` Individual covariates and relative `param_name`.
Each row defines one stratum. If `NA`, then match any value.
Use `.reference` in `param_name` to designate the reference group.
Mutually exclusive with `by`."
rd_risk_params_ <- "named `list` Per-parameter MCMC setups as in `param_name`.
Overrides shared prior/args for parameters specified in the list.
Ignored if automatic definition."

## Run-time knobs.

rd_smpl_steps_ <- "`integer(1)` Number of MCMC iterations per chain."
rd_smpl_burn_ <- "`integer(1)` Number of initial iterations within `steps` to
discard before storing draws."
rd_smpl_thin_ <- "`integer(1)` Keep one iteration every `thin` after burn-in."
rd_smpl_adapt_ <- "`integer(1)` Number of iterations to adapt proposal scales.
They are run before `steps` and are not stored as draws.
If no proposals are adaptable and `burn = 0`, repurposed as burn-in instead."
rd_smpl_chain_ <- "`integer(1)` Number of independent MCMC chains to run."
rd_thread_ <- "`integer(1)` Number of threads for parallel household runs."
rd_progress_ <- "`integer(1)` | `logical(1)` Control progressbar visualization.
One between:
- If `integer`, prints progress every `progress` seconds;
- If `TRUE` prints progress every second (default);
- If `FALSE` or `0`, suppresses progress display."
rd_dots_config_ <- "Advanced configuration."
rd_param_table_ <- "One between
- `data.frame`;
- `posterior::draws`;
- [HousepiPosteriorDraws][housepi::HousepiPosteriorDraws].
"
rd_simul_observ_ <- "`function` Observation process.
Function expects as input the `data.frame` generated by simulation.
If `NULL`, observation process is ignored."
rd_simul_obs_long_ <- "`function` Optional function.
Expects processed `data_wide` and returns a `data_long` data.frame.
When provided, each simulation returns `list(wide = ..., long = ...)`."
rd_natural_log_ <- "Logarithms are natural (base \\eqn{e})."

############################################################################
############################################################################
############################################################################

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

#' @title Attach household wide-format data to a housepi model
#' @description
#' Make an internal copy of the provided household data in wide format,
#' (re)building any definitions that depend on it.
#' @param hhmodel `r rd_model_in_`
#' @param data `r rd_add_wide_`
#' @return `r rd_model_out_`
#' @export
add_data_wide <- function(hhmodel, data) {
  stopifnothousepi(hhmodel)
  if (is.null(data)) stop_h("Missing argument {.arg data}.")
  new_obj <- hhmodel$copy()
  pvt <- new_obj$.__enclos_env__$private
  # drop stale risk params (read pre-change tables) before re-importing data
  pvt$drop_stale_risk_params_()
  pvt$set_data_wide_(assert_household_data(data, NA_character_))
  pvt$rebuild_()
  new_obj
}

#' @title Attach individual-level long-format data to a housepi model
#' @description
#' Make an internal copy of the provided data in long format.
#' @param hhmodel `r rd_model_in_`
#' @param data_long `data.frame` Additional individual long format data.
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
#' @return `data.frame`, one row per metric:
#' \describe{
#'   \item{level}{`individual` or `household`}
#'   \item{metric}{quantity counted in the row}
#'   \item{count}{absolute count}
#'   \item{percent}{percentage of the relevant total}
#' }
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
#' @return `data.frame`, one row per metric:
#' \describe{
#'   \item{metric}{quantity reported in the row}
#'   \item{value}{household size for that metric}
#'   \item{count}{absolute count}
#'   \item{percent}{percentage of individuals}
#' }
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
#' observations, group risks, contact risks, infectivity profile).
#' @details
#' `get_data_wide()` returns the individual-level data the model reads: a column
#' a deterministic observation stands in for carries the derived value, not the
#' imported one. A column an observation produces without standing in for data
#' is not part of the dataset, and is reported by [get_observations()].
#' @param hhmodel `r rd_model_in_`
#' @return `data.frame`
#' @examples
#' hhmodel <- new_housepi_model(flu_1918_frederick, "final_size") |>
#'   add_parameter(import = "template")
#' get_data_wide(hhmodel) |>
#'   head()
#' get_parameters(hhmodel)
#' @rdname get_definitions
#' @export
get_data_wide <- function(hhmodel) {
  stopifnothousepi(hhmodel)
  # UPDATE: drop active getter fields and move logic here
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
get_observations <- function(hhmodel) {
  stopifnothousepi(hhmodel)
  hhmodel$observations
}

#' @rdname get_definitions
#' @export
get_group_risks <- function(hhmodel) {
  stopifnothousepi(hhmodel)
  # NOTE: internal build requires deep clone to avoid in-place
  hhmodel$copy()$group_risks
}

#' @rdname get_definitions
#' @export
get_contact_risks <- function(hhmodel) {
  stopifnothousepi(hhmodel)
  # NOTE: internal build requires deep clone to avoid in-place
  hhmodel$copy()$contact_risks
}

#' @rdname get_definitions
#' @export
get_infectivity_profile <- function(hhmodel) {
  stopifnothousepi(hhmodel)
  hhmodel$infectivity_profile
}

############################################################################
############################################################################
############################################################################

definition_helper_ <- c(
  "Supply a {.arg name} for a single definition.",
  "Provide multiple definitions as {.code name = list()} pairs.",
  "Otherwise import an external list using {.arg import}."
)

# Extra hint for the families shipping a definition template.
template_helper_ <- paste(
  "Load the definitions shipped with the model with",
  '{.code import = "template"}.'
)

# Guard the inline (`...`) bulk form.
assert_bulk_defs_ <- function(dots, has_template = FALSE) {
  if (length(dots) > 0L) return(dots)
  sms <- c("No definitions provided.", definition_helper_)
  if (isTRUE(has_template)) sms <- c(sms, template_helper_)
  stop_h(sms)
}

# Reject `import` combined with single-mode or inline `...` pairs.
assert_no_def_mixing_ <- function(key, dots) {
  if (!is.null(key) || length(dots) > 0L) stop_h(
    "Mixing imported definitions with single and multiple definitions",
    definition_helper_
  )
}

# Build the quosure list for the explicit `import` bulk form.
defs_to_quos_ <- function(import, key, dots, env) {
  assert_no_def_mixing_(key, dots)
  convert_list_to_quosures(import, env)
}

# Recognise the `import` sentinel selecting the definitions shipped with the
# model. Any other character value names a source that does not exist.
is_template_import_ <- function(import) {
  if (!is.character(import)) return(FALSE)
  if (!identical(import, "template")) stop_h(
    "Unknown import source {.val {import}}.",
    template_helper_
  )
  TRUE
}

# Load default templated definitions shipped with the package.
# Reject the sentinel combined with single-mode or inline `...` pairs.
add_template_defs_ <- function(pvt, family, append, key, dots) {
  assert_no_def_mixing_(key, dots)
  holder <- switch(family,
                   parameter = pvt$param_,
                   augmented = pvt$augm_,
                   `infectivity profile` = pvt$infec_profile_)
  if (pvt$generator_$is_test_model()) {
    warn_h("No template {family} definitions available for test model.")
  } else {
    defs <- get_template_def_list(pvt$generator_$source_dir, family)
    holder$add_quo_list(convert_list_to_quosures(defs), append = isTRUE(append))
  }
}


#' @title Add or overwrite parameter definitions
#' @description
#' Add parameter(s) to a housepi model.
#' If a parameter with the same name exists, the former gets overwritten.
#' `r rd_add_modes_`
#' @details
#' `r rd_add_details_`
#' @param hhmodel `r rd_model_in_`
#' @param name `r rd_var_name_`
#' @param value `numeric(1)` Initial value.
#' Required only if constant or random initialization has been disabled.
#' @param const `r rd_var_const_`
#' @param discr `r rd_var_discr_("prior")`
#' @param prior `r rd_var_prior_`
#' @param ... `r rd_var_args_("prior", proposal = TRUE)`
#' In mode 'multiple', the `name = value` definitions.
#' @param min `r rd_var_min_`
#' @param func `r rd_var_func_`
#' @param func_args `r rd_var_funargs_`
#' @param max `r rd_var_max_`
#' @param proposal `r rd_var_proposal_`
#' @param import `r rd_tmpl_import_`
#' @param append `r rd_list_append_`
#' @return `r rd_model_out_`
#' @examples
#' # the single, inline-bulk and list-bulk forms compose freely in a pipe
#' hhmodel <- flu_1918_frederick |>
#'   new_housepi_model("final_size") |>
#'   # single: alpha, with a uniform prior and a lognormal random-walk proposal
#'   add_parameter("alpha", prior = "uniform", proposal = "rw_log") |>
#'   # single: beta, fixed to a constant
#'   add_parameter("beta", value = 0.3, const = TRUE) |>
#'   # inline bulk: several definitions as `name = value` pairs
#'   add_parameter(gamma = 1, hhsize_ref = 4)
#'
#' # list bulk: load a pre-built named list of definitions through `import`
#' my_import <- list(gamma = 1, hhsize_ref = 4)
#' hhmodel <- flu_1918_frederick |>
#'   new_housepi_model("final_size") |>
#'   add_parameter(import = my_import)
#'
#' # `"template"` loads the definitions shipped with the model
#' hhmodel <- flu_1918_frederick |>
#'   new_housepi_model("final_size") |>
#'   add_parameter(import = "template")
#' get_parameters(hhmodel)
#' @export
add_parameter <- function(hhmodel,
                          name = NULL, value = NULL, ...,
                          const = NA, discr = FALSE,
                          prior = NA_character_,
                          func = NULL, func_args = NULL,
                          min = NULL, max = NULL,
                          proposal = NA_character_,
                          import = NULL, append = TRUE) {
  stopifnothousepi(hhmodel)
  new_obj <- hhmodel$copy()
  pvt <- new_obj$.__enclos_env__$private

  #> template mode
  if (is_template_import_(import)) add_template_defs_(
    pvt, "parameter", append, name, rlang::enquos(...)
  )
  #> import mode
  else if (!is.null(import)) pvt$param_$add_quo_list(
    defs_to_quos_(import, name, rlang::enquos(...), rlang::caller_env()),
    append = append
  )
  #> single mode
  else if (!is.null(name)) pvt$param_$add_quo(
    name = name, value = rlang::enquo(value),
    const = const, discr = discr,
    prior = prior, rlang::enquos(...),
    func = rlang::enquo(func), func_args = rlang::enquos(func_args),
    min = rlang::enquo(min), max = rlang::enquo(max),
    proposal = proposal,
    .assert = TRUE
  )
  #> multiple mode
  else pvt$param_$add_quo_list(
    assert_bulk_defs_(rlang::enquos(...), TRUE),
    append = append
  )
  pvt$rebuild_()
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
#' @param discr `r rd_var_discr_("model")`
#' @param func `r rd_var_func_`
#' @param func_args `r rd_var_funargs_`
#' @param model `character(1)` Augmented model distribution.
#' @param ... `r rd_var_args_("model", proposal = TRUE)`
#' In mode 'multiple', the `name = value` definitions.
#' @param min `r rd_aug_min_`
#' @param max `r rd_aug_max_`
#' @param proposal `r rd_var_proposal_`
#' @param import `r rd_tmpl_import_`
#' @param append `r rd_list_append_`
#' @return `r rd_model_out_`
#' @examples
#' hhmodel <- dummy_household |>
#'   new_housepi_model("test") |>
#'   # define an incubation period drawn from a gamma distribution
#'   add_augmented("incub_period",
#'                 model = "gamma", shape = 2, scale = 1, max = 14)
#' @export
add_augmented <- function(hhmodel,
                          name = NULL, value = NULL, ...,
                          status = NA_integer_, const = NA, discr = FALSE,
                          func = NULL, func_args = NULL,
                          model = NA_character_,
                          min = NULL, max = NULL,
                          proposal = NA_character_,
                          import = NULL, append = TRUE) {
  stopifnothousepi(hhmodel)
  new_obj <- hhmodel$copy()
  pvt <- new_obj$.__enclos_env__$private

  #> template mode
  if (is_template_import_(import)) add_template_defs_(
    pvt, "augmented", append, name, rlang::enquos(...)
  )
  #> import mode
  else if (!is.null(import)) pvt$augm_$add_quo_list(
    defs_to_quos_(import, name, rlang::enquos(...), rlang::caller_env()),
    append = append, status = status
  )
  #> single mode
  else if (!is.null(name)) pvt$augm_$add_quo(
    name = name, value = rlang::enquo(value),
    status = status, const = const, discr = discr,
    func = rlang::enquo(func), func_args = rlang::enquos(func_args),
    model = model, rlang::enquos(...),
    min = rlang::enquo(min), max = rlang::enquo(max),
    proposal = proposal,
    .assert = TRUE
  )
  #> multiple mode
  else pvt$augm_$add_quo_list(
    assert_bulk_defs_(rlang::enquos(...), TRUE),
    append = append, status = status
  )
  pvt$rebuild_()
  new_obj
}

#' @title Add an observation
#' @description
#' Add observation associated to data \eqn{x}, at the level set by `level`.
#'
#' At member level, the observation is attached to every infected individual.
#' If `model` is provided, adds to household likelihood during evaluations,
#' and \eqn{x} is sampled from it during simulations.
#'
#' At household level, the observation is a reduction across all members,
#' it's shared by the whole household, and must reduce to a constant value.
#'
#' An observation carrying no `model` produces its own data column, at either
#' level: it is evaluated for every member, and eventually overwrites data
#' column of the same name
#' `r rd_add_modes_`
#' @details
#' `r rd_add_details_`
#' @param hhmodel `r rd_model_in_`
#' @param name `r rd_var_name_`
#' @param model `character(1)` Choice for the observation model distribution.
#' Reserved to member-level observations.
#' @param ... `r rd_var_args_("observation")`
#' In mode 'multiple', the `name = value` definitions.
#' @param level `character(1)` Level the datum is observed at.
#' One between `"member"`, one value per positive status individual, and
#' `"house"`, one value per household reduced over all its members.
#' In mode 'multiple', applies to every entry not setting its own.
#' @param value `r rd_aug_val_`
#' @param status `r rd_aug_status_`
#' @param discr `r rd_var_discr_("model")`
#' @param min `r rd_aug_min_`
#' @param max `r rd_aug_max_`
#' @param import `r rd_import_`
#' @param append `r rd_list_append_`
#' @return `r rd_model_out_`
#' @export
add_observation <- function(hhmodel, name = NULL,
                            model = NA_character_, ...,
                            level = c("member", "house"), value = NULL,
                            status = NA_integer_, discr = FALSE,
                            min = NULL, max = NULL,
                            import = NULL, append = TRUE) {
  stopifnothousepi(hhmodel)
  # bulk default: a per-definition 'level =' overrides it (see add_quo)
  level <- match.arg(level)
  new_obj <- hhmodel$copy()
  pvt <- new_obj$.__enclos_env__$private

  #> import mode
  if (!is.null(import)) pvt$augm_$add_quo_list(
    defs_to_quos_(import, name, rlang::enquos(...), rlang::caller_env()),
    append = append, .kind = level
  )
  #> single mode
  else if (!is.null(name)) pvt$augm_$add_quo(
    name = name, value = rlang::enquo(value), status = status,
    const = NA, discr = discr,
    model = model, rlang::enquos(...),
    min = rlang::enquo(min), max = rlang::enquo(max),
    .assert = TRUE, .kind = level
  )
  #> multiple mode
  else pvt$augm_$add_quo_list(
    assert_bulk_defs_(rlang::enquos(...)),
    append = append, .kind = level
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
#' @param model `character(1)`. Discrete distribution generating the status.
#' @return `r rd_model_out_`
#' @export
add_status <- function(hhmodel, ..., model = "categorical") {
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
    const = NA, discr = FALSE,
    model = model, weights,
    .assert = TRUE, .kind = "status"
  )
  # record the ordered status labels (the dictionary lives here, not in C++)
  pvt$augm_$set_status_levels(names(weights))
  pvt$rebuild_()
  new_obj
}

############################################################################
############################################################################
############################################################################

#' @title Add a stratified group risk definition
#' @description `r rd_risk_desc_("group")`
#' `r rd_add_modes_`
#' @details
#' `r rd_add_details_`
#' @param hhmodel `r rd_model_in_`
#' @param name `r rd_var_name_`
#' @param by `list` Covariates associated to group risk as data headers.
#' Corresponding data columns have to be categorical.
#' @param reference named `list` Covariates associated to reference group.
#' Pass `".none"` to disable reference, giving every risk its own parameter.
#' Default reference group is the one with the largest sample size.
#' @param joint `r rd_risk_joint_`
#' @param strata `r rd_risk_strata_`
#' @param params `r rd_risk_params_`
#' @param prior `r rd_var_prior_`
#' @param ... `r rd_var_args_("prior", proposal = TRUE)`
#' @param min `r rd_var_min_`
#' @param max `r rd_var_max_`
#' @param proposal `r rd_var_proposal_`
#' @param prefix `r rd_grr_prefix_`
#' @param import `r rd_import_`
#' @param append `r rd_list_append_`
#' @references
#' \insertRef{cauchemez_household_2009}{housepi}
#' @return `r rd_model_out_`
#' @export
add_group_risk <-  function(hhmodel,
                            name = NULL, by = NA_character_, ...,
                            reference = NULL, joint = FALSE,
                            strata = NULL, params = NULL,
                            prior = NA_character_,
                            min = NA_real_, max = NA_real_,
                            proposal = NA_character_,
                            prefix = NA_character_,
                            import = NULL, append = TRUE) {
  stopifnothousepi(hhmodel)
  new_obj <- hhmodel$copy()
  pvt <- new_obj$.__enclos_env__$private

  # drop stale parameters and re-derive in `rebuild_()`
  pvt$drop_stale_risk_params_()
  if (is.null(name) && isFALSE(append)) pvt$drop_group_risks_()

  #> import mode
  if (!is.null(import)) pvt$group_risks_$add_quo_list(
    defs_to_quos_(import, name, rlang::enquos(...), rlang::caller_env()),
    append = append
  )
  #> single mode
  else if (!is.null(name)) pvt$group_risks_$add_quo(
    name = name, by = by,
    reference = reference, joint = joint,
    strata = strata, params = params,
    prior = prior, rlang::enquos(...),
    min = rlang::enquo(min), max = rlang::enquo(max),
    proposal = proposal, prefix = prefix
  )
  #> multiple mode
  else pvt$group_risks_$add_quo_list(
    assert_bulk_defs_(rlang::enquos(...)),
    append = append
  )
  pvt$rebuild_()
  new_obj
}

############################################################################
############################################################################
############################################################################

#' @title Add a pairwise contact risk definition
#' @description `r rd_risk_desc_("contact")`
#' @param hhmodel `r rd_model_in_`
#' @param by `character()` Covariates for from/to matching as data headers.
#' Corresponding data columns have to be categorical.
#' @param reference named `list` Covariates associated to reference pair.
#' Names can use `from_`/`to_` prefixes matching `by`,
#' or un-prefixed names (applied to both from and to).
#' @param joint `r rd_risk_joint_`
#' @param symmetric `logical(1)` Treat contact as symmetric.
#' If `TRUE`, contact pairs infector-infectee are treated as equivalent
#' (default).
#' @param strata `r rd_risk_strata_`
#' @param params `r rd_risk_params_`
#' @param prior `r rd_var_prior_`
#' @param ... `r rd_var_args_("prior", proposal = TRUE)`
#' @param min `r rd_var_min_`
#' @param max `r rd_var_max_`
#' @param proposal `r rd_var_proposal_`
#' @param prefix `r rd_risk_prefix_`
#' @param import Named `list` of pre-built definition to import.
#' Cannot be combined with the single definition.
#' @return `r rd_model_out_`
#' @export
add_contact_risk <- function(hhmodel, by = NA_character_, ...,
                             reference = NULL,
                             joint = FALSE, symmetric = TRUE,
                             strata = NULL, params = NULL,
                             prior = NA_character_,
                             min = NA_real_, max = NA_real_,
                             proposal = NA_character_,
                             prefix = NA_character_,
                             import = NULL) {
  stopifnothousepi(hhmodel)
  new_obj <- hhmodel$copy()
  pvt <- new_obj$.__enclos_env__$private
  pvt$assert_contact_risks_supported_()
  # the definition is replaced: its parameters are re-derived by `rebuild_()`
  pvt$drop_stale_risk_params_()

  # contact risk holds a single definition:
  # import a named arg list, or the explicit arguments
  if (is.null(import)) import <- list(
    by = by, reference = reference, joint = joint,
    symmetric = symmetric, strata = strata, params = params,
    prior = prior, rlang::enquos(...),
    min = rlang::enquo(min), max = rlang::enquo(max),
    proposal = proposal, prefix = prefix
  )
  do.call(pvt$contact_risks_$add_quo, import)
  pvt$rebuild_()
  new_obj
}

############################################################################
############################################################################
############################################################################

#' @title Add an infectivity profile definition
#' @description
#' Add definitions for the infectivity profiles.
#' `r rd_add_modes_`
#' @details
#' `r rd_add_details_`
#' @param hhmodel `r rd_model_in_`
#' @param func `r rd_var_name_`
#' @param status `r rd_aug_status_`
#' @param discr `r rd_var_discr_("profile")`
#' @param ... `r rd_var_args_("profile")`
#' In mode 'multiple', the `func` definitions.
#' @param min `r rd_aug_min_`
#' @param max `r rd_aug_max_`
#' @param import `r rd_tmpl_import_`
#' @param append `r rd_list_append_`
#' @return `r rd_model_out_`
#' @examples
#' # `"template"` loads the profile shipped with the model
#' hhmodel <- sarscov2_2021_icofs |>
#'   new_housepi_model("case_ascertained") |>
#'   add_infectivity_profile(import = "template")
#' get_infectivity_profile(hhmodel)
#' @export
add_infectivity_profile <- function(hhmodel,
                                    func = NULL, ...,
                                    status = NA_integer_, discr = FALSE,
                                    min = NULL, max = NULL,
                                    import = NULL, append = TRUE) {
  stopifnothousepi(hhmodel)
  new_obj <- hhmodel$copy()
  pvt <- new_obj$.__enclos_env__$private
  pvt$assert_infec_profile_supported_()
  # impose discrete-time if it's a model choice
  pvt$infec_profile_$discr_model <- pvt$generator_$is_discrete_time

  #> template mode
  if (is_template_import_(import)) add_template_defs_(
    pvt, "infectivity profile", append, func, rlang::enquos(...)
  )
  #> import mode
  else if (!is.null(import)) pvt$infec_profile_$add_quo_list(
    defs_to_quos_(import, func, rlang::enquos(...), rlang::caller_env()),
    append = append
  )
  #> single mode
  else if (!is.null(func)) pvt$infec_profile_$add_quo(
    func = func, status = status, discr = discr,
    rlang::enquos(...),
    min = rlang::enquo(min), max = rlang::enquo(max),
    .assert = TRUE
  )
  #> multiple mode
  else pvt$infec_profile_$add_quo_list(
    assert_bulk_defs_(rlang::enquos(...), has_template = TRUE),
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
#'   add_parameter(import = "template") |>
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
#' @param id_draws Infected members whose public augmented draws are stored.
#' Either a `character()` of individual identifiers,
#' or a `list` with entries `id` (individual identifiers)
#' and/or `house_id` (whole household infected members are stored).
#' When unset, one random eligible household is followed (default).
#' Pass `".none"` to store nothing.
#' @param progress `r rd_progress_`
#' @param ... `r rd_dots_config_`
#' @return An object of class `HousepiPosteriorDraws`.
#' @examples
#' flu_1918_frederick |>
#'   new_housepi_model("final_size") |>
#'   add_parameter(import = "template") |>
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
  id_draws <- pvt$resolve_id_draws_(id_draws)
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
#' @param nrep `integer(1)` Number of replicates to generate.
#' If `param_table` is `NULL`, all replicates share the same parameters.
#' Otherwise, if `nrep` is smaller than number of `param_table` rows,
#' the table is sub-sampled to `nrep` evenly-spaced rows.
#' If `NA` produces a single repetition (default).
#' @param ... `r rd_dots_config_`
#' @return `data.frame` or `list` (when `obs_process_long` is provided)
#' @examples
#' flu_1918_frederick |>
#'   new_housepi_model("final_size") |>
#'   add_parameter(import = "template") |>
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
