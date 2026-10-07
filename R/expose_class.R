## housepi - Epidemiological Analysis of Household Outbreaks.
##
## Copyright © 2024-2026 Institut Pasteur (Paris), CNRS and INSERM.
## Distributed under the terms of the GNU General Public License (GPLv3).
## See the LICENSE file for details.
##

##########################################################################
############################# DOCUMENTATION ##############################
##########################################################################

# nocov start

## Model surface ##

rd_model_in_ <- "`HousepiModel` object."
rd_model_out_ <- paste(
  rd_model_in_,
  "Deep clone of input with new definitions."
)
rd_add_wide_ <- "`data.frame` Household data used for analysis (wide format).
Optional at construction; if omitted, attach it later with
[add_data_wide()][housepi::add_data_wide]."
rd_hhid_ <- "`character(1)` ID of the household of interest."
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
rd_var_discr_ <- function(what, from) {
  sprintf("`logical(1)` Treat %s definition as discrete, if applicable.
          Inferred from %s if `NA` (default).", what, from)
}
rd_obs_discr_from_ <- "the observed data column, or the model requirements"
rd_var_prior_ <- "`character(1)` Prior distribution."
rd_var_val_ <- function(what) {
  sprintf("Value the definition is pinned to in simulations. In inference, it
          pins a constant or deterministic definition, and is ignored when a
          `%s` is given.", what)
}
rd_var_args_ <- function(what, proposal = FALSE, func = FALSE) {
  sink <- if (func) ", or for `func` when the definition is deterministic."
  else "."
  paste0(sprintf("named arguments for the %s distribution", what), sink,
         if (proposal) rd_proposal_opts_)
}
rd_var_min_ <- "`numeric(1)` Truncate to a minimum value.
Relative to the location argument `loc`, when one is given."
rd_var_max_ <- "`numeric(1)` Truncate to a maximum value.
Relative to the location argument `loc`, when one is given."
rd_var_func_ <- "One between:
- `character(1)` Name of latent function;
- rlang expression."
rd_var_proposal_ <- "`character(1)` Proposal used in the Monte Carlo.
Random walks changes current value:
\\describe{
  \\item{`\"rw\"`}{normal increments (default);}
  \\item{`\"rw_log\"`}{lognormal log-increments (positive, non-discrete
    values only);}
  \\item{`\"rw_unif\"`}{uniform increments (bounded).}
}
Independence samplers draw a fresh value from a distribution:
\\describe{
  \\item{`\"forward\"`}{drawn from the variable's own defining distribution
    (its prior for a parameter, its model for a latent variable);}
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
rd_sources_ <- "Either a plain value, the name of another definition,
or an arithmetic expression combining them. A name is looked up among the
latent variables, data columns, observations, group risks and parameters the
model carries: two of them sharing a name cannot be read."
rd_aug_status_ <- "`integer(1)` Individual positive status.
Used to specify status-specific definitions.
If `NA`, associate definition to all infectious individuals (default)."
rd_aug_val_ <- paste("`numeric(1)` | `character(1)` Value the definition is
pinned to in simulations. In inference, it pins a constant or deterministic
definition, and is ignored when the latent is sampled, with a `model` or a
`proposal`.", rd_sources_)
rd_obs_val_ <- paste("`numeric(1)` | `character(1)` Definition of the
observation, depending on `level`:
\\describe{
  \\item{`\"member\"`}{without a `model`, the value computed for every member;
  with a `model`, the datum the model is evaluated at, the like-named data
  column when omitted;}
  \\item{`\"house\"`}{one value for the whole household. The reductions
  `min()`, `max()`, `sum()` and `prod()` fold a per-member value (data column,
  latent or group risk) over the members: their argument is computed for each
  member, and members with no value are left out, e.g. `min(onset_date) + 21`
  or `max(test_date - onset_date)`. A per-member value read outside a
  reduction is an error.}
}", rd_sources_)
rd_aug_min_ <- paste(rd_var_min_, rd_sources_)
rd_aug_max_ <- paste(rd_var_max_, rd_sources_)
rd_lat_flat_bound_ <- "Without a `model`, the bound the latent is sampled within
in inference, reading only data, observations and constant parameters."

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
rd_import_ <- "Named `list` of pre-built definitions to import.
Cannot be combined with single and/or multiple definitions."
rd_tmpl_import_ <- "Named `list` of pre-built definitions to import,
or `\"template\"` for the definitions shipped with the model.
The template is R code evaluated from the model folder, so use it only with
model folders you trust.
Cannot be combined with single and/or multiple definitions."
rd_list_append_ <- "`logical(1)` Append to current definitions.
If `FALSE` drop all pre-existing definitions."

## Risks ##

rd_risk_desc_ <- function(kind) {
  sprintf("Add definitions of a %s risk by given data covariate(s).
Each stratum carries its own absolute risk, unless a reference is chosen:
then every parameter multiplies the baseline risk of the reference stratum.
If no MCMC options are provided, absolute risks are assumed
\\eqn{\\mathrm{uniform}} over their support, and relative ones
\\eqn{\\mathrm{lognormal}[\\mu=0, \\sigma=3]}.
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
Every member, or ordered pair of members for contact risks, must match a row,
and every row must be the best match of at least one of them,
else the definition is refused. So is a table setting no covariate value,
which treats every member alike as the baseline does.
The table is used as written.
A risk read on infected members is checked at every infected status a member
may take, whether or not the data records it.
Rows keyed on different covariates form separate rule groups,
whose risks multiply.
Use `.reference` in `param_name` to designate the reference group;
a rule group without it leaves the risk absolute.
The reserved column `value` pins a stratum to that value in simulation,
while it keeps its prior and is sampled in inference.
The reserved column `const` fixes a pinned stratum in inference too,
so that it carries no prior. The reference stratum takes neither.
Mutually exclusive with `by`."
rd_risk_params_ <- "named `list` Per-parameter MCMC setups as in `param_name`.
Overrides shared prior/args for parameters specified in the list.
A stratum pinned through the `strata` table takes no entry here.
Ignored if automatic definition."
rd_risk_reference_ <- function(cell) {
  sprintf("named `list` Covariates associated to the reference %s,
whose risk becomes the unit that every other parameter multiplies.
Pass `\".auto\"` to take the %s with the largest sample size.
Left empty, each stratum carries its own absolute risk; a covariate the
reference does not name leaves the whole risk absolute.", cell, cell)
}

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
rd_nrep_ <- "`integer(1)` Maximum number of `param_table` rows simulated.
Posterior draws are thinned to their ESS when `thin_ess` is `TRUE`.
Then `nrep` evenly spaced rows are kept (all of them when fewer remain).
If `NA`, every row is simulated (default)."
rd_thin_ess_ <- "`logical(1)` Whether to thin `param_table` to its
bulk effective sample size (ESS).
Ignored if `param_table` has no posterior draws."
rd_simul_observ_ <- "`function` Observation process.
Function expects as input the `data.frame` generated by simulation.
If `NULL`, observation process is ignored."
rd_simul_obs_long_ <- "`function` Optional function.
Expects processed `data_wide` and returns a `data_long` data.frame.
When provided, each simulation returns `list(wide = ..., long = ...)`."
rd_natural_log_ <- "Logarithms are natural (base \\eqn{e})."
rd_global_summary_ <- "Results are summarised at the point estimate and the
credible level held by the global options;
see [set_housepi_summary()]."

## State inspection ##

rd_inspect_state_ <- "The values a `HousepiModel` holds are left by the last
call that ran it:
\\describe{
  \\item{[evaluate_log_lik()]}{the evaluated parameters and latent variables,
  with their log-likelihood contributions.}
  \\item{[run_inference()]}{the final iteration of the last chain.}
  \\item{[run_simulation()]}{the observed data, the parameters of the last
  replicate.}
}"

# nocov end

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

# Ensure object of class HousepiModel.
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
#' An empty cell of a categorical column is missing data, not a category of
#' its own: recode it to a label of your choice to stratify on it.
#' @param hhmodel `r rd_model_in_`
#' @param data `r rd_add_wide_`
#' @return `r rd_model_out_`
#' @examples
#' # construct the model first, attach the data afterwards
#' hhmodel <- new_housepi_model(model = "final_size") |>
#'   add_data_wide(flu_1918_frederick)
#' summarise_house_cases(hhmodel)
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
#' @examples
#' # repeated measures of one individual, matched on `id`
#' measures <- data.frame(id = "hh001m02", day = 1:3,
#'                        viral_load = c(5.2, 6.8, 4.1))
#' hhmodel <- sarscov2_2021_icofs |>
#'   new_housepi_model("case-ascertained") |>
#'   add_data_long(measures)
#' get_data_long(hhmodel)
#' @export
add_data_long <- function(hhmodel, data_long) {
  stopifnothousepi(hhmodel)
  new_obj <- hhmodel$copy()
  pvt <- new_obj$.__enclos_env__$private
  pvt$gen_$import_data_long(data_long)
  pvt$rebuild_()
  new_obj
}

#' @title Summarise observed cases in the loaded household data
#' @description
#' Summarise cases as given from the statuses in the loaded household data.
#' The row `classes`, counts the equivalence classes:
#' households indistinguishable by the model (same size, statuses, risk, etc.).
#' @param hhmodel `r rd_model_in_`
#' @return `data.frame`, one row per metric:
#' \describe{
#'   \item{level}{`individual` or `household`}
#'   \item{metric}{quantity counted in the row}
#'   \item{count}{absolute count}
#'   \item{percent}{percentage of the relevant total, `NA` for `classes`}
#' }
#' @examples
#' hhmodel <- new_housepi_model(flu_1918_frederick, "final_size")
#' summarise_house_cases(hhmodel)
#' @export
summarise_house_cases <- function(hhmodel) {
  stopifnothousepi(hhmodel)
  pvt <- hhmodel$.__enclos_env__$private
  .engine$summarise_house_cases(pvt$gen_$mcmc_ptr)
}

#' @rdname summarise_house_cases
#' @export
summarize_house_cases <- summarise_house_cases

############################################################################
############################################################################
############################################################################

# Keep the definitions whose column `col` matches one of `name`.
# All of them when `name` is `NA`.
subset_definitions_ <- function(defs, name, col = "name") {
  if (length(name) == 1L && is.na(name)) return(defs)
  stopifnot(is.character(name), !anyNA(name))
  unmatched <- setdiff(name, defs[[col]])
  if (length(unmatched)) stop_h(
    "No definition{?s} named {.val {unmatched}}.",
    find_closest_match(unmatched, defs[[col]])
  )
  out <- defs[defs[[col]] %in% name, , drop = FALSE]
  rownames(out) <- NULL
  out
}

#' @title Retrieve currently loaded data and definitions from a housepi model
#' @description
#' Return the definitions currently stored in the housepi model.
#' Each variant returns one category (data, parameters, latent variables,
#' observations, group risks, contact risks, infectivity profile).
#' @details
#' `get_data_wide()` returns the individual-level data the model reads: a column
#' a deterministic observation stands in for carries the derived value, not the
#' imported one. A column an observation produces without standing in for data
#' is not part of the dataset, and is reported by [get_observations()].
#' @param hhmodel `r rd_model_in_`
#' @param name `character()`. Names of the definitions to return.
#' `get_group_risks()` matches the name of the group risk.
#' All of them when `NA` (default).
#' @return `data.frame`
#' @examples
#' hhmodel <- new_housepi_model(flu_1918_frederick, "final_size") |>
#'   add_parameter(import = "template")
#' get_data_wide(hhmodel) |>
#'   head()
#' get_parameters(hhmodel)
#' get_parameters(hhmodel, name = c("alpha", "beta"))
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
get_parameters <- function(hhmodel, name = NA) {
  stopifnothousepi(hhmodel)
  subset_definitions_(hhmodel$parameters, name)
}

#' @rdname get_definitions
#' @export
get_latent <- function(hhmodel, name = NA) {
  stopifnothousepi(hhmodel)
  subset_definitions_(hhmodel$latent, name)
}

#' @rdname get_definitions
#' @export
get_observations <- function(hhmodel, name = NA) {
  stopifnothousepi(hhmodel)
  subset_definitions_(hhmodel$observations, name)
}

#' @rdname get_definitions
#' @export
get_group_risks <- function(hhmodel, name = NA) {
  stopifnothousepi(hhmodel)
  subset_definitions_(hhmodel$group_risks, name, "group_risk")
}

#' @rdname get_definitions
#' @export
get_contact_risks <- function(hhmodel) {
  stopifnothousepi(hhmodel)
  hhmodel$contact_risks
}

#' @rdname get_definitions
#' @export
get_infectivity_profiles <- function(hhmodel) {
  stopifnothousepi(hhmodel)
  hhmodel$infectivity_profiles
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
assert_bulk_defs_ <- function(family, dots, has_template = FALSE) {
  if (length(dots)) return(dots)
  sms <- c("Missing or empty {family} definitions.", definition_helper_)
  if (isTRUE(has_template)) sms <- c(sms, template_helper_)
  stop_h(sms)
}

# Name the formals given in an `add_*` call that only single mode reads.
# `listed` are the formals the list-like modes read besides `...`, `import`
# and `append`.
single_args_ <- function(call, listed) {
  setdiff(names(call)[-1], c("hhmodel", "...", "import", "append", listed))
}

# Reject a `...` definition bound to a formal only single mode reads.
assert_no_reserved_defs_ <- function(singles) {
  if (length(singles)) stop_h(
    paste("A definition can't be named {.or {.arg {singles}}},",
          "{cli::qty(singles)}{?an argument/arguments}",
          "reserved to single mode."),
    "Define it in single mode, or within an {.arg import} list."
  )
}

# Reject `import` combined with single-mode arguments or inline `...` pairs.
assert_no_def_mixing_ <- function(key, dots, singles) {
  if (!is.null(key) || length(dots) || length(singles)) stop_h(
    "Mixing imported definitions with single and multiple definitions",
    definition_helper_
  )
}

# Build the quosure list for the explicit `import` bulk form.
defs_to_quos_ <- function(import, key, dots, singles, env) {
  assert_no_def_mixing_(key, dots, singles)
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
add_template_defs_ <- function(holder, generator, family, append, key, dots,
                               singles) {
  assert_no_def_mixing_(key, dots, singles)
  if (generator$is_test_model()) {
    warn_h("No template {family} definitions available for test model.")
  } else {
    defs <- get_template_def_list(generator$source_dir, family)
    holder$add_quo_list(convert_list_to_quosures(defs), append = isTRUE(append))
  }
}

# Add the definitions of the list-like `add_*` modes: `import = "template"`,
# an `import` list, or the inline `...` pairs of the multiple mode.
# Returns `FALSE` in single mode (`key` given), left to the caller as it
# captures its own arguments; `...` goes to `$add_quo_list()`.
# `call` is the caller's `match.call(expand.dots = FALSE)`, and `listed` the
# formals these modes read on top of `...`, `import` and `append`: any other
# formal the call sets belongs to single mode.
add_listed_defs_ <- function(holder, generator, family, key, import, append,
                             dots, env, call, listed,
                             has_template = FALSE, ...) {
  singles <- single_args_(call, listed)
  if (has_template && is_template_import_(import)) {
    add_template_defs_(holder, generator, family, append, key, dots, singles)
  } else if (length(import)) {
    holder$add_quo_list(defs_to_quos_(import, key, dots, singles, env),
                        append = append, ...)
  } else if (is.null(key)) {
    assert_no_reserved_defs_(singles)
    holder$add_quo_list(assert_bulk_defs_(family, dots, has_template),
                        append = append, ...)
  } else {
    return(FALSE)
  }
  TRUE
}

#' @title Add or overwrite parameter definitions
#' @description
#' Add parameter(s) to a housepi model.
#' If a parameter with the same name exists, the former gets overwritten.
#' `r rd_add_modes_`
#' @param hhmodel `r rd_model_in_`
#' @param name `r rd_var_name_`
#' @param value `numeric(1)` `r rd_var_val_("prior")`
#' @param const `r rd_var_const_`
#' @param discr `r rd_var_discr_("prior", "the model requirements")`
#' @param prior `r rd_var_prior_`
#' @param ... `r rd_var_args_("prior", proposal = TRUE, func = TRUE)`
#' In mode 'multiple', the `name = value` definitions.
#' @param min `r rd_var_min_`
#' @param func `r rd_var_func_`
#' @param max `r rd_var_max_`
#' @param proposal `r rd_var_proposal_`
#' @param import `r rd_tmpl_import_`
#' @param append `r rd_list_append_`
#' @return `r rd_model_out_`
#' @examples
#' # the single, inline-bulk and list-bulk forms compose freely in a pipe
#' hhmodel <- flu_1918_frederick |>
#'   new_housepi_model("final-size") |>
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
#'   new_housepi_model("final-size") |>
#'   add_parameter(import = my_import)
#'
#' # `"template"` loads the definitions shipped with the model
#' hhmodel <- flu_1918_frederick |>
#'   new_housepi_model("final-size") |>
#'   add_parameter(import = "template")
#' get_parameters(hhmodel)
#' @export
add_parameter <- function(hhmodel,
                          name = NULL, value = NULL, ...,
                          const = NA, discr = NA,
                          prior = NA_character_,
                          func = NULL,
                          min = NULL, max = NULL,
                          proposal = NA_character_,
                          import = NULL, append = TRUE) {
  stopifnothousepi(hhmodel)
  new_obj <- hhmodel$copy()
  pvt <- new_obj$.__enclos_env__$private

  #> template, import and multiple modes
  if (!add_listed_defs_(pvt$param_, pvt$gen_, "parameter", name, import,
                        append, rlang::enquos(...), rlang::caller_env(),
                        match.call(expand.dots = FALSE), "name",
                        has_template = TRUE))
    #> single mode
    pvt$param_$add_quo(
      name = name, value = rlang::enquo(value),
      const = const, discr = discr,
      prior = prior, rlang::enquos(...),
      func = rlang::enquo(func),
      min = rlang::enquo(min), max = rlang::enquo(max),
      proposal = proposal,
      .assert = TRUE
    )
  pvt$rebuild_()
  new_obj
}

############################################################################
############################################################################
############################################################################

#' @title Add or update latent variable definitions
#' @description
#' Add latent variables/functions to housepi model.
#' If a latent definitions with the same name and status exists,
#' the former gets overwritten.
#' `r rd_add_modes_`
#' @details
#' A latent with a `model` is sampled from it.
#' A latent the transmission model generates, such as the infection times of
#' [template_case_ascertained], takes no `model`: transmission gives its
#' probability, and inference samples it flat within `min` and `max`.
#' In simulation, transmission sets it, except where it does not, such as the
#' first infection of a household, which takes `value`.
#' @param hhmodel `r rd_model_in_`
#' @param name `r rd_var_name_`
#' @param value `r rd_aug_val_`
#' @param status `r rd_aug_status_`
#' @param const `r rd_var_const_`
#' @param discr `r rd_var_discr_("model", "the model requirements")`
#' @param func `r rd_var_func_`
#' @param model `character(1)` Latent model distribution.
#' Omitted for a latent the transmission model generates (see Details).
#' @param ... `r rd_var_args_("model", proposal = TRUE, func = TRUE)`
#' In mode 'multiple', the `name = value` definitions.
#' @param min `r rd_aug_min_`
#' `r rd_lat_flat_bound_`
#' @param max `r rd_aug_max_`
#' `r rd_lat_flat_bound_`
#' @param proposal `r rd_var_proposal_`
#' @param import `r rd_tmpl_import_`
#' @param append `r rd_list_append_`
#' @return `r rd_model_out_`
#' @examples
#' hhmodel <- dummy_household |>
#'   new_housepi_model("test") |>
#'   # define an incubation period drawn from a gamma distribution
#'   add_latent("incub_period",
#'                 model = "gamma", shape = 2, scale = 1, max = 14)
#' @export
add_latent <- function(hhmodel,
                       name = NULL, value = NULL, ...,
                       status = NA_integer_, const = NA, discr = NA,
                       func = NULL,
                       model = NA_character_,
                       min = NULL, max = NULL,
                       proposal = NA_character_,
                       import = NULL, append = TRUE) {
  stopifnothousepi(hhmodel)
  new_obj <- hhmodel$copy()
  pvt <- new_obj$.__enclos_env__$private

  #> template, import and multiple modes
  if (!add_listed_defs_(pvt$augm_, pvt$gen_, "latent", name, import,
                        append, rlang::enquos(...), rlang::caller_env(),
                        match.call(expand.dots = FALSE), c("name", "status"),
                        has_template = TRUE, status = status))
    #> single mode
    pvt$augm_$add_quo(
      name = name, value = rlang::enquo(value),
      status = status, const = const, discr = discr,
      func = rlang::enquo(func),
      model = model, rlang::enquos(...),
      min = rlang::enquo(min), max = rlang::enquo(max),
      proposal = proposal,
      .assert = TRUE
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
#' An individual whose datum is missing is left out of that likelihood, with a
#' warning: there is nothing to condition on, and the observation constrains
#' neither the individual latent variables nor the parameters. A column holding
#' no value at all is an error.
#'
#' At household level, the observation is a reduction across all members,
#' it's shared by the whole household, and must reduce to a constant value.
#'
#' An observation carrying no `model` produces its own data column:
#' it is evaluated for every member, and overwrites data of the same name.
#'
#' `r rd_add_modes_`
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
#' @param value `r rd_obs_val_`
#' @param func `r rd_var_func_`
#' @param status `r rd_aug_status_`
#' @param discr `r rd_var_discr_("model", rd_obs_discr_from_)`
#' @param min `r rd_aug_min_`
#' @param max `r rd_aug_max_`
#' @param import `r rd_import_`
#' @param append `r rd_list_append_`
#' @return `r rd_model_out_`
#' @examples
#' hhmodel <- sarscov2_2021_icofs |>
#'   new_housepi_model("case-ascertained") |>
#'   add_parameter(import = "template") |>
#'   add_infectivity_profile(import = "template") |>
#'   add_latent(
#'     infec_time = list(min = start_date - 10, max = end_date),
#'     incub_period = list(model = "gamma", mean = 3, sd = 1.4),
#'     onset_time = infec_time + incub_period
#'   ) |>
#'   # member level: the test date is noisy around the onset time
#'   add_observation("test_date", model = "normal",
#'                   mean = onset_time, sd = 1, min = floor(onset_time)) |>
#'   # house level: the earliest positive test among members
#'   add_observation("first_test", value = min(test_date), level = "house")
#' get_observations(hhmodel)
#' @export
add_observation <- function(hhmodel, name = NULL,
                            model = NA_character_, ...,
                            level = c("member", "house"), value = NULL,
                            func = NULL, status = NA_integer_, discr = NA,
                            min = NULL, max = NULL,
                            import = NULL, append = TRUE) {
  stopifnothousepi(hhmodel)
  # bulk default: a per-definition 'level =' overrides it (see add_quo)
  level <- match.arg(level)
  new_obj <- hhmodel$copy()
  pvt <- new_obj$.__enclos_env__$private

  #> import and multiple modes
  if (!add_listed_defs_(pvt$augm_, pvt$gen_, kind2family(level), name,
                        import, append, rlang::enquos(...),
                        rlang::caller_env(), match.call(expand.dots = FALSE),
                        c("name", "status", "level"),
                        status = status, .kind = level))
    #> single mode
    pvt$augm_$add_quo(
      name = name, value = rlang::enquo(value), status = status,
      const = NA, discr = discr, func = rlang::enquo(func),
      model = model, rlang::enquos(...),
      min = rlang::enquo(min), max = rlang::enquo(max),
      .assert = TRUE, .kind = level
    )
  pvt$rebuild_()
  new_obj
}

#' @title Add a status process
#' @description
#' Declare how a member's infected status is generated during simulation, as a
#' categorical distribution over infected statuses. The probabilities are given
#' in declaration order (the \eqn{k}-th probability is that of the \eqn{k}-th
#' infected status); when named, the names become the status labels.
#' During evaluation the status is fixed observed data, so it contributes a
#' constant and the probabilities must be data or numeric constants.
#' @param hhmodel `r rd_model_in_`
#' @param ... Probabilities of the infected statuses, in order. When named, the
#' names label the corresponding statuses: name all of them or none.
#' @return `r rd_model_out_`
#' @examples
#' # simulate infected members as symptomatic (70%) or asymptomatic (30%)
#' simulated <- flu_1918_frederick |>
#'   new_housepi_model("final-size") |>
#'   add_parameter(import = "template") |>
#'   add_status(symptomatic = 0.7, asymptomatic = 0.3) |>
#'   run_simulation()
#' table(simulated$status)
#' @export
add_status <- function(hhmodel, ...) {
  stopifnothousepi(hhmodel)
  new_obj <- hhmodel$copy()
  pvt <- new_obj$.__enclos_env__$private
  weights <- rlang::enquos(...)
  if (!length(weights)) stop_h(
    "No status probability provided.",
    "Provide at least one infected-status probability."
  )
  labels <- rlang::names2(weights)
  if (any(nzchar(labels)) && !all(nzchar(labels))) stop_h(
    "Status labels were given to some probabilities only.",
    "Name every status probability, or none."
  )
  # a single reserved node carries the status process (there is one status)
  pvt$augm_$add_quo(
    name = ".status", value = NULL, status = NA_integer_,
    const = NA, discr = FALSE,
    model = "categorical", weights,
    .assert = TRUE, .kind = "status"
  )
  # record the ordered status labels (the dictionary lives here, not in C++)
  pvt$augm_$set_status_levels(labels)
  pvt$rebuild_()
  new_obj
}

############################################################################
############################################################################
############################################################################

#' @title Add a stratified group risk definition
#' @description `r rd_risk_desc_("group")`
#' `r rd_add_modes_`
#' @param hhmodel `r rd_model_in_`
#' @param name `r rd_var_name_`
#' @param by `list` Covariates associated to group risk as data headers.
#' Corresponding data columns have to be categorical.
#' @param reference `r rd_risk_reference_("group")`
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
#' @examples
#' hhmodel <- new_housepi_model(sarscov2_2021_icofs, "case_ascertained")
#'
#' # single: one susceptibility parameter per age class
#' hhmodel |>
#'   add_group_risk("susceptibility", by = "age") |>
#'   get_group_risks()
#'
#' # inline bulk: relative risks, against a chosen or the most frequent group
#' hhmodel |>
#'   add_group_risk(
#'     infectivity = list(by = "isolated", reference = ".auto"),
#'     susceptibility = list(by = "vacc_status", reference = "unvaccinated")
#'   ) |>
#'   get_parameters()
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

  #> import and multiple modes
  if (!add_listed_defs_(pvt$group_risks_, pvt$gen_,
                        pvt$group_risks_$family, name, import, append,
                        rlang::enquos(...), rlang::caller_env(),
                        match.call(expand.dots = FALSE), "name"))
    #> single mode
    pvt$group_risks_$add_quo(
      name = name, by = by,
      reference = reference, joint = joint,
      strata = strata, params = params,
      prior = prior, rlang::enquos(...),
      min = rlang::enquo(min), max = rlang::enquo(max),
      proposal = proposal, prefix = prefix
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
#' @param reference `r rd_risk_reference_("pair")`
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
#' @examples
#' hhmodel <- new_housepi_model(flu_2014_matsumoto, "final_size")
#'
#' # one parameter per pair of family roles, relative to the most frequent pair
#' hhmodel |>
#'   add_contact_risk("member", reference = ".auto") |>
#'   get_contact_risks()
#'
#' # asymmetric: infector-infectee and infectee-infector pairs differ
#' hhmodel |>
#'   add_contact_risk("member", symmetric = FALSE) |>
#'   get_contact_risks()
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
  pvt$assert_feature_supported_(pvt$gen_$has_contact_risk,
                                "Contact risks")
  # the definition is replaced: its parameters are re-derived by `rebuild_()`
  pvt$drop_stale_risk_params_()

  # contact risk holds a single definition:
  # import a named arg list, or the explicit arguments
  if (!length(import)) import <- list(
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
#' @param hhmodel `r rd_model_in_`
#' @param func `r rd_var_name_`
#' @param status `r rd_aug_status_`
#' @param discr `r rd_var_discr_("profile", "the model time discretization")`
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
#'   new_housepi_model("case-ascertained") |>
#'   add_infectivity_profile(import = "template")
#' get_infectivity_profiles(hhmodel)
#' @export
add_infectivity_profile <- function(hhmodel,
                                    func = NULL, ...,
                                    status = NA_integer_, discr = NA,
                                    min = NULL, max = NULL,
                                    import = NULL, append = TRUE) {
  stopifnothousepi(hhmodel)
  new_obj <- hhmodel$copy()
  pvt <- new_obj$.__enclos_env__$private
  pvt$assert_feature_supported_(pvt$gen_$has_infec_profile,
                                "Infectivity profiles")
  # impose discrete-time if it's a model choice
  pvt$infec_profile_$discr_model <- pvt$gen_$is_discrete_time

  #> template, import and multiple modes
  if (!add_listed_defs_(pvt$infec_profile_, pvt$gen_,
                        "infectivity profile", func, import, append,
                        rlang::enquos(...), rlang::caller_env(),
                        match.call(expand.dots = FALSE), "func",
                        has_template = TRUE))
    #> single mode
    pvt$infec_profile_$add_quo(
      func = func, status = status, discr = discr,
      rlang::enquos(...),
      min = rlang::enquo(min), max = rlang::enquo(max),
      .assert = TRUE
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
#' For models with latent variables
#' this is conditioned on current latent values.
#' @param hhmodel `r rd_model_in_`
#' @param hhid `r rd_hhid_`
#' If `NA` (default) it computes for the whole data.
#' @return `numeric(1)`. The complete-data log-likelihood.
#' `r rd_natural_log_`
#' @examples
#' flu_1918_frederick |>
#'   new_housepi_model("final-size") |>
#'   add_parameter(import = "template") |>
#'   evaluate_log_lik()
#' @export
evaluate_log_lik <- function(hhmodel, hhid = NA_character_) {
  stopifnothousepi(hhmodel)
  pvt <- hhmodel$.__enclos_env__$private
  pvt$build_full_()
  .engine$evaluate_log_lik(pvt$gen_$mcmc_ptr, hhid)
}

#' @title Inspect the augmented data of a housepi model
#' @description
#' Read the dataset of a `HousepiModel` with its latent variables.
#' `r rd_inspect_state_`
#' @param hhmodel `r rd_model_in_`
#' @return `data.frame`, the dataset with its latent variables.
#' @seealso [inspect_house_state()]
#' @examples
#' hhmodel <- flu_1918_frederick |>
#'   new_housepi_model("final-size") |>
#'   add_parameter(import = "template")
#' evaluate_log_lik(hhmodel)
#' inspect_augmented_data(hhmodel) |>
#'   head()
#' @export
inspect_augmented_data <- function(hhmodel) {
  stopifnothousepi(hhmodel)
  gen <- hhmodel$.__enclos_env__$private$gen_
  gen$refactor_data_wide(.engine$export_augmented_data_wide(gen$mcmc_ptr))
}

#' @title Inspect the state of a household in a housepi model
#' @description
#' Read all the variables a `HousepiModel` holds for one household.
#' `r rd_inspect_state_`
#' @param hhmodel `r rd_model_in_`
#' @param hhid `r rd_hhid_`
#' @return `list`, all variables of the household.
#' Entries are named after the homonymous active fields, followed by the
#' latent, observation and transmission log-likelihood contributions.
#' An empty definition is dropped unless the model requires it.
#' @seealso [inspect_augmented_data()]
#' @examples
#' hhmodel <- flu_1918_frederick |>
#'   new_housepi_model("final-size") |>
#'   add_parameter(import = "template")
#' evaluate_log_lik(hhmodel)
#' inspect_house_state(hhmodel, "1")
#' @export
inspect_house_state <- function(hhmodel, hhid) {
  stopifnothousepi(hhmodel)
  gen <- hhmodel$.__enclos_env__$private$gen_
  out <- .engine$export_house_snapshot(gen$mcmc_ptr, as.character(hhid))
  out$data_wide <- gen$refactor_data_wide(out$data_wide)
  ## data long
  if (all(dim(out$data_long) == 0)) out$data_long <- NULL
  else out$data_long <- gen$refactor_data_long(out$data_long)
  ## parameters
  if (!length(gen$param_required) &&
        all(dim(out$parameters) == 0)) out$parameters <- NULL
  ## group risks
  if (!length(gen$group_risks_required) &&
        all(dim(out$group_risks) == 0)) out$group_risks <- NULL
  ## contact risks
  if (!gen$has_contact_risk &&
        all(dim(out$contact_risks) == 0)) out$contact_risks <- NULL
  ## infectivity profiles
  if (!gen$has_infec_profile && !length(out$infectivity_profiles)) {
    out$infectivity_profiles <- NULL
  } else {
    lst <- out$infectivity_profiles
    df <- data.frame(id = lst$id, discr = lst$discr, func = lst$func,
                     min = lst$min, max = lst$max)
    df$func_args <- lst$func_args
    out$infectivity_profiles <- df
  }
  out
}

#' @title Run Bayesian MCMC inference on a housepi model
#' @description
#' Run inference with Bayesian Markov Chain Monte Carlo, by the
#' component-wise Metropolis-Hastings algorithm
#' \insertCite{metropolis_equation_1953,hastings_monte_1970}{housepi}.
#' @param hhmodel `r rd_model_in_`
#' @param steps `r rd_smpl_steps_`
#' @param burn `r rd_smpl_burn_`
#' @param thin `r rd_smpl_thin_`
#' @param adapt `r rd_smpl_adapt_`
#' @param chains `r rd_smpl_chain_`
#' @param threads `r rd_thread_`
#' @param id_draws Infected members whose public latent draws are stored.
#' Either a `character()` of individual identifiers,
#' or a `list` with entries `id` (individual identifiers)
#' and/or `hhid` (whole household infected members are stored).
#' If `NA`, one random eligible household is followed (default).
#' Pass `".none"` to store nothing.
#' @param progress `r rd_progress_`
#' @param warn `logical(1)` Whether to report the warnings the run raises.
#' A quiet run still records the convergence diagnostics it flagged, read back
#' from the `warnings` field of the returned object.
#' @param ... `r rd_dots_config_`
#' @return An object of class `HousepiPosteriorDraws`.
#' @details
#' The run warns on any of: `r paste(diagnostic_flags_, collapse = "; ")`.
#' R-hat and bulk ESS are the rank-normalized ones of
#' \insertCite{vehtari_rank_normalization_2021;textual}{housepi}.
#' @references
#' \insertAllCited{}
#' @examples
#' flu_1918_frederick |>
#'   new_housepi_model("final-size") |>
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
                          warn = TRUE,
                          ...) {
  stopifnothousepi(hhmodel)
  # the whole run is quiet, its build included
  local_quiet_warnings_(!warn)
  ptr <- prepare_inference_(hhmodel, steps, burn, thin, adapt, chains,
                            threads, id_draws, progress, draw_seed_(), ...)
  wrap_inference_(hhmodel, .engine$run_inference_list(list(ptr), 1L))
}

#' Build the engine of `hhmodel` for inference, seeded at `seed`.
#' Shared by [run_inference()] and the replicate sweep.
#' @return The engine pointer.
#' @noRd
prepare_inference_ <- function(hhmodel, steps, burn, thin, adapt, chains,
                               threads, id_draws, progress, seed, ...) {
  assert_advanced_setup_("adapt_window", ...)
  pvt <- hhmodel$.__enclos_env__$private
  id_draws <- resolve_id_draws_(id_draws)
  pvt$assert_posterior_memory_(
    steps = steps, burn = burn, thin = thin,
    chains = chains, id_draws = id_draws
  )
  pvt$build_full_(
    steps = steps, burn = burn, thin = thin, adapt = adapt,
    chains = chains, threads = threads, id_draws = id_draws,
    progress = progress, seed = seed, ...
  )
  pvt$gen_$mcmc_ptr
}

#' Posterior draws of the engine of `hhmodel`, entry `i` of the engine run
#' `ran`, or the error that stopped the run, rendered as any engine error.
#' @noRd
wrap_inference_ <- function(hhmodel, ran, i = 1L) {
  if (!is.null(ran$error)) with_cpp_conditions(stop(ran$error))
  gen <- hhmodel$.__enclos_env__$private$gen_
  HousepiPosteriorDraws$new(gen$raw_name, gen$mcmc_ptr, ran$seconds[[i]])
}

#' @title Simulate a household survey from a housepi model
#' @description
#' Run simulation to generate a household survey.
#' It assumes the same community structure of the loaded data.
#' @param hhmodel `r rd_model_in_`
#' @param obs_process `r rd_simul_observ_`
#' @param obs_process_long `r rd_simul_obs_long_`
#' @param param_table `r rd_param_table_`
#' @param nrep `r rd_nrep_`
#' @param thin_ess `r rd_thin_ess_`
#' Without `param_table`, the number of replicates, which share the pinned
#' parameters and draw the others from their prior. If `NA`, one replicate.
#' Once set, even to 1, the replicates are returned as a list.
#' @param ... `r rd_dots_config_`
#' @return A `data.frame` for a single replicate, when `nrep` is `NA` and
#' `param_table` is `NULL` or a one-row `data.frame`.
#' Otherwise a `list` with one `data.frame` per replicate.
#' With `obs_process_long`, a `list(wide = , long = )` of either.
#' @examples
#' flu_1918_frederick |>
#'   new_housepi_model("final-size") |>
#'   add_parameter(import = "template") |>
#'   run_simulation() |>
#'   head()
#' @export
run_simulation <- function(hhmodel,
                           obs_process = NULL,
                           obs_process_long = NULL,
                           param_table = NULL,
                           nrep = NA_integer_,
                           thin_ess = TRUE,
                           ...) {
  output <- simulate_with_parameters_(
    hhmodel, obs_process = obs_process, obs_process_long = obs_process_long,
    param_table = param_table, thin_ess = thin_ess, nrep = nrep, ...
  )
  output$parameters <- NULL
  # the arguments, not the replicates left after thinning, ask for one frame
  single <- is.null(param_table) || (
    is.data.frame(param_table) && !posterior::is_draws(param_table) &&
      nrow(param_table) == 1L
  )
  if (is.na(nrep) && single) output <- lapply(output, `[[`, 1L)
  if (length(output) == 1) output[[1]] else output
}

#' Simulate from `hhmodel` as [run_simulation()] does, returning the lists of
#' `wide` and `long` surveys, one frame per replicate, and `parameters`: the
#' values each replicate is simulated at, pinned or drawn from the prior, one
#' row per replicate.
#' @keywords internal
#' @noRd
simulate_with_parameters_ <- function(hhmodel,
                                      obs_process = NULL,
                                      obs_process_long = NULL,
                                      param_table = NULL,
                                      nrep = NA_integer_,
                                      thin_ess = TRUE,
                                      ...) {
  stopifnothousepi(hhmodel)
  # simulation specific build
  assert_advanced_setup_("time_step", ...)
  pvt <- hhmodel$.__enclos_env__$private
  pvt$helper_simulation_build_(
    param_table = param_table, nrep = nrep, thin_ess = thin_ess, ...
  )

  # combine observation process to measures and refactor
  obs_func <- pvt$helper_obs_process_(obs_process)
  # run simulations
  output <- .engine$run_simulation(
    pvt$gen_$mcmc_ptr,
    obs_func, NULL,
    obs_process_long
  )

  # drop empty long
  if (!length(output$long)) output$long <- NULL
  output
}
