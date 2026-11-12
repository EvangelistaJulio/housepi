## housepi - Epidemiological Analysis of Household Outbreaks.
## Copyright 2024-2026 Institut Pasteur (Paris), France.
## SPDX-License-Identifier: GPL-3.0-or-later
##

########  Documentation  ########

# nocov start

####  Model surface  ####

rd_model_in_ <- "`HousepiModel` object."
rd_model_out_ <- paste(
  rd_model_in_,
  "Deep clone of input with new definitions."
)
rd_dataset_ <- "`data.frame` Household data, one row per member.
Column `status` holds integer codes (`0` for escaped members), or labels.
Definitions keyed on a status name it the same way."
rd_mdl_choice_ <- "`character(1)` Model to use, one of:
\\describe{
  \\item{`<template>`}{a built-in template name;}
  \\item{`path/to/source`}{folder with implementation file `\"model.cpp\"`;}
  \\item{`path/to/compiled`}{folder with a previously compiled model;}
  \\item{`\"test\"`}{template name for test model;}
  \\item{`NA`}{same as `\"test\"`.}
}"
rd_mdl_astmpl_ <- "`logical(1)`
Whether to store the compiled model as a local template.
If `FALSE`, the model is recompiled in every new session."
rd_mdl_force_ <- "`logical(1)`
Whether to compile even when a cached model exists,
for instance after editing a header under `paths`."
rd_mdl_quiet_ <- "`logical(1)` Whether to hide compilation commands."
rd_mdl_escaped_ <- "`character(1)` Label of status code `0`,
members outside the declared infected statuses,
when the `status` column is character or factor."
rd_mdl_flags_ <- "Named arguments. Additional compilation options:
\\describe{
  \\item{`discr`}{`logical(1)` Whether to compile to discrete time.}
  \\item{`flags`}{`character()` Compilation flags.}
  \\item{`libs`}{`character()` Links to libraries.}
  \\item{`paths`}{`character()` Include paths.}
}"

####  Definition vocabulary  ####

rd_var_name_ <- "`character(1)` Name used by internal calls."
rd_var_const_ <- "`numeric(1)` | `character(1)`
Value the definition is pinned to, in inference and simulation.
Either a number, or a data column, a constant parameter,
or an arithmetic expression of them.
Its prior, proposal and bounds are ignored.
Cannot be combined with `func`."
rd_aug_const_ <- paste(rd_var_const_, "
A latent with a `model` keeps its density at that value,
truncated by `min` and `max`.
A latent the transmission model generates is pinned in inference only:
simulation sets it.")
rd_var_discr_ <- function(what, from) {
  sprintf("`logical(1)` Whether the %s definition is discrete, if applicable.
          If `NA`, inferred from %s.", what, from)
}
rd_obs_discr_from_ <- "the observed data column, or the model requirements"
rd_var_prior_ <- "`character(1)` Prior distribution."
rd_var_args_ <- function(what, proposal = FALSE, func = FALSE) {
  sink <- if (func) ", or for `func` when the definition is deterministic."
  else "."
  paste0(sprintf("Named arguments for the %s distribution", what), sink,
         if (proposal) rd_proposal_opts_)
}
rd_var_min_ <- "`numeric(1)` Truncate to a minimum value.
Relative to the location argument `loc`, when one is given."
rd_var_max_ <- "`numeric(1)` Truncate to a maximum value.
Relative to the location argument `loc`, when one is given."
rd_var_func_ <- "Deterministic definition, one of:
* `character(1)` Name of a latent function, its arguments in `...`;
* a bare name, reading that definition as is;
* an rlang expression."
rd_var_proposal_ <- "`character(1)` Proposal used in the Monte Carlo.
Random walks update the current value:
\\describe{
  \\item{`\"rw\"`}{normal increments;}
  \\item{`\"rw_log\"`}{lognormal log-increments,
    positive non-discrete values only;}
  \\item{`\"rw_unif\"`}{uniform increments, bounded values only.}
}
A discrete variable takes its increments rounded to the nearest whole number.
Independence samplers draw a fresh value from a distribution:
\\describe{
  \\item{`\"forward\"`}{drawn from the variable's own defining distribution,
    its prior for a parameter, its model for a latent variable;}
  \\item{a distribution name}{drawn from that distribution,
    e.g. `proposal = \"normal\"` with `proposal_mean` and `proposal_sd`.}
}
If `NA`, normal increments of a coordinate set by the bounds of the variable:
the value itself when unbounded,
the log of its distance to the bound when bounded on one side,
the logit of its position when bounded on both sides.
The move never leaves those bounds.
A discrete variable takes normal increments of the value itself,
and a move past its bounds is refused.
The latent variables of one member with `NA` move together,
and so do the parameters a reference couples, see [run_inference()].
"
rd_proposal_opts_ <- "
Advanced Monte Carlo options may also be passed through `...`,
each prefixed `proposal_`.
For a random walk, `proposal_scale` sets the scale of its jumps.
If omitted, it is adapted during inference when `adapt > 0`.
If a positive number, it is fixed and never adapted.
For independence sampling from a named distribution,
the distribution parameters are passed the same way,
e.g. `proposal_mean` and `proposal_sd`."
rd_sources_ <- "Either a plain value, the name of another definition,
or an arithmetic expression combining them.
A name is looked up among the latent variables, data columns, observations,
group rates and parameters the model carries.
Two of them sharing a name cannot be read."
rd_aug_status_ <- "`integer(1)` | `character(1)`
Infected status the definition is specific to,
an integer code or a label, as in data column `status`.
If `NA`, the definition applies to every infected member."
rd_obs_func_ <- paste("Deterministic definition of the observation,
taking no `model`, depending on `level`:
\\describe{
  \\item{`\"member\"`}{the value computed for every member;}
  \\item{`\"house\"`}{one value for the whole household.
  The reductions `min()`, `max()`, `sum()` and `prod()` fold over the members
  a per-member value: a data column, a latent or a group rate.
  Their argument is computed for each member, members with no value left out,
  e.g. `min(onset_date) + 21` or `max(test_date - onset_date)`.
  A per-member value read outside a reduction is an error.}
}", rd_sources_)
rd_aug_min_ <- paste(rd_var_min_, rd_sources_)
rd_aug_max_ <- paste(rd_var_max_, rd_sources_)
rd_lat_flat_bound_ <- "Without a `model`,
the bound the latent is sampled within in inference,
reading only data, observations and constant parameters."

####  Definition modes  ####

rd_add_modes_ <- "There are three mutually exclusive definition modes:
\\describe{
  \\item{single}{
  Supply a `name` to add a single definition using the explicit arguments;
  }
  \\item{multiple}{
  Pass pairs `name = list()` (or the shorthand `name = x`)
  for multiple definitions;
  }
  \\item{import}{
  Pass a pre-built named list of definitions through `import`.
  }
}"
rd_dots_multi_ <- function(const = TRUE) {
  paste("In mode 'multiple', the `name = list()` definitions,",
        "or `name = x` for",
        if (const) "`const = x` when `x` is a number, `func = x` otherwise."
        else "`func = x`.")
}
rd_import_ <- "Named `list` of pre-built definitions to import.
Cannot be combined with single or multiple definitions."
rd_tmpl_import_ <- "Named `list` of pre-built definitions to import,
or `\"template\"` for the definitions shipped with the model.
The template is R code evaluated from the model folder,
so use it only with model folders you trust.
Cannot be combined with single or multiple definitions."
rd_list_append_ <- "`logical(1)` Whether to append to the current definitions.
If `FALSE`, the previous definitions are dropped."

####  Rates  ####

rd_rate_desc_ <- function(kind) {
  sprintf("Add definitions of a %s rate by given data covariate(s).
The rate multiplies a transmission rate of the model.
Each stratum carries its own absolute rate, unless a reference is chosen.
Then every parameter is a rate ratio to the reference stratum.
If no MCMC options are provided,
absolute rates are \\eqn{\\mathrm{uniform}} over their support,
and relative ones \\eqn{\\mathrm{lognormal}[\\mu=0, \\sigma=3]}.
Bypass automatic generation by providing a custom pre-built logic table.",
          kind)
}
rd_rate_prefix_ <- "`character(1)` Prefix to names of associated parameters."
rd_grr_prefix_ <- paste(
  rd_rate_prefix_,
  "If `NA`, the rate name."
)
rd_crr_prefix_ <- paste(rd_rate_prefix_, "If `NA`, `contact`.")
rd_rate_joint_ <- "`logical(1)` Whether covariates are joint.
If `FALSE`, they are independent."
rd_rate_strata_ <-
  "`data.frame` Individual covariates and relative `param_name`.
Each row defines one stratum.
A covariate left `NA` matches any value.
Every member, or ordered pair of members for contact rates, must match a row,
and every row must be the best match of at least one of them,
else the definition is refused.
So is a table setting no covariate value,
which treats every member alike as the baseline does.
The table is used as written.
A rate read on infected members is checked at every infected status
a member may take, whether or not the data records it.
Rows keyed on different covariates form separate rule groups,
whose rates multiply.
Use `.reference` in `param_name` to designate the reference group.
A rule group without it leaves the rate absolute.
The reserved numeric column `const` pins a stratum to that value,
in inference and simulation, so that it carries no prior.
The reference stratum takes none.
Mutually exclusive with `by`."
rd_rate_params_ <- "`list` Per-parameter MCMC configurations,
named as in `param_name`.
Overrides shared prior/args for parameters specified in the list.
A stratum pinned through the `strata` table takes no entry here.
Ignored in automatic mode."
rd_rate_reference_ <- function(cell) {
  sprintf("`list` Named covariates of the reference %s,
whose rate becomes the unit that every other parameter multiplies.
Pass `\".auto\"` to take the %s with the largest sample size.
Left empty, each stratum carries its own absolute rate.
A covariate the reference does not name leaves the whole rate absolute.",
          cell, cell)
}

####  Run-time knobs  ####

rd_smpl_steps_ <- "`integer(1)` Number of MCMC steps per chain."
rd_smpl_thin_ <- "`integer(1)` Keep one step every `thin` of `steps`."
rd_smpl_adapt_ <- "`integer(1)` Number of steps to adapt proposal scales
and the joint moves of latent variables and parameters.
Run first and not stored as draws.
If nothing is adaptable, run as burn-in instead."
rd_smpl_chain_ <- "`integer(1)` Number of independent MCMC chains to run."
rd_thread_ <- "`integer(1)` Number of threads,
over which households are spread.
Calls after the same `set.seed()` draw the same replicates
at any `threads`."
rd_progress_ <- "`numeric(1)` | `logical(1)` Control the progress bar,
one of:
\\describe{
  \\item{`numeric`}{seconds between prints, a shorter run printing nothing;}
  \\item{`TRUE`}{one print per second;}
  \\item{`FALSE` or `0`}{no progress bar.}
}
The bar prints as messages, which [suppressMessages()] silences.
In a knitted document, only its closing line prints."
rd_dots_config_ <- "Advanced configuration."
rd_draws_ <- "Parameter values, keyed by parameter name.
Either posterior draws, from [run_inference()] or as `posterior::draws`,
a `data.frame` with one row per draw,
or one value per parameter, as in `c(alpha = 0.1, beta = 0.2)`."
rd_nrep_ <- "`integer(1)` Number of replicates.
With `draws`, at most `nrep` evenly spaced draws are simulated,
picked after thinning to their ESS when `thin_ess` is `TRUE`.
If `NA`, every draw is simulated."
rd_nrep_prior_ <- "Without `draws`, `nrep` replicates are simulated,
each with parameters drawn from their prior."
rd_check_source_ <- function(table = TRUE) {
  paste0(
    "`character(1)` Where the parameter values come from, one of:
\\describe{
  \\item{`\"prior\"`}{drawn from the prior;}
  \\item{`\"posterior\"`}{posterior draws of the model on its attached data",
    if (table) ";}
  \\item{`\"table\"`}{values given by the user,
    or draws of [run_inference()] on another model or other data.}" else ".}",
    "\n}"
  )
}
rd_thin_ess_ <- "`logical(1)`
Whether to thin posterior `draws` to their bulk effective sample size (ESS)."
rd_natural_log_ <- "Logarithms are natural (base \\eqn{e})."
rd_global_summary_ <- "Results are summarised at the point estimate
and credible level of the global options, see [set_housepi_summary()]."
rd_warnings_ <- "`data.frame` Distinct warnings raised,
recorded whether or not they were reported, and counted by the print.
\\describe{
  \\item{message}{Line stating what was flagged;}
  \\item{count}{Times it was raised.}
}"

####  State inspection  ####

rd_inspect_state_ <- "The values a `HousepiModel` holds
are left by the last call that ran it:
\\describe{
  \\item{[evaluate_log_lik()]}{the evaluated parameters and latent variables,
  with their log-likelihood contributions;}
  \\item{[run_inference()]}{the final step of the last chain;}
  \\item{[run_simulation()]}{the observed data,
  the parameters of the last replicate.}
}"

# nocov end

########  Construction  ########

#' @title Create a new household-transmission model instance
#' @description
#' Create a `HousepiModel`, the object inference and simulation run on.
#' @param data `r rd_dataset_`
#' If omitted, attach it later with [add_dataset()].
#' @param model `r rd_mdl_choice_`
#' @param as_template `r rd_mdl_astmpl_`
#' @param force `r rd_mdl_force_`
#' @param quiet `r rd_mdl_quiet_`
#' @param escaped `r rd_mdl_escaped_`
#' @param ... `r rd_mdl_flags_`
#' @return `r rd_model_in_`
#' @examples
#' # build a model from a template using a built-in dataset
#' hh_model <- new_housepi_model(flu_1918_frederick, "final_size")
#' print(hh_model)
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

########  Data  ########

# Assert the object is a `HousepiModel`.
stopifnothousepi <- function(hh_model) {
  if (!inherits(hh_model, "HousepiModel")) stop_h(
    "{.arg hh_model} must be of class {.cls HousepiModel}.",
    "Create one with {.fn new_housepi_model}."
  )
}

#' @title Attach a dataset to a housepi model
#' @description
#' Make an internal copy of the household data, one row per member,
#' (re)building any definitions that depend on it.
#' An empty cell of a categorical column is missing data, not a category.
#' To stratify on it, recode it to a label of your choice.
#' @param hh_model `r rd_model_in_`
#' @param data `r rd_dataset_`
#' @return `r rd_model_out_`
#' @examples
#' # construct the model first, attach the data afterwards
#' hh_model <- new_housepi_model(model = "final_size") |>
#'   add_dataset(flu_1918_frederick)
#' summarise_house_cases(hh_model)
#' @export
add_dataset <- function(hh_model, data) {
  stopifnothousepi(hh_model)
  if (is.null(data)) stop_h("Missing argument {.arg data}.")
  new_obj <- hh_model$copy()
  pvt <- new_obj$.__enclos_env__$private
  pvt$set_dataset_(assert_household_data(data, NA_character_))
  pvt$rebuild_()
  new_obj
}

#' @title Attach a long dataset to a housepi model
#' @description
#' Make an internal copy of member data with any number of rows per member,
#' such as repeated tests over time, matched to the dataset on `id`.
#' @param hh_model `r rd_model_in_`
#' @param data `data.frame` Long dataset, one row per member and entry.
#' @return `r rd_model_out_`
#' @examples
#' # repeated measures of one individual, matched on `id`
#' measures <- data.frame(id = "hh001m02", day = 1:3,
#'                        viral_load = c(5.2, 6.8, 4.1))
#' hh_model <- sarscov2_2021_icofs |>
#'   new_housepi_model("case-ascertained") |>
#'   add_long_dataset(measures)
#' get_long_dataset(hh_model)
#' @export
add_long_dataset <- function(hh_model, data) {
  stopifnothousepi(hh_model)
  new_obj <- hh_model$copy()
  pvt <- new_obj$.__enclos_env__$private
  pvt$set_long_dataset_(data)
  pvt$rebuild_()
  new_obj
}

#' @title Summarise observed cases in the loaded household data
#' @description
#' Summarise cases from the statuses in the loaded household data.
#' The row `classes` counts the equivalence classes,
#' households the model cannot tell apart, e.g. same size, statuses and rates.
#' A sampled latent variable makes every household a class of its own.
#' @param hh_model `r rd_model_in_`
#' @return `data.frame` One row per metric:
#' \describe{
#'   \item{level}{`individual` or `household`;}
#'   \item{metric}{quantity counted in the row;}
#'   \item{count}{absolute count;}
#'   \item{percent}{percentage of the relevant total, `NA` for `classes`.}
#' }
#' @examples
#' hh_model <- new_housepi_model(flu_1918_frederick, "final_size")
#' summarise_house_cases(hh_model)
#' @export
summarise_house_cases <- function(hh_model) {
  stopifnothousepi(hh_model)
  pvt <- hh_model$.__enclos_env__$private
  .engine$summarise_house_cases(pvt$mcmc_ptr_)
}

#' @rdname summarise_house_cases
#' @export
summarize_house_cases <- summarise_house_cases

########  Getters  ########

# Keep the definitions whose column `col` matches one of `name`.
# * `name` NA -> all of them.
subset_definitions_ <- function(defs, name, col = "name") {
  if (length(name) == 1L && is.na(name)) return(defs)
  stopifnot(is.character(name), !anyNA(name))
  unmatched <- setdiff(name, defs[[col]])
  if (length(unmatched)) stop_h(
    "No definition{?s} named {.val {unmatched}}.",
    find_closest_match(unmatched, defs[[col]],
                       otherwise = "Omit {.arg name} to show every definition.")
  )
  out <- defs[defs[[col]] %in% name, , drop = FALSE]
  rownames(out) <- NULL
  out
}

#' @title Retrieve currently loaded data and definitions from a housepi model
#' @description
#' Return the definitions currently stored in the housepi model.
#' Each variant returns one category: dataset, long dataset, parameters,
#' latent variables, observations, status process, group rates,
#' contact rates, infectivity profile or reporting process.
#' @details
#' `get_dataset()` returns the dataset the model reads, one row per member.
#' A column a deterministic observation stands in for
#' carries the computed value, not the imported one.
#' A column an observation produces without standing in for data
#' is not part of the dataset, and is reported by [get_observations()].
#' @param hh_model `r rd_model_in_`
#' @param name `character()` Names of the definitions to return.
#' `get_group_rates()` matches the name of the group rate.
#' If `NA`, all of them.
#' @return `data.frame`
#' `get_status_process()` returns one row per rule of the status process:
#' its `status` and `weight`, the `n` members it applies to,
#' its `rule_group` and its covariates.
#' `get_group_rates()` returns one row per rule of each group rate:
#' its `param_name`, the `n` members it applies to,
#' its `rule_group` and its covariates.
#' `get_contact_rates()` returns one row per rule of the contact rate:
#' its `from_` and `to_` covariates, its `param_name`,
#' the `n` ordered pairs of members it applies to and its `rule_group`.
#' A symmetric contact rate lists the reverse of each rule as a rule of its own.
#' `get_reporting_process()` returns `list(dataset = , long_dataset = )`,
#' see [add_reporting_process()].
#' @examples
#' hh_model <- new_housepi_model(flu_1918_frederick, "final_size") |>
#'   add_parameter(import = "template")
#' get_dataset(hh_model) |>
#'   head()
#' get_parameters(hh_model)
#' get_parameters(hh_model, name = c("alpha", "beta"))
#' @rdname get_definitions
#' @export
get_dataset <- function(hh_model) {
  stopifnothousepi(hh_model)
  pvt <- hh_model$.__enclos_env__$private
  # no data attached yet (deferred until add_dataset() or the wall)
  if (!pvt$has_dataset_()) return(NULL)
  .engine$export_dataset(pvt$mcmc_ptr_) |>
    pvt$refactor_data_("dataset")
}

#' @rdname get_definitions
#' @export
get_long_dataset <- function(hh_model) {
  stopifnothousepi(hh_model)
  pvt <- hh_model$.__enclos_env__$private
  .engine$export_long_dataset(pvt$mcmc_ptr_) |>
    pvt$refactor_data_("long_dataset")
}

#' @rdname get_definitions
#' @export
get_parameters <- function(hh_model, name = NA) {
  stopifnothousepi(hh_model)
  hh_model$.__enclos_env__$private$resolved_param_() |>
    housepi_defin_showcaser() |>
    sort_definitions_() |>
    subset_definitions_(name)
}

#' @rdname get_definitions
#' @export
get_latent <- function(hh_model, name = NA) {
  stopifnothousepi(hh_model)
  hh_model$.__enclos_env__$private$resolved_latent_() |>
    augm_defs_of_kind("latent") |>
    housepi_defin_showcaser() |>
    sort_definitions_() |>
    subset_definitions_(name)
}

#' @rdname get_definitions
#' @export
get_observations <- function(hh_model, name = NA) {
  stopifnothousepi(hh_model)
  resolved <- hh_model$.__enclos_env__$private$resolved_latent_()
  do.call(rbind, lapply(c("member", "house"), function(level) {
    out <- augm_defs_of_kind(resolved, level)
    if (is.null(out)) return(NULL)
    # pinned to its datum, the like-named column
    out[["const"]] <- NULL
    cbind(out, level = level)
  })) |>
    housepi_defin_showcaser() |>
    sort_definitions_() |>
    subset_definitions_(name)
}

#' @rdname get_definitions
#' @export
get_status_process <- function(hh_model) {
  stopifnothousepi(hh_model)
  pvt <- hh_model$.__enclos_env__$private
  dataset <- get_dataset(hh_model)
  levels <- if (!is.null(dataset)) pvt$status_draw_levels_()
  pvt$status_$tidy_table(dataset, levels)
}

#' @rdname get_definitions
#' @export
get_group_rates <- function(hh_model, name = NA) {
  stopifnothousepi(hh_model)
  pvt <- hh_model$.__enclos_env__$private
  pvt$group_rates_ |>
    helper_show_rate_table(setdiff(rate_reserved_cols, rate_pin_cols),
                           pvt$has_dataset_()) |>
    sort_definitions_("group_rate", by_status = FALSE) |>
    subset_definitions_(name, "group_rate")
}

#' @rdname get_definitions
#' @export
get_contact_rates <- function(hh_model) {
  stopifnothousepi(hh_model)
  pvt <- hh_model$.__enclos_env__$private
  pvt$contact_rates_ |>
    helper_show_rate_table(c("param_name", "n", "rule_group"),
                           pvt$has_dataset_()) |>
    sort_definitions_("param_name", by_status = FALSE)
}

#' @rdname get_definitions
#' @export
get_infectivity_profiles <- function(hh_model) {
  stopifnothousepi(hh_model)
  pvt <- hh_model$.__enclos_env__$private
  pvt$resolved_profile_() |>
    pvt$infec_profile_$showcase() |>
    sort_definitions_()
}

#' @rdname get_definitions
#' @export
get_reporting_process <- function(hh_model) {
  stopifnothousepi(hh_model)
  hh_model$.__enclos_env__$private$reporting_process_
}

########  Definitions  ########

definition_helper_ <- c(
  "Supply a {.arg {key_arg}} for a single definition.",
  "Provide multiple definitions as {.code name = list()} pairs.",
  "Otherwise import an external list using {.arg import}."
)

# Extra hint for the families shipping a definition template.
template_helper_ <- paste(
  "Load the definitions shipped with the model with",
  '{.code import = "template"}.'
)

# Stop on missing definitions, followed by the ways to supply them.
# * `single` -> the single definition misses its key `key_arg`;
# * otherwise -> no definition at all.
stop_missing_defs_ <- function(family, key_arg, single, has_template = FALSE) {
  header <- "Missing or empty {family} definitions."
  if (single)
    header <- "Missing {.arg {key_arg}} of the single {family} definition."
  stop_h(header, definition_helper_, if (has_template) template_helper_)
}

# Name the formals given in an `add_*` call that only single mode reads.
# @param listed `character()` Formals the list-like modes read,
# besides `...`, `import` and `append`.
single_args_ <- function(call, listed) {
  setdiff(names(call)[-1], c("hh_model", "...", "import", "append", listed))
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
assert_no_def_mixing_ <- function(key, dots, singles, key_arg) {
  if (!is.null(key) || length(dots) || length(singles)) stop_h(
    "Mixing imported definitions with single and multiple definitions",
    definition_helper_
  )
}

# Whether `import` is the sentinel selecting the definitions the model ships.
# * any other character value -> error, it names no existing source.
is_template_import_ <- function(import) {
  if (!is.character(import)) return(FALSE)
  if (!identical(import, "template")) stop_h(
    "Unknown import source {.val {import}}.",
    template_helper_
  )
  TRUE
}

# Read the templated definitions shipped with the model.
# * test model -> `NULL` with a warning, it ships none.
template_defs_ <- function(compiled, family) {
  if (is_test_model_(compiled$name)) {
    warn_h("No template {family} definitions available for test model.")
    return(NULL)
  }
  get_template_def_list(compiled$source_dir, family)
}

# Load the templated definitions shipped with the package.
add_template_defs_ <- function(holder, compiled, family, append) {
  defs <- template_defs_(compiled, family)
  if (!is.null(defs))
    holder$add_quo_list(convert_list_to_quosures(defs), append = isTRUE(append))
}

# Add the definitions of the list-like `add_*` modes.
# * `import = "template"`, an `import` list, inline `...` pairs -> added
#   through `$add_quo_list()`;
# * `import` with single-mode arguments or inline `...` pairs -> error;
# * single mode, `key` `NA` -> error, as when missing;
# * single mode, `key` given -> `FALSE`, left to the caller,
#   which captures its own arguments, the definitions dropped first
#   when not `append`.
# @param call Caller's `match.call(expand.dots = FALSE)`.
# @param listed `character()` Formals these modes read, the key first,
# besides `...`, `import` and `append`.
# Any other formal the call sets belongs to single mode.
add_listed_defs_ <- function(holder, compiled, family, key, import, append,
                             dots, env, call, listed,
                             has_template = FALSE, ...) {
  singles <- single_args_(call, listed)
  template <- has_template && is_template_import_(import)
  if (length(import)) assert_no_def_mixing_(key, dots, singles, listed[1])
  if (template) {
    add_template_defs_(holder, compiled, family, append)
  } else if (length(import)) {
    holder$add_quo_list(convert_list_to_quosures(import, env),
                        append = append, ...)
  } else if (is.null(key)) {
    if (!length(dots))
      stop_missing_defs_(family, listed[1], length(singles) > 0, has_template)
    assert_no_reserved_defs_(singles)
    holder$add_quo_list(dots, append = append, ...)
  } else {
    if (isTRUE(is.na(key)))
      stop_missing_defs_(family, listed[1], single = TRUE, has_template)
    if (isFALSE(append)) holder$add_quo_list(list(), append = FALSE, ...)
    return(FALSE)
  }
  TRUE
}

#' @title Add or overwrite parameter definitions
#' @description
#' Add parameters to a housepi model.
#' A parameter with the same name as an existing one overwrites it.
#' `r rd_add_modes_`
#' @param hh_model `r rd_model_in_`
#' @param name `r rd_var_name_`
#' @param const `r rd_var_const_`
#' @param discr `r rd_var_discr_("prior", "the model requirements")`
#' @param prior `r rd_var_prior_`
#' @param ... `r rd_var_args_("prior", proposal = TRUE, func = TRUE)`
#' `r rd_dots_multi_()`
#' @param min `r rd_var_min_`
#' @param func `r rd_var_func_`
#' @param max `r rd_var_max_`
#' @param proposal `r rd_var_proposal_`
#' @param import `r rd_tmpl_import_`
#' @param append `r rd_list_append_`
#' @return `r rd_model_out_`
#' @examples
#' # the single, inline-bulk and list-bulk forms compose freely in a pipe
#' hh_model <- flu_1918_frederick |>
#'   new_housepi_model("final-size") |>
#'   # single: alpha, with a uniform prior and a lognormal random-walk proposal
#'   add_parameter("alpha", prior = "uniform", proposal = "rw_log") |>
#'   # single: beta, fixed to a constant
#'   add_parameter("beta", const = 0.3) |>
#'   # inline bulk: several constants as `name = x` pairs
#'   add_parameter(gamma = 1, hh_size_ref = 4)
#'
#' # list bulk: load a pre-built named list of definitions through `import`
#' my_import <- list(gamma = 1, hh_size_ref = 4)
#' hh_model <- flu_1918_frederick |>
#'   new_housepi_model("final-size") |>
#'   add_parameter(import = my_import)
#'
#' # `"template"` loads the definitions shipped with the model
#' hh_model <- flu_1918_frederick |>
#'   new_housepi_model("final-size") |>
#'   add_parameter(import = "template")
#' get_parameters(hh_model)
#' @export
add_parameter <- function(hh_model,
                          name = NULL, ...,
                          const = NULL, discr = NA,
                          prior = NA_character_,
                          func = NULL,
                          min = NULL, max = NULL,
                          proposal = NA_character_,
                          import = NULL, append = TRUE) {
  stopifnothousepi(hh_model)
  new_obj <- hh_model$copy()
  pvt <- new_obj$.__enclos_env__$private

  #> template, import and multiple modes
  if (!add_listed_defs_(pvt$param_, pvt$compiled_, "parameter", name, import,
                        append, rlang::enquos(...), rlang::caller_env(),
                        match.call(expand.dots = FALSE), "name",
                        has_template = TRUE)) {
    #> single mode
    pvt$param_$add_quo(name = name,
                       const = rlang::enquo(const), discr = discr,
                       prior = prior, rlang::enquos(...),
                       func = rlang::enquo(func),
                       min = rlang::enquo(min), max = rlang::enquo(max),
                       proposal = proposal,
                       .assert = TRUE)
  }
  pvt$rebuild_()
  new_obj
}

########  Augmented  ########

#' @title Add or update latent variable definitions
#' @description
#' Add latent variables to a housepi model.
#' A latent with the same name and status as an existing one overwrites it.
#' `r rd_add_modes_`
#' @details
#' A latent with a `model` is sampled from it.
#' A latent the transmission model generates takes no `model`,
#' such as the infection times of [template_case_ascertained].
#' @param hh_model `r rd_model_in_`
#' @param name `r rd_var_name_`
#' @param status `r rd_aug_status_`
#' @param const `r rd_aug_const_`
#' @param discr `r rd_var_discr_("model", "the model requirements")`
#' @param func `r rd_var_func_`
#' @param model `character(1)` Latent model distribution.
#' Omitted for a latent the transmission model generates (see Details).
#' @param ... `r rd_var_args_("model", proposal = TRUE, func = TRUE)`
#' `r rd_dots_multi_()`
#' @param min `r rd_aug_min_` `r rd_lat_flat_bound_`
#' @param max `r rd_aug_max_` `r rd_lat_flat_bound_`
#' @param proposal `r rd_var_proposal_`
#' @param import `r rd_tmpl_import_`
#' @param append `r rd_list_append_`
#' @return `r rd_model_out_`
#' @examples
#' hh_model <- dummy_household |>
#'   new_housepi_model("test") |>
#'   # define an incubation period drawn from a gamma distribution
#'   add_latent("incub_period",
#'                 model = "gamma", shape = 2, scale = 1, max = 14)
#' @export
add_latent <- function(hh_model,
                       name = NULL, ...,
                       status = NA_integer_, const = NULL, discr = NA,
                       func = NULL,
                       model = NA_character_,
                       min = NULL, max = NULL,
                       proposal = NA_character_,
                       import = NULL, append = TRUE) {
  stopifnothousepi(hh_model)
  new_obj <- hh_model$copy()
  pvt <- new_obj$.__enclos_env__$private

  #> template, import and multiple modes
  if (!add_listed_defs_(pvt$augm_, pvt$compiled_, "latent", name, import,
                        append, rlang::enquos(...), rlang::caller_env(),
                        match.call(expand.dots = FALSE), c("name", "status"),
                        has_template = TRUE, status = status)) {
    #> single mode
    pvt$augm_$add_quo(name = name,
                      status = status, const = rlang::enquo(const),
                      discr = discr,
                      func = rlang::enquo(func),
                      model = model, rlang::enquos(...),
                      min = rlang::enquo(min), max = rlang::enquo(max),
                      proposal = proposal,
                      .assert = TRUE)
  }
  pvt$rebuild_()
  new_obj
}

#' @title Add an observation
#' @description
#' Add an observation of data \eqn{x}, at the level set by `level`.
#'
#' At member level, the observation is attached to every infected individual.
#' If `model` is provided, it adds to the household likelihood in evaluations,
#' and \eqn{x} is drawn from it in simulations.
#' An individual whose datum is missing is left out of that likelihood,
#' with a warning.
#' Nothing conditions on it, so the observation constrains neither
#' the individual latent variables nor the parameters.
#' A column holding no value at all is an error.
#'
#' At household level, the observation is a reduction across all members,
#' shared by the whole household, and must reduce to a constant value.
#'
#' An observation carrying no `model` produces its own data column:
#' it is evaluated for every member, and overwrites data of the same name.
#'
#' `r rd_add_modes_`
#' @param hh_model `r rd_model_in_`
#' @param name `r rd_var_name_`
#' @param model `character(1)` Observation model distribution.
#' Reserved to member-level observations.
#' @param ... `r rd_var_args_("observation")`
#' `r rd_dots_multi_(const = FALSE)`
#' @param level `character(1)` Level the datum is observed at, one of:
#' \describe{
#'   \item{`"member"`}{one value per infected member;}
#'   \item{`"house"`}{one value per household, reduced over its members.}
#' }
#' In mode 'multiple', applies to every entry not setting its own.
#' @param func `r rd_obs_func_`
#' @param status `r rd_aug_status_`
#' @param discr `r rd_var_discr_("model", rd_obs_discr_from_)`
#' @param min `r rd_aug_min_`
#' @param max `r rd_aug_max_`
#' @param import `r rd_import_`
#' @param append `r rd_list_append_`
#' @return `r rd_model_out_`
#' @examples
#' hh_model <- sarscov2_2021_icofs |>
#'   new_housepi_model("case-ascertained") |>
#'   add_parameter(import = "template") |>
#'   add_infectivity_profile(import = "template") |>
#'   add_latent(
#'     infec_time = list(min = end_date - 31, max = end_date),
#'     incub_period = list(model = "gamma", mean = 3, sd = 1.4),
#'     onset_time = infec_time + incub_period
#'   ) |>
#'   # member level: the test date is noisy around the onset time
#'   add_observation("test_date", model = "normal",
#'                   mean = onset_time, sd = 1, min = floor(onset_time)) |>
#'   # house level: the earliest positive test among members
#'   add_observation("first_test", func = min(test_date), level = "house")
#' get_observations(hh_model)
#' @export
add_observation <- function(hh_model, name = NULL,
                            model = NA_character_, ...,
                            level = c("member", "house"),
                            func = NULL, status = NA_integer_, discr = NA,
                            min = NULL, max = NULL,
                            import = NULL, append = TRUE) {
  stopifnothousepi(hh_model)
  # bulk default: a per-definition 'level =' overrides it (see add_quo)
  level <- match_choice_(level, c("member", "house"), "level",
                         "add_observation")
  new_obj <- hh_model$copy()
  pvt <- new_obj$.__enclos_env__$private

  #> import and multiple modes
  if (!add_listed_defs_(pvt$augm_, pvt$compiled_, kind2family(level), name,
                        import, append, rlang::enquos(...),
                        rlang::caller_env(), match.call(expand.dots = FALSE),
                        c("name", "status", "level"),
                        status = status, .kind = level)) {
    #> single mode
    pvt$augm_$add_quo(name = name,
                      status = status, discr = discr,
                      func = rlang::enquo(func), model = model,
                      rlang::enquos(...),
                      min = rlang::enquo(min), max = rlang::enquo(max),
                      .assert = TRUE, .kind = level)
  }
  pvt$rebuild_()
  new_obj
}

#' @title Add a status process
#' @description
#' Declare how an infected member's status is drawn during simulation,
#' as a categorical distribution over infected statuses.
#' The probabilities apply to every member when given in `...`,
#' or depend on member covariates when given as a `strata` table.
#' Inference holds the observed status fixed,
#' and does not estimate the probabilities.
#' Statuses are integer codes or labels, as in data column `status`.
#' Data used only as household structure for a simulation can drop the column.
#' @param hh_model `r rd_model_in_`
#' @param ... Probabilities of the infected statuses, in order.
#' When named, the names label the statuses: name all of them or none.
#' Unnamed probabilities draw the integer statuses 1, 2, and so on.
#' @param strata `data.frame` One row per rule:
#' the infected status in column `status`, its non-negative `weight`,
#' and a column per covariate.
#' An `NA` cell matches any value,
#' and a member reads the most specific rule it matches.
#' Rows keyed on different covariates form separate rule groups,
#' whose weights multiply.
#' Every rule group must key on `status`.
#' At every status, each member must match a rule of every rule group,
#' and every rule must apply to some member, else the table is refused.
#' The weights of a member are normalized over the statuses,
#' so they need not sum to 1.
#' `status` holds labels, or integer statuses from 1.
#' Cannot be combined with `...`.
#' @param import Named `list` of arguments, read as `...` or a `strata` entry,
#' or `"template"` for the status process shipped with the model.
#' The template is R code evaluated from the model folder,
#' so use it only with model folders you trust.
#' Cannot be combined with `...` or `strata`.
#' @return `r rd_model_out_`
#' @examples
#' # household structure only: the recorded status codes carry no label
#' households <- subset(flu_1918_frederick, select = -status)
#'
#' # simulate infected members as symptomatic (70%) or asymptomatic (30%)
#' hh_model <- households |>
#'   new_housepi_model("final-size") |>
#'   add_parameter(import = "template") |>
#'   add_status_process(symptomatic = 0.7, asymptomatic = 0.3)
#' get_status_process(hh_model)
#' simulated <- run_simulation(hh_model)
#' table(simulated$status)
#'
#' # the same probabilities, from a pre-built list
#' probs <- list(symptomatic = 0.7, asymptomatic = 0.3)
#' hh_model <- households |>
#'   new_housepi_model("final-size") |>
#'   add_status_process(import = probs)
#'
#' # adults are symptomatic more often than children
#' strata <- data.frame(
#'   age = c("adult", "adult", "child", "child"),
#'   status = c("symptomatic", "asymptomatic", "symptomatic", "asymptomatic"),
#'   weight = c(0.8, 0.2, 0.4, 0.6)
#' )
#' hh_model <- subset(sarscov2_2021_icofs, select = -status) |>
#'   new_housepi_model("final-size") |>
#'   add_status_process(strata = strata)
#' get_status_process(hh_model)
#' @export
add_status_process <- function(hh_model, ..., strata = NULL, import = NULL) {
  stopifnothousepi(hh_model)
  new_obj <- hh_model$copy()
  pvt <- new_obj$.__enclos_env__$private
  weights <- rlang::enquos(...)
  env <- rlang::caller_env()
  if (length(import) && (length(weights) || !is.null(strata))) stop_h(
    "Mixing imported and inline status probabilities.",
    "Provide them in {.arg ...}, in {.arg strata} or in {.arg import}."
  )
  if (is_template_import_(import)) {
    import <- template_defs_(pvt$compiled_, "status process")
    if (is.null(import)) return(new_obj)
  }
  # an imported list holds the arguments
  if (length(import)) {
    strata <- import[["strata"]]
    import <- import[rlang::names2(import) != "strata"]
    weights <- rlang::as_quosures(import, env = env)
  }
  if (!is.null(strata) && length(weights)) stop_h(
    "Mixing a status table with inline status probabilities.",
    "Provide them in {.arg ...} or in {.arg strata}."
  )
  if (is.null(strata)) strata <- status_strata_from_weights_(weights)
  # R keeps the rule table, push_status_process_() sends it coded to the engine
  pvt$status_$set(strata)
  pvt$rebuild_()
  new_obj
}

########  Rates  ########

#' @title Add a stratified group rate definition
#' @description `r rd_rate_desc_("group")`
#' `r rd_add_modes_`
#' @param hh_model `r rd_model_in_`
#' @param name `r rd_var_name_`
#' @param by `list` Data columns of the covariates the group rate stratifies on.
#' The columns must be categorical.
#' @param reference `r rd_rate_reference_("group")`
#' @param joint `r rd_rate_joint_`
#' @param strata `r rd_rate_strata_`
#' @param params `r rd_rate_params_`
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
#' hh_model <- new_housepi_model(sarscov2_2021_icofs, "case_ascertained")
#'
#' # single: one susceptibility parameter per age class
#' hh_model |>
#'   add_group_rate("susceptibility", by = "age") |>
#'   get_group_rates()
#'
#' # inline bulk: rate ratios, against a chosen or the most frequent group
#' hh_model |>
#'   add_group_rate(
#'     infectivity = list(by = "isolated", reference = ".auto"),
#'     susceptibility = list(by = "vacc_status", reference = "unvaccinated")
#'   ) |>
#'   get_parameters()
#' @export
add_group_rate <-  function(hh_model,
                            name = NULL, by = NA_character_, ...,
                            reference = NULL, joint = FALSE,
                            strata = NULL, params = NULL,
                            prior = NA_character_,
                            min = NA_real_, max = NA_real_,
                            proposal = NA_character_,
                            prefix = NA_character_,
                            import = NULL, append = TRUE) {
  stopifnothousepi(hh_model)
  new_obj <- hh_model$copy()
  pvt <- new_obj$.__enclos_env__$private

  # drop stale parameters and re-derive in `rebuild_()`
  pvt$drop_stale_rate_params_()
  if (is.null(name) && isFALSE(append)) pvt$drop_group_rates_()

  #> import and multiple modes
  if (!add_listed_defs_(pvt$group_rates_, pvt$compiled_,
                        pvt$group_rates_$family, name, import, append,
                        rlang::enquos(...), rlang::caller_env(),
                        match.call(expand.dots = FALSE), "name")) {
    #> single mode
    pvt$group_rates_$add_quo(name = name, by = by,
                             reference = reference, joint = joint,
                             strata = strata, params = params,
                             prior = prior, rlang::enquos(...),
                             min = rlang::enquo(min), max = rlang::enquo(max),
                             proposal = proposal, prefix = prefix)
  }
  pvt$rebuild_()
  new_obj
}

#' @title Add a pairwise contact rate definition
#' @description `r rd_rate_desc_("contact")`
#' @param hh_model `r rd_model_in_`
#' @param by `character()` Data columns of the covariates pairs are matched on.
#' The columns must be categorical.
#' @param reference `r rd_rate_reference_("pair")`
#' Names can carry `from_`/`to_` prefixes matching `by`,
#' or no prefix, applying to both ends.
#' @param joint `r rd_rate_joint_`
#' @param symmetric `logical(1)` Whether contact is symmetric.
#' If `TRUE`, infector-infectee and infectee-infector pairs are equivalent.
#' @param strata `r rd_rate_strata_`
#' @param params `r rd_rate_params_`
#' @param prior `r rd_var_prior_`
#' @param ... `r rd_var_args_("prior", proposal = TRUE)`
#' @param min `r rd_var_min_`
#' @param max `r rd_var_max_`
#' @param proposal `r rd_var_proposal_`
#' @param prefix `r rd_crr_prefix_`
#' @param import Named `list` of pre-built definition to import.
#' Cannot be combined with the single definition.
#' @return `r rd_model_out_`
#' @examples
#' hh_model <- new_housepi_model(flu_2014_matsumoto, "final_size")
#'
#' # one parameter per pair of family roles, relative to the most frequent pair
#' hh_model |>
#'   add_contact_rate("member", reference = ".auto") |>
#'   get_contact_rates()
#'
#' # asymmetric: infector-infectee and infectee-infector pairs differ
#' hh_model |>
#'   add_contact_rate("member", symmetric = FALSE) |>
#'   get_contact_rates()
#' @export
add_contact_rate <- function(hh_model, by = NA_character_, ...,
                             reference = NULL,
                             joint = FALSE, symmetric = TRUE,
                             strata = NULL, params = NULL,
                             prior = NA_character_,
                             min = NA_real_, max = NA_real_,
                             proposal = NA_character_,
                             prefix = NA_character_,
                             import = NULL) {
  stopifnothousepi(hh_model)
  new_obj <- hh_model$copy()
  pvt <- new_obj$.__enclos_env__$private
  pvt$assert_feature_supported_(.engine$has_contact_rate(pvt$mcmc_ptr_),
                                "Contact rates")
  # the definition is replaced: its parameters are re-derived by `rebuild_()`
  pvt$drop_stale_rate_params_()

  # single definition, from an imported named list or the explicit arguments
  if (!length(import)) import <- list(
    by = by, reference = reference, joint = joint,
    symmetric = symmetric, strata = strata, params = params,
    prior = prior, rlang::enquos(...),
    min = rlang::enquo(min), max = rlang::enquo(max),
    proposal = proposal, prefix = prefix
  )
  do.call(pvt$contact_rates_$add_quo, import)
  pvt$rebuild_()
  new_obj
}

########  Processes  ########

#' @title Add an infectivity profile definition
#' @description
#' Add infectivity profile definitions to a housepi model.
#' `r rd_add_modes_`
#' @param hh_model `r rd_model_in_`
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
#' hh_model <- sarscov2_2021_icofs |>
#'   new_housepi_model("case-ascertained") |>
#'   add_infectivity_profile(import = "template")
#' get_infectivity_profiles(hh_model)
#' @export
add_infectivity_profile <- function(hh_model,
                                    func = NULL, ...,
                                    status = NA_integer_, discr = NA,
                                    min = NULL, max = NULL,
                                    import = NULL, append = TRUE) {
  stopifnothousepi(hh_model)
  new_obj <- hh_model$copy()
  pvt <- new_obj$.__enclos_env__$private
  pvt$assert_feature_supported_(.engine$has_infec_profile(pvt$mcmc_ptr_),
                                "Infectivity profiles")
  # discrete-time model -> discrete profile
  pvt$infec_profile_$discr_model <- pvt$is_discrete_time_()

  #> template, import and multiple modes
  if (!add_listed_defs_(pvt$infec_profile_, pvt$compiled_,
                        "infectivity profile", func, import, append,
                        rlang::enquos(...), rlang::caller_env(),
                        match.call(expand.dots = FALSE), "func",
                        has_template = TRUE)) {
    #> single mode
    pvt$infec_profile_$add_quo(func = func, status = status, discr = discr,
                               rlang::enquos(...),
                               min = rlang::enquo(min), max = rlang::enquo(max),
                               .assert = TRUE)
  }
  pvt$rebuild_()
  new_obj
}

#' @title Add a reporting process to a housepi model
#' @description
#' Attach the functions turning each simulated household study
#' into the data reported.
#' [run_simulation()], [check_predictive()], [check_performance()] and
#' [check_calibration()] apply them to every replicate.
#' @details
#' Inference does not model the reporting process,
#' which introduces a misspecification.
#' @param hh_model `r rd_model_in_`
#' @param dataset `function` | `NULL`
#' Turns a simulated dataset into the reported dataset.
#' `NULL` removes it.
#' @param long_dataset `function` | `NULL`
#' Turns the reported dataset into a long dataset,
#' which [run_simulation()] returns along with it.
#' `NULL` removes it.
#' @return `r rd_model_out_`
#' @examples
#' # each case goes unreported with probability 0.1
#' miss_cases <- function(data) {
#'   missed <- data$status > 0 & stats::runif(nrow(data)) < 0.1
#'   data$status[missed] <- 0L
#'   data
#' }
#' flu_1918_frederick |>
#'   new_housepi_model("final-size") |>
#'   add_parameter(import = "template") |>
#'   add_reporting_process(dataset = miss_cases) |>
#'   run_simulation() |>
#'   head()
#' @export
add_reporting_process <- function(hh_model, dataset, long_dataset) {
  stopifnothousepi(hh_model)
  if (missing(dataset) && missing(long_dataset)) stop_h(
    "No reporting process given.",
    "Pass a function to {.arg dataset}, {.arg long_dataset} or both."
  )
  new_obj <- hh_model$copy()
  pvt <- new_obj$.__enclos_env__$private
  if (!missing(dataset)) pvt$set_reporting_process_("dataset", dataset)
  if (!missing(long_dataset))
    pvt$set_reporting_process_("long_dataset", long_dataset)
  new_obj
}

########  Evaluation  ########

#' @title Evaluate the log-likelihood of a housepi model at current values
#' @description
#' Compute the complete-data log-likelihood of a `HousepiModel`,
#' transmission plus observation.
#' With latent variables, it is conditioned on their current values.
#' @param hh_model `r rd_model_in_`
#' @param total `logical(1)` If `TRUE`, the sum over households,
#' else one value per household.
#' @return `numeric(1)` The complete-data log-likelihood.
#' With `total = FALSE`, a `numeric` vector named by household ID.
#' `r rd_natural_log_`
#' @examples
#' flu_1918_frederick |>
#'   new_housepi_model("final-size") |>
#'   add_parameter(import = "template") |>
#'   evaluate_log_lik()
#' @export
evaluate_log_lik <- function(hh_model, total = TRUE) {
  stopifnothousepi(hh_model)
  pvt <- hh_model$.__enclos_env__$private
  pvt$build_full_()
  house <- .engine$evaluate_log_lik(pvt$mcmc_ptr_)
  if (total) sum(house) else house
}

#' @title Inspect the augmented data of a housepi model
#' @description
#' Read the dataset of a `HousepiModel` with its latent variables.
#' `r rd_inspect_state_`
#' @param hh_model `r rd_model_in_`
#' @return `data.frame` The dataset with its latent variables.
#' @seealso [inspect_state()]
#' @examples
#' hh_model <- flu_1918_frederick |>
#'   new_housepi_model("final-size") |>
#'   add_parameter(import = "template")
#' evaluate_log_lik(hh_model)
#' inspect_augmented_data(hh_model) |>
#'   head()
#' @export
inspect_augmented_data <- function(hh_model) {
  stopifnothousepi(hh_model)
  pvt <- hh_model$.__enclos_env__$private
  pvt$refactor_data_(.engine$export_augmented_dataset(pvt$mcmc_ptr_),
                     "dataset")
}

#' @title Inspect the state of a housepi model
#' @description
#' Read the parameters a `HousepiModel` holds and how its variables move,
#' and optionally all the variables of one household.
#' `r rd_inspect_state_`
#' @param hh_model `r rd_model_in_`
#' @param hh_id `character(1)` ID of the household of interest.
#' If `NULL`, only the parameters and the sampler are read.
#' @return `list` The parameters and the sampler,
#' and with a household all its variables.
#' Entries are named after their getters, e.g `dataset` for [get_dataset()].
#' Then come the latent, observation and transmission log-likelihoods,
#' `model_quantities`, the named quantities the model exports,
#' see its `template_<model>` help page, and `.sampler`.
#' An empty definition is dropped unless the model requires it.
#' The transmission log-likelihood is `NA` and `model_quantities` absent
#' until the households are evaluated at the current definitions,
#' by [evaluate_log_lik()] or [run_inference()].
#' Each evaluation discards the previous one.
#'
#' `.sampler` reads how the variables move, absent while nothing is sampled.
#' Scales are the defaults until inference adapts them,
#' and a block appears once inference has formed it.
#' It holds two data frames:
#' \describe{
#'   \item{`moves`}{one row per sampled variable, the parameters
#'   and the latent variables of the infected members of the household.
#'   `member` is the member a latent variable belongs to,
#'   `NA` for a parameter.
#'   `proposal` is the declared proposal, `NA` for the default.
#'   `coord` is the coordinate its jumps are drawn in.
#'   `scale` is the scale of the jumps in that coordinate,
#'   the standard deviation of the jump of a variable moving in a block.
#'   For a discrete variable, it is the scale before the jump is rounded
#'   to an integer.
#'   `block` identifies the block the variable moves in,
#'   `NA` if it moves alone;}
#'   \item{`blocks`}{one row per block of variables moving together.
#'   `member` is the member whose latent variables move in it,
#'   `NA` for a block of parameters.
#'   `jump_cov` is the covariance of its joint jump in the coordinates.
#'   `log_scale` is the log of the multiplier of its jumps.}
#' }
#' Values are those of the last chain.
#' @seealso [inspect_augmented_data()]
#' @examples
#' hh_model <- flu_1918_frederick |>
#'   new_housepi_model("final-size") |>
#'   add_parameter(import = "template")
#' evaluate_log_lik(hh_model)
#' inspect_state(hh_model, "1")
#' @export
inspect_state <- function(hh_model, hh_id = NULL) {
  stopifnothousepi(hh_model)
  if (!is.null(hh_id) && (length(hh_id) != 1L || is.na(hh_id))) stop_h(
    "{.arg hh_id} is not a single household ID.",
    "Pass one ID, or {.code NULL} to read the parameters and the sampler only."
  )
  pvt <- hh_model$.__enclos_env__$private
  has_house <- !is.null(hh_id)
  out <- .engine$export_state(
    pvt$mcmc_ptr_, if (has_house) as.character(hh_id) else NA_character_
  )

  if (has_house) {
    out$dataset <- pvt$refactor_data_(out$dataset, "dataset")
    # long dataset
    long <- out$long_dataset
    if (all(dim(long) == 0)) out$long_dataset <- NULL
    else out$long_dataset <- pvt$refactor_data_(long, "long_dataset")
  }
  # parameters
  if (!length(.engine$param_keys(pvt$mcmc_ptr_)) &&
        all(dim(out$parameters) == 0)) out$parameters <- NULL
  if (has_house) {
    # group rates
    if (!length(.engine$group_rate_keys(pvt$mcmc_ptr_)) &&
          all(dim(out$group_rates) == 0)) out$group_rates <- NULL
    # contact rates
    if (!.engine$has_contact_rate(pvt$mcmc_ptr_) &&
          all(dim(out$contact_rates) == 0)) out$contact_rates <- NULL
    # infectivity profiles
    if (!.engine$has_infec_profile(pvt$mcmc_ptr_) &&
          !length(out$infectivity_profiles)) {
      out$infectivity_profiles <- NULL
    } else {
      lst <- out$infectivity_profiles
      df <- data.frame(id = lst$id, discr = lst$discr, func = lst$func,
                       min = lst$min, max = lst$max)
      df$func_args <- lst$func_args
      out$infectivity_profiles <- df
    }
    # model quantities
    if (!length(out$model_quantities)) out$model_quantities <- NULL
  }
  # sampler, absent while nothing is sampled
  smp <- out$.sampler
  out$.sampler <- if (length(smp)) {
    na_empty <- \(x) replace(x, !nzchar(x), NA_character_)
    moves <- as.data.frame(smp$moves)
    cols <- c("member", "proposal", "coord")
    moves[cols] <- lapply(moves[cols], na_empty)
    blocks <- as.data.frame(smp$blocks[names(smp$blocks) != "jump_cov"])
    blocks$member <- na_empty(blocks$member)
    blocks$jump_cov <- smp$blocks$jump_cov
    list(moves = moves, blocks = blocks)
  }
  out
}

#' @title Run Bayesian MCMC inference on a housepi model
#' @description
#' Run inference with Bayesian Markov Chain Monte Carlo,
#' by the component-wise Metropolis-Hastings algorithm
#' \insertCite{metropolis_equation_1953,hastings_monte_1970}{housepi}.
#' The latent variables of one member that declare no proposal
#' move together, by adaptive Metropolis
#' \insertCite{haario_adaptive_2001,andrieu_tutorial_2008}{housepi}.
#' Over the first half of `adapt` they move one at a time.
#' Their joint move then follows the covariance of their values
#' from the second quarter of `adapt` onwards.
#' A latent variable whose bounds or location read another of them,
#' directly or through a function, moves alone.
#' Parameters that declare no proposal move together the same way
#' when a reference couples them:
#' a baseline and each rate ratio of its rule groups
#' matching at least a tenth as many members or pairs as the reference.
#' A model can also request parameters that move together with a baseline,
#' see its `template_<model>` help page.
#' They take as many joint moves per step as they have members.
#' A parameter whose prior or bounds read another of them moves alone,
#' as does a rate matching no member or pair, or whose number is unknown.
#' @param hh_model `r rd_model_in_`
#' @param steps `r rd_smpl_steps_`
#' @param thin `r rd_smpl_thin_`
#' @param adapt `r rd_smpl_adapt_`
#' @param chains `r rd_smpl_chain_`
#' @param threads `integer(1)`
#' Number of threads, shared equally by the chains running at once,
#' each spreading its households over its share.
#' Runs after the same `set.seed()` call draw the same chains
#' at any `threads`.
#' @param id_draws `character()` | `list`
#' Infected members whose public latent draws are stored.
#' A `character()` lists individual identifiers.
#' A `list` takes entries `id`, individual identifiers, and `hh_id`,
#' households whose infected members are all stored.
#' If `NA`, one random eligible household is followed.
#' If `NULL`, nothing is stored.
#' @param progress `r rd_progress_`
#' @param warn `logical(1)` Whether to report the warnings the run raises.
#' A quiet run still records the convergence diagnostics it flagged,
#' read back from the `warnings` field of the returned object.
#' @return An object of class `HousepiPosteriorDraws`.
#' @details
#' With `threads` a multiple of `chains`, every chain runs at once.
#' Otherwise fewer do:
#' 4 chains on 6 threads run 2 at a time on 3 threads each,
#' and 4 chains on 3 threads run one at a time on all 3.
#'
#' The run warns on any of: `r paste(diagnostic_flags_, collapse = "; ")`.
#' R-hat and bulk ESS are the rank-normalized ones of
#' \insertCite{vehtari_rank_normalization_2021;textual}{housepi}.
#' @references
#' \insertAllCited{}
#' @examples
#' flu_1918_frederick |>
#'   new_housepi_model("final-size") |>
#'   add_parameter(import = "template") |>
#'   run_inference(steps = 50, thin = 2,
#'                 chains = 2) |>
#'   get_parameter_draws() |>
#'   posterior::summarize_draws()
#' @export
run_inference <- function(hh_model,
                          steps = 2000L, thin = 10L, adapt = 1000L,
                          chains = 4L, threads = 1L,
                          id_draws = NA,
                          progress = TRUE,
                          warn = TRUE) {
  stopifnothousepi(hh_model)
  # the whole run is quiet, its build included
  local_quiet_warnings_(!warn)
  ptr <- prepare_inference_(hh_model, steps, thin, adapt, chains,
                            threads, id_draws, draw_seed_())
  with_progress_(progress, "sample", 1, function(bar) {
    wrap_inference_(hh_model, .engine$run_inference_list(list(ptr), 1L, bar))
  })
}

#' Build the engine of `hh_model` for inference, seeded at `seed`.
#' Shared by [run_inference()] and the replicate sweep.
#' @return The engine pointer.
#' @noRd
prepare_inference_ <- function(hh_model, steps, thin, adapt, chains,
                               threads, id_draws, seed) {
  settings <- resolve_settings_(steps, thin, adapt, chains)
  id_draws <- resolve_id_draws_(id_draws)
  assert_posterior_memory_(hh_model, settings$total, id_draws)
  pvt <- hh_model$.__enclos_env__$private
  pvt$build_full_(
    steps = settings$steps, thin = thin, adapt = adapt,
    chains = chains, threads = threads, id_draws = id_draws, seed = seed
  )
  pvt$mcmc_ptr_
}

#' Stop when `n_sample` posterior draws of `hh_model` would not fit
#' in the room the memory in use leaves under R's vector memory limit.
#' * no limit -> no check, R raises no error.
#' @noRd
assert_posterior_memory_ <- function(hh_model, n_sample, id_draws) {
  # limit in bytes
  max_bytes <- mem.maxVSize() * 2^20
  if (!is.finite(max_bytes)) return(invisible())
  pvt <- hh_model$.__enclos_env__$private
  # @todo fewer latents and parameters when functions or private
  n_param_cols <- 2L * NROW(pvt$param_$defs)
  n_latent_defs <- NROW(pvt$augm_$defs)
  # structural upper bound: each named or default household at its largest
  # size, status filtering left to the engine
  hh_id <- get_dataset(hh_model)$hh_id
  max_hh_size <- if (!length(hh_id)) 1L else max(table(hh_id))
  n_members_upper <- length(unique(id_draws$id)) +
    length(unique(id_draws$hh_id)) * max_hh_size
  total_cols <- n_param_cols + n_members_upper * n_latent_defs +
    length(unique(hh_id)) + 1L
  # the column-by-sample product exceeds the 32-bit range on large runs
  needed <- 8 * as.double(n_sample) * as.double(total_cols)
  # minor collection overstates the memory in use:
  # * room under it -> holds;
  # * shortfall -> full collection.
  room <- max_bytes - gc(full = FALSE)["Vcells", "used"] * 8
  if (needed <= room) return(invisible())
  room <- max_bytes - gc()["Vcells", "used"] * 8
  if (needed <= room) return(invisible())
  # nolint next: object_usage_linter. read by the cli glue below
  size <- \(b) format(structure(b, class = "object_size"), units = "auto")
  stop_h(
    "Posterior draws would exceed R's memory limit.",
    "They need up to {size(needed)}, and {size(room)} are left.",
    paste("Raise the limit with {.fn mem.maxVSize}, increase {.arg thin},",
          "or reduce {.arg steps} and/or {.arg chains}.")
  )
}

#' Read the posterior draws of entry `i` of engine run `ran` for `hh_model`.
#' * run stopped by an error -> that error, rendered as any engine error.
#' @noRd
wrap_inference_ <- function(hh_model, ran, i = 1L) {
  if (!is.null(ran$error)) with_cpp_conditions(stop(ran$error))
  pvt <- hh_model$.__enclos_env__$private
  HousepiPosteriorDraws$new(pvt$compiled_$name, pvt$mcmc_ptr_,
                            ran$seconds[[i]], pvt$fingerprint_)
}

#' @title Simulate a household survey from a housepi model
#' @description
#' Simulate a household survey on the households of the loaded data.
#' With a reporting process, the generated data go through it,
#' see [add_reporting_process()].
#' @param hh_model `r rd_model_in_`
#' @param draws `r rd_draws_`
#' @param nrep `r rd_nrep_`
#' `r rd_nrep_prior_` If `NA`, one replicate.
#' Once set, even to 1, the replicates are returned as a list.
#' @param thin_ess `r rd_thin_ess_`
#' @param threads `r rd_thread_`
#' @param ... `r rd_dots_config_`
#' @return A `data.frame` for a single replicate,
#' when `nrep` is `NA` and `draws` is `NULL`, named values
#' or a one-row `data.frame`.
#' Otherwise a `list` with one `data.frame` per replicate.
#' With a `long_dataset` reporting process,
#' a `list(dataset = , long_dataset = )` of either.
#' @examples
#' flu_1918_frederick |>
#'   new_housepi_model("final-size") |>
#'   add_parameter(import = "template") |>
#'   run_simulation() |>
#'   head()
#' @export
run_simulation <- function(hh_model,
                           draws = NULL,
                           nrep = NA_integer_,
                           thin_ess = TRUE,
                           threads = 1L,
                           ...) {
  stopifnothousepi(hh_model)
  assert_advanced_settings_("time_step", ...)
  output <- hh_model$.__enclos_env__$private$simulate_(
    draws = draws, nrep = nrep, thin_ess = thin_ess,
    apply_long_dataset = identity, threads = threads, ...
  )
  output$parameters <- NULL
  # the arguments, not the replicates left after thinning, ask for one frame
  single <- is.null(draws) || is_named_values_(draws) || (
    is.data.frame(draws) && !posterior::is_draws(draws) && nrow(draws) == 1L
  )
  if (is.na(nrep) && single) output <- lapply(output, `[[`, 1L)
  if (length(output) == 1) output[[1]] else output
}
