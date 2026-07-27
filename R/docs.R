## housepi - Epidemiological Analysis of Household Outbreaks.
##
## Copyright © 2024-2026 Institut Pasteur (Paris), CNRS and INSERM.
## Distributed under the terms of the GNU General Public License (GPLv3).
## See the LICENSE file for details.
##

rd_model_in_ <- "`HousepiModel` object."
rd_model_out_ <- paste(
  rd_model_in_,
  "Deep clone of input with new definitions."
)

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

############################################################################
############################################################################
############################################################################

rd_var_name_ <- "`character(1)` Name used by internal calls."
rd_var_value_ <- "`numeric(1)` Initial value.
Required only if constant or random initialization has been disabled."
rd_var_const_ <- "`logical(1)` Treat definition as constant.
If `TRUE`, additional Monte Carlo details are ignored."

rd_var_discr_ <- "`logical(1)` Treat definition as discrete.
If applicable, the %s is discretized."
rd_var_prior_ <- "`character(1)` Choice for the prior distribution."
rd_var_model_ <- "`character(1)` Choice for the augmented model distribution."
rd_var_args_ <- "named arguments for the %s distribution."

rd_var_min_ <- "`numeric(1)` Truncate to a minimum value."
rd_var_max_ <- "`numeric(1)` Truncate to a maximum value."

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
# Advanced Monte Carlo knobs supplied through `...`, all prefixed `proposal_`.
rd_proposal_opts_ <- "
Advanced Monte Carlo options may also be passed through `...`, each prefixed
`proposal_`. For a random walk, `proposal_scale` sets the step scale: if
omitted it is adapted during inference when `adapt > 0`; if a positive numeric
it is fixed and never adapted. For independence sampling from a named
distribution, the distribution parameters are supplied the same way (e.g.
`proposal_mean`, `proposal_sd`)."

rd_par_discr_ <- sprintf(rd_var_discr_, "prior")
rd_par_args_ <- paste0(sprintf(rd_var_args_, "prior"), rd_proposal_opts_)

rd_aug_discr_ <- sprintf(rd_var_discr_, "prior")
rd_aug_args_ <- paste0(sprintf(rd_var_args_, "model"), rd_proposal_opts_)

rd_obs_dist_ <- "`character(1)` Choice for the observation model distribution.
Reserved to member-level observations."
rd_obs_level_ <- "`character(1)` Level the datum is observed at.
One between `\"member\"`, one value per positive status individual, and
`\"house\"`, one value per household reduced over all its members.
In mode 'multiple', applies to every entry not setting its own."
rd_obs_args_ <- sprintf(rd_var_args_, "observation")

rd_inf_discr_ <- sprintf(rd_var_discr_, "profile")
rd_inf_args_ <- sprintf(rd_var_args_, "profile")

rd_import_ <- "Named `list` of pre-built definitions to import.
Cannot be combined with single and/or multiple definitions."
rd_list_append_ <- "`logical(1)` Append to current definitions.
If `FALSE` drop all pre-existing definitions."

rd_tmpl_append_ <- "`logical(1)` Append only missing definitions.
If `FALSE` drop all pre-existing definitions."

############################################################################
############################################################################
############################################################################

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
rd_aug_func_ <- "One between:
- `character(1)` Name of augmented function;
- rlang expression."
rd_aug_funargs_ <- "named `list` Arguments of the augmented function."
rd_aug_min_ <- paste(rd_var_min_, rd_types_)
rd_aug_max_ <- paste(rd_var_max_, rd_types_)

############################################################################
############################################################################
############################################################################

rd_risk_prefix_ <- "`character(1)` Prefix to names of associated parameters."
rd_risk_desc_ <- "Add definitions of a %s risk by given data covariate(s).
Risks are relative to a reference class that can be specified:
each parameter multiplies the baseline risk of its reference class.
If no MCMC options are provided, the prior is assumed to be
\\eqn{\\mathrm{lognormal}[\\mu=0, \\sigma=3]} with support \\eqn{[0, 10]}.
Bypass automatic generation by providing a custom pre-built logic table."

rd_grr_desc_ <- sprintf(rd_risk_desc_, "group")
rd_grr_by_ <- "`list` Covariates associated to group risk as data headers.
Corresponding data columns have to be categorical."
rd_grr_prefix_ <- paste(
  rd_risk_prefix_,
  "If `NA` uses the first 3 letters of risk name (default)."
)
rd_grr_refer_ <- "named `list` Covariates associated to reference group.
Pass `\".none\"` to disable reference, giving every risk its own parameter.
Default reference group is the one with the largest sample size."
rd_grr_joint_ <- "`logical(1)` Treat covariates as joint.
If `FALSE`, covariates are independent (default)."
rd_grr_strata_ <- "`data.frame` Individual covariates and relative `param_name`.
Each row defines one stratum. If `NA`, then match any value.
Use `.reference` in `param_name` to designate the reference group.
Mutually exclusive with `by`."
rd_grr_params_ <- "named `list` Per-parameter MCMC setups as in `param_name`.
Overrides shared prior/args for parameters specified in the list.
Ignored if automatic definition."

rd_ctr_desc_ <- sprintf(rd_risk_desc_, "contact")
rd_ctr_by_ <- "`character()` Covariates for from/to matching as data headers.
Corresponding data columns have to be categorical."
rd_ctr_prefix_ <- rd_risk_prefix_
rd_ctr_refer_ <- "named `list` Covariates associated to reference pair.
Names can use `from_`/`to_` prefixes matching `by`,
or un-prefixed names (applied to both from and to)."
rd_ctr_joint_ <- rd_grr_joint_
rd_ctr_symmetric_ <- "`logical(1)` Treat contact as symmetric.
If `TRUE`, contact pairs infector-infectee are treated as equivalent (default)."
rd_ctr_strata_ <- rd_grr_strata_
rd_ctr_params_ <- rd_grr_params_
rd_ctr_import_ <- "Named `list` of pre-built definition to import.
Cannot be combined with the single definition."

############################################################################
############################################################################
############################################################################

rd_progress_ <- "`integer(1)` | `logical(1)` Control progressbar visualization.
One between:
- If `integer`, prints progress every `progress` seconds;
- If `TRUE` prints progress every second (default);
- If `FALSE` or `0`, suppresses progress display."
rd_dots_config_ <- "Advanced configuration."
rd_house_id_ <- "`character(1)` ID of the household of interest."

rd_summ_cases_ <- "`data.frame`, one row per metric:
\\describe{
  \\item{level}{`individual` or `household`}
  \\item{metric}{quantity counted in the row}
  \\item{count}{absolute count}
  \\item{percent}{percentage of the relevant total}
}"
rd_summ_sizes_ <- "`data.frame`, one row per metric:
\\describe{
  \\item{metric}{quantity reported in the row}
  \\item{value}{household size for that metric}
  \\item{count}{absolute count}
  \\item{percent}{percentage of individuals}
}"

rd_add_wide_ <- "`data.frame` Household data used for analysis (wide format).
Optional at construction; if omitted, attach it later with
[add_data_wide()][housepi::add_data_wide]."

rd_simul_nrep_ <- "`integer(1)` Number of replicates to generate.
If `param_table` is `NULL`, all replicates share the same parameters.
Otherwise, if `nrep` is smaller than number of `param_table` rows,
the table is sub-sampled to `nrep` evenly-spaced rows.
If `NA` produces a single repetition (default)."
rd_max_nrep_ <- "`integer(1)` Larger number of replicates.
The parameter table (optionally ESS-thinned) is thinned to `max_nrep` rows."
rd_min_hhsize_ <- "`integer(1)` Minimum household size retained in the outputs.
If `NA`, no lower bound is applied."
rd_max_hhsize_ <- "`integer(1)` Maximum household size retained in the outputs.
If `NA`, no upper bound is applied."
rd_simul_observ_ <- "`function` Observation process.
Function expects as input the `data.frame` generated by simulation.
If `NULL`, observation process is ignored."
rd_simul_obs_long_ <- "`function` Optional function.
Expects processed `data_wide` and returns a `data_long` data.frame.
When provided, each simulation returns `list(wide = ..., long = ...)`."

rd_add_data_long_ <- "`data.frame` Additional individual long format data."

rd_smpl_steps_ <- "`integer(1)` Number of MCMC iterations per chain."
rd_smpl_burn_ <- "`integer(1)` Number of initial iterations within `steps` to
discard before storing draws."
rd_smpl_thin_ <- "`integer(1)` Keep one iteration every `thin` after burn-in."
rd_smpl_adapt_ <- "`integer(1)` Number of iterations to adapt proposal scales.
They are run before `steps` and are not stored as draws.
If no proposals are adaptable and `burn = 0`, repurposed as burn-in instead."
rd_smpl_chain_ <- "`integer(1)` Number of independent MCMC chains to run."

rd_thread_ <- "`integer(1)` Number of threads for parallel household runs."

rd_smpl_iddrw_ <- "Infected members whose public augmented draws are stored.
Either a `character()` of individual identifiers,
or a `list` with entries `id` (individual identifiers)
and/or `house_id` (whole household infected members are stored).
When unset, one random eligible household is followed (default).
Pass `\".none\"` to store nothing."
rd_smpl_cri_level_ <- "`numeric(1)` Level for credible interval.
It must be a value in (0.5, 1.0)."
rd_point_estimate_ <- "`character(1)` Summary statistic used as the predicted
point estimate across simulations.
One between `\"mean\"` and `\"median\"`."

############################################################################
############################################################################
############################################################################

rd_natural_log_ <- "Logarithms are natural (base \\eqn{e})."

rd_pwr_frac_ <- "\\eqn{\\frac{r^5}{{n_\\mathrm{rungs}}^5}}{r^5 / n_rungs^5}"
rd_evi_scheme_ <- sprintf("`character(1)` Scheme used for evidence estimation.
One between
\\describe{
  \\item{`\"adaptive\"`}{uses adaptive thermodynamic ladder (default).}
  \\item{`\"power-fractions\"`}{uses inverse temperatures %s.}
}", rd_pwr_frac_)
rd_evi_rungs_ <- "`integer(1)` Number of rungs in thermodynamic ladder."
rd_evi_method_ <- "`character(1)` Estimator computed from the sampled ladder.
One between
\\describe{
  \\item{`\"stepping-stone\"`}{stepping-stone sampling (default).}
  \\item{`\"thermodynamic\"`}{thermodynamic integration.}
}"
rd_param_table_ <- "One between
- `data.frame`;
- `posterior::draws`;
- [HousepiPosteriorDraws][housepi::HousepiPosteriorDraws].
"
rd_thin_ess_ <- "`logical(1)`. Whether to thin `param_table` to its
bulk effective sample size (ESS).
Ignored if `param_table` has no posterior draws."
rd_measure_dots_ <- "Named arguments. `function`s to compute measures from data.
Each function should return one between `numeric()`, `list` or `data.frame`."

rd_measure_long_ <- paste(
  "`function` | `list` of `function`s.",
  "Applied to `data_long` (derived via `obs_process_long`)."
)
