## housepi - Epidemiological Analysis of Household Outbreaks.
##
## Copyright © 2026 Institut Pasteur (Paris), CNRS and INSERM.
## Distributed under the terms of the GNU General Public License (GPLv3).
## See the LICENSE file for details.
##

rd_model_in_ <- "`HousepiModel` object."
rd_model_out_ <- paste(
  rd_model_in_,
  "Deep clone of input with new definitions."
)

rd_mdl_choice_ <- "`character(1)` Choice for the model.\\cr
One between:\\cr
- template name;\\cr
- path/to/folder with implementation file 'model.cpp';\\cr
- path/to/folder with previously compiled model;
- 'blank': template name for blank model.
- `NA`: fall back to 'blank' (default)."
rd_mdl_astmpl_ <- "`logical(1)` Store dynamic library as a local template.\\cr
If `FALSE` recompile the model at every new session (default)."
rd_mdl_force_ <- "`logical(1)` Compile even if cached model already exists.\\cr
Default is `FALSE`."
rd_mdl_quiet_ <- "`logical(1)` If `TRUE` hide compilation commands.\\cr
Default is `FALSE`."
rd_mdl_flags_ <- "Named arguments. Additional compilation options:\\cr
- `discr` : `logical(1)` If `TRUE` compile to discrete time.\\cr
- `flags` : `character()` Compilation flags.\\cr
- `libs` : `character()` Links to libraries.\\cr
- `paths` : `character()` Include paths."

############################################################################
############################################################################
############################################################################

rd_var_name_ <- "`character(1)` Name used by internal calls."
rd_var_value_ <- "`numeric(1)` Initial value.\\cr
Required only if constant or random initialization has been disabled."
rd_var_const_ <- "`logical(1)` Treat definition as constant.\\cr
If `TRUE`, additional Monte Carlo details are ignored.\\cr
Default is `FALSE`."

rd_var_discr_ <- "`logical(1)` Treat definition as discrete.\\cr
If applicable, the %s is discretized.\\cr
Default is `FALSE`."
rd_var_model_ <- "`character(1)` Choice for the %s distribution."
rd_var_args_ <- "named arguments for the %s distribution."

rd_var_min_ <- "`numeric(1)` Truncate to a minimum value."
rd_var_max_ <- "`numeric(1)` Truncate to a maximum value."

rd_var_proposal_ <- "`character(1)` Proposal used in the Monte Carlo.\\cr
One between:\\cr
- 'sample': sampled from the %s distribution;\\cr
- 'uniform': increments drawn from uniform distribution;\\cr
- 'normal': increments drawn from normal distribution;\\cr
- 'lognormal': log-increments drawn from lognormal distribution (for positive values only);
- `NA` (default): same as 'sample'.
"
rd_proposal_scale_ <- "`numeric(1)` Scale of the proposal (default is `1`).\\cr
Changed automatically during inference if `n_adapt`>0.\\cr
Ignored if proposal is 'sample'."

rd_par_discr_ <- sprintf(rd_var_discr_, "prior")
rd_par_model_ <- sprintf(rd_var_model_, "prior")
rd_par_args_ <- sprintf(rd_var_args_, "prior")
rd_par_proposal_ <- sprintf(rd_var_proposal_, "prior")

rd_aug_discr_ <- sprintf(rd_var_discr_, "model")
rd_aug_model_ <- sprintf(rd_var_model_, "model")
rd_aug_args_ <- sprintf(rd_var_args_, "model")
rd_aug_proposal_ <- sprintf(rd_var_proposal_, "model")

rd_inf_discr_ <- sprintf(rd_var_discr_, "profile")
rd_inf_args_ <- sprintf(rd_var_args_, "profile")

rd_list_args_ <- "Named arguments. Respective list of definition details."
rd_list_append_ <- "`logical(1)` Append to current definitions.\\cr
If `FALSE` drop all pre-existing definitions (default)."

rd_tmpl_append_ <- "`logical(1)` Append only missing definitions.\\cr
If `FALSE` drop all pre-existing definitions (default)."

############################################################################
############################################################################
############################################################################

rd_types_ <- "\\cr
One between:\\cr
- Scalar value;\\cr
- Parameter name (followed by `@P`);\\cr
- Dataset name (followed by `@D`)."

rd_aug_status_ <- "`integer(1)` Individual positive status.\\cr
Used to specify status-specific definitions.
If `NA`, associate definition to all infectious individuals (default)."
rd_aug_val_ <- paste("`numeric(1)` | `character(1)`. Initial value.",
                     rd_types_)
rd_aug_func_ <- "One between:\\cr
- `character(1)` Name of augmented function;\\cr
- rlang expression."
rd_aug_funargs_ <- "named `list` Arguments of the augmented function."
rd_aug_min_ <- paste(rd_var_min_, rd_types_)
rd_aug_max_ <- paste(rd_var_max_, rd_types_)

############################################################################
############################################################################
############################################################################

rd_grr_desc_ <- "Add definitions of a group rate by given data feature(s).\\cr
Group rates are relative to a reference class that can be specified.\\cr
If no MCMC options are provided, the prior is assumed to be\\cr
\\eqn{\\mathrm{lognormal}[\\mu=0, \\sigma=3]} with support \\eqn{[0, 10]}."
rd_grr_by_ <- "`list` Features associated to group rate as data headers.\\cr
Corresponding data columns have to be categorical."
rd_grr_prefix_ <- "`character(1)` Prefix to names of associated parameters.\\cr
If `NA` uses the first 3 letters of rate name (default)."
rd_grr_refer_ <- "named `list` Features associated to reference group.\\cr
Default is the group with the largest sample size."
rd_grr_indep_ <- "`logical(1)` Consider features as independent.\\cr
If `FALSE` features are combined (default)."

############################################################################
############################################################################
############################################################################

rd_progress_ <- "`integer(1)` | `logical(1)` Control progressbar visualization.\\cr
One between:\\cr
- If `integer`, prints progress every `progress` seconds;\\cr
- If `TRUE` prints progress every second (default);\\cr
- If `FALSE` or `0`, suppresses progress display."
rd_house_id_ <- "`character(1)` ID of the household of interest."

rd_add_data_wide_ <- "`data.frame` Household data used for analysis (wide format)."

rd_simul_nsim_ <-  "`integer(1)` Number of simulated outbreaks.\\cr
Ignored if `param_table` is provided, else simulations share parameters."
rd_simul_observ_ <-  "`function` Observation process.\\cr
Function expects as input the `data.frame` generated by simulation.\\cr
If `NULL`, observation process is ignored."
rd_simul_obs_long_ <-  "`function` Optional function.\\cr
Expexts processed `data_wide` and returns a `data_long` data.frame.\\cr
When provided, each simulation returns `list(wide = ..., long = ...)`."

rd_smpl_data_wide_ <- paste(rd_add_data_wide_,
                            "If `NULL`, uses the previously stored data.")

rd_add_data_long_ <- "`data.frame` Additional individual long format data."

rd_smpl_steps_ <- "`integer(1)` Total number of MCMC steps per chain."
rd_smpl_burn_ <- "`integer(1)` Number of initial steps to discard."
rd_smpl_thin_ <- "`integer(1)` Number of steps to skip between samples."
rd_smpl_adapt_ <- "`integer(1)` Number of adaptation steps before sampling."
rd_smpl_chain_ <- "`integer(1)` Number of independent MCMC chains to run."

rd_thread_ <- "`integer(1)` Number of threads for parallel household runs."

rd_smpl_rndin_ <- "`logical(1)` Control random initialization of variables:\\cr
- If `TRUE` initial values are drawn from respective models (default);\\cr
- If `FALSE` uses initial value definitions."
rd_smpl_iddrw_ <- "`character()` Individuals identifiers.\\cr
Store draws of public augmented variables for the chosen individuals."
rd_smpl_cri_level_ <- "`numeric(1)` Level for credible interval\\cr
It must be a value in (0.5, 1.0). Default is 0.95."

############################################################################
############################################################################
############################################################################

rd_evi_scheme_ <- "`character(1)` Scheme used for evidence estimation.\\cr
One between\\cr
- 'power-fractions': uses inverse temperatures \\eqn{\\frac{r^5}{{n_\\mathrm{rungs}}^5}}{r^5 / n_rungs^5} (default).
- 'adaptive': uses adaptive thermodynamic ladder."
rd_evi_rungs_ <- "`integer(1)` Number of rungs in thermodynamic ladder."
rd_param_table_ <- "One between\\cr
- `data.frame`;
- `posterior::draws`;
- [HousepiPosteriorDraws][housepi::HousepiPosteriorDraws].
"
rd_measure_dots_ <- "Named arguments. `function`s to compute respective measures from data.\\cr
Each function should return one between `numeric()`, `list` or `data.frame`."

rd_measure_long_ <- "`function` | `list` of `function`s. Applied to `data_long` (derived via `obs_process_long`)."

rd_ppc_observ_ <- paste(rd_simul_observ_,
                        "This applies only on simulations, not on data.")

############################################################################
############################################################################
############################################################################

rd_hhsampler_size_ <- "`integer(1)` Number of households to sample."
rd_hhsampler_asc_ <- "`logical(1)` Set ascending hierarchy between feature numbers.\\cr
If `FALSE` is descending (default)."
