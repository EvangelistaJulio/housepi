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
# expose_docs.R
# --------------
# Inline calls for shared Roxygen documentation
# between R6 methods and wrappers

rd_model_in_ <- "`HousepiModel` object."
rd_model_out_ <- "A deep clone of the input `HousepiModel` with new definitions."

rd_mdl_choice_ <- "`character(1)` with the choice for the model.\\cr
One between:\\cr
- template name;\\cr
- path/to/folder with implementation files 'model.cpp' and 'model.h';\\cr
- path/to/folder with previously compiled model;
- 'blank': template name for blank model.
- `NA`: fall back to 'blank' (default)."
rd_mdl_astmpl_ <- "`logical(1)` to store the dynamic library as a local template.\\cr
If `FALSE` recompile the model at every new session (default)."
rd_mdl_force_ <- "`logical(1)` to compile even if a cached model already exists.\\cr
Default is `FALSE`."
rd_mdl_quiet_ <- "`logical(1)` to hide compilation commands.\\cr
Default is `FALSE`."
rd_mdl_flags_ <- "Named arguments with additional compilation options:\\cr
- If `discr = TRUE`, uses the model in discrete time, if applicable."

############################################################################
############################################################################
############################################################################

rd_var_name_ <- "`character(1)` name used by internal calls."
rd_var_value_ <- "`numeric(1)` with initial value.\\cr
Required only if constant or random initialization has been disabled."
rd_var_const_ <- "`logical(1)` to treat definition as constant.\\cr
If `TRUE`, additional Monte Carlo details are ignored.\\cr
Default is `FALSE`."

rd_var_discr_ <- "`logical(1)` to treat definition as discrete.\\cr
If applicable, the %s is discretized.\\cr
Default is `FALSE`."
rd_var_model_ <- "`character(1)` with the choice for the %s distribution."
rd_var_args_ <- "named arguments for the %s distribution."

rd_var_min_ <- "`numeric(1)` to truncate to a minimum value."
rd_var_max_ <- "`numeric(1)` to truncate to a maximum value."

rd_var_proposal_ <- "`character(1)` with proposal used in the Monte Carlo.\\cr
One between:\\cr
- 'sample': sampled from the %s distribution;\\cr
- 'uniform': increments drawn from uniform distribution;\\cr
- 'normal': increments drawn from normal distribution;\\cr
- 'lognormal': log-increments drawn from lognormal distribution (for positive values only);
- `NA` (default): same as 'sample'.
"
rd_var_proposal_scale_ <- "`numeric(1)` with scale of the proposal.\\cr 
Changed automatically during inference if `n_adpat`>0.\\cr
Not applicable if proposal is 'sample'. Otherwise scale is `1` by default."

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

rd_list_args_ <- "Named arguments with respective list of definition details."
rd_list_append_ <- "`logical(1)` to append to current definitions.\\cr
If `FALSE` drop all pre-existing definitions (default)."

rd_tmpl_append_ <- "`logical(1)` to append only missing definitions.\\cr
If `FALSE` drop all pre-existing definitions (default)."

############################################################################
############################################################################
############################################################################

rd_types_ <- "\\cr
One between:\\cr
- Scalar value;\\cr
- Parameter name (followed by `@P`);\\cr
- Dataset name (followed by `@D`)."

rd_aug_status_ <- "`integer(1)` with household member positive status.\\cr
It is used to specify status-specific definitions.
If `NA`, associate the definition to all other infectious members (default)."
rd_aug_val_ <- paste("`numeric(1)` or `character(1)` with initial value.",
                     rd_types_)
rd_aug_func_ <- "One between:\\cr
- `character(1)` with the name of the augmented function;\\cr
- rlang expression."
rd_aug_funargs_ <- "named `list` with the arguments of the augmented function."
rd_aug_min_ <- paste(rd_var_min_, rd_types_)
rd_aug_max_ <- paste(rd_var_max_, rd_types_)

############################################################################
############################################################################
############################################################################

rd_grr_by_ <- "`list` with feature(s) associated to group rate as data headers.\\cr
Corresponding data columns have to be categorical."
rd_grr_prefix_ <- "`character(1)` with prefix to parameter names associated to the group rate.\\cr
If `NA` uses the first 3 letters of rate name (default)."
rd_grr_refer_ <- "named `list` to choose reference group providing corresponding features.\\cr
Default is the group with the largest sample size."
rd_grr_indep_ <- "`logical(1)` to consider features as independent.\\cr
If `FALSE` features are combined (default)."

############################################################################
############################################################################
############################################################################

rd_seed_ <- "`integer(1)` with the seed for the random number generator.\\cr
If `NA` uses current R seed (default)."
rd_progress_ <- "`integer(1)` or `logical(1)` to control progress visualization.\\cr
One between:\\cr
- If `integer`, prints progress every `progress` seconds;\\cr
- If `TRUE` prints progress every second (default);\\cr
- If `FALSE` or `0`, suppresses progress display."
rd_house_id_ <- "A `character(1)` with the ID of the household of interest."

rd_add_data_wide_ <- "`data.frame` with the household data used for analysis in wide format."

rd_simul_nsim_ <-  "`integer(1)` with the number of simulated outbreaks."
rd_simul_observ_ <-  "A `function` with the observation process.\\cr
It expects as input the `data.frame` generated by simulation.\\cr
If `NULL`, observation process is ignored."
rd_smpl_data_wide_ <- paste(rd_add_data_wide_,
                            "If `NULL`, uses the previously stored data.")

rd_add_data_long_ <- "`data.frame` with additional individual long format data."

rd_smpl_steps_ <- "`integer(1)`  with total number of MCMC steps per chain."
rd_smpl_burn_ <- "`integer(1)` with number of initial steps to discard."
rd_smpl_thin_ <- "`integer(1)` with number of steps to skip between each saved sample."
rd_smpl_adapt_ <- "`integer(1)` with number of adaptation steps before sampling."
rd_smpl_chain_ <- "`integer(1)` with number of independent MCMC chains to run."
rd_smpl_thread_ <- "`integer(1)` with number of threads to run in parallel."

rd_smpl_rndin_ <- "`logical(1)` to control random initialization of variables:\\cr
- If `TRUE` initial values are drawn from respective distributions (default);\\cr
- If `FALSE` uses initial value definitions."
rd_smpl_iddrw_ <- "`character()` with individuals identifiers.\\cr
Store draws of public augmented variables for the chosen individuals."
rd_smpl_cri_level_ <- "`numeric(1)` with the level for credible interval\\cr
It must be a value in (0.5, 1.0). Default is 0.95."

############################################################################
############################################################################
############################################################################

rd_evi_scheme_ <- "`character(1)` with scheme used for evidence estimation.\\cr
One between\\cr
- 'power-fractions': uses inverse temperatures \\eqn{\\frac{r^5}{{n_\\mathrm{rungs}}^5}}{r^5 / n_rungs^5} (default).
- 'adaptive': uses adaptive thermodynamic ladder."
rd_evi_rungs_ <- "`integer(1)` wiht number of rungs in thermodynamic ladder."

rd_param_table_ <- "One between\\cr
- `data.frame`;
- `posterior::draws`;
- [HousepiPosteriorDraws][housepi::HousepiPosteriorDraws].
"
rd_measure_dots_ <- "Named arguments with `function`s to compute the respective measure from data.\\cr
Each function should return one between `numeric()`, `list` or `data.frame`."

############################################################################
############################################################################
############################################################################

rd_hhsampler_size_ <- "`integer(1) with the number of households to sample."
rd_hhsampler_asc_ <- "`logical(1)` to set ascending hierarchy between feature numbers.\\cr
If `FALSE` is descending (default)."
