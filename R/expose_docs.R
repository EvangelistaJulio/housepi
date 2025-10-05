#      __
#     / /_  ____  __  __ ____ __  ____  __
#    / __ \/ __ \/ / / / ___/ _ \/ __ \/ /
#   / / / / /_/ / /_/ (__  )  __/ /_/ / /
#  /_/ /_/\____/\____/____/\___/ ____/_/
#                             /_/
#
# Copyright (C) 2025 CAMAGLIA Francesco - Institut Pasteur, Paris (France).
#
# expose_docs.R
# --------------
# Inline calls for shared Roxygen documentation
# between R6 methods and wrappers

rd_model_in_ <- "`HousepiModel` object."
rd_model_out_ <- "`HousepiModel`\\cr
A deep clone of the input object with updated definitions."

rd_mdl_choice_ <- "`character(1)` with the choice for the model.\\cr
One between:\\cr
- template name;\\cr
- path/to/folder with implementation files \"model.cpp\" and \"model.h\";\\cr
- path/to/folder with previously compiled model;
- \"blank\": template name for blank model.
- `NA`: fall back to \"blank\" (default)."
rd_mdl_astmpl_ <- "`logical(1)` to store the dynamic library as a local template.\\cr
If `FALSE` recompile the model at every new session (default)."
rd_mdl_force_ <- "`logical(1)` to compile even if a cached model already exists.\\cr
Default is `FALSE`."
rd_mdl_quiet_ <- "`logical(1)` to hide compilation commands and warnings.\\cr
Default is `TRUE`."

############################################################################
############################################################################
############################################################################

rd_var_name_ <- "`character(1)` name used by internal calls."
rd_var_value_ <- "`numeric(1)` with initial value.\\cr
Required only if constant or random initialization has been disabled."
rd_var_const_ <- "`logical(1)` to treat definition as constant (other details are ignored).\\cr
Default is `FALSE`."

rd_var_discr_ <- "`logical(1)` to treat definition as discrete and the %s is discretized (if applicable).\\cr
Default is `FALSE`."
rd_var_model_ <- "`character(1)` wit the choice for the %s distribution."
rd_var_args_ <- "named `numeric()` or `list` with custom arguments for the %s distribution."

rd_var_min_ <- "`numeric(1)` to truncate to a minimum value."
rd_var_max_ <- "`numeric(1)` to truncate to a maximum value."

rd_var_rule_ <- "`character(1)` with the proposer rule used during inference.\\cr
One between:\\cr
- `NA` (default): see \"sample\";\\cr
- \"sample\": sampled from the %s distribution;\\cr
- \"unifrm\": increment is drawn from a uniform distribution;\\cr
- \"normal\": increment is drawn from a normal distribution;\\cr
- \"lognormal\": increment id drawn from a lognormal distribution (for positive values only).
"
rd_var_scale_ <- "`numeric(1)`with the scale of the proposer rule (if applicable).\\cr 
Changed automatically during inference if `n_adpat`>0.\\cr
Default is `1`."

rd_par_discr_ <- sprintf(rd_var_discr_, "prior")
rd_par_model_ <- sprintf(rd_var_model_, "prior")
rd_par_args_ <- sprintf(rd_var_args_, "prior")
rd_par_rule_ <- sprintf(rd_var_rule_, "prior")

rd_aug_discr_ <- sprintf(rd_var_discr_, "model")
rd_aug_model_ <- sprintf(rd_var_model_, "model")
rd_aug_args_ <- sprintf(rd_var_args_, "model")
rd_aug_rule_ <- sprintf(rd_var_rule_, "model")

rd_inf_discr_ <- sprintf(rd_var_discr_, "profile")
rd_inf_args_ <- sprintf(rd_var_args_, "profile")

rd_list_args_ <- "Named arguments with respective list of definition details."
rd_list_append_ <- "`logical(1)` to append to current definitions.\\cr
If `FALSE` drop all pre-existing definitions (default)."
rd_list_from_ <- "named `list` of definitions."

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

rd_aug_status_ <- "`integer(1)` with positive status to specify individual definitions.\\cr
If NA the generic defition is associated to all infectious individuals (default)."
rd_aug_val_ <- paste("`numeric(1)` or `character(1)` with the initial value.", rd_types_)
rd_aug_rndin_ <- "`character(1)` to draw initial value from model distribution.\\cr
Default is `FALSE`."
rd_aug_func_ <- "`character(1)` with the name of the augmented function."
rd_aug_funargs_ <- "named `list` with the arguments of the augmented function."
rd_aug_min_ <- paste(rd_var_min_, rd_types_)
rd_aug_max_ <- paste(rd_var_max_, rd_types_)
rd_aug_expr_ <- "Provide an expression."

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

rd_add_data_wide_ <- "`data.frame` with the household data used for analysis in wide format."

rd_simul_nsim_ <-  "`integer(1)` with the number of simulated outbreaks."
rd_simul_observ_ <-  "A `function` with the observation process.\\cr
It expects as input the `data.frame` generated by simulation.\\cr
If `NULL`, observation process is ignored."
rd_smpl_data_wide_ <- paste(rd_add_data_wide_,
                            "If `NULL`, uses the previously stored data.")

rd_add_data_long_ <- "`data.frame` with the individual data (long format) used for analysis."

rd_smpl_size_ <- "`integer(1)` with the number of posterior draws stored for each chain."
rd_smpl_steps_ <- "`integer(1)`  with the total number of MCMC steps for each chain."
rd_smpl_burn_ <- "`integer(1)` with the number of initial steps to discard as burn-in."
rd_smpl_thin_ <- "`integer(1)` with the thinning interval, i. e. the number of steps to skip between each saved sample."
rd_smpl_adapt_ <- "`integer(1)` with the number of adaptation steps before starting the main sampling."
rd_smpl_chain_ <- "`integer(1)` with the number of independent MCMC chains to run."
rd_smpl_thread_ <- "`integer(1)` with the number of threads to run in parallel."

rd_smpl_rndin_ <- "`logical(1)` to control random initialization of variables:\\cr
- If `TRUE` initial values are randomly drawn from respective distributions (default);\\cr
- If `FALSE` uses initial value definitions."
rd_smpl_iddrw_ <- "`character()` with individuals identifiers.\\cr
Store draws of public augmented variables for the chosen individuals."
rd_smpl_cri_level_ <- "`numeric(1)` with the level for the credible interval in (0.5, 1.0).\\cr
Default is 0.95."

############################################################################
############################################################################
############################################################################

rd_evi_scheme_ <- "`character(1)` with the scheme used for evidence estimation.\\cr
One between\\cr
- \"power-fractions\": uses inverse temperatures \\eqn{\\frac{r^5}{{n_\\mathrm{rungs}}^5}}{r^5 / n_rungs^5}  (default).
- \"adaptive\": uses adaptive thermodynamic ladder."
rd_evi_rungs_ <- "`integer(1)` \\cr
Number of rungs in the thermodynamic ladder."

rd_param_table_ <- "One between\\cr
- `data.frame`;
- `posterior::draws`;
- `HousepiPosteriorDraws`.
"
rd_measure_ <- "A `function` or `list` of `function` used to compute the desired measure from data."
rd_measure_dots_ <- "Named arguments with respective measure `function`.\\cr
If `measure` is provided, arguments are ignored."

############################################################################
############################################################################
############################################################################

rd_hhsampler_size_ <- "`integer(1) with the number of households to sample."
rd_hhsampler_from_ <- "A `list` of `list`s with encoding definitions for each feature."
rd_hhsampler_asc_ <- "A `logical(1)` to set ascending hierarchy between feature numbers.\\cr
If `FALSE` is descending (default)."
