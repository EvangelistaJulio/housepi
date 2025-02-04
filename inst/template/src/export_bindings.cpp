/*
 * export_bindings.cpp
 *
 * Copyright (C) 2025 Francesco Camaglia, Institut Pasteur
 */

//' @name ModelHousepi
//' @title Internal housepi model
//' @description Perform inference and simulation from the model.
//' @field infer Perform MCMC inference of the model parameters. \itemize{
//'   \item Parameter: output_file - The file path where the output will be saved
//'   \item Parameter: progbar - show the progress bar
//' }
//' @field simulate Load a simulated database.
//' @field evaluate Evaluate the loglikelihood according to the current state.
//' @field evaluate_marginal Evaluate the logarithm of the marginal likelihood with the method of power posteriors. \itemize{
//'   \item Parameter: n_rungs
//'   \item Parameter: progbar - show the progress bar
//'   \item Parameter: power_fractions - if true uses power fractions scheme instead of adaptive scheme
//' }
//' @field export_setup Return a List with current setup.
//' @field export_params Return a DataFrame with current parameters.
//' @field export_data Return a DataFrame with current database.
//' @field export_traces Return a DataFrame with MCMC traces.
//' @field export_household_loglike_traces Return a NumericMatrix with household loglikelihood traces.
//' @field export_infec_profile Return a DataFrame values of the infectivity profile.
//' @filed import_setup Import setup from R List.
//' @filed import_infec_profile Import infectivity profile definition from R List.
//' @filed import_params Import parameter definitions from R Data.Frame.
//' @filed import_data Import data definitions from R Data.Frame.
//' @filed import_augm Import augmented variable definitions from R Data.Frame.
//' @field export_household_loglike_traces Return matrix with loglikelihood values.
//' @field show_full_build Show all build messages.
//' @field show_summary Show a summary of the current database.

#include <Rcpp.h>
#include "model.h"

// Expose the Model and MCMC classes to R using Rcpp
RCPP_MODULE(model_housepi_rcpp)
{
  using namespace Rcpp;

  class_<MCMC>("McmcHousepi")
      .constructor()
      .method("infer", &MCMC::infer)
      .method("simulate", &MCMC::simulate)
      .method("evaluate", &MCMC::evaluate)
      .method("evaluate_marginal", &MCMC::evaluate_marginal)
      .method("import_setup", &MCMC::import_setup)
      .method("import_infec_profile", &MCMC::import_infec_profile)
      .method("import_params", &MCMC::import_params)
      .method("import_data", &MCMC::import_data)
      .method("import_augm", &MCMC::import_augm)
      .method("export_setup", &MCMC::export_setup)
      .method("export_params", &MCMC::export_params)
      .method("export_data", &MCMC::export_data)
      .method("export_traces", &MCMC::export_traces)
      .method("export_household_loglike_traces", &MCMC::export_household_loglike_traces)
      .method("export_infec_profile", &MCMC::export_infec_profile)
      .method("show_full_build", &MCMC::show_full_build)
      .method("show_summary", &MCMC::show_summary);

  class_<Model>("ModelHousepi")
      .constructor()
      .derives<MCMC>("McmcHousepi");
}
