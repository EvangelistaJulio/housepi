/*
 * export_bindings.cpp
 *
 * Copyright (C) 2024 Francesco Camaglia, Institut Pasteur
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
//' @field export_loglike_matrix Return matrix with loglikelihood values.
//' @field show_full_build Show all build messages.
//' @field show_summary Show a summary of the current database.
//' @field show_group_classes Show the group classes.

#include <Rcpp.h>
#include "model.h"

// Expose the Model and MCMC classes to R using Rcpp
RCPP_MODULE(model_housepi_rcpp)
{
  using namespace Rcpp;

  class_<MCMC>("McmcHousepi")
      .constructor()
      .method("import_setup", &MCMC::import_setup)
      .method("import_params", &MCMC::import_params)
      .method("import_data", &MCMC::import_data)
      .method("import_augm", &MCMC::import_augm)
      .method("read_setup", &MCMC::read_setup)
      .method("read_params", &MCMC::read_params)
      .method("read_data", &MCMC::read_data)
      .method("read_augm", &MCMC::read_augm)
      .method("infer", &MCMC::infer)
      .method("simulate", &MCMC::simulate)
      .method("evaluate", &MCMC::evaluate)
      .method("evaluate_marginal", &MCMC::evaluate_marginal)
      .method("export_setup", &MCMC::export_setup)
      .method("export_params", &MCMC::export_params)
      .method("export_data", &MCMC::export_data)
      .method("export_traces", &MCMC::export_traces)
      .method("export_loglike_matrix", &MCMC::export_loglike_matrix)
      .method("show_full_build", &MCMC::show_full_build)
      .method("show_summary", &MCMC::show_summary)
      .method("show_group_classes", &MCMC::show_group_classes);

  class_<Model>("ModelHousepi")
      .constructor()
      .derives<MCMC>("McmcHousepi");
}
