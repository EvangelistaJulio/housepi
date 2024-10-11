/*
 * export_bindings.cpp
 *
 * Copyright (C) 2024 Francesco Camaglia, Institut Pasteur
 */

//' @name Modelhousepi12345
//' @title Custom model
//' @description Perform inference and simulation from the model housepi12345.
//' @field new Initialize a household model object.
//' @field build Construct the model configurations and associate the data. \itemize{
//'   \item Parameter: setup_file - setup configuration file
//'   \item Parameter: params_file - parameter definition file
//'   \item Parameter: data_file - database file
//'   \item Parameter: augm_file - augmented data configurations
//'   \item Parameter: extra_file - file path for any extra configurations
//'   \item Parameter: verbose - show build information
//' }
//' @field infer Perform MCMC inference of the model parameters. \itemize{
//'   \item Parameter: output_file - The file path where the output will be saved
//'   \item Parameter: progbar - show the progress bar
//' }
//' @field simulate Load a simulated database.
//' @field evaluate Evaluate the loglikelihood according to the current state.
//' @field marginal Evaluate the logarithm of the marginal likelihood with the method of power posteriors. \itemize{
//'   \item Parameter: num_rungs
//'   \item Parameter: progbar - show the progress bar
//'   \item Parameter: power_fractions - if true uses power fractions scheme instead of adaptive scheme
//' }
//' @field setup Return a List with current setup.
//' @field params Return a DataFrame with current parameters.
//' @field traces Return a DataFrame with MCMC traces.
//' @field loglike_matrix Return matrix with loglikelihood values.
//' @field data Return a DataFrame with current database.
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
      .method("setup", &MCMC::export_setup)
      .method("params", &MCMC::export_params)
      .method("data", &MCMC::export_data)
      .method("traces", &MCMC::export_traces)
      .method("loglike_matrix", &MCMC::export_loglike_matrix)
      .method("show_summary", &MCMC::show_summary)
      .method("show_group_classes", &MCMC::show_group_classes);

  class_<Model>("ModelHousepi")
      .constructor()
      .derives<MCMC>("McmcHousepi");
}
