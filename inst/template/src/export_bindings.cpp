/*
 * export_bin.cpp
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
//' @field setup Return the current setup as a List.
//' @field params Return the current parameters as a DataFrame.
//' @field data Return the current database as a DataFrame.
//' @field show_setup Show the setup configurations.
//' @field show_summary Show a summary of the database.
//' @field show_group_classes Show the group classes.
//' @field show_params Show the parameter configurations.
//' @field show_diagnostics Show the WAIC.
//' @export

#include <Rcpp.h>
#include "model.h"

// Expose the Model and MCMC classes to R using Rcpp
RCPP_MODULE(modelhousepi12345_cpp)
{
 using namespace Rcpp;
 
 class_<MCMC>("MCMChousepi12345")
   .constructor()
   .method("build", &MCMC::build)
   .method("infer", &MCMC::infer)
   .method("simulate", &MCMC::simulate)
   .method("evaluate", &MCMC::evaluate)
   .method("evaluate_marginal", &MCMC::evaluate_marginal)
   .method("setup", &MCMC::export_setup)
   .method("params", &MCMC::export_params)
   .method("data", &MCMC::export_data)
   .method("show_setup", &MCMC::show_setup)
   .method("show_summary", &MCMC::show_summary)
   .method("show_group_classes", &MCMC::show_group_classes)
   .method("show_params", &MCMC::show_parameters)
   .method("show_diagnostics", &MCMC::show_diagnostics);
 
 class_<Model>("Modelhousepi12345")
   .constructor()
   .derives<MCMC>("MCMChousepi12345");
}
