//' @name Model_hemirt12345
//' @title Custom model
//' @description Type 'Model' to see its methods
//' @field new Constructor to initialize a Model object
//' @field build Construct and set up the model \itemize{
//'   \item Parameter: is_infer - boolean flag to indicate whether to perform inference
//'   \item Parameter: is_simulate - boolean flag to indicate whether to perform simulation
//'   \item Parameter: params_file - file path for the model parameters
//'   \item Parameter: data_file - file path for the input data
//'   \item Parameter: augm_file - file path for augmented data
//'   \item Parameter: setup_file - file path for the setup configurations
//'   \item Parameter: extra_file - file path for any extra configurations
//' }
//' @field inference Perform inference on the model \itemize{
//'   \item Parameter: output_file - The file path where the output will be saved
//'   \item Parameter: verbose - boolean flag to suppress automatic prints
//' }
//' @field simulation Perform simulation on the model \itemize{
//'   \item Parameter: output_file - The file path where the output will be saved
//'   \item Parameter: verbose - boolean flag to suppress automatic prints
//' }
//' @field loglikelihood Compute model loglikelihood \itemize{
//'   \item Parameter: verbose - boolean flag to suppress automatic prints
//' }
//' @field display.data Summary of data and model group rates.
//' @field display.results Results of the inference.
//' @export

#include <Rcpp.h>
#include "model.h"
 
// Expose the Model and MCMC classes to R using Rcpp
RCPP_MODULE(model_hemirt12345)
{
 using namespace Rcpp;
 
 class_<MCMC>("MCMC_hemirt12345")
   .constructor("Default constructor")
   .method("build", &MCMC::build, "Construct the model with various configuration files and flags.")
   .method("inference", &MCMC::run_inference, "Perform inference and save results to file.")
   .method("simulation", &MCMC::run_simulation, "Perform simulation and save results to file.")
   .method("loglikelihood", &MCMC::run_compute_loglikelihood, "Perform simulation and save results to file.")
   .method("setup", &MCMC::display_setup, "Display setup.")
   .method("summary", &MCMC::display_summary, "Display summary of the dataset.")
   .method("report", &MCMC::display_report, "Report of model group rates.")
   .method("params", &MCMC::display_parameters, "Display parameters.")
   .method("results", &MCMC::display_thin_results, "Results of the inference.");
 
 class_<Model>("Model_hemirt12345")
   .constructor("Default constructor")
   .derives<MCMC>("MCMC_hemirt12345");
}
 