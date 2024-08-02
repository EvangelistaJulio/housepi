//' @name Modelhousepi12345
//' @title Custom model
//' @description Type 'Model' to see its methods
//' @field new Constructor to initialize a Model object
//' @field build Construct and set up the model \itemize{
//'   \item Parameter: params_file - file path for the model parameters
//'   \item Parameter: data_file - file path for the input data
//'   \item Parameter: augm_file - file path for augmented data
//'   \item Parameter: setup_file - file path for the setup configurations
//'   \item Parameter: extra_file - file path for any extra configurations
//'   \item Parameter: verbose - boolean flag to suppress automatic prints
//' }
//' @field infer Perform inference on the model \itemize{
//'   \item Parameter: output_file - The file path where the output will be saved (temporary)
//'   \item Parameter: progbar - boolean flag to suppress automatic prints
//' }
//' @field simulate Change the current database with a simulation from the model.
//' @field loglikelihood Evaluate the loglikelihood in the current state.
//' @field display.data Summary of data and model group rates.
//' @field display.results Results of the inference.
//' @export

#include <Rcpp.h>
#include "model.h"
 
// Expose the Model and MCMC classes to R using Rcpp
RCPP_EXPOSED_CLASS(Modelhousepi12345)
RCPP_MODULE(modelhousepi12345)
{
 using namespace Rcpp;
 
 class_<MCMC>("MCMChousepi12345")
   .constructor()
   .method("build", &MCMC::build)
   .method("infer", &MCMC::infer)
   .method("simulate", &MCMC::simulate)
   .method("evaluate", &MCMC::evaluate)
   .method("setup", &MCMC::return_setup)
   .method("show_setup", &MCMC::show_setup)
   .method("show_summary", &MCMC::show_summary)
   .method("show_group_classes", &MCMC::show_group_classes)
   .method("show_params", &MCMC::show_parameters)
   .method("show_results", &MCMC::show_results);
 
 class_<Model>("Modelhousepi12345")
   .constructor("Default constructor")
   .derives<MCMC>("MCMChousepi12345");
}
 