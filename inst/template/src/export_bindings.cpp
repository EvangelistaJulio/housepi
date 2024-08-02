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
//'   \item Parameter: output_file - The file path where the output will be saved
//'   \item Parameter: verbose - boolean flag to suppress automatic prints
//' }
//' @field simulate Perform simulation on the model \itemize{
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
RCPP_EXPOSED_CLASS(Modelhousepi12345)
RCPP_MODULE(modelhousepi12345)
{
 using namespace Rcpp;
 
 class_<MCMC>("MCMChousepi12345")
   .constructor()
   .method("build", &MCMC::build)
   .method("infer", &MCMC::run_inference)
   .method("simulate", &MCMC::run_simulation)
   .method("loglikelihood", &MCMC::run_eval_loglikelihood)
   .method("setup", &MCMC::return_setup);
 
 class_<Model>("Modelhousepi12345")
   .constructor("Default constructor")
   .derives<MCMC>("MCMChousepi12345");
}
 