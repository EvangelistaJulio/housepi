//' @name Modelhousepi12345
//' @title Custom model
//' @description Perform inference and simulation from the model housepi12345.
//' @field new Initialize a household model object.
//' @field build Construct and set up the model \itemize{
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
//' @field infer Perform inference on the model \itemize{
//'   \item Parameter: output_file - The file path where the output will be saved (temporary)
//'   \item Parameter: progbar - boolean flag to suppress automatic prints
//' }
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
