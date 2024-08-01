//' @name Modelhousepi12345
//' @title Household modeling class
//' @description Perform inference and simulation from the model housepi12345.
//' @field new Initialize a household model object.
//' @field build Associate definitions and data to the model. \itemize{
//'   \item Parameter: params_file - file path for the model parameters.
//'   \item Parameter: data_file - file path for the input data.
//'   \item Parameter: augm_file - file path for augmented data.
//'   \item Parameter: setup_file - file path for the setup configurations.
//'   \item Parameter: extra_file - file path for any extra configurations.
//'   \item Parameter: verbose - boolean flag to print additional details.
//' }
//' @field infer Perform inference on the model. \itemize{
//'   \item Parameter: output_file - The file path where the output will be saved.
//'   \item Parameter: verbose - boolean flag to print additional details.
//'   \item Parameter: progbar - boolean flag to show the progress bar.
//' }
//' @field simulate Perform simulation on the model. \itemize{
//'   \item Parameter: output_file - The file path where the output will be saved.
//'   \item Parameter: verbose - boolean flag to print additional details.
//' }
//' @field loglikelihood Compute the model loglikelihood. \itemize{
//'   \item Returns: double.
//' }
//' @field setup The setup parameter list. \itemize{
//'   \item Returns: list.
//' }
//' @field show_setup Show choices of the setup.
//' @field show_params Show choices of the parameters.
//' @field show_summary Show a summary of the dataset.
//' @field show_group_classes Show a report of the group classes.
//' @field show_results Show the results of the inference.
//' @export

#include <Rcpp.h>
#include "model.h"

// Expose the Model and MCMC classes to R using Rcpp
RCPP_EXPOSED_CLASS(Modelhousepi12345)
RCPP_MODULE(modelhousepi12345_cpp)
{
	using namespace Rcpp;

	class_<MCMC>("MCMChousepi12345")
		.constructor()
		.method("build", &MCMC::build)
		.method("infer", &MCMC::run_inference)
		.method("simulate", &MCMC::run_simulation)
		.method("loglikelihood", &MCMC::run_eval_loglikelihood)
		.method("setup", &MCMC::return_setup)
		.method("show_setup", &MCMC::show_setup)
		.method("show_summary", &MCMC::show_summary)
		.method("show_group_classes", &MCMC::show_group_classes)
		.method("show_params", &MCMC::show_parameters)
		.method("show_results", &MCMC::show_results);

	class_<Model>("Modelhousepi12345")
		.constructor("Default constructor.")
		.derives<MCMC>("MCMChousepi12345");
}
