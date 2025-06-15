/*      __
 *     / /_  ____  __  __ ____ __  ____  __
 *    / __ \/ __ \/ / / / ___/ _ \/ __ \/ /
 *   / / / / /_/ / /_/ (__  )  __/ /_/ / /
 *  /_/ /_/\____/\____/____/\___/ ____/_/
 *                             /_/
 *
 * Copyright (C) 2025 Francesco Camaglia, Institut Pasteur
 *
 * export_bindings.cpp
 *
 */

#include <Rcpp.h>
#include "model.h"

//' @name ModelHousepi
//' @title Internal housepi model
//' @description Perform inference and simulation from the model.
//' @field init_user_advanced_setup Load user advance setup.
//' @field infer Perform MCMC inference of the model parameters. \itemize{
//'   \item Parameter: output_file - The file path where the output will be saved
//'   \item Parameter: progbar - show the progress bar
//' }
//' @field simulate Load a simulated database.
//' @field evaluate Evaluate the loglikelihood according to the current state.
//' @field estimate_evidence Evaluate the logarithm of the marginal likelihood with the method of power posteriors. \itemize{
//'   \item Parameter: scheme
//'   \item Parameter: n_rungs
//' }
//' @filed import_setup Import setup from R List.
//' @filed import_dataset Import house indexed dataset from R Data.Frame.
//' @filed import_datalist Import individual indexed datalist from R Data.Frame.
//' @field import_multiplicity Imoprt multiplicity definitions from R Data.Frame.
//' @filed import_param Import parameter definitions from R Data.Frame.
//' @filed import_infec_profile Import infectivity profile definition from R List.
//' @filed import_augm Import augmented variable definitions from R Data.Frame.
//' @filed buil_all_rates Assign all contact and group rates.
//' @field export_setup Return a List with current setup.
//' @field export_param Return a DataFrame with current parameters.
//' @field export_dataset Return a DataFrame with current database.
//' @field export_data_report Return a DataFrame with a report of the current database.
//' @field export_house_size_report Return a DataFrame with a report of household sizes.
//' @field export_traces Return a DataFrame with MCMC traces.
//' @field export_house_loglike_mtx Return a NumericMatrix with household loglikelihood traces.
//' @field showcase_infec_profile Return a List with a showcase of the infectivity profile.
//' @field show_full_build Show all build messages.
//' @field show_summary Show a summary of the current database.

RCPP_MODULE(housepi_model_core)
{
    using namespace Rcpp;

    class_<MCMC>("McmcHousepi")
        .constructor()
        .method("init_user_advanced_setup", &MCMC::init_user_advanced_setup)
        .method("infer", &MCMC::infer)
        .method("simulate", &MCMC::simulate)
        .method("evaluate", &MCMC::evaluate)
        .method("estimate_evidence", &MCMC::estimate_evidence)
        .method("import_setup", &MCMC::import_setup)
        .method("import_dataset", &MCMC::import_dataset)
        .method("import_datalist", &MCMC::import_datalist)
        .method("import_multiplicity", &MCMC::import_multiplicity)
        .method("import_param", &MCMC::import_param)
        .method("import_infec_profile", &MCMC::import_infec_profile)
        .method("import_augm", &MCMC::import_augm)
        .method("build_all_rates", &MCMC::build_all_rates)
        .method("export_setup", &MCMC::export_setup)
        .method("export_param", &MCMC::export_param)
        .method("export_dataset", &MCMC::export_dataset)
        .method("export_data_report", &MCMC::export_data_report)
        .method("export_house_size_report", &MCMC::export_house_size_report)
        .method("export_traces", &MCMC::export_traces)
        .method("export_house_loglike_mtx", &MCMC::export_house_loglike_mtx)
        .method("showcase_infec_profile", &MCMC::showcase_infec_profile);

    class_<Model>("ModelHousepi")
        .constructor()
        .derives<MCMC>("McmcHousepi");
}
