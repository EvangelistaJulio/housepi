/* housepi - Epidemiological Analysis of Household Outbreaks.
 *
 * Copyright © 2026 Institut Pasteur (Paris), CNRS and INSERM.
 * Distributed under the terms of the GNU General Public License (GPLv3).
 * See the LICENSE file for details.
 */

#pragma once

#include <string>
#include <vector>

#include "tools/keys.h"

#include "individual.h"
#include "household.h"

#ifdef _HOUSEPI_RCPP
#include <Rcpp.h>
#else
class MCMC;
int standalone(MCMC &model, int argc, char *argv[]);
#endif

#ifndef _HOUSEPI_MODULE
#define _HOUSEPI_MODULE blank
#endif

#define HOUSEPI_DECLARE_EMPTY_KEYS(T)                        \
    namespace                                                \
    {                                                        \
        static const std::vector<std::string> T##_keys = {}; \
    }

#define HOUSEPI_DECLARE_REGISTER_KEYS(T, ...)                                                        \
    namespace                                                                                        \
    {                                                                                                \
        struct T                                                                                     \
        {                                                                                            \
            enum : size_t                                                                            \
            {                                                                                        \
                __VA_ARGS__                                                                          \
            };                                                                                       \
        };                                                                                           \
                                                                                                     \
        static const std::vector<std::string> T##_keys = hsp::key::helper_parse_names(#__VA_ARGS__); \
    }

#define HOUSEPI_MODEL_CONSTRUCTOR()     \
    Model() : MCMC(DataKey_keys,        \
                   HHsizeKey_keys,      \
                   HHsizevecKey_keys,   \
                   HHdoublevecKey_keys, \
                   HHdoublemtxKey_keys, \
                   LongKey_keys,        \
                   ParamKey_keys,       \
                   GroupRateKey_keys,   \
                   AugmKey_keys)        \
    {                                   \
        init_advanced_setup();          \
    }

#define HOUSEPI_OPEN_NAMED_MODULE_NAMESPACE() \
    namespace _HOUSEPI_MODULE                 \
    {
#define HOUSEPI_CLOSE_NAMED_MODULE_NAMESPACE() }

//' @name housepi_compile_flags
//' @title Get Housepi Compilation Flags
//' @description Return the compilation flags used when building the housepi library.
//' @return `character(1)` A string containing the compiled flags.
//' @keywords internal

/// Expose compile flags.
std::string housepi_compile_flags()
{
    return
#ifdef _HOUSEPI_DISCRETE_TIME
        "_HOUSEPI_DISCRETE_TIME:";
#else
        "";
#endif
}

//' @name HousepiRcppModel
//' @title Internal housepi Rcpp model
//' @description Perform inference and simulation from the model.
//' @field run_inference Perform MCMC inference of the model parameters.
//' @field data_wide_for_simulation Dataset without status and augmented entries.
//' @field run_simulation Simulated status and augmented variables.
//' @field run_simulation_for_predict Simulation for the computation of expected measure.
//' @field evaluate Evaluate the log-likelihood according to the current state.
//' @field log_evidence Evaluate the logarithm of the marginal likelihood (aka evidence).
//' Uses the method of power posteriors. \itemize{
//'   \item Parameter: scheme
//'   \item Parameter: n_rungs
//' }
//' @field import_setup Import setup from R List.
//' @field import_data_wide Import house indexed data in wide format from R Data.Frame.
//' @field import_data_long Import individual indexed in long format from R Data.Frame.
//' @field import_multiplicity Imoprt multiplicity definitions from R Data.Frame.
//' @field import_group_rates Imoprt group rates definitions from R Data.Frame.
//' @field import_param Import parameter definitions from R Data.Frame.
//' @field import_infec_profile Import infectivity profile definition from R List.
//' @field import_augm Import augmented variable definitions from R Data.Frame.
//' @field import_param_for_prediction Import parameter tables for expected measures estimation.
//' @field build_contact_rates Assign all contact rates.
//' @field export_setup Return a List with current setup.
//' @field export_infec_profile Return a List with infectivity profile.
//' @field export_data_wide Return a DataFrame with current database.
//' @field summarise_cases Return a DataFrame with a report of the current database.
//' @field summarise_house_sizes Return a DataFrame with a report of household sizes.
//' @field flush_total_log_lik_draws Return a NumericMatrix with the total log-likelihood draws.
//' @field flush_param_draws Return a NumericMatrix with parameters draws.
//' @field flush_param_lprior_draws Return a NumericMatrix with parameter with log-priors draws.
//' @field flush_augm_draws Return a NumericMatrix with augmented draws.
//' @field flush_house_log_lik_draws Return a NumericMatrix with household log-likelihood draws.
//' @field export_house_snapshot Return a List with a snapshot containing full household and model details.
//' @keywords internal

#define STRINGIFY(x) #x
#define TOSTRING(x) STRINGIFY(x)

#ifdef _HOUSEPI_RCPP
#define HOUSEPI_REGISTER_RCPP_MODULE_INTERFACE_IMPL_(module_name)                      \
    RCPP_MODULE(module_name)                                                           \
    {                                                                                  \
        using namespace module_name;                                                   \
                                                                                       \
        Rcpp::function("housepi_compile_flags", &housepi_compile_flags);               \
                                                                                       \
        Rcpp::class_<MCMC>("MCMC." TOSTRING(module_name))                              \
            .method("run_inference", &MCMC::run_inference)                             \
            .method("data_wide_for_simulation", &MCMC::data_wide_for_simulation)       \
            .method("run_simulation", &MCMC::run_simulation)                           \
            .method("run_simulation_for_predict", &MCMC::run_simulation_for_predict)   \
            .method("evaluate_log_lik", &MCMC::evaluate_log_lik)                       \
            .method("log_evidence", &MCMC::log_evidence)                               \
            .method("import_data_wide", &MCMC::import_data_wide)                       \
            .method("import_data_long", &MCMC::import_data_long)                       \
            .method("import_multiplicity", &MCMC::import_multiplicity)                 \
            .method("import_group_rate", &MCMC::import_group_rate)                     \
            .method("import_setup", &MCMC::import_setup)                               \
            .method("import_param", &MCMC::import_param)                               \
            .method("import_infec_profile", &MCMC::import_infec_profile)               \
            .method("import_augm", &MCMC::import_augm)                                 \
            .method("import_param_for_prediction", &MCMC::import_param_for_prediction) \
            .method("build_contact_rates", &MCMC::build_contact_rates)                 \
            .method("export_data_wide", &MCMC::export_data_wide)                       \
            .method("export_data_long", &MCMC::export_data_long)                       \
            .method("summarise_cases", &MCMC::summarise_cases)                         \
            .method("summarise_house_sizes", &MCMC::summarise_house_sizes)             \
            .method("export_setup", &MCMC::export_setup)                               \
            .method("export_augmented_data_wide", &MCMC::export_augmented_data_wide)   \
            .method("export_infec_profile", &MCMC::export_infec_profile)               \
            .method("flush_total_log_lik_draws", &MCMC::flush_total_log_lik_draws)     \
            .method("flush_param_draws", &MCMC::flush_param_draws)                     \
            .method("flush_param_lprior_draws", &MCMC::flush_param_lprior_draws)       \
            .method("flush_augm_draws", &MCMC::flush_augm_draws)                       \
            .method("flush_house_log_lik_draws", &MCMC::flush_house_log_lik_draws)     \
            .method("export_house_snapshot", &MCMC::export_house_snapshot);            \
                                                                                       \
        Rcpp::class_<Model>("Model." TOSTRING(module_name))                            \
            .constructor()                                                             \
            .derives<MCMC>("MCMC." TOSTRING(module_name));                             \
    }
#define HOUSEPI_REGISTER_RCPP_MODULE_INTERFACE() \
    HOUSEPI_REGISTER_RCPP_MODULE_INTERFACE_IMPL_(_HOUSEPI_MODULE)
#else
#define HOUSEPI_REGISTER_RCPP_MODULE_INTERFACE()   \
    int main(int argc, char *argv[])               \
    {                                              \
        using namespace _HOUSEPI_MODULE;           \
        Model user_model;                          \
        return standalone(user_model, argc, argv); \
    }
#endif
