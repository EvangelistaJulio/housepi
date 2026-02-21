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

#define HOUSEPI_DEFINE_MODEL_CONSTRUCTOR()     \
    Model::Model() : MCMC(DataKey_keys,        \
                          HHsizeKey_keys,      \
                          HHsizevecKey_keys,   \
                          HHdoublevecKey_keys, \
                          HHdoublemtxKey_keys, \
                          LongKey_keys,        \
                          ParamKey_keys,       \
                          GroupRateKey_keys,   \
                          AugmKey_keys)        \
    {                                          \
        init_advanced_setup();                 \
    }

#ifdef _HOUSEPI_RCPP
#define HOUSEPI_DECLARE_RCPP_MODULE(housepi_module_name)                               \
    RCPP_MODULE(housepi_module_name)                                                   \
    {                                                                                  \
        Rcpp::function("housepi_compile_flags", &housepi_compile_flags);               \
                                                                                       \
        Rcpp::class_<ExposeBase>("ExposeBase")                                         \
            .method("args", &ExposeBase::args);                                        \
                                                                                       \
        Rcpp::class_<ExposeFunction>("ExposeFunction")                                 \
            .constructor<std::string, SEXP, double, double, bool>()                    \
            .derives<ExposeBase>("ExposeBase")                                         \
            .method("func", &ExposeFunction::func)                                     \
            .method("prim", &ExposeFunction::prim);                                    \
                                                                                       \
        Rcpp::class_<ExposeDistribution>("ExposeDistribution")                         \
            .constructor<std::string, SEXP, double, double, bool>()                    \
            .derives<ExposeBase>("ExposeBase")                                         \
            .method("pdf", &ExposeDistribution::pdf)                                   \
            .method("cdf", &ExposeDistribution::cdf);                                  \
                                                                                       \
        Rcpp::class_<ExposeInfecProfile>("ExposeInfecProfile")                         \
            .constructor<std::string, SEXP, double, double, bool>()                    \
            .derives<ExposeBase>("ExposeBase")                                         \
            .method("inst", &ExposeInfecProfile::inst)                                 \
            .method("cumul", &ExposeInfecProfile::cumul);                              \
                                                                                       \
        Rcpp::class_<MCMC>("MCMC")                                                     \
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
        Rcpp::class_<Model>("Model")                                                   \
            .constructor()                                                             \
            .derives<MCMC>("MCMC");                                                    \
    }
#else
#define HOUSEPI_DECLARE_RCPP_MODULE(housepi_module_name)
#endif
