/* housepi - Epidemiological Analysis of Household Outbreaks.
 *
 * Copyright © 2026 Institut Pasteur (Paris), CNRS and INSERM.
 * Distributed under the terms of the GNU General Public License (GPLv3).
 * See the LICENSE file for details.
 */

#include <memory>
#include <string>
#include <vector>

#include "tools/flags.h"
#include "tools/global.h"
#include "tools/messages.h"
#include "tools/nan.h"

#include "function.h"
#include "distribution.h"
#include "registry.h"
#include "mcmc.h"
#include "model.h"

#include <Rcpp.h>

RCPP_MODULE(housepi_flags)
{
    Rcpp::class_<HousepiRcppCompiledFlags>("HousepiRcppCompiledFlags")
        .constructor<>()
        .method("get_flags", &HousepiRcppCompiledFlags::get_flags);
}

// macro to define Rcpp wrapper classes to export to R
#define MAKE_HOUSEPI_RCPP_WRAPPER(Wrapper, family, set_eval, pdf_call, cdf_call, ...)              \
    class HousepiRcpp##Wrapper                                                                     \
    {                                                                                              \
    public:                                                                                        \
        HousepiRcpp##Wrapper(std::string name_, SEXP args_, double min, double max, bool is_discr) \
        {                                                                                          \
            try                                                                                    \
            {                                                                                      \
                ptr_ = compose_func_registries(name_, family, {__VA_ARGS__}, false);               \
                if (args_ != R_NilValue)                                                           \
                {                                                                                  \
                    Rcpp::List args(args_);                                                        \
                    Rcpp::CharacterVector arg_names;                                               \
                    if (args.hasAttribute("names"))                                                \
                        arg_names = args.names();                                                  \
                    else                                                                           \
                        arg_names = Rcpp::CharacterVector(args.size());                            \
                    for (R_xlen_t i = 0; i < args.size(); ++i)                                     \
                    {                                                                              \
                        const auto this_arg = Rcpp::as<std::string>(arg_names[i]);                 \
                        ptr_->init_attr(this_arg, Attribute(Rcpp::as<double>(args[i])));           \
                    }                                                                              \
                }                                                                                  \
                                                                                                   \
                if (!hsp::na::is(min))                                                             \
                    ptr_->init_attr(hsp::params::MIN, Attribute(min));                             \
                if (!hsp::na::is(max))                                                             \
                    ptr_->init_attr(hsp::params::MAX, Attribute(max));                             \
                                                                                                   \
                if (is_discr)                                                                      \
                    ptr_->init_discretization();                                                   \
                ptr_->init_truncation();                                                           \
                ptr_->assert_is_callable_explicitly();                                             \
            }                                                                                      \
            catch (const std::exception &e)                                                        \
            {                                                                                      \
                hsp::msg::error(                                                                   \
                    "While initializing " family " '" + name_ + "':\n" + e.what());                \
            }                                                                                      \
        }                                                                                          \
        ~HousepiRcpp##Wrapper() = default;                                                         \
        Rcpp::NumericVector pdf_call(Rcpp::NumericVector x) const                                  \
        {                                                                                          \
            ptr_->set_eval##_call();                                                               \
            return call_(x);                                                                       \
        }                                                                                          \
        Rcpp::NumericVector cdf_call(Rcpp::NumericVector x) const                                  \
        {                                                                                          \
            ptr_->set_eval##_primitive();                                                          \
            return call_(x);                                                                       \
        }                                                                                          \
        Rcpp::NumericVector args() const                                                           \
        {                                                                                          \
            Rcpp::NumericVector output(ptr_->n_args());                                            \
            Rcpp::CharacterVector names(ptr_->n_args());                                           \
            for (size_t i = 0; i < ptr_->n_args(); ++i)                                            \
            {                                                                                      \
                output[i] = ptr_->arg(i);                                                          \
                names[i] = ptr_->arg_name(i);                                                      \
            }                                                                                      \
            output.names() = std::move(names);                                                     \
            return output;                                                                         \
        }                                                                                          \
                                                                                                   \
    private:                                                                                       \
        std::unique_ptr<Function> ptr_;                                                            \
                                                                                                   \
        Rcpp::NumericVector call_(Rcpp::NumericVector x) const                                     \
        {                                                                                          \
            try                                                                                    \
            {                                                                                      \
                const R_xlen_t n = x.size();                                                       \
                Rcpp::NumericVector y(n);                                                          \
                for (R_xlen_t i = 0; i < n; ++i)                                                   \
                    y[i] = ptr_->evaluate(x[i]);                                                   \
                return y;                                                                          \
            }                                                                                      \
            catch (const std::exception &e)                                                        \
            {                                                                                      \
                hsp::msg::error(                                                                   \
                    "While calling " family " '" + ptr_->name() + "':\n" + e.what());              \
            }                                                                                      \
        }                                                                                          \
    };

//' @name HousepiRcppFunction
//' @title Rcpp Function from Housepi
//' @description Call housepi registered function.
//' @field initialize\\cr \itemize{
//'   \item Parameter: name `character(1)` with function name.
//'   \item Parameter: args A named `numeric()` with arguments.
//'   \item Parameter: min `numeric(1)`\\cr Minimum value of the function.
//'   \item Parameter: max `numeric(1)`\\cr Maximum value of the function.
//'   \item Parameter: is_discr `logical(1)`\\cr Wheter discretized.
//' }
//' @field func\\cr \itemize{
//'   \item Parameter: x `numeric()`\\cr Input vector.
//' }
//' @return `numeric()`
//' @field prim\\cr \itemize{
//'   \item Parameter: x `numeric()`\\cr Input vector.
//' }
//' @return `numeric()`
//' @keywords internal
RCPP_MODULE(housepi_function)
{
    MAKE_HOUSEPI_RCPP_WRAPPER(Function, "function",
                              set_eval, func, prim,
                              &UserFuncRegistry::inst(),
                              &TranscendentalRegistry::inst(),
                              &ArithmeticRegistry::inst());

    Rcpp::class_<HousepiRcppFunction>("HousepiRcppFunction")
        .constructor<std::string, SEXP, double, double, bool>()
        .method("func", &HousepiRcppFunction::func)
        .method("prim", &HousepiRcppFunction::prim)
        .method("args", &HousepiRcppFunction::args);
}

//' @name HousepiRcppDistribution
//' @title Rcpp Distribution from Housepi
//' @description Call housepi registered distribution.
//' @field initialize\\cr \itemize{
//'   \item Parameter: name `character(1)`\\cr Name of the distribution.
//'   \item Parameter: args `numeric()`\\cr Custom arguments.
//'   \item Parameter: min `numeric(1)`\\cr Minimum value of the distribution.
//'   \item Parameter: max `numeric(1)`\\cr Maximum value of the distribution.
//'   \item Parameter: is_discr `logical(1)`\\cr Wheter discretized.
//' }
//' @field pdf\\cr \itemize{
//'   \item Parameter: x `numeric()`\\cr Input vector.
//' }
//' @return `numeric()`
//' @field cdf\\cr \itemize{
//'   \item Parameter: x `numeric()`\\cr Input vector.
//' }
//' @return `numeric()`
//' @keywords internal
RCPP_MODULE(housepi_distribution)
{
    MAKE_HOUSEPI_RCPP_WRAPPER(Distribution, "distribution",
                              set_eval_as_distr, pdf, cdf,
                              &ContinuosDistrRegistry::inst(),
                              &DiscreteDistrRegistry::inst());

    Rcpp::class_<HousepiRcppDistribution>("HousepiRcppDistribution")
        .constructor<std::string, SEXP, double, double, bool>()
        .method("pdf", &HousepiRcppDistribution::pdf)
        .method("cdf", &HousepiRcppDistribution::cdf)
        .method("args", &HousepiRcppDistribution::args);
}

//' @name HousepiRcppInfecProfile
//' @title Rcpp Infectivity Profile from Housepi
//' @description Call housepi registered infectivity profile.
//' @field initialize\\cr \itemize{
//'   \item Parameter: name `character(1)`\\cr Name of the profile.
//'   \item Parameter: args `numeric()`\\cr Custom arguments.
//'   \item Parameter: min `numeric(1)`\\cr Minimum value of the profile.
//'   \item Parameter: max `numeric(1)`\\cr Maximum value of the profile.
//'   \item Parameter: is_discr `logical(1)`\\cr Wheter discretized.
//' }
//' @field inst\\cr \itemize{
//'   \item Parameter: x `numeric()`\\cr Input vector.
//' }
//' @return `numeric()`
//' @field cumul\\cr \itemize{
//'   \item Parameter: x `numeric()`\\cr Input vector.
//' }
//' @return `numeric()`
//' @keywords internal
RCPP_MODULE(housepi_infec_profile)
{
    MAKE_HOUSEPI_RCPP_WRAPPER(InfecProfile, "infectivity profile",
                              set_eval_as_distr, inst, cumul,
                              &InfecProfileRegistry::inst(),
                              &ContinuosDistrRegistry::inst(),
                              &DiscreteDistrRegistry::inst());

    Rcpp::class_<HousepiRcppInfecProfile>("HousepiRcppInfecProfile")
        .constructor<std::string, SEXP, double, double, bool>()
        .method("inst", &HousepiRcppInfecProfile::inst)
        .method("cumul", &HousepiRcppInfecProfile::cumul)
        .method("args", &HousepiRcppInfecProfile::args);
}

//' @name HousepiRcppModel
//' @title Internal housepi Rcpp model
//' @description Perform inference and simulation from the model.
//' @field infer\\cr Perform MCMC inference of the model parameters. \itemize{
//'   \item Parameter: output_file - The file path where the output will be saved
//'   \item Parameter: progbar - Display the progress bar
//' }
//' @field data_wide_for_simulation\\cr Dataset without status and augmented entries.
//' @field simulate\\cr Simulated status and augmented variables.
//' @field simulate_for_prediction\\cr Simulation for the computation of expected measure.
//' @field evaluate\\cr Evaluate the log-likelihood according to the current state.
//' @field log_evidence\\cr Evaluate the logarithm of the marginal likelihood (aka evidence).
//' Uses the method of power posteriors. \itemize{
//'   \item Parameter: scheme
//'   \item Parameter: n_rungs
//' }
//' @field import_setup\\cr Import setup from R List.
//' @field import_data_wide\\cr Import house indexed data in wide format from R Data.Frame.
//' @field import_data_long\\cr Import individual indexed in long format from R Data.Frame.
//' @field import_multiplicity\\cr Imoprt multiplicity definitions from R Data.Frame.
//' @field import_group_rates\\cr Imoprt group rates definitions from R Data.Frame.
//' @field import_param\\cr Import parameter definitions from R Data.Frame.
//' @field import_infec_profile\\cr Import infectivity profile definition from R List.
//' @field import_augm\\cr Import augmented variable definitions from R Data.Frame.
//' @field import_param_for_prediction\\cr Import parameter tables for expected measures estimation.
//' @field build_all_rates\\cr Assign all contact and group rates.
//' @field export_setup\\cr Return a List with current setup.
//' @field export_param_list_const\\cr Return a CharacterVector with constant parameters.
//' @field export_param_list_sample_standard\\cr Return a CharacterVector with parameters to sample.
//' @field export_param_list_sample_cache\\cr Return a CharacterVector with parameters to sample with cached dependecies.
//' @field export_infec_profile\\cr Return a List with infectivity profile.
//' @field export_data_wide\\cr Return a DataFrame with current database.
//' @field summarise_cases\\cr Return a DataFrame with a report of the current database.
//' @field summarise_house_sizes\\cr Return a DataFrame with a report of household sizes.
//' @field flush_total_log_lik_draws\\cr Return a NumericMatrix with the total log-likelihood draws.
//' @field flush_param_draws\\cr Return a NumericMatrix with parameters draws.
//' @field flush_param_lprior_draws\\cr Return a NumericMatrix with parameter with log-priors draws.
//' @field flush_augm_draws\\cr Return a NumericMatrix with augmented draws.
//' @field flush_house_log_lik_draws\\cr Return a NumericMatrix with household log-likelihood draws.
//' @field export_house_snapshot\\cr Return a List with a snapshot containing full household and model details.
//' @keywords internal
RCPP_MODULE(housepi_model_core)
{
    Rcpp::class_<MCMC>("HousepiRcppMcmc")
        .constructor()
        .method("init_advanced_setup", &MCMC::init_advanced_setup)
        .method("infer", &MCMC::infer)
        .method("data_wide_for_simulation", &MCMC::data_wide_for_simulation)
        .method("simulate", &MCMC::simulate)
        .method("simulate_for_prediction", &MCMC::simulate_for_prediction)
        .method("evaluate_log_lik", &MCMC::evaluate_log_lik)
        .method("log_evidence", &MCMC::log_evidence)
        .method("import_data_wide", &MCMC::import_data_wide)
        .method("import_data_long", &MCMC::import_data_long)
        .method("import_multiplicity", &MCMC::import_multiplicity)
        .method("import_group_rate", &MCMC::import_group_rate)
        .method("import_setup", &MCMC::import_setup)
        .method("import_param", &MCMC::import_param)
        .method("import_infec_profile", &MCMC::import_infec_profile)
        .method("import_augm", &MCMC::import_augm)
        .method("import_param_for_prediction", &MCMC::import_param_for_prediction)
        .method("build_all_rates", &MCMC::build_all_rates)
        .method("export_data_wide", &MCMC::export_data_wide)
        .method("export_data_long", &MCMC::export_data_long)
        .method("summarise_cases", &MCMC::summarise_cases)
        .method("summarise_house_sizes", &MCMC::summarise_house_sizes)
        .method("export_setup", &MCMC::export_setup)
        .method("export_param_list_const", &MCMC::export_param_list_const)
        .method("export_param_list_sample_standard", &MCMC::export_param_list_sample_standard)
        .method("export_param_list_sample_cache", &MCMC::export_param_list_sample_cache)
        .method("export_augm_list_draws", &MCMC::export_augm_list_draws)
        .method("export_augm_list_full", &MCMC::export_augm_list_full)
        .method("export_augm_list_export", &MCMC::export_augm_list_export)
        .method("export_augm_list_discrete", &MCMC::export_augm_list_discrete)
        .method("export_augmented_data_wide", &MCMC::export_augmented_data_wide)
        .method("export_infec_profile", &MCMC::export_infec_profile)
        .method("flush_total_log_lik_draws", &MCMC::flush_total_log_lik_draws)
        .method("flush_param_draws", &MCMC::flush_param_draws)
        .method("flush_param_lprior_draws", &MCMC::flush_param_lprior_draws)
        .method("flush_augm_draws", &MCMC::flush_augm_draws)
        .method("flush_house_log_lik_draws", &MCMC::flush_house_log_lik_draws)
        .method("export_house_snapshot", &MCMC::export_house_snapshot);

    Rcpp::class_<Model>("HousepiRcppModel")
        .constructor()
        .derives<MCMC>("HousepiRcppMcmc");
}
