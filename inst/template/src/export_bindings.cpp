/*
 *        __
 *       / /_  ____  __  __ ____ __  ____  __
 *      / __ \/ __ \/ / / / ___/ _ \/ __ \/ /
 *     / / / / /_/ / /_/ (__  )  __/ /_/ / /
 *    /_/ /_/\____/\____/____/\___/ ____/_/
 *                               /_/
 *
 * Copyright (C) 2025 Francesco Camaglia, Institut Pasteur
 *
 * export_bindings.cpp
 *
 */

#include <Rcpp.h>

#include "tools/flags.h"
#include "tools/na_utils.h"
#include "tools/error.h"
#include "tools/params.h"
#include "function.h"
#include "function_registry.h"
#include "distribution_registry.h"
#include "model.h"

RCPP_MODULE(housepi_flags)
{
    Rcpp::class_<CompiledFlagsHousepi>("CompiledFlagsHousepi")
        .constructor<>()
        .method("get_flags", &CompiledFlagsHousepi::get_flags);
}

// macro to define RCPP wrapper classes to export to R
#define MAKE_HOUSEPI_RCPP_WRAPPER(Wrapper, family, set_eval, pdf_call, cdf_call, ...)          \
    class Wrapper##Housepi                                                                     \
    {                                                                                          \
    public:                                                                                    \
        Wrapper##Housepi(std::string name_, SEXP args_, double min, double max, bool is_discr) \
        {                                                                                      \
            try                                                                                \
            {                                                                                  \
                ptr_ = compose_func_registries(name_, family, {__VA_ARGS__}, false);           \
                if (args_ != R_NilValue)                                                       \
                {                                                                              \
                    Rcpp::List args(args_);                                                    \
                    const Rcpp::CharacterVector arg_names = args.names();                      \
                    for (R_xlen_t i = 0; i < args.size(); ++i)                                 \
                    {                                                                          \
                        const std::string this_arg = Rcpp::as<std::string>(arg_names[i]);      \
                        ptr_->init_attr(this_arg, Attribute(Rcpp::as<double>(args[i])));       \
                    }                                                                          \
                }                                                                              \
                                                                                               \
                if (!hsp::na::is(min))                                                         \
                    ptr_->init_attr(hsp::params::MIN, Attribute(min));                         \
                if (!hsp::na::is(max))                                                         \
                    ptr_->init_attr(hsp::params::MAX, Attribute(max));                         \
                                                                                               \
                if (is_discr)                                                                  \
                    ptr_->init_discretization();                                               \
                ptr_->init_truncation();                                                       \
            }                                                                                  \
            catch (const std::exception &e)                                                    \
            {                                                                                  \
                hsp::msg::error(                                                               \
                    "While initializing " family " '" + name_ + "':\n" + e.what());            \
            }                                                                                  \
        }                                                                                      \
        ~Wrapper##Housepi() = default;                                                         \
        Rcpp::NumericVector pdf_call(Rcpp::NumericVector x) const                              \
        {                                                                                      \
            ptr_->set_eval##_call();                                                           \
            return call_(x);                                                                   \
        }                                                                                      \
        Rcpp::NumericVector cdf_call(Rcpp::NumericVector x) const                              \
        {                                                                                      \
            ptr_->set_eval##_primitive();                                                      \
            return call_(x);                                                                   \
        }                                                                                      \
                                                                                               \
    private:                                                                                   \
        std::unique_ptr<Function> ptr_;                                                        \
                                                                                               \
        Rcpp::NumericVector call_(Rcpp::NumericVector x) const                                 \
        {                                                                                      \
            try                                                                                \
            {                                                                                  \
                const R_xlen_t n = x.size();                                                   \
                Rcpp::NumericVector y(n);                                                      \
                for (R_xlen_t i = 0; i < n; ++i)                                               \
                    y[i] = ptr_->evaluate(x[i]);                                               \
                return y;                                                                      \
            }                                                                                  \
            catch (const std::exception &e)                                                    \
            {                                                                                  \
                hsp::msg::error(                                                               \
                    "While calling '" family "'" + ptr_->name() + "':\n" + e.what());          \
            }                                                                                  \
        }                                                                                      \
    };

//' @name FunctionHousepi
//' @title Function from Housepi
//' @description Call housepi registered function.
//' @field initialize\cr \itemize{
//'   \item Parameter: name `character(1)` with function name.
//'   \item Parameter: args A named `numeric()` with arguments.
//'   \item Parameter: min `numeric(1)`\cr Minimum value of the function.
//'   \item Parameter: max `numeric(1)`\cr Maximum value of the function.
//'   \item Parameter: is_discr `logical(1)`\cr Wheter discretized.
//' }
//' @field func\cr \itemize{
//'   \item Parameter: x `numeric()`\cr Input vector.
//' }
//' @return `numeric()`
//' @field prim\cr \itemize{
//'   \item Parameter: x `numeric()`\cr Input vector.
//' }
//' @return `numeric()`
RCPP_MODULE(housepi_function)
{
    MAKE_HOUSEPI_RCPP_WRAPPER(Function, "function",
                              set_eval, func, prim,
                              &UserFuncRegistry::inst(),
                              &TranscendentalRegistry::inst(),
                              &ArithmeticRegistry::inst());

    Rcpp::class_<FunctionHousepi>("FunctionHousepi")
        .constructor<std::string, SEXP, double, double, bool>()
        .method("func", &FunctionHousepi::func)
        .method("prim", &FunctionHousepi::prim);
}

//' @name DistributionHousepi
//' @title Distribution from Housepi
//' @description Call housepi registered distribution.
//' @field initialize\cr \itemize{
//'   \item Parameter: name `character(1)`\cr Name of the distribution.
//'   \item Parameter: args `numeric()`\cr Custom arguments.
//'   \item Parameter: min `numeric(1)`\cr Minimum value of the distribution.
//'   \item Parameter: max `numeric(1)`\cr Maximum value of the distribution.
//'   \item Parameter: is_discr `logical(1)`\cr Wheter discretized.
//' }
//' @field pdf\cr \itemize{
//'   \item Parameter: x `numeric()`\cr Input vector.
//' }
//' @return `numeric()`
//' @field cdf\cr \itemize{
//'   \item Parameter: x `numeric()`\cr Input vector.
//' }
//' @return `numeric()`
RCPP_MODULE(housepi_distribution)
{
    MAKE_HOUSEPI_RCPP_WRAPPER(Distribution, "distribution",
                              set_eval_as_distr, pdf, cdf,
                              &ContinuosDistrRegistry::inst(),
                              &DiscreteDistrRegistry::inst());

    Rcpp::class_<DistributionHousepi>("DistributionHousepi")
        .constructor<std::string, SEXP, double, double, bool>()
        .method("pdf", &DistributionHousepi::pdf)
        .method("cdf", &DistributionHousepi::cdf);
}

//' @name InfecProfileHousepi
//' @title Infectivity Profile from Housepi
//' @description Call housepi registered infectivity profile.
//' @field initialize\cr \itemize{
//'   \item Parameter: name `character(1)`\cr Name of the profile.
//'   \item Parameter: args `numeric()`\cr Custom arguments.
//'   \item Parameter: min `numeric(1)`\cr Minimum value of the profile.
//'   \item Parameter: max `numeric(1)`\cr Maximum value of the profile.
//'   \item Parameter: is_discr `logical(1)`\cr Wheter discretized.
//' }
//' @field pdf\cr \itemize{
//'   \item Parameter: x `numeric()`\cr Input vector.
//' }
//' @return `numeric()`
//' @field cdf\cr \itemize{
//'   \item Parameter: x `numeric()`\cr Input vector.
//' }
//' @return `numeric()`
RCPP_MODULE(housepi_infec_profile)
{
    MAKE_HOUSEPI_RCPP_WRAPPER(InfecProfile, "infectivity profile",
                              set_eval_as_distr, pdf, cdf,
                              &InfecProfileRegistry::inst(),
                              &ContinuosDistrRegistry::inst(),
                              &DiscreteDistrRegistry::inst());

    Rcpp::class_<InfecProfileHousepi>("InfecProfileHousepi")
        .constructor<std::string, SEXP, double, double, bool>()
        .method("pdf", &InfecProfileHousepi::pdf)
        .method("cdf", &InfecProfileHousepi::cdf);
}

//' @name ModelHousepi
//' @title Internal housepi model
//' @description Perform inference and simulation from the model.
//' @field infer\cr Perform MCMC inference of the model parameters. \itemize{
//'   \item Parameter: output_file - The file path where the output will be saved
//'   \item Parameter: progbar - Display the progress bar
//' }
//' @field data_wide_for_simulation\cr Dataset without status and augmented entries.
//' @field simulate\cr Simulated status and augmented variables.
//' @field evaluate\cr Evaluate the log-likelihood according to the current state.
//' @field log_evidence\cr Evaluate the logarithm of the marginal likelihood (aka evidence).
//' Uses the method of power posteriors. \itemize{
//'   \item Parameter: scheme
//'   \item Parameter: n_rungs
//' }
//' @field import_setup\cr Import setup from R List.
//' @field import_data_wide\cr Import house indexed data in wide format from R Data.Frame.
//' @field import_data_long\cr Import individual indexed in long format from R Data.Frame.
//' @field import_multiplicity\cr Imoprt multiplicity definitions from R Data.Frame.
//' @field import_group_rates\cr Imoprt group rates definitions from R Data.Frame.
//' @field import_param\cr Import parameter definitions from R Data.Frame.
//' @field import_infec_profile\cr Import infectivity profile definition from R List.
//' @field import_augm\cr Import augmented variable definitions from R Data.Frame.
//' @field buil_all_rates\cr Assign all contact and group rates.
//' @field export_setup\cr Return a List with current setup.
//' @field export_param_list_const\cr Return a CharacterVector with constant parameters.
//' @field export_param_list_sample_standard\cr Return a CharacterVector with parameters to sample.
//' @field export_param_list_sample_cache\cr Return a CharacterVector with parameters to sample with cached dependecies.
//' @field export_infec_profile\cr Return a List with infectivity profile.
//' @field export_data_wide\cr Return a DataFrame with current database.
//' @field summarise_cases\cr Return a DataFrame with a report of the current database.
//' @field summarise_house_sizes\cr Return a DataFrame with a report of household sizes.
//' @field export_draws\cr Return a NumericMatrix with MCMC draws.
//' @field export_house_log_lik_draws\cr Return a NumericMatrix with household log-likelihood draws.
//' @field export_house_snapshot\cr Return a List with a snapshot containing full household and model details.

RCPP_MODULE(housepi_model_core)
{
    Rcpp::class_<MCMC>("McmcHousepi")
        .constructor()
        .method("init_advanced_setup", &MCMC::init_advanced_setup)
        .method("clear_memory", &MCMC::clear_memory)
        .method("infer", &MCMC::infer)
        .method("data_wide_for_simulation", &MCMC::data_wide_for_simulation)
        .method("simulate", &MCMC::simulate)
        .method("evaluate_log_lik", &MCMC::evaluate_log_lik)
        .method("log_evidence", &MCMC::log_evidence)
        .method("import_data_wide", &MCMC::import_data_wide)
        .method("import_data_long", &MCMC::import_data_long)
        .method("import_multiplicity", &MCMC::import_multiplicity)
        .method("import_group_rates", &MCMC::import_group_rates)
        .method("import_setup", &MCMC::import_setup)
        .method("import_param", &MCMC::import_param)
        .method("import_infec_profile", &MCMC::import_infec_profile)
        .method("import_augm", &MCMC::import_augm)
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
        .method("export_draws", &MCMC::export_draws)
        .method("export_house_log_lik_draws", &MCMC::export_house_log_lik_draws)
        .method("export_house_snapshot", &MCMC::export_house_snapshot);

    Rcpp::class_<Model>("ModelHousepi")
        .constructor()
        .derives<MCMC>("McmcHousepi");
}
