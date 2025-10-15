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

#include "tools/na_utils.h"
#include "tools/error.h"
#include "model.h"

// macro to define RCPP wrapper classes to export to R
#define MAKE_HOUSEPI_RCPP_WRAPPER_CLASS(WrapperName, type_string, registry_list)                       \
    class WrapperName##Housepi                                                                         \
    {                                                                                                  \
    private:                                                                                           \
        std::shared_ptr<Function> ptr_;                                                                \
                                                                                                       \
    public:                                                                                            \
        WrapperName##Housepi(std::string name_, bool is_cumul, SEXP args_, SEXP range_, bool is_discr) \
        {                                                                                              \
            try                                                                                        \
            {                                                                                          \
                ptr_ = compose_func_registries(name_, type_string, registry_list, is_cumul);           \
                                                                                                       \
                if (args_ != R_NilValue)                                                               \
                {                                                                                      \
                    Rcpp::List args(args_);                                                            \
                    const Rcpp::CharacterVector arg_names = args.names();                              \
                    for (R_xlen_t i = 0; i < args.size(); ++i)                                         \
                    {                                                                                  \
                        const std::string this_arg = Rcpp::as<std::string>(arg_names[i]);              \
                        ptr_->init_attr(this_arg, Rcpp::as<double>(args[i]));                          \
                    }                                                                                  \
                }                                                                                      \
                                                                                                       \
                if (range_ != R_NilValue)                                                              \
                {                                                                                      \
                    Rcpp::List range(range_);                                                          \
                    const Rcpp::CharacterVector range_names = range.names();                           \
                    for (R_xlen_t i = 0; i < range.size(); ++i)                                        \
                        if (!hsp::na::is(Rcpp::as<double>(range[i])))                                 \
                            ptr_->init_attr(Rcpp::as<std::string>(range_names[i]),                     \
                                            Rcpp::as<double>(range[i]));                               \
                }                                                                                      \
                                                                                                       \
                if (is_discr)                                                                          \
                    ptr_->init_discretization();                                                       \
                ptr_->init_truncation();                                                               \
            }                                                                                          \
            catch (const std::exception &e)                                                            \
            {                                                                                          \
                hsp::msg::error(                                                                          \
                    "While initializing " type_string " \"" + name_ + "\":\n" + e.what());             \
            }                                                                                          \
        }                                                                                              \
        ~WrapperName##Housepi() = default;                                                             \
                                                                                                       \
        Rcpp::NumericVector call(Rcpp::NumericVector x)                                                \
        {                                                                                              \
            try                                                                                        \
            {                                                                                          \
                const R_xlen_t n = x.size();                                                           \
                Rcpp::NumericVector y(n);                                                              \
                for (R_xlen_t i = 0; i < n; ++i)                                                       \
                    y[i] = ptr_->evaluate(x[i]);                                                       \
                return y;                                                                              \
            }                                                                                          \
            catch (const std::exception &e)                                                            \
            {                                                                                          \
                hsp::msg::error(                                                                          \
                    "While calling " type_string " \"" + ptr_->name() + "\":\n" + e.what());           \
            }                                                                                          \
            return Rcpp::NumericVector();                                                              \
        }                                                                                              \
    };

//' @name FunctionHousepi
//' @title Function from Housepi
//' @description Call housepi registered function.
//' @field initialize\cr \itemize{
//'   \item Parameter: name `character(1)`\cr Name of the function.
//'   \item Parameter: is_discr `logical(1)`\cr Wheter discrete (not used).
//'   \item Parameter: args `numeric()`\cr Custom arguments.
//'   \item Parameter: range `numeric()`\cr Minimum and maximum range of the function.
//' }
//' @field call\cr \itemize{
//'   \item Parameter: x `numeric()`\cr Input vector.
//' }
//' @return `numeric()`
RCPP_MODULE(housepi_function)
{
#define FUNCTION_REGISTRIES {&UserFuncRegistry::inst(), &TranscendentalRegistry::inst(), &ArithmeticRegistry::inst()}
    MAKE_HOUSEPI_RCPP_WRAPPER_CLASS(Function, "function", FUNCTION_REGISTRIES);

    Rcpp::class_<FunctionHousepi>("FunctionHousepi")
        .constructor<std::string, bool, SEXP, SEXP, bool>()
        .method("call", &FunctionHousepi::call);
}

//' @name DistributionHousepi
//' @title Distribution from Housepi
//' @description Call housepi registered distribution.
//' @field initialize\cr \itemize{
//'   \item Parameter: name `character(1)`\cr Name of the distribution.
//'   \item Parameter: is_discr `logical(1)`\cr Wheter discrete.
//'   \item Parameter: args `numeric()`\cr Custom arguments.
//'   \item Parameter: range `numeric()`\cr Minimum and maximum range of the distribution.
//' }
//' @field call\cr \itemize{
//'   \item Parameter: x `numeric()`\cr Input vector.
//' }
//' @return `numeric()`
RCPP_MODULE(housepi_distribution)
{
#define DISTRIBUTION_REGISTRIES {&ContinuosDistrRegistry::inst(), &DiscreteDistrRegistry::inst()}
    MAKE_HOUSEPI_RCPP_WRAPPER_CLASS(Distribution, "distribution", DISTRIBUTION_REGISTRIES);

    Rcpp::class_<DistributionHousepi>("DistributionHousepi")
        .constructor<std::string, bool, SEXP, SEXP, bool>()
        .method("call", &DistributionHousepi::call);
}

//' @name InfecProfileHousepi
//' @title Infectivity Profile from Housepi
//' @description Call housepi registered infectivity profile.
//' @field initialize\cr \itemize{
//'   \item Parameter: name `character(1)`\cr Name of the profile.
//'   \item Parameter: is_discr `logical(1)`\cr Wheter discrete.
//'   \item Parameter: args `numeric()`\cr Custom arguments.
//'   \item Parameter: range `numeric()`\cr Minimum and maximum range of the profile.
//' }
//' @field call\cr \itemize{
//'   \item Parameter: x `numeric()`\cr Input vector.
//' }
//' @return `numeric()`
#define INFEC_PROFILE_REGISTRIES {&InfecProfileRegistry::inst(), &ContinuosDistrRegistry::inst(), &DiscreteDistrRegistry::inst()}
RCPP_MODULE(housepi_infec_profile)
{
    MAKE_HOUSEPI_RCPP_WRAPPER_CLASS(InfecProfile, "infectivity profile", INFEC_PROFILE_REGISTRIES);

    Rcpp::class_<InfecProfileHousepi>("InfecProfileHousepi")
        .constructor<std::string, bool, SEXP, SEXP, bool>()
        .method("call", &InfecProfileHousepi::call);
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
//' @field evaluate\cr Evaluate the loglikelihood according to the current state.
//' @field logevidence\cr Evaluate the logarithm of the marginal likelihood (aka evidence).
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
//' @field export_house_loglike_draws\cr Return a NumericMatrix with household loglikelihood draws.
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
        .method("evaluate", &MCMC::evaluate)
        .method("logevidence", &MCMC::logevidence)
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
        .method("export_augmented_data_wide", &MCMC::export_augmented_data_wide)
        .method("export_infec_profile", &MCMC::export_infec_profile)
        .method("export_draws", &MCMC::export_draws)
        .method("export_house_loglike_draws", &MCMC::export_house_loglike_draws)
        .method("export_house_snapshot", &MCMC::export_house_snapshot);

    Rcpp::class_<Model>("ModelHousepi")
        .constructor()
        .derives<MCMC>("McmcHousepi");
}
