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

// macro to define RCPP wrapper classes to export to R
#define DEFINE_HOUSEPI_RCPP_WRAPPER_CLASS(WrapperName, type_string, registry_list)               \
    class WrapperName##Housepi                                                                   \
    {                                                                                            \
    private:                                                                                     \
        std::shared_ptr<Function> _ptr;                                                          \
                                                                                                 \
    public:                                                                                      \
        WrapperName##Housepi(std::string name, bool is_discr, Rcpp::List args, Rcpp::List range) \
        {                                                                                        \
            try                                                                                  \
            {                                                                                    \
                _ptr = compose_func_registries(name, type_string, registry_list);                \
                                                                                                 \
                const Rcpp::CharacterVector arg_names = args.names();                            \
                for (R_xlen_t i = 0; i < args.size(); ++i)                                       \
                {                                                                                \
                    bool known_arg = true;                                                       \
                    const std::string this_arg = Rcpp::as<std::string>(arg_names[i]);            \
                    try                                                                          \
                    {                                                                            \
                        _ptr->assert_arg(this_arg);                                              \
                    }                                                                            \
                    catch (const std::exception &e)                                              \
                    {                                                                            \
                        tools::warning(e.what());                                                \
                        known_arg = false;                                                       \
                    }                                                                            \
                    if (known_arg)                                                               \
                        _ptr->init_attr(this_arg, HSP_ATTR_ARGUMENT, Rcpp::as<double>(args[i])); \
                }                                                                                \
                                                                                                 \
                const Rcpp::CharacterVector range_names = range.names();                         \
                for (R_xlen_t i = 0; i < range.size(); ++i)                                      \
                    if (!tools::is_na(Rcpp::as<double>(range[i])))                               \
                        _ptr->init_attr(Rcpp::as<std::string>(range_names[i]),                   \
                                        HSP_ATTR_RANGE,                                          \
                                        Rcpp::as<double>(range[i]));                             \
                                                                                                 \
                _ptr->init_reshaping(is_discr);                                                  \
            }                                                                                    \
            catch (const std::exception &e)                                                      \
            {                                                                                    \
                tools::error(                                                                    \
                    "While initializing " type_string " \"" + name + "\":\n" + e.what());        \
            }                                                                                    \
        }                                                                                        \
        ~WrapperName##Housepi() = default;                                                       \
                                                                                                 \
        Rcpp::NumericVector call(Rcpp::NumericVector x)                                          \
        {                                                                                        \
            try                                                                                  \
            {                                                                                    \
                const R_xlen_t n = x.size();                                                     \
                Rcpp::NumericVector y(n);                                                        \
                for (R_xlen_t i = 0; i < n; ++i)                                                 \
                    y[i] = _ptr->return_value(x[i]);                                             \
                return y;                                                                        \
            }                                                                                    \
            catch (const std::exception &e)                                                      \
            {                                                                                    \
                tools::error(                                                                    \
                    "While calling " type_string " \"" + _ptr->name() + "\":\n" + e.what());     \
            }                                                                                    \
            return Rcpp::NumericVector();                                                        \
        }                                                                                        \
    };

//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//

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
    DEFINE_HOUSEPI_RCPP_WRAPPER_CLASS(Function, "function", FUNCTION_REGISTRIES);

    Rcpp::class_<FunctionHousepi>("FunctionHousepi")
        .constructor<std::string, bool, Rcpp::List, Rcpp::List>()
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
    DEFINE_HOUSEPI_RCPP_WRAPPER_CLASS(Distribution, "distribution", DISTRIBUTION_REGISTRIES);

    Rcpp::class_<DistributionHousepi>("DistributionHousepi")
        .constructor<std::string, bool, Rcpp::List, Rcpp::List>()
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
    DEFINE_HOUSEPI_RCPP_WRAPPER_CLASS(InfecProfile, "infectivity profile", INFEC_PROFILE_REGISTRIES);

    Rcpp::class_<InfecProfileHousepi>("InfecProfileHousepi")
        .constructor<std::string, bool, Rcpp::List, Rcpp::List>()
        .method("call", &InfecProfileHousepi::call);
}

//' @name ModelHousepi
//' @title Internal housepi model
//' @description Perform inference and simulation from the model.
//' @field infer\cr Perform MCMC inference of the model parameters. \itemize{
//'   \item Parameter: output_file - The file path where the output will be saved
//'   \item Parameter: progbar - Display the progress bar
//' }
//' @field dataset_for_simulation\cr Dataset without status and augmented entries.
//' @field simulate\cr Simulated status and augmented variables.
//' @field evaluate\cr Evaluate the loglikelihood according to the current state.
//' @field evidence\cr Evaluate the logarithm of the marginal likelihood with the method of power posteriors. \itemize{
//'   \item Parameter: scheme
//'   \item Parameter: n_rungs
//' }
//' @field import_setup\cr Import setup from R List.
//' @field import_dataset\cr Import house indexed dataset from R Data.Frame.
//' @field import_datalist\cr Import individual indexed datalist from R Data.Frame.
//' @field import_multiplicity\cr Imoprt multiplicity definitions from R Data.Frame.
//' @field import_group_rates\cr Imoprt group rates definitions from R Data.Frame.
//' @field import_param\cr Import parameter definitions from R Data.Frame.
//' @field import_infec_profile\cr Import infectivity profile definition from R List.
//' @field import_augm\cr Import augmented variable definitions from R Data.Frame.
//' @field buil_all_rates\cr Assign all contact and group rates.
//' @field export_setup\cr Return a List with current setup.
//' @field export_param\cr Return a DataFrame with current parameters.
//' @field export_infec_profile\cr Return a List with infectivity profile.
//' @field export_dataset\cr Return a DataFrame with current database.
//' @field export_dataset_report\cr Return a DataFrame with a report of the current database.
//' @field export_house_size_report\cr Return a DataFrame with a report of household sizes.
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
        .method("dataset_for_simulation", &MCMC::dataset_for_simulation)
        .method("simulate", &MCMC::simulate)
        .method("evaluate", &MCMC::evaluate)
        .method("evidence", &MCMC::evidence)
        .method("import_dataset", &MCMC::import_dataset)
        .method("import_datalist", &MCMC::import_datalist)
        .method("import_multiplicity", &MCMC::import_multiplicity)
        .method("import_group_rates", &MCMC::import_group_rates)
        .method("import_setup", &MCMC::import_setup)
        .method("import_param", &MCMC::import_param)
        .method("import_infec_profile", &MCMC::import_infec_profile)
        .method("import_augm", &MCMC::import_augm)
        .method("build_all_rates", &MCMC::build_all_rates)
        .method("export_dataset", &MCMC::export_dataset)
        .method("export_datalist", &MCMC::export_datalist)
        .method("export_dataset_report", &MCMC::export_dataset_report)
        .method("export_house_size_report", &MCMC::export_house_size_report)
        .method("export_setup", &MCMC::export_setup)
        .method("export_param", &MCMC::export_param)
        .method("export_augmented_dataset", &MCMC::export_augmented_dataset)
        .method("export_infec_profile", &MCMC::export_infec_profile)
        .method("export_draws", &MCMC::export_draws)
        .method("export_house_loglike_draws", &MCMC::export_house_loglike_draws)
        .method("export_house_snapshot", &MCMC::export_house_snapshot);

    Rcpp::class_<Model>("ModelHousepi")
        .constructor()
        .derives<MCMC>("McmcHousepi");
}
