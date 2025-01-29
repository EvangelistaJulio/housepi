/*
 * mcmc_build.cpp
 *
 * Copyright (C) 2025 Francesco Camaglia, Institut Pasteur
 */

#include <fstream>
#include <iomanip>

#include "mcmc_build.h"

//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//

/* Auxiliary function to call percentuage. */
inline const std::string _aux_percent_(int num, int den)
{
    return "(" + tools::percentage(static_cast<double>(num) / den, 1) + ")";
}

/////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////  MCMC::Impl : methods (type-defined variables)  ////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

/* Convert a typed std::string to a double */
double MCMC::Impl::type_to_double(const std::string &value,
                                  const std::string &type,
                                  const Household &house, indiv i)
{
    double output = NULL_DOUBLE;

    try
    {
        if (!(tools::is_na(value) && tools::is_na(type)))
        {
            if (!tools::is_na(value) && tools::is_na(type))
                output = tools::stod(value, "argument without type");

            else if (tools::is_na(value) && !tools::is_na(type))
                tools::warning(
                    "Ignoring value type without a paired value.");

            else if (type == TYPE_SCALAR)
                output = tools::stod(value, "argument with scalar type");

            else if (type == TYPE_DATA)
                output = house.data(i, value);

            else if (type == TYPE_PARAM)
                // UPDATE: here it won't be straightforward to drop maps
                if (param_map.at(value).is_constant())
                    output = get_param(value);
                else
                    throw std::invalid_argument(
                        "Parameter to be inferred \"" + value + "\" assigned as constant argument.");

            else
                throw std::invalid_argument(
                    "Unknown variable type specifier \"" + type + "\".");
        }
    }
    catch (const std::exception &e)
    {
        tools::error(e.what());
    }

    return output;
}

/* Convert a std::string with type-dictionary to std::string-to-double */
void MCMC::Impl::assign_augm_type_args(Household &house, indiv i, augmID a,
                                       const str_vec &keys,
                                       const str_vec &vals,
                                       const str_vec &typs)
{
    // check size validity of keys, values and type
    if ((keys.size() != vals.size()) || (typs.size() != vals.size()))
        throw std::invalid_argument(
            "Mismatching definitions in arguments.");

    // loop over all arguments
    for (size_t x = 0; x < keys.size(); ++x)
    {
        std::string key = keys[x];
        std::string value = vals[x];
        std::string type = typs[x];

        if (tools::is_na(value) && tools::is_na(type))
            ; // pass

        else if (!tools::is_na(value) && tools::is_na(type))
            throw std::invalid_argument(
                "Impossible to interpret default value of  \"" + value + "\". Please provide a type.");

        else if (tools::is_na(value) && !tools::is_na(type))
            tools::warning(
                "Found a type without a paired value. Falling back to default.");

        else if (type == TYPE_DATA)
            house.assign_augm_model_scalar_arg(i, a, key, house.data(i, value));

        else if (type == TYPE_PARAM)
        {
            if (param_map.at(value).is_constant())
                house.assign_augm_model_scalar_arg(i, a, key, get_param(value));
            else
                house.assign_augm_model_pointer_arg(i, a, key, get_param_pointer(value));
        }
        else if (type == TYPE_SCALAR)
            house.assign_augm_model_scalar_arg(i, a, key, tools::stod(value, "argument with scalar type"));

        else
            throw std::invalid_argument(
                "Unsopported type \"" + type + "\" found while loading \"" + value + "\" for pointer value.");
    }
}

/* Convert a std::string with type-dictionary to std::string-to-double */
void MCMC::Impl::assign_augm_type_range(Household &house, indiv i, augmID a,
                                        const std::string &min_value,
                                        const std::string &min_type,
                                        const std::string &max_value,
                                        const std::string &max_type)
{
    //-------------------//
    //----- Minimum -----//
    //-------------------//

    if (tools::is_na(min_value) && tools::is_na(min_type))
        ; // pass

    else if (!tools::is_na(min_value) && tools::is_na(min_type))
        throw std::invalid_argument(
            "Impossible to interpret default value of  \"" + min_value + "\". Please provide a type.");

    else if (tools::is_na(min_value) && !tools::is_na(min_type))
        tools::warning(
            "Found a type without a paired value. Falling back to default.");

    else if (min_type == TYPE_DATA)
        house.assign_augm_model_scalar_range_min(i, a, house.data(i, min_value));

    else if (min_value == TYPE_PARAM)
    {
        if (param_map.at(min_value).is_constant())
            house.assign_augm_model_scalar_range_min(i, a, get_param(min_value));
        else
            house.assign_augm_model_pointer_range_min(i, a, get_param_pointer(min_value));
    }

    else if (min_type == TYPE_SCALAR)
        house.assign_augm_model_scalar_range_min(i, a,
                                                 tools::stod(min_value, "scalar type argument"));

    else
        throw std::invalid_argument(
            "Unsopported type \"" + min_type + "\" found while loading \"" + min_value + "\" for pointer value.");

    //-------------------//
    //----- Maximum -----//
    //-------------------//

    if (tools::is_na(max_value) && tools::is_na(max_type))
        ; // pass

    else if (!tools::is_na(max_value) && tools::is_na(max_type))
        throw std::invalid_argument(
            "Impossible to interpret default value of  \"" + max_value + "\". Please provide a type.");

    else if (tools::is_na(max_value) && !tools::is_na(max_type))
        tools::warning(
            "Found a type without a paired value. Falling back to default.");

    else if (max_type == TYPE_DATA)
        house.assign_augm_model_scalar_range_max(i, a, house.data(i, max_value));

    else if (max_type == TYPE_PARAM)
    {
        if (param_map.at(max_value).is_constant())
            house.assign_augm_model_scalar_range_max(i, a, get_param(max_value));
        else
            house.assign_augm_model_pointer_range_max(i, a, get_param_pointer(max_value));
    }

    else if (max_type == TYPE_SCALAR)
        house.assign_augm_model_scalar_range_max(i, a,
                                                 tools::stod(max_value, "argument with scalar type"));

    else
        throw std::invalid_argument(
            "Unsopported type \"" + max_type + "\" found while loading \"" + max_value + "\" for pointer value.");

    //---------------------//
    //----- Reshaping -----//
    //---------------------//

    house.init_augm_model_reshaping(i, a);
}

//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//

/* Use the loaded data to fill-in the respective report. */
void MCMC::Impl::update_report_data()
{
    _report_data_ = {{"n.cases", 0},
                     {"n.miss.stat", 0},
                     {"no.HH.case", 0},
                     {"lst.1.case", 0},
                     {"lst.2.case", 0},
                     {"max.n.HH.cases", 0},
                     {"n.HH.with.max.cases", 0},
                     {"min.HH.size", 1000},
                     {"n.HH.with.min.size", 0},
                     {"max.HH.size", 0},
                     {"n.HH.with.max.size", 0}};

    //===== Set Contact Rates and Perform Database Analysis =====//
    for (const Household &house : house_vector)
    {
        ushort house_size = house.size();
        ushort n_infec = house.n_infectious();

        _report_data_["n.cases"] += n_infec;
        _report_data_["n.miss.stat"] += house.n_missing_status();

        _report_data_["no.HH.case"] += (n_infec == 0);
        _report_data_["lst.1.case"] += (n_infec > 0);
        _report_data_["lst.2.case"] += (n_infec > 1);

        // max number of house cases
        if (_report_data_["max.n.HH.cases"] < n_infec)
        {
            _report_data_["n.HH.with.max.cases"] = 1;
            _report_data_["max.n.HH.cases"] = n_infec;
        }
        else if (_report_data_["max.n.HH.cases"] == n_infec)
        {
            _report_data_["n.HH.with.max.cases"] += 1;
        }

        // minimum house number of members
        if (_report_data_["min.HH.size"] > house_size)
        {
            _report_data_["n.HH.with.min.size"] = 1;
            _report_data_["min.HH.size"] = house_size;
        }
        else if (_report_data_["min.HH.size"] == house_size)
        {
            _report_data_["n.HH.with.min.size"] += 1;
        }

        // maximum house number of members
        if (_report_data_["max.HH.size"] < house_size)
        {
            _report_data_["n.HH.with.max.size"] = 1;
            _report_data_["max.HH.size"] = house_size;
        }
        else if (_report_data_["max.HH.size"] == house_size)
        {
            _report_data_["n.HH.with.max.size"] += 1;
        }
    }
}

/* Return a value from report data. */
size_t MCMC::Impl::report_data(const std::string &key) const
{
    size_t output = 0;
    if (_report_data_.find(key) != _report_data_.end())
        output = _report_data_.at(key);

    else
        throw std::invalid_argument(
            "Unknown summary entry \"" + key + "\".");

    return output;
}

//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////  MCMC::Impl : build augmented  //////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/* Read Model Augmented Variable file. */
void MCMC::Impl::read_augm(const std::string &pathtofile)
{
    try
    {
        tools::read_csv_table(pathtofile, _augm_deftable_, augm_file_headers);
    }
    catch (const std::exception &e)
    {
        tools::error(
            "While reading augmented definitions: " + std::string(e.what()));
    }

    // use loaded definitions to build
    build_augm();
    // store the file name in the report
    _report_pathtofile_["augm"] = pathtofile;
}

/* Assign augmented variables to all individuals.
    UPDATE: this function has to be simplified.
*/
void MCMC::Impl::build_augm()
{
    try
    {
        const dict_str_vec &deftable = _augm_deftable_;
        if (!deftable.empty())
        {
            // check no unused options in deftable
            str_vec entries = {};
            for (const auto &pair : deftable)
                entries.push_back(pair.first);
            tools::check_unused_entries(entries, augm_file_headers);

            // n. of definitions
            size_t table_length = deftable.at("name").size();

            if (table_length > 0)
            {
                // auxiliary variable to load definitions : a -> status lists
                dict_str_vec name_status_check = {};

                //===== Check Augmented Variables lists of usage =====//
                for (size_t k = 0; k < table_length; ++k)
                {
                    std::string a = deftable.at("name")[k];
                    std::string status = deftable.at("status")[k];
                    bool is_init_rand = false;

                    // check that augmID and status are correct
                    if (tools::is_na(a))
                        throw std::invalid_argument(
                            "Augmented variable name can not be undefined.");

                    else if (name_status_check.find(a) == name_status_check.end())
                        // new augmented variable status
                        name_status_check.emplace(a, str_vec{status});

                    else
                    {
                        // check no multiple augmID with same status
                        const str_vec &mylist = name_status_check.at(a);
                        if (std::find(mylist.begin(), mylist.end(), status) == mylist.end())
                            name_status_check.at(a).push_back(status);

                        else
                            throw std::invalid_argument(
                                "Incompatible definitions for augmented variable \"" + a + "\".");
                    }

                    // option for random initialization
                    if (!tools::is_na(deftable.at("init_rand")[k]))
                        is_init_rand = tools::stobool(deftable.at("init_rand")[k]);

                    if ((tools::is_na(deftable.at("model")[k])) && (is_init_rand))
                        throw std::invalid_argument(
                            "Impossible to initialize \"" + a + "\" without a model.");
                }

                //===== Drop previous definitions =====//
                augm_list_full.clear();
                augm_list_to_be_inferred.clear();
                augm_list_model.clear();
                augm_status_map_infer.clear();
                augm_status_map_model.clear();

                //===== Check Univocal Definitions and Allocate New Lists =====//
                for (const auto &pair : name_status_check)
                {
                    const std::string &a = pair.first;
                    const str_vec &mylist = pair.second;

                    for (const std::string &status : mylist)
                    {
                        if (tools::is_na(status))
                        {
                            if (mylist.size() > 1)
                                throw std::invalid_argument(
                                    "Conflicting status-specific and general definitions for augmented \"" + a + "\".");
                        }
                        else
                        {
                            int s = tools::stoi(status, "augmented status");
                            if (s < INFECTED_STATUS)
                            {
                                throw std::invalid_argument(
                                    "Impossible to associate augmented \"" + a + "\" for non infectious status.");
                            }
                            else
                            {
                                augm_status_map_infer.emplace(s, std::unordered_set<std::string>{});
                                augm_status_map_model.emplace(s, std::unordered_set<std::string>{});
                            }
                        }

                        augm_list_full.emplace(a);
                    }
                }

                //===== Case of Only Generic Definitions =====//
                if ((table_length > 0) && (augm_status_map_infer.empty()))
                    _func_augm_model_by_definition_ = &MCMC::Impl::_augm_model_generic_;
                else
                    _func_augm_model_by_definition_ = &MCMC::Impl::_augm_model_by_status_;

                //===========================================//
                //===== Assign new definitions to lists =====//
                //===========================================//

                // note : I didn't manage to collapse this with the previous because of generic status
                for (size_t k = 0; k < table_length; ++k)
                {
                    const std::string &a = deftable.at("name")[k];
                    const std::string &status = deftable.at("status")[k];

                    //----- Augmented Variables to be Inferred -----//
                    if (!tools::is_na(deftable.at("rule")[k]))
                    {
                        augm_list_to_be_inferred.emplace(a);

                        if (tools::is_na(status)) // generic definition
                        {
                            if (!is_augm_generic())
                                for (const auto &pair : augm_status_map_infer)
                                    augm_status_map_infer.at(pair.first).emplace(a);
                        }
                        else // status-specific definition
                        {
                            augm_status_map_infer.at(tools::stoi(status, "augmented status")).emplace(a);
                        }
                    }

                    //----- Augmented Variables with a Model -----//
                    if (!tools::is_na(deftable.at("model")[k]))
                    {
                        augm_list_model.emplace(a);

                        if (tools::is_na(status)) // generic definition
                        {
                            if (!is_augm_generic())
                                for (const auto &pair : augm_status_map_model)
                                    augm_status_map_model.at(pair.first).emplace(a);
                        }
                        else // status-specific definition
                        {
                            augm_status_map_model.at(tools::stoi(status)).emplace(a);
                        }
                    }
                }

                //================================================//
                //===== Finalize assigment in the households =====//
                //================================================//

                if ((!is_augm_infer()) && (!tools::is_na(setup_id_trace)))
                {
                    tools::warning(
                        "Ignoring id trace: no augmented variables to be inferred was provided.");
                    setup_id_trace = NULL_STRING;
                }

                // Allocate memory for augmented variables to each individual
                // Search individual ID for augmented traces follow-up
                dgn_H_star = n_households;
                for (size_t H = 0; H < n_households; ++H)
                {
                    Household &house = house_vector[H];
                    for (indiv i = 0; i < house.size(); ++i)
                    {
                        clear_all_augm(house, i);
                        if (!tools::is_na(setup_id_trace))
                        {
                            if (house.ID(i) == setup_id_trace)
                            {
                                dgn_H_star = H;
                                dgn_i_star = i;
                            }
                        }
                    }
                }

                if (!tools::is_na(setup_id_trace))
                {
                    if (dgn_H_star == n_households)
                    {
                        tools::warning(
                            "Ignoring id trace: no individual \"" + setup_id_trace + "\" was found.");
                        setup_id_trace = NULL_STRING;
                    }
                    else if (house_vector[dgn_H_star].status(dgn_i_star) < INFECTED_STATUS)
                    {
                        tools::warning(
                            "Individual id \"" + setup_id_trace + "\" correspond to a healthy individual.");
                    }
                }
            }
        }
    }
    catch (const std::exception &e)
    {
        tools::error(
            "While building augmented variable definitions: " + std::string(e.what()));
    }
}

/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////  MCMC::Impl : assign augmented variables  //////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

/* Assign augmented variables to single individual. */
void MCMC::Impl::assign_all_augm(Household &house, indiv i)
{
    try
    {
        const dict_str_vec &deftable = _augm_deftable_;

        for (size_t k = 0; k < deftable.at("name").size(); ++k)
        {
            // assign augmented variable of the corresponding status
            const std::string &augm_status = deftable.at("status")[k];

            // assign augmented variable of the corresponding status or undefined status
            if ((tools::is_na(augm_status)) || (house.status(i) == tools::stoi(augm_status, "augmented status")))
            {
                //----- Index -----//
                // |====> UPDATE: to be implemented <===|
                augmID a = deftable.at("name")[k];

                //----- Name -----//
                const std::string &name = deftable.at("name")[k];

                //----- is Discrete -----//
                bool is_discrete = tools::stobool(deftable.at("discrete")[k]);

                //----- Value -----//
                double value = type_to_double(deftable.at("value")[k],
                                              deftable.at("value_type")[k],
                                              house, i);

                //----- Model -----//
                const std::string &model_name = deftable.at("model")[k];

                //----- Proposer Rule -----//
                const std::string &rule_name = deftable.at("rule")[k];
                double rule_scale = type_to_double(deftable.at("rule_scale")[k],
                                                   deftable.at("rule_scale_type")[k],
                                                   house, i);
                double rule_from = type_to_double(deftable.at("rule_from")[k],
                                                  deftable.at("rule_from_type")[k],
                                                  house, i);

                //----- Assign Augmented Variable -----//
                house.assign_augm(i, a, name, value,
                                  is_discrete, model_name,
                                  rule_name, rule_scale, rule_from);

                //----- Arguments -----//
                assign_augm_type_args(house, i, a,
                                      tools::parse_list_args(deftable.at("arg_name")[k]),
                                      tools::parse_list_args(deftable.at("arg_value")[k]),
                                      tools::parse_list_args(deftable.at("arg_type")[k]));

                //----- Range -----//
                assign_augm_type_range(house, i, a,
                                       deftable.at("min")[k],
                                       deftable.at("min_type")[k],
                                       deftable.at("max")[k],
                                       deftable.at("max_type")[k]);

                // random initialization at assignment
                if (tools::stobool(deftable.at("init_rand")[k]))
                    house.augm_model_rand(i, name, get_prng());
            }
        }
    }
    catch (const std::exception &e)
    {
        tools::error(
            "While initializing augmented variables: " + std::string(e.what()));
    }
}

/* Clear all augmented definition of individual. */
void MCMC::Impl::clear_all_augm(Household &house, indiv i)
{
    for (const augmID &a : augm_list_full)
        house.clear_augm(i, a);
}

/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////  MCMC::Impl : state checkers   /////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

/* Check if there are non-constant augmented variable definitions. */
bool MCMC::Impl::is_augm_infer() const { return !augm_list_to_be_inferred.empty(); }

/* Check if there are augmented variable definitions with a model. */
bool MCMC::Impl::is_augm_model() const { return !augm_list_model.empty(); }

//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//

///////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////  MCMC : methods (diagnostics) ///////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/* Init Bayesian diagnostics. */
void MCMC::Impl::init_diagnostics()
{
    ulong sample_size = setup_c_sample_size * setup_n_chains;

    // diagnostics over global loglikelihood
    dgn_count = 0;
    dgn_glob_loglike = dbl_vec(sample_size);

    // initialize  other diagnostic vectors
    if (!is_power_post())
    {
        // diagnostics over parameters
        for (const parID &p : param_to_be_inferred)
            dgn_param_value[p] = dbl_vec(sample_size, std::numeric_limits<double>::quiet_NaN());

        // diagnostics over households
        for (size_t H = 0; H < n_households; ++H)
            house_vector[H].init_diagnostics(sample_size);

        // diagnostics over augmented variables
        if (!tools::is_na(setup_id_trace))
            for (const augmID &a : augm_list_to_be_inferred)
                dgn_augm_value[a] = dbl_vec(sample_size, std::numeric_limits<double>::quiet_NaN());
    }
}

/* Update Bayesian diagnostics */
void MCMC::Impl::update_diagnostics()
{
    // diagnostics for global loglike
    dgn_glob_loglike[dgn_count] = get_global_loglike();

    if (!is_power_post())
    {
        // diagnostics for all parameters
        for (const parID &p : param_to_be_inferred)
            dgn_param_value.at(p)[dgn_count] = get_param(p);

        // diagnostics for all households
        for (size_t H = 0; H < n_households; ++H)
            house_vector[H].update_diagnostics(dgn_count);

        // diagnostics for individual followup
        if (!tools::is_na(setup_id_trace))
            for (const augmID &a : augm_list_to_be_inferred)
                dgn_augm_value.at(a)[dgn_count] = house_vector[dgn_H_star].augm(dgn_i_star, a);
    }

    ++dgn_count;
}

/* Check if augmented variable definition are generic or not (i. e. individual status specific). */
bool MCMC::Impl::is_augm_generic() const
{
    return _func_augm_model_by_definition_ == &MCMC::Impl::_augm_model_generic_;
}

/* Compute the loglikelihood term associated to augmented variables. */
double MCMC::Impl::compute_augm_model_loglike(const Household &house) const
{
    double hh_loglike = 0.;
    // Contribution of model associated to augmented variables
    for (const indiv &i : house.infectious_set())
        (this->*_func_augm_model_by_definition_)(house, i, hh_loglike);

    return hh_loglike;
}

void MCMC::Impl::_augm_model_by_status_(const Household &house, indiv i, double &delta) const
{
    for (const augmID &a : augm_status_map_model.at(house.status(i)))
        delta += house.augm_model_logpdf(i, a);
}

void MCMC::Impl::_augm_model_generic_(const Household &house, indiv i, double &delta) const
{
    for (const augmID &a : augm_list_model)
        delta += house.augm_model_logpdf(i, a);
}

//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//

/* Tell if it's the case of power posteriors. */
bool MCMC::Impl::is_power_post() const { return _n_rungs_ > 0; }

/* The text to print in the Markov chian progress bar. */
const std::string MCMC::Impl::fmt_progtxt(const std::string &txt)
{
    if (is_power_post())
        return "power post " + std::to_string(_rung_) + "/" + std::to_string(_n_rungs_);
    else if (setup_n_chains > 1)
        return txt + " " + std::to_string(_chain_) + "/" + std::to_string(setup_n_chains);
    else
        return txt;
}

//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//

///////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////  MCMC : show summary  ///////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/* Print the report of the database and of the model group rates. */
void MCMC::Impl::show_summary() const
{
    std::ostringstream oss;

    oss << "n. of individuals: " << n_individuals
        << "\n    * n. of cases: " << _report_data_.at("n.cases")
        << " " << _aux_percent_(_report_data_.at("n.cases"), n_individuals);

    if (_report_data_.at("n.miss.stat") > 0)
        oss << "\n    * n. of missing status : " << _report_data_.at("n.miss.stat")
            << " " << _aux_percent_(_report_data_.at("n.miss.stat"), n_individuals);

    oss << "\nn. of households: " << n_households;

    if (_report_data_.at("no.HH.case") > 0)
        oss << "\n    * with no cases: " << _report_data_.at("no.HH.case")
            << " " << _aux_percent_(_report_data_.at("no.HH.case"), n_households);

    oss << "\n    * with at least 1 case: " << _report_data_.at("lst.1.case")
        << " " << _aux_percent_(_report_data_.at("lst.1.case"), n_households);

    if (_report_data_.at("lst.1.case") > 0)
        oss << "\n    * with at least 2 cases: " << _report_data_.at("lst.2.case")
            << " " << _aux_percent_(_report_data_.at("lst.2.case"), n_households);

    oss << "\nmax n. house cases: " << _report_data_.at("max.n.HH.cases")
        << " " << _aux_percent_(_report_data_.at("n.HH.with.max.cases"), n_households)
        << "\nmin n. house members: " << _report_data_.at("min.HH.size")
        << " " << _aux_percent_(_report_data_.at("n.HH.with.min.size"), n_households)
        << "\nmax n. house members: " << _report_data_.at("max.HH.size")
        << " " << _aux_percent_(_report_data_.at("n.HH.with.max.size"), n_households);

    tools::print_stdout(oss.str() + "\n");
}

///////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////  MCMC : show full build  /////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

void MCMC::Impl::show_full_build() const
{
    //----- SETUP -----//
    std::string add_file_info = NULL_STRING;
    if (_report_pathtofile_.find("setup") != _report_pathtofile_.end())
        if (!tools::is_na(_report_pathtofile_.at("setup")))
            add_file_info = "# import: " + _report_pathtofile_.at("setup");

    tools::print_decorator("Setup", add_file_info);
    show_setup();

    //----- PARAMETERS -----//
    add_file_info = NULL_STRING;
    if (_report_pathtofile_.find("params") != _report_pathtofile_.end())
        if (!tools::is_na(_report_pathtofile_.at("params")))
            add_file_info = "# import: " + _report_pathtofile_.at("params");

    tools::print_decorator("Parameters", add_file_info);
    show_all_params();

    //----- DATABASE -----//
    add_file_info = NULL_STRING;
    if (_report_pathtofile_.find("data") != _report_pathtofile_.end())
        if (!tools::is_na(_report_pathtofile_.at("data")))
            add_file_info = "# import: " + _report_pathtofile_.at("data");

    tools::print_decorator("Database", add_file_info);
    show_summary();

    //----- AUGMENTED -----//
    add_file_info = NULL_STRING;
    if (_report_pathtofile_.find("augm") != _report_pathtofile_.end())
        if (!tools::is_na(_report_pathtofile_.at("augm")))
            add_file_info = "# import: " + _report_pathtofile_.at("augm");

    if (!augm_list_full.empty())
        tools::print_decorator("Augmented", add_file_info);
}

//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//

///////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////  MCMC::Impl : print augm variables  ///////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/* Print all augmented variables of each individual */
void MCMC::Impl::print_augm_variables(const std::string &pathtofile) const
{
    if (!tools::is_na(pathtofile))
    {
        std::ofstream augm_ofstream;
        augm_ofstream.open(pathtofile.c_str());

        if (augm_ofstream.is_open())
        {
            //===== Print column headers =====//
            // default headers
            augm_ofstream << "id" << C_DELIM << "house_id" << C_DELIM << "status";

            // data headers
            for (const dataID &d : data_headers)
                augm_ofstream << C_DELIM << d;

            // augmented variables headers
            for (const augmID &a : augm_list_full)
                augm_ofstream << C_DELIM << a;

            augm_ofstream << std::endl;

            //===== Print data and augmented variables =====//
            for (const Household &house : house_vector)
            {
                for (indiv i = 0; i < house.size(); ++i)
                {
                    augm_ofstream << house.ID(i)
                                  << C_DELIM << house.house_ID()
                                  << C_DELIM << house.status(i);
                    for (const dataID &d : data_headers)
                        augm_ofstream << C_DELIM << house.data(i, d);

                    for (const augmID &a : augm_list_full)
                        if (house.status(i) == INFECTED_STATUS)
                            augm_ofstream << C_DELIM << house.augm(i, a);
                        else
                            augm_ofstream << C_DELIM << std::numeric_limits<double>::quiet_NaN();

                    augm_ofstream << std::endl;
                }
            }

            augm_ofstream.close();
        }
        else
        {
            throw std::runtime_error(
                "Unable to open file: " + pathtofile);
        }
    }
}
