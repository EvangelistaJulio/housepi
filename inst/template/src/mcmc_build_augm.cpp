/*      __
 *     / /_  ____  __  __ ____ __  ____  __
 *    / __ \/ __ \/ / / / ___/ _ \/ __ \/ /
 *   / / / / /_/ / /_/ (__  )  __/ /_/ / /
 *  /_/ /_/\____/\____/____/\___/ ____/_/
 *                             /_/
 *
 * Copyright (C) 2025 Francesco Camaglia, Institut Pasteur
 *
 * mcmc_build_augm.cpp
 *
 */

#include "mcmc_build.h"

//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//

/////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////  MCMC::Impl : methods (type-defined variables)  /////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

/* Convert a typed std::string to a double */
double MCMC::Impl::_type_to_double(const std::string &value,
                                   const std::string &type,
                                   const Household &house, indiv i)
{
    double output = HOUSEPI_NA_DOUBLE;

    if (!(tools::is_na(value) && tools::is_na(type)))
    {
        if (!tools::is_na(value) && tools::is_na(type))
            output = tools::stod(value, "argument without type");

        else if (tools::is_na(value) && !tools::is_na(type))
            tools::warning(
                "Ignoring value type without a paired value.");

        else if (tools::stoupper(type) == TYPE_SCALAR)
            output = tools::stod(value, "argument with scalar type");

        else if (tools::stoupper(type) == TYPE_DATA)
            output = house.data(i, value);

        else if (tools::stoupper(type) == TYPE_PARAM)
            // NOTE: here it won't be straightforward to drop maps...
            if (param_map.at(value).is_constant())
                output = get_param(value);
            else
                throw std::invalid_argument(
                    "Parameter to be inferred \"" + value + "\" assigned as constant argument.");

        else
            throw std::invalid_argument(
                "Unknown variable type specifier \"" + type + "\".");
    }

    return output;
}

/* Auxiliary method to assign attributes to augmented function or model. */
void MCMC::Impl::_assign_augm_attributes(Household &house, indiv i, augmID a,
                                         const str_vec &keys,
                                         const str_vec &vals,
                                         const str_vec &typs,
                                         const std::vector<bool> &is_arg,
                                         bool is_var_model)
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
        bool flag = is_arg[x];

        if (tools::is_na(value) && tools::is_na(type))
            ; // pass

        else if (!tools::is_na(value) && tools::is_na(type))
            throw std::invalid_argument(
                "Impossible to interpret default value of  \"" + value + "\". Please provide a type.");

        else if (tools::is_na(value) && !tools::is_na(type))
            tools::warning(
                "Found a type without a paired value. Falling back to default.");

        else if (tools::stoupper(type) == TYPE_DATA)
            house.assign_augm_attr_scalar(i, a, key, flag, house.data(i, value), is_var_model);

        else if (tools::stoupper(type) == TYPE_PARAM)
        {
            if (param_map.at(value).is_constant())
                house.assign_augm_attr_scalar(i, a, key, flag, get_param(value), is_var_model);
            else
            {
                house.assign_augm_attr_pointer(i, a, key, flag, get_param_pointer(value), is_var_model);
                // state parameter cache dependency if variable is a function
                if (is_var_model == VAR_FUNCTION)
                    make_param_cache_dependent(value);
            }
        }

        else if (tools::stoupper(type) == TYPE_SCALAR)
            house.assign_augm_attr_scalar(i, a, key, flag, tools::stod(value, "argument with scalar type"), is_var_model);

        else if (tools::stoupper(type) == TYPE_AUGM)
        {
            if (contains(augm_list_full, value))
            {
                if (a == value)
                    throw std::invalid_argument(
                        "Impossible to set argument equal to its value for \"" + a + "\".");

                else if (house.is_augm_function(i, value))
                    house.assign_augm_attr_function(i, a, key, flag, house.get_augm_function(i, value), is_var_model);

                else if (house.is_augm_constant(i, value))
                    house.assign_augm_attr_scalar(i, a, key, flag, house.augm(i, value), is_var_model);

                else
                    house.assign_augm_attr_pointer(i, a, key, flag, house.get_augm_pointer(i, value), is_var_model);
            }
            else
                throw std::invalid_argument(
                    "Undeclared augmented variable \"" + value + "\".");
        }

        else
            throw std::invalid_argument(
                "Unsopported type \"" + type + "\" found while loading \"" + value + "\" for pointer value.");
    }
}

/* Assign augmented function arguments. */
void MCMC::Impl::_assign_augm_func_args(Household &house, indiv i, augmID a,
                                        const str_vec &keys,
                                        const str_vec &vals,
                                        const str_vec &typs)
{
    _assign_augm_attributes(house, i, a, keys, vals, typs,
                            std::vector<bool>(keys.size(), ATTR_ARGUMENT), VAR_FUNCTION);
}
/* Assign augmented model arguments. */
void MCMC::Impl::_assign_augm_model_args(Household &house, indiv i, augmID a,
                                         const str_vec &keys,
                                         const str_vec &vals,
                                         const str_vec &typs)
{
    _assign_augm_attributes(house, i, a, keys, vals, typs,
                            std::vector<bool>(keys.size(), ATTR_ARGUMENT), VAR_MODEL);
}
/* Assign augmented model range. */
void MCMC::Impl::_assign_augm_model_range(Household &house, indiv i, augmID a,
                                          const str_vec &vals,
                                          const str_vec &typs)
{
    _assign_augm_attributes(house, i, a, {"min", "max"}, vals, typs,
                            std::vector<bool>(2, ATTR_RANGE), VAR_MODEL);
}

//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////  MCMC::Impl : build augmented  ////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/* Produce auxiliary variable to load definitions : a -> status lists */
const dict_str_vec MCMC::Impl::_build_augm_name_status_check()
{
    const dict_str_vec &deftable = _augm_deftable;
    const size_t table_length = deftable.at("name").size();

    dict_str_vec output = {};

    // checka augmented variables lists of usage
    for (size_t k = 0; k < table_length; ++k)
    {
        const std::string &a = deftable.at("name")[k];
        const std::string &status = deftable.at("status")[k];

        // check that augmID and status are correct
        if (tools::is_na(a))
            throw std::invalid_argument(
                "Augmented variable name can not be undefined.");

        else if (output.count(a))
        {
            // check no multiple augmID with same status
            const str_vec &mylist = output.at(a);
            if (!contains(mylist, status))
                output.at(a).push_back(status);
            else
                throw std::invalid_argument(
                    "Incompatible definitions for augmented variable \"" + a + "\".");
        }

        else // new augmented variable status
            output.emplace(a, str_vec{status});
    }

    return output;
}

/* Check univocal definitions and allocate new lists. */
void MCMC::Impl::_build_augm_allocate_univocal_defs(const dict_str_vec &check)
{
    // drop previous definitions
    augm_list_full.clear();
    augm_status_map_infer.clear();
    augm_status_map_model.clear();

    for (const auto &pair : check)
    {
        const std::string &a = pair.first;
        const str_vec &mylist = pair.second;

        for (const std::string &status : mylist)
        {
            if (!tools::is_na(status))
            {
                int s = tools::stoi(status, "augmented status");
                if (s >= INFECTED_STATUS)
                {
                    augm_status_map_infer.emplace(s, std::vector<augmID>{});
                    augm_status_map_model.emplace(s, std::vector<augmID>{});
                }
                else
                    throw std::invalid_argument(
                        "Impossible to associate augmented \"" + a + "\" for non infectious status.");
            }
            else if (mylist.size() > 1)
                throw std::invalid_argument(
                    "Conflicting status-specific and general definitions for augmented \"" + a + "\".");

            augm_list_full.push_back(a);
        }
    }

    // set generic or status-specifc augmented variables
    if (augm_status_map_infer.empty())
        // only generic definitions
        _do_augm_model_by_definition = &MCMC::Impl::__augm_model_generic;
    else
        // status specific definitions
        _do_augm_model_by_definition = &MCMC::Impl::__augm_model_by_status;
}

/* Allocate memory for augmented variables to each individual.
Search individual id for augmented traces follow-up. */
void MCMC::Impl::_build_augm_clear_houses_and_set_id_trace()
{
    // cleaer all augmented definitions
    for (Household &house : house_vector)
        for (indiv i = 0; i < house.size(); ++i)
            clear_all_augm(house, i);

    if (is_augm_sample())
    {
        if (!tools::is_na(setup_id_trace))
        {
            smpl_H_star = n_households;
            for (Household &house : house_vector)
                for (indiv i = 0; i < house.size(); ++i)
                    if (!tools::is_na(setup_id_trace))
                        if (house.id(i) == setup_id_trace)
                        {
                            smpl_H_star = house.house_index();
                            smpl_i_star = i;
                        }

            if (smpl_H_star == n_households)
            {
                tools::warning(
                    "Ignoring id trace: no individual \"" + setup_id_trace + "\" was found.");
                setup_id_trace = HOUSEPI_NA_STRING;
            }
            else if (house_vector[smpl_H_star].status(smpl_i_star) < INFECTED_STATUS)
                tools::warning(
                    "Individual id \"" + setup_id_trace + "\" correspond to a healthy individual.");
        }
    }
    else if (!tools::is_na(setup_id_trace))
    {
        tools::warning(
            "Ignoring id trace: no augmented variables to be inferred was provided.");
        setup_id_trace = HOUSEPI_NA_STRING;
    }
}

/* Assign augmented variables to all individuals. */
void MCMC::Impl::build_augm()
{
    try
    {
        const dict_str_vec &deftable = _augm_deftable;
        if (!deftable.empty())
        {
            // check entries
            str_vec headers = tools::headers_from_deftable(_augm_deftable);
            tools::check_for_missing_entries(headers, augm_file_headers);
            tools::check_for_unused_entries(headers, augm_file_headers);

            // auxiliary variable to load definitions : a -> status lists
            const dict_str_vec name_status_check = _build_augm_name_status_check();
            _build_augm_allocate_univocal_defs(name_status_check);

            //===========================================//
            //===== Assign new definitions to lists =====//
            //===========================================//

            // drop previous definitions
            augm_list_infer.clear();
            augm_list_func.clear();
            augm_list_model.clear();

            // note : I didn't manage to collapse this with the previous because of generic status
            for (size_t k = 0; k < deftable.at("name").size(); ++k)
            {
                const std::string &a = deftable.at("name")[k];
                const std::string &status = deftable.at("status")[k];

                // read option for function
                bool is_function = false;
                if (!tools::is_na(deftable.at("func")[k]))
                    is_function = !tools::is_na(deftable.at("func")[k]);

                // read option for constant variable
                bool is_constant = false;
                if (!tools::is_na(deftable.at("const")[k]))
                    is_constant = tools::stobool(deftable.at("const")[k]);

                // add to list of augmented functions
                if (is_function)
                    augm_list_func.push_back(a);

                // Augmented Variables to be Inferred
                if (!is_constant && !is_function)
                {
                    // check random initialization option
                    if (!tools::is_na(deftable.at("init_rand")[k]))
                        if (tools::stobool(deftable.at("init_rand")[k]))
                            if (tools::is_na(deftable.at("model")[k]))
                                throw std::invalid_argument(
                                    "Impossible to initialize \"" + a + "\" without a model.");

                    // read option for monte carlo proposer
                    if (!tools::is_na(deftable.at("rule")[k]))
                    {
                        if (!contains(augm_list_infer, a))
                            augm_list_infer.push_back(a);

                        // note: impossible to associate rule to missing status
                        if (!tools::is_na(status)) // status-specific definition
                        {
                            int s = tools::stoi(status, "augmented status");
                            if (!contains(augm_status_map_infer.at(s), a))
                                augm_status_map_infer.at(s).push_back(a);
                        }

                        else if (!is_augm_generic()) // generic definition
                            for (const auto &pair : augm_status_map_infer)
                                if (!contains(augm_status_map_infer.at(pair.first), a))
                                    augm_status_map_infer.at(pair.first).push_back(a);
                    }
                }

                // Augmented Variables with a Model
                if (!tools::is_na(deftable.at("model")[k]))
                {
                    if (!contains(augm_list_model, a))
                        augm_list_model.push_back(a);

                    if (!tools::is_na(status)) // status-specific definition
                    {
                        int s = tools::stoi(status, "augmented status");
                        if (!contains(augm_status_map_model.at(s), a))
                            augm_status_map_model.at(s).push_back(a);
                    }

                    else if (!is_augm_generic()) // generic definition
                        for (const auto &pair : augm_status_map_model)
                            if (!contains(augm_status_map_model.at(pair.first), a))
                                augm_status_map_model.at(pair.first).push_back(a);
                }
            }

            // finalize assigment in the households
            _build_augm_clear_houses_and_set_id_trace();
        }
    }
    catch (const std::exception &e)
    {
        tools::error({"While building augmented variable definitions:",
                      std::string(e.what())});
    }
}

/////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////  MCMC::Impl : state checkers   /////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

/* Check if there are ant augmented variable definitions. */
bool MCMC::Impl::is_augm_any() const { return !augm_list_full.empty(); }

/* Check if there are non-constant augmented variable definitions. */
bool MCMC::Impl::is_augm_sample() const { return !augm_list_infer.empty(); }

/* Check if there are augmented variable definitions with a model. */
bool MCMC::Impl::is_augm_model() const { return !augm_list_model.empty(); }

/* Check if augmented variable definition are generic or not (i. e. individual status specific). */
bool MCMC::Impl::is_augm_generic() const { return _do_augm_model_by_definition == &MCMC::Impl::__augm_model_generic; }

/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////  MCMC::Impl : assign augmented variables  /////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

/* Assign augmented variables to single individual. */
void MCMC::Impl::assign_all_augm(Household &house, indiv i)
{
    try
    {
        int status = house.status(i);
        const dict_str_vec &deftab = _augm_deftable;
        std::map<augmID, size_t> indiv_augm_entries;

        // loop over augmented variable definitions to assign placheholders
        for (size_t k = 0; k < deftab.at("name").size(); ++k)
        {
            // assign augmented variable of the corresponding status
            const std::string &augm_status = deftab.at("status")[k];

            // assign augmented variable if undefined status | corresponding status
            if (tools::is_na(augm_status) || (status == tools::stoi(augm_status, "augmented status")))
            {
                augmID a = deftab.at("name")[k]; // UPDATE : index to be implemented

                // assign augmented variable placeholder
                // FIXME: what if constant = T, discrete = T and value != integer?
                house.assign_augm(i, a,
                                  deftab.at("name")[k],
                                  _type_to_double(deftab.at("value")[k], deftab.at("value_at")[k], house, i),
                                  tools::stobool(deftab.at("discr")[k]),
                                  tools::stobool(deftab.at("const")[k]),
                                  deftab.at("func")[k],
                                  deftab.at("model")[k]);

                // add entry to individual augmented variable list
                indiv_augm_entries.emplace(a, k);
            }
        }

        // assign attributes to augmented variables
        for (const auto &pair : indiv_augm_entries)
        {
            const augmID &a = pair.first;
            const size_t &k = pair.second;

            // assign arguments to model
            _assign_augm_model_args(house, i, a,
                                    tools::parse_list_args(deftab.at("arg_key")[k]),
                                    tools::parse_list_args(deftab.at("arg_val")[k]),
                                    tools::parse_list_args(deftab.at("arg_at")[k]));

            // assign range to model
            _assign_augm_model_range(house, i, a,
                                     {deftab.at("min")[k], deftab.at("max")[k]},
                                     {deftab.at("min_at")[k], deftab.at("max_at")[k]});
            house.init_augm_model_reshaping(i, a);

            if (house.is_augm_function(i, a))
                // assign arguments to function
                _assign_augm_func_args(house, i, a,
                                       tools::parse_list_args(deftab.at("farg_key")[k]),
                                       tools::parse_list_args(deftab.at("farg_val")[k]),
                                       tools::parse_list_args(deftab.at("farg_at")[k]));

            else if (!house.is_augm_constant(i, a))
            {
                // assign proposer rule
                house.assign_augm_proposer(i, a,
                                           deftab.at("rule")[k],
                                           _type_to_double(deftab.at("scale")[k], deftab.at("scale_at")[k], house, i),
                                           _type_to_double(deftab.at("from")[k], deftab.at("from_at")[k], house, i));

                // random initialization
                // UPDATE: I would like a more general switch for this
                if (tools::stobool(deftab.at("init_rand")[k]))
                    house.augm_model_rand(i, a, get_prng());
            }
        }
    }
    catch (const std::exception &e)
    {
        tools::error({"While initializing augmented variables:",
                      std::string(e.what())});
    }
}

/* Clear all augmented definition of individual. */
void MCMC::Impl::clear_all_augm(Household &house, indiv i)
{
    for (const augmID &a : augm_list_full)
        house.clear_augm(i, a);
}

//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//

/* Compute the loglikelihood term associated to augmented variables. */
double MCMC::Impl::compute_augm_model_loglike(const Household &house) const
{
    double hh_loglike = 0.;

    // Contribution of model associated to augmented variables
    for (const indiv &i : house.infected_status_list())
        (this->*_do_augm_model_by_definition)(house, i, hh_loglike);

    return hh_loglike;
}

void MCMC::Impl::__augm_model_by_status(const Household &house, indiv i, double &delta) const
{
    for (const augmID &a : augm_status_map_model.at(house.status(i)))
        delta += house.augm_model_logpdf(i, a);
}

void MCMC::Impl::__augm_model_generic(const Household &house, indiv i, double &delta) const
{
    for (const augmID &a : augm_list_model)
        delta += house.augm_model_logpdf(i, a);
}
