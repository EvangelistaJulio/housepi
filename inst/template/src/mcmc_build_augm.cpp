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

    return output;
}

/* Convert a std::string with type-dictionary to std::string-to-double */
void MCMC::Impl::_assign_augm_model(Household &house, indiv i, augmID a,
                                    const str_vec &keys,
                                    const str_vec &vals,
                                    const str_vec &typs,
                                    const std::vector<bool> &is_arg)
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
            house.assign_augm_model_scalar(i, a, key, flag, house.data(i, value));

        else if (tools::stoupper(type) == TYPE_PARAM)
        {
            if (param_map.at(value).is_constant())
                house.assign_augm_model_scalar(i, a, key, flag, get_param(value));
            else
                house.assign_augm_model_pointer(i, a, key, flag, get_param_pointer(value));
        }

        else if (tools::stoupper(type) == TYPE_SCALAR)
            house.assign_augm_model_scalar(i, a, key, flag, tools::stod(value, "argument with scalar type"));

        else if (tools::stoupper(type) == TYPE_AUGM)
        {

            if (augm_list_full.count(value))
            {
                if (a == value)
                    throw std::invalid_argument(
                        "Impossible to set argument equal to its value for \"" + a + "\".");
                else if (house.is_augm_function(i, value))
                    throw std::invalid_argument(
                        "Impossible to set augmented function \"" + value + "\" as argument.");
                else if (house.is_augm_constant(i, value))
                    house.assign_augm_model_scalar(i, a, key, flag, house.augm(i, value));
                else
                    house.assign_augm_model_pointer(i, a, key, flag, house.augm_pointer(i, value));
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
            if (std::find(mylist.begin(), mylist.end(), status) == mylist.end())
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
                    augm_status_map_infer.emplace(s, std::unordered_set<std::string>{});
                    augm_status_map_model.emplace(s, std::unordered_set<std::string>{});
                }
                else
                    throw std::invalid_argument(
                        "Impossible to associate augmented \"" + a + "\" for non infectious status.");
            }
            else if (mylist.size() > 1)
                throw std::invalid_argument(
                    "Conflicting status-specific and general definitions for augmented \"" + a + "\".");

            augm_list_full.emplace(a);
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
                            smpl_H_star = house.index();
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
            augm_infer_list.clear();
            augm_list_model.clear();

            // note : I didn't manage to collapse this with the previous because of generic status
            for (size_t k = 0; k < deftable.at("name").size(); ++k)
            {
                const std::string &a = deftable.at("name")[k];
                const std::string &status = deftable.at("status")[k];

                // read option for constant variable
                bool is_constant = false;
                if (!tools::is_na(deftable.at("constant")[k]))
                    is_constant = tools::stobool(deftable.at("constant")[k]);

                // UPDATE: here I should check if variable is a function
                // Augmented Variables to be Inferred
                if (!is_constant)
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
                        augm_infer_list.emplace(a);

                        if (!tools::is_na(status)) // status-specific definition
                            augm_status_map_infer.at(tools::stoi(status, "augmented status")).emplace(a);

                        else if (!is_augm_generic()) // generic definition
                            for (const auto &pair : augm_status_map_infer)
                                augm_status_map_infer.at(pair.first).emplace(a);
                    }
                }

                //----- Augmented Variables with a Model -----//
                if (!tools::is_na(deftable.at("model")[k]))
                {
                    augm_list_model.emplace(a);

                    if (!tools::is_na(status)) // status-specific definition
                        augm_status_map_model.at(tools::stoi(status)).emplace(a);

                    else if (!is_augm_generic()) // generic definition
                        for (const auto &pair : augm_status_map_model)
                            augm_status_map_model.at(pair.first).emplace(a);
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
bool MCMC::Impl::is_augm_sample() const { return !augm_infer_list.empty(); }

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
        const dict_str_vec &deftab = _augm_deftable;

        // loop over all augmented variable definitions
        for (size_t k = 0; k < deftab.at("name").size(); ++k)
        {
            // assign augmented variable of the corresponding status
            const std::string &augm_status = deftab.at("status")[k];

            // assign augmented variable if undefined status | corresponding status
            if ((tools::is_na(augm_status)) || (house.status(i) == tools::stoi(augm_status, "augmented status")))
            {
                augmID a = deftab.at("name")[k];                // UPDATE : index to be implemented
                const std::string &name = deftab.at("name")[k]; // augmented variable name

                // initial value
                double value = _type_to_double(deftab.at("value")[k], deftab.at("value_at")[k], house, i);
                bool is_constant = tools::stobool(deftab.at("constant")[k]);
                bool is_discrete = tools::stobool(deftab.at("discrete")[k]);
                const std::string &model_name = deftab.at("model")[k];

                // FIXME: what if constant = T, discrete = T and value != int?

                // assign augmented variable
                house.assign_augm(i, a, name, value, is_discrete, model_name);

                // load model arguments and range
                str_vec keys = tools::parse_list_args(deftab.at("arg_name")[k]);
                str_vec vals = tools::parse_list_args(deftab.at("arg_value")[k]);
                str_vec typs = tools::parse_list_args(deftab.at("arg_at")[k]);
                std::vector<bool> is_arg(keys.size(), ARG_DEFAULT);

                keys += {"min", "max"};
                vals += {deftab.at("min")[k], deftab.at("max")[k]};
                typs += {deftab.at("min_at")[k], deftab.at("max_at")[k]};
                is_arg += {ARG_RANGE, ARG_RANGE};

                _assign_augm_model(house, i, a, keys, vals, typs, is_arg);
                house.init_augm_model_reshaping(i, a);

                // non constant definitions
                if (!is_constant)
                {
                    // assign proposer rule
                    const std::string &rule_name = deftab.at("rule")[k];
                    double rule_scale = _type_to_double(deftab.at("scale")[k], deftab.at("scale_at")[k], house, i);
                    double rule_from = _type_to_double(deftab.at("from")[k], deftab.at("from_at")[k], house, i);
                    house.assign_augm_proposer(i, a, rule_name, rule_scale, rule_from);

                    // random initialization
                    // UPDATE: I would like a more general switch for this
                    if (tools::stobool(deftab.at("init_rand")[k]))
                        house.augm_model_rand(i, name, get_prng());
                }
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
