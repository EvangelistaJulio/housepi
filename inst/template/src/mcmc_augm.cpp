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

#include <atomic>

#include "mcmc_impl.h"

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
double MCMC::Impl::_type_to_double(const Individual &member,
                                   const std::string &value,
                                   const std::string &type)
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
            output = member.data(value);

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
void MCMC::Impl::_assign_augm_attributes(Individual &member,
                                         const augmID &a,
                                         const str_vec &keys,
                                         const str_vec &vals,
                                         const str_vec &typs,
                                         const std::vector<bool> &is_arg,
                                         bool is_var_model)
{
    // check size validity of keys, values and type
    if ((keys.size() != vals.size()) || (typs.size() != vals.size()))
        throw std::invalid_argument(
            "Mismatching attributes in definition of augmented variable \"" + a + "\".");

    Variable &xaugm = member.get_augm_variable(a);

    // loop over all attributes
    for (size_t k = 0; k < keys.size(); ++k)
    {
        // inputs
        std::string key = keys[k];
        std::string value = vals[k];
        std::string type = typs[k];
        bool flag = is_arg[k];

        // outputs
        double x = HOUSEPI_NA_DOUBLE;
        const double *y = nullptr;
        std::shared_ptr<Function> f = nullptr;

        if (tools::is_na(value) && tools::is_na(type))
            ; // pass

        else if (!tools::is_na(value) && tools::is_na(type))
            throw std::invalid_argument(
                "Unable to interpret default value of  \"" + value + "\". Please provide a type.");

        else if (tools::is_na(value) && !tools::is_na(type))
            tools::warning(
                "Found a type without a paired value. Falling back to default.");

        else if (tools::stoupper(type) == TYPE_SCALAR)
            x = tools::stod(value, "argument with scalar type");

        else if (tools::stoupper(type) == TYPE_DATA)
            x = member.data(value);

        else if (tools::stoupper(type) == TYPE_PARAM)
        {
            if (param_map.at(value).is_constant())
                x = get_param(value);
            else
            {
                y = get_param_pointer(value);
                if (is_var_model == VAR_FUNCTION)
                    make_param_cache_dependent(value);
            }
        }
        else if (tools::stoupper(type) == TYPE_AUGM)
        {
            if (contains(augm_list_full, value))
            {
                if (a == value)
                    throw std::invalid_argument(
                        "Unable to set argument equal to its value for \"" + a + "\".");
                else
                {
                    Variable &other_augm = member.get_augm_variable(value);

                    if (other_augm.is_function())
                        f = other_augm.get_function();
                    else
                        y = other_augm.get_pointer();
                }
            }
            else
                throw std::invalid_argument(
                    "Undeclared augmented variable \"" + value + "\".");
        }

        else
            throw std::invalid_argument(
                "Unsopported type \"" + type + "\" found while loading \"" + value + "\" for pointer value.");

        // final assignment
        if (!tools::is_na(x)) // scalar attribute
        {
            if (is_var_model == VAR_MODEL)
                xaugm.assign_distr_attr(key, flag, x);
            else
                xaugm.assign_func_attr(key, flag, x);
        }
        else if (y != nullptr) // pointer attribute
        {
            if (is_var_model == VAR_MODEL)
                xaugm.assign_distr_attr(key, flag, y);
            else
                xaugm.assign_func_attr(key, flag, y);
        }
        else if (f != nullptr) // function attribute
        {
            if (is_var_model == VAR_MODEL)
                xaugm.assign_distr_attr(key, flag, f);
            else
                xaugm.assign_func_attr(key, flag, f);
        }
    }
}

/* Assign augmented function arguments. */
void MCMC::Impl::_assign_augm_func_args(Individual &member,
                                        const augmID &a,
                                        const str_vec &keys,
                                        const str_vec &vals,
                                        const str_vec &typs)
{
    _assign_augm_attributes(member, a, keys, vals, typs,
                            std::vector<bool>(keys.size(), ATTR_ARGUMENT), VAR_FUNCTION);
}
/* Assign augmented function range. */
void MCMC::Impl::_assign_augm_func_range(Individual &member,
                                         const augmID &a,
                                         const str_vec &vals,
                                         const str_vec &typs)
{
    _assign_augm_attributes(member, a, {"min", "max"}, vals, typs,
                            std::vector<bool>(2, ATTR_RANGE), VAR_FUNCTION);
}

/* Assign augmented model arguments. */
void MCMC::Impl::_assign_augm_model_args(Individual &member,
                                         const augmID &a,
                                         const str_vec &keys,
                                         const str_vec &vals,
                                         const str_vec &typs)
{
    _assign_augm_attributes(member, a, keys, vals, typs,
                            std::vector<bool>(keys.size(), ATTR_ARGUMENT), VAR_MODEL);
}
/* Assign augmented model range. */
void MCMC::Impl::_assign_augm_model_range(Individual &member,
                                          const augmID &a,
                                          const str_vec &vals,
                                          const str_vec &typs)
{
    _assign_augm_attributes(member, a, {"min", "max"}, vals, typs,
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
                        "Unable to associate augmented \"" + a + "\" for non infectious status.");
            }
            else if (mylist.size() > 1)
                throw std::invalid_argument(
                    "Conflicting status-specific and general definitions for augmented \"" + a + "\".");

            augm_list_full.push_back(a);
        }
    }
}

/* Search individual id for augmented traces follow-up. */
void MCMC::Impl::_build_augm_set_id_trace()
{
    if (is_augm_sample())
    {
        if (!tools::is_na(setup_id_trace))
        {
            for (Household &house : house_vector)
                for (Individual &member : house)
                    if (member.id() == setup_id_trace)
                        smpl_individual_star = &member;

            if (smpl_individual_star == nullptr)
            {
                tools::warning(
                    "Ignoring id trace: no individual \"" + setup_id_trace + "\" was found.");
                setup_id_trace = HOUSEPI_NA_STRING;
            }
            else if (smpl_individual_star->status() < INFECTED_STATUS)
                tools::warning(
                    "Individual id \"" + setup_id_trace + "\" correspond to a healthy individual.");
        }
    }
    else if (!tools::is_na(setup_id_trace))
    {
        tools::warning(
            "Ignoring id trace: no augmented variables to infer.");
        setup_id_trace = HOUSEPI_NA_STRING;
    }
}

/* Assign augmented variables to all individuals. */
void MCMC::Impl::build_augm()
{
    try
    {
        // eventually append infectivity profile definition
        _infec_profile.append_augm_definition();

        const dict_str_vec &deftable = _augm_deftable;

        if (!deftable.empty())
        {
            // check entries
            str_vec headers = tools::headers_from_deftable(deftable);
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
                                    "Unable to initialize \"" + a + "\" without a model.");

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

            // end assigment in the households
            _build_augm_set_id_trace();

            // Allocate memory for augmented variables to each individual.
            for (Household &house : house_vector)
                for (Individual &member : house)
                    clear_all_augm(member);
        }
    }
    catch (const std::exception &e)
    {
        tools::error("While building augmented variable:\n" +
                     std::string(e.what()));
    }
}

/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////  MCMC::Impl : assign augmented variables  /////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

/* Assign augmented variables to single individual. */
void MCMC::Impl::assign_all_augm(Individual &member)
{
    try
    {
        int status = member.status();
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
                member.assign_augm(a,
                                   deftab.at("name")[k],
                                   _type_to_double(member, deftab.at("value")[k], deftab.at("value_at")[k]),
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
            Variable &xaugm = member.get_augm_variable(a);

            // assign arguments to model
            _assign_augm_model_args(member, a,
                                    tools::parse_list_args(deftab.at("arg_key")[k]),
                                    tools::parse_list_args(deftab.at("arg_val")[k]),
                                    tools::parse_list_args(deftab.at("arg_at")[k]));

            // assign range to model
            _assign_augm_model_range(member, a,
                                     {deftab.at("min")[k], deftab.at("max")[k]},
                                     {deftab.at("min_at")[k], deftab.at("max_at")[k]});
            xaugm.init_distr_reshaping();

            if (xaugm.is_function())
            {
                // assign arguments to function
                _assign_augm_func_args(member, a,
                                       tools::parse_list_args(deftab.at("farg_key")[k]),
                                       tools::parse_list_args(deftab.at("farg_val")[k]),
                                       tools::parse_list_args(deftab.at("farg_at")[k]));
                // assign range to function
                _assign_augm_func_range(member, a,
                                        {deftab.at("min")[k], deftab.at("max")[k]},
                                        {deftab.at("min_at")[k], deftab.at("max_at")[k]});
                xaugm.init_func_reshaping();
            }
            else if (!xaugm.is_constant())
            {
                // assign proposer rule
                xaugm.assign_proposer(deftab.at("rule")[k],
                                      _type_to_double(member, deftab.at("scale")[k], deftab.at("scale_at")[k]),
                                      _type_to_double(member, deftab.at("from")[k], deftab.at("from_at")[k]));

                // random initialization
                // UPDATE: I would like a more general switch for this
                if (tools::stobool(deftab.at("init_rand")[k]))
                    xaugm.rand(get_prng());
            }
        }

        // initialize functional augmented model for individuals
        if (is_augm_generic())
            member.init_augm_model(augm_list_model);
        else
            member.init_augm_model(augm_status_map_model[status]);
    }
    catch (const std::exception &e)
    {
        throw std::runtime_error(
            "While initializing augmented variables:\n" + std::string(e.what()));
    }
}

/* Clear all augmented definition of individual. */
void MCMC::Impl::clear_all_augm(Individual &member)
{
    for (const augmID &a : augm_list_full)
        member.clear_augm(a);
}

///////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////  MCMC : init augmented variables  ////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/* Initialize augmented dependecies for all households. */
void MCMC::Impl::_init_augmented_variables()
{
    if (is_augm_any())
    {
        if (is_augm_sample())
        {
            // wheter augmented variable are generic or status specific
            if (is_augm_generic())
                _do_move_augm_indiv = &MCMC::Impl::__move_augm_indiv_generic;
            else
                _do_move_augm_indiv = &MCMC::Impl::__move_augm_indiv_by_status;
        }

        // assign augmented variables
        for (Household &house : house_vector)
            for (Individual &member : house)
                assign_all_augm(member);
    }
}

//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//

/* Compute loglikelihood term associated to augmented models. */
double MCMC::Impl::_compute_augm_model_loglike(Household &house)
{
    double hh_loglike = 0.;
    for (const indiv &i : house.get_infected_status_hhindexes())
        hh_loglike += house.get_member(i).get_augm_model_loglike();
    return hh_loglike;
}

//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////  MCMC : move all augmented variables  /////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/* Move augmented variables from all households. */
void MCMC::Impl::_move_all_augm()
{
    if (is_augm_sample())
    {
        double delta_ll = 0.;
        std::atomic<bool> has_error(false);

#pragma omp parallel for reduction(+ : delta_ll)
        for (size_t H = 0; H < n_houses(); ++H)
        {
            if (has_error.load(std::memory_order_relaxed))
                continue;

            try
            {
                delta_ll += _move_augm_single_house(H);
            }
            catch (const std::exception &e)
            {
#pragma omp critical
                {
                    if (!has_error.exchange(true))
                        _error_message = "While sampling augmented variables for household " +
                                         house_vector[H].house_id() + ":\n" +
                                         std::string(e.what());
                }
            }
        }

        if (!_error_message.empty())
            throw std::runtime_error(_error_message);

        increment_full_loglike(delta_ll);
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////  MCMC : augmented move  //////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/* Move augmented variables of infectious members inside household.
The move is performed one member at a time, one augmented variable at a time.
*/
// UPDATE: smooth the incrementation by reference by just adding...
double MCMC::Impl::_move_augm_single_house(size_t H)
{
    double delta_ll = 0.;
    Household &house = house_vector[H];
    // Loop over infectious residents of the househould
    for (const indiv &i : house.get_infected_status_hhindexes())
        delta_ll += (this->*_do_move_augm_indiv)(house, i);
    return delta_ll;
}

/* Generic move augmented variables one by one. */
double MCMC::Impl::__move_augm_indiv_generic(Household &house, indiv i)
{
    double delta_ll = 0.;
    // Loop over augmented variables of infectious
    for (const augmID &a : augm_list_infer)
        delta_ll += _metropolis_augm_single_move(house, i, a);
    return delta_ll;
}

/* Move augmented variables by status one variable at a time. */
double MCMC::Impl::__move_augm_indiv_by_status(Household &house, indiv i)
{
    double delta_ll = 0.;
    // Loop over augmented variables of infectious
    for (const augmID &a : augm_status_map_infer.at(house.status(i)))
        delta_ll += _metropolis_augm_single_move(house, i, a);
    return delta_ll;
}

/* Take care of a a single augmented Metropolis move. */
double MCMC::Impl::_metropolis_augm_single_move(Household &house, indiv i, const augmID &a)
{
    Variable &xaugm = house.get_member(i).get_augm_variable(a);

    xaugm.propose(get_prng()); // propose new value and buffer old

    // proceed only if the new value is consistent
    if (xaugm.valid_move())
    {
        // update dependencies to perform computations
        house.buffer_all_cache();
        _update_cache_dependecies(house, i);

        double new_obs_ll = _compute_house_obs_loglike(house);
        // UPDATE: I don't need to recompute for all individuals
        double new_augm_value = _compute_augm_model_loglike(house);

        // Metropolis-Hastings algorithm
        const double D_ll = new_obs_ll + new_augm_value - house.get_obs_loglike() - house.get_augm_loglike();
        if (uniform() < std::exp(setup_inv_T * D_ll - xaugm.propose_logratio()))
        {
            house.set_obs_loglike(new_obs_ll);
            house.set_augm_loglike(new_augm_value);
            return D_ll;
        }
        else
            // else recover cached values
            house.recover_all_cache();
    }
    // else recover buffered value
    xaugm.recover_value();

    return 0.;
}
