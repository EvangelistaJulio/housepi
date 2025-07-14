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
double MCMC::Impl::ImplAugmented::type_to_double_(const Individual &member,
                                                  const std::string &value,
                                                  const std::string &type)
{
    double output = HSP_NA_DOUBLE;

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
            if (pImpl_.param.dictionary.at(value).is_constant())
                output = pImpl_.param.get(value);
            else
                throw std::invalid_argument(
                    "Parameter to be inferred \"" + value + "\" assigned as constant argument.");

        else if (tools::stoupper(type) == TYPE_AUGM)
            throw std::invalid_argument(
                "Value initialization with other augmented variable is not supported.");

        else
            throw std::invalid_argument(
                "Unknown variable type specifier \"" + type + "\".");
    }

    return output;
}

/* Auxiliary method to assign attributes to augmented function or model. */
void MCMC::Impl::ImplAugmented::init_attributes_(Individual &member,
                                                   const augm_key &a,
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
        double x = HSP_NA_DOUBLE;
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
            if (pImpl_.param.dictionary.at(value).is_constant())
                x = pImpl_.param.get(value);
            else
            {
                y = pImpl_.param.get_ptr(value);
                if (is_var_model == HSP_VAR_FUNCTION)
                    pImpl_.param.make_cache_dependent(value);
            }
        }
        else if (tools::stoupper(type) == TYPE_AUGM)
        {
            if (list_full % value)
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
                        y = other_augm.get_ptr();
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
            if (is_var_model == HSP_VAR_MODEL)
                xaugm.init_distr_attr(key, flag, x);
            else
                xaugm.init_func_attr(key, flag, x);
        }
        else if (y != nullptr) // pointer attribute
        {
            if (is_var_model == HSP_VAR_MODEL)
                xaugm.init_distr_attr(key, flag, y);
            else
                xaugm.init_func_attr(key, flag, y);
        }
        else if (f != nullptr) // function attribute
        {
            if (is_var_model == HSP_VAR_MODEL)
                xaugm.init_distr_attr(key, flag, f);
            else
                xaugm.init_func_attr(key, flag, f);
        }
    }
}

/* Assign augmented function arguments. */
void MCMC::Impl::ImplAugmented::init_func_args_(Individual &member,
                                                  const augm_key &a,
                                                  const str_vec &keys,
                                                  const str_vec &vals,
                                                  const str_vec &typs)
{
    init_attributes_(member, a, keys, vals, typs,
                       std::vector<bool>(keys.size(), HSP_ATTR_ARGUMENT), HSP_VAR_FUNCTION);
}
/* Assign augmented function range. */
void MCMC::Impl::ImplAugmented::init_func_range_(Individual &member,
                                                   const augm_key &a,
                                                   const str_vec &vals,
                                                   const str_vec &typs)
{
    init_attributes_(member, a, {"min", "max"}, vals, typs,
                       std::vector<bool>(2, HSP_ATTR_RANGE), HSP_VAR_FUNCTION);
}

/* Assign augmented model arguments. */
void MCMC::Impl::ImplAugmented::init_model_args_(Individual &member,
                                                   const augm_key &a,
                                                   const str_vec &keys,
                                                   const str_vec &vals,
                                                   const str_vec &typs)
{
    init_attributes_(member, a, keys, vals, typs,
                       std::vector<bool>(keys.size(), HSP_ATTR_ARGUMENT), HSP_VAR_MODEL);
}
/* Assign augmented model range. */
void MCMC::Impl::ImplAugmented::init_model_range_(Individual &member,
                                                    const augm_key &a,
                                                    const str_vec &vals,
                                                    const str_vec &typs)
{
    init_attributes_(member, a, {"min", "max"}, vals, typs,
                       std::vector<bool>(2, HSP_ATTR_RANGE), HSP_VAR_MODEL);
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
const dict_str_vec MCMC::Impl::ImplAugmented::check_name_status_(const dict_str_vec &deftable)
{
    const size_t table_length = deftable.at("name").size();

    dict_str_vec output = {};

    // checka augmented variables lists of usage
    for (size_t k = 0; k < table_length; ++k)
    {
        const std::string &a = deftable.at("name")[k];
        const std::string &status = deftable.at("status")[k];

        // check that augm_key and status are correct
        if (tools::is_na(a))
            throw std::invalid_argument(
                "Augmented variable name can not be undefined.");

        else if (output.count(a))
        {
            // check no multiple augm_key with same status
            const str_vec &mylist = output.at(a);
            if (mylist % status)
                throw std::invalid_argument(
                    "Incompatible definitions for augmented variable \"" + a + "\".");
            else
                output.at(a).push_back(status);
        }

        else // new augmented variable status
            output.emplace(a, str_vec{status});
    }

    return output;
}

/* Check univocal definitions and allocate new lists. */
void MCMC::Impl::ImplAugmented::allocate_univocal_defs_(const dict_str_vec &check)
{
    // drop previous definitions
    list_full.clear();
    status_map_infer.clear();
    status_map_model.clear();

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
                    status_map_infer.emplace(s, std::vector<augm_key>{});
                    status_map_model.emplace(s, std::vector<augm_key>{});
                }
                else
                    throw std::invalid_argument(
                        "Unable to associate augmented \"" + a + "\" for non infectious status.");
            }
            else if (mylist.size() > 1)
                throw std::invalid_argument(
                    "Conflicting status-specific and general definitions for augmented \"" + a + "\".");

            list_full.push_back(a);
        }
    }
}

/* Assign augmented variables to all individuals. */
void MCMC::Impl::ImplAugmented::build()
{
    try
    {
        // eventually append infectivity profile definition
        pImpl_.infec_profile.append_augm_definition();

        if (!deftable.empty())
        {
            // check entries
            str_vec headers = tools::headers_from_deftable(deftable);
            tools::check_for_missing_entries(headers, file_headers);
            tools::check_for_unused_entries(headers, file_headers);

            // auxiliary variable to load definitions : a -> status lists
            const dict_str_vec name_status_check = check_name_status_(deftable);
            allocate_univocal_defs_(name_status_check);

            //===========================================//
            //===== Assign new definitions to lists =====//
            //===========================================//

            // drop previous definitions
            list_sample.clear();
            list_func.clear();
            list_model.clear();

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
                    list_func.push_back(a);

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
                        if (!(list_sample % a))
                            list_sample.push_back(a);

                        // note: impossible to associate rule to missing status
                        if (!tools::is_na(status)) // status-specific definition
                        {
                            int s = tools::stoi(status, "augmented status");
                            if (!(status_map_infer.at(s) % a))
                                status_map_infer.at(s).push_back(a);
                        }

                        else if (!is_generic()) // generic definition
                            for (const auto &pair : status_map_infer)
                                if (!(status_map_infer.at(pair.first) % a))
                                    status_map_infer.at(pair.first).push_back(a);
                    }
                }

                // Augmented Variables with a Model
                if (!tools::is_na(deftable.at("model")[k]))
                {
                    if (!(list_model % a))
                        list_model.push_back(a);

                    if (!tools::is_na(status)) // status-specific definition
                    {
                        int s = tools::stoi(status, "augmented status");
                        if (!(status_map_model.at(s) % a))
                            status_map_model.at(s).push_back(a);
                    }

                    else if (!is_generic()) // generic definition
                        for (const auto &pair : status_map_model)
                            if (!(status_map_model.at(pair.first) % a))
                                status_map_model.at(pair.first).push_back(a);
                }
            }

            // Allocate memory for augmented variables to each individual.
            for (Household &house : pImpl_.dataset)
                for (Individual &member : house)
                    clear_all(member);

            // end assigment in the households
            set_id_trace_(pImpl_.setup.id_trace);
        }
    }
    catch (const std::exception &e)
    {
        tools::error("While building augmented variable:\n" +
                     std::string(e.what()));
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////  MCMC : Impl : ImplAugmented : sample traces  //////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/* Search individual id for augmented traces follow-up. */
void MCMC::Impl::ImplAugmented::set_id_trace_(str_vec &id_trace)
{
    smpl_individuals.clear();

    if (!id_trace.empty())
    {
        if (is_sample())
        {
            smpl_individuals.clear();
            smpl_individuals.reserve(id_trace.size());

            // check id across all individuals
            str_vec already_stored_id = {};
            for (const Household &house : pImpl_.dataset)
                for (const Individual &member : house)
                    if (id_trace % member.id())
                    {
                        if (already_stored_id % member.id())
                            throw std::invalid_argument(
                                "Multiple individuals share the same id trace \"" + member.id() + "\".");
                        else
                        {
                            already_stored_id.push_back(member.id());
                            smpl_individuals.push_back(&member);
                        }
                    }

            // check all user request are satisfied
            tools::check_for_missing_entries(already_stored_id, id_trace, "trace id");

            if (is_traces())
            {
                for (const Individual *member : smpl_individuals)
                    if (!member->is_infected_status())
                        tools::warning(
                            "Ignornig trace of non infected individual id \"" + member->id() + "\".");
            }
            else
            {
                tools::warning(
                    "Ignoring id trace for individual " + tools::str_list(id_trace));
                id_trace.clear();
            }
        }
        else
        {
            tools::warning(
                "Ignoring id trace: no augmented variables to infer.");
            id_trace.clear();
        }
    }
}

/* Initialize sample for individual augmented variable followup. */
void MCMC::Impl::ImplAugmented::init_traces(ulong sample_size)
{
    if (is_traces())
    {
        smpl_values.clear();
        smpl_values.resize(smpl_individuals.size());
        for (size_t k = 0; k < smpl_individuals.size(); ++k)
        {
            smpl_values[k].clear();
            smpl_values[k].reserve((list_sample + list_func).size());
            for (const augm_key &a : list_sample + list_func)
            {
                smpl_values[k][a].clear();
                smpl_values[k][a].reserve(sample_size);
                smpl_values[k][a].assign(sample_size, HSP_NA_DOUBLE);
            }
        }
    }
}

/* Store sample for individual augmented variable followup */
void MCMC::Impl::ImplAugmented::update_traces(ulong smpl_count)
{
    if (is_traces())
        for (size_t k = 0; k < smpl_individuals.size(); ++k)
        {
            auto &smpl_augm_value_k = smpl_values[k];
            const Individual *member_star = smpl_individuals[k];
            for (const augm_key &a : list_sample + list_func)
                smpl_augm_value_k[a][smpl_count] = member_star->augm(a);
        }
}

/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////  MCMC::Impl : assign augmented variables  /////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

/* Clear all augmented definition of individual. */
void MCMC::Impl::ImplAugmented::clear_all(Individual &member)
{
    for (const augm_key &a : list_full)
        member.clear_augm(a);
}

/* Assign augmented variables to single individual. */
void MCMC::Impl::ImplAugmented::init_all(Individual &member)
{
    try
    {
        int status = member.status();
        const dict_str_vec &deftab = deftable;
        std::map<augm_key, size_t> indiv_augm_entries;

        // loop over augmented variable definitions to assign placheholders
        for (size_t k = 0; k < deftab.at("name").size(); ++k)
        {
            // assign augmented variable of the corresponding status
            const std::string &augm_status = deftab.at("status")[k];

            // assign augmented variable if undefined status | corresponding status
            if (tools::is_na(augm_status) || (status == tools::stoi(augm_status, "augmented status")))
            {
                augm_key a = deftab.at("name")[k]; // UPDATE : index to be implemented

                // assign augmented variable placeholder
                // FIXME: what if constant = T, discrete = T and value != integer?
                member.init_augm(a,
                                   deftab.at("name")[k],
                                   type_to_double_(member, deftab.at("value")[k], deftab.at("value_at")[k]),
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
            const augm_key &a = pair.first;
            const size_t &k = pair.second;
            Variable &xaugm = member.get_augm_variable(a);

            // assign arguments to model
            init_model_args_(member, a,
                               tools::parse_list_args(deftab.at("arg_key")[k]),
                               tools::parse_list_args(deftab.at("arg_val")[k]),
                               tools::parse_list_args(deftab.at("arg_at")[k]));

            // assign range to model
            init_model_range_(member, a,
                                {deftab.at("min")[k], deftab.at("max")[k]},
                                {deftab.at("min_at")[k], deftab.at("max_at")[k]});
            xaugm.init_distr_reshaping();

            if (xaugm.is_function())
            {
                // assign arguments to function
                init_func_args_(member, a,
                                  tools::parse_list_args(deftab.at("farg_key")[k]),
                                  tools::parse_list_args(deftab.at("farg_val")[k]),
                                  tools::parse_list_args(deftab.at("farg_at")[k]));
                // assign range to function
                init_func_range_(member, a,
                                   {deftab.at("min")[k], deftab.at("max")[k]},
                                   {deftab.at("min_at")[k], deftab.at("max_at")[k]});
                xaugm.init_func_reshaping();
            }
            else if (!xaugm.is_constant())
            {
                // assign proposer rule
                xaugm.init_proposer(deftab.at("rule")[k],
                                      type_to_double_(member, deftab.at("scale")[k], deftab.at("scale_at")[k]),
                                      type_to_double_(member, deftab.at("from")[k], deftab.at("from_at")[k]));

                // random initialization
                // UPDATE: I would like a more general switch for this
                if (tools::stobool(deftab.at("init_rand")[k]))
                    xaugm.rand(pImpl_.setup.prng());
            }
        }

        // initialize functional augmented model for individuals
        if (is_generic())
            member.init_augm_model(list_model);
        else
            member.init_augm_model(status_map_model[status]);
    }
    catch (const std::exception &e)
    {
        throw std::runtime_error(
            "While initializing augmented variables:\n" + std::string(e.what()));
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////  MCMC : init augmented variables  ////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/* Initialize augmented dependecies for all households. */
void MCMC::Impl::ImplAugmented::init()
{
    if (is_any())
        for (Household &house : pImpl_.dataset)
            for (Individual &member : house)
                init_all(member);
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
void MCMC::Impl::move_all_augm_()
{
    if (augm.is_sample())
    {
        double delta_ll = 0.;
        std::atomic<bool> has_error(false);

#pragma omp parallel for reduction(+ : delta_ll)
        for (size_t H = 0; H < dataset.n_houses(); ++H)
        {
            Household &this_house = dataset(H);

            if (has_error.load(std::memory_order_relaxed))
                continue;

            try
            {
                delta_ll += move_augm_single_house_(this_house);
            }
            catch (const std::exception &e)
            {
#pragma omp critical
                {
                    if (!has_error.exchange(true))
                        error_message = "While sampling augmented variables for household " +
                                        this_house.house_id() + ":\n" +
                                        std::string(e.what());
                }
            }
        }

        if (!error_message.empty())
            throw std::runtime_error(error_message);

        increment_full_loglike_(delta_ll);
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////  MCMC : augmented move  //////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/* Move augmented variables of infectious members inside household.
The move is performed one member at a time, one augmented variable at a time.
*/
// UPDATE: smooth the incrementation by reference by just adding...
double MCMC::Impl::move_augm_single_house_(Household &house)
{
    double delta_ll = 0.;
    // Loop over infectious residents of the househould
    for (Individual &member_i : house.infected())
        delta_ll += (this->*move_augm_indiv_)(house, member_i);
    return delta_ll;
}

/* Generic move augmented variables one by one. */
double MCMC::Impl::do_move_augm_indiv_generic_(Household &house, Individual &member_i)
{
    double delta_ll = 0.;
    // Loop over augmented variables of infectious
    for (const augm_key &a : augm.list_sample)
        delta_ll += metropolis_augm_single_move_(house, member_i, a);
    return delta_ll;
}

/* Move augmented variables by status one variable at a time. */
double MCMC::Impl::do_move_augm_indiv_by_status_(Household &house, Individual &member_i)
{
    double delta_ll = 0.;
    // Loop over augmented variables of infectious
    for (const augm_key &a : augm.status_map_infer.at(member_i.status()))
        delta_ll += metropolis_augm_single_move_(house, member_i, a);
    return delta_ll;
}

/* Take care of a a single augmented Metropolis move. */
double MCMC::Impl::metropolis_augm_single_move_(Household &house, Individual &member, const augm_key &a)
{
    Variable &xaugm = member.get_augm_variable(a);

    xaugm.propose(setup.prng()); // propose new value and buffer old

    // proceed only if the new value is consistent
    if (xaugm.valid_move())
    {
        // update dependencies to perform computations
        house.buffer_all_cache();
        update_cache_dependecies_(house, member);

        // Metropolis-Hastings algorithm
        double D_ll = -(house.get_obs_loglike() + member.get_augm_model_loglike());
        D_ll += member.compute_and_set_augm_model_loglike();
        D_ll += compute_and_set_household_obs_loglike_(house);

        if (setup.unif() < std::exp(setup.inv_T * D_ll - xaugm.propose_logratio()))
            return D_ll;
        else
        {
            // else recover cached values
            house.recover_buffer_values(true);
            member.recover_augm_model_loglike();
        }
    }
    // else recover buffered value
    xaugm.recover_value();

    return 0.;
}
