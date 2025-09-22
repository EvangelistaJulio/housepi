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
 * mcmc_build_augm.cpp
 *
 */

#include <atomic>

#include "mcmc_impl.h"

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////  MCMC : Impl : Augmented : build  ///////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/// Assign augmented variables to all individuals.
void MCMC::Impl::Augmented::build()
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

            // assert that augmented definitions are univocal
            build_helper_assert_univocal_defs_();
            // build status-specific maps of augmented variables
            build_helper_status_specific_maps_();
            // set individual id for augmented draws follow-up
            build_helper_set_id_draws_(pImpl_.setup.id_draws);

            // allocate memory for augmented variables to each individual
            for (Household &house : pImpl_.dataset)
                for (Individual &member : house)
                    clear_all(member);
        }
    }
    catch (const std::exception &e)
    {
        tools::error("While building augmented variable:\n" +
                     std::string(e.what()));
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////  MCMC : Impl : Augmented : assert univocal definitions  ///////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/// Assert that augmented definitions are univocal.
void MCMC::Impl::Augmented::build_helper_assert_univocal_defs_()
{
    // check augmented variables lists of usage
    dict_str_vec name_status_check = {};
    for (size_t k = 0; k < deftable.at("name").size(); ++k)
    {
        const std::string &a = deftable.at("name")[k];
        const std::string &s = deftable.at("status")[k]; // NOTE: can be NA, i. e. generic

        // check augmented name
        if (tools::is_na(a))
            throw std::invalid_argument(
                "Augmented variable name cannot be undefined.");
        // check that status-specific definition is for infected status
        if (!tools::is_na(s))
            if (!HSP_IS_INFECTED(tools::str::toi(s, "augmented status")))
                throw std::invalid_argument(
                    "Unable to associate augmented \"" + a + "\" for non infectious status " + s + ".");

        // load names and statuses
        name_status_check[a].push_back(s);
    }

    // check that name and status are univocal
    for (const auto &pair : name_status_check)
        if (pair.second.size() != tools::str::unique(pair.second).size())
            throw std::invalid_argument(
                "Incompatible definitions for augmented \"" + pair.first + "\".");
        else if (pair.second.size() > 1 && tools::is_any_na(pair.second))
            throw std::invalid_argument(
                "Conflicting generic and status-specific definitions for augmented \"" + pair.first + "\".");
}

///////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////  MCMC : Impl : Augmented : build function  //////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

void MCMC::Impl::Augmented::build_helper_status_specific_maps_()
{
    // drop previous definitions
    list_sample.clear();
    list_const.clear();
    list_func.clear();
    list_model.clear();
    list_draws.clear();
    status_map_sample.clear();
    status_map_model.clear();

    for (size_t k = 0; k < deftable.at("name").size(); ++k)
    {
        build_helper_variable_with_model_(k);

        if (!tools::is_na(deftable.at("func")[k]))
            build_function_(k);
        else if (tools::str::tobool(deftable.at("const")[k]))
            list_const.push_back(deftable.at("name")[k]);
        else
            build_helper_variable_to_sample_(k);
    }

    // finalize augmented lists
    list_sample = tools::str::unique(list_sample);
    list_const = tools::str::unique(list_const);
    list_func = tools::str::unique(list_func);
    list_model = tools::str::unique(list_model);
    list_full = tools::str::unique(list_sample + list_func + list_const);
    list_export = tools::str::drop_private(list_full);
}

///////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////  MCMC : Impl : Augmented : build functions  /////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/// Deal with the definition of augmented function.
/// UPDATE: improve implementation of constant function definitions
/// NOTE: under development
void MCMC::Impl::Augmented::build_function_(size_t k)
{
    const std::string &func = deftable.at("func")[k];

    if (func == HSP_IDENTITY_FUNC) // optimize 'identity' function
        build_function_identity_helper_(k);
    else
        list_func.push_back(deftable.at("name")[k]);
}

/// Build augmented equality function and change deftable.
void MCMC::Impl::Augmented::build_function_identity_helper_(size_t k)
{
    const std::string &a = deftable.at("name")[k];
    str_vec vals = tools::parse_arg_list(deftable.at("farg_val")[k], HSP_LIST_DELIM);
    str_vec types = tools::parse_arg_list(deftable.at("farg_at")[k], HSP_LIST_DELIM);
    // note: emptiness is checked elsewhere

    if (types.size() != 1 || vals.size() != 1)
        throw std::invalid_argument(
            "Augmented equality function \"" + a + "\" has incorrect arguments.");
    else
    {
        std::string identity_type = tools::str::toup(types)[0];
        bool is_identity_const = true;

        // check the type of requested equality
        if (tools::str::toup(identity_type) == HSP_TYPE_PARAM)
            is_identity_const = pImpl_.param.at(vals[0]).is_constant();
        else if (tools::str::toup(identity_type) == HSP_TYPE_AUGM)
            is_identity_const = false; // UPDATE: use an augmented dependecies tree to see if constant

        if (is_identity_const)
            list_const.push_back(a);
        else
            list_func.push_back(a);
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////  MCMC : Impl : Augmented : build variables  /////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/// Deal with the definition of augmented variable with a model.
void MCMC::Impl::Augmented::build_helper_variable_with_model_(size_t k)
{
    if (!tools::is_na(deftable.at("model")[k]))
    {
        const std::string &a = deftable.at("name")[k];
        const std::string &s = deftable.at("status")[k];

        list_model.push_back(a); // NOTE: make unique elsewhere

        if (tools::is_na(s)) // generic definition
            for (int status : pImpl_.dataset.all_infected_status_list)
                status_map_model[status].push_back(a);
        else // status-specific definition
            status_map_model[tools::str::toi(s)].push_back(a);
    }
}

/// Deal with the definition 'k' of augmented to be sampled.
/// NOTE: this assume definitions are univocal and status correct.
void MCMC::Impl::Augmented::build_helper_variable_to_sample_(size_t k)
{
    const std::string &a = deftable.at("name")[k];
    const std::string &s = deftable.at("status")[k];
    const std::string &value = deftable.at("value")[k];
    const std::string &rule = deftable.at("rule")[k];
    const bool init_rand = tools::str::tobool(deftable.at("init_rand")[k]);

    // check random initialization option
    // UPDATE: check pImpl_.setup.rand_init_params
    // UPDATE: check value
    if (init_rand)
    {
        if (tools::is_na(deftable.at("model")[k]))
            throw std::invalid_argument(
                "Unable to randomly initialize augmented\"" +
                a + "\" without a model.");
    }
    else if (tools::is_na(value))
        throw std::invalid_argument(
            "No initial value provided for augmented \"" + a + "\" " +
            "to be sampled.\nIn alternative, enable random initialization.");

    // read option for monte carlo proposer
    if (!tools::is_na(rule))
    {
        list_sample.push_back(a); // NOTE: make unique elsewhere

        if (tools::is_na(s)) // generic definition
            for (int status : pImpl_.dataset.all_infected_status_list)
                status_map_sample[status].push_back(a);
        else // status-specific definition
            status_map_sample[tools::str::toi(s)].push_back(a);
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////  MCMC : Impl : Augmented : sample draws  ////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/// Set individual id for augmented draws follow-up.
void MCMC::Impl::Augmented::build_helper_set_id_draws_(str_vec &id_draws)
{
    smpl_individuals.clear();

    if (!id_draws.empty())
    {
        if (is_any_sample())
        {
            smpl_individuals.reserve(id_draws.size());

            // check id across all individuals
            str_vec already_stored_id = {};
            for (const Household &house : pImpl_.dataset)
                for (const Individual &member : house)
                    if (id_draws % member.id())
                    {
                        if (already_stored_id % member.id())
                            throw std::invalid_argument(
                                "Multiple individuals share the same id trace \"" + member.id() + "\".");
                        else
                        {
                            already_stored_id.push_back(member.id());
                            smpl_individuals.emplace_back(&member);
                        }
                    }

            // check all user request are satisfied
            tools::check_for_missing_entries(already_stored_id, id_draws, "trace id");

            if (is_any_draw())
            {
                // initialize list with augmented headers for draws
                list_draws = tools::str::unique(tools::str::drop_private(list_sample + list_func));

                for (const auto member : smpl_individuals)
                    if (!member->is_infected_status())
                        tools::warning(
                            "Ignornig trace of non infected individual id \"" + member->id() + "\".");
            }
            else
            {
                // UPDATE: is this even possible?
                tools::warning(
                    "Ignoring id trace for individual(s) " +
                    tools::str::format(id_draws, HSP_MAX_STR_PRINT));
                id_draws.clear();
            }
        }
        else
        {
            tools::warning(
                "Ignoring id trace. No augmented variables to infer.");
            id_draws.clear();
        }
    }
}

/// Initialize sample for individual augmented variable followup.
void MCMC::Impl::Augmented::init_draws(unsigned long sample_size)
{
    if (is_any_draw())
    {
        smpl_values.clear();
        smpl_values.resize(smpl_individuals.size());
        for (size_t k = 0; k < smpl_individuals.size(); ++k)
        {
            smpl_values[k].clear();
            smpl_values[k].reserve(list_draws.size());
            for (const std::string &a : list_draws)
                smpl_values[k][a].assign(sample_size, HSP_NA_DOUBLE);
        }
    }
}

/// Store sample for individual augmented variable followup.
void MCMC::Impl::Augmented::update_draws(unsigned long smpl_count)
{
    if (is_any_draw())
        for (size_t k = 0; k < smpl_individuals.size(); ++k)
        {
            auto &smpl_augm_value_k = smpl_values[k];
            const auto member_star = smpl_individuals[k];
            for (const std::string &a : list_draws)
                smpl_augm_value_k[a][smpl_count] = member_star->augm(a);
        }
}

/// Clears all memory allocated by the Augmented structure.
/// Call this when the data is no longer needed to prevent memory leaks.
void MCMC::Impl::Augmented::clear_draws()
{
    smpl_individuals.clear();
    smpl_individuals.shrink_to_fit();

    for (auto &entry : smpl_values)
    {
        for (auto &pair : entry)
        {
            pair.second.clear();
            pair.second.shrink_to_fit();
        }
        entry.clear();
        entry.rehash(0);
    }

    smpl_values.clear();
    smpl_values.shrink_to_fit();
}

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////  MCMC::Impl : assign augmented variables  //////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/// Clear all augmented definition of individual.
void MCMC::Impl::Augmented::clear_all(Individual &member)
{
    for (const std::string &a : list_full)
        member.clear_augm(a);
}

/// Initialize augmented dependecies for all households.
void MCMC::Impl::Augmented::init()
{
    try
    {
        if (is_any())
            for (Household &house : pImpl_.dataset)
                for (Individual &member : house.infected())
                {
                    // intialize all augmented definitions
                    init_all(member);
                    // initialize augmented model
                    member.init_augm_model_logpdf(status_map_model);
                }
    }
    catch (const std::exception &e)
    {
        throw std::runtime_error(
            "While initializing augmented variables:\n" +
            std::string(e.what()));
    }
}

/// Assign augmented variables to single individual.
void MCMC::Impl::Augmented::init_all(Individual &member)
{
    // store the entries of deftable that are relevant for the individual
    // FIXME: need to be reordered according to dependency tree
    std::vector<size_t> indiv_augm_entries;

    // loop over augmented variable definitions to assign placheholders
    for (size_t k = 0; k < deftable.at("name").size(); ++k)
    {
        const std::string &s = deftable.at("status")[k];

        // assign augmented variable if undefined status or corresponding status
        if (tools::is_na(s) || (member.status() == tools::str::toi(s)))
        {
            // assign augmented variable placeholder
            member.init_augm(deftable.at("name")[k],
                             type_to_double_(member, deftable.at("value")[k], deftable.at("value_at")[k]),
                             tools::str::tobool(deftable.at("const")[k]),
                             tools::str::tobool(deftable.at("discr")[k]),
                             deftable.at("func")[k],
                             tools::str::tobool(deftable.at("fcumul")[k]),
                             deftable.at("model")[k]);

            // add entry to individual augmented variable list
            indiv_augm_entries.push_back(k);
        }
    }

    // assign dependencies to each augmented variable
    // FIXME: this should be ordered according to dependency tree
    for (size_t k : indiv_augm_entries)
        init_single_var_dependencies_(member, k);
    // assign function reshaping and check if inherently constants
    // UPDATE: no need for a second loop if previous follows dependency order
    for (size_t k : indiv_augm_entries)
        member.get_augm_variable(deftable.at("name")[k]).init_func_reshaping();
}

/// Assign model arguments and range. Same for function, otherwise assign proposer.
/// FIXME: if attribute is augmented to be initialized after, it might cause pointer errors.
void MCMC::Impl::Augmented::init_single_var_dependencies_(Individual &member, size_t k)
{
    const std::string a = deftable.at("name")[k];
    Variable &xaugm = member.get_augm_variable(a);

    // assign arguments to model
    if (!tools::is_na(deftable.at("arg_key")[k]))
        init_attribute_list_(HSP_VAR_MODEL, member, a,
                             tools::parse_keys(deftable.at("arg_key")[k], HSP_LIST_DELIM),
                             tools::parse_arg_list(deftable.at("arg_val")[k], HSP_LIST_DELIM),
                             tools::parse_arg_list(deftable.at("arg_at")[k], HSP_LIST_DELIM));

    // assign range to model
    if (!(tools::is_na(deftable.at("min")[k]) && tools::is_na(deftable.at("max")[k])))
        init_attribute_list_(HSP_VAR_MODEL, member, a,
                             HSP_RANGE_KEYS,
                             {deftable.at("min")[k], deftable.at("max")[k]},
                             {deftable.at("min_at")[k], deftable.at("max_at")[k]});
    xaugm.init_distr_reshaping();

    if (xaugm.is_function()) // function
    {
        // assign arguments to function
        if (!tools::is_na(deftable.at("farg_key")[k]))
            init_attribute_list_(HSP_VAR_FUNCTION, member, a,
                                 tools::parse_keys(deftable.at("farg_key")[k], HSP_LIST_DELIM),
                                 tools::parse_arg_list(deftable.at("farg_val")[k], HSP_LIST_DELIM),
                                 tools::parse_arg_list(deftable.at("farg_at")[k], HSP_LIST_DELIM));
        // assign range to function
        if (!(tools::is_na(deftable.at("min")[k]) && tools::is_na(deftable.at("max")[k])))
            init_attribute_list_(HSP_VAR_FUNCTION, member, a,
                                 HSP_RANGE_KEYS,
                                 {deftable.at("min")[k], deftable.at("max")[k]},
                                 {deftable.at("min_at")[k], deftable.at("max_at")[k]});
    }
    else if (xaugm.is_pointer()) // pointer (case of identity function)
    {
        const Attribute &aux_attr = read_single_attribute_(
            member, a,
            tools::parse_arg_list(deftable.at("farg_val")[k], HSP_LIST_DELIM)[0],
            tools::parse_arg_list(deftable.at("farg_at")[k], HSP_LIST_DELIM)[0]);

        if (aux_attr.is_scalar())
            xaugm.init_scalar(aux_attr.get());
        else
            xaugm.init_pointer(aux_attr.get_ptr());
    }
    else if (!xaugm.is_constant()) // variable to be sampled
    {
        // assign proposer rule
        xaugm.init_proposer(deftable.at("rule")[k], tools::str::tod(deftable.at("scale")[k]));

        // random initialization
        // UPDATE: need a more general switch connected to pImpl_.setup.rand_init_params)
        if (tools::str::tobool(deftable.at("init_rand")[k]))
            xaugm.init_rand(pImpl_.setup.prng());
    }
}

/////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////  MCMC::Impl : methods (type-defined variables)  /////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

/// Convert a typed string into a double.
double MCMC::Impl::Augmented::type_to_double_(const Individual &member,
                                              const std::string &value,
                                              const std::string &type)
{
    double output = HSP_NA_DOUBLE;

    if (!(tools::is_na(value) && tools::is_na(type)))
    {
        if (!tools::is_na(value) && tools::is_na(type))
            output = tools::str::tod(value, "argument without type");

        else if (tools::is_na(value) && !tools::is_na(type))
            tools::warning(
                "Ignoring value type without a paired value.");

        else if (tools::str::toup(type) == HSP_TYPE_SCALAR)
            output = tools::str::tod(value, "argument with scalar type");

        else if (tools::str::toup(type) == HSP_TYPE_DATA)
            output = member.data(value);

        else if (tools::str::toup(type) == HSP_TYPE_PARAM)
            if (pImpl_.param.at(value).is_constant())
                output = pImpl_.param.get(value);
            else
                throw std::invalid_argument(
                    "Parameter to be inferred \"" + value + "\" assigned as constant argument.");

        else if (tools::str::toup(type) == HSP_TYPE_AUGM)
            throw std::invalid_argument(
                "Value initialization with other augmented variable is not supported.");

        else
            throw std::invalid_argument(
                "Unknown variable type specifier \"" + type + "\".");
    }

    return output;
}

/// Auxiliary method to assign attributes to augmented function or model.
void MCMC::Impl::Augmented::init_attribute_list_(bool is_var_model,
                                                 Individual &member,
                                                 const std::string &a,
                                                 const str_vec &all_keys,
                                                 const str_vec &all_values,
                                                 const str_vec &all_types)
{
    Variable &xaugm = member.get_augm_variable(a);

    // check size validity of keys, values and type
    const size_t n_keys = all_keys.size();
    if (n_keys > 0)
        if (all_values.size() != n_keys || all_types.size() != n_keys)
            throw std::invalid_argument(
                "Mismatching attributes in definition of " + xaugm.err_name() + ".");

    // loop over all attributes
    for (size_t i = 0; i < all_keys.size(); ++i)
    {
        const std::string key = all_keys[i], value = all_values[i], type = all_types[i];

        // get auxiliary attribute
        const Attribute &aux_attr = read_single_attribute_(member, a, value, type);

        // cache dependecy if an augmented variable depends on a parameter to be inferred
        if (is_var_model == HSP_VAR_FUNCTION)
        {
            if (tools::str::toup(type) == HSP_TYPE_PARAM && aux_attr.is_pointer())
                pImpl_.param.make_cache_dependent(value);

            xaugm.init_func_attr(key, aux_attr);
        }
        else // is_var_model == HSP_VAR_MODEL
            xaugm.init_distr_attr(key, aux_attr);
    }
}

/// Auxiliary method to read typed member values as a single attribute.
Attribute MCMC::Impl::Augmented::read_single_attribute_(Individual &member,
                                                        const std::string &a,
                                                        const std::string &val,
                                                        const std::string &type)
{
    Attribute output;

    if (tools::is_na(val) && tools::is_na(type))
        ; // pass

    else if (!tools::is_na(val) && tools::is_na(type))
        throw std::invalid_argument(
            "Unable to interpret default value of  \"" + val + "\". Please provide a type.");

    else if (tools::is_na(val) && !tools::is_na(type))
        throw std::invalid_argument(
            "Found a type without a paired value. Please correct.");

    else if (tools::str::toup(type) == HSP_TYPE_SCALAR)
        output.init_scalar(tools::str::tod(val, "argument with scalar type"));

    else if (tools::str::toup(type) == HSP_TYPE_DATA)
        output.init_scalar(member.data(val));

    else if (tools::str::toup(type) == HSP_TYPE_PARAM)
    {
        if (pImpl_.param.at(val).is_constant())
            output.init_scalar(pImpl_.param.get(val));
        else
            output.init_pointer(pImpl_.param.get_ptr(val));
    }
    else if (tools::str::toup(type) == HSP_TYPE_AUGM)
    {
        if (!(list_full % val))
            throw std::invalid_argument(
                "Undeclared augmented variable \"" + val + "\".");

        Variable &other_augm = member.get_augm_variable(val);

        // assert_arg no nested definitions
        if (other_augm.is_chain_dependency(member.get_augm_variable(a)))
            throw std::invalid_argument(
                "Found nested definitions for \"" + a + "\" and \"" + val + "\". Please correct.");

        if (other_augm.is_function())
            output.init_function(other_augm.get_function());
        else
            output.init_pointer(other_augm.get_ptr());
    }
    else
        throw std::invalid_argument(
            "Unsopported type \"" + type + "\" found while loading \"" + val + "\" for pointer value.");

    return output;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////  MCMC : move all augmented variables  /////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/// Move augmented variables from all households.
void MCMC::Impl::move_all_augm_()
{
    if (augm.is_any_sample())
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
                        error_message = "While sampling augmented variables for household \"" +
                                        this_house.house_id() + "\":\n" +
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

/// Move augmented variables of infectious members inside household.
/// The move is performed one member at a time, one augmented variable at a time.
/// UPDATE: smooth the incrementation by reference by just adding...
double MCMC::Impl::move_augm_single_house_(Household &house)
{
    double delta_ll = 0.;
    // Loop over infectious residents of the househould
    for (Individual &member_i : house.infected())
        // Loop over augmented variables of infectious
        for (const std::string &a : augm.status_map_sample[member_i.status()])
            delta_ll += metropolis_augm_single_move_(house, member_i, a);
    return delta_ll;
}

/// Take care of a a single augmented Metropolis move.
double MCMC::Impl::metropolis_augm_single_move_(Household &house, Individual &member, const std::string &a)
{
    Variable &xaugm = member.get_augm_variable(a);
    xaugm.propose(setup.prng()); // propose new value and buffer old

    // proceed only if the new value is consistent
    if (xaugm.is_valid_move())
    {
        // update dependencies to perform computations
        house.buffer_all_cache();
        update_cache_dependecies_(house, member);

        // Metropolis-Hastings algorithm
        double D_ll = -(house.get_obs_loglike() + member.get_augm_model_loglike());
        // NOTE: because of dependency between augmented variables
        // I am recomputing everything // UPDATE: to be optimized
        D_ll += member.compute_and_set_augm_model_loglike();
        D_ll += compute_and_set_household_obs_loglike_(house);

        if (setup.unif() < std::exp(setup.inv_T * D_ll - xaugm.get_propose_logratio()))
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
