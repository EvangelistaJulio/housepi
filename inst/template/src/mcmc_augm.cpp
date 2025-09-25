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
/////////////////////////////////  MCMC : Impl : definition helpers  //////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/// Auxiliary method to read typed member values as a single attribute.
Attribute MCMC::Impl::read_attribute_(const DefinitionParser::ASTNode &node,
                                      const Individual &member)
{
    if (node.is_empty())
        return Attribute();

    else if (!node.is_terminal())
    {
        Attribute attr(node.func);
        for (const auto &child : node.children)
            attr.get_function()->init_attr(child.key, read_attribute_(child, member));
        return attr;
    }

    else if (HSP_IS_TYPE_SCALAR(node.type))
        return Attribute(tools::str::tod(node.value, "argument with scalar type"));

    else if (HSP_IS_TYPE_DATA(node.type))
        return Attribute(member.data(node.value));

    else if (HSP_IS_TYPE_PARAM(node.type))
    {
        if (param.at(node.value).is_constant())
            return Attribute(param.get(node.value));
        else
            return Attribute(param.get_ptr(node.value));
    }
    else if (HSP_IS_TYPE_AUGM(node.type))
    {
        if (!(augm.list_full % node.value))
            throw std::invalid_argument(
                "Undeclared augmented variable \"" + node.value + "\".");

        const Variable &other_augm = member.get_augm_variable(node.value);
        if (other_augm.is_function())
            return Attribute(other_augm.get_function());
        else
            return Attribute(other_augm.get_ptr());
    }
    else
        throw std::invalid_argument(
            "Unknown error found while loading \"" + node.str() + "\".");
}
Attribute MCMC::Impl::read_attribute_(const std::string &expr, const Individual &member)
{
    return read_attribute_(DefinitionParser(expr).root, member);
}

/// Convert a typed string into a double.
double MCMC::Impl::read_scalar_(const DefinitionParser::ASTNode &node, const Individual &member)
{
    const Attribute aux_attr_ = read_attribute_(node, member);

    if (aux_attr_.is_empty())
        return HSP_NA_DOUBLE;

    else if (aux_attr_.is_constant())
        return aux_attr_.get();

    else
        throw std::invalid_argument(
            "Unable to interpret definition \"" + node.str() + "\" as a double.");
}
double MCMC::Impl::read_scalar_(const std::string &expr, const Individual &member)
{
    return read_scalar_(DefinitionParser(expr).root, member);
}
void MCMC::Impl::make_param_cache_depen_(const DefinitionParser::ASTNode &node)
{
    // cache dependecy if an augmented function depends on a parameter to be inferred
    if (HSP_IS_TYPE_PARAM(node.type) && !param.at(node.value).is_constant())
        param.make_cache_dependent(node.value);
}

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
            build_helper_set_id_draws_();

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
    for (size_t k = 0; k < deftable["name"].size(); ++k)
    {
        const std::string &a = deftable["name"][k];
        const std::string &s = deftable["status"][k]; // NOTE: can be NA, i. e. generic

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

    // UPDATE
    // Inspect dependency tree, assert correctness and provide ordering

    for (size_t k = 0; k < deftable["name"].size(); ++k)
    {
        build_helper_variable_with_model_(k);

        if (!tools::is_na(deftable["func"][k]))
            build_function_(k);
        else if (tools::str::tobool(deftable["const"][k]))
            list_const.push_back(deftable["name"][k]);
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
    if (DefinitionParser(deftable["func"][k]).is_terminal())
        build_function_identity_helper_(k);
    else
        list_func.push_back(deftable["name"][k]);
}

/// Build augmented equality function and change deftable.
void MCMC::Impl::Augmented::build_function_identity_helper_(size_t k)
{
    const std::string &a = deftable["name"][k];
    const auto parser = DefinitionParser(deftable["func"][k]);

    bool is_identity_const = true;

    // check the type of requested equality
    if (HSP_IS_TYPE_PARAM(parser.root.type))
        is_identity_const = pImpl_.param.at(parser.root.value).is_constant();
    else if (HSP_IS_TYPE_AUGM(parser.root.type))
        is_identity_const = false; // FIXME: use an augmented dependecies tree to see if constant

    if (is_identity_const)
        list_const.push_back(a);
    else
        list_func.push_back(a);
}

///////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////  MCMC : Impl : Augmented : build variables  /////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/// Deal with the definition of augmented variable with a model.
void MCMC::Impl::Augmented::build_helper_variable_with_model_(size_t k)
{
    if (!tools::is_na(deftable["model"][k]))
    {
        const std::string &a = deftable["name"][k];
        const std::string &s = deftable["status"][k];

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
    const std::string &a = deftable["name"][k];
    const std::string &s = deftable["status"][k];
    const std::string &value = deftable["value"][k];
    const std::string &rule = deftable["rule"][k];
    const bool init_rand = tools::str::tobool(deftable["init_rand"][k]);

    // check random initialization option
    // UPDATE: check pImpl_.setup.rand_init_params
    // UPDATE: check value
    if (init_rand)
    {
        if (tools::is_na(deftable["model"][k]))
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
void MCMC::Impl::Augmented::build_helper_set_id_draws_()
{
    str_vec &id_draws = pImpl_.setup.id_draws;
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
                            "Ignornig trace of non infected " + member->err_name() + ".");
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
    for (size_t k = 0; k < deftable["name"].size(); ++k)
    {
        const std::string &s = deftable["status"][k];

        // assign augmented variable if undefined status or corresponding status
        if (tools::is_na(s) || (member.status() == tools::str::toi(s)))
        {
            // assign augmented variable placeholder
            member.init_augm(deftable["name"][k],
                             pImpl_.read_scalar_(deftable["value"][k], member),
                             tools::str::tobool(deftable["const"][k]),
                             tools::str::tobool(deftable["discr"][k]),
                             deftable["func"][k],
                             tools::str::tobool(deftable["fcumul"][k]),
                             deftable["model"][k]);

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
        member.get_augm_variable(deftable["name"][k]).init_func_reshaping();
}

/// Assign model arguments and range. Same for function, otherwise assign proposer.
/// FIXME: if attribute is augmented to be initialized after, it might cause pointer errors.
void MCMC::Impl::Augmented::init_single_var_dependencies_(Individual &member, size_t k)
{
    Variable &xaugm = member.get_augm_variable(deftable["name"][k]);

    // assign arguments to model
    init_model_attributes_(member, k);
    xaugm.init_distr_attr(HSP_PARAM_MIN, pImpl_.read_attribute_(deftable[HSP_PARAM_MIN][k], member));
    xaugm.init_distr_attr(HSP_PARAM_MAX, pImpl_.read_attribute_(deftable[HSP_PARAM_MAX][k], member));
    xaugm.init_distr_reshaping();

    if (xaugm.is_function()) // function
    {
        // assign arguments to function
        init_func_attributes_(member, k);
        xaugm.init_func_attr(HSP_PARAM_MIN, pImpl_.read_attribute_(deftable[HSP_PARAM_MIN][k], member));
        xaugm.init_func_attr(HSP_PARAM_MAX, pImpl_.read_attribute_(deftable[HSP_PARAM_MAX][k], member));
        xaugm.init_func_reshaping();
    }
    else if (xaugm.is_pointer()) // pointer (case of identity function)
    {
        // UPDATE: I think this can be collapsed with above at a given point
        const Attribute &aux_attr = pImpl_.read_attribute_(deftable["func"][k], member);

        if (aux_attr.is_scalar())
            xaugm.init_scalar(aux_attr.get());
        else
            xaugm.init_pointer(aux_attr.get_ptr());
    }
    else if (!xaugm.is_constant()) // variable to be sampled
    {
        // assign proposer rule
        xaugm.init_proposer(deftable["rule"][k], tools::str::tod(deftable["scale"][k]));

        // random initialization
        // UPDATE: need a more general switch connected to pImpl_.setup.rand_init_params)
        if (tools::str::tobool(deftable["init_rand"][k]))
            xaugm.init_rand(pImpl_.setup.prng());
    }
}

/////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////  MCMC::Impl : methods (type-defined variables)  /////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

/// Auxiliary method to assign attributes to augmented model.
void MCMC::Impl::Augmented::init_model_attributes_(Individual &member, size_t k)
{
    Variable &xaugm = member.get_augm_variable(deftable["name"][k]);

    for (const auto &child : DefinitionParser(deftable["model"][k]).root.children)
        xaugm.init_distr_attr(child.key, pImpl_.read_attribute_(child, member));
}

/// Auxiliary method to assign attributes to augmented function.
void MCMC::Impl::Augmented::init_func_attributes_(Individual &member, size_t k)
{
    Variable &xaugm = member.get_augm_variable(deftable["name"][k]);

    for (const auto &child : DefinitionParser(deftable["func"][k]).root.children)
    {
        // get auxiliary attribute
        xaugm.init_func_attr(child.key, pImpl_.read_attribute_(child, member));
        pImpl_.make_param_cache_depen_(child);
    }
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
            if (has_error.load(std::memory_order_relaxed))
                continue;

            Household &this_house = dataset(H);
            try
            {
                delta_ll += move_augm_single_house_(this_house);
            }
            catch (const std::exception &e)
            {
#pragma omp critical
                {
                    if (!has_error.exchange(true))
                        error_message = "While sampling augmented variables for " +
                                        this_house.err_name() + ":\n" +
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
