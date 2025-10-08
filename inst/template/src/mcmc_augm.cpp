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

#include "tools/string.h"
#include "tools/read.h"
#include "tools/error.h"
#include "tools/status.h"
#include "tools/params.h"
#include "tools/sources.h"
#include "mcmc_impl.h"

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////  dependency tree helpers  //////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/// Helper to build augmented definition order and assert for acyclic.
void walk_dependency_tree(const std::string &entry,
                          const u_dict<str_vec> &tree,
                          str_vec &visited,
                          str_vec &recursion_stack,
                          str_vec &ordered)
{
    if (!(tree % entry))
        throw std::invalid_argument(
            "Missing definition for augmented \"" + entry + "\"");
    if (recursion_stack % entry)
        throw std::invalid_argument(
            "Inadmissible nested definitions for augmented \"" + entry + "\"");
    if (visited % entry)
        return;

    recursion_stack.push_back(entry);
    for (auto &dep : tree.at(entry))
        walk_dependency_tree(dep, tree, visited, recursion_stack, ordered);

    recursion_stack = recursion_stack - str_vec({entry});
    visited.push_back(entry);
    ordered.push_back(entry);
}

/// Get topological order of definitions dependencies.
str_vec extract_all_dependencies_order(const u_dict<str_vec> &graph)
{
    // inspect dependency tree for cycles and missing entries
    str_vec visited, recursion_stack, ordered;
    for (const auto &entry : graph)
        if (!(visited % entry.first))
            walk_dependency_tree(entry.first, graph, visited, recursion_stack, ordered);
    return ordered;
}

/// Get topological order of definitions depending on 'start' (start included).
str_vec extract_dependents_order(const std::string &start, const u_dict<str_vec> &graph)
{
    // construct a reversed graph for dependents
    u_dict<str_vec> reverse_graph;
    for (auto &kv : graph)
    {
        const std::string &node = kv.first;
        for (const std::string &dep : kv.second)
            reverse_graph[dep].push_back(node);

        // ensure isolated nodes are present
        if (!(reverse_graph % node))
            reverse_graph[node] = {};
    }

    // walk dependents of start
    str_vec visited, recursion_stack, ordered;
    for (const std::string &dep : reverse_graph.at(start))
        walk_dependency_tree(dep, reverse_graph, visited, recursion_stack, ordered);

    std::reverse(ordered.begin(), ordered.end());
    return str_vec({start}) + ordered;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////  MCMC : Impl : definition helpers  //////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/// Auxiliary method to read typed member values as a single attribute.
Attribute MCMC::Impl::read_attribute_(const DefinitionParser::ASTNode &node,
                                      const Individual &member)
{
    if (node.is_empty()) // empty node
        return Attribute();

    else if (!node.is_terminal()) // function node
    {
        Attribute attr(node.func);
        for (const auto &child : node.children)
            attr.get_function()->init_attr(child.key, read_attribute_(child, member));
        return attr;
    }

    else if (hsp::sources::is_augmented(node.type)) // terminal augmented node
    {
        Attribute attr;
        attr.name(node.value, "augmented");
        const Variable &other_augm = member.get_augm_variable(node.value);
        attr.transfer(other_augm);
        return attr;
    }

    else if (hsp::sources::is_data(node.type)) // terminal data node
    {
        Attribute attr(member.data(node.value));
        attr.name(node.value, "data");
        return attr;
    }

    else if (hsp::sources::is_numeric(node.type)) // terminal numeric node
        return Attribute(hsp::str::tod(node.value, "argument with numeric type"));

    else if (hsp::sources::is_parameter(node.type)) // terminal parameter node
    {
        Attribute attr;
        attr.name(node.value, "parameter");
        if (param.at(node.value).is_constant())
            attr.init_numeric(param.get(node.value));
        else
            attr.init_pointer(param.get_ptr(node.value));
        return attr;
    }

    else
        throw std::invalid_argument(
            "Unknown source type \"" + node.type + "\" found while parsing \"" + node.str() + "\".");
}
Attribute MCMC::Impl::read_attribute_(const std::string &expr, const Individual &member)
{
    return read_attribute_(DefinitionParser(expr).root, member);
}

/// Convert a typed string into a double.
double MCMC::Impl::read_numeric_(const DefinitionParser::ASTNode &node, const Individual &member)
{
    const Attribute aux_attr_ = read_attribute_(node, member);

    if (aux_attr_.is_empty())
        return hsp::na::DOUBLE;
    else if (aux_attr_.is_constant())
        return aux_attr_.get();
    else
        throw std::invalid_argument(
            "Unable to interpret definition \"" + node.str() + "\" as numeric constant.");
}
double MCMC::Impl::read_numeric_(const std::string &expr, const Individual &member)
{
    return read_numeric_(DefinitionParser(expr).root, member);
}
void MCMC::Impl::make_param_cache_depen_(const DefinitionParser::ASTNode &node)
{
    // cache dependecy if an augmented function depends on a parameter to be inferred
    if (hsp::sources::is_parameter(node.type) && !param.at(node.value).is_constant())
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
    // eventually append infectivity profile definition
    pImpl_.infec_profile.append_augm_definition();

    if (!deftable.empty())
    {
        const str_vec file_headers = {"name", "status",
                                      "func", "fcumul",
                                      "value", "const", "discr",
                                      "model", hsp::params::MIN, hsp::params::MAX,
                                      "rule", "scale"};

        // check entries
        str_vec headers = hsp::read::table_headers(deftable);
        hsp::str::assert_no_missing_entries(headers, file_headers);
        hsp::str::assert_no_unused_entries(headers, file_headers);

        helper_build_definition_();
        helper_build_definition_order_();
        helper_build_model_status_maps_();

        // allocate memory for augmented variables to each individual
        for (Household &house : pImpl_.dataset)
            for (Individual &member : house)
                clear_all(member);
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////  MCMC : Impl : Augmented : assert univocal definitions  ///////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/// Assert that augmented definitions are univocal, complete and acyclic.
void MCMC::Impl::Augmented::helper_build_definition_()
{
    // drop previous definitions
    status_to_def_row.clear();

    for (size_t k = 0; k < deftable["name"].size(); ++k)
    {
        // check augmented name
        const std::string &a = deftable["name"][k];
        if (hsp::na::is(a))
            throw std::invalid_argument(
                "Augmented variable name cannot be undefined.");

        // check that status-specific definition is for infected status
        std::string &s = deftable["status"][k];
        if (hsp::na::is(s))
            s = hsp::na::STRING; // enforce internal NA convention

        else if (!hsp::status::is_infected(hsp::str::toi(s, "augmented status")))
            throw std::invalid_argument(
                "Unable to associate augmented \"" + a + "\" for non infectious status " + s + ".");

        // store status to definition row index
        status_to_def_row[s].push_back(k);
    }

    // assert that name and status are univocal
    for (const auto &pair : status_to_def_row)
    {
        str_vec name_to_status = {};
        for (size_t k : pair.second)
            name_to_status.push_back(deftable["name"][k]);

        if (name_to_status.size() != hsp::str::unique(name_to_status).size())
            throw std::invalid_argument(
                "Incompatible definitions for augmented \"" + pair.first + "\".");
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////  MCMC : Impl : Augmented : build variable to sample  ////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/// Deal with the definition 'k' of augmented to be sampled (const = FALSE).
/// NOTE: this assume definitions are univocal and status correct.
void MCMC::Impl::Augmented::helper_build_definition_order_()
{
    // drop previous definitions
    status_map_sample.clear();
    list_full.clear();
    list_export.clear();
    list_draws.clear();

    // NOTE: NA (i. e. generic) does not necessarily need checking.
    if (pImpl_.dataset.all_infected_status_list.empty())
        set_definition_dependency_order_for_(hsp::na::STRING);
    else
        for (int status : pImpl_.dataset.all_infected_status_list)
            set_definition_dependency_order_for_(std::to_string(status));
}

/// Set dependecies ordering for ininitializiation fo given status.
void MCMC::Impl::Augmented::set_definition_dependency_order_for_(const std::string &s)
{
    u_dict<str_vec> tree;       // populate augmented dependency tree
    u_dict<size_t> name_to_row; // ... and reverse dictionary name-to-row-index
    for (size_t k : get_unsorted_definition_rows_by_(s))
    {
        str_vec xdepen = {};
        for (const std::string &c : str_vec({"value", "model", "func", hsp::params::MIN, hsp::params::MAX}))
            xdepen += DefinitionParser(deftable[c][k]).list_augmented_depen();
        tree[deftable["name"][k]] = hsp::str::unique(xdepen);
        name_to_row[deftable["name"][k]] = k;
    }

    // check dependecies and translate order from name to deftable row index
    try
    {
        status_to_def_row[s] = {}; // NOTE:recycle to provide order
        for (const std::string &a : extract_all_dependencies_order(tree))
            status_to_def_row[s].push_back(name_to_row[a]);
    }
    catch (const std::exception &e)
    {
        throw std::invalid_argument(std::string(e.what()) + " for status " + s + ".");
    }

    // define usage lists respecting definition order
    for (size_t k : status_to_def_row[s])
        set_sample_map_for_(s, k);

    // finalize augmented lists
    list_full = hsp::str::unique(list_full + list_export + list_draws); // NOTE: must be first!
    list_export = hsp::str::drop_private(hsp::str::unique(list_export + list_draws));
    list_draws = hsp::str::drop_private(hsp::str::unique(list_draws));
}

/// Extract status-specific definitions row indexes.
std::vector<size_t> MCMC::Impl::Augmented::get_unsorted_definition_rows_by_(const std::string &s) const
{
    // load the list of status-specific definition rows
    if (status_to_def_row % s)
    {
        std::vector<size_t> allrows = status_to_def_row.at(s);

        // extract names
        str_vec name_to_status = {};
        for (size_t k : allrows)
            name_to_status.push_back(deftable.find("name")->second[k]);

        // add generic definitions to current list (if not generic already)
        if (!hsp::na::is(s))
            if (status_to_def_row % hsp::na::STRING)
                for (size_t k : status_to_def_row.at(hsp::na::STRING))
                    // skip generic if status-specific already exist
                    if (!(name_to_status % deftable.find("name")->second[k]))
                        allrows.push_back(k);

        return allrows;
    }

    else if (status_to_def_row % hsp::na::STRING)
        // use generic definitions if status-specific don't exist
        return status_to_def_row.at(hsp::na::STRING);

    else
        throw std::invalid_argument(
            "No definitions compatible with augmented status " + s + ".");
}

/// Set sample map and relative usage list.
void MCMC::Impl::Augmented::set_sample_map_for_(const std::string &s, size_t k)
{
    const std::string &a = deftable["name"][k];

    if (hsp::str::tobool(deftable["const"][k]))
        list_export.push_back(a);

    else if (!hsp::na::is(deftable["func"][k]))
        set_function_or_identity_for_(k);

    else if (hsp::na::is(deftable["value"][k]) && hsp::na::is(deftable["model"][k]))
        throw std::invalid_argument(
            "Unable to initialize augmented\"" + a + "\" for status " + s + " without value nor model");

    else if (!hsp::na::is(deftable["rule"][k]))
    {
        // read option for monte carlo proposer
        list_draws.push_back(a);
        if (!hsp::na::is(s))
            status_map_sample[hsp::str::toi(s)].push_back(a);
    }
    else
        throw std::runtime_error(
            "Ambigous non-constant definition of augmented \"" + a +
            "\" for status " + s + ". No rule was provided.");
}

///////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////  MCMC : Impl : Augmented : build functions  /////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/// Build augmented function and inspect for constant.
void MCMC::Impl::Augmented::set_function_or_identity_for_(size_t k)
{
    const std::string &a = deftable["name"][k];
    const DefinitionParser parser(deftable["func"][k]);

    if (parser.is_terminal())
    {
        bool is_identity_const = true;

        // check the type of requested equality
        if (hsp::sources::is_parameter(parser.root.type))
            is_identity_const = pImpl_.param.at(parser.root.value).is_constant();

        else if (hsp::sources::is_augmented(parser.root.type))
            is_identity_const = false; // FIXME: use an augmented dependecies tree to see if constant

        if (is_identity_const)
            list_export.push_back(a);
        else
            list_draws.push_back(a);
    }
    else if (parser.is_implicit_callable())
        list_draws.push_back(a);
    else
        list_full.push_back(a);
}

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////  MCMC : Impl : Augmented : build variables models  /////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/// Deal with the definition of augmented variable with a model.
/// NOTE: here there's no need to initialize to generic since only used for posterior sampling.
void MCMC::Impl::Augmented::helper_build_model_status_maps_()
{
    // drop previous definitions
    status_map_model.clear();

    for (int status : pImpl_.dataset.all_infected_status_list)
    {
        u_dict<str_vec> tree;       // populate augmented dependency tree
        u_dict<size_t> name_to_row; // ... and reverse dictionary name to index
        str_vec without_model;      // list of augmented without models
        for (const size_t k : get_unsorted_definition_rows_by_(std::to_string(status)))
        {
            const std::string &a = deftable["name"][k];
            str_vec xdepen = {};
            for (const std::string &c : str_vec({"value", "model", "func", hsp::params::MIN, hsp::params::MAX}))
                xdepen += DefinitionParser(deftable[c][k]).list_augmented_depen();

            if (hsp::na::is(deftable["model"][k]))
                without_model.push_back(a);
            else
                tree[a] = hsp::str::unique(xdepen);

            name_to_row[a] = k;
        }

        // prune dependency tree from augmented without model
        for (auto &pair : tree)
            pair.second = pair.second - without_model;

        // assign order of dependents for augmented model computation
        for (const auto &pair : tree)
            status_map_model[status].emplace(pair.first, extract_dependents_order(pair.first, tree));
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////  MCMC : Impl : Augmented : assign augmented variables  ////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/// Clear all augmented definition of individual.
void MCMC::Impl::Augmented::clear_all(Individual &member)
{
    for (const std::string &a : list_full)
        member.clear_augm(a);
}

/// Initialize augmented dependecies for evaluations.
/// NOTE: this must not be called before simulation which build on demand.
void MCMC::Impl::Augmented::init()
{
    if (is_any())
        for (Household &house : pImpl_.dataset)
            for (Individual &member : house.infected())
                try
                {
                    // intialize all augmented definitions
                    build_all(member);
                    // initialize augmented model
                    member.init_augm_model_logratios(status_map_model[member.status()]);
                }
                catch (const std::exception &e)
                {
                    throw std::runtime_error(
                        "While initializing " + section_name_ +
                        " for " + member.err_name() + ":\n" +
                        std::string(e.what()));
                }

    helper_init_id_draws_();
}

/// Assign all augmented variables to single individual.
void MCMC::Impl::Augmented::build_all(Individual &member)
{
    // loop over augmented variable definitions to assign placheholders
    for (size_t k : init_get_definition_order_(member))
        helper_init_attributes_(member, k);
}

/// Assign augmented variables only if missing (used for simulations).
void MCMC::Impl::Augmented::build_missing(Individual &member)
{
    // loop over augmented variable definitions to assign placheholders
    for (size_t k : init_get_definition_order_(member))
        if (member.get_augm_variable(deftable["name"][k]).is_empty())
            helper_init_attributes_(member, k);
}

/// Get the definition order for current member according to its status.
std::vector<size_t> MCMC::Impl::Augmented::init_get_definition_order_(Individual &member)
{
    const std::string s = std::to_string(member.status());

    // NOTE: if status hasn't been observed in dataset, build on demand
    if (!(pImpl_.dataset.all_infected_status_list % member.status()))
        if (member.is_infected_status())
        {
            set_definition_dependency_order_for_(s);
            pImpl_.dataset.all_infected_status_list.push_back(member.status());
        }

    if (status_to_def_row % s)
        return status_to_def_row[s];
    else if (status_to_def_row % hsp::na::STRING)
        return status_to_def_row[hsp::na::STRING];
    else if (member.is_infected_status())
        throw std::invalid_argument(
            "No augmented definitions for " + member.err_name() + " with status " + s + ".");
    else
        throw std::invalid_argument(
            "Missing generic augmented definition for non-infectious " + member.err_name() + ".");
}

/// Assign model arguments and range. Same for function, otherwise assign proposer.
void MCMC::Impl::Augmented::helper_init_attributes_(Individual &member, size_t k)
{
    const std::string &a = deftable["name"][k];
    member.init_augm(a,
                     pImpl_.read_numeric_(deftable["value"][k], member),
                     hsp::str::tobool(deftable["const"][k]),
                     hsp::str::tobool(deftable["discr"][k]),
                     deftable["func"][k],
                     hsp::str::tobool(deftable["fcumul"][k]),
                     deftable["model"][k]);
    Variable &xaugm = member.get_augm_variable(a);

    if (xaugm.is_function()) // function
        init_func_attributes_(member, k);

    else if (xaugm.is_pointer()) // pointer (case of identity function)
        xaugm.transfer(pImpl_.read_attribute_(deftable["func"][k], member));

    if (!xaugm.is_constant()) // variable to be sampled
        xaugm.init_proposer(deftable["rule"][k], hsp::str::tod(deftable["scale"][k]));

    // assign arguments to model
    init_model_attributes_(member, k);
}

/// Auxiliary method to assign attributes to augmented model.
void MCMC::Impl::Augmented::init_model_attributes_(Individual &member, size_t k)
{
    Variable &xaugm = member.get_augm_variable(deftable["name"][k]);

    for (const auto &child : DefinitionParser(deftable["model"][k]).root.children)
        xaugm.init_distr_attr(child.key, pImpl_.read_attribute_(child, member));

    for (const std::string &key : {hsp::params::MIN, hsp::params::MAX})
        xaugm.init_distr_attr(key, pImpl_.read_attribute_(deftable[key][k], member));
    xaugm.init_distr_reshaping();
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

    for (const std::string &key : {hsp::params::MIN, hsp::params::MAX})
        xaugm.init_func_attr(key, pImpl_.read_attribute_(deftable[key][k], member));
    xaugm.init_func_reshaping();
}

/////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////  MCMC : Impl : Augmented : user reset  /////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

/// Reset augmented initial values for given individual.
void MCMC::Impl::Augmented::user_reset(Individual &member)
{
    for (const std::string &a : status_map_sample[member.status()])
        member.get_augm_variable(a).user_reset(pImpl_.setup.prng());
}

/// Initialize augmented initial values for all individuals.
void MCMC::Impl::Augmented::user_reset()
{
    if (is_any())
        for (Household &house : pImpl_.dataset)
            for (Individual &member : house.infected())
                user_reset(member);
}

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////  MCMC : Impl : Augmented : sample draws  ////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/// Set individual id for augmented draws follow-up.
void MCMC::Impl::Augmented::helper_init_id_draws_()
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
            hsp::str::assert_no_missing_entries(already_stored_id, id_draws, "trace id");

            if (is_any_draw())
            {
                for (const auto member : smpl_individuals)
                    if (!member->is_infected_status())
                        hsp::warning(
                            "Ignornig draws of non infected " + member->err_name() + ".");
            }
            else
            {
                // UPDATE: is this even possible?
                hsp::warning(
                    "Ignoring id draws for individual(s) " + hsp::str::format(id_draws));
                id_draws.clear();
            }
        }
        else
        {
            hsp::warning(
                "Ignoring id draws. No augmented variables to infer.");
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
                smpl_values[k][a].assign(sample_size, hsp::na::DOUBLE);
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
        // deep clean pre-existing draws
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
double MCMC::Impl::move_augm_single_house_(Household &house)
{
    double delta_ll = 0.;
    // loop over infectious residents of the househould
    for (Individual &member_i : house.infected())
        // loop over augmented variables of infectious
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
        double delta_ll = -house.get_obs_loglike();
        delta_ll += member.compute_and_set_augm_model_logratio(a);
        delta_ll += compute_and_set_household_obs_loglike_(house);

        if (setup.unif() < std::exp(setup.inv_T * delta_ll - xaugm.get_propose_logratio()))
            return delta_ll;
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
