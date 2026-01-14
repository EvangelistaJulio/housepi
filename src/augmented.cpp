/* housepi - Epidemiological Analysis of Household Outbreaks.
 *
 * Copyright © 2026 Institut Pasteur (Paris), CNRS and INSERM.
 * Distributed under the terms of the GNU General Public License (GPLv3).
 * See the LICENSE file for details.
 */

#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

#include "tools/global.h"
#include "tools/messages.h"
#include "tools/nan.h"
#include "tools/operators.h"
#include "tools/parser.h"
#include "tools/random.h"
#include "tools/read.h"
#include "tools/string.h"
#include "tools/types.h"

#include "variable.h"
#include "individual.h"
#include "household.h"
#include "data.h"
#include "parameter.h"
#include "infec_profile.h"
#include "augmented.h"

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////  dependency tree helpers  //////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

namespace
{
    using std::string;
    using std::vector;

    /// Helper to build augmented definition order and assert for acyclic.
    void walk_dependency_tree(const string &entry,
                              const hsp::str_table &tree,
                              vector<string> &visited,
                              vector<string> &recursion_stack,
                              vector<string> &ordered)
    {
        if (!tree.count(entry))
            throw std::invalid_argument(
                "Could not find a definition for augmented '" + entry + "'.");
        if (recursion_stack % entry)
            throw std::invalid_argument(
                "Found inadmissible cyclic definitions for augmented '" + entry + "'");
        if (visited % entry)
            return;

        recursion_stack.push_back(entry);
        for (auto &dep : tree.at(entry))
            walk_dependency_tree(dep, tree, visited, recursion_stack, ordered);

        recursion_stack = recursion_stack - vector<string>({entry});
        visited.push_back(entry);
        ordered.push_back(entry);
    }

    /// Get topological order of definitions dependencies.
    vector<string> extract_all_dependencies_order(const hsp::str_table &graph)
    {
        // inspect dependency tree for cycles and missing entries
        vector<string> visited, recursion_stack, ordered;
        for (const auto &[a, vec] : graph)
            if (!(visited % a))
                walk_dependency_tree(a, graph, visited, recursion_stack, ordered);
        return ordered;
    }

    /// Get topological order of definitions depending on 'start' (start included).
    vector<string> extract_dependents_order(const string &start, const hsp::str_table &graph)
    {
        // construct a reversed graph for dependents
        hsp::str_table reverse_graph;
        for (const auto &[node, dependents] : graph)
        {
            for (const auto &dep : dependents)
                reverse_graph[dep].push_back(node);

            // ensure isolated nodes are present
            if (!reverse_graph.count(node))
                reverse_graph[node] = {};
        }

        // walk dependents of start
        vector<string> visited, recursion_stack, ordered;
        for (const auto &dep : reverse_graph.at(start))
            walk_dependency_tree(dep, reverse_graph, visited, recursion_stack, ordered);

        std::reverse(ordered.begin(), ordered.end());
        return vector<string>({start}) + ordered;
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////  MCMC : Impl : Augmented : build  ///////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

Augmented::Augmented(std::shared_ptr<DataWide> dataptr,
                     std::shared_ptr<Parameter> paramptr,
                     std::shared_ptr<InfecProfile> infprofptr)
    : pDataWide_(dataptr),
      pParameter_(paramptr),
      pInfecProfile_(infprofptr) {}
Augmented::~Augmented() = default;

/// Assign augmented variables to all individuals.
void Augmented::build(const hsp::str_table &input)
{
    this->deftable_ = input;
    helper_assert_no_reserved_();

    // eventually append infectivity profile definitions
    pInfecProfile_->transfer_augm_definitions(this->deftable_);

    // skip if no definitions are present
    if (deftable_.empty())
        return;

    // append discrete infectivity profile definitions if any
    if (pInfecProfile_->is_augmented() && pInfecProfile_->is_discrete_requested())
        for (const auto a : {hsp::infprof::INSTANT, hsp::infprof::CUMUL})
            user_list_discrete.emplace_back(a);

    // check entries
    const vector<string> file_headers = {"name", "status",
                                         "func", "fcumul",
                                         "value", "const", "discr",
                                         "model", hsp::params::MIN, hsp::params::MAX,
                                         "proposal", "proposal_scale"};
    const auto headers = hsp::extract_keys(deftable_);
    hsp::str::assert_no_missing_entries(headers, file_headers);
    hsp::str::assert_no_unused_entries(headers, file_headers);

    helper_build_definition_by_status_();
    helper_build_definition_order_();
    helper_build_model_status_maps_();
    user_list_discrete = hsp::str::unique(user_list_discrete);

    // allocate memory for augmented variables to each individual
    for (auto &house : *pDataWide_)
        for (auto &member : house)
            clear_all(member);
}

// check augmented name
void Augmented::helper_assert_no_reserved_()
{
    if (deftable_.empty() || !deftable_.count("name"))
        return;

    for (const auto &a : deftable_["name"])
    {
        if (hsp::na::is(a))
            throw std::invalid_argument(
                "Found an augmented variable with undefined name.");
        if (a == hsp::infprof::INSTANT || a == hsp::infprof::CUMUL)
            throw std::invalid_argument(
                "Tried to use the reserved augmented name '" + a + "'.\nPlease, choose a different name.");
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////  MCMC : Impl : Augmented : assert univocal definitions  ///////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/// Assert that augmented status/generic definitions are non degenerate.
/// Also, assert that status-specific definitions are for infected status only.
void Augmented::helper_build_definition_by_status_()
{
    status_to_def_row.clear(); // drop previous definitions

    for (size_t k = 0; k < deftable_["name"].size(); ++k)
    {
        string s = deftable_["status"][k];
        if (hsp::na::is(s))
            s = ""; // enforce internal NA convention

        else if (!hsp::status::is_infected(hsp::str::toi(s, "augmented status")))
            throw std::invalid_argument("Impossible to associate augmented '" + deftable_["name"][k] +
                                        "' for non infectious status " + s + ".");

        status_to_def_row[s].push_back(k); // store status to definition row index
    }

    // assert that name and status are univocal
    for (const auto &[s, idxs] : status_to_def_row)
    {
        vector<string> name_to_status = {};
        for (size_t k : idxs)
            name_to_status.push_back(deftable_["name"][k]);

        auto degenerate = name_to_status - hsp::str::unique(name_to_status);
        if (!degenerate.empty())
            throw std::invalid_argument(
                "Incompatible definitions found for augmented " + hsp::str::format(degenerate));
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////  MCMC : Impl : Augmented : build variable to sample  ////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/// Deal with the definition 'k' of augmented to be sampled (const = FALSE).
/// NOTE: this assume definitions are univocal and status correct.
void Augmented::helper_build_definition_order_()
{
    // drop previous definitions
    status_map_sample.clear();
    list_full.clear();
    list_export.clear();
    list_draws.clear();

    // NOTE: NA (i. e. generic) does not necessarily need checking.
    if (pDataWide_->all_infected_status_list.empty())
        set_definition_dependency_order_for_("");
    else
        for (int status : pDataWide_->all_infected_status_list)
            set_definition_dependency_order_for_(std::to_string(status));

    // finalize augmented lists
    list_full = hsp::str::unique(list_full + list_export + list_draws); // NOTE: must be first!
    list_export = hsp::str::drop_private(hsp::str::unique(list_export + list_draws));
    list_draws = hsp::str::drop_private(hsp::str::unique(list_draws));
}

/// Set dependecies ordering for ininitializiation fo given status.
/// At the same time assert that definition are complet and acyclic.
void Augmented::set_definition_dependency_order_for_(const string &s)
{
    hsp::str_table tree;           // populate augmented dependency tree
    hsp::dict<size_t> name_to_row; // ... and reverse dictionary name-to-row-index
    for (size_t k : get_unsorted_definition_rows_(s, true))
    {
        vector<string> xdepen = {};
        for (const string c : {"value", "model", "func", hsp::params::MIN, hsp::params::MAX})
            xdepen += hsp::parse::Definition(deftable_[c][k]).list_augmented_depen();
        tree[deftable_["name"][k]] = hsp::str::unique(xdepen);
        name_to_row[deftable_["name"][k]] = k;
    }

    // check dependecies and translate order from name to deftable row index
    try
    {
        status_to_def_row[s] = {}; // NOTE: recycle from build to provide order

        // Dropping spurios private definition for no public definition depends on it
        for (const auto &a : extract_all_dependencies_order(tree))
            if (!hsp::str::is_private(a) ||
                !hsp::str::drop_private(extract_dependents_order(a, tree)).empty())
                status_to_def_row[s].push_back(name_to_row[a]);
    }
    catch (const std::exception &e)
    {
        throw std::invalid_argument(string(e.what()) + " for status " + s + ".");
    }

    // define usage lists respecting definition order
    for (size_t k : status_to_def_row[s])
        set_sample_map_for_(s, k);
}

/// Extract status-specific definitions row indexes.
vector<size_t> Augmented::get_unsorted_definition_rows_(const string &s,
                                                        bool is_pre_process) const
{
    // load the list of status-specific definition rows
    if (status_to_def_row.count(s))
    {
        vector<size_t> allrows = status_to_def_row.at(s);

        // extract names
        vector<string> name_to_status = {};
        for (size_t k : allrows)
            name_to_status.push_back(deftable_.find("name")->second[k]);

        // add generic definitions to current list (if not generic already)
        if (is_pre_process && !hsp::na::is(s) && status_to_def_row.count(""))
            for (size_t k : status_to_def_row.at(""))
                // skip generic if status-specific already exist
                if (!(name_to_status % deftable_.find("name")->second[k]))
                    allrows.push_back(k);

        return allrows;
    }

    if (status_to_def_row.count(""))
        // use generic definitions if status-specific don't exist
        return status_to_def_row.at("");

    throw std::invalid_argument(
        "No definitions were compatible with augmented status " + s + ".");
}

/// Set sample map and relative usage list.
void Augmented::set_sample_map_for_(const string &s, size_t k)
{
    const auto &a = deftable_["name"][k];

    if (hsp::str::tobool(deftable_["const"][k]))
        list_export.push_back(a);

    else if (!hsp::na::is(deftable_["func"][k]))
        set_function_or_identity_for_(k);

    else if (hsp::na::is(deftable_["value"][k]) && hsp::na::is(deftable_["model"][k]))
        throw std::invalid_argument(
            "Impossible to initialize augmented '" + a + "' for status " + s + " without a value or model.");

    else if (hsp::na::is(deftable_["proposal"][k]))
        throw std::runtime_error(
            "No proposal was provided for non-constant definition of augmented '" + a +
            "' for status " + s + ".");
    else
    {
        // read option for monte carlo proposer
        list_draws.push_back(a);
        if (!hsp::na::is(s))
            status_map_sample[hsp::str::toi(s)].push_back(a);
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////  MCMC : Impl : Augmented : build functions  /////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/// Build augmented function and inspect for constant.
void Augmented::set_function_or_identity_for_(size_t k)
{
    const auto &a = deftable_["name"][k];
    const hsp::parse::Definition parser(deftable_["func"][k]);

    if (parser.is_terminal())
    {
        bool is_identity_const = true;

        // check the type of requested equality
        if (hsp::sources::is_parameter(parser.root.type))
            is_identity_const = pParameter_->at(parser.root.value).is_constant();

        else if (hsp::sources::is_augmented(parser.root.type))
            is_identity_const = false; // UPDATE: use an augmented dependecies tree to see if constant

        if (is_identity_const)
            list_export.push_back(a);
        else
            list_draws.push_back(a);
    }
    else if (parser.is_callable_implicitly())
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
void Augmented::helper_build_model_status_maps_()
{
    // drop previous definitions
    status_map_model.clear();

    for (int status : pDataWide_->all_infected_status_list)
    {
        hsp::str_table tree;          // populate augmented dependency tree
        vector<string> without_model; // list of augmented without models
        for (const size_t k : get_unsorted_definition_rows_(std::to_string(status), false))
        {
            const auto &a = deftable_["name"][k];
            vector<string> xdepen = {};
            for (const auto c : {"value", "model", "func", hsp::params::MIN, hsp::params::MAX})
                xdepen += hsp::parse::Definition(deftable_[c][k]).list_augmented_depen();

            if (hsp::na::is(deftable_["model"][k]))
                without_model.push_back(a);
            else
                tree[a] = hsp::str::unique(xdepen);
        }

        // prune dependency tree from augmented without model
        for (auto &[key, vec] : tree)
            vec = vec - without_model;

        // assign order of dependents for augmented model computation
        for (const auto &[key, vec] : tree)
            status_map_model[status].emplace(key, extract_dependents_order(key, tree));
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////  MCMC : Impl : Augmented : assign augmented variables  ////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/// Clear all augmented definition of individual.
void Augmented::clear_all(Individual &member)
{
    for (const auto &a : list_full)
        member.clear_augm(a);
}

/// Initialize augmented dependecies for evaluations.
/// NOTE: this must not be called before simulation which build on demand.
void Augmented::init()
{
    if (!is_any())
        return;

    for (auto &house : *pDataWide_)
        for (auto &member : house.infected())
        {
            init_all(member);
            init_model_logratios(member);
        }
}

/// Assign all augmented variables to single individual.
void Augmented::init_all(Individual &member)
{
    // loop over augmented variable definitions to assign placheholders
    for (size_t k : init_get_definition_order_(member))
        helper_init_attributes_(member, k);
}

/// Assign augmented variables only if missing (used for simulations).
void Augmented::init_missing(Individual &member)
{
    // loop over augmented variable definitions to assign placheholders
    for (size_t k : init_get_definition_order_(member))
        if (member.get_augm_variable(deftable_["name"][k]).is_empty())
            helper_init_attributes_(member, k);
}

// Assign all augmented models to single individual.
void Augmented::init_model_logratios(Individual &member)
{
    try
    {
        member.init_augm_model_logratios(status_map_model[member.status()]);
    }
    catch (const std::exception &e)
    {
        throw std::runtime_error(
            "While initializing augmented model for " + member.err_name() + ":\n" +
            string(e.what()));
    }
}

/// Get the definition order for current member according to its status.
vector<size_t> Augmented::init_get_definition_order_(Individual &member)
{
    const string s = std::to_string(member.status());

    // NOTE: if status hasn't been observed in data, build on demand
    if (member.is_infected_status() && !(pDataWide_->all_infected_status_list % member.status()))
    {
        set_definition_dependency_order_for_(s);
        pDataWide_->all_infected_status_list.push_back(member.status());
    }

    if (status_to_def_row.count(s))
        return status_to_def_row[s];
    if (status_to_def_row.count(""))
        return status_to_def_row[""];
    if (member.is_infected_status())
        throw std::invalid_argument(
            "No augmented definitions were found for " + member.err_name() + " with status " + s + ".");
    throw std::invalid_argument(
        "Generic augmented definition for non-infectious " + member.err_name() + " was missing.");
}

/// Assign model arguments and range. Same for function, otherwise assign proposer.
void Augmented::helper_init_attributes_(Individual &member, size_t k)
{
    const auto &a = deftable_["name"][k];

    try
    {
        member.init_augm(a, pParameter_->read_numeric(deftable_["value"][k], &member).get(),
                         hsp::str::tobool(deftable_["const"][k]),
                         deftable_["func"][k], hsp::str::tobool(deftable_["fcumul"][k]),
                         deftable_["model"][k]);
        auto &xaugm = member.get_augm_variable(a);

        init_model_attributes_(member, k); // initialize model
        init_func_attributes_(member, k);  // initialize function

        if (xaugm.is_pointer()) // pointer (case of identity function)
            xaugm.transfer(pParameter_->read_attribute(deftable_["func"][k], &member));

        // assign proposer if variable to be sampled
        if (!xaugm.is_constant())
            xaugm.init_proposer(deftable_["proposal"][k],
                                hsp::str::tod(deftable_["proposal_scale"][k]));

        // check if satisfies discrete requirement, if any
        if (!user_list_discrete.empty() && (user_list_discrete % a) && !xaugm.is_discrete())
            throw std::invalid_argument(
                "Definition of " + xaugm.err_name() +
                " did not return a discrete object as requested by the model.");
    }
    catch (const std::exception &e)
    {
        throw std::runtime_error("While initializing augmented '" + a +
                                 "' for " + member.err_name() + ":\n" +
                                 string(e.what()));
    }
}

/// Auxiliary method to assign attributes to augmented model.
void Augmented::init_model_attributes_(Individual &member, size_t k)
{
    auto &xaugm = member.get_augm_variable(deftable_["name"][k]);

    if (!xaugm.is_model())
        return;

    for (const auto &child : hsp::parse::Definition(deftable_["model"][k]).root.children)
        xaugm.init_distr_attr(child.key, pParameter_->read_attribute(child, &member));

    for (const auto key : {hsp::params::MIN, hsp::params::MAX})
        xaugm.init_distr_attr(key, pParameter_->read_attribute(deftable_[key][k], &member));

    xaugm.init_distr_reshaping(hsp::str::tobool(deftable_["discr"][k]));
}

/// Auxiliary method to assign attributes to augmented function.
void Augmented::init_func_attributes_(Individual &member, size_t k)
{
    auto &xaugm = member.get_augm_variable(deftable_["name"][k]);

    if (!xaugm.is_function())
        return;

    // initialize function
    for (const auto &child : hsp::parse::Definition(deftable_["func"][k]).root.children)
    {
        xaugm.init_func_attr(child.key, pParameter_->read_attribute(child, &member));
        pParameter_->make_cache_dependent(child);
    }

    // initialize range
    for (const auto key : {hsp::params::MIN, hsp::params::MAX})
        xaugm.init_func_attr(key, pParameter_->read_attribute(deftable_[key][k], &member));

    // initialize discretization and truncation
    xaugm.init_func_reshaping(hsp::str::tobool(deftable_["discr"][k]));
}

/////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////  MCMC : Impl : Augmented : user reset  /////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

/// Reset augmented initial values for given individual.
void Augmented::user_reset(hsp::rnd_gen &prng, Individual &member)
{
    for (const auto &a : status_map_sample[member.status()])
        member.get_augm_variable(a).user_reset(prng);
}

/// Initialize augmented initial values for all individuals.
void Augmented::user_reset(hsp::rnd_gen &prng)
{
    if (!is_any())
        return;

    for (auto &house : *pDataWide_)
        for (auto &member : house.infected())
            user_reset(prng, member);
}

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////  MCMC : Impl : Augmented : sample draws  ////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/// Set individual id for augmented draws follow-up.
void Augmented::init_individual_sample(vector<string> &id_draws)
{
    smpl_individuals.clear();

    if (id_draws.empty())
        return;

    if (!is_any_sample())
    {
        hsp::msg::warning(
            "Ignored id draws. No augmented variables to sample were defined.");
        id_draws.clear();
        return;
    }

    smpl_individuals.reserve(id_draws.size());
    vector<string> already_stored_id = {};
    for (const auto &house : *pDataWide_)
        for (const auto &member : house)
            if (id_draws % member.id())
            {
                // check id across all individuals
                if (already_stored_id % member.id())
                    throw std::invalid_argument(
                        "Multiple individuals shared the same id '" + member.id() + "'.");
                already_stored_id.push_back(member.id());

                // check status of chosen individual
                if (member.is_infected_status())
                    smpl_individuals.emplace_back(&member);
                else
                    hsp::msg::warning(
                        "Ignored draws of non infected " + member.err_name() + ".");
            }

    // check all user request are satisfied
    hsp::str::assert_no_missing_entries(already_stored_id, id_draws, "trace id");
}

/// Initialize sample for individual augmented variable followup.
void Augmented::init_draws(unsigned long sample_size)
{
    if (!is_any_draw())
        return;

    smpl_values.clear();
    smpl_values.resize(smpl_individuals.size());
    for (size_t k = 0; k < smpl_individuals.size(); ++k)
    {
        smpl_values[k].clear();
        smpl_values[k].reserve(list_draws.size());
        for (const auto &a : list_draws)
            smpl_values[k][a].assign(sample_size, hsp::na::DOUBLE);
    }
}

/// Store sample for individual augmented variable followup.
void Augmented::store_draws(unsigned long smpl_count)
{
    if (!is_any_draw())
        return;

    for (size_t k = 0; k < smpl_individuals.size(); ++k)
    {
        auto &smpl_augm_value_k = smpl_values[k];
        const auto member_star = smpl_individuals[k];
        for (const auto &a : list_draws)
            smpl_augm_value_k[a][smpl_count] = member_star->aug(a);
    }
}
