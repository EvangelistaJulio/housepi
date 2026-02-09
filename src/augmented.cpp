/* housepi - Epidemiological Analysis of Household Outbreaks.
 *
 * Copyright © 2026 Institut Pasteur (Paris), CNRS and INSERM.
 * Distributed under the terms of the GNU General Public License (GPLv3).
 * See the LICENSE file for details.
 */

#include <algorithm>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

#include "tools/keys.h"
#include "tools/loglike.h"
#include "tools/messages.h"
#include "tools/nan.h"
#include "tools/operators.h"
#include "tools/parser.h"
#include "tools/random.h"
#include "tools/read.h"
#include "tools/string.h"
#include "tools/types.h"

#include "augmented.h"
#include "data_wide.h"
#include "function.h"
#include "household.h"
#include "individual.h"
#include "infec_profile.h"
#include "parameter.h"
#include "variable.h"

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////  Dependency Tree Helpers  //////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

namespace
{
    /// Check augmented names and assert no reserved names are used.
    void helper_assert_augm_input(const hsp::str_table &input)
    {
        const std::vector<std::string> file_headers = {"name", "status",
                                                       "func", "fcumul",
                                                       "value", "const", "discr",
                                                       "model", "min", "max",
                                                       "proposal", "proposal_scale"};
        const std::vector<std::string> headers = hsp::read::extract_keys(input);
        hsp::str::assert_no_missing_entries(headers, file_headers);
        hsp::str::assert_no_unused_entries(headers, file_headers);

        for (const std::string &augmstr : input.at("name"))
        {
            if (hsp::na::is(augmstr))
                throw std::invalid_argument(
                    "Found an augmented variable with undefined name.");
            if (augmstr == hsp::infprof::AUGM_FUNC)
                throw std::invalid_argument(
                    "Tried to use the reserved augmented name '" +
                    augmstr + "'.\nPlease, choose a different one.");
        }
    }

    /// Helper to build augmented definition order and assert for acyclic.
    void walk_dependency_tree(const std::string &entry,
                              const hsp::str_table &tree,
                              std::vector<std::string> &visited,
                              std::vector<std::string> &recursion_stack,
                              std::vector<std::string> &ordered)
    {
        if (!tree.count(entry))
            throw std::invalid_argument(
                "Could not find a definition for augmented '" + entry + "'.");

        if (std::find(recursion_stack.begin(), recursion_stack.end(), entry) != recursion_stack.end())
            throw std::invalid_argument(
                "Found inadmissible cyclic definitions for augmented '" + entry + "'");

        if (std::find(visited.begin(), visited.end(), entry) != visited.end())
            return;

        recursion_stack.push_back(entry);
        for (const std::string &dep : tree.at(entry))
            walk_dependency_tree(dep, tree, visited, recursion_stack, ordered);

        recursion_stack = recursion_stack - std::vector<std::string>({entry});
        visited.push_back(entry);
        ordered.push_back(entry);
    }

    /// Get topological order of definitions dependencies.
    std::vector<std::string> extract_all_dependencies_order(const hsp::str_table &graph)
    {
        // inspect dependency tree for cycles and missing entries
        std::vector<std::string> visited, recursion_stack, ordered;
        for (const auto &[a, vec] : graph)
            if (std::find(visited.begin(), visited.end(), a) == visited.end())
                walk_dependency_tree(a, graph, visited, recursion_stack, ordered);
        return ordered;
    }

    /// Get topological order of definitions depending on 'start' (start included).
    std::vector<std::string> extract_dependents_order(const std::string &start, const hsp::str_table &graph)
    {
        // construct a reversed graph for dependents
        hsp::str_table reverse_graph;
        for (const auto &[node, dependents] : graph)
        {
            for (const std::string &dep : dependents)
                reverse_graph[dep].push_back(node);

            // ensure isolated nodes are present
            if (!reverse_graph.count(node))
                reverse_graph[node] = {};
        }

        // walk dependents of start
        std::vector<std::string> visited, recursion_stack, ordered;
        for (const std::string &dep : reverse_graph.at(start))
            walk_dependency_tree(dep, reverse_graph, visited, recursion_stack, ordered);

        std::reverse(ordered.begin(), ordered.end());
        return std::vector<std::string>({start}) + ordered;
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////  Augmented : build  //////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

Augmented::Augmented(const std::vector<std::string> &keys,
                     DataWide *dataptr,
                     Parameter *paramptr,
                     InfecProfile *infprofptr)
    : pRegistry(std::make_unique<hsp::key::Registry>(keys)),
      pDataWide_(dataptr),
      pParameter_(paramptr),
      pInfecProfile_(infprofptr) {}

Augmented::~Augmented() = default;

void Augmented::build(const hsp::str_table &input)
{
    pRegistry->clear();

    if (input.empty() && !pInfecProfile_->is_augmented())
        return;

    helper_assert_augm_input(input);      // check entries
    deftable_ = input;                    // store user definitions
    helper_build_append_infec_profile_(); // add infectivity profile
    helper_build_definition_by_status_(); // process definition by status
    helper_build_definition_order_();     // proicess definition order

    // allocate memory for augmented variables to each individual
    for (Household &house : *pDataWide_)
        for (Individual &member : house)
            clear_all(member);
}

void Augmented::helper_build_append_infec_profile_()
{
    if (!pInfecProfile_->is_augmented())
        return;

    for (const std::string &col : hsp::read::extract_keys(deftable_))
        deftable_[col] += pInfecProfile_->augm_profile_definitions_.at(col);

    pRegistry->add(hsp::infprof::AUGM_FUNC);
    pInfecProfile_->augm_key = pRegistry->from_string(hsp::infprof::AUGM_FUNC);

    // append discrete infectivity profile definitions if any
    if (pInfecProfile_->is_discrete_requested())
        user_discrete_set.emplace(pInfecProfile_->augm_key);
}

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////  Augmented : assert univocal definitions  //////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

void Augmented::helper_build_definition_by_status_()
{
    status_to_def_row.clear();

    for (size_t k = 0; k < deftable_["name"].size(); ++k)
    {
        // fill registry with optional augmented definitions
        const std::string augmstr = deftable_["name"][k];
        if (!pRegistry->count(augmstr))
            pRegistry->add(augmstr);

        // check associated status
        int status = hsp::status::MISSING;
        const std::string &s = deftable_["status"][k];
        if (!hsp::na::is(s))
        {
            status = hsp::str::toi(s, "augmented status");
            if (!hsp::status::is_infected(status))
                throw std::invalid_argument(
                    "Impossible to associate augmented '" + augmstr +
                    "' for non infectious status " + s + ".");
        }
        status_to_def_row[status].push_back(k); // store status to definition row index
    }

    // assert that name and status are univocal
    for (const auto &[status, idxs] : status_to_def_row)
    {
        std::vector<std::string> name_to_status = {};
        for (size_t k : idxs)
            name_to_status.push_back(deftable_["name"][k]);

        const std::vector<std::string> degenerate = name_to_status - hsp::str::unique(name_to_status);
        if (!degenerate.empty())
            throw std::invalid_argument(
                "Incompatible definitions found for augmented " + hsp::str::format(degenerate));
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////  Augmented : build variable to sample  ///////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

void Augmented::helper_build_definition_order_()
{
    status_to_sample.clear();
    list_export.clear();
    list_draws.clear();

    // NOTE: NA (i. e. generic) does not necessarily need checking.
    if (pDataWide_->all_infected_statuses.empty())
        set_definition_dependency_order_for_(hsp::status::MISSING);
    else
        for (int status : pDataWide_->all_infected_statuses)
            set_definition_dependency_order_for_(status);
}

void Augmented::set_definition_dependency_order_for_(int status)
{
    // columns presenting definitions possibly containing sources
    const auto columns = {"value", "model", "func", "min", "max"};

    hsp::str_table tree;           // populate augmented dependency tree
    hsp::dict<size_t> name_to_row; // ... and reverse dictionary name-to-row-index
    for (size_t k : get_unsorted_definition_rows_(status, true))
    {
        std::vector<std::string> xdepen = {};
        for (const auto c : columns)
            xdepen += hsp::parse::Definition(deftable_[c][k]).list_augmented_depen();
        tree[deftable_["name"][k]] = hsp::str::unique(xdepen);
        name_to_row[deftable_["name"][k]] = k;
    }

    // check dependecies and translate order from name to deftable row index
    try
    {
        status_to_def_row[status].clear(); // NOTE: recycled from build to store order

        // Ignore spurios private definition if no public definition depends on it
        // NOTE: ensure private infectivity is stored.
        for (const std::string &augmstr : extract_all_dependencies_order(tree))
            if (augmstr == hsp::infprof::AUGM_FUNC ||
                !hsp::str::is_private(augmstr) ||
                !hsp::str::drop_private(extract_dependents_order(augmstr, tree)).empty())
                status_to_def_row[status].push_back(name_to_row.at(augmstr));
    }
    catch (const std::exception &e)
    {
        throw std::invalid_argument(
            std::string(e.what()) + " for status " + std::to_string(status) + ".");
    }

    // define usage lists respecting definition order
    for (size_t k : status_to_def_row.at(status))
        set_sample_map_for_(status, k);
}

std::vector<size_t> Augmented::get_unsorted_definition_rows_(int status, bool is_pre_process) const
{
    // load the list of status-specific definition rows
    if (status_to_def_row.count(status))
    {
        std::vector<size_t> allrows = status_to_def_row.at(status);

        // extract names
        std::vector<std::string> name_to_status = {};
        for (size_t k : allrows)
            name_to_status.push_back(deftable_.find("name")->second[k]);

        // add generic definitions to current list (if not generic already)
        // skip generic if a status-specific already exist
        if (is_pre_process &&
            hsp::status::is_infected(status) &&
            status_to_def_row.count(hsp::status::MISSING))
            for (size_t k : status_to_def_row.at(hsp::status::MISSING))
                if (std::find(name_to_status.begin(), name_to_status.end(), deftable_.find("name")->second[k]) == name_to_status.end())
                    allrows.push_back(k);

        return allrows;
    }

    if (status_to_def_row.count(hsp::status::MISSING))
        // use generic definitions if status-specific don't exist
        return status_to_def_row.at(hsp::status::MISSING);

    throw std::invalid_argument(
        "No definitions were compatible with augmented status " +
        std::to_string(status) + ".");
}

void Augmented::set_sample_map_for_(int status, size_t k)
{
    const std::string &augmstr = deftable_["name"][k];
    size_t a = pRegistry->from_string(augmstr);

    if (hsp::str::tobool(deftable_["const"][k]))
    {
        if (!hsp::str::is_private(augmstr))
            list_export.emplace(a);
    }
    else if (!hsp::na::is(deftable_["func"][k]))
    {
        set_function_or_identity_for_(k);
    }
    else if (hsp::na::is(deftable_["value"][k]) && hsp::na::is(deftable_["model"][k]))
    {
        throw std::invalid_argument(
            "Impossible to initialize augmented '" + augmstr +
            "' for status " + std::to_string(status) + " without a value or model.");
    }
    else if (hsp::na::is(deftable_["proposal"][k]))
    {
        throw std::runtime_error(
            "No proposal was provided for non-constant definition of augmented '" +
            augmstr + "' for status " + std::to_string(status) + ".");
    }
    else
    {
        if (!hsp::str::is_private(augmstr))
        {
            list_draws.emplace(a);
            list_export.emplace(a);
        }

        if (hsp::status::is_infected(status))
            status_to_sample[status].push_back(a);
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////  Augmented : build functions  ////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

void Augmented::set_function_or_identity_for_(size_t k)
{
    const std::string &augmstr = deftable_["name"][k];
    size_t a = pRegistry->from_string(augmstr);
    const hsp::parse::Definition parser(deftable_["func"][k]);

    if (parser.is_terminal())
    {
        bool is_identity_const = true;

        // check the type of requested equality
        if (hsp::parse::is_parameter(parser.root.type))
            is_identity_const = pParameter_->at(pParameter_->pRegistry->from_string(parser.root.value)).is_constant();

        else if (hsp::parse::is_augmented(parser.root.type))
            is_identity_const = false; // UPDATE: use an augmented dependecies tree to see if constant

        if (!hsp::str::is_private(augmstr))
        {
            list_export.emplace(a);
            if (!is_identity_const)
                list_draws.emplace(a);
        }
    }
    else if (parser.is_callable_implicitly("x") && !hsp::str::is_private(augmstr))
    {
        list_draws.emplace(a);
        list_export.emplace(a);
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////  Augmented : assign augmented variables  ///////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

void Augmented::clear_all(Individual &member)
{
    size_t n = pRegistry->size();
    member.augm_vec_.resize(n);
    for (size_t a = 0; a < n; ++a)
    {
        member.augm_vec_[a].name(pRegistry->as_string(a));
        member.augm_vec_[a].set_err_type("augmented");
    }
}

void Augmented::init()
{
    if (pRegistry->size() == 0)
        return;

    init_model_list_();

    for (Household &house : *pDataWide_)
        for (size_t i : house.infected_hhindexes())
        {
            Individual &member = house.member(i);
            init_missing(member, true);
            member.augm_idx_status_ = status_to_index.at(member.status());
        }
}

void Augmented::init_missing(Individual &member, bool force)
{
    // loop over augmented variable definitions to assign placheholders
    for (size_t k : init_get_definition_order_(member))
        if (force || member.augm_vec_[pRegistry->from_string(deftable_["name"][k])].is_empty())
            helper_init_attributes_(member, k);
}

std::vector<size_t> Augmented::init_get_definition_order_(Individual &member)
{
    int status = member.status();

    if (member.is_infected_status() &&
        pDataWide_->all_infected_statuses.count(status) == 0)
    {
        set_definition_dependency_order_for_(status);
        pDataWide_->all_infected_statuses.emplace(status);
    }

    if (status_to_def_row.count(status))
        return status_to_def_row[status];
    if (status_to_def_row.count(hsp::status::MISSING))
        return status_to_def_row[hsp::status::MISSING];

    if (member.is_infected_status())
        throw std::invalid_argument(
            "No augmented definitions were found for " +
            member.err_name() + " with status " + std::to_string(status) + ".");

    throw std::invalid_argument(
        "Generic augmented definition for non-infectious " +
        member.err_name() + " was missing.");
}

void Augmented::helper_init_attributes_(Individual &member, size_t k)
{
    const std::string &augmstr = deftable_["name"][k];
    size_t a = pRegistry->from_string(augmstr);

    try
    {
        Variable &xaugm = member.augm_vec_[a];
        xaugm.clear();

        if (hsp::str::tobool(deftable_["const"][k]))
        {
            xaugm.init_numeric(
                pParameter_->read_numeric(deftable_["value"][k], &member).value());
            xaugm.enforce_user_const(); // enforce user constant choice
        }
        else
            init_func_attributes_(member, k);

        // initialize model
        init_model_attributes_(member, k);

        // assign proposal if variable to be sampled
        // NOTE: not necessarily according to a model
        if (!xaugm.is_constant())
            xaugm.init_proposal(deftable_["proposal"][k],
                                hsp::str::tod(deftable_["proposal_scale"][k]));

        // check if satisfies discrete requirement, if any
        if (user_discrete_set.count(a) &&
            !xaugm.is_discrete())
            throw std::invalid_argument(
                "Definition of " + xaugm.err_name() +
                " did not return a discrete object as requested by model.");
    }
    catch (const std::exception &e)
    {
        throw std::runtime_error(
            "While initializing augmented '" + augmstr +
            "' for " + member.err_name() + ":\n" + std::string(e.what()));
    }
}

void Augmented::init_model_attributes_(Individual &member, size_t k)
{
    const std::string &augmstr = deftable_["name"][k];
    Variable &xaugm = member.augm_vec_[pRegistry->from_string(augmstr)];
    hsp::parse::Definition parser(deftable_["model"][k]);
    xaugm.init_distr(parser.root.func);

    if (!xaugm.is_model())
        return;

    for (const auto &child : parser.root.children)
        xaugm.init_distr_attr(child.key, pParameter_->read_attribute(child, &member));

    for (const auto key : {"min", "max"})
        xaugm.init_distr_attr(key, pParameter_->read_attribute(deftable_[key][k], &member));

    xaugm.init_distr_reshaping(hsp::str::tobool(deftable_["discr"][k]));
}

void Augmented::init_func_attributes_(Individual &member, size_t k)
{
    const std::string &augmstr = deftable_["name"][k];
    Variable &xaugm = member.augm_vec_[pRegistry->from_string(augmstr)];
    hsp::parse::Definition parser(deftable_["func"][k]);

    // pointer mode
    if (parser.is_terminal())
    {
        xaugm.transfer(pParameter_->read_attribute(parser.root, &member));
        return;
    }

    // function mode
    bool is_primitive = hsp::str::tobool(deftable_["fcumul"][k]);
    xaugm.init_function(compose_func_registries(
        parser.root.func, "function",
        {"user", "transcendental", "arithmetic",
         "continuos", "discrete", "infecprofile"},
        is_primitive));

    // attribute
    for (const auto &child : parser.root.children)
    {
        xaugm.init_func_attr(child.key, pParameter_->read_attribute(child, &member));
        pParameter_->make_cache_dependent(child);
    }

    // range
    for (const auto key : {"min", "max"})
        xaugm.init_func_attr(key, pParameter_->read_attribute(deftable_[key][k], &member));

    // discretization and truncation
    xaugm.init_func_reshaping(hsp::str::tobool(deftable_["discr"][k]));
}

/////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////  Augmented : user reset  ////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

void Augmented::user_reset(hsp::rnd_gen &prng, Individual &member)
{
    for (size_t a : status_to_sample[member.status()])
        member.augm_vec_[a].user_reset(prng);

    for (size_t a = 0; a < pRegistry->size(); ++a)
        member.augm_vec_[a].optimize();
}

void Augmented::user_reset(hsp::rnd_gen &prng)
{
    if (pRegistry->size() == 0)
        return;

    for (Household &house : *pDataWide_)
        for (size_t i : house.infected_hhindexes())
            user_reset(prng, house.member(i));
}

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////  Augmented : build augmented models  ////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

void Augmented::init_model_list_()
{
    status_to_index.clear();
    model_blocks_.clear();
    model_blocks_.resize(pDataWide_->all_infected_statuses.size());

    size_t i = 0;
    for (int status : pDataWide_->all_infected_statuses)
    {
        ModelBlock_ &block = model_blocks_[i];
        block.depen.resize(pRegistry->size());

        // associate blocks for ordered model lists
        const hsp::str_table &ordered_models = get_ordered_models_(status);
        if (!ordered_models.empty())
        {
            is_any_model_ = true;
            block.models = pRegistry->from_string(hsp::read::extract_keys(ordered_models));
            for (const auto &[augmstr, others] : ordered_models)
                block.depen[pRegistry->from_string(augmstr)] = pRegistry->from_string(others);
        }
        // update status to index map
        status_to_index.emplace(status, i++);
    }
}

hsp::str_table Augmented::get_ordered_models_(int status) const
{
    // columns presenting definitions possibly containing sources
    const auto columns = {"value", "model", "func", "min", "max"};

    hsp::str_table tree;                    // populate augmented dependency tree
    hsp::str_table ordered_models;          // store hierarchical order of models
    std::vector<std::string> without_model; // list of augmented without models
    for (size_t k : get_unsorted_definition_rows_(status, false))
    {
        const std::string &aumgstr = deftable_.at("name")[k];
        std::vector<std::string> xdepen = {};
        for (auto c : columns)
            xdepen += hsp::parse::Definition(deftable_.at(c)[k]).list_augmented_depen();

        if (hsp::na::is(deftable_.at("model")[k]))
            without_model.push_back(aumgstr);
        else
            tree[aumgstr] = hsp::str::unique(xdepen);
    }

    // prune dependency tree from augmented without model
    for (auto &[key, vec] : tree)
        vec = vec - without_model;

    // assign order of dependents for augmented model computation
    for (const auto &[key, vec] : tree)
        ordered_models.emplace(key, extract_dependents_order(key, tree));

    return ordered_models;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////  Augmented : log-likelihood  //////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

double Augmented::materialize_augm_model_log_lik(Individual &member)
{
    const ModelBlock_ &block = model_blocks_[member.augm_idx_status_];
    if (block.models.empty())
        return 0.;

    double new_ll = 0.0;
    for (size_t a : block.models)
        new_ll += member.augm_vec_[a].get_distr_logpdf();
    member.augm_log_lik_.set(new_ll);

    return new_ll;
}

double Augmented::materialize_augm_model_logratio(Individual &member, size_t a)
{
    const ModelBlock_ &block = model_blocks_[member.augm_idx_status_];
    if (block.depen[a].empty())
        return 0.;

    double new_ll = member.augm_log_lik_.value(); // store old value in buffer
    for (size_t b : block.depen[a])             // depen includes 'a'
        new_ll += member.augm_vec_[b].get_distr_logratio();

    member.augm_log_lik_.set(new_ll);

    return member.augm_log_lik_.logratio();
}

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////  Augmented : sample draws  ///////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

void Augmented::init_individual_sample(std::vector<std::string> &id_draws)
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
    std::vector<std::string> existing_ids = {};
    for (const Household &house : *pDataWide_)
        for (const Individual &member : house)
            if (std::find(id_draws.begin(), id_draws.end(), member.id()) != id_draws.end())
            {
                // check id across all individuals
                if (std::find(existing_ids.begin(), existing_ids.end(), member.id()) != existing_ids.end())
                    throw std::invalid_argument(
                        "Multiple individuals shared the same id '" + member.id() + "'.");
                existing_ids.push_back(member.id());

                // check status of chosen individual
                if (member.is_infected_status())
                    smpl_individuals.emplace_back(&member);
                else
                    hsp::msg::warning(
                        "Ignored draws of non infected " + member.err_name() + ".");
            }

    // check all user request are satisfied
    hsp::str::assert_no_missing_entries(existing_ids, id_draws, "draws id");
}

void Augmented::init_draws(unsigned long sample_size)
{
    if (!is_any_draw())
        return;

    size_t n = list_draws.size();
    smpl_values.clear();
    smpl_values.resize(smpl_individuals.size());
    for (size_t i = 0; i < smpl_individuals.size(); ++i)
    {
        auto &member_matrix = smpl_values[i];
        member_matrix.reserve(n);
        for (size_t j = 0; j < n; ++j)
            member_matrix.emplace_back(sample_size, hsp::na::DOUBLE);
    }
}

void Augmented::store_draws(unsigned long smpl_count)
{
    if (!is_any_draw())
        return;

    for (size_t i = 0; i < smpl_individuals.size(); ++i)
    {
        const Individual *member_ptr = smpl_individuals[i];
        auto &member_matrix = smpl_values[i];

        size_t c = 0;
        for (size_t a : list_draws)
            member_matrix[c++][smpl_count] = member_ptr->augm_vec_[a].value();
    }
}
