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
#include "augm_container.h"
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
    using namespace hsp;

    /// Check augmented names and assert no reserved names are used.
    void helper_assert_augm_input(const str_table &input)
    {
        const std::vector<std::string> file_headers = {"name", "status",
                                                       "func", "fcumul",
                                                       "value", "const", "discr",
                                                       "model", "min", "max",
                                                       "proposal", "proposal_scale"};
        const std::vector<std::string> headers = read::extract_keys(input);
        str::assert_no_missing_entries(headers, file_headers);
        str::assert_no_unused_entries(headers, file_headers);

        for (const std::string &augmstr : input.at("name"))
        {
            if (na::is(augmstr))
                throw std::invalid_argument(
                    "Found an augmented variable with undefined name.");
            if (augmstr == infprof::AUGM_FUNC)
                throw std::invalid_argument(
                    "Tried to use the reserved augmented name '" +
                    augmstr + "'.\nPlease, choose a different one.");
        }
    }

    /// Helper to build augmented definition order and assert for acyclic.
    void walk_dependency_tree(const std::string &entry,
                              const str_table &tree,
                              std::vector<std::string> &visited,
                              std::vector<std::string> &recursion_stack,
                              std::vector<std::string> &ordered)
    {
        if (!tree.count(entry))
            throw std::invalid_argument(
                "Could not find a definition for augmented '" + entry + "'");

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
    std::vector<std::string> extract_all_dependencies_order(const str_table &graph)
    {
        // inspect dependency tree for cycles and missing entries
        std::vector<std::string> visited, recursion_stack, ordered;
        for (const auto &[a, vec] : graph)
            if (std::find(visited.begin(), visited.end(), a) == visited.end())
                walk_dependency_tree(a, graph, visited, recursion_stack, ordered);
        return ordered;
    }

    /// Get topological order of definitions depending on 'start' (start included).
    std::vector<std::string> extract_dependents_order(const std::string &start, const str_table &graph)
    {
        // construct a reversed graph for dependents
        str_table reverse_graph;
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
////////////////////////////////////  *** Augmented : class ***  //////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

namespace hsp
{
    Augmented::Augmented(DataWide *dataptr,
                         Parameter *paramptr,
                         InfecProfile *infprofptr)
        : pDataWide_(dataptr),
          pParameter_(paramptr),
          pInfecProfile_(infprofptr) {}
    Augmented::~Augmented() = default;

    void Augmented::init_registry(const std::vector<std::string> &keys)
    {
        pRegistry = std::make_unique<key::Registry>(keys);
    }

    void Augmented::build(const str_table &input)
    {
        // drop previous definitions
        is_any_sample_ = false;
        is_any_model_ = false;

        pRegistry->clear();
        map_status_to_def_row_.clear();
        map_status_to_block_.clear();
        status_blocks_.clear();
        list_export.clear();

        // skip if no augmented definitions
        if (input.empty() && !pInfecProfile_->is_augmented())
            return;

        helper_assert_augm_input(input);      // check entries
        deftable_ = input;                    // store user definitions
        helper_build_append_infec_profile_(); // add infectivity profile
        helper_build_definition_by_status_(); // process definition by status
        helper_build_allocation_order_();     // allocation and order
    }

    void Augmented::helper_build_append_infec_profile_()
    {
        if (!pInfecProfile_->is_augmented())
            return;

        for (const std::string &col : read::extract_keys(deftable_))
            deftable_[col] += pInfecProfile_->augm_definitions.at(col);

        pRegistry->add(infprof::AUGM_FUNC);
        pInfecProfile_->augm_key = pRegistry->from_string(infprof::AUGM_FUNC);

        // append discrete infectivity profile definitions if any
        if (pInfecProfile_->is_discrete_requested())
            user_discrete_set_.emplace(pInfecProfile_->augm_key);
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////  Augmented : assert univocal definitions  //////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////

    void Augmented::helper_build_definition_by_status_()
    {
        for (size_t k = 0; k < deftable_["name"].size(); ++k)
        {
            // fill registry with optional augmented definitions
            const std::string augmstr = deftable_["name"][k];
            pRegistry->add(augmstr);

            // check associated status
            int status = status::MISSING;
            const std::string &s = deftable_["status"][k];
            if (!na::is(s))
            {
                status = str::toi(s, "augmented status");
                if (!status::is_infected(status))
                    throw std::invalid_argument(
                        "Impossible to associate augmented '" + augmstr +
                        "' for non infectious status " + s + ".");
            }
            map_status_to_def_row_[status].push_back(k); // store status to definition row index
        }

        // assert that name and status are univocal
        for (const auto &[status, idxs] : map_status_to_def_row_)
        {
            std::vector<std::string> name_to_status = {};
            for (size_t k : idxs)
                name_to_status.push_back(deftable_["name"][k]);

            const std::vector<std::string> degenerate = name_to_status - str::unique(name_to_status);
            if (!degenerate.empty())
                throw std::invalid_argument(
                    "Incompatible definitions found for augmented " + str::format(degenerate));
        }
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////  Augmented : build definition order  ////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////

    void Augmented::helper_build_allocation_order_()
    {
        std::unordered_set<int> all_infected_statuses = {};

        size_t n = pRegistry->size();
        for (Household &house : *pDataWide_)
            for (Individual &member : house)
            {
                // UPDATE: this could be reserved during member construction
                member.augm_vec_.resize(n);
                for (size_t a = 0; a < n; ++a)
                {
                    member.augm_vec_[a].name(pRegistry->as_string(a));
                    member.augm_vec_[a].set_err_type("augmented");
                }
                if (member.is_infected_status())
                    all_infected_statuses.emplace(member.status());
            }

        // first, load status specific user definitions
        // NOTE: NA (i. e. generic) does not necessarily need checking
        for (auto &[status, vec] : map_status_to_def_row_)
            if (status::is_infected(status))
                set_definition_dependency_order_for_(status);

        // associate observed status to corresponding status block
        // if not defined specifically, then associate to generic
        for (int status : all_infected_statuses)
            if (!map_status_to_block_.count(status))
            {
                if (!map_status_to_block_.count(status::MISSING))
                    set_definition_dependency_order_for_(status::MISSING);
                map_status_to_block_.emplace(status, map_status_to_block_.at(status::MISSING));
            }

        // if no observed infected status or status-specific definitions
        // then load generic definitition
        if (all_infected_statuses.empty() ||
            (map_status_to_def_row_.size() == 1 &&
             status::is_missing(map_status_to_def_row_.begin()->first)))
            set_definition_dependency_order_for_(status::MISSING);
    }

    void Augmented::set_definition_dependency_order_for_(int status)
    {
        // skip if already passed
        if (map_status_to_block_.count(status))
            return;
        // throw if not defined
        if (!map_status_to_def_row_.count(status))
            throw std::invalid_argument(
                "No definitions were compatible with augmented status " +
                std::to_string(status) + ".");

        // columns presenting definitions possibly containing sources
        const auto columns = {"value", "model", "func", "min", "max"};

        str_table tree;           // populate augmented dependency tree
        dict<size_t> name_to_row; // ... and reverse dictionary name-to-row-index
        for (size_t k : get_unsorted_definition_rows_(status))
        {
            std::vector<std::string> xdepen = {};
            for (const auto c : columns)
                xdepen += parse::Definition(deftable_[c][k]).list_augmented_depen();
            tree[deftable_["name"][k]] = str::unique(xdepen);
            name_to_row[deftable_["name"][k]] = k;
        }

        // check dependecies and translate order from name to deftable row index
        try
        {
            // NOTE: this method is just called once per status
            map_status_to_block_.emplace(status, status_blocks_.size());
            status_blocks_.emplace_back(status);

            // Ignore spurios private definition if no public definition depends on it
            // NOTE: ensure private infectivity is stored.
            for (const std::string &augmstr : extract_all_dependencies_order(tree))
                if (augmstr == infprof::AUGM_FUNC ||
                    !str::is_private(augmstr) ||
                    !str::drop_private(extract_dependents_order(augmstr, tree)).empty())
                    status_blocks_.back().order.emplace_back(name_to_row[augmstr]);

            // define usage lists respecting definition order
            for (size_t k : status_blocks_.back().order)
                set_sample_list_for_(status, k);
        }
        catch (const std::exception &e)
        {
            throw std::invalid_argument(
                std::string(e.what()) + " for status " + std::to_string(status) + ".");
        }
    }

    std::vector<size_t> Augmented::get_unsorted_definition_rows_(int status) const
    {
        // load the list of status-specific definition rows
        std::vector<size_t> allrows = map_status_to_def_row_.at(status);

        // extract names
        std::vector<std::string> name_to_status = {};
        for (size_t k : allrows)
            name_to_status.push_back(deftable_.find("name")->second[k]);

        // add generic definitions to current list (if not generic already)
        // skip generic if a status-specific already exist
        if (status::is_infected(status) &&
            map_status_to_def_row_.count(status::MISSING))
            for (size_t k : map_status_to_def_row_.at(status::MISSING))
                if (std::find(name_to_status.begin(), name_to_status.end(), deftable_.find("name")->second[k]) == name_to_status.end())
                    allrows.push_back(k);

        return allrows;
    }

    void Augmented::set_sample_list_for_(int status, size_t k)
    {
        const std::string &augmstr = deftable_["name"][k];
        size_t a = pRegistry->from_string(augmstr);

        bool is_to_export = true;
        bool is_to_draw = true;

        if (str::tobool(deftable_["const"][k]))
        {
            is_to_draw = false;
        }
        else if (!na::is(deftable_["func"][k]))
        {
            const parse::Definition parser(deftable_["func"][k]);
            if (parser.is_terminal())
            {
                // check the type of requested equality
                if (parse::is_parameter(parser.root.type))
                    is_to_draw = !pParameter_->at(pParameter_->pRegistry->from_string(parser.root.value)).is_constant();
                // UPDATE: can I use dependecies tree to see if constant function?
            }
            else if (!parser.is_callable_implicitly("x"))
            {
                is_to_export = false;
            }
        }
        else if (na::is(deftable_["value"][k]) && na::is(deftable_["model"][k]))
        {
            throw std::invalid_argument(
                "Impossible to initialize augmented '" + augmstr + "'  without a value or model");
        }
        else if (na::is(deftable_["proposal"][k]))
        {
            throw std::runtime_error(
                "No proposal was provided for non-constant definition of augmented '" +
                augmstr + "'");
        }
        else
        {
            is_any_sample_ = true;
            status_blocks_[map_status_to_block_.at(status)].sample.emplace_back(a);
        }

        if (!str::is_private(augmstr) && is_to_export)
        {
            if (std::find(list_export.begin(), list_export.end(), a) == list_export.end())
                list_export.emplace_back(a);

            if (is_to_draw)
                status_blocks_[map_status_to_block_.at(status)].draws.emplace_back(a);
        }
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////  Augmented : initialize augmented  //////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////

    void Augmented::init()
    {
        if (pRegistry->size() == 0)
            return;

        helper_init_model_list_();

        for (Household &house : *pDataWide_)
            for (size_t i : house.infected_hhindexes())
                init_missing(house.member(i), true);
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////  Augmented : init augmented models  /////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////

    void Augmented::helper_init_model_list_()
    {
        for (StatusBlock_ &block : status_blocks_)
        {
            block.depen.resize(pRegistry->size());

            // associate blocks for ordered model lists
            const str_table &ordered_models = get_ordered_models_(block);
            if (!ordered_models.empty())
            {
                is_any_model_ = true;
                block.models = pRegistry->from_string(read::extract_keys(ordered_models));

                for (const auto &[augmstr, others] : ordered_models)
                    block.depen[pRegistry->from_string(augmstr)] = pRegistry->from_string(others);
            }
        }
    }

    str_table Augmented::get_ordered_models_(const StatusBlock_ &block) const
    {
        // columns presenting definitions possibly containing sources
        const auto columns = {"value", "model", "func", "min", "max"};

        str_table tree;                         // populate augmented dependency tree
        str_table ordered_models;               // store hierarchical order of models
        std::vector<std::string> without_model; // list of augmented without models
        for (size_t k : block.order)
        {
            const std::string &aumgstr = deftable_.at("name")[k];
            std::vector<std::string> xdepen = {};
            for (auto c : columns)
                xdepen += parse::Definition(deftable_.at(c)[k]).list_augmented_depen();

            if (na::is(deftable_.at("model")[k]))
                without_model.push_back(aumgstr);
            else
                tree[aumgstr] = str::unique(xdepen);
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
    /////////////////////////////////  Augmented : init missing attribute  ////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////

    void Augmented::init_missing(Individual &member, bool force)
    {
        // loop over augmented variable definitions to assign placheholders
        for (size_t k : helper_materialize_definition_order_(member))
            if (force || member.augm_vec_[pRegistry->from_string(deftable_["name"][k])].is_empty())
                helper_init_attributes_(member, k);
    }

    std::vector<size_t> Augmented::helper_materialize_definition_order_(Individual &member)
    {
        int status = member.status();

        // if known status, then load its definition
        const auto it = map_status_to_block_.find(status);
        if (it != map_status_to_block_.end())
        {
            const size_t i = it->second;
            member.augm_idx_ = i;
            return status_blocks_[i].order;
        }

        // otherwise if it exists fallback to generic definition
        const auto generic_it = map_status_to_block_.find(status::MISSING);
        if (generic_it != map_status_to_block_.end())
        {
            const size_t i = generic_it->second;
            map_status_to_block_.emplace(status, i);
            member.augm_idx_ = i;
            return status_blocks_[i].order;
        }

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

            if (str::tobool(deftable_["const"][k]))
            {
                xaugm.init_numeric(
                    pParameter_->read_numeric(deftable_["value"][k], &member).value());
                xaugm.enforce_user_const(); // enforce user constant choice
            }
            else
                helper_init_func_attributes_(member, k);

            // initialize model
            helper_init_model_attributes_(member, k);

            // assign proposal if variable to be sampled
            // NOTE: not necessarily according to a model
            if (!xaugm.is_constant())
                xaugm.init_proposal(deftable_["proposal"][k],
                                    str::tod(deftable_["proposal_scale"][k]));

            // check if satisfies discrete requirement, if any
            if (user_discrete_set_.count(a) && !xaugm.is_discrete())
                throw std::invalid_argument(
                    "Definition of " + xaugm.err_name() +
                    " did not return a discrete object as requested by model.");
        }
        catch (const std::exception &e)
        {
            throw std::runtime_error(
                "While initializing augmented '" + augmstr + "' for " +
                member.err_name() + " (status " + std::to_string(member.status()) + "):\n" +
                std::string(e.what()));
        }
    }

    void Augmented::helper_init_model_attributes_(Individual &member, size_t k)
    {
        const std::string &augmstr = deftable_["name"][k];
        Variable &xaugm = member.augm_vec_[pRegistry->from_string(augmstr)];
        parse::Definition parser(deftable_["model"][k]);
        xaugm.init_distr(parser.root.func);

        if (!xaugm.is_model())
            return;

        for (const auto &child : parser.root.children)
            xaugm.init_distr_attr(child.key, pParameter_->read_attribute(child, &member));

        for (const auto key : {"min", "max"})
            xaugm.init_distr_attr(key, pParameter_->read_attribute(deftable_[key][k], &member));

        xaugm.init_distr_reshaping(str::tobool(deftable_["discr"][k]));
    }

    void Augmented::helper_init_func_attributes_(Individual &member, size_t k)
    {
        const std::string &augmstr = deftable_["name"][k];
        Variable &xaugm = member.augm_vec_[pRegistry->from_string(augmstr)];
        parse::Definition parser(deftable_["func"][k]);

        // pointer mode
        if (parser.is_terminal())
        {
            xaugm.transfer(pParameter_->read_attribute(parser.root, &member));
            return;
        }

        // function mode
        bool is_primitive = str::tobool(deftable_["fcumul"][k]);
        xaugm.init_function(compose_func_registries(
            parser.root.func, "function",
            {"custom", "transcendental", "arithmetic",
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
        xaugm.init_func_reshaping(str::tobool(deftable_["discr"][k]));
    }

    /////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////  Augmented : user reset  ////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////

    void Augmented::user_reset(rnd_gen &prng, AugmContainer &member)
    {
        for (size_t a : status_blocks_[member.augm_idx_].sample)
            member.augm_vec_[a].user_reset(prng);

        for (size_t a = 0; a < pRegistry->size(); ++a)
            member.augm_vec_[a].optimize();
    }

    void Augmented::user_reset(rnd_gen &prng)
    {
        if (pRegistry->size() == 0)
            return;

        for (Household &house : *pDataWide_)
            for (size_t i : house.infected_hhindexes())
                user_reset(prng, house.member(i));
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////  Augmented : log-likelihood  //////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////

    const std::vector<size_t> &Augmented::list_sample(const AugmContainer &member) const
    {
        return status_blocks_[member.augm_idx_].sample;
    }
    const std::vector<size_t> &Augmented::list_draws(const AugmContainer &member) const
    {
        return status_blocks_[member.augm_idx_].draws;
    }

    double Augmented::commit_augm_model_log_lik(AugmContainer &member)
    {
        const StatusBlock_ &block = status_blocks_[member.augm_idx_];
        if (block.models.empty())
            return 0.;
        double new_ll = 0.0;
        for (size_t a : block.models)
            new_ll += member.augm_vec_[a].get_distr_log_pdf();
        member.augm_log_lik_.set(new_ll);
        return new_ll;
    }

    double Augmented::commit_augm_model_log_ratio(AugmContainer &member, size_t a)
    {
        const StatusBlock_ &block = status_blocks_[member.augm_idx_];
        if (block.depen[a].empty())
            return 0.;
        double delta_ll = 0;
        // NOTE: vector of dependencies 'depen' includes also 'a'
        for (size_t b : block.depen[a])
            delta_ll += member.augm_vec_[b].get_distr_log_ratio();
        member.augm_log_lik_.increment(delta_ll);
        return delta_ll;
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////  Augmented : sample draws  ///////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////

    void Augmented::init_individual_sample(std::vector<std::string> &id_draws)
    {
        smpl_individuals.clear();

        if (id_draws.empty())
            return;

        if (!is_any_sample())
        {
            msg::warning(
                "Ignored id draws. No augmented variables to sample were defined.");
            id_draws.clear();
            return;
        }

        smpl_individuals.reserve(id_draws.size());
        std::vector<std::string> existing_ids = {};
        for (const Household &house : *pDataWide_)
            for (const Individual &member : house)
            {
                if (std::find(id_draws.begin(), id_draws.end(), member.id()) == id_draws.end())
                    continue;

                // unique IDs are granted by class DataWide load
                existing_ids.push_back(member.id());

                // check status of chosen individual
                if (member.is_infected_status())
                    smpl_individuals.emplace_back(&member);
                else
                    msg::warning(
                        "Ignored draws of non infected " + member.err_name() + ".");
            }
        smpl_individuals.shrink_to_fit();

        // check all user request are satisfied
        str::assert_no_missing_entries(existing_ids, id_draws, "draws id");
    }

    void Augmented::init_draws(unsigned long sample_size)
    {
        if (!is_any_draw())
            return;

        smpl_values.clear();
        smpl_values.resize(smpl_individuals.size());
        for (size_t i = 0; i < smpl_individuals.size(); ++i)
        {
            auto &member_matrix = smpl_values[i];
            size_t n = list_draws(*smpl_individuals[i]).size();
            member_matrix.reserve(n);
            for (size_t j = 0; j < n; ++j)
                member_matrix.emplace_back(sample_size, na::DOUBLE);
        }
    }

    void Augmented::store_draws(unsigned long smpl_count)
    {
        if (!is_any_draw())
            return;

        for (size_t i = 0; i < smpl_individuals.size(); ++i)
        {
            const AugmContainer &member_ptr = *smpl_individuals[i];
            auto &member_matrix = smpl_values[i];

            size_t c = 0;
            for (size_t a : list_draws(member_ptr))
            {
                const Variable &xaugm = member_ptr.augm_vec_[a];
                member_matrix[c++][smpl_count] = xaugm.value();
            }
        }
    }
}
