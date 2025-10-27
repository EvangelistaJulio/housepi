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
#include "variable.h"
#include "household.h"
#include "distribution_registry.h"
#include "mcmc_impl.h"

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////  AugmentedMap : setter  ////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/// A string containing a summary of the augmented map.
std::string AugmentedMap::str_augm_summary() const
{
    std::string output;
    for (const auto &pair : augm_map_)
    {
        output += "- " + pair.first + " : ";
        if (pair.second.is_function() && !pair.second.get_function()->is_implicit_callable())
            output += pair.second.get_function()->str_summary();
        else
            output += pair.second.str_summary();
        output += "\n";
    }
    return output;
}

/// Associate functional methods for augmented model pdfs.
void AugmentedMap::init_augm_model_logratios(const u_dict<str_vec> &models_ordered)
{
    // >> clean pre-existing
    for (auto &pair : model_logratios_)
    {
        pair.second.clear();
        pair.second.shrink_to_fit();
    }
    model_logratios_.clear();
    model_logratios_.rehash(0);
    // <<

    // associate models if applicable
    for (const auto &pair : models_ordered)
    {
        const std::string &a = pair.first;
        model_logratios_[a] = {};

        // repeat for all dependent models
        for (const std::string &b : pair.second)
        {
            const Variable &xaugm = augm_map_[b];
            model_logratios_[a].emplace_back([&xaugm]()
                                             { return xaugm.get_distr_logratio(); });
        }
    }
}

/// Value of log-likelihood associated to augmented variable models.
/// Store the log-likelihood of individual augmented variables.
double AugmentedMap::compute_and_set_augm_model_log_lik()
{
    buffer_ = augm_log_lik_; // store old value in buffer
    augm_log_lik_ = 0.0;
    for (const auto &pair : model_logratios_)
        augm_log_lik_ += augm_map_[pair.first].get_distr_logpdf();
    assert_log_lik_not_NaN();
    return augm_log_lik_;
}

/// Change in log-likelihood associated to given augmented.
double AugmentedMap::compute_and_set_augm_model_logratio(const std::string &a)
{
    buffer_ = augm_log_lik_; // store old value in buffer
    double delta_ll = 0.;
    for (const std::function<double()> &f : model_logratios_[a])
        delta_ll += f();
    augm_log_lik_ += delta_ll;
    assert_log_lik_not_NaN(a);
    return delta_ll;
}

/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////  MCMC : Impl : InfecProfile : assign template  ///////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

/// Set model template while constructiing model.
void MCMC::Impl::InfecProfile::init_template(u_dict<std::string> deflist)
{
    // if not previous build, proceed with template
    // if no template, assume no need for an infectivity profile
    if (!hsp::na::is(profile_->name()) ||
        !augm_infec_profile_deflist_.empty() ||
        deflist.empty())
        return;

    try
    {
        if (!deflist.count("func"))
            throw std::invalid_argument(
                "A function name was missing for infectivity profile template.");

        const str_vec file_headers = {"func", "discr", hsp::params::MAX};

        for (const std::string &key : file_headers)
            if (!deflist.count(key))
                deflist.emplace(key, hsp::na::STRING);

        // discretization
        deflist["discr"] = "false";
        if (is_discrete_requested_)
            deflist["discr"] = "true";
        else if (deflist.count("discr") && !hsp::na::is(deflist.at("discr")))
            deflist["discr"] = deflist.at("discr");

        build(deflist);

        // forbid augmented infectivity profile for template
        if (is_augmented())
            throw std::invalid_argument(
                "Tried to define an infectivity profile template depending on augmented variables.");

        // UPDATE: maybe also not depending on parameters?
    }
    catch (const std::exception &e)
    {
        throw std::runtime_error(
            "While building infectivity profile template:\n" + std::string(e.what()));
    }
}
void MCMC::init_infec_profile_template(u_dict<std::string> x) { pImpl_->infec_profile.init_template(x); }

/// Store values of the infectivity profile.
void MCMC::Impl::InfecProfile::init_lookup_()
{
    inst_ = &InfecProfile::do_inst_lookup_;
    cumul_ = &InfecProfile::do_cumul_lookup_;

    if (hsp::na::is(profile_->get_max()))
        size_lut_ = hsp::infprof::LOOKUP_SIZE;
    else
        size_lut_ = profile_->get_max() + 1;

    lookup_inst_.resize(size_lut_);
    lookup_cumul_.resize(size_lut_);
    lookup_buffer_inst_.resize(size_lut_);
    lookup_buffer_cumul_.resize(size_lut_);

    update(true);
}

/// Clear infectivity profile lookup table.
void MCMC::Impl::InfecProfile::clear_lookup_()
{
    size_lut_ = 0;

    lookup_inst_.clear();
    lookup_inst_.shrink_to_fit();
    lookup_cumul_.clear();
    lookup_cumul_.shrink_to_fit();
    lookup_buffer_inst_.clear();
    lookup_buffer_inst_.shrink_to_fit();
    lookup_buffer_cumul_.clear();
    lookup_buffer_cumul_.shrink_to_fit();
}

/// If infectivity profile is a lookup table, store its values.
/// UPDATE: this should not be recoverd if constant.
void MCMC::Impl::InfecProfile::update(bool is_cache_depen)
{
    if (!is_cache_depen || !is_lookup())
        return;

    // move current lookup table to buffer
    lookup_inst_.swap(lookup_buffer_inst_);
    lookup_cumul_.swap(lookup_buffer_cumul_);

    // Set infectivity to 0 on the same day
    lookup_inst_[0] = 0.;
    lookup_cumul_[0] = 0.;

    for (size_t i = 1; i < size_lut_; ++i)
    {
        double x = do_inst_default_(nullptr, i);
        lookup_inst_[i] = x;
        lookup_cumul_[i] = x + lookup_cumul_[i - 1];
    }
}

/// If infectivity profile is a lookup table, recover its buffered values.
/// UPDATE: this should not be recoverd if constant.
void MCMC::Impl::InfecProfile::recover(bool is_cache_depen)
{
    if (!is_cache_depen && !is_lookup())
        return;

    lookup_inst_.swap(lookup_buffer_inst_);
    lookup_cumul_.swap(lookup_buffer_cumul_);
}

/// Build infectivity profile distribution .
void MCMC::Impl::InfecProfile::build(const u_dict<std::string> &deflist)
{
    augm_infec_profile_deflist_.clear();

    if (deflist.empty())
        return;

    const str_vec file_headers = {"func", "discr", hsp::params::MAX};

    // check entries
    hsp::str::assert_no_unused_entries(hsp::read::list_headers(deflist), file_headers);

    // check definition list has a function
    if (!deflist.count("func"))
        throw std::invalid_argument(
            "A function name choice was missing in the infectivity profile definition.");
    else if (DefinitionParser(deflist.at("func")).is_na())
        throw std::invalid_argument(
            "Provided NA as a function name in the infectivity profile definition.");

    // check if definition is generic or augmented
    bool is_augm = DefinitionParser(deflist.at("func")).is_augmented_depen() ||
                   DefinitionParser(deflist.at(hsp::params::MAX)).is_augmented_depen();

    if (is_augm)
        build_augmented_profile_(deflist);
    else
        build_generic_profile_(deflist);
}

/// Build for generic infectivity profile.
void MCMC::Impl::InfecProfile::build_generic_profile_(const u_dict<std::string> &deflist)
{
    inst_ = &InfecProfile::do_inst_default_;
    cumul_ = &InfecProfile::do_cumul_default_;

    DefinitionParser parser(DefinitionParser(deflist.at("func")));
    profile_ = compose_distr_registries(parser.root.func,
                                        "infectivity profile",
                                        {&InfecProfileRegistry::inst(),
                                         &ContinuosDistrRegistry::inst(),
                                         &DiscreteDistrRegistry::inst()});
    for (const auto &node : parser.root.children)
    {
        profile_->init_attr(node.key, pImpl_.read_attribute_(node));
        pImpl_.make_param_cache_depen_(node);
    }

    // discretization
    if (deflist.count("discr") && hsp::str::tobool(deflist.at("discr")))
        profile_->init_discretization();
    else if (is_discrete_requested_)
        hsp::msg::error(
            "Definition of infectivity profile did not return a discrete object as requested by the model.");

    // range minimum
    if (profile_->is_discrete())
        profile_->init_attr(hsp::params::MIN, pImpl_.read_attribute_("1|N"));
    else
        profile_->init_attr(hsp::params::MIN, pImpl_.read_attribute_("0|N"));

    // range maximum
    if (deflist.count(hsp::params::MAX))
        profile_->init_attr(hsp::params::MAX, pImpl_.read_attribute_(deflist.at(hsp::params::MAX)));
    // range truncation
    profile_->init_truncation();

    // if discrete and generic, initializes lookup table
    if (profile_->is_discrete())
        init_lookup_();
    else
        clear_lookup_();
}

/// Build for individual augmented function infectivity profile.
void MCMC::Impl::InfecProfile::build_augmented_profile_(const u_dict<std::string> &deflist)
{
    // use generic build to check profile name at this stage
    // NOTE: this also control for discretization
    // UPDATE: not very clean, maybe refactor
    build_generic_profile_({{"func", DefinitionParser(deflist.at("func")).root.func},
                            {"discr", deflist.at("discr")}});

    // drop generic definition
    clear_lookup_();
    profile_ = std::make_unique<Distribution>();
    inst_ = &InfecProfile::do_inst_augm_;
    cumul_ = &InfecProfile::do_cumul_augm_;

    // initialize for instantaneous and cumulative profiles
    bool is_discrete = hsp::str::tobool(deflist.at("discr"));

    // match infectivity profile and augmented variable definition
    augm_infec_profile_deflist_ = {
        {"name", str_vec({hsp::infprof::INSTANTANEOUS, hsp::infprof::CUMULATIVE})},
        {"status", str_vec(2, hsp::na::STRING)}, // UPDATE
        {"func", str_vec(2, deflist.at("func"))},
        {"fcumul", {"false", "true"}},
        {"value", str_vec(2, hsp::na::STRING)},
        {"const", str_vec(2, "false")},
        {"discr", str_vec(2, hsp::str::format(is_discrete))},
        {"model", str_vec(2, hsp::na::STRING)},
        {hsp::params::MIN, str_vec(2, (is_discrete) ? "1|N" : "0|N")},
        {hsp::params::MAX, str_vec(2, deflist.at(hsp::params::MAX))},
        {"rule", str_vec(2, hsp::na::STRING)},
        {"scale", str_vec(2, hsp::na::STRING)}};
}

// Eventually append infectivity profile definitions to augmented.
void MCMC::Impl::InfecProfile::transfer_augm_definitions(u_dict<str_vec> &augm_deftable) const
{
    if (!is_augmented())
        return;

    for (const auto &pair : augm_infec_profile_deflist_)
        augm_deftable.at(pair.first) += pair.second;

    if (is_discrete_requested_)
        for (const std::string &a : {hsp::infprof::INSTANTANEOUS, hsp::infprof::CUMULATIVE})
            pImpl_.augm.user_list_discrete.emplace_back(a);
}

/// Request discrete infectivity profile.
void MCMC::request_infec_profile_discrete() { pImpl_->infec_profile.request_discrete(); }

/////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////  MCMC : Impl : InfecProfile : getter methods  ///////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

/// Lookup method for discrete instantaneous infectivity profile.
double MCMC::Impl::InfecProfile::do_inst_lookup_(const Individual *, double t)
{
    if (profile_->is_valid(t) && t < size_lut_)
        return lookup_inst_.at(static_cast<size_t>(t));
    else
        return 0;
}
/// Lookup method for discrete cumulative infectivity profile.
double MCMC::Impl::InfecProfile::do_cumul_lookup_(const Individual *, double t)
{
    if (profile_->is_below(t))
        return 0.;
    else if ((t >= size_lut_) || (profile_->is_above(t)))
        return lookup_cumul_.at(size_lut_ - 1);
    else
        return lookup_cumul_.at(static_cast<size_t>(t));
}

double MCMC::infec_profile_inst(const Individual &i, double t) const { return pImpl_->infec_profile.inst(&i, t); }
double MCMC::infec_profile_cumul(const Individual &i, double t) const { return pImpl_->infec_profile.cumul(&i, t); }

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
    if (!tree.count(entry))
        throw std::invalid_argument(
            "Could not find a definition for augmented '" + entry + "'.");
    if (recursion_stack % entry)
        throw std::invalid_argument(
            "Found inadmissible nested definitions for augmented '" + entry + "'");
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
        if (!reverse_graph.count(node))
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
const Attribute &MCMC::Impl::read_attribute_(const DefinitionParser::ASTNode &node,
                                             const Individual &member) const
{
    storage_ = Attribute(); // reset storage to default empty

    if (!node.is_terminal()) // function node
    {
        Attribute temporary(node.func);
        for (const auto &child : node.children)
            temporary.get_function()->init_attr(child.key, read_attribute_(child, member));
        storage_ = std::move(temporary);
    }
    else if (hsp::sources::is_data(node.type)) // terminal data node
    {
        storage_ = Attribute(member.data(node.value));
        storage_.set_name(node.value, "data");
    }
    else if (hsp::sources::is_numeric(node.type)) // terminal numeric node
        storage_ = Attribute(hsp::str::tod(node.value, "argument with numeric type"));

    else if (hsp::sources::is_augmented(node.type)) // augmented terminal node
        return member.get_augm_variable(node.value);

    else if (hsp::sources::is_group_rate(node.type)) // group rate terminal node
        return member.get_group_rate_attribute(node.value);

    else if (hsp::sources::is_parameter(node.type)) // terminal parameter node
        return param.at(node.value);

    else if (!node.is_empty())
        throw std::invalid_argument(
            "Unrecognized source type '" + node.type + "' found while parsing '" + node.str() + "'.");

    return storage_;
}
const Attribute &MCMC::Impl::read_attribute_(const std::string &expr, const Individual &member) const
{
    return read_attribute_(DefinitionParser(expr).root, member);
}

/// Convert a typed string into a double.
const Attribute &MCMC::Impl::read_numeric_(const DefinitionParser::ASTNode &node, const Individual &member) const
{
    const Attribute &aux_attr_ = read_attribute_(node, member);
    if (!aux_attr_.is_empty() && !aux_attr_.is_constant())
        throw std::invalid_argument(
            "Definition '" + node.str() + "' could not be interpreted as a numeric constant.");
    return aux_attr_;
}
const Attribute &MCMC::Impl::read_numeric_(const std::string &expr, const Individual &member) const
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
void MCMC::Impl::Augmented::build(const u_dict<str_vec> &deftable)
{
    this->deftable_ = deftable;
    helper_assert_no_reserved_();

    // eventually append infectivity profile definitions
    pImpl_.infec_profile.transfer_augm_definitions(this->deftable_);

    if (deftable_.empty())
        return;

    // check entries
    const str_vec file_headers = {"name", "status",
                                  "func", "fcumul",
                                  "value", "const", "discr",
                                  "model", hsp::params::MIN, hsp::params::MAX,
                                  "rule", "scale"};
    const str_vec headers = hsp::read::table_headers(deftable_);
    hsp::str::assert_no_missing_entries(headers, file_headers);
    hsp::str::assert_no_unused_entries(headers, file_headers);

    helper_build_definition_();
    helper_build_definition_order_();
    helper_build_model_status_maps_();

    // allocate memory for augmented variables to each individual
    for (Household &house : pImpl_.data_wide)
        for (Individual &member : house)
            clear_all(member);
}

// check augmented name
void MCMC::Impl::Augmented::helper_assert_no_reserved_()
{
    if (deftable_.empty() || !deftable_.count("name"))
        return;

    for (const std::string &a : deftable_["name"])
        if (hsp::na::is(a))
            throw std::invalid_argument(
                "Found an augmented variable with undefined name.");
        else if (a == hsp::infprof::INSTANTANEOUS || a == hsp::infprof::CUMULATIVE)
            throw std::invalid_argument(
                "Tried to use the reserved augmented name '" + a + "'.\nPlease, choose a different name.");
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

    // check that status-specific definition is for infected status
    for (size_t k = 0; k < deftable_["name"].size(); ++k)
    {
        std::string s = deftable_["status"][k];
        if (hsp::na::is(s))
            s = hsp::na::STRING; // enforce internal NA convention

        else if (!hsp::status::is_infected(hsp::str::toi(s, "augmented status")))
            throw std::invalid_argument(
                "Impossible to associate augmented '" + deftable_["name"][k] + "' for non infectious status " + s + ".");

        // store status to definition row index
        status_to_def_row[s].push_back(k);
    }

    // assert that name and status are univocal
    for (const auto &pair : status_to_def_row)
    {
        str_vec name_to_status = {};
        for (size_t k : pair.second)
            name_to_status.push_back(deftable_["name"][k]);

        if (name_to_status.size() != hsp::str::unique(name_to_status).size())
            throw std::invalid_argument(
                "Incompatible definitions found for augmented '" + pair.first + "'.");
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
    if (pImpl_.data_wide.all_infected_status_list.empty())
        set_definition_dependency_order_for_(hsp::na::STRING);
    else
        for (int status : pImpl_.data_wide.all_infected_status_list)
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
            xdepen += DefinitionParser(deftable_[c][k]).list_augmented_depen();
        tree[deftable_["name"][k]] = hsp::str::unique(xdepen);
        name_to_row[deftable_["name"][k]] = k;
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
    if (status_to_def_row.count(s))
    {
        std::vector<size_t> allrows = status_to_def_row.at(s);

        // extract names
        str_vec name_to_status = {};
        for (size_t k : allrows)
            name_to_status.push_back(deftable_.find("name")->second[k]);

        // add generic definitions to current list (if not generic already)
        if (!hsp::na::is(s) && status_to_def_row.count(hsp::na::STRING))
            for (size_t k : status_to_def_row.at(hsp::na::STRING))
                // skip generic if status-specific already exist
                if (!(name_to_status % deftable_.find("name")->second[k]))
                    allrows.push_back(k);

        return allrows;
    }

    else if (status_to_def_row.count(hsp::na::STRING))
        // use generic definitions if status-specific don't exist
        return status_to_def_row.at(hsp::na::STRING);

    else
        throw std::invalid_argument(
            "No definitions were compatible with augmented status " + s + ".");
}

/// Set sample map and relative usage list.
void MCMC::Impl::Augmented::set_sample_map_for_(const std::string &s, size_t k)
{
    const std::string &a = deftable_["name"][k];

    if (hsp::str::tobool(deftable_["const"][k]))
        list_export.push_back(a);

    else if (!hsp::na::is(deftable_["func"][k]))
        set_function_or_identity_for_(k);

    else if (hsp::na::is(deftable_["value"][k]) && hsp::na::is(deftable_["model"][k]))
        throw std::invalid_argument(
            "Impossible to initialize augmented '" + a + "' for status " + s + " without a value or model.");

    else if (!hsp::na::is(deftable_["rule"][k]))
    {
        // read option for monte carlo proposer
        list_draws.push_back(a);
        if (!hsp::na::is(s))
            status_map_sample[hsp::str::toi(s)].push_back(a);
    }
    else
        throw std::runtime_error(
            "No rule was provided for non-constant definition of augmented '" + a +
            "' for status " + s + ".");
}

///////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////  MCMC : Impl : Augmented : build functions  /////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/// Build augmented function and inspect for constant.
void MCMC::Impl::Augmented::set_function_or_identity_for_(size_t k)
{
    const std::string &a = deftable_["name"][k];
    const DefinitionParser parser(deftable_["func"][k]);

    if (parser.is_terminal())
    {
        bool is_identity_const = true;

        // check the type of requested equality
        if (hsp::sources::is_parameter(parser.root.type))
            is_identity_const = pImpl_.param.at(parser.root.value).is_constant();

        else if (hsp::sources::is_augmented(parser.root.type))
            is_identity_const = false; // UPDATE: use an augmented dependecies tree to see if constant

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

    for (int status : pImpl_.data_wide.all_infected_status_list)
    {
        u_dict<str_vec> tree;       // populate augmented dependency tree
        u_dict<size_t> name_to_row; // ... and reverse dictionary name to index
        str_vec without_model;      // list of augmented without models
        for (const size_t k : get_unsorted_definition_rows_by_(std::to_string(status)))
        {
            const std::string &a = deftable_["name"][k];
            str_vec xdepen = {};
            for (const std::string &c : str_vec({"value", "model", "func", hsp::params::MIN, hsp::params::MAX}))
                xdepen += DefinitionParser(deftable_[c][k]).list_augmented_depen();

            if (hsp::na::is(deftable_["model"][k]))
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
void MCMC::clear_all_augm(Individual &member) { pImpl_->augm.clear_all(member); }

/// Initialize augmented dependecies for evaluations.
/// NOTE: this must not be called before simulation which build on demand.
void MCMC::Impl::Augmented::init()
{
    helper_init_id_draws_();

    if (!is_any())
        return;

    for (Household &house : pImpl_.data_wide)
        for (Individual &member : house.infected())
        {
            init_all(member);
            init_model_logratios(member);
        }
}

/// Assign all augmented variables to single individual.
void MCMC::Impl::Augmented::init_all(Individual &member)
{
    // loop over augmented variable definitions to assign placheholders
    for (size_t k : init_get_definition_order_(member))
        helper_init_attributes_(member, k);
}

/// Assign augmented variables only if missing (used for simulations).
void MCMC::Impl::Augmented::init_missing(Individual &member)
{
    // loop over augmented variable definitions to assign placheholders
    for (size_t k : init_get_definition_order_(member))
        if (member.get_augm_variable(deftable_["name"][k]).is_empty())
            helper_init_attributes_(member, k);
}
void MCMC::init_missing_augm(Individual &member)
{
    pImpl_->augm.init_missing(member);
    pImpl_->augm.user_reset(member);
}

// Assign all augmented models to single individual.
void MCMC::Impl::Augmented::init_model_logratios(Individual &member)
{
    try
    {
        member.init_augm_model_logratios(status_map_model[member.status()]);
    }
    catch (const std::exception &e)
    {
        throw std::runtime_error(
            "While initializing " + section_name_ + " model " +
            " for " + member.err_name() + ":\n" +
            std::string(e.what()));
    }
}

/// Get the definition order for current member according to its status.
std::vector<size_t> MCMC::Impl::Augmented::init_get_definition_order_(Individual &member)
{
    const std::string s = std::to_string(member.status());

    // NOTE: if status hasn't been observed in data, build on demand
    if (member.is_infected_status() && !(pImpl_.data_wide.all_infected_status_list % member.status()))
    {
        set_definition_dependency_order_for_(s);
        pImpl_.data_wide.all_infected_status_list.push_back(member.status());
    }

    if (status_to_def_row.count(s))
        return status_to_def_row[s];
    else if (status_to_def_row.count(hsp::na::STRING))
        return status_to_def_row[hsp::na::STRING];
    else if (member.is_infected_status())
        throw std::invalid_argument(
            "No augmented definitions were found for " + member.err_name() + " with status " + s + ".");
    else
        throw std::invalid_argument(
            "Generic augmented definition for non-infectious " + member.err_name() + " was missing.");
}

/// Assign model arguments and range. Same for function, otherwise assign proposer.
void MCMC::Impl::Augmented::helper_init_attributes_(Individual &member, size_t k)
{
    const std::string &a = deftable_["name"][k];

    try
    {
        member.init_augm(a, pImpl_.read_numeric_(deftable_["value"][k], member).get(),
                         hsp::str::tobool(deftable_["const"][k]),
                         deftable_["func"][k], hsp::str::tobool(deftable_["fcumul"][k]),
                         deftable_["model"][k]);
        Variable &xaugm = member.get_augm_variable(a);

        // assign arguments to model and associate discretization
        init_model_attributes_(member, k);
        if (xaugm.is_function()) // function
            init_func_attributes_(member, k);
        else if (xaugm.is_pointer()) // pointer (case of identity function)
            xaugm.transfer(pImpl_.read_attribute_(deftable_["func"][k], member));

        // assign proposer if variable to be sampled
        if (!xaugm.is_constant())
            xaugm.init_proposer(deftable_["rule"][k], hsp::str::tod(deftable_["scale"][k]));

        // check if satisfies discrete requirement, if any
        if (!user_list_discrete.empty() && (user_list_discrete % a) && !xaugm.is_discrete())
            throw std::invalid_argument(
                "Definition of " + xaugm.err_name() +
                " did not return a discrete object as requested by the model.");
    }
    catch (const std::exception &e)
    {
        throw std::runtime_error(
            "While initializing " + section_name_ + " '" + a +
            "' for " + member.err_name() + ":\n" +
            std::string(e.what()));
    }
}
void MCMC::request_augm_discrete(const std::string &x) { pImpl_->augm.user_list_discrete.emplace_back(x); }

/// Auxiliary method to assign attributes to augmented model.
void MCMC::Impl::Augmented::init_model_attributes_(Individual &member, size_t k)
{
    Variable &xaugm = member.get_augm_variable(deftable_["name"][k]);

    for (const auto &child : DefinitionParser(deftable_["model"][k]).root.children)
        xaugm.init_distr_attr(child.key, pImpl_.read_attribute_(child, member));

    for (const std::string &key : {hsp::params::MIN, hsp::params::MAX})
        xaugm.init_distr_attr(key, pImpl_.read_attribute_(deftable_[key][k], member));

    xaugm.init_distr_reshaping(hsp::str::tobool(deftable_["discr"][k]));
}

/// Auxiliary method to assign attributes to augmented function.
void MCMC::Impl::Augmented::init_func_attributes_(Individual &member, size_t k)
{
    Variable &xaugm = member.get_augm_variable(deftable_["name"][k]);

    // initialize function
    for (const auto &child : DefinitionParser(deftable_["func"][k]).root.children)
    {
        xaugm.init_func_attr(child.key, pImpl_.read_attribute_(child, member));
        pImpl_.make_param_cache_depen_(child);
    }

    // initialize range
    for (const std::string &key : {hsp::params::MIN, hsp::params::MAX})
        xaugm.init_func_attr(key, pImpl_.read_attribute_(deftable_[key][k], member));

    // initialize discretization and truncation
    xaugm.init_func_reshaping(hsp::str::tobool(deftable_["discr"][k]));
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
    if (!is_any())
        return;

    for (Household &house : pImpl_.data_wide)
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

    if (id_draws.empty())
        return;

    if (is_any_sample())
    {
        smpl_individuals.reserve(id_draws.size());

        // check id across all individuals
        str_vec already_stored_id = {};
        for (const Household &house : pImpl_.data_wide)
            for (const Individual &member : house)
                if (id_draws % member.id())
                {
                    if (already_stored_id % member.id())
                        throw std::invalid_argument(
                            "Multiple individuals shared the same id '" + member.id() + "'.");
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
                    hsp::msg::warning(
                        "Ignored draws of non infected " + member->err_name() + ".");
        }
        else
        {
            // UPDATE: is this even possible?
            hsp::msg::warning(
                "Ignored id draws for individual(s) " + hsp::str::format(id_draws));
            id_draws.clear();
        }
    }
    else
    {
        hsp::msg::warning(
            "Ignored id draws. No augmented variables to sample were defined.");
        id_draws.clear();
    }
}

/// Initialize sample for individual augmented variable followup.
void MCMC::Impl::Augmented::init_draws(unsigned long sample_size)
{
    if (!is_any_draw())
        return;

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

/// Store sample for individual augmented variable followup.
void MCMC::Impl::Augmented::update_draws(unsigned long smpl_count)
{
    if (!is_any_draw())
        return;

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
    double delta_ll = 0.;
    std::atomic<bool> has_error(false);

#pragma omp parallel for reduction(+ : delta_ll)
    for (size_t H = 0; H < data_wide.report_house.total; ++H)
    {
        if (has_error.load(std::memory_order_relaxed))
            continue;

        Household &this_house = data_wide(H);
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

    increment_total_log_lik_(delta_ll);
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
        double delta_ll = -house.get_obs_log_lik();
        delta_ll += member.compute_and_set_augm_model_logratio(a);
        delta_ll += compute_and_set_household_obs_log_lik_(house);

        if (setup.unif() < std::exp(setup.inv_T * delta_ll - xaugm.get_propose_logratio()))
            return delta_ll;
        else
        {
            // else recover cached values
            house.recover_buffer_values(true);
            member.recover_augm_model_log_lik();
        }
    }
    // else recover buffered value
    xaugm.recover_value();

    return 0.;
}
