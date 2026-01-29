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
#include "tools/messages.h"
#include "tools/operators.h"
#include "tools/parser.h"
#include "tools/random.h"
#include "tools/read.h"
#include "tools/string.h"
#include "tools/types.h"

#include "attribute.h"
#include "augmented_map.h"
#include "data.h"
#include "distribution.h"
#include "function.h"
#include "group_rate.h"
#include "infec_profile.h"
#include "parameter.h"
#include "variable.h"

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////  *** AugmentedContainer : class ***  //////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/// A string containing a summary of the augmented map.
std::string AugmentedContainer::str_augm_summary() const
{
    std::string output;
    for (const Variable &var : augm_vec_)
    {
        output += "- " + std::string(var.name()) + " : ";
        if (var.is_function() && !var.get_function()->is_callable_implicitly())
            output += var.get_function()->str_summary();
        else
            output += var.str_summary();
        output += "\n";
    }
    return output;
}

/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////  *** Parameter : struct ***  ////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

Parameter::Parameter(const std::vector<std::string> &keys)
    : pRegistry(std::make_shared<hsp::key::Registry>(keys)) {}

/// Assign parameter definitions to dictionary.
void Parameter::build(const hsp::str_table &input)
{
    pRegistry->clear();

    if (input.empty())
        return;

    // check entries
    const std::vector<std::string> file_headers = {"name", "value", "const", "discr",
                                                   "prior", hsp::params::MIN, hsp::params::MAX,
                                                   "proposal", "proposal_scale"};
    const std::vector<std::string> headers = hsp::read::extract_keys(input);
    hsp::str::assert_no_missing_entries(headers, file_headers);
    hsp::str::assert_no_unused_entries(headers, file_headers);

    // assert that parameter is not already defined
    if (input.at("name").size() != hsp::str::unique(input.at("name")).size())
        throw std::invalid_argument(
            "Found multiple definitions for the same parameters.");

    clear_lists_();
    // load parameters from the row with options
    for (size_t k = 0; k < input.at("name").size(); ++k)
    {
        const std::string &paramstr = input.at("name")[k]; // parameter name

        if (!pRegistry->count(paramstr))
        {
            pRegistry->add(paramstr);
            param_vec_.emplace_back();
            param_vec_.back().name(paramstr, "parameter");
        }

        if (hsp::str::tobool(input.at("const")[k]))
            helper_build_constant_(input, k);
        else
            helper_build_to_sample_(input, k);
    }

    // initialize all parameters to standard
    list_sample_standard = list_sample;
    for (const size_t &p : list_sample_standard)
        if (at(p).is_adaptable())
            list_adaptable_sample_standard.push_back(p);
        else
            list_non_adapt_sample_standard.push_back(p);
}

/// Build parameter dictionaries for expected measure on posterior sample.
void Parameter::build_for_expec(const hsp::dbl_table &input)
{
    if (input.empty())
        return;

    const std::vector<std::string> headers = hsp::read::extract_keys(input);
    // assert that parameter is not already defined
    if (headers.size() != hsp::str::unique(headers).size())
        throw std::invalid_argument(
            "Found multiple definitions for the same parameters.");

    clear_lists_();

    // update list of constant parameters to all parameters
    list_const = pRegistry->from_string(pRegistry->get_all());

    // store as sample values
    smpl_values.clear();

    // reset definition dictionary to all constants NAs.
    for (const std::string &paramstr : pRegistry->get_all())
        if (std::find(headers.begin(), headers.end(), paramstr) == headers.end())
            smpl_values.emplace_back();
        else
            smpl_values.emplace_back(input.at(paramstr));
}

/// Drop current lists of parameters.
void Parameter::clear_lists_()
{
    size_t n = pRegistry->size();
    param_vec_.clear();
    param_vec_.reserve(n);
    param_vec_.assign(n, Variable("parameter"));
    for (size_t p = 0; p < n; ++p)
        param_vec_[p].name(pRegistry->as_string(p));

    list_const.clear();
    list_sample.clear();
    list_sample_standard.clear();
    list_adaptable_sample_standard.clear();
    list_non_adapt_sample_standard.clear();
    list_sample_cached.clear();
    list_adaptable_sample_cached.clear();
    list_non_adapt_sample_cached.clear();

    smpl_values.clear();
    smpl_priors.clear();
}

/// Helper to build constant parameter.
void Parameter::helper_build_constant_(const hsp::str_table &input, size_t k)
{
    const std::string &paramstr = input.at("name")[k]; // parameter name

    if (hsp::na::is(input.at("value")[k]))
        throw std::invalid_argument(
            "Constant parameter '" + paramstr + "' had undefined value.");

    // assign constant parameter to the set of variables
    const size_t p = pRegistry->from_string(paramstr);
    param_vec_[p].transfer(read_numeric(input.at("value")[k]));
    // check if discrete as requested
    if (user_discrete_set.count(p) && !param_vec_[p].is_discrete())
        throw std::invalid_argument("Definition of " + param_vec_[p].err_name() +
                                    " did not return a discrete object as requested by model.");

    // update the lists of constant parameters
    list_const.push_back(p);
}

/// Helper to build parameter to sample.
void Parameter::helper_build_to_sample_(const hsp::str_table &input, size_t k)
{
    const std::string &paramstr = input.at("name")[k]; // parameter name

    // assert parameter to sample has a valid value initialization
    if (hsp::na::is(input.at("value")[k]) && hsp::na::is(input.at("prior")[k]))
        throw std::invalid_argument(
            "Non-constant parameter '" + paramstr + "' had undefined prior.");

    // create a new parameter Variable
    const size_t p = pRegistry->from_string(paramstr);
    param_vec_[p].transfer(read_numeric(input.at("value")[k]));
    Variable &new_param = at(p);

    // assign a prior distribution to the parameter
    hsp::parse::Definition parser(input.at("prior")[k]);
    new_param.init_distr(parser.root.func);

    // load prior arguments
    for (const auto &child : parser.root.children)
        new_param.init_distr_attr(child.key, read_numeric(child));

    // load parameter range
    new_param.init_distr_attr(hsp::params::MIN, read_numeric(input.at(hsp::params::MIN)[k]));
    new_param.init_distr_attr(hsp::params::MAX, read_numeric(input.at(hsp::params::MAX)[k]));
    new_param.init_distr_reshaping(hsp::str::tobool(input.at("discr")[k]));

    // load proposal name
    std::string proposal_name = "sample";
    if (!hsp::na::is(input.at("proposal")[k]))
        proposal_name = input.at("proposal")[k];

    // load proposal scale
    double proposal_scale = 1.;
    if (!hsp::na::is(input.at("proposal_scale")[k]))
        proposal_scale = hsp::str::tod(input.at("proposal_scale")[k]);

    // assign the proposer to the parameter
    new_param.init_proposer(proposal_name, proposal_scale);
    // check if discrete as requested
    if (user_discrete_set.count(p) && !new_param.is_discrete())
        throw std::invalid_argument("Definition of " + new_param.err_name() +
                                    " did not return a discrete object as requested by model.");

    // update the lists of parameters to sample
    list_sample.push_back(p);
}

/// Adjust lists of parameters to be inferred making parameter as dependent on cache.
void Parameter::make_cache_dependent(size_t p)
{
    // remove a parameter from standard if it was defined as such
    if (std::find(list_sample_standard.begin(), list_sample_standard.end(), p) != list_sample_standard.end())
    {
        const std::vector<size_t> vec({p});
        list_sample_standard = list_sample_standard - vec;
        list_adaptable_sample_standard = list_adaptable_sample_standard - vec;
        list_non_adapt_sample_standard = list_non_adapt_sample_standard - vec;
    }

    // add the parameter as cache if not already present
    if (std::find(list_sample_cached.begin(), list_sample_cached.end(), p) == list_sample_cached.end())
    {
        list_sample_cached.push_back(p);
        if (at(p).is_adaptable())
            list_adaptable_sample_cached.push_back(p);
        else
            list_non_adapt_sample_cached.push_back(p);
    }
}

/// Load parameters values of the n-th posterrior sample.
void Parameter::load_sample_values(size_t n)
{
    for (const std::string &paramstr : pRegistry->get_all())
    {
        const size_t p = pRegistry->from_string(paramstr);
        const std::vector<double> &vec = smpl_values[p];
        if (!vec.empty())
            at(p).set(vec[n]);
    }
}

void Parameter::make_cache_dependent(const hsp::parse::Entry &node)
{
    // cache dependecy if an augmented function depends on a parameter to be inferred
    if (hsp::parse::is_parameter(node.type))
    {
        const size_t p = pRegistry->from_string(node.value);
        if (!at(p).is_constant())
            make_cache_dependent(p);
    }
}

/// Randomly initialize values of the parameters to be inferred.
void Parameter::user_reset(hsp::rnd_gen &prng)
{
    for (const size_t &p : list_sample)
        at(p).user_reset(prng);
}

/// Propose a valid move (return true), otherwise update refusal (false).
bool Parameter::propose_valid_move(size_t p, hsp::rnd_gen &prng)
{
    Variable &this_par = at(p);

    // set new value of the parameter and store the old to buffer
    this_par.propose(prng);

    if (this_par.is_valid_move())
        return true;

    // recover buffer if not valid move
    this_par.recover();
    return false;
}

/// Adpat the proposal of model parameters.
void Parameter::adapt_proposal_of_sampled(unsigned long step, unsigned long tot)
{
    const double ratio = static_cast<double>(step) / static_cast<double>(tot);

    // adpat all standard model parameters
    for (const size_t p : list_adaptable_sample_standard)
        at(p).adapt_proposal(ratio);

    // adapt all parameters with cached dependecies
    for (const size_t p : list_adaptable_sample_cached)
        at(p).adapt_proposal(ratio);
}

/// Initialize posterior sample of parameters and respective priors.
void Parameter::init_draws(unsigned long sample_size)
{
    const size_t n = list_sample.size();
    smpl_values.clear();
    smpl_priors.clear();
    smpl_values.reserve(n);
    smpl_priors.reserve(n);
    for (size_t i = 0; i < n; ++i)
    {
        smpl_values.emplace_back(sample_size, hsp::na::DOUBLE);
        smpl_priors.emplace_back(sample_size, hsp::na::DOUBLE);
    }
}

/// Store current sample of parameters and respective priors.
void Parameter::store_draws(unsigned long smpl_count)
{
    size_t i = 0;
    for (const size_t p : list_sample)
    {
        const Variable &this_par = at(p);
        smpl_values[i][smpl_count] = this_par.get();
        smpl_priors[i][smpl_count] = this_par.get_distr_logpdf();
        ++i;
    }
}

/// Auxiliary method to read a definition as a single attribute.
const Attribute &Parameter::read_attribute(const hsp::parse::Entry &node, const Individual *member_ptr) const
{
    storage_ = Attribute(); // reset auxiliary storage to default empty

    if (!node.is_terminal()) // function node
    {
        Attribute temporary(node.func);
        for (const auto &child : node.children)
            temporary.get_function()->init_attr(child.key, read_attribute(child, member_ptr));
        storage_ = std::move(temporary);
    }
    else if (hsp::parse::is_numeric(node.type)) // terminal numeric node
    {
        storage_ = Attribute(hsp::str::tod(node.value, "argument with numeric type"));
    }
    else if (hsp::parse::is_parameter(node.type)) // terminal parameter node
    {
        return at(pRegistry->from_string(node.value));
    }
    else if (!node.is_empty())
    {
        if (member_ptr == nullptr)
            throw std::invalid_argument(
                "Not supported individidual specific definition found while parsing '" + node.str() + "'.");

        const Individual &member = *member_ptr;
        if (hsp::parse::is_data(node.type)) // terminal data node
        {
            storage_ = Attribute(member.get_data_wide_vec()[pDataWideRegistry_->from_string(node.value)]);
            storage_.name(node.value, "data");
        }
        else if (hsp::parse::is_augmented(node.type)) // augmented terminal node
        {
            return member.get_augm_variable(pAugmRegistry_->from_string(node.value));
        }
        else if (hsp::parse::is_group_rate(node.type)) // group rate terminal node
        {
            return member.get_group_rate_attribute(pGroupRateRegistry_->from_string(node.value));
        }
        else
            throw std::invalid_argument(
                "Unrecognized source type '" + node.type + "' found while parsing '" + node.str() + "'.");
    }

    return storage_;
}

/// Convert a typed string into a double.
const Attribute &Parameter::read_numeric(const hsp::parse::Entry &node, const Individual *member_ptr) const
{
    const Attribute &aux_attr_ = read_attribute(node, member_ptr);
    if (!aux_attr_.is_empty() && !aux_attr_.is_constant())
        throw std::invalid_argument(
            "Definition '" + node.str() + "' could not be interpreted as a numeric constant.");
    return aux_attr_;
}

/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////  *** GroupRate : struct ***  ////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

GroupRate::GroupRate(const std::vector<std::string> &keys,
                     std::shared_ptr<DataWide> dataptr,
                     std::shared_ptr<Parameter> paramptr)
    : pRegistry(std::make_shared<hsp::key::Registry>(keys)),
      pDataWide_(dataptr),
      pParameter_(paramptr) {}
GroupRate::~GroupRate() = default;

/// Update model report and set group rates associate to new individual.
/// NOTE: this has to run even if input is empty to initialize existing rates to reference.
void GroupRate::build(const hsp::str_table &input)
{
    pRegistry->clear();

    size_t n_rows = 0;
    std::vector<std::string> data_headers = {};
    if (!input.empty())
    {
        data_headers = helper_assert_definitions_(input);
        n_rows = input.begin()->second.size();
    }
    else if (pRegistry->size() == 0)
        return;

    // loop over all household members
    for (Household &house : *pDataWide_)
        for (Individual &member : house)
        {
            hsp::str_table these_rates = {};

            for (size_t r = 0; r < n_rows; ++r)
            {
                // check if all conditions are satisfied
                bool is_this_class = true;
                for (const std::string &d : data_headers)
                    if (!hsp::na::is(input.at(d)[r]))
                        is_this_class &= member.data_vec_[pDataWide_->pRegistry->from_string(d)] == hsp::str::tod(input.at(d)[r]);

                // if it's the case, add the parameter name
                if (is_this_class)
                    these_rates[input.at("group_rate")[r]].push_back(input.at("param_name")[r]);
            }

            helper_init_all_rates_(member, these_rates);
        }
}

/// Verify that group rates definitions are correct.
std::vector<std::string> GroupRate::helper_assert_definitions_(const hsp::str_table &input)
{
    const std::vector<std::string> file_headers = {"group_rate", "param_name"};

    // check that tempalte headers are present
    std::vector<std::string> data_headers = hsp::read::extract_keys(input);
    hsp::str::assert_no_missing_entries(data_headers, file_headers);

    // remove template headers to keep only data
    data_headers = data_headers - file_headers;
    // check that required data columns are present
    hsp::str::assert_no_unused_entries(data_headers, pDataWide_->headers);

    // assert the group rate names are allowed
    for (const std::string &groupstr : input.at("group_rate"))
    {
        if (groupstr == hsp::rates::CONTACT || groupstr == hsp::rates::MISSING)
            throw std::invalid_argument(
                "Tried to use the reserved group rate name '" + groupstr + "'.\nPlease, choose a different name.");
        if (hsp::na::is(groupstr))
            throw std::invalid_argument(
                "Found a missing group rate name in the definitions.");
        if (!pRegistry->count(groupstr))
            pRegistry->add(groupstr);
    }

    return data_headers;
}

/// Convert to a dictionary of parameter pointers and pass the pointers to the individual.
void GroupRate::helper_init_all_rates_(Individual &member, const hsp::str_table &these_rates)
{
    // initialize group rates to reference
    member.group_rate_vec_.clear();
    member.group_rate_vec_.reserve(pRegistry->size());
    for (size_t g = 0; g < pRegistry->size(); ++g)
    {
        member.group_rate_vec_.emplace_back(Attribute(hsp::rates::VALUE_REFERENCE));
        member.group_rate_vec_.back().name(pRegistry->as_string(g), "group rate");
    }

    // assign custom definitions
    for (const auto &[groupstr, vec] : these_rates)
    {
        // extract non-reference parameter pointers
        std::vector<const Attribute *> param_pointers = {};
        for (const std::string &paramstr : vec)
        {
            if (hsp::na::is(paramstr) || paramstr == hsp::rates::MISSING)
                throw std::invalid_argument("A group rate definition '" + groupstr +
                                            "' was missing for " + member.err_name() + ".");
            if (paramstr != hsp::rates::REFERENCE)
                param_pointers.emplace_back(
                    pParameter_->at(pParameter_->pRegistry->from_string(paramstr)).get_ptr());
        }

        if (!param_pointers.empty())
        {
            size_t g = pRegistry->from_string(groupstr);
            member.group_rate_vec_[g].transfer(Attribute(param_pointers));
        }
    }
}

/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////  *** InfecProfile : struct ***  //////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

/// Infectivity profile constructor.
InfecProfile::InfecProfile(std::shared_ptr<Parameter> paramptr)
    : pParameter_(paramptr) {}
/// Infectivity profile deconstructor.
InfecProfile::~InfecProfile() = default;

/// Set model template while constructiing model.
void InfecProfile::init_template(hsp::str_table input)
{
    // if not previous build, proceed with template
    // if no template, assume no need for an infectivity profile
    if (generic_profile != nullptr ||
        !augm_profile_definitions_.empty() ||
        input.empty())
        return;

    try
    {
        if (!input.count("func"))
            throw std::invalid_argument(
                "A function name was missing for infectivity profile template.");

        const size_t n_rows = hsp::read::csv_n_rows(input);
        const auto file_headers = {"func", "status", "discr", hsp::params::MAX};

        for (const auto &key : file_headers)
            if (!input.count(key))
                input.emplace(key, std::vector<std::string>(n_rows, ""));

        // discretization
        input["discr"] = std::vector<std::string>(n_rows, "false");
        if (is_discrete_requested_)
            input["discr"] = std::vector<std::string>(n_rows, "true");
        else if (input.count("discr"))
        {
            for (size_t r = 0; r < n_rows; ++r)
                if (!hsp::na::is(input.at("discr")[r]))
                    input["discr"][r] = input.at("discr")[r];
        }

        build(input); // run build as if external input

        // forbid augmented infectivity profile for template
        if (is_augmented())
            throw std::invalid_argument(
                "Impossible to define infectivity profile template with augmented dependencies.");

        // UPDATE: maybe also not depending on parameters?
    }
    catch (const std::exception &e)
    {
        throw std::runtime_error(
            "While building infectivity profile template:\n" + std::string(e.what()));
    }
}

/// Store values of the infectivity profile.
void InfecProfile::init_lookup_()
{
    inst_ = &InfecProfile::do_inst_lookup_;
    cumul_ = &InfecProfile::do_cumul_lookup_;

    if (generic_profile->get_upper().is_empty())
        size_lut_ = hsp::infprof::LOOKUP_SIZE;
    else
        size_lut_ = static_cast<size_t>(generic_profile->upper()) + 1;

    lookup_inst_.resize(size_lut_);
    lookup_cumul_.resize(size_lut_);
    lookup_buffer_inst_.resize(size_lut_);
    lookup_buffer_cumul_.resize(size_lut_);

    update();
}

/// Clear infectivity profile lookup table.
void InfecProfile::clear_lookup_()
{
    size_lut_ = 0;

    lookup_inst_.clear();
    lookup_cumul_.clear();
    lookup_buffer_inst_.clear();
    lookup_buffer_cumul_.clear();
}

/// If infectivity profile is a lookup table, store its values.
/// UPDATE: this should not be recoverd if constant.
void InfecProfile::update()
{
    if (!is_lookup())
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
void InfecProfile::recover()
{
    if (!is_lookup())
        return;

    lookup_inst_.swap(lookup_buffer_inst_);
    lookup_cumul_.swap(lookup_buffer_cumul_);
}

/// Build infectivity profile distribution.
/// NOTE: at the moment I am ignoring minimum
void InfecProfile::build(const hsp::str_table &input)
{
    augm_profile_definitions_.clear();

    if (input.empty())
        return;

    // check all entries are legitimiate
    hsp::str::assert_no_unused_entries(
        hsp::read::extract_keys(input),
        {"func", "status", "discr", hsp::params::MIN, hsp::params::MAX});

    // check definition list has a function
    if (!input.count("func"))
        throw std::invalid_argument(
            "A function name choice was missing in infectivity profile definition.");
    for (const std::string &in_name : input.at("func"))
        if (hsp::parse::Definition(in_name).is_na())
            throw std::invalid_argument(
                "Provided NA as a function name in infectivity profile definition.");

    // check if definition is generic or augmented
    if (hsp::read::csv_n_rows(input) > 1 ||
        hsp::parse::Definition(input.at("func")[0]).is_augmented_depen() ||
        hsp::parse::Definition(input.at(hsp::params::MAX)[0]).is_augmented_depen() ||
        !hsp::na::is(input.at("status")[0]))
        build_augmented_profile_(input);
    else
        build_generic_profile_(input);
}

/// Build for generic infectivity profile.
/// NOTE: no status dependence allowed here.
void InfecProfile::build_generic_profile_(const hsp::str_table &input)
{
    inst_ = &InfecProfile::do_inst_default_;
    cumul_ = &InfecProfile::do_cumul_default_;

    generic_profile = helper_build_profile_({{"func", {input.at("func")[0]}},
                                             {"discr", {input.count("discr") ? input.at("discr")[0] : "false"}}});

    // range minimum
    if (input.count(hsp::params::MIN))
        generic_profile->init_attr(hsp::params::MIN, pParameter_->read_attribute(input.at(hsp::params::MIN)[0]));
    else
    {
        const std::string default_min_expr = generic_profile->is_discrete() ? "1|N" : "0|N";
        generic_profile->init_attr(hsp::params::MIN, pParameter_->read_attribute(default_min_expr));
    }

    // range maximum
    if (input.count(hsp::params::MAX))
        generic_profile->init_attr(hsp::params::MAX, pParameter_->read_attribute(input.at(hsp::params::MAX)[0]));
    // range truncation
    generic_profile->init_truncation();

    // if discrete and generic, initializes lookup table
    clear_lookup_();
    if (generic_profile->is_discrete())
        init_lookup_();
}

/// Build for individual augmented function infectivity profile.
void InfecProfile::build_augmented_profile_(const hsp::str_table &input)
{
    // use helper to check profile name and discretization
    for (size_t r = 0; r < hsp::read::csv_n_rows(input); ++r)
        helper_build_profile_({{"func", {hsp::parse::Definition(input.at("func")[r]).root.func}},
                               {"discr", {input.at("discr")[r]}}});
    // drop generic definition
    generic_profile = nullptr;
    clear_lookup_();

    inst_ = &InfecProfile::do_inst_augm_;
    cumul_ = &InfecProfile::do_cumul_augm_;
    augm_profile_definitions_.clear();

    for (size_t r = 0; r < hsp::read::csv_n_rows(input); ++r)
    {
        // initialize for instantaneous and cumulative profiles
        bool is_discrete = hsp::str::tobool(input.at("discr")[r]);

        // match infectivity profile and augmented variable definition
        augm_profile_definitions_["func"].push_back(hsp::infprof::INSTANT);
        augm_profile_definitions_["status"].push_back(input.at("status")[r]);
        augm_profile_definitions_["func"].push_back(input.at("func")[r]);
        augm_profile_definitions_["fcumul"].push_back("false");
        augm_profile_definitions_["value"].push_back("");
        augm_profile_definitions_["const"].push_back("false");
        augm_profile_definitions_["discr"].push_back(hsp::str::format(is_discrete));
        augm_profile_definitions_["model"].push_back("");
        augm_profile_definitions_[hsp::params::MIN].push_back(is_discrete ? "1|N" : "0|N");
        augm_profile_definitions_[hsp::params::MAX].push_back(input.at(hsp::params::MAX)[r]);
        augm_profile_definitions_["proposal"].push_back("");
        augm_profile_definitions_["proposal_scale"].push_back("");
    }
}

std::unique_ptr<Distribution> InfecProfile::helper_build_profile_(const hsp::str_table &input)
{
    hsp::parse::Definition parser(hsp::parse::Definition(input.at("func")[0]));
    auto placeholder = compose_distr_registries(parser.root.func,
                                                "infectivity profile",
                                                {"infecprofile", "continuos", "discrete"});
    for (const auto &node : parser.root.children)
    {
        placeholder->init_attr(node.key, pParameter_->read_attribute(node));
        pParameter_->make_cache_dependent(node);
    }

    // discretization
    if (input.count("discr") && hsp::str::tobool(input.at("discr")[0]))
        placeholder->init_discretization();
    else if (is_discrete_requested_)
        hsp::msg::error(
            "Definition of infectivity profile did not return a discrete object as requested by the model.");
    return placeholder;
}

/////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////  MCMC : Impl : InfecProfile : getter methods  ///////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

/// Default method for instantaneous infectivity profile.
double InfecProfile::do_inst_default_(const Individual *, double t)
{
    return generic_profile->pdf(t);
}
/// Lookup method for discrete instantaneous infectivity profile.
double InfecProfile::do_inst_lookup_(const Individual *, double t)
{
    if (generic_profile->is_valid(t) && t < static_cast<double>(size_lut_))
        return lookup_inst_.at(static_cast<size_t>(t));
    return 0;
}
/// Augmented method for instantaneous infectivity profile.
double InfecProfile::do_inst_augm_(const Individual *member, double t)
{
    return member->get_augm_variable(augm_infecprof_key).get_function()->pdf(t);
}

/// Default method for cumulative infectivity profile.
double InfecProfile::do_cumul_default_(const Individual *, double t)
{
    return generic_profile->cdf(t);
}
/// Lookup method for discrete cumulative infectivity profile.
double InfecProfile::do_cumul_lookup_(const Individual *, double t)
{
    if (generic_profile->is_below(t))
        return 0.;
    if (t >= static_cast<double>(size_lut_) || generic_profile->is_above(t))
        return lookup_cumul_.at(size_lut_ - 1);
    return lookup_cumul_.at(static_cast<size_t>(t));
}
/// Augmented method for cumulative infectivity profile.
double InfecProfile::do_cumul_augm_(const Individual *member, double t)
{
    return member->get_augm_variable(augm_infecprof_key).get_function()->cdf(t);
}
