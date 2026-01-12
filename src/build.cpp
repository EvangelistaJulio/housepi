/* housepi - Epidemiological Analysis of Household Outbreaks.
 *
 * Copyright © 2026 Institut Pasteur (Paris), CNRS and INSERM.
 * Distributed under the terms of the GNU General Public License (GPLv3).
 * See the LICENSE file for details.
 */

#include "tools/global.h"
#include "tools/messages.h"
#include "tools/operators.h"
#include "tools/random.h"
#include "tools/read.h"
#include "tools/string.h"
#include "tools/types.h"

#include "attribute.h"
#include "function.h"
#include "function_registry.h"
#include "distribution.h"
#include "distribution_registry.h"
#include "variable.h"
#include "augmented_map.h"
#include "data.h"
#include "parameter.h"
#include "group_rate.h"
#include "infec_profile.h"

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////  *** AugmentedMap : class ***  ////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/// A string containing a summary of the augmented map.
std::string AugmentedMap::str_augm_summary() const
{
    std::string output;
    for (const auto &pair : augm_map_)
    {
        output += "- " + pair.first + " : ";
        if (pair.second.is_function() && !pair.second.get_function()->is_callable_implicitly())
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
    // clean pre-existing
    for (auto &pair : model_logratios_)
        pair.second.clear();
    model_logratios_.clear();

    // associate models if applicable
    for (const auto &pair : models_ordered)
    {
        const std::string &a = pair.first;
        auto &vec = model_logratios_[a];
        vec.clear();
        vec.reserve(pair.second.size());

        // repeat for all dependent models
        for (const std::string &b : pair.second)
            vec.emplace_back(&augm_map_[b]);
    }
}

/// Value of log-likelihood associated to augmented variable models.
/// Store the log-likelihood of individual augmented variables.
double AugmentedMap::materialize_augm_model_log_lik()
{
    buffer_ = augm_log_lik_; // store old value in buffer
    augm_log_lik_ = 0.0;
    for (const auto &pair : model_logratios_)
        augm_log_lik_ += augm_map_[pair.first].get_distr_logpdf();
    assert_log_lik_not_NaN();
    return augm_log_lik_;
}

/// Change in log-likelihood associated to given augmented.
double AugmentedMap::materialize_augm_model_logratio(const std::string &a)
{
    buffer_ = augm_log_lik_; // store old value in buffer
    double delta_ll = 0.;
    for (const Variable *xaugm : model_logratios_[a])
        delta_ll += xaugm->get_distr_logratio();
    augm_log_lik_ += delta_ll;
    assert_log_lik_not_NaN(a);
    return delta_ll;
}

/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////  *** Parameter : struct ***  ////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

/// Assign parameter definitions to dictionary.
void Parameter::build(const u_dict<str_vec> &input)
{
    if (input.empty())
        return;

    // check entries
    const str_vec file_headers = {"name", "value", "const", "discr",
                                  "prior", hsp::params::MIN, hsp::params::MAX,
                                  "proposal", "proposal_scale"};
    const str_vec headers = hsp::read::extract_keys(input);
    hsp::str::assert_no_missing_entries(headers, file_headers);
    hsp::str::assert_no_unused_entries(headers, file_headers);
    // assert that parameter is not already defined
    if (input.at("name").size() != hsp::str::unique(input.at("name")).size())
        throw std::invalid_argument(
            "Found multiple definitions for the same parameters.");

    dictionary_.clear();
    clear_lists_();

    // load parameters from the row with options
    for (size_t k = 0; k < input.at("name").size(); ++k)
        if (hsp::str::tobool(input.at("const")[k]))
            build_constant_parameter_(input, k);
        else
            build_parameter_to_sample_(input, k);

    // initialize all parameters to standard
    list_sample_standard = list_sample;
    for (const std::string &p : list_sample_standard)
        if (at(p).is_adaptable())
            list_adaptable_sample_standard.push_back(p);
        else
            list_non_adapt_sample_standard.push_back(p);
}

/// Build parameter dictionaries for expected measure on posterior sample.
void Parameter::build_for_expec(const u_dict<dbl_vec> &input)
{
    if (input.empty())
        return;

    const str_vec headers = hsp::read::extract_keys(input);
    // assert that parameter is not already defined
    if (headers.size() != hsp::str::unique(headers).size())
        throw std::invalid_argument(
            "Found multiple definitions for the same parameters.");

    clear_lists_();

    // store values as a posterior sample
    smpl_values = input;

    // reset definition dictionary to all constants NAs.
    for (const std::string &p : headers)
        if (dictionary_.count(p))
            at(p).transfer(Attribute(hsp::na::DOUBLE));
        else
            dictionary_.emplace(p, Variable(p, "parameter", hsp::na::DOUBLE));
    // update list of constant parameters
    list_const = hsp::read::extract_keys(dictionary_);
}

/// Drop current lists of parameters.
void Parameter::clear_lists_()
{
    list_const.clear();
    list_sample.clear();
    list_sample_standard.clear();
    list_adaptable_sample_standard.clear();
    list_non_adapt_sample_standard.clear();
    list_sample_cached.clear();
    list_adaptable_sample_cached.clear();
    list_non_adapt_sample_cached.clear();
}

/// Helper to build constant parameter.
void Parameter::build_constant_parameter_(const u_dict<str_vec> &input, size_t k)
{
    const std::string &p = input.at("name")[k]; // parameter name

    if (hsp::na::is(input.at("value")[k]))
        throw std::invalid_argument(
            "Constant parameter '" + p + "' had undefined value.");

    // assign constant parameter to the set of variables
    double value = read_numeric(input.at("value")[k]).get();
    dictionary_.emplace(p, Variable(p, "parameter", value));
    // update the lists of constant parameters
    list_const.push_back(p);
}

/// Helper to build parameter to sample.
void Parameter::build_parameter_to_sample_(const u_dict<str_vec> &input, size_t k)
{
    const std::string &p = input.at("name")[k]; // parameter name

    // assert parameter to sample has a valid value initialization
    if (hsp::na::is(input.at("value")[k]) && hsp::na::is(input.at("prior")[k]))
        throw std::invalid_argument(
            "Non-constant parameter '" + p + "' had undefined prior.");

    // load proposal name
    std::string proposal_name = "sample";
    if (!hsp::na::is(input.at("proposal")[k]))
        proposal_name = input.at("proposal")[k];

    // load proposal scale
    double proposal_scale = 1.;
    if (!hsp::na::is(input.at("proposal_scale")[k]))
        proposal_scale = hsp::str::tod(input.at("proposal_scale")[k]);

    // create a new parameter Variable
    dictionary_.emplace(
        p, Variable(p, "parameter", read_numeric(input.at("value")[k]).get(), false));
    Variable &new_param = at(p);

    // assign a prior distribution to the parameter
    hsp::parse::Definition parser(input.at("prior")[k]);
    new_param.init_distr(parser);
    // load prior arguments
    for (const auto &child : parser.root.children)
        new_param.init_distr_attr(child.key, read_numeric(child));
    // load parameter range
    new_param.init_distr_attr(hsp::params::MIN, read_numeric(input.at(hsp::params::MIN)[k]));
    new_param.init_distr_attr(hsp::params::MAX, read_numeric(input.at(hsp::params::MAX)[k]));
    new_param.init_distr_reshaping(hsp::str::tobool(input.at("discr")[k]));

    // assign the proposer to the parameter
    new_param.init_proposer(proposal_name, proposal_scale);
    // update the lists of parameters to sample
    list_sample.push_back(p);
}

/// Adjust lists of parameters to be inferred making parameter as dependent on cache.
void Parameter::make_cache_dependent(const std::string &p)
{
    // remove a parameter from standard if it was defined as such
    if (list_sample_standard % p)
    {
        list_sample_standard = list_sample_standard - str_vec({p});
        list_adaptable_sample_standard = list_adaptable_sample_standard - str_vec({p});
        list_non_adapt_sample_standard = list_non_adapt_sample_standard - str_vec({p});
    }

    // add the parameter as cache if not alread present
    if (!(list_sample_cached % p))
    {
        list_sample_cached.push_back(p);
        if (at(p).is_adaptable())
            list_adaptable_sample_cached.push_back(p);
        else
            list_non_adapt_sample_cached.push_back(p);
    }
}

void Parameter::make_cache_dependent(const hsp::parse::Entry &node)
{
    // cache dependecy if an augmented function depends on a parameter to be inferred
    if (hsp::sources::is_parameter(node.type) && !at(node.value).is_constant())
        make_cache_dependent(node.value);
}

/// Randomly initialize values of the parameters to be inferred.
void Parameter::user_reset(rnd_gen &prng)
{
    for (const std::string &p : list_sample)
        at(p).user_reset(prng);
}

/// Load parameters values of the n-th posterrior sample.
void Parameter::load_sample_values(size_t n)
{
    for (const auto &pair : smpl_values)
        at(pair.first).set(pair.second[n]);
}

/// Initialize posterior sample of parameters and respective priors.
void Parameter::init_draws(unsigned long sample_size)
{
    smpl_values.clear();
    smpl_values.reserve(list_sample.size());
    smpl_priors.clear();
    smpl_priors.reserve(list_sample.size());
    for (const std::string &p : list_sample)
    {
        smpl_values[p].clear();
        smpl_values[p].reserve(sample_size);
        smpl_values[p].assign(sample_size, hsp::na::DOUBLE);
        smpl_priors[p].clear();
        smpl_priors[p].reserve(sample_size);
        smpl_priors[p].assign(sample_size, hsp::na::DOUBLE);
    }
}

/// Store current sample of parameters and respective priors.
void Parameter::store_draws(unsigned long smpl_count)
{
    for (const std::string &p : list_sample)
    {
        const Variable &this_par = this->at(p);
        smpl_values[p][smpl_count] = this_par.get();
        smpl_priors[p][smpl_count] = this_par.get_distr_logpdf();
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
    else if (hsp::sources::is_numeric(node.type)) // terminal numeric node
    {
        storage_ = Attribute(hsp::str::tod(node.value, "argument with numeric type"));
    }
    else if (hsp::sources::is_parameter(node.type)) // terminal parameter node
    {
        return at(node.value);
    }
    else if (!node.is_empty())
    {
        if (member_ptr == nullptr)
            throw std::invalid_argument(
                "Not supported individidual specific definition found while parsing '" + node.str() + "'.");

        const Individual &member = *member_ptr;
        if (hsp::sources::is_data(node.type)) // terminal data node
        {
            storage_ = Attribute(member.data(node.value));
            storage_.set_name(node.value, "data");
        }
        else if (hsp::sources::is_augmented(node.type)) // augmented terminal node
        {
            return member.get_augm_variable(node.value);
        }
        else if (hsp::sources::is_group_rate(node.type)) // group rate terminal node
        {
            return member.get_group_rate_attribute(node.value);
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

GroupRate::GroupRate(std::shared_ptr<DataWide> dataptr,
                     std::shared_ptr<Parameter> paramptr)
    : pDataWide_(dataptr),
      pParameter_(paramptr) {}
GroupRate::~GroupRate() = default;

/// Update model report and set group rates associate to new individual.
void GroupRate::build(const u_dict<str_vec> &input)
{
    if (input.empty())
        return;

    const str_vec data_headers(helper_assert_definitions_(input));

    // loop over all household members
    const size_t n_rows = input.at("group_rate").size();
    for (Household &house : *pDataWide_)
        for (Individual &member : house)
        {
            u_dict<str_vec> these_rates = {};

            // obtain group rates for given member
            for (size_t r = 0; r < n_rows; ++r)
            {
                // check if all conditions are satisfied
                bool is_this_class = true;
                for (const std::string &d : data_headers)
                    if (!hsp::na::is(input.at(d)[r]))
                        is_this_class &= (member.data(d) == hsp::str::tod(input.at(d)[r]));

                // if it's the case, add the parameter name
                if (is_this_class)
                    these_rates[input.at("group_rate")[r]].push_back(input.at("param_name")[r]);
            }

            // assign rates to the individual
            helper_init_all_rates_(member, these_rates);
        }
}

/// Verify that group rates definitions are correct.
str_vec GroupRate::helper_assert_definitions_(const u_dict<str_vec> &input)
{
    // check that tempalte headers are present
    str_vec data_headers = hsp::read::extract_keys(input);
    hsp::str::assert_no_missing_entries(data_headers, file_headers_);

    // remove template headers to keep only data
    data_headers = data_headers - file_headers_;
    // check that required data columns are present
    hsp::str::assert_no_unused_entries(data_headers, pDataWide_->headers);

    // assert the group rate names are allowed
    headers.clear();
    for (const std::string &g : input.at("group_rate"))
    {
        if (g == hsp::rates::CONTACT || g == hsp::rates::MISSING)
            throw std::invalid_argument(
                "Tried to use the reserved group rate name '" + g + "'.\nPlease, choose a different name.");
        if (hsp::na::is(g))
            throw std::invalid_argument(
                "Found a missing group rate name in the definitions.");
        headers.push_back(g);
    }
    headers = hsp::str::unique(headers);

    return data_headers;
}

/// Convert to a dictionary of parameter pointers and pass the pointers to the individual.
void GroupRate::helper_init_all_rates_(Individual &member, const u_dict<str_vec> &these_rates)
{
    for (const auto &pair : these_rates)
    {
        const std::string &g = pair.first;

        // extract non-reference parameter pointers
        std::vector<const Attribute *> vec_par_ptr = {};
        for (const std::string &p : pair.second)
        {
            if (hsp::na::is(p) || p == hsp::rates::MISSING)
                throw std::invalid_argument(
                    "A group rate defintion '" + g + "' was missing for " + member.err_name() + ".");
            if (p != hsp::rates::REFERENCE)
                vec_par_ptr.emplace_back(pParameter_->at(p).get_ptr());
        }

        if (vec_par_ptr.empty()) // stick to reference
            member.init_group_rate_to_reference(g);

        else if (vec_par_ptr.size() == 1) // use simple pointer
            member.init_group_rate(g, Attribute(vec_par_ptr[0]));

        else // use variadic product
        {
            /// UPDATE: I could define a call to deal with this cases
            /// to avoid including headers of function and registries...
            std::shared_ptr<Function> fun = compose_func_registries(
                "prod", "function", {&ArithmeticRegistry::inst()}, false);

            // initialize attributes to pointers
            fun->init_attr(hsp::params::X, Attribute(vec_par_ptr[0]));
            for (size_t i = 1; i < vec_par_ptr.size(); ++i)
                fun->init_attr(hsp::params::X + std::to_string(i + 1), Attribute(vec_par_ptr[i]));

            member.init_group_rate(g, Attribute(fun));
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
    : generic_profile(std::make_unique<Distribution>()),
      pParameter_(paramptr) {}
/// Infectivity profile deconstructor.
InfecProfile::~InfecProfile() = default;

/// Set model template while constructiing model.
void InfecProfile::init_template(u_dict<str_vec> input)
{
    // if not previous build, proceed with template
    // if no template, assume no need for an infectivity profile
    if (!hsp::na::is(generic_profile->name()) ||
        !augm_profile_definitions_.empty() ||
        input.empty())
        return;

    try
    {
        if (!input.count("func"))
            throw std::invalid_argument(
                "A function name was missing for infectivity profile template.");

        const size_t n_rows = hsp::read::csv_n_rows(input);
        const str_vec file_headers = {"func", "status", "discr", hsp::params::MAX};

        for (const std::string &key : file_headers)
            if (!input.count(key))
                input.emplace(key, str_vec(n_rows, hsp::na::STRING));

        // discretization
        input["discr"] = str_vec(n_rows, "false");
        if (is_discrete_requested_)
            input["discr"] = str_vec(n_rows, "true");
        else if (input.count("discr"))
        {
            for (size_t r = 0; r < n_rows; ++r)
                if (!hsp::na::is(input.at("discr")[r]))
                    input["discr"][r] = input.at("discr")[r];
        }

        build(input);

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
        size_lut_ = generic_profile->upper() + 1;

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
void InfecProfile::build(const u_dict<str_vec> &input)
{
    augm_profile_definitions_.clear();

    if (input.empty())
        return;

    // check entries
    const str_vec file_headers = {"func", "status", "discr",
                                  hsp::params::MIN, hsp::params::MAX};
    hsp::str::assert_no_unused_entries(hsp::read::extract_keys(input), file_headers);
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
void InfecProfile::build_generic_profile_(const u_dict<str_vec> &input)
{
    inst_ = &InfecProfile::do_inst_default_;
    cumul_ = &InfecProfile::do_cumul_default_;

    generic_profile = helper_build_profile_({{"func", input.at("func")[0]},
                                             {"discr", input.count("discr") ? input.at("discr")[0] : "false"}});

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
    if (generic_profile->is_discrete())
        init_lookup_();
    else
        clear_lookup_();
}

/// Build for individual augmented function infectivity profile.
void InfecProfile::build_augmented_profile_(const u_dict<str_vec> &input)
{
    // use helper to check profile name and discretization
    for (size_t r = 0; r < hsp::read::csv_n_rows(input); ++r)
        helper_build_profile_({{"func", hsp::parse::Definition(input.at("func")[r]).root.func},
                               {"discr", input.at("discr")[r]}});
    // drop generic definition
    generic_profile = std::make_unique<Distribution>();
    clear_lookup_();

    inst_ = &InfecProfile::do_inst_augm_;
    cumul_ = &InfecProfile::do_cumul_augm_;
    augm_profile_definitions_.clear();

    for (size_t r = 0; r < hsp::read::csv_n_rows(input); ++r)
    {
        // initialize for instantaneous and cumulative profiles
        bool is_discrete = hsp::str::tobool(input.at("discr")[r]);

        // match infectivity profile and augmented variable definition
        augm_profile_definitions_["func"] += {hsp::infprof::INSTANT, hsp::infprof::CUMUL};
        augm_profile_definitions_["status"] += str_vec(2, input.at("status")[r]);
        augm_profile_definitions_["func"] += str_vec(2, input.at("func")[r]);
        augm_profile_definitions_["fcumul"] += {"false", "true"};
        augm_profile_definitions_["value"] += str_vec(2, hsp::na::STRING);
        augm_profile_definitions_["const"] += str_vec(2, "false");
        augm_profile_definitions_["discr"] += str_vec(2, hsp::str::format(is_discrete));
        augm_profile_definitions_["model"] += str_vec(2, hsp::na::STRING);
        augm_profile_definitions_[hsp::params::MIN] += str_vec(2, (is_discrete) ? "1|N" : "0|N");
        augm_profile_definitions_[hsp::params::MAX] += str_vec(2, input.at(hsp::params::MAX)[r]);
        augm_profile_definitions_["proposal"] += str_vec(2, hsp::na::STRING);
        augm_profile_definitions_["proposal_scale"] += str_vec(2, hsp::na::STRING);
    }
}

std::unique_ptr<Distribution> InfecProfile::helper_build_profile_(const u_dict<std::string> &input)
{
    hsp::parse::Definition parser(hsp::parse::Definition(input.at("func")));
    auto placeholder = compose_distr_registries(parser.root.func,
                                                "infectivity profile",
                                                {&InfecProfileRegistry::inst(),
                                                 &ContinuosDistrRegistry::inst(),
                                                 &DiscreteDistrRegistry::inst()});
    for (const auto &node : parser.root.children)
    {
        placeholder->init_attr(node.key, pParameter_->read_attribute(node));
        pParameter_->make_cache_dependent(node);
    }

    // discretization
    if (input.count("discr") && hsp::str::tobool(input.at("discr")))
        placeholder->init_discretization();
    else if (is_discrete_requested_)
        hsp::msg::error(
            "Definition of infectivity profile did not return a discrete object as requested by the model.");
    return placeholder;
}

// Eventually append infectivity profile definitions to augmented.
void InfecProfile::transfer_augm_definitions(u_dict<str_vec> &augm_deftable) const
{
    if (!is_augmented())
        return;

    for (const auto &pair : augm_profile_definitions_)
        augm_deftable.at(pair.first) += pair.second;
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
    if (generic_profile->is_valid(t) && t < size_lut_)
        return lookup_inst_.at(static_cast<size_t>(t));
    return 0;
}
/// Augmented method for instantaneous infectivity profile.
double InfecProfile::do_inst_augm_(const Individual *i, double t)
{
    return i->get_augm_variable(hsp::infprof::INSTANT)(t);
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
    if (t >= size_lut_ || generic_profile->is_above(t))
        return lookup_cumul_.at(size_lut_ - 1);
    return lookup_cumul_.at(static_cast<size_t>(t));
}
/// Augmented method for cumulative infectivity profile.
double InfecProfile::do_cumul_augm_(const Individual *i, double t)
{
    return i->get_augm_variable(hsp::infprof::CUMUL)(t);
}
