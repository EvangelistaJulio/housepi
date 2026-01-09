/* housepi - Epidemiological Analysis of Household Outbreaks.
 *
 * Copyright © 2026 Institut Pasteur (Paris), CNRS and INSERM.
 * Distributed under the terms of the GNU General Public License (GPLv3).
 * See the LICENSE file for details.
 */

#include <chrono>

#include "tools/types.h"
#include "tools/na_utils.h"
#include "tools/string.h"
#include "tools/read.h"
#include "tools/messages.h"
#include "tools/global.h"

#include "function.h"
#include "function_registry.h"
#include "variable.h"
#include "data.h"
#include "data_reports.h"
#include "mcmc_impl.h"

MCMC::Impl::Impl(std::function<void()> func_init_inference,
                 std::function<void()> func_init_simulation,
                 std::function<double(Household &)> func_household_log_lik,
                 std::function<void(Household &)> func_simulate_house,
                 std::function<void(Household &, Individual &)> func_update_cache)
    : atomic_error(std::make_unique<hsp::atomic::Error>()),
      setup(std::make_shared<Setup>()),
      data_wide(std::make_shared<DataWide>()),
      data_long(std::make_shared<DataLong>(data_wide)),
      param(*this), 
      group_rates(*this), 
      infec_profile(*this), 
      augm(*this),
      callback_init_inference_(std::move(func_init_inference)),
      callback_init_simulation_(std::move(func_init_simulation)),
      callback_house_log_lik_(std::move(func_household_log_lik)),
      callback_simulate_house_(std::move(func_simulate_house)),
      callback_update_cache_(std::move(func_update_cache)) {}
MCMC::Impl::~Impl() = default;

/////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////  MCMC : cache attributes  ///////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

/// Extra integer vectors for all households.
void MCMC::cache_house_integer_vector(const std::string &e)
{
    for (Household &house : pImpl_->data_wide->restricted())
        house.init_cache_integer_vector(e);
}
/// Extra size_t vectors for all households.
void MCMC::cache_house_size_t_vector(const std::string &e)
{
    for (Household &house : pImpl_->data_wide->restricted())
        house.init_cache_size_t_vector(e);
}
/// Extra double vectors for all households.
void MCMC::cache_house_double_vector(const std::string &e)
{
    for (Household &house : pImpl_->data_wide->restricted())
        house.init_cache_double_vector(e);
}
/// Extra double square matrix for all households.
void MCMC::cache_house_double_matrix(const std::string &e)
{
    for (Household &house : pImpl_->data_wide->restricted())
        house.init_cache_double_matrix(e);
}

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////  MCMC : Impl : GroupRates : build  ///////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/// Update model report and set group rates associate to new individual.
void MCMC::Impl::GroupRates::build(const u_dict<str_vec> &input)
{
    if (input.empty())
        return;

    const str_vec data_headers(helper_assert_definitions_(input));

    // loop over all household members
    const size_t n_rows = input.at("group_rate").size();
    for (Household &house : *pImpl_.data_wide)
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
str_vec MCMC::Impl::GroupRates::helper_assert_definitions_(const u_dict<str_vec> &input)
{
    // check that tempalte headers are present
    str_vec data_headers = hsp::read::extract_keys(input);
    hsp::str::assert_no_missing_entries(data_headers, file_headers_);

    // remove template headers to keep only data
    data_headers = data_headers - file_headers_;
    // check that required data columns are present
    hsp::str::assert_no_unused_entries(data_headers, pImpl_.data_wide->headers);

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
void MCMC::Impl::GroupRates::helper_init_all_rates_(Individual &member,
                                                    const u_dict<str_vec> &these_rates)
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
                vec_par_ptr.emplace_back(pImpl_.param.at(p).get_ptr());
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

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////  MCMC : Impl : Parameter : build  ////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/// Assign parameter definitions to dictionary.
void MCMC::Impl::Parameter::build(const u_dict<str_vec> &input)
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
void MCMC::Impl::Parameter::build_for_expec(const u_dict<dbl_vec> &input)
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
void MCMC::Impl::Parameter::clear_lists_()
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
void MCMC::Impl::Parameter::build_constant_parameter_(const u_dict<str_vec> &input, size_t k)
{
    const std::string &p = input.at("name")[k];

    if (hsp::na::is(input.at("value")[k]))
        throw std::invalid_argument(
            "Constant parameter '" + p + "' had undefined value.");

    // assign constant parameter to the set of variables
    double value = pImpl_.read_numeric_(input.at("value")[k]).get();
    dictionary_.emplace(p, Variable(p, "parameter", value));
    // update the lists of constant parameters
    list_const.push_back(p);
}

/// Helper to build parameter to sample.
void MCMC::Impl::Parameter::build_parameter_to_sample_(const u_dict<str_vec> &input, size_t k)
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
        p, Variable(p, "parameter", pImpl_.read_numeric_(input.at("value")[k]).get(), false));
    Variable &new_param = at(p);

    // assign a prior distribution to the parameter
    hsp::parse::Definition parser(input.at("prior")[k]);
    new_param.init_distr(parser);
    // load prior arguments
    for (const auto &child : parser.root.children)
        new_param.init_distr_attr(child.key, pImpl_.read_numeric_(child));
    // load parameter range
    new_param.init_distr_attr(hsp::params::MIN, pImpl_.read_numeric_(input.at(hsp::params::MIN)[k]));
    new_param.init_distr_attr(hsp::params::MAX, pImpl_.read_numeric_(input.at(hsp::params::MAX)[k]));
    new_param.init_distr_reshaping(hsp::str::tobool(input.at("discr")[k]));

    // assign the proposer to the parameter
    new_param.init_proposer(proposal_name, proposal_scale);
    // update the lists of parameters to sample
    list_sample.push_back(p);
}

/// Adjust lists of parameters to be inferred making parameter as dependent on cache.
void MCMC::Impl::Parameter::make_cache_dependent(const std::string &p)
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

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////  MCMC : constructor  /////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

MCMC::MCMC()
    : pImpl_(std::make_unique<Impl>(
          [this]()
          { this->init_model_for_inference(); },
          [this]()
          { this->init_model_for_simulation(); },
          [this](Household &h) -> double
          { return this->compute_household_obs_log_lik(h); },
          [this](Household &h)
          { return this->simulate_single_household(h); },
          [this](Household &h, Individual &i)
          { this->update_cache_dependecies(h, i); })) {}
MCMC::~MCMC() = default;

///////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////  MCMC : general methods  /////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

double MCMC::uniform() { return pImpl_->setup->unif(); }
double MCMC::par(const std::string &p) const { return pImpl_->param.get(p); }

size_t MCMC::get_max_hh_cases() const { return pImpl_->data_wide->report_hh_size->max_cases.first; }
void MCMC::assert_no_missing_status() const
{
    if (pImpl_->data_wide->report_indiv->missing > 0)
        hsp::msg::error(
            "Individual with missing status are incompatible with current model.");
}

///////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////  MCMC : contact rates ///////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/// Build all group rates and contact rates and set rate per group. // UPDATE
void MCMC::build_all_rates()
{
    try
    {
        for (Household &house : *pImpl_->data_wide)
            init_contact_rates_(house);
    }
    catch (const std::exception &e)
    {
        hsp::msg::error("While initializing contact rates:\n" +
                        std::string(e.what()));
    }
}

/// Set contact rate according to user assignement. // UPDATE
void MCMC::init_contact_rates_(Household &house)
{
    // loop over all couples
    for (const Individual &member_k : house)     // from : member_k
        for (const Individual &member_i : house) // -> to : member_i
            if (member_k != member_i)
            {
                _pair_contact_rate = hsp::na::STRING;
                init_household_contact_rates(house, member_k, member_i);

                if (!hsp::na::is(_pair_contact_rate))
                {
                    // default reference rate
                    Attribute y(hsp::rates::VALUE_REFERENCE);

                    if (_pair_contact_rate == hsp::rates::MISSING)
                        throw std::invalid_argument(
                            "Missing contact rate for individuals '" +
                            member_k.id() + "' and '" + member_i.id() +
                            "' (" + house.err_name() + ").");

                    else if (_pair_contact_rate != hsp::rates::REFERENCE)
                        y = pImpl_->param.at(_pair_contact_rate); // parameter rate

                    house.init_contact_rates(static_cast<size_t>(member_k), static_cast<size_t>(member_i), y);
                }
            }
}
