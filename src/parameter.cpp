/* housepi - Epidemiological Analysis of Household Outbreaks.
 *
 * Copyright © 2026 Institut Pasteur (Paris), CNRS and INSERM.
 * Distributed under the terms of the GNU General Public License (GPLv3).
 * See the LICENSE file for details.
 */

#include <algorithm>
#include <stdexcept>
#include <string>
#include <vector>

#include "tools/dependency.h"
#include "tools/keys.h"
#include "tools/messages.h"
#include "tools/nan.h"
#include "tools/operators.h"
#include "tools/parser.h"
#include "tools/random.h"
#include "tools/read.h"
#include "tools/string.h"
#include "tools/types.h"

#include "parameter.h"
#include "data_wide.h"
#include "function.h"
#include "group_rate.h"
#include "individual.h"
#include "variable.h"

/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////  *** Parameter : struct ***  ////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

namespace housepi
{
    void Parameter::build(const str_table &input)
    {
        khold_.clear_optional();

        if (input.empty())
            return;

        // check entries
        const std::vector<std::string> file_headers = {"name", "value", "const", "discr",
                                                       "func",
                                                       "prior", "min", "max",
                                                       "proposal", "proposal_scale"};
        const std::vector<std::string> headers = read::extract_keys(input);
        str::assert_no_missing_entries(headers, file_headers);
        str::assert_no_unused_entries(headers, file_headers);

        // assert that parameter is not already defined
        if (input.at("name").size() != str::unique(input.at("name")).size())
            throw std::invalid_argument(
                "Found multiple definitions for the same parameter.");

        // build dependency graph and validate function parameters
        bool has_func_params = false;
        dict<size_t> name_to_row;
        for (size_t k = 0; k < input.at("name").size(); ++k)
        {
            const std::string &pstr = input.at("name")[k];
            name_to_row[pstr] = k;

            // validate: function parameters can only depend on other parameters
            if (!na::is(input.at("func")[k]))
            {
                has_func_params = true;
                parse::Definition parser(input.at("func")[k]);
                if (parser.is_augmented_depen())
                    throw std::invalid_argument(
                        "Parameter function '" + pstr +
                        "' cannot depend on augmented variables.");
                if (parser.is_data_wide_depen())
                    throw std::invalid_argument(
                        "Parameter function '" + pstr +
                        "' cannot depend on data variables.");
            }
        }

        // topological sort only when function parameters exist (preserves input order otherwise)
        std::vector<std::string> ordered;
        if (has_func_params)
        {
            // use input order as iteration base for stable topological sort
            str_table dep_graph;
            for (size_t k = 0; k < input.at("name").size(); ++k)
            {
                const std::string &pstr = input.at("name")[k];
                std::vector<std::string> deps;
                for (const auto col : {"func", "min", "max"})
                    deps += parse::Definition(input.at(col)[k]).list_parameter_depen();
                dep_graph[pstr] = str::unique(deps);
            }
            ordered = dependency::order_dependencies(dep_graph);
        }
        else
        {
            ordered = input.at("name");
        }

        init_param_vec_();
        clear_lists_();
        // load parameters in dependency order
        for (const std::string &pstr : ordered)
        {
            size_t k = name_to_row.at(pstr);

            if (!khold_.count(pstr))
            {
                khold_.add(pstr);
                param_vec_.emplace_back();
                param_vec_.back().name(pstr);
                param_vec_.back().set_err_type("parameter");
            }

            if (!na::is(input.at("func")[k]))
                helper_build_function_(input, k);
            else if (str::tobool(input.at("const")[k]))
                helper_build_constant_(input, k);
            else
                helper_build_to_sample_(input, k);
        }

        // initialize all parameters to standard
        list_sample_standard = list_sample;
        for (size_t p : list_sample_standard)
            if (at(p).is_adaptable())
                list_adaptable_sample_standard.push_back(p);
            else
                list_non_adapt_sample_standard.push_back(p);
    }

    void Parameter::build_for_predic(const dbl_table &input)
    {
        if (input.empty())
            return;

        const std::vector<std::string> headers = read::extract_keys(input);
        // assert that parameter is not already defined
        if (headers.size() != str::unique(headers).size())
            throw std::invalid_argument(
                "Found multiple definitions for the same parameters.");

        clear_lists_();
        simulation_values.clear();
        // update list of constant parameters to all parameters
        list_const = khold_.from_string(khold_.get_all());

        // reset definition dictionary to all constants NAs.
        for (const std::string &pstr : khold_.get_all())
            if (std::find(headers.begin(), headers.end(), pstr) == headers.end())
                simulation_values.emplace_back();
            else
                simulation_values.emplace_back(input.at(pstr));
    }

    void Parameter::init_param_vec_()
    {
        size_t n = khold_.size();
        param_vec_.clear();
        param_vec_.resize(n);
        for (size_t p = 0; p < n; ++p)
        {
            param_vec_[p].name(khold_.as_string(p));
            param_vec_[p].set_err_type("parameter");
        }
    }

    void Parameter::clear_lists_()
    {
        list_const.clear();
        list_func.clear();
        list_sample.clear();
        list_sample_standard.clear();
        list_adaptable_sample_standard.clear();
        list_non_adapt_sample_standard.clear();
        list_sample_cached.clear();
        list_adaptable_sample_cached.clear();
        list_non_adapt_sample_cached.clear();

        simulation_values.clear();
    }

    void Parameter::helper_build_constant_(const str_table &input, size_t k)
    {
        const std::string &pstr = input.at("name")[k]; // parameter name

        if (na::is(input.at("value")[k]))
            throw std::invalid_argument(
                "Constant parameter '" + pstr + "' had undefined value.");

        // assign constant parameter to the set of variables
        size_t p = khold_.from_string(pstr);
        param_vec_[p].transfer(read_numeric(input.at("value")[k]));
        // check if discrete as requested
        if (user_discrete_set_.count(p) && !param_vec_[p].is_discrete())
            throw std::invalid_argument(
                "Definition of " + param_vec_[p].err_name() +
                " did not return a discrete object as requested by model.");

        // update the lists of constant parameters
        list_const.push_back(p);
    }

    void Parameter::helper_build_to_sample_(const str_table &input, size_t k)
    {
        const std::string &pstr = input.at("name")[k]; // parameter name

        // assert parameter to sample has a valid value initialization
        if (na::is(input.at("value")[k]) && na::is(input.at("prior")[k]))
            throw std::invalid_argument(
                "Non-constant parameter '" + pstr + "' had undefined prior.");

        // create a new parameter Variable
        size_t p = khold_.from_string(pstr);
        Variable &xparam = at(p);
        xparam.transfer(read_numeric(input.at("value")[k]));

        // assign a prior distribution to the parameter
        parse::Definition parser(input.at("prior")[k]);
        xparam.init_distr(parser.root.func);

        // load prior arguments
        for (const auto &child : parser.root.children)
            xparam.init_distr_attr(child.key, read_numeric(child));
        // load parameter range
        for (const auto key : {"min", "max"})
            xparam.init_distr_attr(key, read_numeric(input.at(key)[k]));
        xparam.init_distr_reshaping(str::tobool(input.at("discr")[k]));

        // assign proposal
        xparam.init_proposal(input.at("proposal")[k],
                             str::tod(input.at("proposal_scale")[k]));

        // check if discrete as requested
        if (user_discrete_set_.count(p) && !xparam.is_discrete())
            throw std::invalid_argument(
                "Definition of " + xparam.err_name() +
                " did not return a discrete object as requested by model.");

        // update the lists of parameters to sample
        list_sample.push_back(p);
    }

    void Parameter::helper_build_function_(const str_table &input, size_t k)
    {
        const std::string &pstr = input.at("name")[k];
        size_t p = khold_.from_string(pstr);
        Variable &xparam = at(p);

        parse::Definition parser(input.at("func")[k]);

        // pointer mode (terminal = alias for another parameter)
        if (parser.is_terminal())
        {
            xparam.transfer(read_attribute(parser.root));
            list_const.push_back(p);
            list_func.push_back(p);
            return;
        }

        // function mode
        xparam.init_function(compose_func_registries(
            parser.root.func, "function",
            {"custom", "transcendental", "arithmetic"}));

        // bind function arguments
        for (const auto &child : parser.root.children)
            xparam.init_func_attr(child.key, read_attribute(child));

        // range
        for (const auto key : {"min", "max"})
            xparam.init_func_attr(key, read_attribute(input.at(key)[k]));

        // discretization
        xparam.init_func_reshaping(str::tobool(input.at("discr")[k]));

        // check if discrete as requested
        if (user_discrete_set_.count(p) && !xparam.is_discrete())
            throw std::invalid_argument(
                "Definition of " + xparam.err_name() +
                " did not return a discrete object as requested by model.");

        // classify: deterministic, not sampled
        list_const.push_back(p);
        list_func.push_back(p);
    }

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

    void Parameter::user_reset(rnd_gen &prng)
    {
        for (size_t p : list_sample)
            at(p).user_reset(prng);

        for (size_t p = 0; p < khold_.size(); ++p)
            at(p).optimize();
    }

    void Parameter::load_sample_values(size_t r)
    {
        for (size_t p = 0; p < khold_.size(); ++p)
        {
            const std::vector<double> &vec = simulation_values[p];
            if (!vec.empty())
                at(p).set(vec[r]);
        }
    }

    void Parameter::make_cache_dependent(const parse::Entry &node)
    {
        // cache dependecy if an augmented function depends on a parameter to be inferred
        if (parse::is_parameter(node.type))
        {
            size_t p = khold_.from_string(node.value);
            if (!at(p).is_constant())
                make_cache_dependent(p);
        }
    }

    bool Parameter::propose_valid_move(size_t p, rnd_gen &prng)
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

    void Parameter::adapt_proposal_of_sampled(unsigned long step, unsigned long tot)
    {
        const double ratio = static_cast<double>(step) / static_cast<double>(tot);

        // adpat all standard model parameters
        for (size_t p : list_adaptable_sample_standard)
            at(p).adapt_proposal(ratio);

        // adapt all parameters with cached dependecies
        for (size_t p : list_adaptable_sample_cached)
            at(p).adapt_proposal(ratio);
    }

    const Attribute &Parameter::read_attribute(
        const parse::Entry &node,
        const Individual *member_ptr) const
    {
        storage_ = Attribute(); // reset auxiliary storage to default empty

        if (!node.is_terminal())
        {
            helper_store_read_function_(node, member_ptr);
        }
        else if (parse::is_numeric(node.type))
        {
            storage_.init_numeric(str::tod(node.value, "argument with numeric type"));
        }
        else if (parse::is_parameter(node.type))
        {
            return at(khold_.from_string(node.value));
        }
        else if (!node.is_empty())
        {
            if (!member_ptr)
            {
                std::string sms = "Individual specific definitions as " + node.err_str() +
                                  " are not supported here.";
                // add tip if augmented (which does not require '@')
                if (parse::is_augmented(node.type) && node.is_terminal())
                    sms += "\nDid you forget rlang injection operator `!!`?";
                throw std::invalid_argument(sms);
            }

            if (parse::is_data_wide(node.type))
            {
                size_t d = pDataWideKeys_->from_string(node.value);
                storage_.init_numeric(member_ptr->data(d));
                storage_.name(node.value);
                storage_.set_err_type("data");
            }
            else if (parse::is_augmented(node.type))
            {
                return member_ptr->get_augm_variable(pAugmKeys_->from_string(node.value));
            }
            else if (parse::is_group_rate(node.type))
            {
                return member_ptr->get_group_rate_attribute(pGroupRateKeys_->from_string(node.value));
            }
            else
                throw std::invalid_argument(
                    "Unrecognized source type '" + node.type +
                    "' found while parsing " + node.err_str() + ".");
        }

        return storage_;
    }

    void Parameter::helper_store_read_function_(
        const parse::Entry &node,
        const Individual *member_ptr) const
    {
        Attribute parent;
        parent.init_function(compose_func_registries(
            node.func, "function",
            {"custom", "transcendental", "arithmetic",
             "continuos", "discrete", "infecprofile"}));

        for (const auto &child : node.children)
            parent.init_func_attr(child.key, read_attribute(child, member_ptr));

        storage_ = parent;
    }

    const Attribute &Parameter::read_numeric(
        const parse::Entry &node,
        const Individual *member_ptr) const
    {
        storage_ = read_attribute(node, member_ptr);
        if (!storage_.is_numeric() && !storage_.is_empty())
            throw std::invalid_argument(
                "Definition " + node.err_str() +
                " could not be interpreted as a numeric constant.");
        return storage_;
    }
}
