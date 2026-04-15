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

#include "tools/keys.h"
#include "tools/parser.h"
#include "tools/read.h"
#include "tools/string.h"
#include "tools/types.h"

#include "infec_profile.h"
#include "individual.h"
#include "data_wide.h"
#include "distribution.h"
#include "function.h"
#include "parameter.h"

/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////  *** InfecProfile : struct ***  //////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

namespace housepi
{
    namespace infprof
    {
        void LookUpTable::init(const Distribution *ptr, size_t n)
        {
            generic_ptr_ = ptr;
            size_ = n;
            value_.assign(n, 0);
            old_.assign(n, 0);
            materialize_();
        }
        void LookUpTable::clear() noexcept
        {
            value_.clear();
            old_.clear();
        }
        void LookUpTable::recover() noexcept
        {
            if (generic_ptr_->is_constant())
                return;
            std::swap(value_, old_);
        }
        void LookUpTable::update()
        {
            if (generic_ptr_->is_constant())
                return;
            std::swap(value_, old_);
            materialize_();
        }

        void LookUpTable::Table_::clear() noexcept
        {
            inst.clear();
            cumul.clear();
        }
        void LookUpTable::Table_::assign(size_t n, double v)
        {
            inst.assign(n, v);
            cumul.assign(n, v);
        }

        void LookUpTable::materialize_()
        {
            for (size_t i = 1; i < size_; ++i)
            {
                double x = generic_ptr_->pdf(i);
                value_.inst[i] = x;
                value_.cumul[i] = x + value_.cumul[i - 1];
            }
        }
    }

    InfecProfile::InfecProfile() = default;
    InfecProfile::~InfecProfile() = default;

    void InfecProfile::init_template(Parameter &param, str_table input)
    {
        // if not previous build, proceed with template
        // if no template, assume no need for an infectivity profile
        if (generic || !augm_definitions.empty() || input.empty())
            return;

        try
        {
            if (!input.count("func"))
                throw std::invalid_argument(
                    "A function name was missing for infectivity profile template.");

            size_t n_rows = read::csv_n_rows(input);
            const auto file_headers = {"func", "status", "discr", "max"};

            for (const auto &key : file_headers)
                if (!input.count(key))
                    input.emplace(key, std::vector<std::string>(n_rows, ""));

            // discretization
            input["discr"] = std::vector<std::string>(n_rows, "false");
            if (is_discrete_requested_)
                input["discr"] = std::vector<std::string>(n_rows, "true");

            else if (input.count("discr"))
                for (size_t r = 0; r < n_rows; ++r)
                    if (!na::is(input.at("discr")[r]))
                        input["discr"][r] = input.at("discr")[r];

            build(param, input); // run build as if external input

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

    void InfecProfile::init_template(Parameter &param, dict<std::string> input)
    {
        str_table new_input = {};
        for (const auto &[key, val] : input)
            new_input[key] = {val};
        init_template(param, std::move(new_input));
    }

    void InfecProfile::build(Parameter &param, const str_table &input)
    {
        augm_definitions.clear();

        if (input.empty())
            return;

        // check all entries are legitimiate
        str::assert_no_unused_entries(
            read::extract_keys(input),
            {"func", "status", "discr", "min", "max"});

        // check definition list has a function
        if (!input.count("func"))
            throw std::invalid_argument(
                "A function name choice was missing in infectivity profile definition.");
        for (const std::string &in_name : input.at("func"))
            if (parse::Definition(in_name).is_na())
                throw std::invalid_argument(
                    "Provided NA as a function name in infectivity profile definition.");

        // check if definition is generic or augmented
        if (read::csv_n_rows(input) > 1 ||
            parse::Definition(input.at("func")[0]).is_augmented_depen() ||
            parse::Definition(input.at("max")[0]).is_augmented_depen() ||
            !na::is(input.at("status")[0]))
            build_augmented_(param, input);
        else
            build_generic_(param, input);
    }

    void InfecProfile::build_generic_(Parameter &param, const str_table &input)
    {
        mode_ = Mode::Generic;
        generic = helper_build_(param,
            {{"func", {input.at("func")[0]}},
             {"discr", {input.count("discr") ? input.at("discr")[0] : "false"}}});
        lookup_.clear();
        augm_definitions.clear();

        // range minimum
        std::string min_expr = generic->is_discrete() ? "1|N" : "0|N"; // default
        if (input.count("min"))
            min_expr = input.at("min")[0];
        generic->init_attr("min", param.read_attribute(min_expr));

        // range maximum
        if (input.count("max"))
            generic->init_attr("max", param.read_attribute(input.at("max")[0]));

        // range truncation
        generic->init_truncation();

        // if discrete and generic, initializes lookup table
        if (generic->is_discrete())
        {
            mode_ = Mode::LookUp;
            size_t n = infprof::LOOKUP_SIZE;
            if (!generic->upper().is_empty())
                n = 1 + generic->upper().value();
            lookup_.init(generic.get(), n);
        }
    }

    void InfecProfile::build_augmented_(Parameter &param, const str_table &input)
    {
        // drop generic definition
        mode_ = Mode::Augmented;
        generic = nullptr;
        lookup_.clear();
        augm_definitions.clear();

        // use helper to check profile name and discretization
        for (size_t r = 0; r < read::csv_n_rows(input); ++r)
            helper_build_(param, {{"func", {parse::Definition(input.at("func")[r]).root.func}},
                                  {"discr", {input.at("discr")[r]}}});

        // match infectivity profile and augmented variable definition
        for (size_t r = 0; r < read::csv_n_rows(input); ++r)
        {
            bool is_discrete = str::tobool(input.at("discr")[r]);

            augm_definitions["name"].push_back(infprof::AUGM_FUNC);
            augm_definitions["status"].push_back(input.at("status")[r]);
            augm_definitions["func"].push_back(input.at("func")[r]);
            augm_definitions["fcumul"].push_back("false");
            augm_definitions["value"].push_back("");
            augm_definitions["const"].push_back("false");
            augm_definitions["discr"].push_back(str::format(is_discrete));
            augm_definitions["model"].push_back("");
            augm_definitions["min"].push_back(is_discrete ? "1|N" : "0|N");
            augm_definitions["max"].push_back(input.at("max")[r]);
            augm_definitions["proposal"].push_back("");
            augm_definitions["proposal_scale"].push_back("");
        }
    }

    std::unique_ptr<Distribution> InfecProfile::helper_build_(Parameter &param, const str_table &input)
    {
        parse::Definition parser(parse::Definition(input.at("func")[0]));
        auto output = compose_distr_registries(parser.root.func,
                                               "infectivity profile",
                                               {"infecprofile", "continuos", "discrete"});
        for (const auto &node : parser.root.children)
        {
            output->init_attr(node.key, param.read_attribute(node));
            param.make_cache_dependent(node);
        }

        // discretization
        if (input.count("discr") && str::tobool(input.at("discr")[0]))
            output->init_discretization();
        else if (is_discrete_requested_)
            throw std::invalid_argument(
                "Definition of infectivity profile did not return a discrete object as requested by the model.");
        return output;
    }
}
