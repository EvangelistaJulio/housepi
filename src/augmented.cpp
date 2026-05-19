/* housepi - Epidemiological Analysis of Household Outbreaks.
 *
 * Copyright © 2024-2026 Institut Pasteur (Paris), CNRS and INSERM.
 * Distributed under the terms of the GNU General Public License (GPLv3).
 * See the LICENSE file for details.
 */

#include <algorithm>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

#include "tools/dependency.h"
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
#include "mcmc_base.h"
#include "data_wide.h"
#include "distribution.h"
#include "function.h"
#include "household.h"
#include "individual.h"
#include "infec_profile.h"
#include "parameter.h"
#include "variable.h"

namespace housepi
{
    /////////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////  *** InfecProfile : struct ***  //////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////

    // UPDATE: this could be just forward declared maybe
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

        if (is_discrete_requested_ ||
            (input.count("discr") && str::tobool(input.at("discr")[0])))
            output->init_discretization();
        return output;
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////  *** Augmented : class ***  //////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////

    Augmented::Augmented(MCMCBase *mcmc)
        : mcmc_(mcmc) {}
    Augmented::~Augmented() = default;

    void Augmented::build(InfecProfile &infec_profile, const str_table &input)
    {
        // drop previous definitions
        is_any_sample_ = false;
        is_any_model_ = false;

        khold_.clear_optional();
        map_status_to_def_row_.clear();
        map_status_to_block_.clear();
        status_blocks_.clear();
        list_export.clear();
        indiv_augm_.clear();

        // skip if no augmented definitions
        if (input.empty() && !infec_profile.is_augmented())
            return;

        helper_build_allocate_storage_();                  // allocate and plug-in
        helper_build_assert_augm_input_(input);            // check entries
        helper_build_append_infec_profile_(infec_profile); // add infectivity profile
        helper_build_definition_by_status_();              // process definition by status
        helper_build_set_allocation_order_();              // allocation and order
    }

    void Augmented::helper_build_allocate_storage_()
    {
        indiv_augm_.resize(mcmc_->data_wide().report_indiv.total);
        for (Household &house : mcmc_->data_wide())
            for (Individual &member : house)
                member.paugm_ = &indiv_augm_[member.index_];
    }

    /// Check augmented names and assert no reserved names are used.
    void Augmented::helper_build_assert_augm_input_(const str_table &input)
    {
        const std::vector<std::string> file_headers = {"name", "status",
                                                       "func", "fcumul",
                                                       "value", "const", "discr",
                                                       "model", "min", "max",
                                                       "proposal", "proposal_scale"};
        const std::vector<std::string> headers = read::extract_keys(input);
        str::assert_no_missing_entries(headers, file_headers);
        str::assert_no_unused_entries(headers, file_headers);

        for (const std::string &astr : input.at("name"))
        {
            if (na::is(astr))
                throw std::invalid_argument(
                    "Found an augmented variable with undefined name.");
            if (astr == infprof::AUGM_FUNC)
                throw std::invalid_argument(
                    "Tried to use the reserved augmented name '" +
                    astr + "'.\nPlease, choose a different one.");
        }

        deftable_ = input; // store user definitions
    }

    void Augmented::helper_build_append_infec_profile_(InfecProfile &infec_profile)
    {
        if (!infec_profile.is_augmented())
            return;

        for (const std::string &col : read::extract_keys(deftable_))
            deftable_[col] += infec_profile.augm_definitions.at(col);

        khold_.add(infprof::AUGM_FUNC);
        infec_profile.augm_key = khold_.from_string(infprof::AUGM_FUNC);

        // append discrete infectivity profile definitions if any
        if (infec_profile.is_discrete_requested())
            user_discrete_set_.emplace(infec_profile.augm_key);
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////  Augmented : assert univocal definitions  //////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////

    void Augmented::helper_build_definition_by_status_()
    {
        for (size_t k = 0; k < deftable_["name"].size(); ++k)
        {
            // fill registry with optional augmented definitions
            const std::string astr = deftable_["name"][k];
            khold_.add(astr);

            // check associated status
            int status = status::MISSING;
            const std::string &s = deftable_["status"][k];
            if (!na::is(s))
            {
                status = str::toi(s, "augmented status");
                if (!status::is_infected(status))
                    throw std::invalid_argument(
                        "Impossible to associate augmented '" + astr +
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

    void Augmented::helper_build_set_allocation_order_()
    {
        std::unordered_set<int> all_infected_statuses = {};

        // NOTE: storage already allocated and pointers set by helper_build_allocate_storage_()
        size_t n = khold_.size();
        for (Household &house : mcmc_->data_wide())
            for (Individual &member : house)
            {
                Indiv::Augm &iaugm = *member.paugm_;
                iaugm.vec.resize(n);
                for (size_t a = 0; a < n; ++a)
                {
                    iaugm.vec[a].name(khold_.as_string(a));
                    iaugm.vec[a].set_err_type("augmented");
                }
                if (member.is_infected_status())
                    all_infected_statuses.emplace(member.status());
            }

        // first, load status specific user definitions
        // NOTE: NA (i. e. generic) does not necessarily need checking
        for (auto &[status, vec] : map_status_to_def_row_)
            if (status::is_infected(status))
                set_definition_dependency_order_for_(status, status);

        // associate observed status to corresponding status block
        // if not defined specifically, then associate to generic
        for (int status : all_infected_statuses)
            if (!map_status_to_block_.count(status))
            {
                if (!map_status_to_block_.count(status::MISSING))
                    set_definition_dependency_order_for_(status::MISSING, status);
                map_status_to_block_.emplace(status, map_status_to_block_.at(status::MISSING));
            }

        // if no observed infected status or status-specific definitions
        // then load generic definitition
        if (all_infected_statuses.empty() ||
            (map_status_to_def_row_.size() == 1 &&
             status::is_missing(map_status_to_def_row_.begin()->first)))
            set_definition_dependency_order_for_(status::MISSING, status::MISSING);
    }

    void Augmented::set_definition_dependency_order_for_(int status, int err_status)
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
            // NOTE: but ensure private infectivity is stored.
            for (const std::string &astr : dependency::order_dependencies(tree))
                if (astr == infprof::AUGM_FUNC ||
                    !str::is_private(astr) ||
                    !str::drop_private(dependency::order_conditioned_on(astr, tree)).empty())
                    status_blocks_.back().order.emplace_back(name_to_row[astr]);

            // define usage lists respecting definition order
            set_sample_list_(status);
        }
        catch (const std::exception &e)
        {
            std::string sms = e.what();
            if (status::is_missing(err_status))
                sms += "\nCheck your generic definitions.";
            else
                sms += " for status " + std::to_string(err_status) + ".";
            throw std::invalid_argument(sms);
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

    ///////////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////  Augmented : sample lists  //////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////

    void Augmented::set_sample_list_(int status)
    {
        StatusBlock_ &block = status_blocks_[map_status_to_block_.at(status)];
        for (size_t k : block.order)
            set_sample_list_for_single_defininition_(status, k);

        // NOTE: by construction of draws list, I am forced to drop private at the end
        list_export = khold_.from_string(str::drop_private(khold_.as_string(list_export)));
        block.draws = khold_.from_string(str::drop_private(khold_.as_string(block.draws)));
    }

    void Augmented::set_sample_list_for_single_defininition_(int status, size_t k)
    {
        const std::string &astr = deftable_["name"][k];
        size_t a = khold_.from_string(astr);

        bool is_to_export = true;
        bool is_to_draw = true;

        if (str::tobool(deftable_["const"][k]))
        {
            is_to_draw = false;
        }
        else if (!(na::is(deftable_["func"][k]) && na::is(deftable_["value"][k])))
        {
            auto [is_to_export_, is_to_draw_] = set_sample_list_for_function_(status, k);
            is_to_export = is_to_export_;
            is_to_draw = is_to_draw_;
        }
        else if (na::is(deftable_["value"][k]) && na::is(deftable_["model"][k]))
        {
            throw std::invalid_argument(
                "Impossible to initialize augmented '" + astr + "'  without a value or model");
        }
        else if (na::is(deftable_["proposal"][k]))
        {
            throw std::runtime_error(
                "No proposal was provided for non-constant definition of augmented '" +
                astr + "'");
        }
        else
        {
            is_any_sample_ = true;
            status_blocks_[map_status_to_block_.at(status)].sample.emplace_back(a);
        }

        // nothing to do if not to export
        if (!is_to_export)
            return;

        if (std::find(list_export.begin(), list_export.end(), a) == list_export.end())
            list_export.emplace_back(a);

        if (is_to_draw)
            status_blocks_[map_status_to_block_.at(status)].draws.emplace_back(a);
    }

    std::pair<bool, bool> Augmented::set_sample_list_for_function_(int status, size_t k)
    {
        bool is_to_export = true;
        bool is_to_draw = true;

        // NOTE: value definition has priority over function
        const std::string &str_val = deftable_["value"][k];
        const parse::Definition parser(na::is(str_val) ? deftable_["func"][k] : str_val);

        // Set potentially frequency data
        mcmc_->data_wide().set_frequency_covariates(parser.list_data_wide_depen());

        // NOTE: do not export if not callable (otherwise crashes)
        if (!parser.is_callable_implicitly("x"))
            is_to_export = false;

        bool all_consts = true; // start assuming all constants
        for (const std::string &pstr : parser.list_parameter_depen())
        {
            size_t p = mcmc_->param.khold().from_string(pstr);
            all_consts &= mcmc_->param.at(p).is_constant();
        }

        // inspect status specific augmented definitions and check if to draw (i. e. not const)
        // NOTE: taking advatange of the order!
        for (const std::string &bstr : parser.list_augmented_depen())
        {
            const std::vector<size_t> &v = status_blocks_[map_status_to_block_.at(status)].draws;
            if (std::find(v.begin(), v.end(), khold_.from_string(bstr)) != v.end())
                all_consts = false;
        }

        // to draw if not all constant dependencies
        is_to_draw = !all_consts;

        return {is_to_export, is_to_draw};
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////  Augmented : initialize augmented  //////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////

    void Augmented::init()
    {
        if (khold_.empty())
            return;

        helper_init_model_list_();

        // define augmented variables for all augmented individuals in ordefr
        for (Household &house : mcmc_->data_wide())
            for (size_t i : house.infected_hhindexes())
            {
                Individual &member = house.member(i);
                for (size_t k : helper_materialize_definition_order_(member))
                    helper_init_attributes_(member, k);
            }
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////  Augmented : init augmented models  /////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////

    void Augmented::helper_init_model_list_()
    {
        for (StatusBlock_ &block : status_blocks_)
        {
            block.depen.resize(khold_.size());

            // associate blocks for ordered model lists
            const str_table &ordered_models = get_ordered_models_(block);
            if (!ordered_models.empty())
            {
                is_any_model_ = true;
                block.models = khold_.from_string(read::extract_keys(ordered_models));

                for (const auto &[astr, others] : ordered_models)
                    block.depen[khold_.from_string(astr)] = khold_.from_string(others);
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
            ordered_models.emplace(key, dependency::order_conditioned_on(key, tree));

        return ordered_models;
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////  Augmented : init missing attribute  ////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////

    void Augmented::clear_all()
    {
        if (khold_.empty())
            return;

        for (Household &house : mcmc_->data_wide())
            for (size_t i : house.infected_hhindexes())
                clear_definitions(house.member(i));
    }

    void Augmented::clear_definitions(Individual &member)
    {
        for (Variable &xaugm : member.paugm_->vec)
            xaugm.clear();
    }

    void Augmented::init_missing_definitions(Individual &member)
    {
        for (size_t k : helper_materialize_definition_order_(member))
            if (member.paugm_->vec[khold_.from_string(deftable_["name"][k])].is_empty())
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
            member.paugm_->idx = i;
            return status_blocks_[i].order;
        }

        // otherwise if it exists fallback to generic definition
        const auto generic_it = map_status_to_block_.find(status::MISSING);
        if (generic_it != map_status_to_block_.end())
        {
            const size_t i = generic_it->second;
            map_status_to_block_.emplace(status, i);
            member.paugm_->idx = i;
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
        const std::string &astr = deftable_["name"][k];
        size_t a = khold_.from_string(astr);

        try
        {
            Variable &xaugm = member.paugm_->vec[a];

            if (str::tobool(deftable_["const"][k]))
            {
                xaugm.init_numeric(
                    mcmc_->param.read_numeric(deftable_["value"][k], &member).value());
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
                "While initializing augmented '" + astr + "' for " +
                member.err_name() + " (status " + std::to_string(member.status()) + "):\n" +
                std::string(e.what()));
        }
    }

    void Augmented::helper_init_model_attributes_(Individual &member, size_t k)
    {
        const std::string &astr = deftable_["name"][k];
        Variable &xaugm = member.paugm_->vec[khold_.from_string(astr)];
        parse::Definition parser(deftable_["model"][k]);

        xaugm.init_distr(parser.root.func);

        if (!xaugm.has_model())
            return;

        for (const auto &child : parser.root.children)
            xaugm.init_distr_attr(child.key, mcmc_->param.read_attribute(child, &member));

        for (const auto key : {"min", "max"})
            xaugm.init_distr_attr(key, mcmc_->param.read_attribute(deftable_[key][k], &member));

        xaugm.init_distr_reshaping(str::tobool(deftable_["discr"][k]));
    }

    void Augmented::helper_init_func_attributes_(Individual &member, size_t k)
    {
        const std::string &astr = deftable_["name"][k];
        Variable &xaugm = member.paugm_->vec[khold_.from_string(astr)];
        std::string str_val = deftable_["value"][k];
        parse::Definition parser(na::is(str_val) ? deftable_["func"][k] : str_val);

        // pointer mode
        if (parser.is_terminal())
        {
            xaugm.transfer(mcmc_->param.read_attribute(parser.root, &member));
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
            xaugm.init_func_attr(child.key, mcmc_->param.read_attribute(child, &member));
            mcmc_->param.make_cache_dependent(child);
        }

        // range
        for (const auto key : {"min", "max"})
            xaugm.init_func_attr(key, mcmc_->param.read_attribute(deftable_[key][k], &member));

        // discretization and truncation
        xaugm.init_func_reshaping(str::tobool(deftable_["discr"][k]));
    }

    /////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////  Augmented : user reset  ////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////

    void Augmented::user_reset(rnd_gen &prng, Individual &member)
    {
        for (size_t a : list_sample(member))
            member.paugm_->vec[a].user_reset(prng);

        for (Variable &xaugm : member.paugm_->vec)
            xaugm.optimize();
    }

    void Augmented::user_reset(rnd_gen &prng)
    {
        if (khold_.empty())
            return;

        for (Household &house : mcmc_->data_wide())
            for (size_t i : house.infected_hhindexes())
                user_reset(prng, house.member(i));
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////  Augmented : log-likelihood  //////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////

    const std::vector<size_t> &Augmented::list_sample(const Individual &member) const
    {
        return status_blocks_[member.paugm_->idx].sample;
    }
    const std::vector<size_t> &Augmented::list_draws(const Individual &member) const
    {
        return status_blocks_[member.paugm_->idx].draws;
    }
    const std::vector<size_t> &Augmented::list_models(const Individual &member) const
    {
        return status_blocks_[member.paugm_->idx].models;
    }

    std::string format_exception_for(Individual &member, const std::exception &e)
    {
        std::string sms = e.what();
        if (!sms.empty() && sms.back() == '.')
            sms.pop_back();
        return sms + " for " + member.err_name() + ".";
    }

    double Augmented::commit_augm_model_log_lik(Individual &member)
    {
        try
        {
            const std::vector<size_t> &models = list_models(member);

            if (models.empty())
                return 0.;
            double new_ll = 0.0;
            for (size_t a : models)
                new_ll += member.paugm_->vec[a].get_distr_log_pdf();
            member.paugm_->log_lik.set(new_ll);
            return new_ll;
        }
        catch (const std::exception &e)
        {
            throw std::runtime_error(format_exception_for(member, e));
        }
    }

    double Augmented::commit_augm_model_log_ratio(Individual &member, size_t a)
    {
        try
        {
            const StatusBlock_ &block = status_blocks_[member.paugm_->idx];
            if (block.depen[a].empty())
                return 0.;
            double delta_ll = 0;
            // NOTE: vector of dependencies 'depen' includes also 'a'
            for (size_t b : block.depen[a])
                delta_ll += member.paugm_->vec[b].get_distr_log_ratio();
            member.paugm_->log_lik.increment(delta_ll);
            return delta_ll;
        }
        catch (const std::exception &e)
        {
            throw std::runtime_error(format_exception_for(member, e));
        }
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////  Augmented : sample draws  ///////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////

    std::vector<const Individual *> Augmented::init_individual_sample(
        const std::vector<std::string> &id_draws) const
    {
        std::vector<const Individual *> output;

        if (id_draws.empty())
            return output;

        if (!is_any_sample())
        {
            msg::warning(
                "Ignored id draws. No augmented variables to sample were defined.");
            return output;
        }

        output.reserve(id_draws.size());
        std::vector<std::string> existing_ids = {};
        for (const Household &house : mcmc_->data_wide())
            for (const Individual &member : house)
            {
                if (std::find(id_draws.begin(), id_draws.end(), member.id()) == id_draws.end())
                    continue;

                // unique IDs are granted by class DataWide load
                existing_ids.push_back(member.id());

                // check status of chosen individual
                if (member.is_infected_status())
                    output.emplace_back(&member);
                else
                    msg::warning(
                        "Ignored draws of non infected " + member.err_name() + ".");
            }
        output.shrink_to_fit();

        // check all user request are satisfied
        str::assert_no_missing_entries(existing_ids, id_draws, "draws id");
        return output;
    }
}
