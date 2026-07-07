/* housepi - Epidemiological Analysis of Household Outbreaks.
 *
 * Copyright © 2024-2026 Institut Pasteur (Paris), CNRS and INSERM.
 * Distributed under the terms of the GNU General Public License (GPLv3).
 * See the LICENSE file for details.
 */

#include <algorithm>
#include <cmath>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

#include "tools/dependency.h"
#include "tools/keys.h"
#include "tools/logproba.h"
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
#include "translator.h"

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

    void InfecProfile::init_template(Translator &translator, str_table input)
    {
        // if not previous build, proceed with template
        // if no template, assume no need for an infectivity profile
        if (generic || !augm_definitions.empty() || input.empty())
            return;

        try
        {
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

            build(translator, input); // run build as if external input

            // forbid augmented infectivity profile for template
            if (is_augmented())
                msg::stop(
                    "Impossible to define infectivity profile template with augmented dependencies.");

            // UPDATE: maybe also not depending on parameters?
        }
        catch (const std::exception &e)
        {
            msg::stop(
                std::string(e.what()) +
                "\nWhile building the infectivity profile template.");
        }
    }

    void InfecProfile::init_template(Translator &translator, dict<std::string> input)
    {
        str_table new_input = {};
        for (const auto &[key, val] : input)
            new_input[key] = {val};
        init_template(translator, std::move(new_input));
    }

    void InfecProfile::build(Translator &translator, const str_table &input)
    {
        augm_definitions.clear();

        if (input.empty())
            return;

        // check if definition is generic or augmented
        if (read::csv_n_rows(input) > 1 ||
            parse::Definition(input.at("func")[0]).is_augmented_depen() ||
            parse::Definition(input.at("max")[0]).is_augmented_depen() ||
            !na::is(input.at("status")[0]))
            build_augmented_(translator, input);
        else
            build_generic_(translator, input);
    }

    void InfecProfile::build_generic_(Translator &translator, const str_table &input)
    {
        mode_ = Mode::Generic;
        generic = helper_build_(translator,
                                {{"func", {input.at("func")[0]}},
                                 {"discr", {input.count("discr") ? input.at("discr")[0] : "false"}}});
        lookup_.clear();
        augm_definitions.clear();

        // range minimum
        std::string min_expr = generic->is_discrete() ? "1|N" : "0|N"; // default
        if (input.count("min"))
            min_expr = input.at("min")[0];
        generic->init_attr("min", translator.to_attribute(min_expr));

        // range maximum
        if (input.count("max"))
            generic->init_attr("max", translator.to_attribute(input.at("max")[0]));

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

    void InfecProfile::build_augmented_(Translator &translator, const str_table &input)
    {
        // drop generic definition
        mode_ = Mode::Augmented;
        generic = nullptr;
        lookup_.clear();
        augm_definitions.clear();

        // use helper to check profile name and discretization
        for (size_t r = 0; r < read::csv_n_rows(input); ++r)
            helper_build_(translator, {{"func", {parse::Definition(input.at("func")[r]).root.func}},
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
            augm_definitions["prior"].push_back("");
            augm_definitions["min"].push_back(is_discrete ? "1|N" : "0|N");
            augm_definitions["max"].push_back(input.at("max")[r]);
            augm_definitions["proposal"].push_back("");
            augm_definitions["proposal_scale"].push_back("");
            augm_definitions["kind"].push_back("augm");
        }
    }

    std::unique_ptr<Distribution> InfecProfile::helper_build_(Translator &translator, const str_table &input)
    {
        parse::Definition parser(parse::Definition(input.at("func")[0]));
        auto output = compose_distr_registries(parser.root.func,
                                               "infectivity profile",
                                               {"infecprofile", "continuos", "discrete"});
        for (const auto &node : parser.root.children)
        {
            output->init_attr(node.key, translator.to_attribute(node));
            translator.make_cache_dependent(node);
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

    /// Map a definition's kind label to its typed Kind.
    static Augmented::Kind kind_from_string_(const std::string &s)
    {
        if (s == "obs")
            return Augmented::Kind::Obs;
        if (s == "house")
            return Augmented::Kind::House;
        if (s == "status")
            return Augmented::Kind::Status;
        // if (s == "augm")
        return Augmented::Kind::Augm;
    }

    Augmented::Augmented(MCMCBase *mcmc)
        : mcmc_(mcmc) {}
    Augmented::~Augmented() = default;

    void Augmented::build(InfecProfile &infec_profile, const str_table &input)
    {
        // drop previous definitions
        is_any_sample_ = false;
        is_any_prior_ = false;
        is_any_obs_ = false;
        prior_param_depen_.clear();
        obs_param_depen_.clear();

        khold_.clear_optional();
        map_status_to_def_row_.clear();
        map_status_to_block_.clear();
        status_blocks_.clear();
        list_export.clear();
        list_obs_export.clear();
        list_house_export.clear();
        varia.clear();
        house_varia.clear();
        status_obs_def_ = na::INDEX;
        status_obs_key_ = na::INDEX;

        // skip if no augmented definitions
        if (input.empty() && !infec_profile.is_augmented())
            return;

        helper_build_assert_input_(input);    // check entries
        helper_build_append_(infec_profile);  // add infectivity profile
        helper_build_definition_by_status_(); // process definition by status
        helper_build_set_allocation_order_(); // allocation and order
        helper_init_observation_list_();      // observation nodes
        helper_init_house_obs_list_();        // household-level observations
        assert_status_obs_consistent_();      // assert status observation
        assert_no_obs_entry_collision_();     // no definition reads an obs entry column
    }

    /// Check augmented names and assert no reserved names are used.
    void Augmented::helper_build_assert_input_(const str_table &input)
    {
        // skip if no definitions
        if (input.empty())
            return;

        for (const std::string &astr : input.at("name"))
        {
            if (astr == infprof::AUGM_FUNC)
                msg::stop(
                    "Tried to use the reserved augmented name " +
                    msg::val(astr) + ".\nPlease, choose a different one.");
        }

        deftable_.clear();
        deftable_.append(input);
    }

    bool Augmented::is_obs_node(const std::string &name) const noexcept
    {
        for (size_t k = 0; k < deftable_.c_name.size(); ++k)
            if (deftable_.c_kind[k] == Kind::Obs && deftable_.c_name[k] == name)
                return true;
        return false;
    }

    void Augmented::set_prior_param_depen_()
    {
        // initialize list of parameter dependent augmented/observation
        prior_param_depen_.assign(mcmc_->param.khold().size(), false);
        obs_param_depen_.assign(mcmc_->param.khold().size(), false);

        // augmented nodes that (transitively) feed an observation, from the obs
        // inverse-dependency map: a parameter read by such a deterministic node
        // moves the observation through the node's value, not only the prior
        std::vector<bool> feeds_obs(khold_.size(), false);
        for (const StatusBlock_ &block : status_blocks_)
            for (size_t a = 0; a < block.obs_depen.size(); ++a)
                if (!block.obs_depen[a].empty())
                    feeds_obs[a] = true;

        for (size_t k = 0; k < deftable_.c_name.size(); ++k)
        {
            if (deftable_.c_name[k] == infprof::AUGM_FUNC ||
                deftable_.c_kind[k] == Kind::House ||
                deftable_.c_kind[k] == Kind::Status)
                continue;

            const bool is_obs = deftable_.c_kind[k] == Kind::Obs;
            // a deterministic node forwards its parameters to the observation it feeds
            const bool to_obs = is_obs ||
                (na::is(deftable_.c_prior[k]) && feeds_obs[khold(deftable_.c_name[k])]);

            for (const std::vector<std::string> *col : deftable_.source_cols())
                for (const std::string &pstr : parse::Definition((*col)[k]).list_parameter_depen())
                {
                    const size_t p = mcmc_->param.khold(pstr);
                    if (!is_obs)
                        prior_param_depen_[p] = true;
                    if (to_obs)
                        obs_param_depen_[p] = true;
                }
        }
    }

    void Augmented::DefTable_::clear()
    {
        for (auto col : {&c_name, &c_status, &c_func, &c_fcumul,
                         &c_value, &c_const, &c_discr, &c_prior,
                         &c_min, &c_max, &c_proposal, &c_scale})
            col->clear();
        c_kind.clear();
    }

    void Augmented::DefTable_::append(const str_table &input)
    {
        c_name += input.at("name");
        c_status += input.at("status");
        c_func += input.at("func");
        c_fcumul += input.at("fcumul");
        c_value += input.at("value");
        c_const += input.at("const");
        c_discr += input.at("discr");
        c_prior += input.at("prior");
        c_min += input.at("min");
        c_max += input.at("max");
        c_proposal += input.at("proposal");
        c_scale += input.at("proposal_scale");
        for (const std::string &s : input.at("kind"))
            c_kind.push_back(kind_from_string_(s));
    }

    void Augmented::helper_build_append_(InfecProfile &infec_profile)
    {
        if (!infec_profile.is_augmented())
            return;

        deftable_.append(infec_profile.augm_definitions);
        khold_.add(infprof::AUGM_FUNC);
        infec_profile.augm_key = khold(infprof::AUGM_FUNC);
        infec_profile.bind_augm(this); // for runtime per-member profile resolution

        // append discrete infectivity profile definitions if any
        if (infec_profile.is_discrete_requested())
            user_discrete_set_.emplace(infec_profile.augm_key);
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////  Augmented : assert univocal definitions  //////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////

    void Augmented::helper_build_definition_by_status_()
    {
        for (size_t k = 0; k < deftable_.c_name.size(); ++k)
        {
            // household-level observations are not per-member
            // (see helper_init_house_obs_list_ / init_house_for_simulation).
            if (deftable_.c_kind[k] == Kind::House)
                continue;

            // fill registry with optional augmented definitions
            const std::string astr = deftable_.c_name[k];
            khold_.add(astr);

            // status observation is per-member, but out of status block
            if (deftable_.c_kind[k] == Kind::Status)
            {
                status_obs_def_ = k;
                status_obs_key_ = khold(astr);
                continue;
            }

            // check associated status
            int status = status::MISSING;
            const std::string &s = deftable_.c_status[k];
            if (!na::is(s))
                status = str::toi(s, "augmented status");
            map_status_to_def_row_[status].push_back(k); // store status to definition row index
        }

        // NOTE: assertion of univocal name-status is done in R
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////  Augmented : build definition order  ////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////

    void Augmented::helper_build_set_allocation_order_()
    {
        std::unordered_set<int> all_infected_statuses = {};

        // allocate storage
        varia.resize(mcmc_->data_wide().report_indiv.total);

        size_t n = khold_.size();
        for (Household &house : mcmc_->data_wide())
            for (Individual &member : house)
            {
                IndivBlock &iaugm = varia_of_(member);
                iaugm.vec.resize(n);
                for (size_t a = 0; a < n; ++a)
                    iaugm.vec[a].name(
                        khold(a),
                        is_obs_node(khold(a)) ? ErrType::Observation : ErrType::Augmented);

                if (member.is_infected_status())
                    all_infected_statuses.emplace(member.status());
            }

        // no per-member definition to place into status blocks: a model carrying
        // only household-level observations needs no status block (they are handled
        // by helper_init_house_obs_list_). Skip block generation so such a model can
        // be built before its per-member augmented definitions are declared; a
        // genuinely required augmented variable is still enforced at the run wall.
        if (map_status_to_def_row_.empty())
            return;

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
            msg::stop(
                "No definitions were compatible with augmented status " +
                std::to_string(status) + ".");

        str_table tree;           // populate augmented dependency tree
        dict<size_t> name_to_row; // ... and reverse dictionary name-to-row-index
        for (size_t k : get_unsorted_definition_rows_(status))
        {
            std::vector<std::string> xdepen = {};
            for (const std::vector<std::string> *col : deftable_.source_cols())
                xdepen += parse::Definition((*col)[k]).list_augmented_depen();
            tree[deftable_.c_name[k]] = str::unique(xdepen);
            name_to_row[deftable_.c_name[k]] = k;
        }

        // check dependecies and translate order from name to deftable row index
        try
        {
            // NOTE: this method is just called once per status
            map_status_to_block_.emplace(status, status_blocks_.size());
            status_blocks_.emplace_back(status);

            // Ignore spurios private definition if no kept definition depends on it.
            // NOTE: infectivity profile (private) can depend on private
            auto is_kept_consumer = [](const std::string &c)
            { return !str::is_private(c) || c == infprof::AUGM_FUNC; };
            for (const std::string &astr : dependency::order_dependencies(tree))
            {
                const auto consumers = dependency::order_conditioned_on(astr, tree);
                if (astr == infprof::AUGM_FUNC || !str::is_private(astr) ||
                    std::any_of(consumers.begin(), consumers.end(), is_kept_consumer))
                    status_blocks_.back().order.emplace_back(name_to_row[astr]);
            }

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
            msg::stop(sms);
        }
    }

    std::vector<size_t> Augmented::get_unsorted_definition_rows_(int status) const
    {
        // load the list of status-specific definition rows
        std::vector<size_t> allrows = map_status_to_def_row_.at(status);

        // extract names
        std::vector<std::string> name_to_status = {};
        for (size_t k : allrows)
            name_to_status.push_back(deftable_.c_name[k]);

        // add generic definitions to current list (if not generic already)
        // skip generic if a status-specific already exist
        if (status::is_infected(status) &&
            map_status_to_def_row_.count(status::MISSING))
            for (size_t k : map_status_to_def_row_.at(status::MISSING))
                if (!contains(name_to_status, deftable_.c_name[k]))
                    allrows.push_back(k);

        return allrows;
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////  Augmented : sample lists  //////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////

    void Augmented::set_sample_list_(int status)
    {
        StatusBlock_ &block = block_at_(status);
        for (size_t k : block.order)
            set_sample_list_for_single_defininition_(status, k);

        // NOTE: by construction of draws list, I am forced to drop private at the end
        list_export = khold(str::drop_private(khold(list_export)));
        block.draws = khold(str::drop_private(khold(block.draws)));
    }

    void Augmented::set_sample_list_for_single_defininition_(int status, size_t k)
    {
        const std::string &astr = deftable_.c_name[k];
        size_t a = khold(astr);

        bool is_to_export = true;
        bool is_to_draw = true;

        if (str::tobool(deftable_.c_const[k]))
        {
            is_to_draw = false;
        }
        else if (!(na::is(deftable_.c_func[k]) && na::is(deftable_.c_value[k])))
        {
            auto [is_to_export_, is_to_draw_] = set_sample_list_for_function_(status, k);
            is_to_export = is_to_export_;
            is_to_draw = is_to_draw_;
        }
        else if (na::is(deftable_.c_value[k]) && na::is(deftable_.c_prior[k]))
        {
            msg::stop(
                "Impossible to initialize augmented " + msg::val(astr) + " without a value or prior");
        }
        else if (na::is(deftable_.c_proposal[k]))
        {
            msg::stop(
                "No proposal was provided for non-constant definition of augmented " +
                msg::val(astr));
        }
        else
        {
            is_any_sample_ = true;
            block_at_(status).sample.emplace_back(a);

            // adaptive if proposal is parametric and scale is NA
            // NOTE: user defined proposal_scale freezes the scale
            if (deftable_.c_proposal[k] != "sample" &&
                na::is(deftable_.c_scale[k]))
            {
                is_any_adaptive_ = true;
                block_at_(status).adaptive.emplace_back(a);
            }
        }

        // NOTE: observation exports are listed separately
        if (deftable_.c_kind[k] == Kind::Obs)
            is_to_export = false;

        // nothing to do if not to export
        if (!is_to_export)
            return;

        if (!contains(list_export, a))
            list_export.emplace_back(a);

        if (is_to_draw)
            block_at_(status).draws.emplace_back(a);
    }

    std::pair<bool, bool> Augmented::set_sample_list_for_function_(int status, size_t k)
    {
        bool is_to_export = true;
        bool is_to_draw = true;

        // NOTE: value definition has priority over function
        const std::string &str_val = deftable_.c_value[k];
        const parse::Definition parser(na::is(str_val) ? deftable_.c_func[k] : str_val);

        mcmc_->data_wide().set_eq_class_data_names(parser.list_data_wide_depen());

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
            if (contains(block_at_(status).draws, khold(bstr)))
                all_consts = false;

        // to draw if not all constant dependencies
        is_to_draw = !all_consts;

        return {is_to_export, is_to_draw};
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////  Augmented : initialize augmented  //////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////

    void Augmented::init_for_evaluation()
    {
        if (khold_.empty())
            return;

        // raise if deterministic observation node for evaluation
        for (size_t k = 0; k < deftable_.c_name.size(); ++k)
            if (deftable_.c_kind[k] == Kind::Obs && na::is(deftable_.c_prior[k]))
                msg::stop(
                    "Deterministic observation node " + msg::val(deftable_.c_name[k]) +
                    " cannot be used for evaluation, only simulation.");

        helper_init_prior_list_();

        // define augmented variables for all augmented individuals in ordefr
        for (Household &house : mcmc_->data_wide())
            for (size_t i : house.infected_hhindexes())
            {
                Individual &member = house.member(i);
                for (size_t k : helper_materialize_definition_order_(member))
                    helper_init_attributes_(member, k);
            }

        // parameters that may move the augmented log-prior (params now registered)
        set_prior_param_depen_();
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////  Augmented : init augmented priors  /////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////

    void Augmented::helper_init_prior_list_()
    {
        for (StatusBlock_ &block : status_blocks_)
        {
            block.depen.resize(khold_.size());

            // associate blocks for ordered prior lists
            const str_table &ordered_priors = get_ordered_priors_(block);
            if (!ordered_priors.empty())
            {
                is_any_prior_ = true;
                block.priors = khold(read::extract_keys(ordered_priors));

                for (const auto &[astr, others] : ordered_priors)
                    block.depen[khold(astr)] = khold(others);
            }
        }
    }

    str_table Augmented::get_ordered_priors_(const StatusBlock_ &block) const
    {
        str_table tree;                         // populate augmented dependency tree
        str_table ordered_priors;               // store hierarchical order of priors
        std::vector<std::string> without_prior; // list of augmented without priors
        for (size_t k : block.order)
        {
            const std::string &aumgstr = deftable_.c_name[k];
            std::vector<std::string> xdepen = {};
            for (const std::vector<std::string> *col : deftable_.source_cols())
                xdepen += parse::Definition((*col)[k]).list_augmented_depen();

            // NOTE: observations hold distributions in c_prior but are NOT priors
            if (na::is(deftable_.c_prior[k]) || deftable_.c_kind[k] == Kind::Obs)
                without_prior.push_back(aumgstr);
            else
                tree[aumgstr] = str::unique(xdepen);
        }

        // prune dependency tree from augmented without prior
        for (auto &[key, vec] : tree)
            vec = vec - without_prior;

        // assign the ordered consumers of each prior for augmented prior computation
        for (const auto &[key, vec] : tree)
            ordered_priors.emplace(key, dependency::order_conditioned_on(key, tree));

        return ordered_priors;
    }

    void Augmented::helper_init_observation_list_()
    {
        for (StatusBlock_ &block : status_blocks_)
        {
            block.obs_depen.resize(khold_.size());

            // direct @A dependency of every block definition (name -> its deps).
            // An observation reaches its sampled ancestors through deterministic
            // nodes, so the inverse map must follow the transitive closure (as the
            // augmented-prior map does), else moving such an ancestor would skip
            // the observation ratio.
            str_table tree;
            for (size_t k : block.order)
            {
                std::vector<std::string> xdepen = {};
                for (const std::vector<std::string> *col : deftable_.source_cols())
                    xdepen += parse::Definition((*col)[k]).list_augmented_depen();
                tree[deftable_.c_name[k]] = str::unique(xdepen);
            }

            for (size_t k : block.order)
            {
                if (deftable_.c_kind[k] != Kind::Obs)
                    continue;

                const size_t obs_key = khold(deftable_.c_name[k]);

                // unique observation keys to return as simulated columns
                if (!contains(list_obs_export, obs_key))
                    list_obs_export.push_back(obs_key);

                // a deterministic observation node carries no distribution
                if (na::is(deftable_.c_prior[k]))
                    continue;

                is_any_obs_ = true;
                block.obs.push_back(obs_key);

                // register the observation under every augmented it transitively
                // reads, so moving any (sampled) ancestor recomputes it
                std::vector<std::string> visited, recursion_stack, depen;
                dependency::walk_tree(deftable_.c_name[k], tree, visited, recursion_stack, depen);
                for (const std::string &lstr : depen)
                    if (lstr != deftable_.c_name[k])
                        block.obs_depen[khold(lstr)].push_back(obs_key);
            }
        }
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////  Augmented : household observations  /////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////

    void Augmented::helper_init_house_obs_list_()
    {
        for (size_t k = 0; k < deftable_.c_name.size(); ++k)
            if (deftable_.c_kind[k] == Kind::House)
                list_house_export.push_back(k);
    }

    void Augmented::init_house_for_simulation()
    {
        DataWide &data_wide = mcmc_->data_wide();

        // reject an unsatisfiable cycle (augmented variable reading a live house observation)
        assert_no_live_house_obs_cycle_();
        // reject a group/contact rate keyed on a generated observation column
        assert_no_obs_rate_collision_();

        house_varia.clear();
        if (list_house_export.empty())
            return;

        // household observations are first-class Variables: build the reduction for
        // every house once and point each standing-in data column at it. Live
        // reductions (over generated observations) are rebuilt as they are drawn.
        const size_t n_houses = data_wide.report_house.total;
        house_varia.assign(n_houses, std::vector<Variable>(list_house_export.size()));
        for (size_t H = 0; H < n_houses; ++H)
        {
            Household &house = data_wide.at(H);
            for (size_t b = 0; b < list_house_export.size(); ++b)
                helper_init_house_attributes_(house, b);
            wire_house_obs_cells_(house);
        }
    }

    void Augmented::wire_house_obs_cells_(Household &house)
    {
        DataWide &data_wide = mcmc_->data_wide();
        for (size_t b = 0; b < list_house_export.size(); ++b)
        {
            const std::string &name = house_obs_name(b);
            if (!data_wide.khold().count(name))
                continue; // export-only household observation: no standing-in column
            const size_t d = data_wide.khold(name);
            for (const Individual &member : house)
                data_wide.wire_cell(d, member, &house_varia[Index::of(house)][b]);
        }
    }

    void Augmented::reset_obs_data_for_simulation(Household &house)
    {
        // clear every member's observation Variables before re-sampling: a member
        // infected in a previous repetition must not leak a stale value into this
        // one. The live household reductions point at these slots and follow them.
        for (Individual &member : house)
            for (size_t a : list_obs_export)
                varia_of_(member).vec[a].clear();
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////  Augmented : data-bound household observations (simulation)  /////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////

    bool Augmented::house_obs_is_live_(size_t k) const
    {
        // live iff a source references an observation node's data column, i.e. a
        // column generated during simulation rather than fixed input data.
        std::vector<std::string> obs_names;
        for (size_t a : list_obs_export)
            obs_names.push_back(khold(a));
        for (const std::vector<std::string> *col : deftable_.source_cols())
        {
            const std::string &expr = (*col)[k];
            if (na::is(expr))
                continue;
            for (const std::string &d : parse::Definition(expr).list_data_wide_depen())
                if (contains(obs_names, d))
                    return true;
        }
        return false;
    }

    void Augmented::assert_no_live_house_obs_cycle_() const
    {
        const DataWide &data_wide = mcmc_->data_wide();
        std::vector<std::string> live_names;
        for (size_t b = 0; b < list_house_export.size(); ++b)
        {
            const size_t k = list_house_export[b];
            const std::string &name = deftable_.c_name[k];
            if (data_wide.khold().count(name) && house_obs_is_live_(k))
                live_names.push_back(name);
        }
        if (live_names.empty())
            return;

        for (size_t k = 0; k < deftable_.c_name.size(); ++k)
        {
            if (deftable_.c_kind[k] != Kind::Augm)
                continue;
            for (const std::vector<std::string> *col : deftable_.source_cols())
            {
                const std::string &expr = (*col)[k];
                if (na::is(expr))
                    continue;
                for (const std::string &dstr : parse::Definition(expr).list_data_wide_depen())
                    if (contains(live_names, dstr))
                        msg::stop(
                            "Cyclic dependency for simulation: augmented " + msg::val(deftable_.c_name[k]) +
                            " reads data column " + msg::val(dstr) + ", a household observation derived from a "
                                                                       "generated observation that depends on augmented variables.\n"
                                                                       "Define " +
                            msg::val(deftable_.c_name[k]) + " without referencing " + msg::val(dstr) +
                            " (an anchor-free distribution; the likelihood is invariant to the time origin).");
            }
        }
    }

    void Augmented::assert_no_obs_entry_collision_() const
    {
        if (list_obs_export.empty())
            return;

        // entry-point data column of each single observation node
        const std::vector<std::string> obs_names = khold(list_obs_export);

        // augmented dependency tree (name -> @A deps), accumulated across status rows
        str_table tree;
        for (size_t k = 0; k < deftable_.c_name.size(); ++k)
        {
            if (deftable_.c_kind[k] != Kind::Augm &&
                deftable_.c_kind[k] != Kind::Obs)
                continue;
            std::vector<std::string> &deps = tree[deftable_.c_name[k]];
            for (const std::vector<std::string> *col : deftable_.source_cols())
                deps += parse::Definition((*col)[k]).list_augmented_depen();
            deps = str::unique(deps);
        }

        const auto data_obs_refs = [&](const std::string &expr)
        {
            std::vector<std::string> hits;
            if (na::is(expr))
                return hits;
            for (const std::string &dstr : parse::Definition(expr).list_data_wide_depen())
                if (contains(obs_names, dstr))
                    hits.push_back(dstr);
            return hits;
        };

        for (size_t k = 0; k < deftable_.c_name.size(); ++k)
        {
            if (deftable_.c_kind[k] != Kind::Augm && deftable_.c_kind[k] != Kind::Obs)
                continue;
            const std::string &name = deftable_.c_name[k];

            // Rule 1: a sampled quantity's prior (distribution + support) must not read an
            // observation. A prior may only depend on parameters and other latents.
            if (!na::is(deftable_.c_prior[k]))
                for (const std::string *col : {&deftable_.c_prior[k], &deftable_.c_min[k], &deftable_.c_max[k]})
                    for (const std::string &dstr : data_obs_refs(*col))
                        msg::stop(
                            "Definition collision: the prior of " + msg::val(name) + " reads observation " + msg::val(dstr) +
                            ".\nA prior cannot depend on an observed or generated value.");

            // Rule 2: a deterministic definition may read an observation (a downstream
            // readout) unless it closes a loop: the observation transitively depends on it.
            for (const std::string *col : {&deftable_.c_value[k], &deftable_.c_func[k]})
                for (const std::string &dstr : data_obs_refs(*col))
                {
                    if (dstr == name)
                        continue; // observation node reading its own entry column
                    // loop iff the observation 'd' is a (transitive) consumer of 'name'
                    const auto consumers = dependency::order_conditioned_on(name, tree);
                    if (contains(consumers, dstr))
                        msg::stop(
                            "Cyclic definition: " + msg::val(name) + " reads observation " + msg::val(dstr) +
                            ", but " + msg::val(dstr) + " depends on " + msg::val(name) +
                            ".\nBreak the loop by defining one of them independently of the other.");
                }
        }
    }

    void Augmented::assert_no_obs_rate_collision_() const
    {
        if (list_obs_export.empty())
            return;

        // entry-point data column of each single observation node
        std::vector<std::string> obs_names = khold(list_obs_export);

        const auto reject = [&](const std::string &origin, const std::string &d)
        {
            msg::stop(
                "Definition collision: a " + origin + " is keyed on observation " + msg::val(d) +
                ".\nThe " + origin + " feeds the transmission model, which the observation is "
                                     "generated from, so it cannot be assigned coherently during simulation. Key the " +
                origin + " on a column that is not an observation node.");
        };

        for (const std::string &d : mcmc_->group_rate.data_wide_depen())
            if (contains(obs_names, d))
                reject("group rate", d);
        for (const std::string &d : mcmc_->contact_rate.data_wide_depen())
            if (contains(obs_names, d))
                reject("contact rate", d);
    }

    void Augmented::helper_init_house_attributes_(Household &house, size_t j)
    {
        const size_t k = list_house_export[j];
        Variable &hvar = house_varia[Index::of(house)][j];
        build_house_obs_variable_(hvar, house, k);
    }

    void Augmented::build_house_obs_variable_(Variable &hvar, Household &house, size_t k)
    {
        const std::string &astr = deftable_.c_name[k];
        hvar.name(astr, ErrType::HouseholdObservation);

        try
        {
            const std::string &str_val = deftable_.c_value[k];
            parse::Definition parser(na::is(str_val) ? deftable_.c_func[k] : str_val);

            // bare reference (no reduction): bind to the first member, which for a
            // singleton household is the household value itself.
            if (parser.is_terminal())
            {
                hvar.transfer(mcmc_->translator.to_attribute(parser.root, &house.member(0)));
                return;
            }

            // reduction: compose the (variadic) function once
            hvar.init_function(compose_func_registries(
                parser.root.func, "function",
                {"custom", "transcendental", "arithmetic",
                 "continuos", "discrete", "infecprofile"},
                str::tobool(deftable_.c_fcumul[k])));

            // bind operands: a data (@D) source is spread one live operand per member;
            // any other source is bound once at household level. Every member is bound
            // (the function folds na.rm-style at eval, skipping NA), so a generated
            // observation drawn later is tracked with no rebuild. Only an immutable NA
            // *constant* (a data cell with no value) is skipped at build, since the
            // function rejects NA constant args and it can never acquire a value. An
            // undrawn observation is empty (not constant), so it is bound live.
            // The first operand seeds the reduction ('x'), the rest are appended.
            size_t op = 0;
            const auto append_operand = [&](const Attribute &attr)
            {
                const std::string key = op == 0 ? "x" : "x" + std::to_string(op);
                ++op;
                hvar.init_func_attr(key, attr);
            };
            for (const auto &child : parser.root.children)
                if (parse::is_data_wide(child.type))
                    for (const Individual &member : house)
                    {
                        const Attribute &attr = mcmc_->translator.to_attribute(child, &member);
                        if (attr.is_constant() && na::is(attr.value()))
                            continue;
                        append_operand(attr);
                    }
                else
                    append_operand(mcmc_->translator.to_attribute(child, &house.member(0)));

            // every operand was NA: the reduction is empty, yield NA (na.rm convention)
            if (op == 0)
            {
                hvar.init_numeric(na::DOUBLE);
                return;
            }

            // range and reshaping
            hvar.init_func_attr("min", mcmc_->translator.to_attribute(deftable_.c_min[k], &house.member(0)));
            hvar.init_func_attr("max", mcmc_->translator.to_attribute(deftable_.c_max[k], &house.member(0)));
            hvar.init_func_reshaping(str::tobool(deftable_.c_discr[k]));
        }
        catch (const std::exception &e)
        {
            msg::stop(
                std::string(e.what()) +
                "\nFor household observation " + msg::val(astr) + " of " +
                house.err_name() + ".");
        }
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////  Augmented : data-bound household observations  //////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////

    bool Augmented::house_obs_is_pure_data_(size_t k) const
    {
        // every source (value/func/min/max) must reference data or numerics only:
        // augmented, parameter and group-rate references are not constant data.
        for (const std::vector<std::string> *col : deftable_.source_cols())
        {
            const std::string &expr = (*col)[k];
            if (na::is(expr))
                continue;
            parse::Definition def(expr);
            if (def.is_augmented_depen() || def.is_parameter_depen() ||
                entry_has_group_rate_(def.root))
                return false;
        }
        return true;
    }

    bool Augmented::entry_has_group_rate_(const parse::Entry &node)
    {
        if (parse::is_group_rate(node.type))
            return true;
        for (const parse::Entry &child : node.children)
            if (entry_has_group_rate_(child))
                return true;
        return false;
    }

    void Augmented::assert_status_obs_consistent_()
    {
        if (!is_any_status_obs())
            return;

        const size_t k = status_obs_def_;
        const std::string &astr = deftable_.c_name[k];

        // status must be weights must be constant data: a status fixed in the data contributes
        // a constant to the likelihood, dropped during inference.
        // NOTE: outcome k resolves to infected status (INFECTED + k); the matching
        // augmented/observation block is selected at simulation (status-specific if
        // declared, otherwise the generic one), so no block count is enforced here.
        if (!house_obs_is_pure_data_(k))
            msg::stop(
                "Status observation " + msg::val(astr) + " must depend only on data or numeric constants.");
    }

    void Augmented::materialize_data_bound_house_obs(DataWide &data_wide)
    {
        if (list_house_export.empty())
            return;

        // build the household-observation Variables (inference home, shared with the
        // simulation init) and point each standing-in data column's read-surface at
        // its Variable, so @D readers resolve the derived value without a data copy.
        const size_t n_houses = data_wide.report_house.total;
        house_varia.assign(n_houses, std::vector<Variable>(list_house_export.size()));

        for (size_t b = 0; b < list_house_export.size(); ++b)
        {
            const size_t k = list_house_export[b];
            const std::string &astr = deftable_.c_name[k];
            // export-only household observation: its name is not a data column
            if (!data_wide.khold().count(astr))
                continue;

            // a household observation standing in for required data must be a pure
            // function of data: a latent/parameter-dependent one has no fixed value
            // during inference (it would impose an indicator constraint).
            if (!house_obs_is_pure_data_(k))
                msg::stop(
                    "Data-bound household observation " + msg::val(astr) + "must be a pure function of data.\n"
                    "Live observations are only meaningful in simulation.");

            const size_t d = data_wide.khold(astr);
            for (size_t H = 0; H < n_houses; ++H)
            {
                Household &house = data_wide.at(H);
                Variable &hvar = house_varia[H][b];
                build_house_obs_variable_(hvar, house, k);
                const double v = hvar.value();

                if (data_wide.requires_discrete(d) && v != std::floor(v))
                    msg::stop(
                        "Household observation " + msg::val(astr) + " for " + house.err_name() +
                        " was expected to be discrete-valued, as requested by model.");
                if (data_wide.requires_positive(d) && !(v > 0.))
                    msg::stop(
                        "Household observation " + msg::val(astr) + " for " + house.err_name() +
                        " was expected to be greater than zero, as requested by model.");

                // enforce inference-optimization
                hvar.optimize();

                for (const Individual &member : house)
                    data_wide.wire_cell(d, member, &hvar);
            }
        }
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
        for (Variable &avar : varia_of_(member).vec)
            avar.clear();
    }

    void Augmented::init_definitions_for_simulation(Individual &member)
    {
        for (size_t k : helper_materialize_definition_order_(member))
            if (varia_of_(member).vec[khold(deftable_.c_name[k])].is_empty())
                helper_init_attributes_(member, k, true);
    }

    std::vector<size_t> Augmented::helper_materialize_definition_order_(Individual &member)
    {
        int status = member.status();

        // if known status, then load its definition
        const auto it = map_status_to_block_.find(status);
        if (it != map_status_to_block_.end())
        {
            const size_t i = it->second;
            varia_of_(member).idx = i;
            return status_blocks_[i].order;
        }

        // otherwise if it exists fallback to generic definition
        const auto generic_it = map_status_to_block_.find(status::MISSING);
        if (generic_it != map_status_to_block_.end())
        {
            const size_t i = generic_it->second;
            map_status_to_block_.emplace(status, i);
            varia_of_(member).idx = i;
            return status_blocks_[i].order;
        }

        if (member.is_infected_status())
            msg::stop(
                "No augmented definitions were found for " +
                member.err_name() + " with status " + std::to_string(status) + ".");

        msg::stop(
            "Generic augmented definition for non-infectious " +
            member.err_name() + " was missing.");
    }

    void Augmented::helper_init_attributes_(Individual &member, size_t k, bool simulation)
    {
        const std::string &astr = deftable_.c_name[k];
        size_t a = khold(astr);

        try
        {
            Variable &avar = varia_of_(member).vec[a];

            if (str::tobool(deftable_.c_const[k]))
            {
                // placeholder for observation node if generated by simulation
                if (simulation &&
                    (deftable_.c_kind[k] == Kind::Obs || deftable_.c_kind[k] == Kind::Status))
                    avar.init_numeric(0.);
                else
                    avar.init_numeric(
                        mcmc_->translator.to_numeric(deftable_.c_value[k], &member).value());
                avar.enforce_user_const(); // enforce user constant choice
            }
            else
                helper_init_func_attributes_(member, k);

            // initialize prior
            helper_init_prior_attributes_(member, k);

            // assign proposal if variable to be sampled
            if (!avar.is_constant())
                avar.init_proposal(deftable_.c_proposal[k], str::tod(deftable_.c_scale[k]));

            // check if satisfies discrete requirement, if any
            if (user_discrete_set_.count(a) && !avar.is_discrete())
                msg::stop(
                    "Definition of " + avar.err_name() +
                    " did not return a discrete object as requested by prior.");
        }
        catch (const std::exception &e)
        {
            msg::stop(
                std::string(e.what()) +
                "\nFor augmented " + msg::val(astr) + " of " +
                member.err_name() + " (status " + std::to_string(member.status()) + ").");
        }
    }

    void Augmented::helper_init_prior_attributes_(Individual &member, size_t k)
    {
        const std::string &astr = deftable_.c_name[k];
        Variable &avar = varia_of_(member).vec[khold(astr)];
        parse::Definition parser(deftable_.c_prior[k]);

        avar.init_prior(parser.root.func);

        if (!avar.has_prior())
            return;
        for (const auto &child : parser.root.children)
            avar.init_prior_attr(child.key, mcmc_->translator.to_attribute(child, &member));
        avar.init_prior_attr("min", mcmc_->translator.to_attribute(deftable_.c_min[k], &member));
        avar.init_prior_attr("max", mcmc_->translator.to_attribute(deftable_.c_max[k], &member));
        avar.init_prior_reshaping(str::tobool(deftable_.c_discr[k]));
    }

    void Augmented::helper_init_func_attributes_(Individual &member, size_t k)
    {
        const std::string &astr = deftable_.c_name[k];
        Variable &avar = varia_of_(member).vec[khold(astr)];
        std::string str_val = deftable_.c_value[k];
        parse::Definition parser(na::is(str_val) ? deftable_.c_func[k] : str_val);

        // pointer mode
        if (parser.is_terminal())
        {
            avar.transfer(mcmc_->translator.to_attribute(parser.root, &member));
            return;
        }

        // function mode
        bool is_primitive = str::tobool(deftable_.c_fcumul[k]);
        avar.init_function(compose_func_registries(
            parser.root.func, "function",
            {"custom", "transcendental", "arithmetic",
             "continuos", "discrete", "infecprofile"},
            is_primitive));

        // attribute
        for (const auto &child : parser.root.children)
        {
            avar.init_func_attr(child.key, mcmc_->translator.to_attribute(child, &member));
            mcmc_->translator.make_cache_dependent(child);
        }

        // range
        avar.init_func_attr("min", mcmc_->translator.to_attribute(deftable_.c_min[k], &member));
        avar.init_func_attr("max", mcmc_->translator.to_attribute(deftable_.c_max[k], &member));

        // discretization and truncation
        avar.init_func_reshaping(str::tobool(deftable_.c_discr[k]));
    }

    /////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////  Augmented : user reset  ////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////

    void Augmented::sample_obs(Individual &member, rnd_gen &prng)
    {
        // draw each observation into its Variable; @D consumers read it straight from
        // there (the data column, if any, is never overwritten).
        for (size_t a : list_obs(member))
            varia_of_(member).vec[a].draw_value(prng);
    }

    void Augmented::simulate_infected_status(Household &house, Individual &member, rnd_gen &prng)
    {
        int new_status = status::INFECTED;

        if (is_any_status_obs())
        {
            Variable &svar = varia_of_(member).vec[status_obs_key_];
            // lazy materialize status observation on first use
            if (svar.is_empty())
                helper_init_attributes_(member, status_obs_def_, true);
            svar.draw_value(prng);

            // NOTE: I am assuming value draws are integer in [0, max - 1]
            new_status += static_cast<int>(svar.value());
        }

        house.set_new_status_(member, new_status);
    }

    void Augmented::restart(rnd_gen &prng)
    {
        if (!is_any_sample())
            return;

        for (Household &house : mcmc_->data_wide())
            for (size_t i : house.infected_hhindexes())
                restart(house.member(i), prng);
    }

    [[noreturn]] void Augmented::throw_obs_init_error(const Individual &member) const
    {
        const StatusBlock_ &block = block_of_(member);

        for (size_t k : block.order)
        {
            if (deftable_.c_kind[k] != Kind::Obs)
                continue;

            const Variable &ovar = varia_of_(member).vec[khold(deftable_.c_name[k])];
            if (std::isfinite(ovar.log_obs()))
                continue;

            // domain: loc + [min, max]
            std::vector<std::string> augm_dep, param_dep;
            for (const parse::Entry &child : parse::Definition(deftable_.c_prior[k]).root.children)
                if (child.key == "loc")
                {
                    child.append_augmented_depen(augm_dep);
                    child.append_parameter_depen(param_dep);
                }
            for (const std::string &col : {deftable_.c_min[k], deftable_.c_max[k]})
            {
                parse::Definition def(col);
                augm_dep += def.list_augmented_depen();
                param_dep += def.list_parameter_depen();
            }

            if (augm_dep.size() > 0)
                msg::stop(
                    "Exceeded attempts to initialize augmented variables of " +
                    member.err_name() + " due to observations.\n" +
                    "Widen augmented prior of " + str::format(str::unique(augm_dep), 3, false) +
                    ", revise the observation domain of " + msg::val(deftable_.c_name[k]) +
                    ", or check the observed individual data.");

            // non-constant parameter in domain cannot be fixed by per-individual resampling
            for (const std::string &pstr : str::unique(param_dep))
            {
                size_t p = mcmc_->param.khold().from_string(pstr);
                if (!mcmc_->param.at(p).is_constant())
                    msg::stop(
                        "Failed to initialize due to observations of " + member.err_name() + ".\n" +
                        "Try to remove free parameter dependencies in observation domain of " + msg::val(deftable_.c_name[k]) + ".");
            }
        }

        // otherwise observed data are just unreachable
        msg::stop(
            "Failed to initialize due to observations of " + member.err_name() + ".\n" +
            "The observed values lie outside of respective domain.\n" +
            "Check the individual data or revise the domain.");
    }

    void Augmented::optimize()
    {
        if (khold_.empty())
            return;

        for (Household &house : mcmc_->data_wide())
            for (size_t i : house.infected_hhindexes())
            {
                Individual &member = house.member(i);
                for (Variable &avar : varia_of_(member).vec)
                    avar.optimize();
            }
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////  Augmented : log-likelihood  //////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////

    double Augmented::commit_augm_log_prior(Individual &member)
    {
        const std::vector<size_t> &priors = list_priors(member);

        if (priors.empty())
            return 0.;
        double new_ll = 0.0;
        for (size_t a : priors)
            new_ll += varia_of_(member).vec[a].log_prior();
        varia_of_(member).log_prior.set(new_ll);
        return new_ll;
    }

    double Augmented::commit_augm_log_prior_ratio(Individual &member, size_t a)
    {
        const StatusBlock_ &block = block_of_(member);
        if (block.depen[a].empty())
            return 0.;
        double delta_ll = 0;
        // NOTE: vector of dependencies 'depen' includes also 'a'
        for (size_t b : block.depen[a])
            delta_ll += varia_of_(member).vec[b].prior_log_ratio();
        varia_of_(member).log_prior.increment(delta_ll);
        return delta_ll;
    }

    double Augmented::commit_obs_log_lik(Individual &member)
    {
        if (list_obs(member).empty())
            return 0.;
        const double new_ll = eval_log_obs(member);
        varia_of_(member).log_obs.set(new_ll);
        return new_ll;
    }

    double Augmented::commit_obs_log_lik_ratio(Individual &member, size_t a)
    {
        const StatusBlock_ &block = block_of_(member);
        // skip unless the moved augmented 'a' parameterizes some observation node
        if (block.obs_depen[a].empty())
            return 0.;

        // log-proba must be recomputed (observation value is just a datum)
        double new_ll = 0.0;
        for (size_t b : block.obs)
            new_ll += varia_of_(member).vec[b].log_obs();
        const double delta_ll = new_ll - varia_of_(member).log_obs.value();
        varia_of_(member).log_obs.set(new_ll);
        return delta_ll;
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
            msg::warn(
                "Ignored" + msg::arg("id_draws") + ": no augmented variables to sample were defined.");
            return output;
        }

        output.reserve(id_draws.size());
        std::vector<std::string> existing_ids = {};
        for (const Household &house : mcmc_->data_wide())
            for (const Individual &member : house)
            {
                if (!contains(id_draws, member.id()))
                    continue;

                // unique IDs are granted by class DataWide load
                existing_ids.push_back(member.id());

                // check status of chosen individual
                if (member.is_infected_status())
                    output.emplace_back(&member);
                else
                    msg::warn(
                        "Ignored draws of non infected " + member.err_name() + ".");
            }
        output.shrink_to_fit();

        return output;
    }
}
