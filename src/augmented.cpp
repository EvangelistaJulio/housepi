/* housepi - Epidemiological Analysis of Household Outbreaks.
 *
 * Copyright © 2024-2026 Institut Pasteur (Paris), CNRS and INSERM.
 * Distributed under the terms of the GNU General Public License (GPLv3).
 * See the LICENSE file for details.
 */

#include <algorithm>
#include <cmath>
#include <memory>
#include <random>
#include <stdexcept>
#include <string>
#include <unordered_set>
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
#include "registry.h"
#include "variable.h"
#include "translator.h"

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

    /////////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////  *** InfecProfile : struct ***  //////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////

    InfecProfile::InfecProfile() = default;
    InfecProfile::~InfecProfile() = default;

    void InfecProfile::build(Translator &translator, const str_table &input)
    {
        augm_definitions.clear();

        if (input.empty())
            return;

        // check if definition is generic or augmented
        if (read::csv_n_rows(input) > 1 ||
            parse::Definition(input.at("func")[0]).consumes_augmented() ||
            parse::Definition(input.at("min")[0]).consumes_augmented() ||
            parse::Definition(input.at("max")[0]).consumes_augmented() ||
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
                                 {"discr", {input.at("discr")[0]}}});
        lookup_.clear();
        augm_definitions.clear();

        // range minimum
        generic->init_attr("min", translator.to_attribute(input.at("min")[0]));

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
            // a discrete-time model discretizes whatever profile it is handed
            bool discr = is_discr_requested_ || str::tobool(input.at("discr")[r]);

            augm_definitions["name"].push_back(infprof::AUGM_FUNC);
            augm_definitions["status"].push_back(input.at("status")[r]);
            augm_definitions["func"].push_back(input.at("func")[r]);
            augm_definitions["fcumul"].push_back("false");
            augm_definitions["value"].push_back("");
            augm_definitions["const"].push_back("false");
            augm_definitions["discr"].push_back(str::format(discr));
            augm_definitions["model"].push_back("");
            augm_definitions["min"].push_back(input.at("min")[r]);
            augm_definitions["max"].push_back(input.at("max")[r]);
            augm_definitions["proposal"].push_back("");
            augm_definitions["kind"].push_back("augm");
        }
    }

    std::unique_ptr<Distribution> InfecProfile::helper_build_(Translator &translator, const str_table &input)
    {
        parse::Definition parser(parse::Definition(input.at("func")[0]));

        std::unique_ptr<Distribution> output = compose_distr_registries(
            parser.root.func,
            ErrType::InfecProfile,
            {Registry::InfecProfile, Registry::Continuos, Registry::Discrete});

        for (const auto &node : parser.root.children)
        {
            output->init_attr(node.key, translator.to_attribute(node));
            translator.mark_param_produces_transm_cache(node);
        }

        if (is_discr_requested_ ||
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

    /// Registries searched when composing an augmented / household-observation
    /// reduction function (shared by per-member and per-household builders).
    static const std::vector<Registry> AUGM_FUNC_REGISTRIES = {
        Registry::Custom, Registry::Transcendental, Registry::Arithmetic,
        Registry::Continuos, Registry::Discrete, Registry::InfecProfile};

    /// Error-message family of a definition row's level.
    static ErrType kind_to_err_type(Augmented::Kind kind)
    {
        switch (kind)
        {
        case Augmented::Kind::Member:
            return ErrType::Observation;
        case Augmented::Kind::House:
            return ErrType::HouseholdObservation;
        default:
            return ErrType::Augmented;
        }
    }

    /// Map a definition's kind label to its typed Kind.
    static Augmented::Kind string_to_kind(const std::string &s)
    {
        if (s == "member")
            return Augmented::Kind::Member;
        if (s == "house")
            return Augmented::Kind::House;
        if (s == "status")
            return Augmented::Kind::Status;
        // if (s == "augm")
        return Augmented::Kind::Augm;
    }

    void Augmented::DefTable_::clear()
    {
        for (auto col : {&c_name, &c_status, &c_func, &c_fcumul,
                         &c_value, &c_const, &c_discr, &c_model,
                         &c_min, &c_max, &c_proposal})
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
        c_model += input.at("model");
        c_min += input.at("min");
        c_max += input.at("max");
        c_proposal += input.at("proposal");
        for (const std::string &s : input.at("kind"))
            c_kind.push_back(string_to_kind(s));
    }

    void Augmented::build(InfecProfile &infec_profile, const str_table &input)
    {
        build_restore_internal_();            // reset internal state
        build_assert_input_(input);           // check entries
        build_append_profile_(infec_profile); // add infectivity profile

        // skip if no augmented definitions
        if (deftable_.empty())
            return;

        build_definition_by_status_();         // process definition by status
        build_set_allocation_order_();         // allocation and order
        init_surrogate_list_();                // surrogate columns
        init_observation_list_();              // observation nodes
        init_transm_producers_();              // keys reaching the transmission side
        assert_status_obs_consistent_();       // assert status observation
        assert_no_obs_entry_collision_();      // no definition reads an obs entry column
        assert_no_sim_override_latent_read_(); // simulation override reads data/params only
    }

    void Augmented::build_restore_internal_()
    {
        is_any_sample_ = false;
        is_any_augm_model_ = false;
        is_any_obs_model_ = false;
        param_produces_augm_ll_.clear();
        param_produces_obs_ll_.clear();

        deftable_.clear();
        khold_.clear_optional();
        map_status_to_def_row_.clear();
        map_status_to_block_.clear();
        status_blocks_.clear();
        list_augm_export.clear();
        list_obs_member.clear();
        list_obs_house.clear();
        list_surrogate.clear();
        varia.clear();
        house_varia.clear();
        status_obs_def_ = na::INDEX;
        status_obs_key_ = na::INDEX;
    }

    void Augmented::build_assert_input_(const str_table &input)
    {
        if (input.empty())
            return;

        for (const std::string &astr : input.at("name"))
            if (astr == infprof::AUGM_FUNC)
                msg::stop(
                    "Tried to use the reserved augmented name " +
                    msg::val(astr) + ".\nPlease, choose a different one.");

        deftable_.append(input);
    }

    void Augmented::build_append_profile_(InfecProfile &infec_profile)
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

    // UPDATE: these flags are not precise if multi-status!
    bool Augmented::is_member_obs(const std::string &name) const noexcept
    {
        for (size_t k = 0; k < deftable_.c_name.size(); ++k)
            if (deftable_.c_kind[k] == Kind::Member && deftable_.c_name[k] == name)
                return true;
        return false;
    }

    bool Augmented::is_house_obs(const std::string &name) const noexcept
    {
        for (size_t k = 0; k < deftable_.c_name.size(); ++k)
            if (deftable_.c_kind[k] == Kind::House && deftable_.c_name[k] == name)
                return true;
        return false;
    }

    void Augmented::build_definition_by_status_()
    {
        for (size_t k = 0; k < deftable_.c_name.size(); ++k)
        {
            // skip household-level observations, not per-member
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

    void Augmented::build_set_allocation_order_()
    {
        std::unordered_set<int> all_infected_statuses = {};

        // allocate storage
        varia.resize(mcmc_->data_wide().report_indiv.total);

        size_t n = khold_.size();
        for (const Household &house : mcmc_->data_wide())
            for (const Individual &indiv : house)
            {
                IndivBlock &iaugm = varia_of_(indiv);
                iaugm.vec.resize(n);
                for (size_t a = 0; a < n; ++a)
                    if (is_member_obs(khold(a)))
                        iaugm.vec[a].name(khold(a), ErrType::Observation);
                    else
                        iaugm.vec[a].name(khold(a), ErrType::Augmented);

                if (indiv.is_infected_status())
                    all_infected_statuses.emplace(indiv.status());
            }

        // skip if no per-member definition to place into status blocks
        // NOTE: household-level observations need no status block
        if (map_status_to_def_row_.empty())
            return;

        // first, load status specific user definitions
        // NOTE: NA (i. e. generic) does not necessarily need checking
        for (auto &[status, vec] : map_status_to_def_row_)
            if (status::is_infected(status))
                set_definition_order_for_(status, status);

        // associate observed status to corresponding status block
        // if not defined specifically, then associate to generic
        for (int status : all_infected_statuses)
            if (!map_status_to_block_.count(status))
            {
                if (!map_status_to_block_.count(status::MISSING))
                    set_definition_order_for_(status::MISSING, status);
                map_status_to_block_.emplace(status, map_status_to_block_.at(status::MISSING));
            }

        // if no observed infected status or status-specific definitions
        // then load generic definitition
        if (all_infected_statuses.empty() ||
            (map_status_to_def_row_.size() == 1 &&
             status::is_missing(map_status_to_def_row_.begin()->first)))
            set_definition_order_for_(status::MISSING, status::MISSING);
    }

    /// Report the observed status left without definitions.
    static std::string undeclared_status_message_(int err_status)
    {
        if (status::is_missing(err_status))
            return "No generic augmented definition was declared.\n"
                   "Provide a definition without " +
                   msg::arg("status") + " to apply it to every infected status.";

        return "No augmented definition was declared for status " + std::to_string(err_status) + ".\n" +
               "Try to associate a status-specific definition.\n" +
               "In alternative, drop " + msg::arg("status") + " from existing definitions to make them generic.";
    }

    void Augmented::set_definition_order_for_(int status, int err_status)
    {
        // skip if already passed
        if (map_status_to_block_.count(status))
            return;
        // throw if not defined
        if (!map_status_to_def_row_.count(status))
            msg::stop(undeclared_status_message_(err_status));

        const std::vector<size_t> rows = get_unsorted_definition_rows_(status);
        dpdcy::prod_graph pgraph = build_producer_graph_(rows);
        dict<size_t> name_to_row; // reverse dictionary name-to-row-index
        for (size_t k : rows)
            name_to_row[deftable_.c_name[k]] = k;

        assert_producers_declared_(pgraph, err_status);

        // assert acyclicity and translate order from name to deftable row index
        try
        {
            // NOTE: this method is just called once per status
            map_status_to_block_.emplace(status, status_blocks_.size());
            status_blocks_.emplace_back(status);

            // ignore spurious hidden definition if no kept definition consumes it
            // NOTE: infectivity profile (hidden) can consume hidden
            auto is_kept_consumer = [](const std::string &c)
            { return !str::is_hidden(c) || c == infprof::AUGM_FUNC; };

            for (const std::string &astr : dpdcy::order_producers_first(pgraph, ErrType::Augmented))
            {
                const auto consumers = dpdcy::transitive_consumers(astr, pgraph, ErrType::Augmented);
                if (astr == infprof::AUGM_FUNC || !str::is_hidden(astr) ||
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
                sms += "\nCheck status-specific definitions for status " + std::to_string(err_status) + ".";
            msg::stop(sms);
        }
    }

    dpdcy::prod_graph Augmented::build_producer_graph_(const std::vector<size_t> &rows) const
    {
        dpdcy::prod_graph pgraph;
        for (size_t k : rows)
        {
            auto &producers = pgraph[deftable_.c_name[k]];
            for (const DefTable_::col *col : deftable_.source_cols())
                for (const std::string &astr : parse::Definition((*col)[k]).list_augmented_producers())
                    producers.insert(astr);
        }
        return pgraph;
    }

    void Augmented::assert_producers_declared_(const dpdcy::prod_graph &pgraph,
                                               int status) const
    {
        for (const auto &node : pgraph)
            for (const std::string &astr : node.second)
            {
                if (pgraph.count(astr))
                    continue;

                std::string sms;

                // declared for another status and report accordingly
                if (contains(deftable_.c_name, astr))
                {
                    if (status::is_missing(status))
                        sms = "No generic augmented definition was found for " + msg::val(astr) + ".";
                    else
                        sms = "No augmented definition was found for " + msg::val(astr) +
                              " for status " + std::to_string(status) + ".";
                    msg::stop(sms);
                }

                // declared elsewhere and forgot sourcing
                sms = "Could not find an augmented definition for " + msg::val(astr) + ".";
                if (mcmc_->data_wide().khold().count(astr))
                    sms += "\nDid you mean " + msg::code(astr + "@D") + "?";
                else if (mcmc_->group_risk.khold().count(astr))
                    sms += "\nDid you mean " + msg::code(astr + "@G") + "?";
                else if (mcmc_->param.khold().count(astr))
                    sms += "\nDid you mean " + msg::code(astr + "@P") + "?";
                else // true augmented definition missing
                    sms += "\nDeclare it with " + msg::fn("add_augmented") + ".";
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
            set_sample_list_for_single_definition_(status, k);

        // NOTE: to construct draws list, I drop hidden at the end
        list_augm_export = khold(str::drop_hidden(khold(list_augm_export)));
        block.draws = khold(str::drop_hidden(khold(block.draws)));
    }

    void Augmented::set_sample_list_for_single_definition_(int status, size_t k)
    {
        const std::string &astr = deftable_.c_name[k];
        size_t a = khold(astr);

        bool is_to_export = true;
        bool is_to_draw = true;

        if (str::tobool(deftable_.c_const[k]))
        {
            is_to_draw = false;
        }
        // deterministic definition: a function, or a value with no model. A value
        // paired with a model is a sampled variable (the value is a simulation-only
        // override, see helper_init_attributes_), so it falls through to the sample branch.
        else if (!na::is(deftable_.c_func[k]) ||
                 (!na::is(deftable_.c_value[k]) && na::is(deftable_.c_model[k])))
        {
            auto [is_to_export_, is_to_draw_] = set_sample_list_for_function_(status, k);
            is_to_export = is_to_export_;
            is_to_draw = is_to_draw_;
        }
        else if (na::is(deftable_.c_value[k]) && na::is(deftable_.c_model[k]))
        {
            msg::stop(
                "Impossible to initialize augmented " + msg::val(astr) + " without a value or model");
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

            // adaptive only for a random-walk kernel carrying no pinned scale
            // NOTE: a user-defined scale (folded into the proposal call) freezes it
            parse::Definition pdef(deftable_.c_proposal[k]);
            const std::string &pf = pdef.root.func;
            bool has_scale = false;
            for (const parse::Entry &child : pdef.root.children)
                if (child.key == "scale")
                    has_scale = true;

            bool is_random_walk = (pf == "rw" || pf == "rw_log" || pf == "rw_unif");
            if (is_random_walk && !has_scale)
            {
                is_any_adaptive_ = true;
                block_at_(status).adaptive.emplace_back(a);
            }
        }

        // NOTE: observation exports are listed separately
        if (deftable_.c_kind[k] == Kind::Member)
            is_to_export = false;

        // nothing to do if not to export
        if (!is_to_export)
            return;

        if (!contains(list_augm_export, a))
            list_augm_export.emplace_back(a);

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

        // raw data are equivalence-class covariates if not-produced
        std::vector<std::string> covariates;
        for (const std::string &dstr : parser.list_data_wide_producers())
            if (!is_member_obs(dstr) && !is_house_obs(dstr))
                covariates.push_back(dstr);
        mcmc_->data_wide().set_eq_class_data_names(covariates);

        // NOTE: do not export if not callable (otherwise crashes)
        if (!parser.is_callable_implicitly("x"))
            is_to_export = false;

        bool all_consts = true; // start assuming all constants
        for (const std::string &pstr : parser.list_parameter_producers())
            all_consts &= mcmc_->param.at(mcmc_->param.khold(pstr)).is_constant();

        // inspect status specific augmented definitions and check if to draw (i. e. not const)
        // NOTE: taking advatange of the order!
        for (const std::string &bstr : parser.list_augmented_producers())
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

        materialize_surrogate_columns(false);
        init_model_list_();

        // define augmented variables for all augmented individuals in order
        for (const Household &house : mcmc_->data_wide())
            for (const Individual &indiv : house.infected_indivs())
                for (size_t k : helper_materialize_definition_order_(indiv))
                    helper_init_attributes_(indiv, k);

        // parameters that may move augmented/observation models
        set_param_producer_lists_();
    }

    void Augmented::set_param_producer_lists_()
    {
        // mark parameter as producers ofaugmented/observation models
        param_produces_augm_ll_.assign(mcmc_->param.khold().size(), false);
        param_produces_obs_ll_.assign(mcmc_->param.khold().size(), false);

        std::vector<bool> augm_produces_obs(khold_.size(), false);
        for (const StatusBlock_ &block : status_blocks_)
            for (size_t a = 0; a < block.obs_consumers.size(); ++a)
                if (!block.obs_consumers[a].empty())
                    augm_produces_obs[a] = true;

        for (size_t k = 0; k < deftable_.c_name.size(); ++k)
        {
            if (deftable_.c_name[k] == infprof::AUGM_FUNC ||
                deftable_.c_kind[k] == Kind::House ||
                deftable_.c_kind[k] == Kind::Status)
                continue;

            const bool is_obs = deftable_.c_kind[k] == Kind::Member; // else Kind::Augm
            // a deterministic augmented forwards its parameters to the observation it produces
            const bool to_obs = is_obs ||
                                (na::is(deftable_.c_model[k]) &&
                                 augm_produces_obs[khold(deftable_.c_name[k])]);

            for (const DefTable_::col *col : deftable_.source_cols())
                for (size_t p : mcmc_->param.khold(parse::Definition((*col)[k]).list_parameter_producers()))
                {
                    if (!is_obs)
                        param_produces_augm_ll_[p] = true;
                    if (to_obs)
                        param_produces_obs_ll_[p] = true;
                }
        }
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////  Augmented : init augmented models  /////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////

    void Augmented::init_model_list_()
    {
        for (StatusBlock_ &block : status_blocks_)
        {
            block.consumers.resize(khold_.size());

            // index the consumer lists by their producer key
            const str_table &model_consumers = get_model_consumers_(block);
            if (!model_consumers.empty())
            {
                is_any_augm_model_ = true;
                block.augm_models = khold(read::extract_keys(model_consumers));

                for (const auto &[astr, consumers] : model_consumers)
                    block.consumers[khold(astr)] = khold(consumers);
            }
        }
    }

    str_table Augmented::get_model_consumers_(const StatusBlock_ &block) const
    {
        const dpdcy::prod_graph pgraph = build_producer_graph_(block.order);

        std::vector<std::string> without_model;
        for (size_t k : block.order)
            if (na::is(deftable_.c_model[k]) || deftable_.c_kind[k] == Kind::Member)
                without_model.push_back(deftable_.c_name[k]);

        // walk the whole graph, so a chain through a modelless node keeps its edge,
        // then keep only the models whose density a move can shift
        // NOTE: an unshifted consumer cancels in the ratio, a missing one does not
        str_table model_consumers;
        for (const auto &[key, producers] : pgraph)
            if (!contains(without_model, key))
                model_consumers.emplace(key,
                                        dpdcy::transitive_consumers(key, pgraph, ErrType::Augmented) -
                                            without_model);

        return model_consumers;
    }

    void Augmented::init_observation_list_()
    {
        for (StatusBlock_ &block : status_blocks_)
        {
            block.obs_consumers.resize(khold_.size());

            // An observation reaches its sampled producers through deterministic
            // nodes, so the consumer map must follow the transitive closure (as the
            // augmented-model map does), else moving such a producer would skip
            // the observation ratio.
            dpdcy::prod_graph pgraph = build_producer_graph_(block.order);

            for (size_t k : block.order)
            {
                if (deftable_.c_kind[k] != Kind::Member)
                    continue;

                // a surrogate column is data, not drawn
                bool skip = false;
                for (const Surrogate &surr : list_surrogate)
                    if (surr.row == k)
                        skip = surr.is_static;
                if (skip)
                    continue;

                const size_t obs_key = khold(deftable_.c_name[k]);

                // unique observation keys to return as simulated columns
                if (!contains(list_obs_member, obs_key))
                    list_obs_member.push_back(obs_key);

                // a deterministic observation node carries no distribution
                if (na::is(deftable_.c_model[k]))
                    continue;

                is_any_obs_model_ = true;
                block.obs_models.push_back(obs_key);

                // list the observation under each of its producers,
                // so moving any  of them recomputes it
                for (const std::string &astr : dpdcy::transitive_producers(deftable_.c_name[k], pgraph, ErrType::Augmented))
                    if (astr != deftable_.c_name[k]) // transitive_producers returns 'start' last
                        block.obs_consumers[khold(astr)].push_back(obs_key);
            }
        }
    }

    void Augmented::init_transm_producers_()
    {
        // transmission reads only model-declared keys
        // and eventually augmented infectivity profile
        std::vector<std::string> transm_astr;
        for (size_t a : user_transm_set_)
            transm_astr.push_back(khold(a));
        transm_astr.push_back(infprof::AUGM_FUNC);

        for (StatusBlock_ &block : status_blocks_)
        {
            block.produces_transm.assign(khold_.size(), false);
            // declaring no transmission-key skip computation
            if (user_transm_set_.empty())
                continue;

            // whatever a node consumes, transmission consumes
            const dpdcy::prod_graph pgraph = build_producer_graph_(block.order);
            for (const std::string &node : transm_astr)
                if (pgraph.count(node)) // a root needs no definition in every block
                    for (const std::string &pstr :
                         dpdcy::transitive_producers(node, pgraph, ErrType::Augmented))
                        block.produces_transm[khold(pstr)] = true;
        }
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////  Augmented : surrogate columns  //////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////

    void Augmented::init_surrogate_list_()
    {
        // look for surrogate column definitions
        // i.e. deterministic definition producing standing for a data column
        dict<size_t> name_to_row;
        for (size_t k = 0; k < deftable_.c_name.size(); ++k)
            // household reduction or member-level definition carrying no model
            if (deftable_.c_kind[k] == Kind::House ||
                (deftable_.c_kind[k] == Kind::Member && na::is(deftable_.c_model[k])))
                name_to_row[deftable_.c_name[k]] = k;

        if (name_to_row.empty())
            return;

        // order producer-before-consumer surrogate columns
        dpdcy::prod_graph pgraph;
        for (const auto &[nm, k] : name_to_row)
        {
            auto &producers = pgraph[nm];
            for (const DefTable_::col *col : deftable_.source_cols())
                for (const std::string &dstr : parse::Definition((*col)[k]).list_data_wide_producers())
                    if (dstr != nm && name_to_row.count(dstr))
                        producers.insert(dstr);
        }

        // extract member-level live observations (with a model)
        std::unordered_set<std::string> live_names;
        for (size_t k = 0; k < deftable_.c_name.size(); ++k)
            if (deftable_.c_kind[k] == Kind::Member && !na::is(deftable_.c_model[k]))
                live_names.insert(deftable_.c_name[k]);

        // build list of surrogate data columns in order
        for (const std::string &nm : dpdcy::order_producers_first(pgraph, ErrType::Observation))
        {
            Surrogate surr;
            surr.row = name_to_row.at(nm);
            surr.is_house = deftable_.c_kind[surr.row] == Kind::House;
            surr.is_static = def_is_pure_data_(surr.row);

            // update house-level observations
            if (surr.is_house)
            {
                surr.slot = list_obs_house.size();
                list_obs_house.push_back(surr.row);
            }

            // check if depends on live, then it's live
            if (surr.is_static)
                for (const DefTable_::col *col : deftable_.source_cols())
                    for (const std::string &dstr : parse::Definition((*col)[surr.row]).list_data_wide_producers())
                        if (live_names.count(dstr))
                            surr.is_static = false;
            if (!surr.is_static)
                live_names.insert(nm);

            list_surrogate.push_back(std::move(surr));
        }
    }

    std::vector<std::string> Augmented::get_produced_cols()
    {
        std::vector<std::string> produced_cols = khold(list_obs_member);
        for (const Surrogate &surr : list_surrogate)
            produced_cols.push_back(deftable_.c_name[surr.row]);
        return produced_cols;
    }

    void Augmented::init_house_for_simulation()
    {
        DataWide &data_wide = mcmc_->data_wide();

        assert_no_live_house_obs_cycle_();
        data_wide.reset_readcells(get_produced_cols());

        // reserve household-level observations
        house_varia.clear();
        if (!list_obs_house.empty())
            house_varia.assign(data_wide.report_house.total,
                               std::vector<Variable>(list_obs_house.size()));

        for (const Household &house : data_wide)
        {
            // wire stochastic member-level observations to data read-cell
            for (const Individual &indiv : house)
                for (size_t a : list_obs_member)
                    data_wide.wire_readcell(data_wide.khold(khold(a)), indiv,
                                            &varia_of_(indiv).vec[a]);

            // build and wire surrogate columns in dependency order
            for (const Surrogate &surr : list_surrogate)
                // defer live member observation (drawn with the member)
                if (surr.is_house || surr.is_static)
                    build_surrogate_column_(surr, house, false);
        }
    }

    void Augmented::build_surrogate_column_(const Surrogate &surr,
                                            const Household &house,
                                            bool for_evaluation)
    {
        DataWide &data_wide = mcmc_->data_wide();
        const std::string &name = deftable_.c_name[surr.row];
        const size_t d = data_wide.khold(name);

        // household-level
        if (surr.is_house)
        {
            Variable &hvar = house_varia[Index::of(house)][surr.slot];
            build_definition_variable_(hvar, surr.row, house.at(0), &house);
            if (for_evaluation)
            {
                assert_surrogate_datum_requirements_(d, hvar.value(), name, house);
                hvar.finalize(); // enforce inference-finalization
            }
            for (const Individual &indiv : house)
                data_wide.wire_readcell(d, indiv, &hvar);
            return;
        }

        // member level
        const size_t a = khold(name);
        for (const Individual &indiv : house)
        {
            Variable &avar = varia_of_(indiv).vec[a];
            build_definition_variable_(avar, surr.row, indiv, nullptr);
            if (for_evaluation)
            {
                assert_surrogate_datum_requirements_(d, avar.value(), name, house);
                avar.finalize(); // enforce inference-finalization
            }
            data_wide.wire_readcell(d, indiv, &avar);
        }
    }

    void Augmented::assert_surrogate_datum_requirements_(
        size_t d, double v, const std::string &name, const Household &house) const
    {
        const DataWide &data_wide = mcmc_->data_wide();
        if (data_wide.requires_discrete(d) && v != std::floor(v))
            msg::stop(
                "Observation " + msg::val(name) + " for " + data_wide.err_name(house) +
                " was expected to be discrete-valued, as requested by model.");
        if (data_wide.requires_positive(d) && !(v > 0.))
            msg::stop(
                "Observation " + msg::val(name) + " for " + data_wide.err_name(house) +
                " was expected to be greater than zero, as requested by model.");
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////  Augmented : definition collisions and variable build  ////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////

    /// Data names (@D) referenced by 'expr' that are also in 'names'.
    std::vector<std::string> referred_data(std::string_view expr,
                                           const std::vector<std::string> &names)
    {
        parse::Definition def(expr);
        std::vector<std::string> output = {};
        if (def.is_na())
            return output;
        for (const std::string &dstr : def.list_data_wide_producers())
            if (contains(names, dstr))
                output.push_back(dstr);
        return output;
    }

    void Augmented::assert_no_live_house_obs_cycle_() const
    {
        std::vector<std::string> live_names;
        for (const Surrogate &surr : list_surrogate)
            if (surr.is_house && !surr.is_static)
                live_names.push_back(deftable_.c_name[surr.row]);
        if (live_names.empty())
            return;

        for (size_t k = 0; k < deftable_.c_name.size(); ++k)
        {
            if (deftable_.c_kind[k] != Kind::Augm)
                continue;
            for (const DefTable_::col *col : deftable_.source_cols())
                for (const std::string &dstr : referred_data((*col)[k], live_names))
                    msg::stop(
                        "Found cyclic dependency in simulation.\nAugmented " +
                        msg::val(deftable_.c_name[k]) + " reads generated household observation " + msg::val(dstr) +
                        ".\nDefine " + msg::val(deftable_.c_name[k]) +
                        " without referencing " + msg::val(dstr) + ", or viceversa.");
        }
    }

    void Augmented::assert_no_obs_entry_collision_() const
    {
        if (list_obs_member.empty())
            return;

        // entry-point data column of each single observation node
        const std::vector<std::string> obs_names = khold(list_obs_member);

        // producer graph (name -> the @A names it reads), over augmented + obs rows
        std::vector<size_t> rows;
        for (size_t k = 0; k < deftable_.c_name.size(); ++k)
            if (deftable_.c_kind[k] == Kind::Augm || deftable_.c_kind[k] == Kind::Member)
                rows.push_back(k);
        const dpdcy::prod_graph pgraph = build_producer_graph_(rows);

        for (size_t k = 0; k < deftable_.c_name.size(); ++k)
        {
            if (deftable_.c_kind[k] != Kind::Augm && deftable_.c_kind[k] != Kind::Member)
                continue;
            const std::string &name = deftable_.c_name[k];

            // The two rules below partition DefTable_::source_cols() by role: the
            // model support {model,min,max} vs the deterministic value {value,func}.
            // Keep them in sync with source_cols() if that set ever changes.

            // Rule 1: a sampled quantity's model (distribution + support) must not read an
            // observation. Its producers may only be parameters and other latents.
            if (!na::is(deftable_.c_model[k]))
                for (const std::string *col : {&deftable_.c_model[k], &deftable_.c_min[k], &deftable_.c_max[k]})
                    for (const std::string &dstr : referred_data(*col, obs_names))
                        msg::stop(
                            "Definition collision: the model of " + msg::val(name) + " reads observation " + msg::val(dstr) +
                            ".\nAn augmented model cannot depend on an observed or generated value.");

            // Rule 2: a deterministic definition may read an observation (a downstream
            // readout) unless it closes a loop, i.e. the observation is also its consumer.
            for (const std::string *col : {&deftable_.c_value[k], &deftable_.c_func[k]})
                for (const std::string &dstr : referred_data(*col, obs_names))
                {
                    if (dstr == name)
                        continue; // observation node reading its own entry column
                    // loop iff the observation 'd' is a (transitive) consumer of 'name'
                    const auto consumers = dpdcy::transitive_consumers(name, pgraph, ErrType::Augmented);
                    if (contains(consumers, dstr))
                        msg::stop(
                            "Cyclic definition: " + msg::val(name) + " reads observation " + msg::val(dstr) +
                            ", but " + msg::val(dstr) + " depends on " + msg::val(name) +
                            ".\nBreak the loop by defining one of them independently of the other.");
                }
        }
    }

    void Augmented::assert_no_sim_override_latent_read_() const
    {
        // a simulation-only override (a value paired with a model) is applied when the
        // variable is first drawn in simulation, before any other latent, so it may depend
        // on data columns and parameters only -- never on augmented variables or observations
        // (which are not yet available at that point). Unlike the obs-entry collision check,
        // this must run even when no observation is defined.
        const std::vector<std::string> obs_names = khold(list_obs_member);
        for (size_t k = 0; k < deftable_.c_name.size(); ++k)
        {
            if (deftable_.c_kind[k] != Kind::Augm ||
                na::is(deftable_.c_model[k]) || na::is(deftable_.c_value[k]))
                continue;
            const std::string &name = deftable_.c_name[k];
            if (parse::Definition(deftable_.c_value[k]).consumes_augmented())
                msg::stop(
                    "Simulation override of " + msg::val(name) +
                    " reads an augmented variable.\nA value paired with a model is applied "
                    "before other latents are drawn; depend on data or parameters only.");
            for (const std::string &dstr : referred_data(deftable_.c_value[k], obs_names))
                msg::stop(
                    "Simulation override of " + msg::val(name) + " reads observation " + msg::val(dstr) +
                    ".\nA value paired with a model is applied during simulation, before " +
                    msg::val(dstr) + " is generated; depend on data or parameters only.");
        }
    }

    void Augmented::build_definition_variable_(Variable &var,
                                               size_t k,
                                               const Individual &anchor,
                                               const Household *house)
    {
        const Kind kind = deftable_.c_kind[k];
        const std::string &astr = deftable_.c_name[k];
        var.name(astr, kind_to_err_type(kind));

        try
        {
            // a deterministic definition is driven by its value, anything else by its
            // function: a value paired with a model is a simulation-only override,
            // applied before this
            const bool has_model = !na::is(deftable_.c_model[k]);
            const std::string &value = deftable_.c_value[k];
            parse::Definition parser(has_model || na::is(value) ? deftable_.c_func[k] : value);

            // pointer mode: a bare reference points at what it names
            if (parser.is_terminal())
            {
                var.transfer(mcmc_->translator.to_attribute(parser.root, &anchor));
                return;
            }

            // function mode: compose it once
            var.init_function(compose_func_registries(
                parser.root.func, ErrType::Function, AUGM_FUNC_REGISTRIES,
                str::tobool(deftable_.c_fcumul[k])));

            // a surrogate member-observation materializes a data column for every member:
            // where the cells it reads are empty it has nothing to stand for
            const bool is_surrogate = kind == Kind::Member && !has_model;

            bool is_defined = true;
            if (kind == Kind::House)
                is_defined = bind_house_reduction_(var, parser.root, *house, anchor);
            else
                for (const parse::Entry &child : parser.root.children)
                {
                    const Attribute &attr = mcmc_->translator.to_attribute(child, &anchor);
                    if (is_surrogate && attr.is_constant() && na::is(attr.value()))
                    {
                        is_defined = false;
                        break;
                    }
                    var.init_func_attr(child.key, attr);
                    mcmc_->translator.mark_param_produces_transm_cache(child);
                }

            // no operand carries a value: the definition yields NA here
            if (!is_defined)
            {
                var.init_numeric(na::DOUBLE);
                return;
            }

            // range and reshaping
            var.init_func_attr("min", mcmc_->translator.to_attribute(deftable_.c_min[k], &anchor));
            var.init_func_attr("max", mcmc_->translator.to_attribute(deftable_.c_max[k], &anchor));
            var.init_func_reshaping(str::tobool(deftable_.c_discr[k]));
        }
        catch (const std::exception &e)
        {
            const DataWide &data_wide = mcmc_->data_wide();
            msg::stop(
                std::string(e.what()) + "\nCheck " + var.err_name() + " for " +
                (house ? data_wide.err_name(*house) : data_wide.err_name(anchor)) + ".");
        }
    }

    bool Augmented::bind_house_reduction_(Variable &var,
                                          const parse::Entry &root,
                                          const Household &house,
                                          const Individual &anchor)
    {
        size_t n_operand = 0;
        const auto append = [&](const Attribute &attr)
        {
            var.init_func_attr(n_operand ? "x" + std::to_string(n_operand) : "x", attr);
            ++n_operand;
        };

        for (const parse::Entry &child : root.children)
            // a per-member data cell unfolds over the house; a household observation is
            // already one value for the whole house and enters the reduction once
            if (parse::is_data_wide(child.type) && !is_house_obs(child.value))
                for (const Individual &indiv : house)
                {
                    const Attribute &attr = mcmc_->translator.to_attribute(child, &indiv);
                    if (attr.is_constant() && na::is(attr.value()))
                        continue; // an empty cell drops out of the reduction
                    append(attr);
                }
            else
                append(mcmc_->translator.to_attribute(child, &anchor));

        return n_operand > 0;
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////  Augmented : data-bound surrogate columns  ///////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////

    bool Augmented::def_is_pure_data_(size_t k) const
    {
        for (const DefTable_::col *col : deftable_.source_cols())
        {
            const std::string &expr = (*col)[k];
            if (na::is(expr))
                continue;
            parse::Definition def(expr);
            if (def.consumes_augmented() || def.consumes_parameter() ||
                entry_has_group_risk_(def.root))
                return false;
        }
        return true;
    }

    bool Augmented::entry_has_group_risk_(const parse::Entry &node)
    {
        if (parse::is_group_risk(node.type))
            return true;
        for (const parse::Entry &child : node.children)
            if (entry_has_group_risk_(child))
                return true;
        return false;
    }

    void Augmented::assert_status_obs_consistent_()
    {
        if (!is_any_status_obs())
            return;

        const size_t k = status_obs_def_;
        const std::string &astr = deftable_.c_name[k];

        // a status observation's weights must be constant data: a status fixed in the
        // data contributes a constant to the likelihood, dropped during inference.
        // NOTE: outcome k resolves to infected status (INFECTED + k); the matching
        // augmented/observation block is selected at simulation (status-specific if
        // declared, otherwise the generic one), so no block count is enforced here.
        if (!def_is_pure_data_(k))
            msg::stop(
                "Status observation " + msg::val(astr) + " must depend only on data or numeric constants.");
    }

    void Augmented::materialize_surrogate_columns(bool for_export)
    {
        if (list_surrogate.empty())
            return;

        DataWide &data_wide = mcmc_->data_wide();

        // build the surrogate Variables (inference home, shared with the simulation init)
        // and point each standing-in data column's read-surface at its Variable, so @D
        // readers resolve the surrogate value without a data copy.
        house_varia.assign(data_wide.report_house.total,
                           std::vector<Variable>(list_obs_house.size()));

        for (const Surrogate &surr : list_surrogate)
        {
            const std::string &astr = deftable_.c_name[surr.row];
            // a surrogate column must resolve from data: a latent, parameter or group-risk
            // source has no fixed value during inference (it would impose an indicator
            // constraint).
            if (!def_is_pure_data_(surr.row))
            {
                // exporting the data view is not inferring: the column has no derived
                // value to report and keeps its datum, the run wall still refuses it
                if (for_export)
                    continue;
                msg::stop(
                    "Cannot infer from " +
                    err_type_name(kind_to_err_type(deftable_.c_kind[surr.row])) + " " +
                    msg::val(astr) + ": it must resolve from data alone.\n"
                                     "A definition reading a latent variable or a parameter is only "
                                     "meaningful in simulation.");
            }

            for (Household &house : data_wide)
                build_surrogate_column_(surr, house, true);
        }
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////  Augmented : init missing attribute  ////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////

    void Augmented::clear_all()
    {
        if (khold_.empty())
            return;

        for (const Household &house : mcmc_->data_wide())
            for (const Individual &indiv : house.infected_indivs())
                clear_definitions(indiv);
    }

    void Augmented::clear_definitions(const Individual &indiv)
    {
        for (Variable &avar : varia_of_(indiv).vec)
            avar.clear();
    }

    void Augmented::init_definitions_for_simulation(const Individual &indiv)
    {
        for (size_t k : helper_materialize_definition_order_(indiv))
            if (varia_of_(indiv).vec[khold(deftable_.c_name[k])].is_empty())
                helper_init_attributes_(indiv, k, true);
    }

    std::vector<size_t> Augmented::helper_materialize_definition_order_(const Individual &indiv)
    {
        int status = indiv.status();

        // if known status, then load its definition
        const auto it = map_status_to_block_.find(status);
        if (it != map_status_to_block_.end())
        {
            const size_t i = it->second;
            varia_of_(indiv).idx = i;
            return status_blocks_[i].order;
        }

        // otherwise if it exists fallback to generic definition
        const auto generic_it = map_status_to_block_.find(status::MISSING);
        if (generic_it != map_status_to_block_.end())
        {
            const size_t i = generic_it->second;
            map_status_to_block_.emplace(status, i);
            varia_of_(indiv).idx = i;
            return status_blocks_[i].order;
        }

        if (indiv.is_infected_status())
            msg::stop(
                "No augmented definitions were found for " +
                mcmc_->data_wide().err_name(indiv) + " with status " + std::to_string(status) + ".");

        msg::stop(
            "Generic augmented definition for non-infectious " +
            mcmc_->data_wide().err_name(indiv) + " was missing.");
    }

    void Augmented::helper_init_attributes_(const Individual &indiv, size_t k, bool simulation)
    {
        const std::string &astr = deftable_.c_name[k];
        size_t a = khold(astr);

        try
        {
            Variable &avar = varia_of_(indiv).vec[a];

            if (str::tobool(deftable_.c_const[k]))
            {
                // placeholder for observation node if generated by simulation
                if (simulation &&
                    (deftable_.c_kind[k] == Kind::Member || deftable_.c_kind[k] == Kind::Status))
                    avar.init_numeric(0.);
                else
                    avar.init_numeric(
                        mcmc_->translator.to_numeric(deftable_.c_value[k], &indiv).value());
                avar.enforce_user_const(); // enforce user constant choice
            }
            // simulation-only override: a sampled augmented (model present) carrying a value is
            // pinned to that data-bound value in generative simulation (the inference model is a
            // nuisance law, unfit to generate). to_numeric evaluates the constant value eagerly --
            // the index reads its infection time before any lazy pass, so a live function would
            // still read 0 -- and the early return skips the model so restart leaves it untouched.
            else if (simulation && !na::is(deftable_.c_value[k]) && !na::is(deftable_.c_model[k]))
            {
                avar.init_numeric(
                    mcmc_->translator.to_numeric(deftable_.c_value[k], &indiv).value());
                return;
            }
            else
                build_definition_variable_(avar, k, indiv, nullptr);

            // initialize model attributes
            helper_init_distr_attributes_(indiv, k);

            // assign proposal if variable to be sampled
            // (arguments are folded into the proposal call, if any); simulation
            // draws from the model and never proposes, so its data-bound args
            // (e.g. an observation column) need not resolve here
            if (!simulation && !avar.is_constant())
            {
                parse::Definition prop(deftable_.c_proposal[k]);
                dict<Attribute> prop_args;
                for (const auto &child : prop.root.children)
                    // per-member context lets proposal args fold data cells
                    prop_args.emplace(child.key, mcmc_->translator.to_numeric(child, &indiv));
                avar.init_proposal(prop.root.func, prop_args);
            }

            // check if satisfies discrete requirement, if any
            if (user_discrete_set_.count(a) && !avar.is_discrete())
                msg::stop(
                    "Definition of " + avar.err_name() +
                    " did not return a discrete object as requested by model.");
        }
        catch (const std::exception &e)
        {
            msg::stop(
                std::string(e.what()) +
                "\nCheck augmented " + msg::val(astr) + " for " +
                mcmc_->data_wide().err_name(indiv) + " (status " + std::to_string(indiv.status()) + ").");
        }
    }

    void Augmented::helper_init_distr_attributes_(const Individual &indiv, size_t k)
    {
        const std::string &astr = deftable_.c_name[k];
        Variable &avar = varia_of_(indiv).vec[khold(astr)];
        parse::Definition parser(deftable_.c_model[k]);

        avar.init_distr(parser.root.func);

        if (!avar.has_distr())
            return;
        for (const auto &child : parser.root.children)
            avar.init_distr_attr(child.key, mcmc_->translator.to_attribute(child, &indiv));
        avar.init_distr_attr("min", mcmc_->translator.to_attribute(deftable_.c_min[k], &indiv));
        avar.init_distr_attr("max", mcmc_->translator.to_attribute(deftable_.c_max[k], &indiv));
        avar.init_distr_reshaping(str::tobool(deftable_.c_discr[k]));
    }

    /////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////  Augmented : user reset  ////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////

    void Augmented::sample_obs(const Individual &indiv, rnd_gen &prng)
    {
        // draw each observation into its Variable; @D consumers read it straight from
        // there (the data column, if any, is never overwritten).
        for (size_t a : list_obs_models(indiv))
            varia_of_(indiv).vec[a].draw_value(prng);
    }

    void Augmented::simulate_infected_status(Household &house,
                                             Individual &indiv,
                                             rnd_gen &prng)
    {
        StatusLevel new_status = status::INFECTED;

        if (is_any_status_obs())
        {
            Variable &svar = varia_of_(indiv).vec[status_obs_key_];
            // lazy materialize status observation on first use
            if (svar.is_empty())
                helper_init_attributes_(indiv, status_obs_def_, true);
            svar.draw_value(prng);

            // NOTE: I am assuming value draws are integer in [0, max - 1]
            new_status += static_cast<StatusLevel>(svar.value());
        }

        house.set_new_status_(indiv, new_status);
    }

    void Augmented::restart(rnd_gen &prng)
    {
        if (!is_any_sample())
            return;

        for (const Household &house : mcmc_->data_wide())
            for (const Individual &indiv : house.infected_indivs())
                restart(indiv, prng);
    }

    [[noreturn]] void Augmented::throw_obs_init_error(const Individual &indiv) const
    {
        const StatusBlock_ &block = block_of_(indiv);

        for (size_t k : block.order)
        {
            if (deftable_.c_kind[k] != Kind::Member)
                continue;

            const Variable &ovar = varia_of_(indiv).vec[khold(deftable_.c_name[k])];
            if (std::isfinite(ovar.guarded_log_pdf()))
                continue;

            // domain: loc + [min, max]
            std::vector<std::string> augm_producers, param_producers;
            for (const parse::Entry &child : parse::Definition(deftable_.c_model[k]).root.children)
                if (child.key == "loc")
                {
                    child.append_augmented_producers(augm_producers);
                    child.append_parameter_producers(param_producers);
                }
            for (const std::string &col : {deftable_.c_min[k], deftable_.c_max[k]})
            {
                parse::Definition def(col);
                augm_producers += def.list_augmented_producers();
                param_producers += def.list_parameter_producers();
            }

            if (augm_producers.size() > 0)
                msg::stop(
                    "Exceeded attempts to initialize augmented variables of " +
                    mcmc_->data_wide().err_name(indiv) + " due to observations.\n" +
                    "Widen augmented model of " + str::format(str::unique(augm_producers), false, 3) +
                    ", revise the observation domain of " + msg::val(deftable_.c_name[k]) +
                    ", or check the observed individual data.");

            // non-constant parameter in domain cannot be fixed by per-individual resampling
            for (const std::string &pstr : str::unique(param_producers))
                if (!mcmc_->param.at(mcmc_->param.khold(pstr)).is_constant())
                    msg::stop(
                        "Failed to initialize due to observations of " + mcmc_->data_wide().err_name(indiv) + ".\n" +
                        "Try to remove free parameter dependencies in observation domain of " + msg::val(deftable_.c_name[k]) + ".");
        }

        // otherwise observed data are just unreachable
        msg::stop(
            "Failed to initialize due to observations of " + mcmc_->data_wide().err_name(indiv) + ".\n" +
            "The observed values lie outside of respective domain.\n" +
            "Check the individual data or revise the domain.");
    }

    void Augmented::finalize()
    {
        if (khold_.empty())
            return;

        for (const Household &house : mcmc_->data_wide())
            for (const Individual &indiv : house.infected_indivs())
                for (Variable &avar : varia_of_(indiv).vec)
                    avar.finalize();
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////  Augmented : log-likelihood  //////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////

    double Augmented::commit_augm_log_lik(const Individual &indiv)
    {
        const std::vector<size_t> &models = list_augm_models(indiv);

        if (models.empty())
            return 0.;
        double new_ll = 0.0;
        // guarded: a moved parameter (e.g. a param feeding a model loc/min/max) may push a
        // value out of its shifted support; that must read as -Inf, not a stale finite density.
        for (size_t a : models)
            new_ll += varia_of_(indiv).vec[a].guarded_log_pdf();
        varia_of_(indiv).augm_log_lik.set(new_ll);
        return new_ll;
    }

    double Augmented::commit_augm_log_lik_ratio(const Individual &indiv, size_t a)
    {
        const StatusBlock_ &block = block_of_(indiv);
        if (block.consumers[a].empty())
            return 0.;
        double delta_ll = 0;
        // NOTE: consumers[a] includes 'a' itself, so model of 'a' is counted
        for (size_t b : block.consumers[a])
            delta_ll += varia_of_(indiv).vec[b].log_ratio();
        varia_of_(indiv).augm_log_lik.increment(delta_ll);
        return delta_ll;
    }

    double Augmented::sum_guarded_augm_log_lik(const Individual &indiv, size_t a) const
    {
        double ll = 0.;
        // NOTE: consumers[a] includes 'a' itself, so model of 'a' is counted
        for (size_t b : block_of_(indiv).consumers[a])
            ll += varia_of_(indiv).vec[b].guarded_log_pdf();
        return ll;
    }

    double Augmented::commit_augm_log_lik_ratio_with_consumers(const Individual &indiv,
                                                               size_t a,
                                                               double L_old)
    {
        // 'a' has already been proposed: its model-consumers now read the moved parameter live.
        // Compare their guarded model sum against the pre-move snapshot L_old.
        const double delta_ll = sum_guarded_augm_log_lik(indiv, a) - L_old;
        varia_of_(indiv).augm_log_lik.increment(delta_ll);
        return delta_ll;
    }

    double Augmented::commit_obs_log_lik(const Individual &indiv)
    {
        if (list_obs_models(indiv).empty())
            return 0.;
        const double new_ll = eval_log_obs(indiv);
        varia_of_(indiv).obs_log_lik.set(new_ll);
        return new_ll;
    }

    double Augmented::commit_obs_log_lik_ratio(const Individual &indiv, size_t a)
    {
        const StatusBlock_ &block = block_of_(indiv);
        // skip unless the moved augmented 'a' parameterizes some observation node
        if (block.obs_consumers[a].empty())
            return 0.;

        // log-proba must be recomputed (observation value is just a datum)
        double new_ll = 0.0;
        for (size_t b : block.obs_models)
            new_ll += varia_of_(indiv).vec[b].guarded_log_pdf();
        const double delta_ll = new_ll - varia_of_(indiv).obs_log_lik.value();
        varia_of_(indiv).obs_log_lik.set(new_ll);
        return delta_ll;
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////  Augmented : sample draws  ///////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////

    namespace
    {
        struct AugmDraws_
        {
            std::vector<std::string> ids;
            std::vector<std::string> house_ids;
            std::vector<const Individual *> ptrs;

            void store(const DataWide &data_wide, const Household &house, const Individual &indiv)
            {
                ids.emplace_back(data_wide.id(indiv));
                house_ids.emplace_back(data_wide.house_id(house));
                ptrs.emplace_back(&indiv);
            }
        };
    }

    std::vector<const Individual *> Augmented::init_individual_sample(
        const std::vector<std::string> &id_draws,
        const std::vector<std::string> &house_draws,
        rnd_gen &prng) const
    {
        if (id_draws.size() == 1 && id_draws[0] == ".none" && house_draws.empty())
            return {};

        if (!is_any_sample())
        {
            if (!id_draws.empty() || !house_draws.empty())
                msg::warn(
                    "Ignored" + msg::arg("id_draws") +
                    ": no augmented variables to sample were defined.");
            return {};
        }

        bool is_default = id_draws.size() == 0 && house_draws.empty();

        AugmDraws_ seen;
        std::vector<const Household *> eligible;

        // explicit individual selections
        // store eligible households for default
        for (const Household &house : mcmc_->data_wide())
            if (house.n_infected_status() > 0)
            {
                if (is_default)
                {
                    eligible.push_back(&house);
                    continue;
                }

                if (contains(house_draws, mcmc_->data_wide().house_id(house)))
                    for (const Individual &indiv : house.infected_indivs())
                        seen.store(mcmc_->data_wide(), house, indiv);
                else
                    for (const Individual &indiv : house.infected_indivs())
                        if (contains(id_draws, mcmc_->data_wide().id(indiv)))
                            seen.store(mcmc_->data_wide(), house, indiv);
            }

        // raise warnings for draws not found in the data
        std::vector<std::string> missed_id = id_draws - seen.ids;
        std::vector<std::string> missed_house = house_draws - seen.house_ids;
        if (!missed_id.empty())
            msg::warn(
                "Ignored draws of non infected " + msg::qty(missed_id) +
                "individual{?s} " + str::format(missed_id));
        if (!missed_house.empty())
            msg::warn(
                "Ignored draws for " + msg::qty(missed_house) +
                "household{?s} with no infected members " + str::format(missed_house));

        // default: follow one random household with at least one infected
        if (!eligible.empty())
        {
            std::uniform_int_distribution<size_t> pick(0, eligible.size() - 1);
            const Household &house = *eligible[pick(prng)];
            for (const Individual &indiv : house.infected_indivs())
                seen.store(mcmc_->data_wide(), house, indiv);
        }

        return seen.ptrs;
    }
}
