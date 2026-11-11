/* housepi - Epidemiological Analysis of Household Outbreaks.
 * Copyright 2024-2026 Institut Pasteur (Paris), France.
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#include <algorithm>
#include <atomic>
#include <chrono>
#include <cmath>
#include <exception>
#include <functional>
#include <memory>
#include <random>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

#include "tools/keys.h"
#include "tools/logproba.h"
#include "tools/messages.h"
#include "tools/nan.h"
#include "tools/parallel.h"
#include "tools/parser.h"
#include "tools/random.h"
#include "tools/string.h"
#include "tools/vector.h"
#include "augmented.h"
#include "data.h"
#include "data/dataset.h"
#include "data/household.h"
#include "data/individual.h"
#include "data/long_dataset.h"
#include "distribution.h"
#include "function.h"
#include "infec_profile.h"
#include "mcmc.h"
#include "mcmc/base.h"
#include "mcmc/settings.h"
#include "model.h"
#include "parameter.h"
#include "rates.h"
#include "registry/functions.h"
#include "registry/models.h"
#include "translator.h"
#include "variable.h"

namespace housepi
{
    //////////////////////////////  registry::RegModels  //////////////////////////////

    namespace registry
    {
        RegModels &RegModels::instance()
        {
            static RegModels models;
            return models;
        }
    }

    //////////////////////////////  Translator  //////////////////////////////

    const Attribute &Translator::to_attribute(const parse::Entry &node,
                                              const Individual *indiv_ptr,
                                              const Household *house_ptr) const
    {
        storage_ = Attribute();

        if (!node.is_terminal())
        {
            helper_store_read_function_(node, indiv_ptr, house_ptr);
        }
        else if (parse::is_numeric(node.type))
        {
            storage_.init_numeric(str::tod(node.value, "argument with numeric type"));
        }
        else if (parse::is_parameter(node.type))
        {
            return param_.at(param_.consult(node.value));
        }
        else if (!node.is_empty())
        {
            // a household value reads what holds one value for the whole house
            if (!indiv_ptr && house_ptr)
            {
                if (reads_member_cell_(node))
                    msg::stop(
                        "Per-member " + node.err_str() + " is read outside a reduction.\n" +
                        "Fold it over the household members with min(), max(), sum() or prod(), e.g. " +
                        msg::code("max(" + node.value + ")") + ".");
                indiv_ptr = &house_ptr->at(0);
            }

            if (!indiv_ptr)
            {
                std::string sms = "Individual specific definitions as " + node.err_str() +
                                  " are not supported here.";
                // add tip if the reference resolved to a latent variable
                if (parse::is_latent(node.type) && node.is_terminal())
                    sms += "\nDid you forget rlang injection operator " + msg::code("!!" + node.value) + "?";
                msg::stop(sms);
            }

            if (parse::is_dataset(node.type))
            {
                // wire the column's read-cell: its producer points the cell at itself
                // the reader needs no rule of its own
                storage_.init_pointer(
                    &data_.at(data_.dataset().consult(node.value), *indiv_ptr));
                storage_.name(node.value, ErrType::Data);
            }
            else if (parse::is_latent(node.type))
            {
                return augm_.at(augm_.consult(node.value), *indiv_ptr);
            }
            else if (parse::is_group_rate(node.type))
            {
                return group_rate_.at(group_rate_.consult(node.value), *indiv_ptr);
            }
        }

        return storage_;
    }

    void Translator::helper_store_read_function_(const parse::Entry &node,
                                                 const Individual *indiv_ptr,
                                                 const Household *house_ptr) const
    {
        std::unique_ptr<Function> func = registry::create_function(
            node.func, ErrType::Function, registry::func_distr_search_order(false));
        if (is_house_fold(node, house_ptr))
            func->set_eval_as_fold();

        Attribute parent;
        parent.init_function(std::move(func));
        const bool is_defined = bind_arguments(parent, node, indiv_ptr, house_ptr);

        storage_ = Attribute();
        if (is_defined)
            storage_ = parent;
        else
            storage_.init_numeric(na::DOUBLE);
    }

    namespace
    {
        /// Whether `f` is an explicit reduction, folding household members.
        bool is_reduction_(std::string_view f)
        {
            return f == "sum" || f == "prod" || f == "min" || f == "max";
        }
    }

    bool Translator::reads_member_cell_(const parse::Entry &node) const
    {
        if (node.is_terminal())
            return (parse::is_dataset(node.type) && !augm_.def().has_row(node.value, augm::Kind::House)) ||
                   parse::is_latent(node.type) ||
                   parse::is_group_rate(node.type);
        // an inner reduction folds its own operands
        if (is_reduction_(node.func))
            return false;
        return reads_member_cell_(node.children);
    }

    bool Translator::reads_member_cell_(const std::vector<parse::Entry> &nodes) const
    {
        return std::any_of(nodes.begin(), nodes.end(),
                           [this](const parse::Entry &child)
                           { return reads_member_cell_(child); });
    }

    bool Translator::is_house_fold(const parse::Entry &node,
                                   const Household *house_ptr) const
    {
        return house_ptr && !node.is_terminal() && is_reduction_(node.func) &&
               reads_member_cell_(node.children);
    }

    bool Translator::bind_arguments(Attribute &parent,
                                    const parse::Entry &node,
                                    const Individual *indiv_ptr,
                                    const Household *house_ptr) const
    {
        const bool is_fold = is_house_fold(node, house_ptr);

        // operands keep their own key, unless the fold makes the call variadic
        size_t n_operand = 0;
        const auto append = [&](const std::string &key, const Attribute &attr)
        {
            ++n_operand;

            parent.init_func_attr(!is_fold ? key : n_operand > 1 ? "x" + std::to_string(n_operand)
                                                                 : "x",
                                  attr);
        };
        // within a household, an operand carrying no value:
        // * drops out of a fold;
        // * leaves any other call undefined.
        const auto is_valueless = [&](const Attribute &attr)
        { return house_ptr && attr.is_constant() && attr.is_na(); };

        for (const parse::Entry &child : node.children)
            if (is_fold && reads_member_cell_(child))
                for (const Individual &indiv : *house_ptr)
                {
                    const Attribute &attr = to_attribute(child, &indiv, house_ptr);
                    if (!is_valueless(attr))
                        append(child.key, attr);
                }
            else
            {
                const Attribute &attr = to_attribute(child, indiv_ptr, house_ptr);
                if (!is_valueless(attr))
                    append(child.key, attr);
                else if (!is_fold)
                    return false;
            }

        return !is_fold || n_operand > 0;
    }

    const Attribute &Translator::to_attribute(const std::string &expr,
                                              const Individual *indiv_ptr) const
    {
        return to_attribute(parse::Definition(expr).root, indiv_ptr);
    }

    const Attribute &Translator::to_numeric(const parse::Entry &node,
                                            const Individual *indiv_ptr) const
    {
        storage_ = to_attribute(node, indiv_ptr);
        if (!storage_.is_empty() && !storage_.is_constant())
            msg::stop(
                "Definition " + node.err_str() +
                " could not be interpreted as a numeric constant.");
        return storage_;
    }

    dict<Attribute> Translator::to_numeric_args(const parse::Entry &node,
                                                const Individual *indiv_ptr) const
    {
        dict<Attribute> args;
        for (const parse::Entry &child : node.children)
            args.emplace(child.key, to_numeric(child, indiv_ptr));
        return args;
    }

    void Translator::bind_distr(Variable &var,
                                const parse::Entry &model,
                                const parse::Entry &min,
                                const parse::Entry &max,
                                bool discr,
                                const Individual *indiv_ptr,
                                bool numeric) const
    {
        const auto read = [&](const parse::Entry &node) -> const Attribute &
        { return numeric ? to_numeric(node, indiv_ptr) : to_attribute(node, indiv_ptr); };
        for (const parse::Entry &child : model.children)
            var.init_distr_attr(child.key, read(child));
        var.init_distr_attr("min", read(min));
        var.init_distr_attr("max", read(max));
        var.reshape_distribution(discr);
    }

    void Translator::mark_param_produces_transm_cache(const parse::Entry &node) const
    {
        if (!parse::is_parameter(node.type))
            return;
        size_t p = param_.consult(node.value);
        if (!param_.at(p).is_constant())
            param_.mark_produces_transm_cache(p);
    }

    //////////////////////////////  Posterior  //////////////////////////////

    namespace
    {
        /// Individuals selected for latent draws, with their ids and household ids.
        struct LatentDraws_
        {
            std::vector<std::string> ids;
            std::vector<std::string> hh_ids;
            std::vector<const Individual *> ptrs;

            void store(const data::Dataset &dataset,
                       const Household &house,
                       const Individual &indiv)
            {
                ids.emplace_back(dataset.id(indiv));
                hh_ids.emplace_back(dataset.hh_id(house));
                ptrs.emplace_back(&indiv);
            }
        };

        /// Resolve individual and household ids for latent draws, and return tracked pointers.
        /// @note Draw lists are assumed to hold unique ids.
        std::vector<const Individual *> helper_init_individual_sample(const Augmented &augm,
                                                                      const data::Dataset &dataset,
                                                                      const std::vector<std::string> &id_draws,
                                                                      const std::vector<std::string> &house_draws,
                                                                      RndGen &prng)
        {
            if (id_draws.size() == 1 && id_draws[0] == ".none" && house_draws.empty())
                return {};

            if (!augm.scheme().is_any_sample())
            {
                if (!id_draws.empty() || !house_draws.empty())
                    msg::warn(
                        "Ignored " + msg::arg("id_draws") +
                        ": no latent variables to sample were defined.");
                return {};
            }

            bool is_default = id_draws.size() == 0 && house_draws.empty();

            LatentDraws_ seen;
            std::vector<std::string> known_id, known_house;
            std::vector<const Household *> eligible;

            // collect explicit selections and the identifiers they name, or eligible households for the default
            for (const Household &house : dataset)
            {
                if (is_default)
                {
                    if (house.n_infected_status() > 0)
                        eligible.push_back(&house);
                    continue;
                }

                const bool is_house = vect::contains(house_draws, dataset.hh_id(house));
                if (is_house)
                    known_house.emplace_back(dataset.hh_id(house));
                for (const Individual &indiv : house)
                {
                    const bool is_indiv = vect::contains(id_draws, dataset.id(indiv));
                    if (is_indiv)
                        known_id.emplace_back(dataset.id(indiv));
                    if ((is_house || is_indiv) && indiv.is_infected_status())
                        seen.store(dataset, house, indiv);
                }
            }

            // a draw naming no member or household stops, a known one without infection warns
            const std::vector<std::string> unknown =
                vect::concatenate(vect::difference(id_draws, known_id), vect::difference(house_draws, known_house));
            if (!unknown.empty())
                msg::stop(
                    "Unknown " + msg::arg("id_draws") + " " + msg::qty(unknown) +
                    "identifier{?s} " + str::format(unknown));
            const std::vector<std::string> missed_id = vect::difference(id_draws, seen.ids);
            const std::vector<std::string> missed_house = vect::difference(house_draws, seen.hh_ids);
            if (!missed_id.empty())
                msg::warn(
                    "Ignored draws of non infected " + msg::qty(missed_id) +
                    "individual{?s} " + str::format(missed_id));
            if (!missed_house.empty())
                msg::warn(
                    "Ignored draws for " + msg::qty(missed_house) +
                    "household{?s} with no infected members " + str::format(missed_house));

            // default: follow one random household with an infected member
            if (!eligible.empty())
            {
                std::uniform_int_distribution<size_t> pick(0, eligible.size() - 1);
                const Household &house = *eligible[pick(prng)];
                for (const Individual &indiv : house.infected_indivs())
                    seen.store(dataset, house, indiv);
            }

            return seen.ptrs;
        }
    }

    /// Initialize the parameter draws of the posterior.
    void Posterior::init_param_draws_(const Parameter &param,
                                      Settings::Step sample_size)
    {
        for (size_t p : param.scheme().sample)
        {
            param_values.add(param.consult(p), sample_size);
            param_priors.add("param_lprior_" + param.consult(p), sample_size);
        }
        // function-derived parameters: value trace only, no prior
        for (size_t p : param.scheme().func)
            param_values.add(param.consult(p), sample_size);
    }

    void Posterior::init_house_draws_(const Data &data,
                                      Settings::Step sample_size)
    {
        house_log_lik.cols.assign(data.n_eq_classes(), std::vector<double>(sample_size, na::DOUBLE));

        for (const Household &house : data.dataset())
        {
            house_log_lik.names.push_back(data.dataset().hh_id(house));
            house_log_lik.source.push_back(data.state.eq_class(house));
        }
    }

    /// Initialize the latent draws of the posterior.
    void Posterior::init_latent_draws_(const Augmented &augm,
                                       const data::Dataset &dataset,
                                       const Settings &settings,
                                       Settings::Step sample_size)
    {
        tracked_indivs = helper_init_individual_sample(
            augm, dataset, settings.id_draws, settings.house_draws, settings.prng());

        for (const Individual *indiv : tracked_indivs)
            for (size_t a : augm.list_draws(*indiv))
                latent_values.add(dataset.id(*indiv) + "_" + augm.consult(a), sample_size);
    }

    void Posterior::init(const Workspace &chain, const Settings &settings)
    {
        param_values = {};
        param_priors = {};
        latent_values = {};
        house_log_lik = {};
        tracked_indivs.clear();

        const Settings::Step sample_size = settings.c_sample_size * settings.n_chains;

        init_param_draws_(chain.param, sample_size);
        init_house_draws_(chain.data, sample_size);
        init_latent_draws_(chain.augm, chain.data.dataset(), settings, sample_size);
    }

    void Posterior::store(const Workspace &chain,
                          const Settings &settings,
                          Settings::Step local)
    {
        const size_t col = chain.id * settings.c_sample_size + local;

        // parameter values and priors
        const Parameter &param = chain.param;
        size_t i = 0;
        for (size_t p : param.scheme().sample)
        {
            const Variable &var = param.at(p);
            param_values.cols[i][col] = var.value();
            param_priors.cols[i++][col] = var.log_pdf();
        }
        // function-derived parameters: store value, no prior
        for (size_t p : param.scheme().func)
            param_values.cols[i++][col] = param.at(p).value();

        // pointwise household log-likelihood
        const Data &data = chain.data;
        for (size_t R = 0; R < data.n_eq_classes(); ++R)
            house_log_lik.cols[R][col] = chain.house_log_lik(data.eq_class_rep(R));

        // latent variable draws, member by member
        const Augmented &augm = chain.augm;
        size_t c = 0;
        for (const Individual *indiv : tracked_indivs)
            for (size_t a : augm.list_draws(*indiv))
                latent_values.cols[c++][col] = chain.aug(a, *indiv);
    }

    //////////////////////////////  Workspace  //////////////////////////////

    Workspace::Workspace(std::unique_ptr<ModelBase> model_ptr,
                         const Definitions &defs,
                         const Settings &settings,
                         const Data *shared)
        : model(std::move(model_ptr)),
          defs(defs),
          data(*model, shared),
          param(model->param_keys(), defs.param),
          rates(model->group_rate_keys()),
          infec_profile(augm),
          augm(model->latent_keys(),
               defs.latent,
               data,
               param,
               infec_profile,
               defs.status_process,
               translator),
          translator(param, augm, data, rates.group())
    {
        model->bind(*this, settings);
        model->init_model_requirements();
    }
    Workspace::~Workspace() = default;

    void Workspace::build(Definitions::Slot slot)
    {
        using Slot = Definitions::Slot;
        // a rebuilt definition invalidates the last evaluation
        data.state.set_evaluated(false);
        switch (slot)
        {
        case Slot::Param:
            param.build(translator);
            break;
        case Slot::GroupRate:
            data.set_eq_class_keys(Data::Source::GroupRate,
                                   rates.build_group(data.dataset(), param, defs.group_rate));
            break;
        case Slot::ContactRate:
            data.set_eq_class_keys(Data::Source::ContactRate,
                                   rates.build_contact(data.dataset(), param, defs.contact_rate));
            break;
        case Slot::InfecProfile:
            infec_profile.build(translator, defs.infec_profile);
            break;
        case Slot::Latent:
            augm.build();
            break;
        case Slot::ParamPins:
            param.build_for_predic(defs.param_pins);
            break;
        case Slot::LatentPins:
            augm.pin_values(defs.latent_pins);
            break;
        }
    }

    void Workspace::init_equivalence_class_compression()
    {
        // one class per household when:
        // * a latent is sampled;
        // * the model declares latents, no definition provided yet.
        data.init_eq_classes(augm.scheme().is_any_sample() || (augm.def().empty() && augm.size() > 0));
        data.state.allocate_transm_log_lik();
    }

    void Workspace::assert_valid_infectivity_profile() const
    {
        infec_profile.assert_valid();
        if (!infec_profile.is_latent())
            return;

        // check per-member profiles
        for (const Household &house : data.dataset())
            for (const Individual &indiv : house.infected_indivs())
            {
                const Variable &avar = augm.at(infec_profile.latent_key, indiv);
                // @todo Check that `avar` holds a function.
                infec_profile.assert_zero_at_origin(*avar.function());
            }
    }

    //////  Workspace : Getter  //////

    std::vector<std::string> Workspace::get_static_data_names() const
    {
        // exclude exported latent and observation-node columns
        // an observation node is generated by the simulation: the generated value wins over static data
        // a deterministic observation is produced too, as a column of its own: it never surfaces as static data
        std::vector<std::string> produced = vect::concatenate(
            augm.consult(augm.scheme().list_latent_export), augm.consult(augm.scheme().list_obs_member));
        for (const augm::DeterministicObs &dobs : augm.scheme().list_deterministic_obs)
            if (dobs.is_own_column())
                produced.push_back(augm.def()[dobs.row].name);
        return vect::difference(data.dataset().all(), produced);
    }

    std::vector<std::string> Workspace::get_col_names() const
    {
        std::vector<std::string> output = {"id", "hh_id"};
        vect::append(output, get_static_data_names());
        vect::append(output, {"status"});

        output.reserve(output.size() +
                       augm.scheme().list_latent_export.size() +
                       augm.scheme().list_obs_member.size() +
                       augm.scheme().list_deterministic_obs.size());

        for (auto idx : augm.scheme().list_latent_export)
            output.push_back(augm.consult(idx));

        for (auto idx : augm.scheme().list_obs_member)
            output.push_back(augm.consult(idx));

        for (const augm::DeterministicObs &dobs : augm.scheme().list_deterministic_obs)
            if (dobs.is_own_column())
                output.push_back(augm.def()[dobs.row].name);

        return output;
    }

    std::vector<size_t> Workspace::get_row_offsets() const
    {
        size_t n_houses = data.dataset().report_house.total;
        std::vector<size_t> output(n_houses + 1, 0);
        for (size_t H = 0; H < n_houses; ++H)
            output[H + 1] = output[H] + data.dataset().at(H).size();
        return output;
    }

    //////////////////////////////  MCMCBase  //////////////////////////////

    MCMCBase::MCMCBase(std::unique_ptr<ModelBase> model)
    {
        workspaces_.push_back(std::make_unique<Workspace>(std::move(model), defs, settings));
    }
    MCMCBase::~MCMCBase() = default;

    //////////////////////////////  MCMC  //////////////////////////////

    MCMC::MCMC(std::string_view model_name)
        : MCMCBase(registry::RegModels::instance().create(model_name)),
          model_name_(model_name) {}

    void MCMC::detach_data_store()
    {
        Workspace &chain = this->chain();
        if (!chain.data.is_shared())
            return;
        Data fresh(*chain.model);
        // model requests outlive the data
        fresh.dataset().copy_requests(chain.data.dataset());
        chain.data.share(fresh);
    }

    //////  Evaluation  //////

    std::vector<double> MCMC::evaluate_log_lik()
    {
        Workspace &chain = this->chain();
        // using '0 chains' as a state variable
        settings.n_chains = 0;
        init_for_evaluation_(chain);
        chain.param.load_sample_values(0);
        chain.param.draw_sample(settings.prng());
        chain.augm.draw_sample(settings.prng());
        commit_total_log_posterior_(chain, all_terms_(chain));
        chain.data.state.set_evaluated(true);

        // each household reads the terms of its class representative
        std::vector<double> output;
        output.reserve(chain.data.dataset().report_house.total);
        for (const Household &house : chain.data.dataset())
            output.push_back(chain.house_log_lik(chain.data.eq_class_rep(house)));
        return output;
    }

    //////  Inference  //////

    void MCMC::prepare_inference()
    {
        Workspace &chain = this->chain();
        init_monte_carlo_(chain);
        init_for_evaluation_(chain);
        init_param_routes_(chain);
        posterior.init(chain, settings);
        settings.n_done = 0;
        build_workspaces_();
    }

    double MCMC::sample()
    {
        const auto start = std::chrono::steady_clock::now();
        const size_t n_workspaces = workspaces_.size();
        std::atomic<size_t> n_completed{0};
        // workspace w runs chains w, w + n_workspaces, ... one after another
        const parallel::Failure failure = parallel::for_units(
            n_workspaces, [&](size_t w)
            {
                Workspace &chain = *workspaces_[w];
                for (Settings::Chain id = w; id < settings.n_chains; id += n_workspaces)
                {
                    chain.id = id;
                    init_monte_carlo_chain_(chain);
                    init_guarded_commit_total_(chain);
                    run_chain_(chain);
                }
                n_completed.fetch_add(1, std::memory_order_relaxed); },
            n_workspaces, settings.budget.inner, *settings.stop_flag);

        // keep the workspace that failed
        // otherwise the one that started the highest chain
        const auto kept = failure.error
                              ? workspaces_.begin() + failure.index
                              : std::max_element(workspaces_.begin(), workspaces_.end(),
                                                 [](const auto &a, const auto &b)
                                                 { return a->id < b->id; });
        std::swap(workspaces_.front(), *kept);
        workspaces_.resize(1);
        chain().id = 0;
        if (failure.error)
            std::rethrow_exception(failure.error);
        // a workspace stopped by another engine or the user
        if (n_completed.load() < n_workspaces)
            throw parallel::Cancelled{};
        chain().data.state.set_evaluated(true);
        return std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
    }

    //////  Workspaces  //////

    void MCMC::build_workspaces_()
    {
        workspaces_.resize(1);
        // the first workspace raised every warning a build raises
        msg::WarningBuffer muted;
        while (workspaces_.size() < settings.budget.outer)
            workspaces_.push_back(rebuild_workspace_());
    }

    std::unique_ptr<Workspace> MCMC::rebuild_workspace_()
    {
        using Slot = Definitions::Slot;
        auto output = std::make_unique<Workspace>(
            registry::RegModels::instance().copy(model_name_, *chain().model),
            defs, settings, &chain().data);
        auto slot_order = {Slot::Param, Slot::GroupRate, Slot::ContactRate,
                           Slot::InfecProfile, Slot::Latent,
                           Slot::ParamPins, Slot::LatentPins};
        for (Slot slot : slot_order)
            output->build(slot);
        init_for_evaluation_(*output);
        init_param_routes_(*output);
        return output;
    }

    //////  Setup  //////

    void MCMC::init_monte_carlo_(Workspace &chain)
    {
        // check there is something to sample
        if (!chain.param.is_any_sample() && !chain.augm.scheme().is_any_sample())
            msg::stop(
                "Nothing to sample.\nCheck your parameter and latent definitions.");

        // whether adaptation is required
        bool is_adaptive = !chain.param.scheme().adaptive.empty() ||
                           chain.augm.scheme().is_any_adaptive();
        // initialize settings for mcmc
        settings.init_mcmc_configuration(is_adaptive);
    }

    void MCMC::init_for_evaluation_(Workspace &chain)
    {
        data::Dataset &dataset = chain.data.dataset();
        dataset.assert_not_all_missing_status();

        // materialize deterministic observations before anything reads the data store
        chain.data.reset_readcells(chain.augm.get_produced_cols());
        chain.augm.init_for_evaluation();
        dataset.assert_required_everyone(chain.data.state);
        // a sampled latent -> one class per household
        chain.data.init_eq_classes(chain.augm.scheme().is_any_sample());
        chain.data.state.allocate_cache(dataset);
        chain.data.state.allocate_transm_log_lik();

        chain.param.finalize();
        chain.rates.finalize();
        chain.augm.finalize();
        chain.infec_profile.finalize();
        chain.data.finalize();

        chain.rates.assert_identifiability(chain.param);
        chain.assert_valid_infectivity_profile();

        chain.model->init_model_for_evaluation();
        for (size_t R = 0; R < chain.data.n_eq_classes(); ++R)
            chain.model->init_house_for_evaluation(chain.data.eq_class_rep(R));
    }

    void MCMC::init_param_routes_(Workspace &chain)
    {
        const data::Dataset &dataset = chain.data.dataset();
        const augm::Scheme &as = chain.augm.scheme();
        chain.routes.assign(chain.param.size(), {});

        const std::vector<std::vector<HouseIdx>> houses_of = chain.rates.list_houses_of_param(chain.param.size());
        for (size_t p : chain.param.scheme().sample)
        {
            ParamRoute &route = chain.routes[p];
            route.terms.transm = vect::contains(chain.param.scheme().sample_cached, p);
            route.terms.latent = as.is_any_latent_model() && as.param_produces_latent_log_lik(p);
            route.terms.obs = as.is_any_obs_model() && as.param_produces_obs_log_lik(p);

            // read by rates alone -> the classes reading it
            if (route.terms.latent || route.terms.obs || houses_of[p].empty() || !chain.param.is_rate_only(p))
                continue;
            for (HouseIdx H : houses_of[p])
                route.classes.push_back(chain.data.state.eq_class(dataset.at(H)));
            vect::sort_unique(route.classes);
            // every class reads it -> no gain
            if (route.classes.size() == chain.data.n_eq_classes())
                route.classes.clear();
        }
        init_param_blocks_(chain);
    }

    void MCMC::init_param_blocks_(Workspace &chain)
    {
        chain.pblocks.clear();
        chain.in_pblock.assign(chain.param.size(), 0);
        for (std::vector<size_t> &keys : chain.rates.list_param_rw_blocks(chain.param))
        {
            ParamBlock block;
            block.keys = std::move(keys);
            bool is_sparse = true;
            for (size_t p : block.keys)
            {
                chain.in_pblock[p] = 1;
                const ParamRoute &route = chain.routes[p];
                block.route.terms |= route.terms;
                is_sparse = is_sparse && !route.classes.empty();
                block.route.classes.insert(block.route.classes.end(), route.classes.begin(), route.classes.end());
            }
            // a member recomputing every class -> the block does
            vect::sort_unique(block.route.classes);
            if (!is_sparse || block.route.classes.size() == chain.data.n_eq_classes())
                block.route.classes.clear();
            chain.pblocks.push_back(std::move(block));
        }
    }

    //////  Log-posterior  //////

    void MCMC::commit_total_log_posterior_(Workspace &chain, const LogLikTerms &terms)
    {
        if (terms.transm)
        {
            chain.data.state.save_cache();
            chain.infec_profile.update_cache();
        }

        // class representative, its terms weighted by class size
        LogLik ll = parallel::for_reduce<LogLik>(
            chain.data.n_eq_classes(),
            [&](size_t R) -> LogLik
            {
                return commit_house_log_posterior_(chain, chain.data.eq_class_rep(R), terms) * chain.data.eq_class_size(R);
            },
            settings.budget.inner);
        chain.tot_log_lik.set(ll, terms);
    }

    LogLik MCMC::commit_house_log_posterior_(Workspace &chain,
                                             const Household &house,
                                             const LogLikTerms &terms)
    {
        LogLik ll;

        // re-compute transmission cached quantities before the commit
        if (terms.transm)
            chain.model->update_house_transmission_cache(house);

        // transmission
        ll.transm = commit_transm_log_lik_(chain, house);
        // observation
        if (terms.obs)
            for (const Individual &indiv : house.infected_indivs())
                ll.obs += chain.augm.commit_obs_log_lik(indiv);
        // latent
        if (terms.latent)
            for (const Individual &indiv : house.infected_indivs())
                ll.latent += chain.augm.commit_latent_log_lik(indiv);
        return ll;
    }

    void MCMC::recover_total_log_lik_(Workspace &chain, const LogLikTerms &terms)
    {
        chain.tot_log_lik.recover(terms);

        // recover old log-likelihood values for households
        for (size_t R = 0; R < chain.data.n_eq_classes(); ++R)
        {
            const Household &house = chain.data.eq_class_rep(R);
            recover_transm_log_lik_(chain, house);

            // recover from commit_total_log_posterior_
            if (terms.latent)
                for (const Individual &indiv : house.infected_indivs())
                    recover_latent_log_lik_(chain, indiv);
            if (terms.obs)
                for (const Individual &indiv : house.infected_indivs())
                    recover_obs_log_lik_(chain, indiv);
        }

        if (terms.transm)
        {
            chain.data.state.recover_cache();
            chain.infec_profile.recover();
        }
    }

    //////  Chain  //////

    void MCMC::init_monte_carlo_chain_(Workspace &chain)
    {
        chain.param.restart(settings.prng(chain.id));
        chain.augm.restart(settings.prng(chain.id));

        // re-initialize the adaptable proposal of every variable, latent block factors from default scales
        for (size_t p : chain.param.scheme().adaptive)
            chain.param.at(p).init_proposal();
        if (chain.augm.scheme().is_any_sample())
            for (const Household &house : chain.data.dataset())
                for (const Individual &indiv : house.infected_indivs())
                {
                    for (size_t a : chain.augm.list_adaptive(indiv))
                        chain.augm.at(a, indiv).init_proposal();
                    if (block_rw::Factor *factor = chain.augm.block_factor(indiv))
                    {
                        std::vector<double> scales = {};
                        std::vector<bool> lattice = {};
                        for (size_t a : chain.augm.list_block(indiv))
                        {
                            scales.push_back(chain.augm.at(a, indiv).block_scale());
                            lattice.push_back(chain.augm.at(a, indiv).block_on_lattice());
                        }
                        factor->init(scales, std::move(lattice));
                    }
                }
    }

    void MCMC::init_guarded_commit_total_(Workspace &chain)
    {
        chain.data.state.save_cache();
        chain.infec_profile.update_cache();

        // single serial pass
        LogLik ll;
        for (size_t R = 0; R < chain.data.n_eq_classes(); ++R)
            ll += init_guarded_commit_house_(chain, chain.data.eq_class_rep(R)) * chain.data.eq_class_size(R);

        chain.tot_log_lik.set(ll, all_terms_(chain));
    }

    LogLik MCMC::init_guarded_commit_house_(Workspace &chain, Household &house)
    {
        LogLik ll;

        // observation and latent models
        // redraw a member's latents until its observation and latent log-likelihoods are finite
        if (chain.augm.scheme().is_any_obs_model() ||
            chain.augm.scheme().is_any_latent_model())
            for (const Individual &indiv : house.infected_indivs())
            {
                if (chain.augm.list_obs_models(indiv).empty() &&
                    chain.augm.list_latent_models(indiv).empty())
                    continue;

                size_t attempt = 0;
                while (!std::isfinite(chain.augm.eval_log_obs(indiv) +
                                      chain.augm.eval_log_latent(indiv)))
                {
                    ++attempt;
                    if (!chain.augm.scheme().is_any_sample() ||
                        attempt >= N_INIT_ATTEMPTS)
                    {
                        if (!std::isfinite(chain.augm.eval_log_obs(indiv)))
                            chain.augm.throw_obs_init_error(indiv);
                        chain.augm.throw_latent_init_error(indiv);
                    }
                    chain.augm.restart(indiv, settings.prng(chain.id));
                }
                ll.obs += chain.augm.commit_obs_log_lik(indiv);
                ll.latent += chain.augm.commit_latent_log_lik(indiv);
            }

        // transmission, on the latents the observation and latent models accepted
        chain.model->update_house_transmission_cache(house);
        ll.transm = commit_transm_log_lik_(chain, house);
        if (!std::isfinite(ll.transm))
            msg::stop(
                "The transmission probability of " + chain.data.dataset().str(house) +
                " is numerically zero at the initial values.\n"
                "Check the parameters and latents the transmission model reads.");

        return ll;
    }

    void MCMC::run_chain_(Workspace &chain)
    {
        //////  Adaptation  //////

        // latent blocks move one latent at a time over the warm-up, its second half shapes their factors
        const auto n_warm = static_cast<Settings::Step>(block_rw::WARMUP * settings.n_adapt);
        for (ParamBlock &block : chain.pblocks)
        {
            std::vector<double> scales = {};
            std::vector<bool> lattice = {};
            for (size_t p : block.keys)
            {
                scales.push_back(chain.param.at(p).block_scale());
                lattice.push_back(chain.param.at(p).block_on_lattice());
            }
            block.factor.init(scales, std::move(lattice));
        }
        // moment windows after the warm-up, doubling from WINDOW steps, the last running to the end
        Settings::Step window = block_rw::WINDOW, window_end = n_warm;
        const auto next_window = [&]
        {
            window_end += window;
            window *= 2;
            if (window_end + window > settings.n_adapt)
                window_end = settings.n_adapt;
        };
        next_window();
        for (Settings::Step i = 0; i < settings.n_adapt; ++i)
        {
            chain.is_warm = i < n_warm;
            // drop the moments of the first warm-up half and of each window, where a chain may still drift
            if ((n_warm > 0 && i == n_warm / 2) || i == window_end)
            {
                for (block_rw::Factor &factor : chain.augm.state.block_factors)
                    factor.forget();
                for (ParamBlock &block : chain.pblocks)
                    block.factor.forget();
                if (i == window_end)
                    next_window();
            }
            // gains falling linearly to 0 at the end of the adaptation, the block one from the end of the warm-up
            chain.gain = 1. - static_cast<double>(i) / settings.n_adapt;
            chain.block_gain = chain.is_warm ? 0. : 1. - static_cast<double>(i - n_warm) / (settings.n_adapt - n_warm);
            move_all_<true>(chain);
        }

        //////  Burn-in  //////

        for (Settings::Step i = 0; i < settings.n_burn; ++i)
            move_all_<false>(chain);

        //////  Posterior sampling  //////

        for (Settings::Step i = 0; i < settings.c_sample_size; ++i)
        {
            for (Settings::Step t = 1; t <= settings.n_thin; ++t)
                move_all_<false>(chain);
            posterior.store(chain, settings, i);
        }
    }

    template <bool Adapt>
    void MCMC::move_all_(Workspace &chain)
    {
        // move latent variables to sample, if any
        move_all_latent_<Adapt>(chain);

        // move standard model parameters, then those with cached dependencies
        // block members move together after the warm-up, one move per member
        const bool is_block = !chain.pblocks.empty() && !(Adapt && chain.is_warm);
        for (size_t p : chain.param.scheme().sample_standard)
            if (!(is_block && chain.in_pblock[p]))
                move_single_parameter_<Adapt>(chain, p);
        for (size_t p : chain.param.scheme().sample_cached)
            if (!(is_block && chain.in_pblock[p]))
                move_single_parameter_<Adapt>(chain, p);
        for (ParamBlock &block : chain.pblocks)
            if (is_block)
                for (size_t q = 0; q < block.keys.size(); ++q)
                    move_param_block_<Adapt>(chain, block);
            else if (Adapt)
                observe_param_block_(chain, block);
        settings.count_step();
    }

    //////  Parameter moves  //////

    template <bool Adapt>
    void MCMC::move_single_parameter_(Workspace &chain, size_t p)
    {
        Variable &this_par = chain.param.at(p);
        const ParamRoute &route = chain.routes[p];
        // set new value of the parameter and store the old to buffer
        this_par.propose(settings.prng(chain.id));

        // a null move is accepted without evaluation
        if (this_par.is_null_move())
        {
            if constexpr (Adapt)
                this_par.adapt_proposal(1., chain.gain);
            return;
        }

        // proceed only if the new value is consistent, else refuse the move
        if (this_par.is_valid_move())
        {
            commit_param_move_(chain, route);

            // Metropolis-Hastings acceptance log-ratio
            const double r = chain.tot_log_lik.ratio(route.terms) + this_par.get_move_log_ratio();

            if constexpr (Adapt)
                this_par.adapt_proposal(accept_probability(r), chain.gain);
            if (settings.exponential(chain.id) > -r)
                return;

            // else move is refused and log-likelihood recovered
            recover_param_move_(chain, route);
        }
        else if constexpr (Adapt)
            this_par.adapt_proposal(0., chain.gain);
        this_par.restore();
    }

    template <bool Adapt>
    void MCMC::move_param_block_(Workspace &chain, ParamBlock &block)
    {
        if constexpr (Adapt)
            block.factor.reshape();
        block.factor.propose(settings.prng(chain.id));

        // a null jump is accepted without evaluation
        if (block.factor.is_null_jump())
        {
            if constexpr (Adapt)
            {
                block.factor.adapt(1., chain.block_gain);
                observe_param_block_(chain, block);
            }
            return;
        }

        // move every member by its jump, the log-ratio starting from their Jacobians
        double r = 0.;
        for (size_t i = 0; i < block.keys.size(); ++i)
            r -= chain.param.at(block.keys[i]).block_jump(block.factor.jump[i]);

        // a jump rounded onto a bound -> refused before the likelihood
        const bool is_evaluated = r > -inf::DOUBLE;
        bool is_accepted = false;
        if (is_evaluated)
        {
            for (size_t p : block.keys)
                r += chain.param.at(p).log_ratio();
            commit_param_move_(chain, block.route);
            r += chain.tot_log_lik.ratio(block.route.terms);
            is_accepted = settings.exponential(chain.id) > -r;
        }
        if constexpr (Adapt)
            block.factor.adapt(is_evaluated ? accept_probability(r) : 0., chain.block_gain);
        if (!is_accepted)
        {
            if (is_evaluated)
                recover_param_move_(chain, block.route);
            for (size_t p : block.keys)
                chain.param.at(p).restore();
        }
        if constexpr (Adapt)
            observe_param_block_(chain, block);
    }

    //////  Latent moves  //////

    template <bool Adapt>
    void MCMC::move_all_latent_(Workspace &chain)
    {
        if (!chain.augm.scheme().is_any_sample())
            return;

        // run per-house latent moves, and accumulate accepted log-density deltas
        LogLik ll = parallel::for_reduce<LogLik>(
            chain.data.dataset().report_house.total,
            [&](size_t H) -> LogLik
            { return move_latent_single_house_<Adapt>(chain, H); },
            settings.budget.inner);
        chain.tot_log_lik.increment(ll);
    }

    template <bool Adapt>
    LogLik MCMC::move_latent_single_house_(Workspace &chain, size_t H)
    {
        const Household &house = chain.data.dataset().at(H);

        LogLik ll;
        for (const Individual &indiv : house.infected_indivs())
        {
            // a latent block moves together, one latent at a time over the warm-up
            block_rw::Factor *factor = chain.augm.block_factor(indiv);
            const bool is_block = factor && !(Adapt && chain.is_warm);
            for (size_t a : is_block ? chain.augm.list_sample_alone(indiv) : chain.augm.list_sample(indiv))
                ll += move_latent_single_variable_<Adapt>(chain, house, indiv, a);
            if (is_block)
                ll += move_latent_block_<Adapt>(chain, house, indiv, *factor);
            else if (Adapt && factor)
                observe_latent_block_(chain, indiv, *factor);
        }
        return ll;
    }

    template <bool Adapt>
    LogLik MCMC::move_latent_single_variable_(Workspace &chain,
                                              const Household &house,
                                              const Individual &indiv,
                                              size_t a)
    {
        Variable &avar = chain.augm.at(a, indiv);

        // track model consumers, whose densities shift when `a` moves
        const bool has_consumers = chain.augm.has_model_consumers(indiv, a);
        double cons_ll_old = 0.;
        if (has_consumers)
            cons_ll_old = chain.augm.sum_guarded_latent_log_lik(indiv, a);

        // whether `a` produces transmission, or its cache
        const bool produces_transm = chain.augm.produces_transm(indiv, a);

        // propose new value and save old
        avar.propose(settings.prng(chain.id));

        // a null move is accepted without evaluation
        if (avar.is_null_move())
        {
            if constexpr (Adapt)
                avar.adapt_proposal(1., chain.gain);
            return {0., 0., 0.};
        }

        if (avar.is_valid_move())
        {
            // recompute the cached quantities consuming this member
            if (produces_transm)
                update_member_transm_cache_(chain, house, indiv);

            const double d_latent =
                has_consumers
                    ? chain.augm.commit_latent_log_lik_ratio_with_consumers(indiv, a, cons_ll_old)
                    : chain.augm.commit_latent_log_lik_ratio(indiv, a);
            const double d_obs = chain.augm.commit_obs_log_lik_ratio(indiv, a);
            const double d_transm = produces_transm ? commit_transm_log_lik_ratio_(chain, house) : 0.;

            const double r = d_obs + d_transm + d_latent - avar.proposal_log_ratio();
            if constexpr (Adapt)
                avar.adapt_proposal(accept_probability(r), chain.gain);
            if (settings.exponential(chain.id) > -r)
                return {d_obs, d_transm, d_latent};

            recover_member_move_(chain, house, indiv, {chain.augm.produces_obs(indiv, a), produces_transm, true});
        }
        else if constexpr (Adapt)
            avar.adapt_proposal(0., chain.gain);
        avar.restore();

        return {0., 0., 0.};
    }

    template <bool Adapt>
    LogLik MCMC::move_latent_block_(Workspace &chain,
                                    const Household &house,
                                    const Individual &indiv,
                                    block_rw::Factor &factor)
    {
        const std::vector<size_t> &keys = chain.augm.list_block(indiv);
        const LogLikTerms &terms = chain.augm.block_terms(indiv);
        if constexpr (Adapt)
            factor.reshape();
        factor.propose(settings.prng(chain.id));

        // a null jump is accepted without evaluation
        if (factor.is_null_jump())
        {
            if constexpr (Adapt)
            {
                factor.adapt(1., chain.block_gain);
                observe_latent_block_(chain, indiv, factor);
            }
            return {};
        }

        // move every latent by its jump, the log-ratio starting from their Jacobians
        double r = 0.;
        for (size_t i = 0; i < keys.size(); ++i)
            r -= chain.augm.at(keys[i], indiv).block_jump(factor.jump[i]);

        // a jump rounded onto a bound -> refused before the likelihood
        LogLik ll;
        const bool is_evaluated = r > -inf::DOUBLE;
        if (is_evaluated)
        {
            if (terms.transm)
                update_member_transm_cache_(chain, house, indiv);

            // member log-likelihoods recomputed whole, model consumers included
            if (terms.latent)
                ll.latent = commit_member_latent_log_lik_ratio_(chain, indiv);
            if (terms.obs)
                ll.obs = commit_member_obs_log_lik_ratio_(chain, indiv);
            if (terms.transm)
                ll.transm = commit_transm_log_lik_ratio_(chain, house);
            r += ll.value();
        }
        const bool is_accepted = is_evaluated && settings.exponential(chain.id) > -r;
        if constexpr (Adapt)
            factor.adapt(is_evaluated ? accept_probability(r) : 0., chain.block_gain);
        if (!is_accepted)
        {
            if (is_evaluated)
                recover_member_move_(chain, house, indiv, terms);
            for (size_t a : keys)
                chain.augm.at(a, indiv).restore();
            ll = {};
        }
        if constexpr (Adapt)
            observe_latent_block_(chain, indiv, factor);
        return ll;
    }

    //////////////////////////////  run_inference  //////////////////////////////

    namespace
    {
        /// Bind every engine to a stop flag until destroyed.
        struct StopFlagBinding
        {
            StopFlagBinding(const std::vector<MCMC *> &e,
                            std::atomic<bool> &flag)
                : mcmc_ptrs_(e) { bind_(&flag); }
            ~StopFlagBinding() { bind_(nullptr); }

        private:
            const std::vector<MCMC *> &mcmc_ptrs_;

            void bind_(std::atomic<bool> *flag) const
            {
                for (MCMC *engine : mcmc_ptrs_)
                    engine->settings.stop_flag = flag;
            }
        };

        /// Sample each engine as a unit, writing its seconds.
        /// * unit -> on one thread, the engine dispatches its chains and their households.
        parallel::Failure sample_engines(const std::vector<MCMC *> &engines,
                                         std::vector<double> &seconds,
                                         size_t workers,
                                         std::atomic<bool> &stop_flag)
        {
            return parallel::for_units(
                engines.size(), [&](size_t i)
                { seconds[i] = engines[i]->sample(); },
                workers, 1, stop_flag);
        }

        /// Return the sum over engines of the fraction of steps done.
        double fraction_done(const std::vector<MCMC *> &engines)
        {
            double done = 0.;
            for (const MCMC *engine : engines)
                done += engine->settings.fraction_done();
            return done;
        }
    }

    InferenceRun run_inference(const std::vector<MCMC *> &mcmc_vec,
                               size_t workers,
                               const std::function<void(double)> &monitor)
    {
        InferenceRun run{std::vector<double>(mcmc_vec.size(), na::DOUBLE), {}};

        // prepare in order on the calling thread, the first error ends the run
        for (size_t i = 0; i < mcmc_vec.size(); ++i)
            try
            {
                mcmc_vec[i]->prepare_inference();
            }
            catch (const std::exception &)
            {
                run.failure = {std::current_exception(), i};
                return run;
            }

        // the first error stops every engine at its next step
        std::atomic<bool> stop_flag{false};
        const StopFlagBinding binding(mcmc_vec, stop_flag);

        run.failure = parallel::monitored(
            [&]
            { return sample_engines(mcmc_vec, run.seconds, workers, stop_flag); },
            [&]
            { monitor(fraction_done(mcmc_vec)); },
            stop_flag);
        return run;
    }
}
