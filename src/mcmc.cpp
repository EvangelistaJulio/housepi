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
#include "registry.h"
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

        std::unique_ptr<ModelBase> RegModels::create(std::string_view name) const
        {
            auto it = factory_map_.find(name);
            // # nocov start
            if (it == factory_map_.end())
                msg::stop(
                    "Model registration failed: unknown model " +
                    msg::val(std::string(name)) + ".\n");
            // # nocov end
            return it->second();
        }
    }

    //////////////////////////////  Translator  //////////////////////////////

    const Attribute &Translator::to_attribute(const parse::Entry &node,
                                              const Individual *indiv_ptr,
                                              const Household *house_ptr) const
    {
        storage_ = Attribute(); // reset auxiliary storage to default empty

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

    void Posterior::init(MCMCBase &mcmc)
    {
        param_values = {};
        param_priors = {};
        latent_values = {};
        house_log_lik = {};
        tracked_indivs.clear();

        const Settings &settings = mcmc.settings;
        const Settings::Step sample_size = settings.c_sample_size * settings.n_chains;

        init_param_draws_(mcmc.param, sample_size);
        init_house_draws_(mcmc.data, sample_size);
        init_latent_draws_(mcmc.augm, mcmc.data.dataset(), settings, sample_size);
    }

    void Posterior::store(MCMCBase &mcmc, Settings::Chain chain_id, Settings::Step local)
    {
        const size_t col = chain_id * mcmc.settings.c_sample_size + local;

        // parameter values and priors
        const Parameter &param = mcmc.param;
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

        // pointwise household log-likelihood: transmission + observation + latent
        // with a latent or observation model, every household is its own equivalence class
        const Data &data = mcmc.data;
        const bool with_obs = mcmc.augm.scheme().is_any_obs_model();
        const bool with_model = mcmc.augm.scheme().is_any_latent_model();
        for (size_t R = 0; R < data.n_eq_classes(); ++R)
        {
            const Household &rep = data.eq_class_rep(R);
            double value = mcmc.data.state.transm_log_lik(rep).value();
            if (with_obs)
                value += mcmc.house_obs_log_lik(rep);
            if (with_model)
                value += mcmc.house_latent_log_lik(rep);
            house_log_lik.cols[R][col] = value;
        }

        // latent variable draws, member by member
        const Augmented &augm = mcmc.augm;
        size_t c = 0;
        for (const Individual *indiv : tracked_indivs)
            for (size_t a : augm.list_draws(*indiv))
                latent_values.cols[c++][col] = mcmc.aug(a, *indiv);
    }

    //////////////////////////////  MCMCBase  //////////////////////////////

    MCMCBase::MCMCBase(const ModelBase &model)
        : data(model),
          param(model.param_keys()),
          rates(model.group_rate_keys()),
          augm(this, model.latent_keys()),
          translator(param, augm, data, rates.group()) {}
    MCMCBase::~MCMCBase() = default;

    void MCMCBase::init_equivalence_class_compression()
    {
        // singleton equivalence classes when sampling latents, with latent models or member observation nodes
        bool force_singleton = augm.scheme().is_singleton_eq_class();
        // false may also mean definitions not provided yet: if required, force singleton classes
        if (!force_singleton && augm.scheme().list_latent_export.size() == 0)
            force_singleton = augm.size() > 0;

        data.init_eq_classes(force_singleton);
        data.state.allocate_transm_log_lik();
    }

    namespace infprof
    {
        namespace
        {
            /// Throw if the profile transmits at generation time 0.
            /// Only a determinate profile is checked: an argument resolved per member leaves it undetermined.
            void assert_zero_at_origin(const Function &f, bool is_discr_requested)
            {
                // skip for hand-written infectivity profiles
                if (!is_discr_requested && registry::instance(registry::Kind::InfecProfile).contains(f.name()))
                    return;

                const Domain &domain = f.domain();
                // relocated domain, not support(): the dense rule of support() would relax a floor at zero
                const Codomain reach = domain.relocated();
                // the bounds keep every attainable draw past the origin
                if (reach.is_positive())
                    return;

                std::string sms = "Infectivity profile " + msg::val(f.name());

                if (na::is(domain.delocate(0.)))
                {
                    // the origin reads an undrawn value the bounds do not keep away from it
                    sms += " is not kept away from generation " +
                           std::string(is_discr_requested ? "day" : "time") + " 0.\n" +
                           msg::alert("Its support is " + reach.str() + ", so it can put infectivity on the " +
                                      (is_discr_requested ? "day" : "instant") + " of infection.");
                }
                else
                {
                    const double y = f.pdf(0.);
                    if (std::abs(y) < 1.e-10 || na::is(y))
                        return;

                    sms += " is non-zero at generation ";
                    if (is_discr_requested)
                        sms += "day 0 (F[0] = " + str::format(y) + ").\n" +
                               "A discrete-time profile admits no transmission on the day of infection.";
                    else
                        sms += "time 0 (f(0) = " + str::format(y) + ").\n" +
                               "A continuous-time profile admits no transmission at the instant of infection.";
                }

                msg::stop(
                    sms +
                    "\nConsider truncating the distribution with " + msg::arg("min") +
                    ", or shifting with " + msg::arg("loc") + ".");
            }
        }
    }

    void MCMCBase::assert_valid_infectivity_profile() const
    {
        if (!infec_profile.generic && !infec_profile.is_latent())
            return;

        bool is_discr = infec_profile.is_required_discrete();

        if (infec_profile.generic)
            infprof::assert_zero_at_origin(*infec_profile.generic, is_discr);
        else
            for (const Household &house : data.dataset())
                for (const Individual &indiv : house.infected_indivs())
                {
                    // check status-specific profiles
                    const Variable &avar = augm.at(infec_profile.latent_key, indiv);
                    // @todo Check that `avar` holds a function.
                    infprof::assert_zero_at_origin(*avar.function(), is_discr);
                }
    }

    //////  MCMCBase : Getter  //////

    std::vector<std::string> MCMCBase::get_static_data_names() const
    {
        // exclude exported latent and observation-node columns
        // an observation node is generated by the simulation: the generated value wins over static data
        // a surrogate column is produced too, as a column of its own: it never surfaces as static data
        std::vector<std::string> produced = vect::concatenate(
            augm.consult(augm.scheme().list_latent_export), augm.consult(augm.scheme().list_obs_member));
        for (const augm::Surrogate &surr : augm.scheme().list_surrogate)
            if (surr.is_own_column())
                produced.push_back(augm.def()[surr.row].name);
        return vect::difference(data.dataset().all(), produced);
    }

    std::vector<std::string> MCMCBase::get_col_names() const
    {
        std::vector<std::string> output = {"id", "hh_id"};
        vect::append(output, get_static_data_names());
        vect::append(output, {"status"});

        output.reserve(output.size() +
                       augm.scheme().list_latent_export.size() +
                       augm.scheme().list_obs_member.size() +
                       augm.scheme().list_surrogate.size());

        for (auto idx : augm.scheme().list_latent_export)
            output.push_back(augm.consult(idx));

        for (auto idx : augm.scheme().list_obs_member)
            output.push_back(augm.consult(idx));

        for (const augm::Surrogate &surr : augm.scheme().list_surrogate)
            if (surr.is_own_column())
                output.push_back(augm.def()[surr.row].name);

        return output;
    }

    std::vector<size_t> MCMCBase::get_row_offsets() const
    {
        size_t n_houses = data.dataset().report_house.total;
        std::vector<size_t> output(n_houses + 1, 0);
        for (size_t H = 0; H < n_houses; ++H)
            output[H + 1] = output[H] + data.dataset().at(H).size();
        return output;
    }

    //////////////////////////////  MCMC  //////////////////////////////

    MCMC::MCMC(std::string_view model_name)
        : MCMC(registry::RegModels::instance().create(model_name)) {}

    MCMC::MCMC(std::unique_ptr<ModelBase> model_ptr)
        : MCMCBase(*model_ptr),
          model(std::move(model_ptr))
    {
        model->bind_solver(this);
        model->init_model_requirements();
    }

    void MCMC::detach_data_store()
    {
        if (!data.is_shared())
            return;
        Data fresh(*model);
        fresh.dataset().copy_requests(data.dataset()); // model requests outlive the data
        data.share(fresh);
    }

    //////  Evaluation  //////

    double MCMC::evaluate_log_lik(std::string_view hh_id)
    {
        settings.n_chains = 0; // using '0 chains' as a state variable
        init_for_evaluation_();
        param.draw_sample(settings.prng());
        augm.draw_sample(settings.prng());

        // complete-data log-likelihood at current state: observation + transmission + latent model
        if (na::is(hh_id))
        {
            // evaluation across all households
            commit_total_log_posterior_(all_terms_());
            data.state.set_evaluated(true);
            return tot_log_lik.value();
        }

        // single household evaluation
        const Household &house = data.dataset().at(hh_id);
        infec_profile.update_cache();
        const double output = commit_house_log_posterior_(house, all_terms_()).value();
        data.state.set_evaluated(house);
        return output;
    }

    //////  Inference  //////

    void MCMC::prepare_inference()
    {
        init_monte_carlo_();
        init_for_evaluation_();
        init_param_routes_();
        posterior.init(*this);
        settings.n_done = 0;
    }

    double MCMC::sample()
    {
        const auto start = std::chrono::steady_clock::now();

        // run multiple MCMC chains
        for (Settings::Chain chain_id = 0; chain_id < settings.n_chains; ++chain_id)
        {
            init_monte_carlo_chain_();
            init_guarded_commit_total_();
            run_chain_(chain_id);
        }
        data.state.set_evaluated(true);
        return std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
    }

    //////  Setup  //////

    void MCMC::init_monte_carlo_()
    {
        // check there is something to sample
        if (!param.is_any_sample() && !augm.scheme().is_any_sample())
            msg::stop(
                "Nothing to sample.\nCheck your parameter and latent definitions.");

        // whether adaptation is required
        bool is_adaptive = !param.scheme().adaptive.empty() || augm.scheme().is_any_adaptive();
        // initialize settings for mcmc
        settings.init_mcmc_configuration(is_adaptive);
    }

    void MCMC::init_for_evaluation_()
    {
        data::Dataset &dataset = this->data.dataset();
        dataset.assert_not_all_missing_status();

        // materialize surrogate data columns before anything reads the data store
        data.reset_readcells(augm.get_produced_cols());
        augm.init_for_evaluation();
        dataset.assert_required_everyone(data.state);
        data.init_eq_classes(augm.scheme().is_singleton_eq_class());
        data.state.allocate_cache(dataset);
        data.state.allocate_transm_log_lik();

        param.finalize();
        rates.finalize();
        augm.finalize();
        infec_profile.finalize();
        data.finalize();

        rates.assert_identifiability(param);
        assert_valid_infectivity_profile();

        model->init_model_for_evaluation();
        for (size_t R = 0; R < data.n_eq_classes(); ++R)
            model->init_house_for_evaluation(data.eq_class_rep(R));
    }

    void MCMC::init_param_routes_()
    {
        const data::Dataset &dataset = data.dataset();
        const augm::Scheme &as = augm.scheme();
        routes_.assign(param.size(), {});

        const std::vector<std::vector<HouseIdx>> houses_of = rates.list_houses_of_param(param.size());
        for (size_t p : param.scheme().sample)
        {
            ParamRoute &route = routes_[p];
            route.terms.transm = vect::contains(param.scheme().sample_cached, p);
            route.terms.latent = as.is_any_latent_model() && as.param_produces_latent_log_lik(p);
            route.terms.obs = as.is_any_obs_model() && as.param_produces_obs_log_lik(p);

            // read by rates alone -> the classes reading it
            if (route.terms.latent || route.terms.obs || houses_of[p].empty() || !param.is_rate_only(p))
                continue;
            for (HouseIdx H : houses_of[p])
                route.classes.push_back(data.state.eq_class(dataset.at(H)));
            vect::sort_unique(route.classes);
            // every class reads it -> no gain
            if (route.classes.size() == data.n_eq_classes())
                route.classes.clear();
        }
        init_param_blocks_();
    }

    void MCMC::init_param_blocks_()
    {
        pblocks_.clear();
        in_pblock_.assign(param.size(), 0);
        for (std::vector<size_t> &keys : rates.list_param_blocks(param))
        {
            ParamBlock block;
            block.keys = std::move(keys);
            bool is_sparse = true;
            for (size_t p : block.keys)
            {
                in_pblock_[p] = 1;
                const ParamRoute &route = routes_[p];
                block.route.terms |= route.terms;
                is_sparse = is_sparse && !route.classes.empty();
                block.route.classes.insert(block.route.classes.end(), route.classes.begin(), route.classes.end());
            }
            // a member recomputing every class -> the block does
            vect::sort_unique(block.route.classes);
            if (!is_sparse || block.route.classes.size() == data.n_eq_classes())
                block.route.classes.clear();
            pblocks_.push_back(std::move(block));
        }
    }

    //////  Log-posterior  //////

    void MCMC::commit_total_log_posterior_(const LogLikTerms &terms)
    {
        if (terms.transm)
        {
            data.state.save_cache();
            infec_profile.update_cache();
        }

        // with a latent or observation model, equivalence classes coincide with households
        LogLik ll = parallel::for_reduce<LogLik>(
            data.n_eq_classes(),
            [&](size_t R) -> LogLik
            {
                return commit_house_log_posterior_(data.eq_class_rep(R), terms) * data.eq_class_size(R);
            },
            settings.budget.inner);
        tot_log_lik.set(ll, terms);
    }

    LogLik MCMC::commit_house_log_posterior_(const Household &house, const LogLikTerms &terms)
    {
        LogLik ll;

        // re-compute transmission cached quantities before the commit
        if (terms.transm)
            model->update_house_transmission_cache(house);

        // transmission
        ll.transm = commit_transm_log_lik_(house);
        // observation
        if (terms.obs)
            for (const Individual &indiv : house.infected_indivs())
                ll.obs += augm.commit_obs_log_lik(indiv);
        // latent
        if (terms.latent)
            for (const Individual &indiv : house.infected_indivs())
                ll.latent += augm.commit_latent_log_lik(indiv);
        return ll;
    }

    void MCMC::recover_total_log_lik_(const LogLikTerms &terms)
    {
        tot_log_lik.recover(terms);

        // recover old log-likelihood values for households
        for (size_t R = 0; R < data.n_eq_classes(); ++R)
        {
            const Household &house = data.eq_class_rep(R);
            recover_transm_log_lik_(house);

            // recover from commit_total_log_posterior_
            if (terms.latent)
                for (const Individual &indiv : house.infected_indivs())
                    recover_latent_log_lik_(indiv);
            if (terms.obs)
                for (const Individual &indiv : house.infected_indivs())
                    recover_obs_log_lik_(indiv);
        }

        if (terms.transm)
        {
            data.state.recover_cache();
            infec_profile.recover();
        }
    }

    //////  Chain  //////

    void MCMC::init_monte_carlo_chain_()
    {
        param.restart(settings.prng()); // reset parameters
        augm.restart(settings.prng());  // reset latent

        // re-initialize the adaptable proposal of every variable, latent block factors from default scales
        for (size_t p : param.scheme().adaptive)
            param.at(p).init_proposal();
        if (augm.scheme().is_any_sample())
            for (const Household &house : data.dataset())
                for (const Individual &indiv : house.infected_indivs())
                {
                    for (size_t a : augm.list_adaptive(indiv))
                        augm.at(a, indiv).init_proposal();
                    if (block_rw::Factor *factor = augm.block_factor(indiv))
                    {
                        std::vector<double> scales = {};
                        for (size_t a : augm.list_block(indiv))
                            scales.push_back(augm.at(a, indiv).block_scale());
                        factor->init(scales);
                    }
                }
    }

    void MCMC::init_guarded_commit_total_()
    {
        data.state.save_cache();
        infec_profile.update_cache();

        // single serial pass
        LogLik ll;
        for (size_t R = 0; R < data.n_eq_classes(); ++R)
            ll += init_guarded_commit_house_(data.eq_class_rep(R)) * data.eq_class_size(R);

        tot_log_lik.set(ll, all_terms_());
    }

    LogLik MCMC::init_guarded_commit_house_(Household &house)
    {
        LogLik ll;

        // observation
        // resample each individual variable upon failure
        if (augm.scheme().is_any_obs_model())
            for (const Individual &indiv : house.infected_indivs())
            {
                if (augm.list_obs_models(indiv).empty())
                    continue;

                size_t attempt = 0;
                while (!std::isfinite(augm.eval_log_obs(indiv)))
                {
                    ++attempt;
                    if (!augm.scheme().is_any_sample() || attempt >= N_INIT_ATTEMPTS)
                        augm.throw_obs_init_error(indiv);
                    augm.restart(indiv, settings.prng());
                }
                ll.obs += augm.commit_obs_log_lik(indiv);
            }

        // transmission, on the latents the observations accepted
        model->update_house_transmission_cache(house);
        ll.transm = commit_transm_log_lik_(house);
        if (!std::isfinite(ll.transm))
            msg::stop(
                "Initial transmission log-likelihood returned infinity for " +
                data.dataset().str(house) + ".");

        // latent
        if (augm.scheme().is_any_latent_model())
            for (const Individual &indiv : house.infected_indivs())
            {
                if (augm.list_latent_models(indiv).empty())
                    continue;

                ll.latent += augm.commit_latent_log_lik(indiv);
                if (!std::isfinite(ll.latent))
                    msg::stop(
                        "Initial latent log-likelihood returned infinity for " +
                        data.dataset().str(indiv) + ".");
            }

        return ll;
    }

    void MCMC::run_chain_(Settings::Chain chain_id)
    {
        //////  Adaptation  //////

        // latent blocks move one latent at a time over the warm-up, its second half shapes their factors
        const auto n_warm = static_cast<Settings::Step>(block_rw::WARMUP * settings.n_adapt);
        for (ParamBlock &block : pblocks_)
        {
            std::vector<double> scales = {};
            for (size_t p : block.keys)
                scales.push_back(param.at(p).block_scale());
            block.factor.init(scales);
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
            is_warm_ = i < n_warm;
            // drop the moments of the first warm-up half and of each window, where a chain may still drift
            if ((n_warm > 0 && i == n_warm / 2) || i == window_end)
            {
                for (block_rw::Factor &factor : augm.state.block_factors)
                    factor.forget();
                for (ParamBlock &block : pblocks_)
                    block.factor.forget();
                if (i == window_end)
                    next_window();
            }
            // gains falling linearly to 0 at the end of the adaptation, the block one from the end of the warm-up
            gain_ = 1. - static_cast<double>(i) / settings.n_adapt;
            block_gain_ = is_warm_ ? 0. : 1. - static_cast<double>(i - n_warm) / (settings.n_adapt - n_warm);
            move_all_<true>();
        }

        //////  Burn-in  //////

        for (Settings::Step i = 0; i < settings.n_burn; ++i)
            move_all_<false>();

        //////  Posterior sampling  //////

        for (Settings::Step i = 0; i < settings.c_sample_size; ++i)
        {
            for (Settings::Step t = 1; t <= settings.n_thin; ++t)
                move_all_<false>();
            posterior.store(*this, chain_id, i);
        }
    }

    template <bool Adapt>
    void MCMC::move_all_()
    {
        // move latent variables to sample, if any
        move_all_latent_<Adapt>();

        // move standard model parameters, then those with cached dependencies
        // block members move together after the warm-up, one move per member
        const bool is_block = !pblocks_.empty() && !(Adapt && is_warm_);
        for (size_t p : param.scheme().sample_standard)
            if (!(is_block && in_pblock_[p]))
                move_single_parameter_<Adapt>(p);
        for (size_t p : param.scheme().sample_cached)
            if (!(is_block && in_pblock_[p]))
                move_single_parameter_<Adapt>(p);
        for (ParamBlock &block : pblocks_)
            if (is_block)
                for (size_t q = 0; q < block.keys.size(); ++q)
                    move_param_block_<Adapt>(block);
            else if (Adapt)
                observe_param_block_(block);
        settings.count_step();
    }

    //////  Parameter moves  //////

    template <bool Adapt>
    void MCMC::move_single_parameter_(size_t p)
    {
        Variable &this_par = param.at(p);
        const ParamRoute &route = routes_[p];
        // set new value of the parameter and store the old to buffer
        this_par.propose(settings.prng());

        // proceed only if the new value is consistent, else refuse the move
        if (this_par.is_valid_move())
        {
            commit_param_move_(route);

            // Metropolis-Hastings acceptance log-ratio
            const double r = tot_log_lik.ratio(route.terms) + this_par.get_move_log_ratio();

            if constexpr (Adapt)
                this_par.adapt_proposal(accept_probability(r), gain_);
            if (settings.exponential() > -r)
                return;

            // else move is refused and log-likelihood recovered
            recover_param_move_(route);
        }
        else if constexpr (Adapt)
            this_par.adapt_proposal(0., gain_);
        this_par.restore();
    }

    template <bool Adapt>
    void MCMC::move_param_block_(ParamBlock &block)
    {
        if constexpr (Adapt)
            block.factor.reshape();
        block.factor.propose(settings.prng());

        // move every member by its jump, the log-ratio starting from their Jacobians
        double r = 0.;
        for (size_t i = 0; i < block.keys.size(); ++i)
            r -= param.at(block.keys[i]).block_jump(block.factor.jump[i]);

        // a jump rounded onto a bound -> refused before the likelihood
        const bool is_evaluated = r > -inf::DOUBLE;
        bool is_accepted = false;
        if (is_evaluated)
        {
            for (size_t p : block.keys)
                r += param.at(p).log_ratio();
            commit_param_move_(block.route);
            r += tot_log_lik.ratio(block.route.terms);
            is_accepted = settings.exponential() > -r;
        }
        if constexpr (Adapt)
            block.factor.adapt(is_evaluated ? accept_probability(r) : 0., block_gain_);
        if (!is_accepted)
        {
            if (is_evaluated)
                recover_param_move_(block.route);
            for (size_t p : block.keys)
                param.at(p).restore();
        }
        if constexpr (Adapt)
            observe_param_block_(block);
    }

    //////  Latent moves  //////

    template <bool Adapt>
    void MCMC::move_all_latent_()
    {
        if (!augm.scheme().is_any_sample())
            return;

        // run per-house latent moves, and accumulate accepted log-density deltas
        LogLik ll = parallel::for_reduce<LogLik>(
            data.dataset().report_house.total,
            [&](size_t H) -> LogLik
            { return move_latent_single_house_<Adapt>(H); },
            settings.budget.inner);
        tot_log_lik.increment(ll);
    }

    template <bool Adapt>
    LogLik MCMC::move_latent_single_house_(size_t H)
    {
        const Household &house = data.dataset().at(H);

        LogLik ll;
        for (const Individual &indiv : house.infected_indivs())
        {
            // a latent block moves together, one latent at a time over the warm-up
            block_rw::Factor *factor = augm.block_factor(indiv);
            const bool is_block = factor && !(Adapt && is_warm_);
            for (size_t a : is_block ? augm.list_sample_alone(indiv) : augm.list_sample(indiv))
                ll += move_latent_single_variable_<Adapt>(house, indiv, a);
            if (is_block)
                ll += move_latent_block_<Adapt>(house, indiv, *factor);
            else if (Adapt && factor)
                observe_latent_block_(indiv, *factor);
        }
        return ll;
    }

    template <bool Adapt>
    LogLik MCMC::move_latent_single_variable_(const Household &house,
                                              const Individual &indiv,
                                              size_t a)
    {
        Variable &avar = augm.at(a, indiv);

        // track model consumers, whose densities shift when `a` moves
        const bool has_consumers = augm.has_model_consumers(indiv, a);
        double cons_ll_old = 0.;
        if (has_consumers)
            cons_ll_old = augm.sum_guarded_latent_log_lik(indiv, a);

        // whether `a` produces transmission, or its cache
        const bool produces_transm = augm.produces_transm(indiv, a);

        avar.propose(settings.prng()); // propose new value and save old
        if (avar.is_valid_move())      // proceed only if new value is valid
        {
            // recompute the cached quantities consuming this member
            if (produces_transm)
                update_member_transm_cache_(house, indiv);

            const double d_latent =
                has_consumers
                    ? augm.commit_latent_log_lik_ratio_with_consumers(indiv, a, cons_ll_old)
                    : augm.commit_latent_log_lik_ratio(indiv, a);
            const double d_obs = augm.commit_obs_log_lik_ratio(indiv, a);
            const double d_transm = produces_transm ? commit_transm_log_lik_ratio_(house) : 0.;

            const double r = d_obs + d_transm + d_latent - avar.proposal_log_ratio();
            if constexpr (Adapt)
                avar.adapt_proposal(accept_probability(r), gain_);
            if (settings.exponential() > -r)
                return {d_obs, d_transm, d_latent};

            recover_member_move_(house, indiv, {augm.produces_obs(indiv, a), produces_transm, true});
        }
        else if constexpr (Adapt)
            avar.adapt_proposal(0., gain_);
        avar.restore();

        return {0., 0., 0.};
    }

    template <bool Adapt>
    LogLik MCMC::move_latent_block_(const Household &house,
                                    const Individual &indiv,
                                    block_rw::Factor &factor)
    {
        const std::vector<size_t> &keys = augm.list_block(indiv);
        const LogLikTerms &terms = augm.block_terms(indiv);
        if constexpr (Adapt)
            factor.reshape();
        factor.propose(settings.prng());

        // move every latent by its jump, the log-ratio starting from their Jacobians
        double r = 0.;
        for (size_t i = 0; i < keys.size(); ++i)
            r -= augm.at(keys[i], indiv).block_jump(factor.jump[i]);

        // a jump rounded onto a bound -> refused before the likelihood
        LogLik ll;
        const bool is_evaluated = r > -inf::DOUBLE;
        if (is_evaluated)
        {
            if (terms.transm)
                update_member_transm_cache_(house, indiv);

            // member log-likelihoods recomputed whole, model consumers included
            if (terms.latent)
                ll.latent = commit_member_latent_log_lik_ratio_(indiv);
            if (terms.obs)
                ll.obs = commit_member_obs_log_lik_ratio_(indiv);
            if (terms.transm)
                ll.transm = commit_transm_log_lik_ratio_(house);
            r += ll.value();
        }
        const bool is_accepted = is_evaluated && settings.exponential() > -r;
        if constexpr (Adapt)
            factor.adapt(is_evaluated ? accept_probability(r) : 0., block_gain_);
        if (!is_accepted)
        {
            if (is_evaluated)
                recover_member_move_(house, indiv, terms);
            for (size_t a : keys)
                augm.at(a, indiv).restore();
            ll = {};
        }
        if constexpr (Adapt)
            observe_latent_block_(indiv, factor);
        return ll;
    }

    //////////////////////////////  run_inference  //////////////////////////////

    namespace
    {
        /// Bind every engine to a stop flag until destroyed.
        struct StopFlagBinding
        {
            StopFlagBinding(const std::vector<MCMC *> &e,
                            const std::atomic<bool> &flag)
                : mcmc_ptrs_(e) { bind_(&flag); }
            ~StopFlagBinding() { bind_(nullptr); }

        private:
            const std::vector<MCMC *> &mcmc_ptrs_;

            void bind_(const std::atomic<bool> *flag) const
            {
                for (MCMC *engine : mcmc_ptrs_)
                    engine->settings.stop_flag = flag;
            }
        };

        /// Sample each engine as a unit, writing its seconds.
        parallel::Failure sample_engines(const std::vector<MCMC *> &engines,
                                         std::vector<double> &seconds,
                                         size_t workers,
                                         size_t inner,
                                         std::atomic<bool> &stop_flag)
        {
            return parallel::for_units(
                engines.size(), [&](size_t i)
                { seconds[i] = engines[i]->sample(); },
                workers, inner, stop_flag);
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
        size_t inner = 1;
        for (size_t i = 0; i < mcmc_vec.size(); ++i)
            try
            {
                mcmc_vec[i]->prepare_inference();
                inner = std::max(inner, mcmc_vec[i]->settings.budget.inner);
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
            { return sample_engines(mcmc_vec, run.seconds, workers, inner, stop_flag); },
            [&]
            { monitor(fraction_done(mcmc_vec)); },
            stop_flag);
        return run;
    }
}
