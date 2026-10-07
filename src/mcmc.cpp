/* housepi - Epidemiological Analysis of Household Outbreaks.
 *
 * Copyright © 2024-2026 Institut Pasteur (Paris), CNRS and INSERM.
 * Distributed under the terms of the GNU General Public License (GPLv3).
 * See the LICENSE file for details.
 */

#include <algorithm>
#include <atomic>
#include <chrono>
#include <cmath>
#include <exception>
#include <memory>
#include <random>
#include <stdexcept>
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
#include "data/household.h"
#include "data/individual.h"
#include "data/long.h"
#include "data/wide.h"
#include "distribution.h"
#include "function.h"
#include "infec_profile.h"
#include "mcmc.h"
#include "mcmc/base.h"
#include "model.h"
#include "parameter.h"
#include "registry.h"
#include "risks.h"
#include "mcmc/setup.h"
#include "translator.h"
#include "variable.h"

namespace housepi
{
    //////////////////////////////  registry::RegModels  //////////////////////////////

    namespace registry
    {
        RegModels &RegModels::instance()
        {
            static RegModels _inst;
            return _inst;
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
            return param_.at(param_.khold(node.value));
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

            if (parse::is_data_wide(node.type))
            {
                // wire the column's read-cell: its producer points the cell at itself
                // the reader needs no rule of its own
                storage_.init_pointer(
                    &pdata_->at(pdata_->wide().khold(node.value), *indiv_ptr));
                storage_.name(node.value, ErrType::Data);
            }
            else if (parse::is_latent(node.type))
            {
                return augm_.at(augm_.khold(node.value), *indiv_ptr);
            }
            else if (parse::is_group_risk(node.type))
            {
                return pgroup_risk_->at(pgroup_risk_->khold(node.value), *indiv_ptr);
            }
        }

        return storage_;
    }

    void Translator::helper_store_read_function_(const parse::Entry &node,
                                                 const Individual *indiv_ptr,
                                                 const Household *house_ptr) const
    {
        std::unique_ptr<Function> func = registry::create_function(
            node.func, ErrType::Function, registry::func_distr_kinds(false));
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

    /// Whether `f` is an explicit reduction, folding household members.
    bool is_reduction_(std::string_view f)
    {
        return f == "sum" || f == "prod" || f == "min" || f == "max";
    }

    bool Translator::reads_member_cell_(const parse::Entry &node) const
    {
        if (node.is_terminal())
            return (parse::is_data_wide(node.type) && !augm_.def().has_row(node.value, augm::Kind::House)) ||
                   parse::is_latent(node.type) ||
                   parse::is_group_risk(node.type);
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
        // - drops out of a fold;
        // - leaves any other call undefined.
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
    const Attribute &Translator::to_numeric(
        const std::string &expr,
        const Individual *indiv_ptr) const
    {
        return to_numeric(parse::Definition(expr).root, indiv_ptr);
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
        size_t p = param_.khold(node.value);
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
            std::vector<std::string> hhids;
            std::vector<const Individual *> ptrs;

            void store(const data::Wide &data_wide,
                       const Household &house,
                       const Individual &indiv)
            {
                ids.emplace_back(data_wide.id(indiv));
                hhids.emplace_back(data_wide.hhid(house));
                ptrs.emplace_back(&indiv);
            }
        };

        /// Resolve individual and household ids for latent draws, and return tracked pointers.
        /// @note Draw lists are assumed to hold unique ids.
        std::vector<const Individual *> helper_init_individual_sample(const Augmented &augm,
                                                                      const data::Wide &data_wide,
                                                                      const std::vector<std::string> &id_draws,
                                                                      const std::vector<std::string> &house_draws,
                                                                      rnd_gen &prng)
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
            std::vector<const Household *> eligible;

            // collect explicit selections, or eligible households for the default
            for (const Household &house : data_wide)
            {
                if (house.n_infected_status() == 0)
                    continue;

                if (is_default)
                {
                    eligible.push_back(&house);
                    continue;
                }

                if (vect::contains(house_draws, data_wide.hhid(house)))
                    for (const Individual &indiv : house.infected_indivs())
                        seen.store(data_wide, house, indiv);
                else
                    for (const Individual &indiv : house.infected_indivs())
                        if (vect::contains(id_draws, data_wide.id(indiv)))
                            seen.store(data_wide, house, indiv);
            }

            // warn on draws not found in the data
            std::vector<std::string> missed_id = vect::difference(id_draws, seen.ids);
            std::vector<std::string> missed_house = vect::difference(house_draws, seen.hhids);
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
                    seen.store(data_wide, house, indiv);
            }

            return seen.ptrs;
        }
    }

    /// Initialize the parameter draws of the posterior.
    void Posterior::init_param_draws_(const Parameter &param,
                                      unsigned long sample_size)
    {
        for (size_t p : param.scheme().sample)
        {
            param_values.add(param.khold(p), sample_size);
            param_priors.add("param_lprior_" + param.khold(p), sample_size);
        }
        // function-derived parameters: value trace only, no prior
        for (size_t p : param.scheme().func)
            param_values.add(param.khold(p), sample_size);
    }

    void Posterior::init_house_draws_(const data::Wide &data_wide,
                                      unsigned long sample_size)
    {
        house_log_lik.cols.assign(data_wide.eq_class_size(), std::vector<double>(sample_size, na::DOUBLE));

        size_t H = 0;
        for (const Household &house : data_wide)
        {
            house_log_lik.names.push_back(data_wide.hhid(house));
            house_log_lik.source.push_back(data_wide.get_eq_class(H++));
        }
    }

    /// Initialize the latent draws of the posterior.
    void Posterior::init_latent_draws_(const Augmented &augm,
                                       const data::Wide &data_wide,
                                       const Setup &setup,
                                       unsigned long sample_size)
    {
        tracked_indivs = helper_init_individual_sample(
            augm, data_wide, setup.id_draws, setup.house_draws, setup.prng());

        for (const Individual *indiv : tracked_indivs)
            for (size_t a : augm.list_draws(*indiv))
                latent_values.add(data_wide.id(*indiv) + "_" + augm.khold(a), sample_size);
    }

    void Posterior::init(MCMCBase &mcmc)
    {
        param_values = {};
        param_priors = {};
        latent_values = {};
        house_log_lik = {};
        tracked_indivs.clear();

        const Setup &setup = mcmc.setup;
        sample_size_ = setup.c_sample_size * setup.n_chains;

        init_param_draws_(mcmc.param, sample_size_);
        init_house_draws_(mcmc.data.wide(), sample_size_);
        init_latent_draws_(mcmc.augm, mcmc.data.wide(), setup, sample_size_);
    }

    void Posterior::store(MCMCBase &mcmc, unsigned long chain_id, unsigned long local)
    {
        const size_t col = chain_id * mcmc.setup.c_sample_size + local;

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
        const data::Wide &data_wide = mcmc.data.wide();
        const bool with_obs = mcmc.augm.scheme().is_any_obs_model();
        const bool with_model = mcmc.augm.scheme().is_any_latent_model();
        for (size_t R = 0; R < data_wide.eq_class_size(); ++R)
        {
            const Household &head = data_wide.eq_class_head(R);
            double value = mcmc.data.state.transm_log_lik(data_wide.eq_class_head(R)).value();
            if (with_obs)
                value += mcmc.house_obs_log_lik(head);
            if (with_model)
                value += mcmc.house_latent_log_lik(head);
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

    MCMCBase::MCMCBase()
        : augm(this), translator(param, augm) {}
    MCMCBase::~MCMCBase() = default;

    void MCMCBase::init_key_holder(
        const std::vector<std::string> &data_wide_keys,
        const std::vector<std::string> &hh_member_keys,
        const std::vector<std::string> &hh_member_vec_keys,
        const std::vector<std::string> &hh_double_vec_keys,
        const std::vector<std::string> &hh_double_mtx_keys,
        const std::vector<std::string> &data_long_keys,
        const std::vector<std::string> &param_keys,
        const std::vector<std::string> &group_risk_keys,
        const std::vector<std::string> &latent_keys)
    {
        data.wide().init_key_holder(data_wide_keys);
        data.wide().init_cache_key_holders(hh_member_keys,
                                           hh_member_vec_keys,
                                           hh_double_vec_keys,
                                           hh_double_mtx_keys);
        data.long_format().init_key_holder(data_long_keys);
        param.init_key_holder(param_keys);
        risks.init_key_holder(group_risk_keys);
        augm.init_key_holder(latent_keys);

        translator.bind_data(&data);
        translator.bind_group_risk(&risks.group());
    }

    void MCMCBase::init_equivalence_class_compression()
    {
        // unique equivalence classes when sampling latents, with latent models or member observation nodes
        bool force_unique = augm.scheme().is_unique_eq_class();
        // false may also mean definitions not provided yet: if required, force unique classes
        if (!force_unique && augm.scheme().list_latent_export.size() == 0)
            force_unique = augm.khold().size() > 0;

        data.wide().init_eq_classes(force_unique, data.state);
        data.state.allocate_transm_log_lik(data.wide());
    }

    namespace infprof
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

    void MCMCBase::assert_valid_infectivity_profile() const
    {
        if (!infec_profile.generic && !infec_profile.is_latent())
            return;

        bool is_discr = infec_profile.is_required_discrete();

        if (infec_profile.generic)
            infprof::assert_zero_at_origin(*infec_profile.generic, is_discr);
        else
            for (const Household &house : data.wide())
                for (const Individual &indiv : house.infected_indivs())
                {
                    // check status-specific profiles
                    const Variable &avar = augm.at(infec_profile.latent_key, indiv);
                    // @todo Check that `avar` holds a function.
                    infprof::assert_zero_at_origin(*avar.function(), is_discr);
                }
    }

    //////  MCMCBase : Monte Carlo  //////

    void MCMCBase::init_monte_carlo_()
    {
        // check there is something to sample
        if (!param.is_any_sample() && !augm.scheme().is_any_sample())
            msg::stop(
                "Nothing to sample.\nCheck your parameter and latent definitions.");

        // whether adaptation is required
        bool is_adaptive = !param.scheme().adaptive.empty() || augm.scheme().is_any_adaptive();
        // initialize setup for mcmc
        setup.init_mcmc_configuration(is_adaptive);
    }

    void MCMCBase::init_monte_carlo_chain_()
    {
        setup.prog_set_text();       // set progress bar text
        param.restart(setup.prng()); // reset parameters
        augm.restart(setup.prng());  // reset latent

        // re-initialize the adaptable proposal of every variable
        for (size_t p : param.scheme().adaptive)
            param.at(p).init_proposal();
        if (augm.scheme().is_any_adaptive())
            for (const Household &house : data.wide())
                for (const Individual &indiv : house.infected_indivs())
                    for (size_t a : augm.list_adaptive(indiv))
                        augm.at(a, indiv).init_proposal();
    }

    void MCMCBase::recover_total_log_lik_(bool with_cache,
                                          bool with_latent_ll,
                                          bool with_obs_ll)
    {
        tot_transm_log_lik.recover();
        if (with_obs_ll)
            tot_obs_log_lik.recover();
        if (with_latent_ll)
            tot_latent_log_lik.recover();

        // recover old log-likelihood values for households
        for (size_t R = 0; R < data.wide().eq_class_size(); ++R)
        {
            const Household &house = data.wide().eq_class_head(R);
            recover_transm_log_lik_(house);

            // recover from commit_total_log_posterior_
            if (with_latent_ll)
                for (const Individual &indiv : house.infected_indivs())
                    recover_latent_log_lik_(indiv);
            if (with_obs_ll)
                for (const Individual &indiv : house.infected_indivs())
                    recover_obs_log_lik_(indiv);
        }

        if (with_cache)
        {
            data.state.recover_cache();
            infec_profile.recover();
        }
    }

    std::vector<std::string> MCMCBase::get_static_data_names() const
    {
        // exclude exported latent and observation-node columns
        // an observation node is generated by the simulation: the generated value wins over static data
        // a surrogate column is produced too, as a column of its own: it never surfaces as static data
        std::vector<std::string> produced = vect::concatenate(
            augm.khold(augm.scheme().list_latent_export), augm.khold(augm.scheme().list_obs_member));
        for (const augm::Surrogate &surr : augm.scheme().list_surrogate)
            if (surr.is_own_column())
                produced.push_back(augm.def()[surr.row].name);
        return vect::difference(data.wide().khold().all(), produced);
    }

    std::vector<size_t> MCMCBase::get_row_offsets() const
    {
        size_t n_houses = data.wide().report_house.total;
        std::vector<size_t> output(n_houses + 1, 0);
        for (size_t H = 0; H < n_houses; ++H)
            output[H + 1] = output[H] + data.wide().at(H).size();
        return output;
    }

    double MCMCBase::house_obs_log_lik(const Household &house) const
    {
        double obs = 0.;
        for (const Individual &indiv : house.infected_indivs())
            obs += obs_log_lik_(indiv).value();
        return obs;
    }

    double MCMCBase::house_latent_log_lik(const Household &house) const
    {
        double model = 0.;
        for (const Individual &indiv : house.infected_indivs())
            model += latent_log_lik_(indiv).value();
        return model;
    }

    std::vector<std::string> MCMCBase::get_col_names() const
    {
        std::vector<std::string> output = {"id", "hhid"};
        vect::append(output, get_static_data_names());
        vect::append(output, {"status"});

        output.reserve(output.size() +
                       augm.scheme().list_latent_export.size() +
                       augm.scheme().list_obs_member.size() +
                       augm.scheme().list_surrogate.size());

        for (auto idx : augm.scheme().list_latent_export)
            output.push_back(augm.khold(idx));

        for (auto idx : augm.scheme().list_obs_member)
            output.push_back(augm.khold(idx));

        for (const augm::Surrogate &surr : augm.scheme().list_surrogate)
            if (surr.is_own_column())
                output.push_back(augm.def()[surr.row].name);

        return output;
    }

    //////////////////////////////  MCMC  //////////////////////////////

    MCMC::MCMC(std::string_view model_name)
        : model(registry::RegModels::instance().create(model_name))
    {
        init_key_holder(model->key_Data(),
                        model->key_HHmember(),
                        model->key_HHmembervec(),
                        model->key_HHdoublevec(),
                        model->key_HHdoublemtx(),
                        model->key_Long(),
                        model->key_Param(),
                        model->key_GroupRisk(),
                        model->key_Latent());
        model->bind_solver(this);
        model->init_model_requirements();
    }

    void MCMC::detach_data_store()
    {
        if (!data.is_shared())
            return;
        Data fresh;
        fresh.wide().init_key_holder(model->key_Data());
        fresh.wide().copy_requests(data.wide()); // model requests outlive the data
        fresh.wide().init_cache_key_holders(
            model->key_HHmember(),
            model->key_HHmembervec(),
            model->key_HHdoublevec(),
            model->key_HHdoublemtx());
        fresh.long_format().init_key_holder(model->key_Long());
        data.share(fresh);
    }

    //////  Evaluation  //////

    double MCMC::evaluate_log_lik(std::string_view hhid)
    {
        setup.n_chains = 0; // using '0 chains' as a state variable
        init_for_evaluation_();
        param.draw_sample(setup.prng());
        augm.draw_sample(setup.prng());

        // evaluation at current state
        LogLikAccumulator acc;
        if (na::is(hhid))
        {
            // evaluation across all households
            commit_total_log_posterior_(true);
            acc.obs_ll = tot_obs_log_lik.value();
            acc.transm_ll = tot_transm_log_lik.value();
            acc.latent_ll = tot_latent_log_lik.value();
        }
        else
        {
            // single household evaluation
            infec_profile.update_cache();
            acc = commit_house_log_posterior_(
                data.wide().at(hhid),
                true, augm.scheme().is_any_latent_model(), augm.scheme().is_any_obs_model());
        }

        // complete-data household log-likelihood: observation + transmission + latent model
        return acc.obs_ll + acc.transm_ll + acc.latent_ll;
    }

    void MCMC::init_for_evaluation_()
    {
        data::Wide &data_wide = this->data.wide();
        data_wide.assert_not_all_missing_status();

        // materialize surrogate data columns before anything reads the data store
        data.reset_readcells(augm.get_produced_cols());
        augm.init_for_evaluation();
        data_wide.assert_required_everyone(data.state);
        data_wide.init_eq_classes(augm.scheme().is_unique_eq_class(), data.state);
        data.state.allocate_cache(data_wide);
        data.state.allocate_transm_log_lik(data_wide);

        param.finalize();
        risks.finalize();
        augm.finalize();
        infec_profile.finalize();
        data.finalize();

        risks.assert_identifiability(param);
        assert_valid_infectivity_profile();

        model->init_model_for_evaluation();
        for (size_t R = 0; R < data_wide.eq_class_size(); ++R)
            model->init_house_for_evaluation(data_wide.eq_class_head(R));
    }

    void MCMC::commit_total_log_posterior_(bool with_cache,
                                           bool with_latent_ll,
                                           bool with_obs_ll)
    {
        data::Wide &data_wide = this->data.wide();

        if (with_cache)
        {
            data.state.save_cache();
            infec_profile.update_cache();
        }

        // with a latent or observation model, equivalence classes coincide with households
        const bool do_latent = with_latent_ll && augm.scheme().is_any_latent_model();
        const bool do_obs = with_obs_ll && augm.scheme().is_any_obs_model();
        LogLikAccumulator acc = parallel::for_reduce<LogLikAccumulator>(
            data_wide.eq_class_size(),
            [&](size_t R) -> LogLikAccumulator
            {
                double w = data_wide.eq_class_weight(R);
                LogLikAccumulator a = commit_house_log_posterior_(
                    data_wide.eq_class_head(R),
                    with_cache, do_latent, do_obs);
                return {w * a.obs_ll, w * a.transm_ll, w * a.latent_ll};
            },
            setup.budget.inner);

        tot_transm_log_lik.set(acc.transm_ll);
        if (do_obs)
            tot_obs_log_lik.set(acc.obs_ll);
        if (do_latent)
            tot_latent_log_lik.set(acc.latent_ll);
    }

    LogLikAccumulator MCMC::commit_house_log_posterior_(const Household &house,
                                                        bool with_cache,
                                                        bool with_latent_ll,
                                                        bool with_obs_ll)
    {
        LogLikAccumulator acc;

        // re-compute transmission cached quantities before the commit
        if (with_cache)
            model->update_house_transmission_cache(house);

        // transmission
        acc.transm_ll = commit_transm_log_lik_(house);
        // observation
        if (with_obs_ll)
            for (const Individual &indiv : house.infected_indivs())
                acc.obs_ll += augm.commit_obs_log_lik(indiv);
        // latent
        if (with_latent_ll)
            for (const Individual &indiv : house.infected_indivs())
                acc.latent_ll += augm.commit_latent_log_lik(indiv);
        return acc;
    }

    //////  Metropolis-Hastings  //////

    InferenceRun run_inference(const std::vector<MCMC *> &mcmc_vec, size_t workers)
    {
        const size_t n = mcmc_vec.size();
        InferenceRun run{std::vector<double>(n, na::DOUBLE), {}};

        // prepare in order on the calling thread, the first error ends the run
        size_t inner = 1;
        for (size_t i = 0; i < n; ++i)
            try
            {
                mcmc_vec[i]->prepare_inference();
                inner = std::max(inner, mcmc_vec[i]->setup.budget.inner);
            }
            catch (const std::exception &)
            {
                run.failure = {std::current_exception(), i};
                return run;
            }

        // the first error stops every engine at its next iteration
        std::atomic<bool> stop_flag{false};
        for (MCMC *engine : mcmc_vec)
            engine->setup.stop_flag = &stop_flag;
        run.failure = parallel::for_units(
            n, [&](size_t i)
            { run.seconds[i] = mcmc_vec[i]->sample(); },
            workers, inner, stop_flag);
        for (MCMC *engine : mcmc_vec)
            engine->setup.stop_flag = nullptr;

        // close bars in order, on the calling thread
        for (MCMC *engine : mcmc_vec)
        {
            if (run.failure.error)
                engine->setup.prog_crash();
            else
                engine->setup.prog_close();
        }
        return run;
    }

    void MCMC::prepare_inference()
    {
        init_monte_carlo_();
        init_for_evaluation_();
        posterior.init(*this);
        setup.prog_open_for_mcmc_sample();
    }

    double MCMC::sample()
    {
        const auto start = std::chrono::steady_clock::now();
        setup.prog_set_text(); // refresh label

        // run multiple MCMC chains
        for (unsigned long chain_id = 0; chain_id < setup.n_chains; ++chain_id)
        {
            init_chain_();
            run_chain_(chain_id);
        }
        return std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
    }

    void MCMC::init_guarded_commit_total_()
    {
        data.state.save_cache();
        infec_profile.update_cache();

        // single serial pass
        LogLikAccumulator acc;
        for (size_t R = 0; R < data.wide().eq_class_size(); ++R)
        {
            double w = data.wide().eq_class_weight(R);
            LogLikAccumulator a = init_guarded_commit_house_(data.wide().eq_class_head(R));
            acc += {w * a.obs_ll, w * a.transm_ll, w * a.latent_ll};
        }

        tot_transm_log_lik.set(acc.transm_ll);
        if (augm.scheme().is_any_obs_model())
            tot_obs_log_lik.set(acc.obs_ll);
        if (augm.scheme().is_any_latent_model())
            tot_latent_log_lik.set(acc.latent_ll);
    }

    LogLikAccumulator MCMC::init_guarded_commit_house_(Household &house)
    {
        LogLikAccumulator acc;

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
                    if (!augm.scheme().is_any_sample() || attempt >= N_INIT_ATTEMPTS_)
                        augm.throw_obs_init_error(indiv);
                    augm.restart(indiv, setup.prng());
                }
                acc.obs_ll += augm.commit_obs_log_lik(indiv);
            }

        // transmission, on the latents the observations accepted
        model->update_house_transmission_cache(house);
        acc.transm_ll = commit_transm_log_lik_(house);
        if (data.state.transm_log_lik(house).is_inf())
            msg::stop(
                "Initial transmission log-likelihood returned infinity for " +
                data.wide().str(house) + ".");

        // latent
        if (augm.scheme().is_any_latent_model())
            for (const Individual &indiv : house.infected_indivs())
            {
                if (augm.list_latent_models(indiv).empty())
                    continue;

                acc.latent_ll += augm.commit_latent_log_lik(indiv);
                if (latent_log_lik_(indiv).is_inf())
                    msg::stop(
                        "Initial latent log-likelihood returned infinity for " +
                        data.wide().str(indiv) + ".");
            }

        return acc;
    }

    void MCMC::run_mc_finite_adaptation_()
    {
        for (unsigned long i = 0; i < setup.n_adapt; ++i)
        {
            move_all_<true>();

            // fraction of the adaptation budget spent, for vanishing adaptation
            const double ratio = static_cast<double>(i) / setup.n_adapt;

            for (size_t p : param.scheme().adaptive)
                param.at(p).adapt_proposal(ratio);
            if (augm.scheme().is_any_adaptive())
                for (const Household &house : data.wide())
                    for (const Individual &indiv : house.infected_indivs())
                        for (size_t a : augm.list_adaptive(indiv))
                            augm.at(a, indiv).adapt_proposal(ratio);
        }
    }

    void MCMC::run_mc_posterior_sampling_(unsigned long chain_id)
    {
        for (unsigned long i = 0; i < setup.n_burn; ++i)
            move_all_<false>();

        unsigned long local = 0;
        for (unsigned long i = setup.n_burn; i < setup.n_steps; i += setup.n_thin)
        {
            for (unsigned long t = 1; t <= setup.n_thin; ++t)
                move_all_<false>();
            posterior.store(*this, chain_id, local++);
        }
    }

    template <bool Adapt>
    void MCMC::move_all_()
    {
        // move latent variables to sample, if any
        move_all_latent_<Adapt>();

        // move standard model parameters
        for (size_t p : param.scheme().sample_standard)
            move_single_parameter_<Adapt>(p, false);
        // move parameters with cached dependencies
        for (size_t p : param.scheme().sample_cached)
            move_single_parameter_<Adapt>(p, true);

        // update progress bar
        setup.prog_update();
    }

    template <bool Adapt>
    void MCMC::move_single_parameter_(size_t p, bool with_cache)
    {
        Variable &this_par = param.at(p);
        // set new value of the parameter and store the old to buffer
        this_par.propose(setup.prng());

        // proceed only if the new value is consistent, else refuse the move
        if (this_par.is_valid_move())
        {
            const bool with_latent_ll = augm.scheme().is_any_latent_model() &&
                                        augm.scheme().param_produces_latent_log_lik(p);

            //////  Metropolis-Hastings  //////

            const bool with_obs_ll = augm.scheme().is_any_obs_model() &&
                                     augm.scheme().param_produces_obs_log_lik(p);

            commit_total_log_posterior_(with_cache, with_latent_ll, with_obs_ll);

            // Metropolis-Hastings acceptance log-ratio
            double r = tot_transm_log_lik.ratio();
            if (with_obs_ll)
                r += tot_obs_log_lik.ratio();
            if (with_latent_ll)
                r += tot_latent_log_lik.ratio();
            r += this_par.get_move_log_ratio();

            if (setup.exponential() > -r)
            {
                if constexpr (Adapt)
                    this_par.commit_proposal_move(setup.adapt_smoothing);
                return;
            }

            // else move is refused and log-likelihood recovered
            recover_total_log_lik_(with_cache, with_latent_ll, with_obs_ll);
        }
        this_par.restore();
        if constexpr (Adapt)
            this_par.refuse_proposal_move(setup.adapt_smoothing);
    }

    template <bool Adapt>
    void MCMC::move_all_latent_()
    {
        if (!augm.scheme().is_any_sample())
            return;

        // run per-house latent moves, and accumulate accepted log-density deltas
        LogLikAccumulator acc = parallel::for_reduce<LogLikAccumulator>(
            data.wide().report_house.total,
            [&](size_t H) -> LogLikAccumulator
            { return move_latent_single_house_<Adapt>(H); },
            setup.budget.inner);
        if (augm.scheme().is_any_obs_model())
            tot_obs_log_lik.increment(acc.obs_ll);
        tot_transm_log_lik.increment(acc.transm_ll);
        tot_latent_log_lik.increment(acc.latent_ll);
    }

    template <bool Adapt>
    LogLikAccumulator MCMC::move_latent_single_house_(size_t H)
    {
        const Household &house = data.wide().at(H);

        LogLikAccumulator acc;
        for (const Individual &indiv : house.infected_indivs())
            for (size_t a : augm.list_sample(indiv))
                acc += move_latent_single_variable_<Adapt>(house, indiv, a);
        return acc;
    }

    template <bool Adapt>
    LogLikAccumulator MCMC::move_latent_single_variable_(const Household &house,
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

        avar.propose(setup.prng()); // propose new value and save old
        if (avar.is_valid_move())   // proceed only if new value is valid
        {
            // recompute the cached quantities consuming this member
            if (produces_transm)
            {
                data.state.save_member_cache(house, indiv);
                model->update_member_transmission_cache(house, indiv);
            }

            //////  Metropolis-Hastings  //////

            const double d_latent = has_consumers
                                        ? augm.commit_latent_log_lik_ratio_with_consumers(indiv, a, cons_ll_old)
                                        : commit_latent_log_lik_ratio(indiv, a);

            const double d_obs = commit_obs_log_lik_ratio(indiv, a);
            const double d_transm = produces_transm ? commit_transm_log_lik_ratio_(house) : 0.;

            const double r = d_obs + d_transm + d_latent - avar.proposal_log_ratio();
            if (setup.exponential() > -r)
            {
                if constexpr (Adapt)
                    avar.commit_proposal_move(setup.adapt_smoothing);
                return {d_obs, d_transm, d_latent};
            }

            if (augm.produces_obs(indiv, a))
                recover_obs_log_lik_(indiv);
            if (produces_transm)
                recover_transm_log_lik_(house);
            recover_latent_log_lik_(indiv);

            if (produces_transm)
                data.state.recover_member_cache(house, indiv);
        }
        avar.restore();
        if constexpr (Adapt)
            avar.refuse_proposal_move(setup.adapt_smoothing);

        return {0., 0., 0.};
    }
}
