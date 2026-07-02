/* housepi - Epidemiological Analysis of Household Outbreaks.
 *
 * Copyright © 2024-2026 Institut Pasteur (Paris), CNRS and INSERM.
 * Distributed under the terms of the GNU General Public License (GPLv3).
 * See the LICENSE file for details.
 */

#include <cmath>
#include <limits>
#include <memory>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

#include "tools/keys.h"
#include "tools/logproba.h"
#include "tools/messages.h"
#include "tools/nan.h"
#include "tools/operators.h"
#include "tools/parallel.h"
#include "tools/parser.h"
#include "tools/string.h"
#include "tools/types.h"

#include "augmented.h"
#include "contact_rate.h"
#include "data_wide.h"
#include "data_long.h"
#include "distribution.h"
#include "function.h"
#include "group_rate.h"
#include "household.h"
#include "individual.h"
#include "infec_profile.h"
#include "mcmc.h"
#include "mcmc_base.h"
#include "model.h"
#include "parameter.h"
#include "registry.h"
#include "setup.h"
#include "translator.h"

namespace housepi
{
    /////////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////// *** ModelRegistry : class ***  //////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////

    ModelRegistry &ModelRegistry::instance()
    {
        static ModelRegistry _inst;
        return _inst;
    }

    std::unique_ptr<ModelBase> ModelRegistry::create(std::string_view name) const
    {
        auto it = factory_map_.find(name);
        if (it == factory_map_.end())
            msg::stop(
                "Model registration failed: unknown model " + msg::val(std::string(name)) + ".\n");
        return it->second();
    }

    /////////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////  *** Translator : struct ***  //////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////

    const Attribute &Translator::to_attribute(
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
            return param_.at(param_.khold(node.value));
        }
        else if (!node.is_empty())
        {
            if (!member_ptr)
            {
                std::string sms = "Individual specific definitions as " + node.err_str() +
                                  " are not supported here.";
                // add tip if augmented (which does not require '@')
                if (parse::is_augmented(node.type) && node.is_terminal())
                    sms += "\nDid you forget rlang injection operator " + msg::code("!!" + node.value) + "?";
                msg::stop(sms);
            }

            if (parse::is_data_wide(node.type))
            {
                if (!pdata_wide_->khold().count(node.value) && augm_.is_obs_node(node.value))
                {
                    // if requested data has no data column, but is individual-observation
                    // wire individual-observation (see Augmented::varia)
                    size_t a = augm_.khold(node.value);
                    storage_.init_pointer(&augm_.at(a, *member_ptr));
                }
                else
                {
                    // otherwise wire data column's read-cell live, raw or observation
                    size_t d = pdata_wide_->khold(node.value);
                    storage_.init_pointer(&pdata_wide_->cell(d, *member_ptr));
                }
                storage_.name(node.value, ErrType::Data);
            }
            else if (parse::is_augmented(node.type))
            {
                return augm_.at(augm_.khold(node.value), *member_ptr);
            }
            else if (parse::is_group_rate(node.type))
            {
                return pgroup_rate_->at(pgroup_rate_->khold(node.value), *member_ptr);
            }
#ifndef NDEBUG
            else
                msg::stop(
                    "Unrecognized source type " + msg::val(node.type) +
                    " found while parsing " + node.err_str() + ".");
#endif
        }

        return storage_;
    }

    void Translator::helper_store_read_function_(
        const parse::Entry &node,
        const Individual *member_ptr) const
    {
        Attribute parent;
        parent.init_function(compose_func_registries(
            node.func, "function",
            {"custom", "transcendental", "arithmetic",
             "continuos", "discrete", "infecprofile"}));

        for (const auto &child : node.children)
            parent.init_func_attr(child.key, to_attribute(child, member_ptr));

        storage_ = parent;
    }
    const Attribute &Translator::to_attribute(
        const std::string &expr,
        const Individual *member_ptr) const
    {
        return to_attribute(parse::Definition(expr).root, member_ptr);
    }

    const Attribute &Translator::to_numeric(
        const parse::Entry &node,
        const Individual *member_ptr) const
    {
        storage_ = to_attribute(node, member_ptr);
        if (!storage_.is_empty() && !storage_.is_constant())
            msg::stop(
                "Definition " + node.err_str() +
                " could not be interpreted as a numeric constant.");
        return storage_;
    }
    const Attribute &Translator::to_numeric(
        const std::string &expr,
        const Individual *member_ptr) const
    {
        return to_numeric(parse::Definition(expr).root, member_ptr);
    }

    void Translator::make_cache_dependent(const parse::Entry &node) const
    {
        if (!parse::is_parameter(node.type))
            return;
        size_t p = param_.khold(node.value);
        if (!param_.at(p).is_constant())
            param_.make_cache_dependent(p);
    }

    /////////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////  *** Posterior : struct ***  ///////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////

    namespace
    {
        // Help initializing parameter draws,
        void helper_init_param_draws(
            Posterior *post,
            const Parameter &param,
            unsigned long sample_size)
        {
            size_t n_sample = param.list_sample.size();
            size_t n_func = param.list_func.size();
            post->param_values.reserve(n_sample + n_func);
            post->param_priors.reserve(n_sample);
            post->param_names.reserve(n_sample + n_func);
            post->param_prior_names.reserve(n_sample);
            for (size_t p : param.list_sample)
            {
                post->param_values.emplace_back(sample_size, na::DOUBLE);
                post->param_priors.emplace_back(sample_size, na::DOUBLE);
                post->param_names.push_back(param.khold(p));
                post->param_prior_names.push_back("param_lprior_" + param.khold(p));
            }
            // Function-derived parameters: value trace only (no prior).
            for (size_t p : param.list_func)
            {
                post->param_values.emplace_back(sample_size, na::DOUBLE);
                post->param_names.push_back(param.khold(p));
            }
        }

        /// Help initializing household transmission log-likelihood draws.
        void helper_init_house_draws(
            Posterior *post,
            const DataWide &data_wide,
            unsigned long sample_size)
        {
            post->house_log_lik.reserve(data_wide.eq_class_size());
            for (size_t R = 0; R < data_wide.eq_class_size(); ++R)
                post->house_log_lik.emplace_back(sample_size, na::DOUBLE);

            size_t n_houses = data_wide.report_house.total;
            post->house_names.reserve(n_houses);
            for (Household house : data_wide)
                post->house_names.push_back(house.house_id());
        }

        // Help initializing augmented draws,
        void helper_init_augm_draws(
            Posterior *post,
            const Augmented &augm,
            const std::vector<std::string> &id_draws,
            unsigned long sample_size)
        {
            post->tracked_members = augm.init_individual_sample(id_draws);

            if (post->tracked_members.empty())
                return;

            post->augm_values.resize(post->tracked_members.size());
            for (size_t i = 0; i < post->tracked_members.size(); ++i)
            {
                const Individual &member = *post->tracked_members[i];
                const auto &draws_list = augm.list_draws(member);
                auto &member_matrix = post->augm_values[i];
                member_matrix.reserve(draws_list.size());
                for (size_t a : draws_list)
                {
                    member_matrix.emplace_back(sample_size, na::DOUBLE);
                    post->augm_col_names.push_back(member.id() + "_" + augm.khold(a));
                }
            }
        }
    }

    void Posterior::init(MCMCBase &mcmc)
    {
        count_ = 0;
        is_power_post_ = false;
        tot_log_lik.clear();
        tot_augm_log_prior.clear();
        param_values.clear();
        param_priors.clear();
        house_log_lik.clear();
        augm_values.clear();
        tracked_members.clear();
        param_names.clear();
        param_prior_names.clear();
        augm_col_names.clear();
        house_names.clear();

        const Setup &setup = mcmc.setup;
        unsigned long sample_size = setup.c_sample_size * setup.n_chains;
        is_power_post_ = setup.is_power_post();

        // power posteriors only track the joint log-likelihood total for the
        // thermodynamic ladder; ordinary inference never reads it, so it is the
        // single quantity stored in that mode (no pointwise / parameter draws).
        if (is_power_post_)
        {
            tot_log_lik.reserve(1);
            tot_log_lik.emplace_back(sample_size, na::DOUBLE);
            return;
        }

        // allocate the augmented-variable prior total only if any is defined
        if (mcmc.augm.is_any_prior())
            tot_augm_log_prior.emplace_back(sample_size, na::DOUBLE);

        helper_init_param_draws(this, mcmc.param, sample_size);
        helper_init_house_draws(this, mcmc.data_wide(), sample_size);
        helper_init_augm_draws(this, mcmc.augm, setup.id_draws, sample_size);
    }

    void Posterior::store(MCMCBase &mcmc)
    {
        // power posteriors store only the joint (tempered) log-likelihood that
        // feeds the thermodynamic ladder: transmission plus observation
        if (is_power_post_)
        {
            double value = mcmc.tot_transm_log_lik.value();
            if (mcmc.augm.is_any_obs())
                value += mcmc.tot_obs_log_lik.value();
            tot_log_lik[0][count_++] = value;
            return;
        }

        // parameter values and priors
        const Parameter &param = mcmc.param;
        size_t i = 0;
        for (size_t p : param.list_sample)
        {
            const Variable &v = param.at(p);
            param_values[i][count_] = v.value();
            param_priors[i++][count_] = v.log_prior();
        }
        // function-derived parameters: store value, no prior
        for (size_t p : param.list_func)
            param_values[i++][count_] = param.at(p).value();

        // pointwise per-household log-likelihood: transmission plus, when an obs
        // node exists, the per-household observation contribution summed live from
        // per-member log_obs (kept current by the move machinery). Obs nodes force
        // unique eq-classes, so each class is one household.
        const DataWide &data_wide = mcmc.data_wide();
        const bool with_obs = mcmc.augm.is_any_obs();
        for (size_t R = 0; R < data_wide.eq_class_size(); ++R)
        {
            double value = data_wide.house_transm_log_lik(R);
            if (with_obs)
                value += mcmc.house_log_obs(data_wide.eq_class_head(R));
            house_log_lik[R][count_] = value;
        }

        // augmented variable draws
        const Augmented &augm = mcmc.augm;
        for (size_t j = 0; j < tracked_members.size(); ++j)
        {
            const Individual &member = *tracked_members[j];
            auto &member_matrix = augm_values[j];
            size_t c = 0;
            for (size_t a : augm.list_draws(member))
                member_matrix[c++][count_] = mcmc.aug(a, member);
        }

        // total augmented-variable prior, stored separately when present (feeds lprior)
        if (!tot_augm_log_prior.empty())
            tot_augm_log_prior[0][count_] = mcmc.tot_augm_log_prior.value();
        // update sample count
        ++count_;
    }

    /////////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////  *** MCMCBase : class ***  ////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////

    MCMCBase::MCMCBase()
        : pDataStore(std::make_shared<DataStore>()),
          augm(this), translator(param, augm) {}
    MCMCBase::~MCMCBase() = default;

    void MCMCBase::set_data_store(std::shared_ptr<DataStore> ds)
    {
        pDataStore = std::move(ds);
        translator.bind_data_wide(&pDataStore->data_wide);
    }

    void MCMCBase::init_key_holder(
        const std::vector<std::string> &data_wide_keys,
        const std::vector<std::string> &hh_size_keys,
        const std::vector<std::string> &hh_size_vec_keys,
        const std::vector<std::string> &hh_double_vec_keys,
        const std::vector<std::string> &hh_double_mtx_keys,
        const std::vector<std::string> &data_long_keys,
        const std::vector<std::string> &param_keys,
        const std::vector<std::string> &group_rate_keys,
        const std::vector<std::string> &augm_keys)
    {
        data_wide().init_key_holder(data_wide_keys);
        data_wide().init_cache_key_holders(hh_size_keys,
                                           hh_size_vec_keys,
                                           hh_double_vec_keys,
                                           hh_double_mtx_keys);
        data_long().init_key_holder(data_long_keys);
        param.init_key_holder(param_keys);
        group_rate.init_key_holder(group_rate_keys);
        augm.init_key_holder(augm_keys);

        translator.bind_data_wide(&data_wide());
        translator.bind_group_rate(&group_rate);
    }

    /// Wheter hosuehold members i and j are equivalent in terms of rates.
    /// Two members are in the same if group-rate attributes and
    /// every contact-rate involving one matches the corresponding one involving the
    /// other, both outgoing and incoming, across all other members.
    bool MCMCBase::is_same_rate_partition(const Household &house, size_t i, size_t j) const
    {
        if (i == j)
            return true;

        const Individual &member_i = house.member(i);
        const Individual &member_j = house.member(j);

        // check all group rates
        if (!group_rate.same_group(member_i, member_j))
            return false;

        // check all contact rates
        if (!contact_rate.same_contacts(house, i, j))
            return false;

        return true;
    }

    void MCMCBase::init_equivalence_class_compression()
    {
        // NOTE: equivalence compression is incompatible with sampling augmented data,
        // per-household augmented priors and per-member observation nodes
        bool force_unique = augm.is_any_sample() || augm.is_any_prior() || augm.is_any_obs();
        // however, false could also mean definitions not provided:
        // check if required and in such case set unique preventively
        if (!force_unique && augm.list_export.size() == 0)
            force_unique = augm.khold().size() > 0;

        data_wide().init_eq_classes(force_unique);
        data_wide().allocate_transm_log_lik();
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////  MCMCBase : monte carlo  /////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////

    void MCMCBase::init_monte_carlo_()
    {
        // check if there is actually something to sample
        if (!param.is_any_sample() && !augm.is_any_sample())
            msg::stop(
                "Nothing to sample.\nCheck your parameter and augmented definitions.");

        // whether adaptation is required
        bool is_adaptive = !param.list_adaptive.empty() || augm.is_any_adaptive();
        // initialize setup for mcmc
        setup.init_mcmc_configuration(is_adaptive);

        // set moving average smoothing factor for every Variable
        Variable::set_adapt_smoothing(1. / setup.n_adapt_window);
    }

    void MCMCBase::init_monte_carlo_chain_()
    {
        setup.prog_set_text();       // set progress bar text
        param.restart(setup.prng()); // reset parameters
        augm.restart(setup.prng());  // reset augmented

        // re-initialize the acceptance estimate for every variable
        for (size_t p : param.list_adaptive)
            param.at(p).reset_accept_rate();

        if (augm.is_any_adaptive())
            for (Household &house : data_wide())
                for (size_t i : house.infected_hhindexes())
                {
                    Individual &member = house.member(i);
                    for (size_t a : augm.list_adaptive(member))
                        get_augm_varia(a, member).reset_accept_rate();
                }
    }

    void MCMCBase::recover_total_log_density_(bool is_cache_depen, bool with_augm_prior,
                                              bool with_obs)
    {
        // recover old value of full log-like
        // observation / augmented priors were only touched by commit when the matching
        // flag was set: otherwise recovering would roll them back to a stale value.
        if (with_obs)
            tot_obs_log_lik.recover();
        tot_transm_log_lik.recover();
        if (with_augm_prior)
            tot_augm_log_prior.recover();

        // recover old log-likelihood values for households
        for (size_t R = 0; R < data_wide().eq_class_size(); ++R)
        {
            Household &house = data_wide().eq_class_head(R);
            data_wide().transm_log_lik(house).recover();

            // recover from commit_total_log_posterior_
            if (with_augm_prior)
                for (size_t i : house.infected_hhindexes())
                    augm.varia[Index::of(house.member(i))].log_prior.recover();
            if (with_obs)
                for (size_t i : house.infected_hhindexes())
                    augm.varia[Index::of(house.member(i))].log_obs.recover();
        }

        if (is_cache_depen)
        {
            data_wide().recover_old_cached();
            infec_profile.recover();
        }
    }

    std::vector<std::string> MCMCBase::get_static_data_names() const
    {
        // exclude exported augmented and observation-node columns: the latter are
        // generated by the simulation, so the generated value wins over static data.
        // data-bound household observations are likewise produced (exported as house
        // columns), so they must not also surface as a static data column.
        std::vector<std::string> produced =
            augm.khold(augm.list_export) + augm.khold(augm.list_obs_export);
        for (size_t b = 0; b < augm.list_house_export.size(); ++b)
            if (data_wide().khold().count(augm.house_obs_name(b)))
                produced.push_back(augm.house_obs_name(b));
        return data_wide().khold().all() - produced;
    }

    std::vector<size_t> MCMCBase::get_row_offsets() const
    {
        size_t n_houses = data_wide().report_house.total;
        std::vector<size_t> output(n_houses + 1, 0);
        for (size_t H = 0; H < n_houses; ++H)
            output[H + 1] = output[H] + data_wide().at(H).size();
        return output;
    }

    double MCMCBase::house_log_obs(const Household &house) const
    {
        double obs = 0.;
        for (size_t i : house.infected_hhindexes())
            obs += augm.varia[Index::of(house.member(i))].log_obs.value();
        return obs;
    }

    std::vector<std::string> MCMCBase::get_col_names() const
    {
        std::vector<std::string> output = {"id", "house_id"};
        output += get_static_data_names();
        output += {"status"};

        output.reserve(output.size() +
                       augm.list_export.size() +
                       augm.list_obs_export.size() +
                       augm.list_house_export.size());

        for (auto idx : augm.list_export)
            output.push_back(augm.khold(idx));

        for (auto idx : augm.list_obs_export)
            output.push_back(augm.khold(idx));

        for (size_t c = 0; c < augm.list_house_export.size(); ++c)
            output.push_back(augm.house_obs_name(c));

        return output;
    }

    /////////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////  *** MCMC : class ***  /////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////

    MCMC::MCMC(std::string_view model_name)
        : model(ModelRegistry::instance().create(model_name))
    {
        init_key_holder(model->get_DataKey(),
                        model->get_HHuintKey(),
                        model->get_HHuintvecKey(),
                        model->get_HHdoublevecKey(),
                        model->get_HHdoublemtxKey(),
                        model->get_LongKey(),
                        model->get_ParamKey(),
                        model->get_GroupRateKey(),
                        model->get_AugmKey());
        model->bind_solver(this);
        model->init_advanced_setup();
    }

    void MCMC::detach_data_store()
    {
        if (pDataStore.use_count() <= 1)
            return;
        auto new_ds = std::make_shared<DataStore>();
        new_ds->data_wide.init_key_holder(model->get_DataKey());
        new_ds->data_wide.init_cache_key_holders(
            model->get_HHuintKey(),
            model->get_HHuintvecKey(),
            model->get_HHdoublevecKey(),
            model->get_HHdoublemtxKey());
        new_ds->data_long.init_key_holder(model->get_LongKey());
        set_data_store(std::move(new_ds));
    }

    //////////////////////////
    //////  EVALUATION  //////
    //////////////////////////

    double MCMC::evaluate_log_lik(std::string_view house_id)
    {
        setup.n_chains = 0; // using '0 chains' as a state variable
        init_for_evaluation_();
        param.restart(setup.prng());
        augm.restart(setup.prng());

        // evaluation at current state
        LogProbaAccumulator acc;
        if (na::is(house_id))
        {
            // evaluation across all households
            commit_total_log_posterior_(true);
            acc.obs = tot_obs_log_lik.value();
            acc.transm = tot_transm_log_lik.value();
        }
        else
        {
            // single household evaluation
            infec_profile.update(); // eventually recompute infectivity profile
            acc = commit_house_log_posterior_(data_wide().at(house_id), true,
                                              augm.is_any_prior(), augm.is_any_obs());
        }

        // complete-data household log-likelihood (joint over augmented variables):
        // observation + transmission
        return acc.obs + acc.transm;
    }

    void MCMC::init_for_evaluation_()
    {
        data_wide().assert_not_all_missing_status();
        // materialize household observations standing in for required data columns
        // before any @D reference (augmented bounds, model) reads the data store.
        augm.materialize_data_bound_house_obs(data_wide());
        augm.init_for_evaluation();
        data_wide().init_eq_classes(
            augm.is_any_sample() || augm.is_any_prior() || augm.is_any_obs());
        data_wide().allocate_all_cache();
        data_wide().allocate_transm_log_lik();

        param.optimize();
        group_rate.optimize();
        contact_rate.optimize();
        augm.optimize();
        infec_profile.optimize();
        data_wide().optimize_cells();

        model->init_model_for_evaluation();
        for (size_t R = 0; R < data_wide().eq_class_size(); ++R)
            model->init_house_for_evaluation(data_wide().eq_class_head(R));
    }

    void MCMC::commit_total_log_posterior_(bool is_cache_depen, bool with_augm_prior, bool with_obs)
    {
        if (is_cache_depen)
        {
            data_wide().save_all_cached();
            // NOTE: save infectivity profile if lookup tables are active
            infec_profile.update();
        }

        // NOTE: eq-class coincide with all households when obs/prior exist
        const bool do_prior = with_augm_prior && augm.is_any_prior();
        const bool do_obs = with_obs && augm.is_any_obs();
        LogProbaAccumulator acc = parallel::for_reduce<LogProbaAccumulator>(
            data_wide().eq_class_size(),
            [&](size_t R) -> LogProbaAccumulator
            {
                double w = data_wide().eq_class_weight(R);
                LogProbaAccumulator a = commit_house_log_posterior_(
                    data_wide().eq_class_head(R), is_cache_depen, do_prior, do_obs);
                return {w * a.obs, w * a.transm, w * a.aprior};
            });
        if (do_obs)
            tot_obs_log_lik.set(acc.obs);
        tot_transm_log_lik.set(acc.transm);
        if (do_prior)
            tot_augm_log_prior.set(acc.aprior);
    }

    LogProbaAccumulator MCMC::commit_house_log_posterior_(Household &house, bool is_cache_depen,
                                                          bool with_augm_prior, bool with_obs)
    {
        LogProbaAccumulator acc;
        if (with_obs)
            for (size_t i : house.infected_hhindexes())
                acc.obs += augm.commit_obs_log_lik(house.member(i));

        // save all household cache dependencies
        if (is_cache_depen)
            model->update_house_cache_dependecies(house);
        acc.transm = commit_transm_log_lik_(house);
        if (with_augm_prior)
            for (size_t i : house.infected_hhindexes())
                acc.aprior += augm.commit_augm_log_prior(house.member(i));
        return acc;
    }

    double MCMC::commit_total_augm_log_prior_ratio_()
    {
        double acc = parallel::for_reduce<double>(
            data_wide().eq_class_size(),
            [&](size_t R) -> double
            {
                Household &house = data_wide().eq_class_head(R);
                double aprior = 0.;
                for (size_t i : house.infected_hhindexes())
                    aprior += augm.commit_augm_log_prior(house.member(i));
                return data_wide().eq_class_weight(R) * aprior;
            });
        tot_augm_log_prior.set(acc);
        return tot_augm_log_prior.ratio();
    }

    void MCMC::recover_total_augm_log_prior_()
    {
        tot_augm_log_prior.recover();
        for (size_t R = 0; R < data_wide().eq_class_size(); ++R)
        {
            Household &house = data_wide().eq_class_head(R);
            for (size_t i : house.infected_hhindexes())
                augm.varia[Index::of(house.member(i))].log_prior.recover();
        }
    }

    ///////////////////////////////////
    //////  METROPOLIS-HASTINGS  //////
    ///////////////////////////////////

    void MCMC::run_inference()
    {
        init_monte_carlo_();
        init_for_evaluation_();
        setup.prog_open_for_mcmc_sample();
        run_monte_carlo_();
        setup.prog_close();
    }

    void MCMC::run_monte_carlo_()
    {
        setup.prog_set_text(); // refresh label
        posterior.init(*this); // initialize posterior sample

        // run multiple MCMC chains
        for (setup.chain_id = 0; setup.chain_id < setup.n_chains; ++setup.chain_id)
        {
            init_monte_carlo_chain_();
            init_guarded_commit_total_();
            if (setup.inv_temp == 0.0)
            {
                run_mc_finite_adaptation_<true>();
                run_mc_posterior_sampling_<true>();
            }
            else
            {
                run_mc_finite_adaptation_<false>();
                run_mc_posterior_sampling_<false>();
            }
        }
    }

    void MCMC::init_guarded_commit_total_()
    {
        // cache prologue
        data_wide().save_all_cached();
        infec_profile.update();

        // single serial pass
        LogProbaAccumulator acc;
        for (size_t R = 0; R < data_wide().eq_class_size(); ++R)
        {
            double w = data_wide().eq_class_weight(R);
            LogProbaAccumulator a = init_guarded_commit_house_(data_wide().eq_class_head(R));
            acc += {w * a.obs, w * a.transm, w * a.aprior};
        }

        if (augm.is_any_obs())
            tot_obs_log_lik.set(acc.obs);
        tot_transm_log_lik.set(acc.transm);
        if (augm.is_any_prior())
            tot_augm_log_prior.set(acc.aprior);
    }

    LogProbaAccumulator MCMC::init_guarded_commit_house_(Household &house)
    {
        LogProbaAccumulator acc;

        // Observation (individual):
        // resample each member variable upon failure
        if (augm.is_any_obs())
            for (size_t i : house.infected_hhindexes())
            {
                Individual &member = house.member(i);
                if (augm.list_obs(member).empty())
                    continue;

                size_t attempt = 0;
                while (!std::isfinite(augm.eval_log_obs(member)))
                {
                    ++attempt;
                    if (!augm.is_any_sample() || attempt >= N_INIT_ATTEMPTS_)
                        augm.throw_obs_init_error(member);
                    augm.restart(member, setup.prng());
                }
                acc.obs += augm.commit_obs_log_lik(member);
            }

        // Transmission:
        model->update_house_cache_dependecies(house);
        acc.transm = commit_transm_log_lik_(house);
        if (data_wide().transm_log_lik(house).is_inf())
            msg::stop(
                "Initial transmission log-likelihood returned infinity for " +
                house.err_name() + ".");

        // Augmented Prior:
        if (augm.is_any_prior())
            for (size_t i : house.infected_hhindexes())
            {
                Individual &member = house.member(i);
                acc.aprior += augm.commit_augm_log_prior(member);
                if (augm.varia[Index::of(member)].log_prior.is_inf())
                    msg::stop(
                        "Initial augmented log-prior returned infinity for " +
                        member.err_name() + ".");
            }

        return acc;
    }

    template <bool Beta0>
    void MCMC::run_mc_finite_adaptation_()
    {
        for (unsigned long i = 0; i < setup.n_adapt; ++i)
        {
            move_all_<Beta0>();

            // compute vanishing ratio
            const double ratio = static_cast<double>(i + 1) / setup.n_adapt;
            // adapt all adaptive parameters
            for (size_t p : param.list_adaptive)
                param.at(p).adapt_proposal(ratio);
            // adapt all adaptive augmented variables
            if (augm.is_any_adaptive())
                for (Household &house : data_wide())
                    for (size_t i : house.infected_hhindexes())
                    {
                        Individual &member_i = house.member(i);
                        for (size_t a : augm.list_adaptive(member_i))
                            get_augm_varia(a, member_i).adapt_proposal(ratio);
                    }
        }
    }

    template <bool Beta0>
    void MCMC::run_mc_posterior_sampling_()
    {
        for (unsigned long i = 0; i < setup.n_burn; ++i)
            move_all_<Beta0>();

        for (unsigned long i = setup.n_burn; i < setup.n_steps; i += setup.n_thin)
        {
            for (unsigned long t = 1; t <= setup.n_thin; ++t)
                move_all_<Beta0>();
            // NOTE: recomputing likelihood, not tracked at inverse temperature 0
            if constexpr (Beta0)
                commit_total_log_posterior_(true, false, augm.is_any_obs());
            posterior.store(*this);
        }
    }

    template <bool Beta0>
    void MCMC::move_all_()
    {
        // move augmented variables to sample, if any
        move_all_augm_<Beta0>();

        // move standard model parameters
        for (size_t p : param.list_sample_standard)
            move_single_parameter_<Beta0>(p, false);
        // move parameters with cached dependecies
        for (size_t p : param.list_sample_cached)
            move_single_parameter_<Beta0>(p, true);

        // update progress bar
        setup.prog_update();
    }

    template <bool Beta0>
    void MCMC::move_single_parameter_(size_t p, bool is_cache_depen)
    {
        Variable &this_par = param.at(p);
        // set new value of the parameter and store the old to buffer
        this_par.propose(setup.prng());

        // proceed only if the new value is consistent
        // otherwise automatically refuse move
        if (this_par.is_valid_move())
        {
            // re-evaluate the augmented log-prior only if the parameter affects it
            const bool with_augm_prior = augm.is_any_prior() && augm.param_affects_augm_prior(p);

            ////////////////////
            ///  METROPOLIS  ///
            ////////////////////

            if constexpr (Beta0)
            {
                // at inverse temperature 0 the likelihood totals are not tracked
                // they are refreshed from scratch before each stored draw
                double r = this_par.get_move_log_ratio();
                if (with_augm_prior)
                    r += commit_total_augm_log_prior_ratio_();
                if (setup.exponential() > -r)
                    return;
                if (with_augm_prior)
                    recover_total_augm_log_prior_();
            }
            else
            {
                // re-evaluate the observation log-like only if the parameter affects it
                const bool with_obs = augm.is_any_obs() && augm.param_affects_obs_log_lik(p);

                // if parameter is cache dependent, update cache
                commit_total_log_posterior_(is_cache_depen, with_augm_prior, with_obs);

                const double r = setup.inv_temp * (tot_transm_log_lik.ratio() +
                                                   (with_obs ? tot_obs_log_lik.ratio() : 0.)) +
                                 (with_augm_prior ? tot_augm_log_prior.ratio() : 0.) +
                                 this_par.get_move_log_ratio();
                if (setup.exponential() > -r)
                    return;
                // else move is refused and log-likelihood recovered
                recover_total_log_density_(is_cache_depen, with_augm_prior, with_obs);
            }
        }
        // recover old value of parameter
        this_par.recover();
    }

    template <bool Beta0>
    void MCMC::move_all_augm_()
    {
        if (!augm.is_any_sample())
            return;

        // run per-house augmented moves, and accumulate accepted log-density deltas
        LogProbaAccumulator acc = parallel::for_reduce<LogProbaAccumulator>(
            data_wide().report_house.total,
            [&](size_t H) -> LogProbaAccumulator
            { return move_augm_single_house_<Beta0>(H); });
        // at inverse temperature 0 the likelihood totals are not tracked
        // they are refreshed from scratch before each stored draw
        if constexpr (!Beta0)
        {
            if (augm.is_any_obs())
                tot_obs_log_lik.increment(acc.obs);
            tot_transm_log_lik.increment(acc.transm);
        }
        tot_augm_log_prior.increment(acc.aprior);
    }

    template <bool Beta0>
    LogProbaAccumulator MCMC::move_augm_single_house_(size_t H)
    {
        Household &house = data_wide().at(H);

        LogProbaAccumulator acc;
        for (size_t i : house.infected_hhindexes())
        {
            Individual &member_i = house.member(i);
            for (size_t a : augm.list_sample(member_i))
                acc += move_augm_single_variable_<Beta0>(house, member_i, a);
        }
        return acc;
    }

    template <bool Beta0>
    LogProbaAccumulator MCMC::move_augm_single_variable_(Household &house, Individual &member, size_t a)
    {
        Variable &avar = get_augm_varia(a, member);

        avar.propose(setup.prng()); // propose new value and save old
        if (avar.is_valid_move())   // proceed only if new value is valid
        {
            size_t i = member.hhindex(); // household-relative index

            // update cache dependencies to perform computations
            data_wide().save_member_cache(house, i);
            model->update_member_cache_dependecies(house, member);

            ////////////////////
            ///  METROPOLIS  ///
            ////////////////////

            const double d_prior = commit_augm_log_prior_ratio(member, a);

            // NOTE: likelihood terms are tempered, the augmented prior is not.
            if constexpr (Beta0)
            {
                // at inverse temperature 0 the likelihood totals are not tracked
                // they are refreshed from scratch before each stored draw
                if (setup.exponential() > avar.proposal_log_ratio() - d_prior)
                    return {0., 0., d_prior};
            }
            else
            {
                const double d_obs = commit_augm_log_obs_ratio(member, a);
                const double d_transm = commit_transm_log_lik_ratio_(house);

                const double r = setup.inv_temp * (d_obs + d_transm) + d_prior - avar.proposal_log_ratio();
                if (setup.exponential() > -r)
                    return {d_obs, d_transm, d_prior};

                // recover log-likelihood values
                if (augm.affects_obs(member, a))
                    augm.varia[Index::of(member)].log_obs.recover();
                data_wide().transm_log_lik(house).recover();
            }
            // recover log-prior value
            augm.varia[Index::of(member)].log_prior.recover();
            // recover cached values
            data_wide().recover_member_old_cached(house, i);
        }
        avar.recover(); // recover old value

        return {0., 0., 0.};
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////  MCMC : estimate log_evidence  ////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////

    namespace
    {
        /// Welford online's algorithm for mean and variance.
        /// Based on @cite welford_note_1962
        /// A single non-finite draw makes the rung mean and variance undefined and
        /// is reported as {-Inf, NaN}: at inverse temperature 0 the chain legitimately
        /// visits -Inf log-likelihood states, which the caller uses to fall back to a
        /// non-adaptive ladder (the adaptive tangent placement is undefined there).
        std::pair<double, double> compute_stats(const std::vector<double> &data)
        {
            size_t n = 0;
            double avg = 0., sqr_diff = 0.;
            for (double v : data)
            {
                if (!std::isfinite(v))
                {
                    n = 0;
                    break;
                }
                const double delta = v - avg;
                avg += delta / ++n;
                sqr_diff += delta * (v - avg);
            }

            if (n == 0)
                return {-std::numeric_limits<double>::infinity(),
                        std::numeric_limits<double>::quiet_NaN()};
            return {avg, n > 1 ? sqr_diff / static_cast<double>(n - 1) : 0.};
        }
    }

    PowerLadder::PowerLadder(size_t n_rungs)
    {
        size_t n = n_rungs + 1;
        inv_temp.reserve(n);
        mean_log_lik.reserve(n);
        var_log_lik.reserve(n);
        draws.reserve(n);
    }

    void PowerLadder::push(double inv_t, std::vector<double> row)
    {
        const auto [avg, var] = compute_stats(row);
        inv_temp.push_back(inv_t);
        mean_log_lik.push_back(avg);
        var_log_lik.push_back(var);
        draws.push_back(std::move(row));
    }

    void PowerLadder::insert(size_t i, double inv_t, std::vector<double> row)
    {
        const auto [avg, var] = compute_stats(row);
        inv_temp.insert(inv_temp.begin() + i, inv_t);
        mean_log_lik.insert(mean_log_lik.begin() + i, avg);
        var_log_lik.insert(var_log_lik.begin() + i, var);
        draws.insert(draws.begin() + i, std::move(row));
    }

    PowerLadder MCMC::run_power_posteriors(std::string_view scheme, unsigned long n_rungs)
    {
        // check number of rungs for the thermodynamic ladder
        if (n_rungs < 2)
            msg::stop(
                "N. of rungs was not greater than 1.");

        setup.n_rungs = n_rungs;          // NOTE: set before progress bar
        init_monte_carlo_();              // initialize mcmc sampling
        init_for_evaluation_();           // initialize for evaluation
        setup.prog_open_for_power_post(); // print progress bar
        // UPDATE: always under estimate progress-bar time if starting from itm = 0

        //////////////////////////
        ///  POWER-POSTERIORS  ///
        //////////////////////////

        PowerLadder ladder(n_rungs);
        power_post_zero_(ladder);

        // NOTE: adaptive is undefined if -Int at inverse-temperature = 0, fall back
        ladder.scheme = scheme;
        if (scheme == "adaptive" && !std::isfinite(ladder.mean_log_lik.front()))
        {
            msg::warn(
                "Adaptive scheme is undefined at inverse temperature 0 with -Inf log-likelihood.\n"
                "Falling back to the power-fractions scheme.");
            ladder.scheme = "power-fractions";
        }

        if (ladder.scheme == "adaptive")
            power_post_adaptive_(ladder);
        else
            power_post_fractions_(ladder);

        setup.prog_close(); // close progress bar

        // restore setup state variables
        setup.inv_temp = 1.0;
        setup.rung_id = 0;
        setup.n_rungs = 0;

        return ladder;
    }

    void MCMC::power_post_zero_(PowerLadder &ladder)
    {
        setup.rung_id = 0;
        setup.inv_temp = 0.0;
        run_monte_carlo_();
        ladder.push(setup.inv_temp, posterior.tot_log_lik[0]);
    }

    void MCMC::power_post_adaptive_(PowerLadder &ladder)
    {
        // inverse temperature 1
        setup.rung_id = 1;
        setup.inv_temp = 1.0;
        run_monte_carlo_();
        ladder.push(setup.inv_temp, posterior.tot_log_lik[0]);

        // symmetrized sum_i DKL between powerpost[i] and [i+1]
        // equals to sum_i Delta_inv_temp_i * (mean_log_lik[i+1] - mean_log_lik[i])
        // subdivide interval with the 'worst' contribution
        while (setup.rung_id < setup.n_rungs)
        {
            ++setup.rung_id;

            const auto &itm = ladder.inv_temp;
            const auto &avg = ladder.mean_log_lik;
            const auto &var = ladder.var_log_lik;

            size_t k = 0;
            double worst = -1.;
            for (size_t i = 0; i + 1 < itm.size(); ++i)
            {
                const double s = std::fabs((itm[i + 1] - itm[i]) * (avg[i + 1] - avg[i]));
                if (s > worst)
                {
                    worst = s;
                    k = i;
                }
            }

            // place new inverse temperature within [itm[k], itm[k + 1]]
            double new_t;
            if (!std::isfinite(avg[k]) || !std::isfinite(avg[k + 1]) ||
                !std::isfinite(var[k]) || !std::isfinite(var[k + 1]) ||
                avg[k + 1] <= avg[k])
                // midpoint if non-monotone (Monte Carlo error)
                // or non-finite statistic (e.g. itm = 0)
                // tangent formula is undefined
                new_t = 0.5 * (itm[k] + itm[k + 1]);
            else
            {
                // intersection of tangents at interval ends
                new_t = (avg[k + 1] - itm[k + 1] * var[k + 1] - avg[k] + itm[k] * var[k]) / (var[k] - var[k + 1]);
                // variance-weighted average if intersection falls outside interval
                if (new_t <= itm[k] || new_t >= itm[k + 1])
                    new_t = itm[k] + var[k + 1] * (itm[k + 1] - itm[k]) / (var[k] + var[k + 1]);
            }

            // power posterior sampling
            setup.inv_temp = new_t;
            run_monte_carlo_();
            ladder.insert(k + 1, setup.inv_temp, posterior.tot_log_lik[0]);
        }
    }

    void MCMC::power_post_fractions_(PowerLadder &ladder)
    {
        for (setup.rung_id = 1; setup.rung_id <= setup.n_rungs; ++setup.rung_id)
        {
            setup.inv_temp = std::pow(static_cast<double>(setup.rung_id) / static_cast<double>(setup.n_rungs), 5.);
            run_monte_carlo_();
            ladder.push(setup.inv_temp, posterior.tot_log_lik[0]);
        }
    }
}
