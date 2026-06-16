/* housepi - Epidemiological Analysis of Household Outbreaks.
 *
 * Copyright © 2024-2026 Institut Pasteur (Paris), CNRS and INSERM.
 * Distributed under the terms of the GNU General Public License (GPLv3).
 * See the LICENSE file for details.
 */

#include <memory>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

#include "tools/keys.h"
#include "tools/logdensity.h"
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
            throw std::invalid_argument(
                "Model registration failed: unknown model '" + std::string(name) + "'.\n");
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
                    sms += "\nDid you forget rlang injection operator '!!'?";
                throw std::invalid_argument(sms);
            }

            if (parse::is_data_wide(node.type))
            {
                size_t d = pdata_wide_->khold(node.value);
                storage_.init_numeric(pdata_wide_->data(d, *member_ptr));
                storage_.name(node.value);
                storage_.set_err_type("data");
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
                throw std::invalid_argument(
                    "Unrecognized source type '" + node.type +
                    "' found while parsing " + node.err_str() + ".");
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
            throw std::invalid_argument(
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
            post->house_transm_log_lik.reserve(data_wide.eq_class_size());
            for (size_t R = 0; R < data_wide.eq_class_size(); ++R)
                post->house_transm_log_lik.emplace_back(sample_size, na::DOUBLE);

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
            post->augm_individuals = augm.init_individual_sample(id_draws);

            if (post->augm_individuals.empty())
                return;

            post->augm_values.resize(post->augm_individuals.size());
            for (size_t i = 0; i < post->augm_individuals.size(); ++i)
            {
                const Individual &member = *post->augm_individuals[i];
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
        tot_transm_log_lik.clear();
        tot_augm_log_prior.clear();
        param_values.clear();
        param_priors.clear();
        house_transm_log_lik.clear();
        house_obs_log_lik.clear();
        augm_values.clear();
        augm_individuals.clear();
        param_names.clear();
        param_prior_names.clear();
        augm_col_names.clear();
        house_names.clear();

        const Setup &setup = mcmc.setup;
        unsigned long sample_size = setup.c_sample_size * setup.n_chains;
        is_power_post_ = setup.is_power_post();

        // always allocate the transmission log-likelihood total (internal use)
        tot_transm_log_lik.reserve(1);
        tot_transm_log_lik.emplace_back(sample_size, na::DOUBLE);

        // allocate the augmented-variable prior total only if any is defined
        if (mcmc.augm.is_any_prior())
            tot_augm_log_prior.emplace_back(sample_size, na::DOUBLE);

        /* --FUTURE-LOG-EVIDENCE--
        // UPDATE: I would like to store also power posterior sample.
        // NOTE: it will depend on chosen rungs in that case.
        if (is_power_post_)
            return;
        */

        helper_init_param_draws(this, mcmc.param, sample_size);
        helper_init_house_draws(this, mcmc.data_wide(), sample_size);
        helper_init_augm_draws(this, mcmc.augm, setup.id_draws, sample_size);

        // pointwise observation likelihood, per eq-class, only when an obs node exists.
        // obs nodes force unique eq-classes, so each class is one household.
        if (mcmc.augm.is_any_obs())
        {
            house_obs_log_lik.reserve(mcmc.data_wide().eq_class_size());
            for (size_t R = 0; R < mcmc.data_wide().eq_class_size(); ++R)
                house_obs_log_lik.emplace_back(sample_size, na::DOUBLE);
        }
    }

    void Posterior::store(MCMCBase &mcmc)
    {
        if (!is_power_post_)
        {
            // parameter values and priors
            const Parameter &param = mcmc.param;
            size_t i = 0;
            for (size_t p : param.list_sample)
            {
                const Variable &v = param.at(p);
                param_values[i][count_] = v.value();
                param_priors[i++][count_] = v.get_log_prior();
            }
            // function-derived parameters: store value, no prior
            for (size_t p : param.list_func)
                param_values[i++][count_] = param.at(p).value();

            // pointwise per-household transmission log-likelihood
            const DataWide &data_wide = mcmc.data_wide();
            for (size_t R = 0; R < data_wide.eq_class_size(); ++R)
                house_transm_log_lik[R][count_] = data_wide.house_transm_log_lik(R);

            // pointwise per-household observation log-likelihood (when an obs node exists);
            // summed live from per-member log_obs (kept current by the move machinery)
            if (!house_obs_log_lik.empty())
                for (size_t R = 0; R < data_wide.eq_class_size(); ++R)
                    house_obs_log_lik[R][count_] = mcmc.house_log_obs(data_wide.eq_class_head(R));

            // augmented variable draws
            const Augmented &augm = mcmc.augm;
            for (size_t j = 0; j < augm_individuals.size(); ++j)
            {
                const Individual &member = *augm_individuals[j];
                auto &member_matrix = augm_values[j];
                size_t c = 0;
                for (size_t a : augm.list_draws(member))
                    member_matrix[c++][count_] = mcmc.aug(a, member);
            }
        }

        // total (tempered) transmission log-likelihood — internal (standalone / log-evidence)
        tot_transm_log_lik[0][count_] = mcmc.tot_transm_log_lik.value();
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
            throw std::length_error(
                "Nothing to sample. Check your parameter and augmented definitions.");

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
        param.restart(setup.prng()); // reset parameter user definitions
        augm.restart(setup.prng());  // reset augmented user definitions

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

    void MCMCBase::assert_finite_log_lik() const
    {
        if (!tot_obs_log_lik.is_inf() && !tot_transm_log_lik.is_inf() &&
            !tot_augm_log_prior.is_inf())
            return;

        for (size_t R = 0; R < data_wide().eq_class_size(); ++R)
        {
            const Household &house = data_wide().eq_class_head(R);

            if (data_wide().transm_log_lik(house).is_inf())
                throw std::runtime_error(
                    "Initial transmission log-likelihood returned infinity for " +
                    house.err_name() + ".");

            if (augm.is_any_sample())
                for (size_t i : house.infected_hhindexes())
                    if (augm.varia[house.member(i).index_].log_prior.is_inf())
                        throw std::runtime_error(
                            "Initial augmented log-prior returned infinity for " +
                            house.member(i).err_name() + ".");

            if (augm.is_any_obs())
                for (size_t i : house.infected_hhindexes())
                    if (augm.varia[house.member(i).index_].log_obs.is_inf())
                        throw std::runtime_error(
                            "Initial observation log-likelihood returned infinity for " +
                            house.member(i).err_name() + ".");
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

            // recover from commit_total_log_density_
            if (with_augm_prior)
                for (size_t i : house.infected_hhindexes())
                    augm.varia[house.member(i).index_].log_prior.recover();
            if (with_obs)
                for (size_t i : house.infected_hhindexes())
                    augm.varia[house.member(i).index_].log_obs.recover();
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
        return data_wide().khold().all() - augm.khold(augm.list_export) -
               augm.khold(augm.list_obs_export);
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
            obs += augm.varia[house.member(i).index_].log_obs.value();
        return obs;
    }

    std::vector<std::string> MCMCBase::get_col_names() const
    {
        std::vector<std::string> static_data_names = get_static_data_names();
        const size_t n_shared = static_data_names.size();
        const size_t n_augm_cols = augm.list_export.size();
        const size_t n_obs_cols = augm.list_obs_export.size();

        std::vector<std::string> output(2 + n_shared + 1 + n_augm_cols + n_obs_cols);

        output[0] = "id";
        output[1] = "house_id";
        for (size_t i = 0; i < n_shared; ++i)
            output[2 + i] = static_data_names[i];

        output[2 + n_shared] = "status";
        for (size_t c = 0; c < n_augm_cols; ++c)
            output[3 + n_shared + c] = augm.khold(augm.list_export[c]);
        // observation-node (generated data) columns follow the augmented ones
        for (size_t c = 0; c < n_obs_cols; ++c)
            output[3 + n_shared + n_augm_cols + c] = augm.khold(augm.list_obs_export[c]);

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

        LogDensityAccumulator acc;
        if (na::is(house_id))
        {
            // full evaluation: observation and transmission over all households
            commit_total_log_density_(true);
            acc.obs = tot_obs_log_lik.value();
            acc.transm = tot_transm_log_lik.value();
        }
        else
        {
            // single household: its observation and transmission contribution
            infec_profile.update(); // eventually recompute infectivity profile
            acc = commit_house_log_density_(data_wide().at(house_id), true,
                                            augm.is_any_prior(), augm.is_any_obs());
        }

        // complete-data household log-likelihood (joint over augmented variables):
        // observation + transmission. Priors are excluded.
        return acc.obs + acc.transm;
    }

    void MCMC::init_for_evaluation_()
    {
        data_wide().assert_not_all_missing_status();
        augm.init_for_evaluation();
        data_wide().init_eq_classes(
            augm.is_any_sample() || augm.is_any_prior() || augm.is_any_obs());
        data_wide().allocate_all_cache();
        data_wide().allocate_transm_log_lik();

        param.optimize();
        group_rate.optimize();
        contact_rate.optimize();
        augm.optimize();

        model->init_model_for_evaluation();
        for (size_t R = 0; R < data_wide().eq_class_size(); ++R)
            model->init_house_for_evaluation(data_wide().eq_class_head(R));
    }

    double MCMC::commit_transm_log_lik_(Household &house, bool is_cache_depen)
    {
        try
        {
            // save all household cache dependencies
            if (is_cache_depen)
                model->update_house_cache_dependecies(house);

            // compute and set household observation log-likelihood
            // NOTE: augmented-variable priors are accumulated separately
            return commit_transm_log_lik_(house);
        }
        catch (const std::exception &e)
        {
            throw std::runtime_error("While computing log-likelihood for " +
                                     house.err_name() + ":\n" + std::string(e.what()));
        }
    }

    void MCMC::commit_total_log_density_(bool is_cache_depen, bool with_augm_prior, bool with_obs)
    {
        // save infectivity profile if lookup tables are active
        if (is_cache_depen)
        {
            data_wide().save_all_cached();
            infec_profile.update();
        }

        // NOTE: classes are unique whenever augmented are sampled, carry a prior or
        // carry an observation node, so per-member sums are weight-1 over singletons
        const bool do_prior = with_augm_prior && augm.is_any_prior();
        const bool do_obs = with_obs && augm.is_any_obs();
        LogDensityAccumulator acc = parallel::for_reduce<LogDensityAccumulator>(
            data_wide().eq_class_size(),
            [&](size_t R) -> LogDensityAccumulator
            {
                double w = data_wide().eq_class_weight(R);
                LogDensityAccumulator a = commit_house_log_density_(
                    data_wide().eq_class_head(R), is_cache_depen, do_prior, do_obs);
                // NOTE: aggregate-init order matches LogDensityAccumulator: obs, transm, prior
                return {w * a.obs, w * a.transm, w * a.aprior};
            });
        if (do_obs)
            tot_obs_log_lik.set(acc.obs);
        tot_transm_log_lik.set(acc.transm);
        if (do_prior)
            tot_augm_log_prior.set(acc.aprior);
    }

    LogDensityAccumulator MCMC::commit_house_log_density_(Household &house, bool is_cache_depen,
                                                          bool with_augm_prior, bool with_obs)
    {
        LogDensityAccumulator acc;
        if (with_obs)
            for (size_t i : house.infected_hhindexes())
                acc.obs += augm.commit_log_lik_obs(house.member(i));
        acc.transm = commit_transm_log_lik_(house, is_cache_depen);
        if (with_augm_prior)
            for (size_t i : house.infected_hhindexes())
                acc.aprior += augm.commit_augm_log_prior(house.member(i));
        return acc;
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
        posterior.init(*this); // initialize posterior sample

        // run multiple MCMC chains&
        for (setup.chain_id = 0; setup.chain_id < setup.n_chains; ++setup.chain_id)
        {
            try
            {
                init_monte_carlo_chain_();
                // NOTE: init all cached with flag true
                commit_total_log_density_(true);
                assert_finite_log_lik();
            }
            catch (const std::exception &e)
            {
                throw std::runtime_error(
                    "While initializing Markov chain:\n" + std::string(e.what()));
            }

            run_mc_finite_adaptation_();
            run_mc_posterior_sampling_();
        }
    }

    void MCMC::run_mc_finite_adaptation_()
    {
        for (unsigned long i = 0; i < setup.n_adapt; ++i)
        {
            move_all_();

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

    void MCMC::run_mc_posterior_sampling_()
    {
        for (unsigned long i = 0; i < setup.n_burn; ++i)
            move_all_();

        for (unsigned long i = setup.n_burn; i < setup.n_steps; i += setup.n_thin)
        {
            for (unsigned long t = 1; t <= setup.n_thin; ++t)
                move_all_();
            posterior.store(*this);
        }
    }

    void MCMC::move_all_()
    {
        // move augmented variables to sample, if any
        move_all_augm_();

        // move standard model parameters
        for (size_t p : param.list_sample_standard)
            move_single_parameter_(p, false);
        // move parameters with cached dependecies
        for (size_t p : param.list_sample_cached)
            move_single_parameter_(p, true);

        // update progress bar
        setup.prog_update();
    }

    void MCMC::move_single_parameter_(size_t p, bool is_cache_depen)
    {
        Variable &this_par = param.at(p);
        // set new value of the parameter and store the old to buffer
        this_par.propose(setup.prng());

        // proceed only if the new value is consistent
        // otherwise automatically refuse move
        if (this_par.is_valid_move())
        {
            // re-evaluate the augmented log-prior / observation only if this
            // parameter feeds the respective term
            const bool with_augm_prior =
                augm.is_any_prior() && augm.param_affects_augm_prior(p);
            const bool with_obs =
                augm.is_any_obs() && augm.param_affects_obs_log_lik(p);

            // if parameter is cache dependent, update cache
            commit_total_log_density_(is_cache_depen, with_augm_prior, with_obs);

            // perform Metropolis-Hastings to accept or not move
            // NOTE: likelihoods (transmission and observation) are tempered; priors
            // (parameter and augmented-variable) enter untempered.
            const double r = setup.inv_T * tot_transm_log_lik.log_ratio() +
                             (with_augm_prior ? tot_augm_log_prior.log_ratio() : 0.) +
                             (with_obs ? setup.inv_T * tot_obs_log_lik.log_ratio() : 0.) +
                             this_par.get_move_log_ratio();
            if (setup.exponential() > -r)
                return;
            // else move is refused and log-likelihood recovered
            recover_total_log_density_(is_cache_depen, with_augm_prior, with_obs);
        }
        // recover old value of parameter
        this_par.recover();
    }

    void MCMC::move_all_augm_()
    {
        if (!augm.is_any_sample())
            return;

        // run per-house augmented moves, and accumulate accepted log-density deltas
        LogDensityAccumulator acc = parallel::for_reduce<LogDensityAccumulator>(
            data_wide().report_house.total,
            [&](size_t H) -> LogDensityAccumulator
            { return move_augm_single_house_(H); });
        if (augm.is_any_obs())
            tot_obs_log_lik.increment(acc.obs);
        tot_transm_log_lik.increment(acc.transm);
        tot_augm_log_prior.increment(acc.aprior);
    }

    LogDensityAccumulator MCMC::move_augm_single_house_(size_t H)
    {
        Household &house = data_wide().at(H);

        try
        {
            LogDensityAccumulator acc;
            for (size_t i : house.infected_hhindexes())
            {
                Individual &member_i = house.member(i);
                for (size_t a : augm.list_sample(member_i))
                    acc += move_augm_single_variable_(house, member_i, a);
            }
            return acc;
        }
        catch (const std::exception &e)
        {
            throw std::runtime_error("While sampling augmented variables for " +
                                     house.err_name() + ":\n" + std::string(e.what()));
        }
    }

    LogDensityAccumulator MCMC::move_augm_single_variable_(Household &house, Individual &member, size_t a)
    {
        Variable &avar = get_augm_varia(a, member);

        avar.propose(setup.prng()); // propose new value and save old

        // proceed only if the new value is consistent
        if (avar.is_valid_move())
        {
            size_t i = member.hhindex();
            // update dependencies to perform computations
            data_wide().save_member_cache(house, i);
            model->update_member_cache_dependecies(house, member);

            // Metropolis-Hastings. Moving this latent changes the transmission
            // likelihood, its own augmented prior, and the observation nodes it
            // parameterizes (obs_depen). Likelihood terms are tempered, prior is not.
            // only touch the observation density when this latent parameterizes an
            // obs node, so a rejected move does not roll log_obs back to a stale value.
            const bool with_obs = augm.is_any_obs() && augm.affects_obs(member, a);
            const double d_transm = commit_transm_log_lik_ratio_(house);
            const double d_prior = commit_augm_log_prior_ratio(member, a); 
            // FIXME: am I double checking 'affects_obs' basically?
            const double d_obs = with_obs ? commit_augm_log_obs_ratio(member, a) : 0.;
            if (setup.exponential() >
                avar.get_proposal_log_ratio() - setup.inv_T * (d_obs + d_transm) - d_prior)
                // aggregate-init order matches LogDensityAccumulator: obs, transm, prior
                return {d_obs, d_transm, d_prior};

            // else recover cached values
            augm.varia[member.index_].log_prior.recover();
            if (with_obs)
                augm.varia[member.index_].log_obs.recover();
            data_wide().transm_log_lik(house).recover();
            data_wide().recover_member_old_cached(house, i);
        }
        // else recover old value
        avar.recover();

        return {0., 0., 0.};
    }

    /* --FUTURE-LOG-EVIDENCE--
    ///////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////  MCMC : estimate log_evidence  ////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////

    namespace
    {
        /// Welford online's algorithm for mean and variance.
        /// Based on @cite welford_note_1962
        std::pair<double, double> compute_stats(const std::vector<double> &data)
        {
            size_t i = 0;
            size_t N = data.size();
            double x, delta, delta2;
            double sqr_diff = 0.;

            double avg = 0.;
            while (i < N)
            {
                x = data[i++];
                delta = x - avg;
                avg += delta / i;
                delta2 = x - avg;
                sqr_diff += delta * delta2;
            }

            return {avg, sqr_diff / static_cast<double>(N - 1.)};
        }
    }

    std::vector<std::vector<double>> MCMC::power_posteriors(std::string_view scheme, unsigned long n_rungs)
    {
        // check scheme choice
        const std::vector<std::string> options = {"adaptive", "power-fractions"};
        if (!contains(options, scheme))
            throw std::invalid_argument(
                "Invoked invalid power posterior scheme '" + std::string(scheme) + "'." +
                str::find_closest(scheme, options));

        // check number of rungs for the thermodynamic ladder
        if (n_rungs < 2)
            throw std::invalid_argument(
                "N. of rungs was not greater than 1.");

        setup.n_rungs = n_rungs;          // NOTE: set before progress bar
        double old_inv_T = setup.inv_T;   // save user setup
        init_monte_carlo_();              // initialize mcmc sampling
        init_for_evaluation_();           // initialize for evaluation
        setup.prog_open_for_power_post(); // print progress bar

        // thermodynamic ladder estimates
        std::vector<double> itm = {}; // inverse temperature
        std::vector<double> avg = {}; // mean log-likelihood
        std::vector<double> var = {}; // variance log-likelihood

        // run power posteriors
        power_post_zero_(itm, avg, var);
        if (scheme == "adaptive")
            power_post_adaptive_(itm, avg, var);
        else if (scheme == "power-fractions")
            power_post_fractions_(itm, avg, var);

        setup.prog_close();         // close progress bar
        setup.inv_T = old_inv_T; // restore user setup

        return {itm, avg, var};
    }

    void MCMC::power_post_zero_(std::vector<double> &itm, std::vector<double> &avg, std::vector<double> &var)
    {
        // UPDATE: should I use previous inference for itm = 1 if available?
        // UPDATE: always under-estimate time if starting from itm = 0

        setup.rung_id = 0;
        setup.inv_T = 0.0;
        run_monte_carlo_();
        const auto [stats_avg, stats_var] = compute_stats(posterior.tot_transm_log_lik);

        itm.push_back(setup.inv_T);
        avg.push_back(stats_avg);
        var.push_back(stats_var);
    }

    void MCMC::power_post_adaptive_(std::vector<double> &itm, std::vector<double> &avg, std::vector<double> &var)
    {
        // Inverse Temperature 1
        {
            setup.rung_id = 1;
            setup.inv_T = 1.0;
            run_monte_carlo_();
            const auto [stats_avg, stats_var] = compute_stats(posterior.tot_transm_log_lik);

            itm.push_back(setup.inv_T);
            avg.push_back(stats_avg);
            var.push_back(stats_var);
        }

        // Inverse Temperatures (0, 1)
        size_t k = 0; // index to interval lower bound
        bool is_right = false;
        while (setup.rung_id < setup.n_rungs)
        {
            ++setup.rung_id;

            // adaptive choice of inverse temperature
            // intersection between tangents
            double new_t = (avg[k + 1] - itm[k + 1] * var[k + 1] - avg[k] + itm[k] * var[k]) / (var[k] - var[k + 1]);
            // otherwise weighted average
            if (new_t <= itm[k] || new_t >= itm[k + 1])
                new_t = itm[k] + var[k + 1] * (itm[k + 1] - itm[k]) / (var[k] + var[k + 1]);

            // power posterior sampling
            setup.inv_T = new_t;
            run_monte_carlo_();
            const auto [stats_avg, stats_var] = compute_stats(posterior.tot_transm_log_lik);

            // find larger rectangular discretization error
            is_right = (itm[k + 1] - setup.inv_T) * (avg[k + 1] - stats_avg) > (new_t - itm[k]) * (stats_avg - avg[k]);

            // insert in-between for all lists
            itm.insert(itm.begin() + k + 1, setup.inv_T);
            avg.insert(avg.begin() + k + 1, stats_avg);
            var.insert(var.begin() + k + 1, stats_var);

            // if right interval is chosen, index points to inserted as lower bound
            k += (is_right);
        }
    }

    void MCMC::power_post_fractions_(std::vector<double> &itm, std::vector<double> &avg, std::vector<double> &var)
    {
        for (setup.rung_id = 1; setup.rung_id <= setup.n_rungs; ++setup.rung_id)
        {
            // inverse tempearture = (rung_id / #rungs)^5
            setup.inv_T = std::pow(static_cast<double>(setup.rung_id) / static_cast<double>(setup.n_rungs), 5.);
            run_monte_carlo_();

            const auto [stats_avg, stats_var] = compute_stats(posterior.tot_transm_log_lik);
            itm.push_back(setup.inv_T);
            avg.push_back(stats_avg);
            var.push_back(stats_var);
        }
    }
    */
}
