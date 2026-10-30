/* housepi - Epidemiological Analysis of Household Outbreaks.
 * Copyright 2024-2026 Institut Pasteur (Paris), France.
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#pragma once

#include <functional>
#include <memory>
#include <string>
#include <string_view>
#include <vector>

#include "tools/logproba.h"
#include "tools/nan.h"
#include "tools/parallel.h"
#include "augmented.h"
#include "data/dataset.h"
#include "data/long_dataset.h"
#include "infec_profile.h"
#include "mcmc/base.h"
#include "mcmc/settings.h"
#include "model.h"
#include "parameter.h"
#include "rates.h"

//////////////////////////////  MCMC  //////////////////////////////

namespace housepi
{
    class Individual;
    class Household;

    /// Maximum resamples per member, at initialization.
    static constexpr size_t N_INIT_ATTEMPTS = 100;

    /// Engine running a built-in or custom model: evaluation, inference and simulation.
    class MCMC : public MCMCBase
    {
    public:
        std::unique_ptr<ModelBase> model; /// Transmission model.

        //////  Constructor  //////

        /// Create the engine of built-in or custom model `model_name`.
        explicit MCMC(std::string_view model_name);

        /// Detach a shared data store before mutation (copy-on-write).
        void detach_data_store();

        //////  Evaluation  //////

        /// Initialize the model and evaluate the complete-data log-likelihood, joint over latents.
        /// * sum of observation, transmission and latent terms;
        /// * missing `hh_id` -> every household, flagged evaluated.
        double evaluate_log_lik(std::string_view hh_id);

        //////  Inference  //////

        /// Assert definitions and allocate the posterior, before sample().
        /// @note Call on the main thread: it may raise warnings.
        void prepare_inference();
        /// Run every chain: initialization, adaptation, then posterior sampling.
        /// * completion -> every household flagged evaluated, at the last state.
        /// @return Seconds spent sampling.
        /// @note Call after prepare_inference(). Off the main thread, it calls no R.
        double sample();

        //////  Getter  //////

        /// Number of parameter blocks.
        size_t n_param_blocks() const noexcept { return pblocks_.size(); }
        /// Parameters of parameter block `b`, baseline first.
        const std::vector<size_t> &param_block_keys(size_t b) const { return pblocks_[b].keys; }
        /// Factor of parameter block `b`.
        const block_rw::Factor &param_block_factor(size_t b) const { return pblocks_[b].factor; }

    private:
        /// Log-likelihood terms a parameter move recomputes.
        struct ParamRoute
        {
            LogLikTerms terms = {}; /// Terms the move shifts.
            /// Equivalence classes reading the parameter.
            /// * read by more than rates, or by every class -> empty, every class recomputed.
            std::vector<size_t> classes = {};
        };
        /// Parameters moving together, one block per baseline and the rates coupled to it.
        struct ParamBlock
        {
            std::vector<size_t> keys = {}; /// Parameters, baseline first.
            ParamRoute route = {};         /// Union of the routes of the members.
            block_rw::Factor factor = {};  /// Jump factor and running moments of the coordinates.
        };

        /// Create the engine of model `model_ptr`.
        explicit MCMC(std::unique_ptr<ModelBase> model_ptr);

        //////  Setup  //////

        /// Assert configuration and initialize Monte Carlo.
        void init_monte_carlo_();
        /// Assert definitions, then initialize for evaluation.
        void init_for_evaluation_();
        /// Set the route of every sampled parameter.
        void init_param_routes_();
        /// Set the parameter blocks, after the routes.
        void init_param_blocks_();

        //////  Log-posterior  //////

        /// Return every term the model has: the cache, and the latent and observation terms it declares.
        LogLikTerms all_terms_() const noexcept
        {
            return {augm.scheme().is_any_obs_model(), true, augm.scheme().is_any_latent_model()};
        }
        /// Compute the total log-posterior into `tot_log_lik`, storing every old value.
        /// * `terms`: the cache to update and the latent and observation terms to re-evaluate.
        void commit_total_log_posterior_(const LogLikTerms &terms);
        /// Commit the posterior contribution of `house`.
        LogLik commit_house_log_posterior_(const Household &house, const LogLikTerms &terms);
        /// Compute and store the transmission log-likelihood of `house`.
        double commit_transm_log_lik_(const Household &house)
        {
            double output = model->compute_house_transm_log_lik(house);
            data.state.transm_log_lik(house).set(output);
            return output;
        }
        /// Compute and store the transmission log-likelihood of `house`.
        /// @return Log-ratio against the stored value.
        double commit_transm_log_lik_ratio_(const Household &house)
        {
            commit_transm_log_lik_(house);
            return data.state.transm_log_lik(house).ratio();
        }
        /// Recover previous total log-likelihood.
        /// @param terms Terms the refused move recomputed.
        void recover_total_log_lik_(const LogLikTerms &terms);
        /// Recover previous transmission log-likelihood of `house`, after a rejected move.
        void recover_transm_log_lik_(const Household &house) noexcept
        {
            data.state.transm_log_lik(house).recover();
        }
        /// Recover previous observation log-likelihood of `m`, after a rejected move.
        void recover_obs_log_lik_(const Individual &m) noexcept
        {
            augm.state.of(m).obs_log_lik.recover();
        }
        /// Recover previous latent log-likelihood of `m`, after a rejected move.
        void recover_latent_log_lik_(const Individual &m) noexcept
        {
            augm.state.of(m).latent_log_lik.recover();
        }

        //////  Chain  //////

        /// Initialize one Markov chain.
        void init_monte_carlo_chain_();
        /// Compute the total log-posterior at initialization, guarded.
        /// * non-finite contribution -> recovered in place, or error.
        void init_guarded_commit_total_();
        /// Compute the guarded contribution of `house`.
        LogLik init_guarded_commit_house_(Household &house);
        /// Run chain `chain_id`: adaptation or burn-in, then sampling.
        /// * draws -> stored in the column block of chain `chain_id`.
        /// @note Call after init_chain_(). It calls no R and raises no warning.
        void run_chain_(Settings::Chain chain_id);
        /// Move every parameter and latent.
        /// * `Adapt` -> moves credited to the proposals, read by adaptation only.
        template <bool Adapt>
        void move_all_();

        //////  Parameter moves  //////

        /// Attempt a move of parameter `p`.
        /// * read by rates alone -> only the classes reading it are recomputed.
        template <bool Adapt>
        void move_single_parameter_(size_t p);
        /// Attempt a move of parameter block `block`, adapting its factor.
        /// * log-posterior terms of the union of the routes -> recomputed once.
        template <bool Adapt>
        void move_param_block_(ParamBlock &block);
        /// Add the coordinates of the members of `block` to the moments of its factor.
        void observe_param_block_(ParamBlock &block)
        {
            for (size_t i = 0; i < block.keys.size(); ++i)
                block.factor.z[i] = param.at(block.keys[i]).block_coord();
            block.factor.observe();
        }
        /// Commit the log-likelihood terms of a move along `route`.
        /// * classes -> the total transmission moves by their delta.
        void commit_param_move_(const ParamRoute &route)
        {
            if (route.classes.empty())
                return commit_total_log_posterior_(route.terms);

            const double delta = parallel::for_reduce<double>(
                route.classes.size(),
                [&](size_t k) -> double
                {
                    const Household &house = data.eq_class_rep(route.classes[k]);
                    return commit_transm_log_lik_ratio_(house) * data.eq_class_size(route.classes[k]);
                },
                settings.budget.inner);
            tot_log_lik.transm.increment(delta);
        }
        /// Recover the log-likelihood terms of a refused move along `route`.
        void recover_param_move_(const ParamRoute &route)
        {
            if (route.classes.empty())
                return recover_total_log_lik_(route.terms);

            tot_log_lik.transm.recover();
            for (size_t R : route.classes)
                recover_transm_log_lik_(data.eq_class_rep(R));
        }

        //////  Latent moves  //////

        /// Move the latents of every infected member.
        template <bool Adapt>
        void move_all_latent_();
        /// Move the latents of infected members of household `H`, one member at a time.
        /// * latent block -> moved together after the other latents, one at a time over the warm-up.
        /// @return Accepted deltas of transmission and latent-model log-likelihoods.
        template <bool Adapt>
        LogLik move_latent_single_house_(size_t H);
        /// Move latent `a` of member `indiv`.
        /// @return Accepted deltas of transmission and latent-model log-likelihoods.
        template <bool Adapt>
        LogLik move_latent_single_variable_(const Household &house,
                                            const Individual &indiv,
                                            size_t a);
        /// Attempt a move of the latent block of member `indiv`, adapting `factor`.
        /// * member latent and observation log-likelihoods -> recomputed whole;
        /// * transmission -> one commit, if a block latent produces it.
        /// @return Accepted deltas of observation, transmission and latent-model log-likelihoods.
        template <bool Adapt>
        LogLik move_latent_block_(const Household &house,
                                  const Individual &indiv,
                                  block_rw::Factor &factor);
        /// Add the coordinates of the latent block of member `indiv` to the running moments of `factor`.
        void observe_latent_block_(const Individual &indiv, block_rw::Factor &factor)
        {
            const std::vector<size_t> &keys = augm.list_block(indiv);
            for (size_t i = 0; i < keys.size(); ++i)
                factor.z[i] = augm.at(keys[i], indiv).block_coord();
            factor.observe();
        }
        /// Save the transmission cache of member `indiv`, then update it.
        void update_member_transm_cache_(const Household &house, const Individual &indiv)
        {
            data.state.save_member_cache(house, indiv);
            model->update_member_transmission_cache(house, indiv);
        }
        /// Compute and store the latent-model log-likelihood of member `indiv`.
        /// @return Log-ratio against the stored value.
        double commit_member_latent_log_lik_ratio_(const Individual &indiv)
        {
            augm.commit_latent_log_lik(indiv);
            return latent_log_lik_(indiv).ratio();
        }
        /// Compute and store the observation log-likelihood of member `indiv`.
        /// @return Log-ratio against the stored value.
        double commit_member_obs_log_lik_ratio_(const Individual &indiv)
        {
            augm.commit_obs_log_lik(indiv);
            return obs_log_lik_(indiv).ratio();
        }
        /// Recover terms `terms` of a refused move of member `indiv`.
        /// * transmission -> its member cache too.
        void recover_member_move_(const Household &house, const Individual &indiv, const LogLikTerms &terms)
        {
            if (terms.obs)
                recover_obs_log_lik_(indiv);
            if (terms.transm)
            {
                recover_transm_log_lik_(house);
                data.state.recover_member_cache(house, indiv);
            }
            if (terms.latent)
                recover_latent_log_lik_(indiv);
        }

        std::vector<ParamRoute> routes_ = {};  /// Route of each sampled parameter, by key.
        std::vector<ParamBlock> pblocks_ = {}; /// Blocks of the run, see Rates::list_param_blocks().
        std::vector<char> in_pblock_ = {};     /// Whether each parameter belongs to a block.
        bool is_warm_ = false;                 /// Whether the adaptation warm-up runs, block members moving one at a time.
        double gain_ = 0.;                     /// Robbins-Monro gain of the scalar proposals at this step.
        double block_gain_ = 0.;               /// Robbins-Monro gain of the block factors at this step.
    };

    /// Outcome of run_inference().
    struct InferenceRun
    {
        std::vector<double> seconds; /// Seconds spent sampling, NA where the engine did not complete.
        parallel::Failure failure;   /// Error that stopped the run, with the index of its engine.
    };

    /// Run inference on every engine, `workers` of them at once.
    /// Component-wise Metropolis-Hastings algorithm.
    /// Based on @cite metropolis_equation_1953 and @cite hastings_monte_1970
    ///
    /// * preparation -> in order, on the calling thread;
    /// * sampling -> one engine per unit, in an arena of the largest `inner` among engines;
    /// * sampling on a new thread, `monitor` called on the calling thread every 250 ms;
    /// * error while preparing -> nothing sampled;
    /// * error while sampling, or `monitor` throwing -> every engine stops at its next step.
    /// @param monitor Callback taking the sum over engines of the fraction of steps done.
    /// @note An exception thrown by `monitor` escapes once every engine has stopped.
    InferenceRun run_inference(const std::vector<MCMC *> &engine_list, size_t workers,
                               const std::function<void(double)> &monitor);
}
