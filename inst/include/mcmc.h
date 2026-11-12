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
        //////  Constructor  //////

        /// Create the engine of built-in or custom model `model_name`.
        explicit MCMC(std::string_view model_name);

        /// Detach a shared data store before mutation (copy-on-write).
        void detach_data_store();

        //////  Evaluation  //////

        /// Initialize the model and evaluate the complete-data log-likelihood of every household.
        /// * sum of observation, transmission and latent terms, as stored in the posterior draws;
        /// * every household flagged evaluated.
        /// @return One value per household, in dataset order.
        std::vector<double> evaluate_log_lik();

        //////  Inference  //////

        /// Assert definitions and allocate the posterior, before sample().
        /// * one workspace per concurrent chain, those after the first rebuilt from definitions, warnings dropped.
        /// @note Call on the main thread: it may raise warnings.
        void prepare_inference();
        /// Run every chain: initialization, adaptation, then posterior sampling.
        /// * workspace `w` of `W` -> chains `w`, `w + W`, ... one after another, every workspace at once;
        /// * each workspace -> its households in an arena of `budget.inner` slots;
        /// * one workspace -> on the calling thread;
        /// * any exit -> one workspace kept as the first, back at chain 0, the others freed;
        /// * chain failing -> its workspace kept, the error rethrown;
        /// * otherwise -> workspace that started the highest chain kept;
        /// * workspace stopped by another engine or the user -> throws parallel::Cancelled;
        /// * completion -> every household flagged evaluated.
        /// @return Seconds spent sampling.
        /// @note Call through run_inference(), after prepare_inference(). Off the main thread, it calls no R.
        double sample();

        //////  Getter  //////

        /// Number of parameter blocks.
        size_t n_param_blocks() const noexcept { return chain().pblocks.size(); }
        /// Parameters of parameter block `b`, baseline first.
        const std::vector<size_t> &param_block_keys(size_t b) const { return chain().pblocks[b].keys; }
        /// Factor of parameter block `b`.
        const block_rw::Factor &param_block_factor(size_t b) const { return chain().pblocks[b].factor; }

    private:
        /// Name the model is registered under.
        std::string model_name_;

        //////  Workspaces  //////

        /// Build one workspace per concurrent chain, those after the first rebuilt from definitions.
        void build_workspaces_();
        /// Build a workspace from definitions, sharing the imported tables of the first workspace.
        /// * model copied from the first workspace, after its evaluation setup;
        /// * model requests, every slot in `Slot` order, then evaluation setup and parameter routes.
        std::unique_ptr<Workspace> rebuild_workspace_();

        //////  Setup  //////

        /// Assert configuration and initialize Monte Carlo.
        void init_monte_carlo_(Workspace &chain);
        /// Assert definitions, then initialize for evaluation.
        void init_for_evaluation_(Workspace &chain);
        /// Set the route of every sampled parameter.
        void init_param_routes_(Workspace &chain);
        /// Set the parameter blocks, after the routes.
        void init_param_blocks_(Workspace &chain);

        //////  Log-posterior  //////

        /// Return every term the model has: the cache, and the latent and observation terms it declares.
        LogLikTerms all_terms_(const Workspace &chain) const noexcept
        {
            return {chain.augm.scheme().is_any_obs_model(), true, chain.augm.scheme().is_any_latent_model()};
        }
        /// Compute the total log-posterior into `tot_log_lik`, storing every old value.
        /// * `terms`: the cache to update and the latent and observation terms to re-evaluate.
        void commit_total_log_posterior_(Workspace &chain, const LogLikTerms &terms);
        /// Commit the posterior contribution of `house`.
        LogLik commit_house_log_posterior_(Workspace &chain,
                                           const Household &house,
                                           const LogLikTerms &terms);
        /// Compute and store the transmission log-likelihood of `house`.
        double commit_transm_log_lik_(Workspace &chain, const Household &house)
        {
            double output = chain.model->compute_house_transm_log_lik(house);
            chain.data.state.transm_log_lik(house).set(output);
            return output;
        }
        /// Compute and store the transmission log-likelihood of `house`.
        /// @return Log-ratio against the stored value.
        double commit_transm_log_lik_ratio_(Workspace &chain, const Household &house)
        {
            commit_transm_log_lik_(chain, house);
            return chain.data.state.transm_log_lik(house).ratio();
        }
        /// Recover previous total log-likelihood.
        /// @param terms Terms the refused move recomputed.
        void recover_total_log_lik_(Workspace &chain, const LogLikTerms &terms);
        /// Recover previous transmission log-likelihood of `house`, after a rejected move.
        void recover_transm_log_lik_(Workspace &chain, const Household &house) noexcept
        {
            chain.data.state.transm_log_lik(house).recover();
        }
        /// Recover previous observation log-likelihood of `m`, after a rejected move.
        void recover_obs_log_lik_(Workspace &chain, const Individual &m) noexcept
        {
            chain.augm.state.of(m).obs_log_lik.recover();
        }
        /// Recover previous latent log-likelihood of `m`, after a rejected move.
        void recover_latent_log_lik_(Workspace &chain, const Individual &m) noexcept
        {
            chain.augm.state.of(m).latent_log_lik.recover();
        }

        //////  Chain  //////

        /// Initialize one Markov chain.
        void init_monte_carlo_chain_(Workspace &chain);
        /// Compute the total log-posterior at initialization, guarded.
        /// * non-finite contribution -> recovered in place, or error.
        void init_guarded_commit_total_(Workspace &chain);
        /// Compute the guarded contribution of `house`.
        LogLik init_guarded_commit_house_(Workspace &chain, Household &house);
        /// Run chain `chain.id`: adaptation/burn-in, then sampling.
        /// * draws -> stored in the column block of chain `chain.id`.
        /// @note Call after init_guarded_commit_total_(). It calls no R and raises no warning.
        void run_chain_(Workspace &chain);
        /// Move every parameter and latent.
        /// * `Adapt` -> moves credited to the proposals, read by adaptation only.
        template <bool Adapt>
        void move_all_(Workspace &chain);

        //////  Parameter moves  //////

        /// Attempt a move of parameter `p`.
        /// * read by rates alone -> only the classes reading it are recomputed.
        template <bool Adapt>
        void move_single_parameter_(Workspace &chain, size_t p);
        /// Attempt a move of parameter block `block`, adapting its factor.
        /// * log-posterior terms of the union of the routes -> recomputed once.
        template <bool Adapt>
        void move_param_block_(Workspace &chain, ParamBlock &block);
        /// Add the coordinates of the members of `block` to the moments of its factor.
        void observe_param_block_(Workspace &chain, ParamBlock &block)
        {
            for (size_t i = 0; i < block.keys.size(); ++i)
                block.factor.z[i] = chain.param.at(block.keys[i]).block_coord();
            block.factor.observe();
        }
        /// Commit the log-likelihood terms of a move along `route`.
        /// * classes -> the total transmission moves by their delta.
        void commit_param_move_(Workspace &chain, const ParamRoute &route)
        {
            if (route.classes.empty())
                return commit_total_log_posterior_(chain, route.terms);

            const double delta = parallel::for_reduce<double>(
                route.classes.size(),
                [&](size_t k) -> double
                {
                    const Household &house = chain.data.eq_class_rep(route.classes[k]);
                    return commit_transm_log_lik_ratio_(chain, house) * chain.data.eq_class_size(route.classes[k]);
                },
                settings.budget.inner);
            chain.tot_log_lik.transm.increment(delta);
        }
        /// Recover the log-likelihood terms of a refused move along `route`.
        void recover_param_move_(Workspace &chain, const ParamRoute &route)
        {
            if (route.classes.empty())
                return recover_total_log_lik_(chain, route.terms);

            chain.tot_log_lik.transm.recover();
            for (size_t R : route.classes)
                recover_transm_log_lik_(chain, chain.data.eq_class_rep(R));
        }

        //////  Latent moves  //////

        /// Move the latents of every infected member.
        template <bool Adapt>
        void move_all_latent_(Workspace &chain);
        /// Move the latents of infected members of household `H`, one member at a time.
        /// * latent block -> moved together after the other latents, one latent at a time over the warm-up;
        /// * draws -> household generator, restarted at (seed, chain, step, `H`).
        /// @return Accepted deltas of transmission and latent-model log-likelihoods.
        template <bool Adapt>
        LogLik move_latent_single_house_(Workspace &chain, size_t H);
        /// Move latent `a` of member `indiv`, drawing from `prng`.
        /// @return Accepted deltas of transmission and latent-model log-likelihoods.
        template <bool Adapt>
        LogLik move_latent_single_variable_(Workspace &chain,
                                            const Household &house,
                                            const Individual &indiv,
                                            size_t a,
                                            RndGen &prng);
        /// Attempt a move of the latent block of member `indiv`, drawing from `prng`, adapting `factor`.
        /// * member latent and observation log-likelihoods -> recomputed whole;
        /// * transmission -> one commit, if a block latent produces it.
        /// @return Accepted deltas of observation, transmission and latent-model log-likelihoods.
        template <bool Adapt>
        LogLik move_latent_block_(Workspace &chain,
                                  const Household &house,
                                  const Individual &indiv,
                                  block_rw::Factor &factor,
                                  RndGen &prng);
        /// Add the coordinates of the latent block of member `indiv` to the running moments of `factor`.
        void observe_latent_block_(Workspace &chain, const Individual &indiv, block_rw::Factor &factor)
        {
            const std::vector<size_t> &keys = chain.augm.list_block(indiv);
            for (size_t i = 0; i < keys.size(); ++i)
                factor.z[i] = chain.augm.at(keys[i], indiv).block_coord();
            factor.observe();
        }
        /// Save the transmission cache of member `indiv`, then update it.
        void update_member_transm_cache_(Workspace &chain, const Household &house, const Individual &indiv)
        {
            chain.data.state.save_member_cache(house, indiv);
            chain.model->update_member_transmission_cache(house, indiv);
        }
        /// Compute and store the latent-model log-likelihood of member `indiv`.
        /// @return Log-ratio against the stored value.
        double commit_member_latent_log_lik_ratio_(Workspace &chain, const Individual &indiv)
        {
            chain.augm.commit_latent_log_lik(indiv);
            return chain.latent_log_lik(indiv).ratio();
        }
        /// Compute and store the observation log-likelihood of member `indiv`.
        /// @return Log-ratio against the stored value.
        double commit_member_obs_log_lik_ratio_(Workspace &chain, const Individual &indiv)
        {
            chain.augm.commit_obs_log_lik(indiv);
            return chain.obs_log_lik(indiv).ratio();
        }
        /// Recover terms `terms` of a refused move of member `indiv`.
        /// * transmission -> its member cache too.
        void recover_member_move_(Workspace &chain,
                                  const Household &house,
                                  const Individual &indiv,
                                  const LogLikTerms &terms)
        {
            if (terms.obs)
                recover_obs_log_lik_(chain, indiv);
            if (terms.transm)
            {
                recover_transm_log_lik_(chain, house);
                chain.data.state.recover_member_cache(house, indiv);
            }
            if (terms.latent)
                recover_latent_log_lik_(chain, indiv);
        }
    };

    /// Outcome of run_inference().
    struct InferenceRun
    {
        /// Seconds spent sampling, NA where the engine did not complete.
        std::vector<double> seconds;
        /// Error that stopped the run, with the index of its engine.
        parallel::Failure failure;
    };

    /// Run inference on every engine, `workers` of them at once.
    /// Component-wise Metropolis-Hastings algorithm.
    /// Based on @cite metropolis_equation_1953 and @cite hastings_monte_1970
    ///
    /// * preparation -> in order, on the calling thread;
    /// * sampling -> one engine per unit, each dispatching its chains, see MCMC::sample();
    /// * sampling on a new thread, `monitor` called on the calling thread every 250 ms;
    /// * error while preparing -> nothing sampled;
    /// * error while sampling, or `monitor` throwing -> every engine stops at its next step.
    /// @param monitor Callback taking the sum over engines of the fraction of steps done.
    /// @note An exception thrown by `monitor` escapes once every engine has stopped.
    InferenceRun run_inference(const std::vector<MCMC *> &engine_list, size_t workers,
                               const std::function<void(double)> &monitor);
}
