/* housepi - Epidemiological Analysis of Household Outbreaks.
 * Copyright 2024-2026 Institut Pasteur (Paris), France.
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#pragma once

#include <memory>
#include <string>
#include <vector>

#include "tools/logproba.h"
#include "tools/nan.h"
#include "augmented.h"
#include "data.h"
#include "infec_profile.h"
#include "mcmc/posterior.h"
#include "mcmc/settings.h"
#include "parameter.h"
#include "proposal.h"
#include "rates.h"
#include "status_process.h"
#include "translator.h"

//////////////////////////////  ChainState  //////////////////////////////

namespace housepi
{
    class ModelBase;

    /// Log-likelihood terms a parameter move recomputes.
    struct ParamRoute
    {
        /// Terms the move shifts.
        LogLikTerms terms = {};
        /// Equivalence classes reading the parameter.
        /// * read by more than rates, or by every class -> empty, every class recomputed.
        std::vector<size_t> classes = {};
    };

    /// Parameters moving together, one block per baseline and the rates coupled to it.
    struct ParamBlock
    {
        /// Parameters, baseline first.
        std::vector<size_t> keys = {};
        /// Union of the routes of the members.
        ParamRoute route = {};
        /// Jump factor and running moments of the coordinates.
        block_rw::Factor factor = {};
    };

    /// State of one Markov chain: modules built from definitions, log-likelihood totals, adaptation.
    /// * modules point into each other -> never copied, rebuilt.
    struct ChainState
    {
        /// Register the keys declared by `model` in every key holder.
        explicit ChainState(const ModelBase &model);
        ChainState(const ChainState &) = delete;
        ChainState &operator=(const ChainState &) = delete;

        //////  Modules  //////

        /// Imported tables, shared across engines, and run state.
        Data data;
        Parameter param;
        Rates rates;
        InfecProfile infec_profile;
        StatusProcess status_process;
        Augmented augm;
        Translator translator;

        //////  Totals and adaptation  //////

        /// Log-likelihood totals over households.
        LogLikTotal tot_log_lik;
        /// Route of each sampled parameter, by key.
        std::vector<ParamRoute> routes = {};
        /// Blocks of the run, see Rates::list_param_rw_blocks().
        std::vector<ParamBlock> pblocks = {};
        /// Whether each parameter belongs to a block.
        std::vector<char> in_pblock = {};
        /// Whether the adaptation warm-up runs, block members moving one at a time.
        bool is_warm = false;
        /// Robbins-Monro gain of the scalar proposals at this step.
        double gain = 0.;
        /// Robbins-Monro gain of the block factors at this step.
        double block_gain = 0.;

        //////  Setter  //////

        /// Initialize equivalence-class compression.
        /// * latent definitions required, not imported yet -> one class per household.
        void init_equivalence_class_compression();

        /// Throw if the generic profile transmits at generation time 0.
        /// * status-specific profile -> checked per member, with its variable.
        void assert_valid_infectivity_profile() const;

        //////  Getter  //////

        /// Latent `a` of member `m`.
        double aug(size_t a, const Individual &m) const { return augm.at(a, m).value(); }

        /// Data column names shared across simulations.
        std::vector<std::string> get_static_data_names() const;
        /// Column names of the per-replicate shallow list.
        std::vector<std::string> get_col_names() const;
        /// Per-household row offsets, so parallel simulations write disjoint rows.
        std::vector<size_t> get_row_offsets() const;

        //////  Log-likelihood  //////

        /// Observation-model log-likelihood of member `m`.
        const LogProba &obs_log_lik(const Individual &m) const noexcept { return augm.state.of(m).obs_log_lik; }
        /// Latent-model log-likelihood of member `m`.
        const LogProba &latent_log_lik(const Individual &m) const noexcept { return augm.state.of(m).latent_log_lik; }
        /// Sum the observation log-likelihood over infected members of `house`.
        double house_obs_log_lik(const Household &house) const
        {
            double obs = 0.;
            for (const Individual &indiv : house.infected_indivs())
                obs += obs_log_lik(indiv).value();
            return obs;
        }
        /// Sum the latent log-likelihood over infected members of `house`.
        double house_latent_log_lik(const Household &house) const
        {
            double model = 0.;
            for (const Individual &indiv : house.infected_indivs())
                model += latent_log_lik(indiv).value();
            return model;
        }
        /// Sum the transmission, observation and latent log-likelihood stored for `house`.
        /// * observation, latent -> only when the model declares them.
        double house_log_lik(const Household &house) const
        {
            double value = data.state.transm_log_lik(house).value();
            if (augm.scheme().is_any_obs_model())
                value += house_obs_log_lik(house);
            if (augm.scheme().is_any_latent_model())
                value += house_latent_log_lik(house);
            return value;
        }
    };
}

//////////////////////////////  MCMCBase  //////////////////////////////

namespace housepi
{
    class MCMCBase
    {
    public:
        Settings settings;
        Posterior posterior;

        //////  Constructor  //////

        /// Build the chain state over the keys declared by `model`.
        explicit MCMCBase(const ModelBase &model);
        virtual ~MCMCBase();

        //////  Getter  //////

        /// State of the chain.
        ChainState &chain() noexcept { return *chain_; }
        const ChainState &chain() const noexcept { return *chain_; }

    private:
        /// State of the chain.
        std::unique_ptr<ChainState> chain_;
    };
}
