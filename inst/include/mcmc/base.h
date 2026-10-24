/* housepi - Epidemiological Analysis of Household Outbreaks.
 * Copyright 2024-2026 Institut Pasteur (Paris), France.
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#pragma once

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
#include "risks.h"
#include "status_process.h"
#include "translator.h"

//////////////////////////////  MCMCBase  //////////////////////////////

namespace housepi
{
    class ModelBase;

    class MCMCBase
    {
    public:
        Settings settings;
        Posterior posterior;
        LogLikTotal tot_log_lik; /// Log-likelihood totals over households.

        //////  Key holders  //////

        Data data; /// Imported tables.
        Parameter param;
        Risks risks;
        InfecProfile infec_profile;
        StatusProcess status_process;
        Augmented augm;
        Translator translator;

        //////  Constructor  //////

        /// Register the keys declared by `model` in every key holder.
        explicit MCMCBase(const ModelBase &model);
        virtual ~MCMCBase();

        //////  Data  //////

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
        /// Sum the observation log-likelihood over infected members of `house`.
        double house_obs_log_lik(const Household &house) const
        {
            double obs = 0.;
            for (const Individual &indiv : house.infected_indivs())
                obs += obs_log_lik_(indiv).value();
            return obs;
        }
        /// Sum the latent log-likelihood over infected members of `house`.
        double house_latent_log_lik(const Household &house) const
        {
            double model = 0.;
            for (const Individual &indiv : house.infected_indivs())
                model += latent_log_lik_(indiv).value();
            return model;
        }

    protected:
        //////  Log-likelihood  //////

        /// Observation-model log-likelihood of individual `m`.
        const LogProba &obs_log_lik_(const Individual &m) const noexcept
        {
            return augm.state.of(m).obs_log_lik;
        }
        /// Latent-model log-likelihood of individual `m`.
        const LogProba &latent_log_lik_(const Individual &m) const noexcept
        {
            return augm.state.of(m).latent_log_lik;
        }
    };
}
