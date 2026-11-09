/* housepi - Epidemiological Analysis of Household Outbreaks.
 * Copyright 2024-2026 Institut Pasteur (Paris), France.
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#pragma once

#include <string>
#include <utility>
#include <vector>

#include "tools/nan.h"
#include "tools/vector.h"
#include "mcmc/settings.h"

namespace housepi
{
    class Individual;
    struct ChainState;
    struct Data;
    struct Parameter;
    struct Augmented;
    namespace data
    {
        struct Dataset;
    }

    /// Store posterior MCMC draws of every chain.
    struct Posterior
    {
        using headers = std::vector<std::string>;

        /// Draw columns of one table, shaped at initialization.
        /// * table column `k` reads stored column `source[k]`, or `k` when `source` is empty.
        struct Draws
        {
            /// Draws per stored column.
            matrix cols;
            /// Name of each table column.
            headers names;
            /// Stored column of each table column.
            std::vector<size_t> source;

            /// Append stored column `name` of `n_draws` NA draws.
            void add(std::string name, size_t n_draws)
            {
                cols.emplace_back(n_draws, na::DOUBLE);
                names.push_back(std::move(name));
            }
        };

        Draws param_values;
        /// Log-prior of each sampled parameter.
        Draws param_priors;
        /// Latent values of tracked members, member by member.
        Draws latent_values;
        /// Pointwise log-likelihood: observation + transmission + latent model.
        /// * one stored column per equivalence class, `source` giving the class of each household.
        Draws house_log_lik;
        std::vector<const Individual *> tracked_indivs;

        /// Initialize posterior sample of `chain`.
        void init(const ChainState &chain, const Settings &settings);
        /// Store the values of `chain` as draw `local` of chain `chain_id`, in its own column block.
        void store(const ChainState &chain,
                   const Settings &settings,
                   Settings::Chain chain_id,
                   Settings::Step local);

    private:
        /// Initialize the parameter draws of the posterior.
        void init_param_draws_(const Parameter &param,
                               Settings::Step sample_size);
        /// Initialize the household transmission log-likelihood draws of the posterior.
        void init_house_draws_(const Data &data,
                               Settings::Step sample_size);
        /// Initialize the latent draws of the posterior.
        void init_latent_draws_(const Augmented &augm,
                                const data::Dataset &dataset,
                                const Settings &settings,
                                Settings::Step sample_size);
    };
}
