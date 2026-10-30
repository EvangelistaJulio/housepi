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
    class MCMCBase;
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
            matrix cols;                /// Draws per stored column.
            headers names;              /// Name of each table column.
            std::vector<size_t> source; /// Stored column of each table column.

            /// Append stored column `name` of `n_draws` NA draws.
            void add(std::string name, size_t n_draws)
            {
                cols.emplace_back(n_draws, na::DOUBLE);
                names.push_back(std::move(name));
            }
        };

        Draws param_values;  /// Parameter values.
        Draws param_priors;  /// Log-prior of each sampled parameter.
        Draws latent_values; /// Latent values of tracked members, member by member.
        /// Pointwise log-likelihood: observation + transmission + latent model.
        /// * one stored column per equivalence class, `source` giving the class of each household.
        Draws house_log_lik;
        std::vector<const Individual *> tracked_indivs; /// Members to track.

        /// Initialize posterior sample.
        void init(MCMCBase &mcmc);
        /// Store current values as draw `local` of chain `chain_id`, in its own column block.
        void store(MCMCBase &mcmc,
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
