/* housepi - Epidemiological Analysis of Household Outbreaks.
 *
 * Copyright © 2024-2026 Institut Pasteur (Paris), CNRS and INSERM.
 * Distributed under the terms of the GNU General Public License (GPLv3).
 * See the LICENSE file for details.
 */

#pragma once

#include <string>
#include <vector>

#include "tools/nan.h"

namespace housepi
{
    class Individual;
    class MCMCBase;

    /// Centralised storage for posterior MCMC draws.
    struct Posterior
    {
        using matrix = std::vector<std::vector<double>>;
        /// house_log_lik = transmission + observation + latent model
        matrix house_log_lik;                            // pointwise log-likelihood
        matrix param_values;                             // parameter values
        matrix param_priors;                             // single parameter priors
        std::vector<matrix> latent_values;                 // selected latent values
        std::vector<const Individual *> tracked_indivs; // members to track

        using headers = std::vector<std::string>;
        /// Column name metadata captured at initialization
        headers param_names;
        headers param_prior_names;
        headers latent_col_names;
        headers house_names;

        /// Initialize posterior sample.
        void init(MCMCBase &mcmc);
        /// Update posterior sample storing current values.
        /// Chain `chain_id` writes its own column block, `local` being the
        /// index of the draw within that chain.
        void store(MCMCBase &mcmc, unsigned long chain_id, unsigned long local);

        /// Return the number of draws of the whole sample, over all chains.
        unsigned long sample_size() const noexcept { return sample_size_; }

    private:
        unsigned long sample_size_ = 0; /// Draws of the whole sample
    };
}
