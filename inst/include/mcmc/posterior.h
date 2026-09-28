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
    struct Parameter;
    struct Augmented;
    namespace data
    {
        struct Wide;
    }
    struct Setup;

    /// Store posterior MCMC draws of every chain.
    struct Posterior
    {
        using matrix = std::vector<std::vector<double>>;
        /// Pointwise log-likelihood: transmission, observation and latent model.
        matrix house_log_lik;
        matrix param_values;                            /// Parameter values.
        matrix param_priors;                            /// Log-prior of each parameter.
        std::vector<matrix> latent_values;              /// Latent values of tracked members.
        std::vector<const Individual *> tracked_indivs; /// Members to track.

        //////  Column names, captured at initialization  //////

        using headers = std::vector<std::string>;
        headers param_names;       /// Parameter names.
        headers param_prior_names; /// Parameter prior names.
        headers latent_col_names;  /// Latent column names.
        headers house_names;       /// Household names.

        /// Initialize posterior sample.
        void init(MCMCBase &mcmc);
        /// Store current values as draw `local` of chain `chain_id`, in its own column block.
        void store(MCMCBase &mcmc, unsigned long chain_id, unsigned long local);

        /// Return the number of draws of the whole sample, over all chains.
        unsigned long sample_size() const noexcept { return sample_size_; }

    private:
        unsigned long sample_size_ = 0; /// Draws of the whole sample.

        /// Initialize the parameter draws of the posterior.
        void init_param_draws_(const Parameter &param,
                               unsigned long sample_size);
        /// Initialize the household transmission log-likelihood draws of the posterior.
        void init_house_draws_(const data::Wide &data_wide,
                               unsigned long sample_size);
        /// Initialize the latent draws of the posterior.
        void init_latent_draws_(const Augmented &augm,
                                const data::Wide &data_wide,
                                const Setup &setup,
                                unsigned long sample_size);
    };
}
