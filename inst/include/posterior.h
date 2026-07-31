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
        /// NOTE: for simplicity tot_* are matrices of 1 column
        /// NOTE: tot_log_lik is filled only under power posteriors as the tempered
        /// transmission + observation; ordinary inference stores the pointwise
        /// house_log_lik = transmission + observation + augmented model
        matrix tot_log_lik;                              // joint tempered log-likelihood
        matrix house_log_lik;                            // pointwise log-likelihood
        matrix param_values;                             // parameter values
        matrix param_priors;                             // single parameter priors
        std::vector<matrix> augm_values;                 // selected augmented values
        std::vector<const Individual *> tracked_indivs; // members to track

        using headers = std::vector<std::string>;
        /// Column name metadata captured at initialization
        headers param_names;
        headers param_prior_names;
        headers augm_col_names;
        headers house_names;

        /// Initialize posterior sample.
        void init(MCMCBase &mcmc);
        /// Update posterior sample storing current values.
        void store(MCMCBase &mcmc);

        /// Return sample counter.
        unsigned long count() const noexcept { return count_; }

    private:
        unsigned long count_ = 0;    /// Internal count
        bool is_power_post_ = false; /// Wheter using power posteriors
    };
}
