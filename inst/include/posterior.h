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
        /// Return sample counter.
        unsigned long count() const noexcept { return count_; }

        /// Numeric storage
        std::vector<std::vector<double>> total_log_lik; // NOTE: external vector of size 1 for simplicity
        std::vector<std::vector<double>> param_values;
        std::vector<std::vector<double>> param_priors;
        std::vector<std::vector<double>> house_log_lik;
        std::vector<const Individual *> augm_individuals;
        std::vector<std::vector<std::vector<double>>> augm_values;

        /// Column name metadata captured at initialization
        const std::vector<std::string> total_log_lik_names = {"log_lik"};
        std::vector<std::string> param_names;
        std::vector<std::string> param_prior_names;
        std::vector<std::string> augm_col_names;
        std::vector<std::string> house_names;

        /// Initialize posterior sample.
        void init(MCMCBase &mcmc);
        /// Update posterior sample storing current values.
        void store(MCMCBase &mcmc);

    private:
        unsigned long count_ = 0;    /// Internal count
        bool is_power_post_ = false; /// Wheter using power posteriors
    };
}
