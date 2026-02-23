/* housepi - Epidemiological Analysis of Household Outbreaks.
 *
 * Copyright © 2026 Institut Pasteur (Paris), CNRS and INSERM.
 * Distributed under the terms of the GNU General Public License (GPLv3).
 * See the LICENSE file for details.
 */

#pragma once

#include <vector>

namespace bax
{
    /// Functions for the generations lookup table.
    void init_generations_lookuptable(unsigned max_N);

    unsigned long get_num_generations(unsigned n_positive);
    std::vector<unsigned> get_generation_vector(unsigned n_positive, unsigned long i);
    unsigned long get_generation_multipl(unsigned n_positive, unsigned long i);

    /// Recursive likelihood formulation lookuptable.
    void init_escape_patterns_lookuptable(unsigned max_N, bool clean_combinations = true);

    /// Heterogeneous chain binomial log-likelihood for ordered probabilities (seropositives, seronegatives).
    /// @param N_star n. of seropositves
    /// @param commu_esc_rate ordered community escape rates
    /// @param p2p_esc_rate ordered person-to-person escape rates
    double chain_binomial_heter_log_lik(unsigned N_star,
                                        const std::vector<double> &commu_esc_rate,
                                        const std::vector<std::vector<double>> &p2p_esc_rate);
}
