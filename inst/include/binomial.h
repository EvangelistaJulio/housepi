/* housepi - Epidemiological Analysis of Household Outbreaks.
 *
 * Copyright © 2026 Institut Pasteur (Paris), CNRS and INSERM.
 * Distributed under the terms of the GNU General Public License (GPLv3).
 * See the LICENSE file for details.
 */

#pragma once

#include <vector>

namespace hsp::bino
{
    /// Functions for the generations lookup table.
    void init_generations_lookuptable(unsigned max_N);

    /// Homogeneous chain binomial log-likelihood.
    /// @param n_infec n. of infected
    /// @param house_size n. of household members
    /// @param commu_esc_rate community escape rate
    /// @param p2p_esc_rate person-to-person escape rate
    double chain_binom_homogeneous_log_lik(
        unsigned n_infec,
        unsigned house_size,
        double commu_escape_rate,
        double p2p_escape_rate);

    /// Recursive likelihood formulation lookuptable.
    void init_escape_patterns_lookuptable(unsigned max_N, bool clean_combinations = true);

    /// Heterogeneous chain binomial log-likelihood for ordered probabilities (seropositives, seronegatives).
    /// @param n_infec n. of infected
    /// @param commu_esc_rate ordered community escape rates
    /// @param p2p_esc_rate ordered person-to-person escape rates
    double chain_binom_heterogeneous_log_lik(
        unsigned n_infec,
        const std::vector<double> &commu_esc_rate,
        const std::vector<std::vector<double>> &p2p_esc_rate);

    /// Grouped chain binomial log-likelihood exploiting type degeneracy.
    /// @param n_infec n. of infected
    /// @param type_counts count of each type among infected
    /// @param commu_esc_rate_type community escape rate per type
    /// @param p2p_esc_rate_type person-to-person escape rate by type pair
    /// @param commu_esc_rate per-individual community escape rates (all members)
    /// @param p2p_esc_rate per-individual p2p escape rates (all members)
    double chain_binom_grouped_log_lik(
        unsigned n_infec,
        const std::vector<unsigned> &type_counts,
        const std::vector<double> &commu_esc_rate_type,
        const std::vector<std::vector<double>> &p2p_esc_rate_type,
        const std::vector<double> &commu_esc_rate,
        const std::vector<std::vector<double>> &p2p_esc_rate);
}
