/* housepi - Epidemiological Analysis of Household Outbreaks.
 *
 * Copyright © 2026 Institut Pasteur (Paris), CNRS and INSERM.
 * Distributed under the terms of the GNU General Public License (GPLv3).
 * See the LICENSE file for details.
 */

#pragma once

#include <vector>

namespace housepi::bino
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

    /// Grouped chain binomial log-likelihood.
    ///
    /// Exploits type degeneracy: when multiple seropositives share the same type, the algorithm enumerates
    /// distinct type-composition sequences per generation instead of all individual
    /// permutations. Multinomial coefficients account for within-type interchangeability.
    ///
    /// @param n_infec number of infected members
    /// @param n_types number of types
    /// @param type_counts counts of in-house seropositives
    /// @param type_indexes indexes of household member representative for type
    /// @param commu_esc_rate per-individual community escape rates
    /// @param p2p_esc_rate per-individual p2p escape rates
    double chain_binom_grouped_log_lik(
        size_t n_infec, size_t n_types,
        const std::vector<size_t> &type_counts,
        const std::vector<size_t> &type_indexes,
        const std::vector<double> &commu_esc_rate,
        const std::vector<std::vector<double>> &p2p_esc_rate);
}
