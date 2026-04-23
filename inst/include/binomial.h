/* housepi - Epidemiological Analysis of Household Outbreaks.
 *
 * Copyright © 2026 Institut Pasteur (Paris), CNRS and INSERM.
 * Distributed under the terms of the GNU General Public License (GPLv3).
 * See the LICENSE file for details.
 */

#pragma once

#include <string>
#include <vector>

namespace housepi::bino
{
    /// Functions for the generations lookup table.
    void init_generations_lookuptable(unsigned max_N);

    /// Chain binomial log-likelihood dispatcher.
    ///
    /// Routes to the cheapest algorithm consistent with the partition:
    /// `n_parts == 1` with uniform rates → closed-form homogeneous;
    /// `n_parts == n_infec` → bitmask DP heterogeneous;
    /// otherwise → multitype final-size DP (grouped).
    double chain_binom_log_lik(
        size_t n_infec, size_t n_parts,
        const std::vector<size_t> &part_counts,
        const std::vector<size_t> &part_indexes,
        const std::vector<double> &commu_esc_rate,
        const std::vector<std::vector<double>> &p2p_esc_rate);

    /// Symbolic chain-binomial polynomial for debugging.
    /// @param max_n_infec number of seropositives (< 8)
    /// @return formula string with tokens B<i> (community escape of infected i, 1-indexed)
    ///         and Q<i>_<j> (p2p escape from infector i to infectee j, 1-indexed; i=0 = community)
    std::string auxiliary_chain_binomial_debugger(unsigned max_n_infec);
}
