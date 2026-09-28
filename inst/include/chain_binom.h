/* housepi - Epidemiological Analysis of Household Outbreaks.
 *
 * Copyright © 2024-2026 Institut Pasteur (Paris), CNRS and INSERM.
 * Distributed under the terms of the GNU General Public License (GPLv3).
 * See the LICENSE file for details.
 */

#pragma once

#include <string>
#include <vector>

#include "data/individual.h"

namespace housepi::chain_binom
{
    /// Clear all lookup tables between successive evaluations.
    void clear_lookuptables();

    /// Ensure homogeneous generation tables are available up to `n` infected.
    void ensure_homogeneous(MemberIdx n);

    /// Build and store the partition tree of equivalence class `idx`.
    /// Only the first `n_parts` entries of `part_counts` are read.
    void ensure_partitioned(HouseIdx idx,
                            const std::vector<MemberIdx> &part_counts,
                            MemberIdx n_parts);

    /// Compute the chain-binomial log-likelihood with the cheapest algorithm for the partition.
    ///
    /// - `n_parts == 1`: homogeneous;
    /// - `n_parts == n_infec`: bitmask DP, heterogeneous;
    /// - otherwise: multitype partition tree, grouped.
    double log_likelihood(HouseIdx dp_idx,
                          MemberIdx n_infec,
                          MemberIdx n_parts,
                          const std::vector<MemberIdx> &part_indexes,
                          const std::vector<double> &commu_infec_rate,
                          const std::vector<std::vector<double>> &p2p_infec_rate);

    /// Return the symbolic chain-binomial polynomial, for debugging.
    ///
    /// - `B<i>`: community escape of infected `i`, 1-indexed;
    /// - `Q<i>_<j>`: escape of infectee `j` from infector `i`, 1-indexed, `i = 0` for community.
    /// @param max_n_infec Number of seropositives, below 8.
    std::string aux_debugger(unsigned max_n_infec);
}
