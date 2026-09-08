/* housepi - Epidemiological Analysis of Household Outbreaks.
 *
 * Copyright © 2024-2026 Institut Pasteur (Paris), CNRS and INSERM.
 * Distributed under the terms of the GNU General Public License (GPLv3).
 * See the LICENSE file for details.
 */

#pragma once

#include <string>
#include <vector>

#include "tools/types.h"

namespace housepi::chain_binom
{
    /// Clear all lookup tables between successive evaluations.
    void clear_lookuptables();

    /// Ensure homogeneous generation tables are available up to n infected.
    void ensure_homogeneous(MemberIdx n);

    /// Pre-build and store the partition tree for equivalence class idx.
    /// Only the first n_parts entries of part_counts are read.
    void ensure_partitioned(HouseIdx idx,
                            const std::vector<MemberIdx> &part_counts,
                            MemberIdx n_parts);

    /// Chain binomial log-likelihood dispatcher.
    ///
    /// Routes to the cheapest algorithm consistent with the partition:
    /// `n_parts == 1` → homogeneous;
    /// `n_parts == n_infec` → bitmask DP heterogeneous;
    /// otherwise → multitype partition tree (grouped).
    double log_likelihood(HouseIdx dp_idx,
                          MemberIdx n_infec,
                          MemberIdx n_parts,
                          const std::vector<MemberIdx> &part_indexes,
                          const std::vector<double> &commu_infec_rate,
                          const std::vector<std::vector<double>> &p2p_infec_rate);

    /// Symbolic chain-binomial polynomial for debugging.
    /// @param max_n_infec number of seropositives (< 8)
    /// @return formula string with tokens B<i> (community escape of infected i, 1-indexed)
    ///         and Q<i>_<j> (p2p escape from infector i to infectee j, 1-indexed; i=0 = community)
    std::string aux_debugger(unsigned max_n_infec);
}
