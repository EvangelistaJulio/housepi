/* housepi - Epidemiological Analysis of Household Outbreaks.
 *
 * Copyright © 2024-2026 Institut Pasteur (Paris), CNRS and INSERM.
 * Distributed under the terms of the GNU General Public License (GPLv3).
 * See the LICENSE file for details.
 */

#pragma once

#include <cstdint>
#include <string>
#include <utility>
#include <vector>

#include "data/individual.h"

namespace housepi::chain_binom
{
    using gentype = std::vector<MemberIdx>; /// Count vector, one entry per generation or per part.

    //////////////////////////////  chain_binom::PartitionTree  //////////////////////////////

    /// Hold the multitype partition tree of a chain binomial.
    /// Based on @cite ball_unified_1986 and @cite addy_generalized_1991
    /// @todo Flatten the tree.
    struct PartitionTree
    {
        /// Sub-vectors of the root, in a linear extension of the componentwise order.
        std::vector<gentype> sub_vecs;
        /// Per sub-vector `t`: each smaller `w`, as its index and log C(t, w).
        std::vector<std::vector<std::pair<std::uint32_t, double>>> ancestors;
    };

    //////////////////////////////  chain_binom::LookUpTables  //////////////////////////////

    /// Hold lookup tables of one model, shared by its chains.
    ///
    /// - written while the model builds, serially;
    /// - read-only during evaluation.
    struct LookUpTables
    {
        //////  Transmission generations  //////

        std::vector<unsigned long> num_generations;               /// Number of generations, per infected count.
        std::vector<std::vector<unsigned long>> gen_multiplicity; /// Multiplicity of each generation.
        std::vector<std::vector<gentype>> generations;            /// Count-compressed generations, per infected count.

        //////  Cached partition trees  //////

        std::vector<PartitionTree> partition_trees; /// Partition trees, per equivalence class.

        //////  Setter  //////

        /// Ensure homogeneous generation tables are available up to `n` infected.
        void ensure_homogeneous(MemberIdx n);

        /// Build and store the partition tree of equivalence class `idx`.
        /// Only the first `n_parts` entries of `part_counts` are read.
        void ensure_partitioned(HouseIdx idx,
                                const std::vector<MemberIdx> &part_counts,
                                MemberIdx n_parts);
    };

    /// Compute the chain-binomial log-likelihood with the cheapest algorithm for the partition.
    ///
    /// - `n_parts == 1`: homogeneous;
    /// - `n_parts == n_infec`: bitmask DP, heterogeneous;
    /// - otherwise: multitype partition tree, grouped.
    double log_likelihood(const LookUpTables &lut,
                          HouseIdx dp_idx,
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
