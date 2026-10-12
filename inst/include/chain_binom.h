/* housepi - Epidemiological Analysis of Household Outbreaks.
 * Copyright 2024-2026 Institut Pasteur (Paris), France.
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#pragma once

#include <cstdint>
#include <map>
#include <utility>
#include <vector>

#include "data/individual.h"

namespace housepi::chain_binom
{
    using gentype = std::vector<MemberIdx>; /// Count vector, one entry per part.
    using ShapeIdx = std::uint32_t;         /// Index of a shape (part-count vector) in LookUpTables.

    /// Cap to the within-household risk partitions.
    ///
    /// * heterogeneous -> 3^K steps per evaluation, costs 0.5 s at K = 15;
    /// * partitioned -> at least 2^K tree nodes.
    constexpr MemberIdx max_parts = 15;

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

    /// Hold partition trees of one model, shared by its chains.
    ///
    /// * written while the model builds, serially;
    /// * read-only during evaluation.
    struct LookUpTables
    {
    public:
        //////  Setter  //////

        /// Return the shape index of `part_counts`, building its partition tree if absent.
        /// Only the first `n_parts` entries of `part_counts` are read.
        ShapeIdx ensure_partitioned(const std::vector<MemberIdx> &part_counts,
                                    MemberIdx n_parts);

    private:
        //////  Cached partition trees  //////

        std::vector<PartitionTree> partition_trees_; /// Partition trees, one per shape.
        std::map<gentype, ShapeIdx> shapes_;         /// Shape index, per part-count vector.

        /// Read the tables.
        friend double log_likelihood(const LookUpTables &lut,
                                     MemberIdx n_infec,
                                     MemberIdx n_parts,
                                     ShapeIdx shape_idx,
                                     const std::vector<MemberIdx> &part_indexes,
                                     const std::vector<double> &commu_infec_rate,
                                     const std::vector<std::vector<double>> &p2p_infec_rate);
    };

    /// Compute the chain-binomial log-likelihood with the cheapest algorithm for the partition.
    ///
    /// * `n_parts == 1` -> generation sizes, homogeneous;
    /// * `n_parts == n_infec` -> bitmask DP, heterogeneous;
    /// * otherwise -> multitype partition tree, grouped.
    /// @param shape_idx Shape index from LookUpTables::ensure_partitioned(), read only when grouped.
    double log_likelihood(const LookUpTables &lut,
                          MemberIdx n_infec,
                          MemberIdx n_parts,
                          ShapeIdx shape_idx,
                          const std::vector<MemberIdx> &part_indexes,
                          const std::vector<double> &commu_infec_rate,
                          const std::vector<std::vector<double>> &p2p_infec_rate);
}
