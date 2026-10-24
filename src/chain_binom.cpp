/* housepi - Epidemiological Analysis of Household Outbreaks.
 * Copyright 2024-2026 Institut Pasteur (Paris), France.
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <map>
#include <numeric>

#include "tools/math.h"
#include "tools/messages.h"
#include "tools/nan.h"
#include "chain_binom.h"

//////////////////////////////  chain_binom  //////////////////////////////

namespace housepi::chain_binom
{
    using std::vector;

    //////////////////////////////  Homogeneous log-likelihood  //////////////////////////////

    namespace
    {
        /// Compute the homogeneous chain-binomial log-likelihood of `n` infected members.
        ///
        /// Model based on @cite longini_household_1982
        /// Run a Reed-Frost chain seeded by community infections at generation 0.
        /// Final size ignores timing, so seeding at generation 0 is exact, after @cite ball_unified_1986
        /// Sum over generation sizes `c_0, c_1, ...` adding up to `n`.
        /// Sequences are grouped by state, so the sum takes O(n^3) products:
        /// * `s`: members infected before the last generation;
        /// * `l`: members of the last generation.
        /// @note Probabilities below 1e-308 underflow, returning -Inf.
        double log_lik_homogeneous_(Scratch &scratch,
                                    MemberIdx n,
                                    double commu_infec_rate,
                                    double p2p_infec_rate)
        {
            const double commu_infec_proba = -math::expm1(-commu_infec_rate);
            const double commu_esc_proba = std::exp(-commu_infec_rate);
            const double p2p_esc_proba = std::exp(-p2p_infec_rate);

            // probability of infection by a generation of n_prev infectors, per generation size
            auto &p2p_infec_proba = scratch.p2p_infec_proba;
            p2p_infec_proba.assign(n + 1, 0.);
            for (unsigned n_prev = 1; n_prev <= n; ++n_prev)
                p2p_infec_proba[n_prev] = -math::expm1(-(n_prev * p2p_infec_rate));

            // probability of escaping the community, per number of members
            auto &commu_esc_pow = scratch.commu_esc_pow;
            commu_esc_pow.assign(n + 1, 1.);
            for (unsigned k = 1; k <= n; ++k)
                commu_esc_pow[k] = commu_esc_pow[k - 1] * commu_esc_proba;

            // probability of reaching state (s, l), at s * (n + 1) + l
            auto &state_proba = scratch.state_proba;
            state_proba.assign((n + 1) * (n + 1), 0.);

            // generation (0) : C(n, c) ways for c members infected by community, the others escaped
            double ways = 1.;
            for (unsigned c = 1; c <= n; ++c)
            {
                ways *= commu_infec_proba * (n - c + 1) / c;
                state_proba[c] = ways * commu_esc_pow[n - c];
            }

            // generation (i) : C(n - s - l, c) ways for c members infected by generation i - 1
            double total = 0., p2p_esc_earlier = 1.; // escape from the s members before generation i - 1
            for (unsigned s = 0; s < n; ++s, p2p_esc_earlier *= p2p_esc_proba)
                for (unsigned l = 1; s + l <= n; ++l)
                {
                    const double prev = state_proba[s * (n + 1) + l];
                    if (s + l == n)
                    {
                        total += prev;
                        continue;
                    }
                    const double member_proba = p2p_infec_proba[l] * p2p_esc_earlier;
                    const unsigned n_left = n - s - l;
                    ways = prev;
                    for (unsigned c = 1; c <= n_left; ++c)
                    {
                        ways *= member_proba * (n_left - c + 1) / c;
                        state_proba[(s + l) * (n + 1) + c] += ways;
                    }
                }

            return std::log(total);
        }

        //////////////////////////////  Heterogeneous log-likelihood  //////////////////////////////

        /// Compute the heterogeneous chain-binomial log-likelihood over subset bitmasks.
        /// Based on @cite addy_generalized_1991
        double log_lik_heterogeneous_(Scratch &scratch,
                                      MemberIdx n_infec,
                                      const std::vector<double> &commu_infec_rate,
                                      const matrix &p2p_infec_rate)
        {
            const size_t Nsets = size_t(1) << n_infec; // number of subsets R, 2^n_infec

            // log-escape probabilities q of every subset R, built incrementally
            auto &log_q = scratch.log_q;
            log_q.assign(Nsets * n_infec, 0.);

            // empty R: community escape only, O(n)
            for (size_t j = 0; j < n_infec; ++j)
                log_q[j] = -commu_infec_rate[j];

            // non-empty R: O(2^n * n)
            for (size_t R = 1; R < Nsets; ++R)
            {
                const size_t low = R & (~R + 1);                            // isolate lowest set bit of R
                const size_t i = static_cast<size_t>(__builtin_ctzll(low)); // index of that bit (= some member of R)
                const size_t Rprev = R ^ low;                               // R with that member removed
                for (size_t j = 0; j < n_infec; ++j)
                    log_q[R * n_infec + j] = log_q[Rprev * n_infec + j] - p2p_infec_rate[i][j];
            }

            auto &log_pi = scratch.log_pi;
            log_pi.assign(Nsets, 0.); // log pi(empty) = 0

            for (size_t S = 1; S < Nsets; ++S)
            {
                double sum = 0.;
                // enumerate all proper subsets of S in decreasing order
                for (size_t R = (S - 1) & S;; R = (R - 1) & S)
                {
                    const size_t SmR = S & ~R; // bitmask for S \ R
                    double expo = 0.;
                    for (size_t j = 0; j < n_infec; ++j)
                        if ((SmR >> j) & size_t(1)) // j in S \ R
                            expo += log_q[R * n_infec + j];
                    sum += std::exp(log_pi[R] + expo);
                    if (R == 0) // last subset
                        break;
                }
                log_pi[S] = (sum < 1.0) ? std::log1p(-sum) : -inf::DOUBLE;
            }

            return log_pi[Nsets - 1];
        }

        //////////////////////////////  chain_binom::PartitionTree : build and evaluate  //////////////////////////////

        /// Build the partition tree for the infected partition held in c_root[0, K).
        PartitionTree build_partition_tree_(const MemberIdx *c_root,
                                            MemberIdx K)
        {
            size_t total = 1;
            for (MemberIdx s = 0; s < K; ++s)
                total *= (c_root[s] + 1);

            // enumerate all v <= c_root via odometer increment
            vector<gentype> sub_vecs;
            sub_vecs.reserve(total);
            {
                gentype cur(K, 0);
                while (true)
                {
                    sub_vecs.push_back(cur);
                    MemberIdx s = 0;
                    while (s < K && cur[s] == c_root[s])
                        cur[s++] = 0;
                    if (s == K)
                        break;
                    ++cur[s];
                }
            }

            // linear extension of the componentwise order: |v|_1 ascending, lex as tiebreak
            std::sort(sub_vecs.begin(), sub_vecs.end(),
                      [](const gentype &a, const gentype &b)
                      {
                          unsigned sa = 0, sb = 0;
                          for (MemberIdx x : a)
                              sa += x;
                          for (MemberIdx x : b)
                              sb += x;
                          return sa != sb ? sa < sb : a < b;
                      });

            std::map<gentype, SubVecIdx> idx_of;
            for (SubVecIdx i = 0; i < sub_vecs.size(); ++i)
                idx_of[sub_vecs[i]] = i;

            PartitionTree dp;
            dp.sub_vecs = std::move(sub_vecs);
            dp.ancestors.resize(dp.sub_vecs.size());

            for (size_t ti = 1; ti < dp.sub_vecs.size(); ++ti)
            {
                const auto &t = dp.sub_vecs[ti];
                gentype w(K, 0);
                while (true)
                {
                    if (w != t)
                    {
                        double log_C = 0.;
                        for (MemberIdx s = 0; s < K; ++s)
                            log_C += std::lgamma(t[s] + 1) - std::lgamma(w[s] + 1) -
                                     std::lgamma(t[s] - w[s] + 1);
                        dp.ancestors[ti].emplace_back(idx_of[w], log_C);
                    }
                    MemberIdx s = 0;
                    while (s < K && w[s] == t[s])
                        w[s++] = 0;
                    if (s == K)
                        break;
                    ++w[s];
                }
            }

            return dp;
        }
    }

    ShapeIdx LookUpTables::ensure_partitioned(const std::vector<MemberIdx> &part_counts,
                                              MemberIdx n_parts)
    {
        auto [it, is_new] = shapes_.try_emplace(
            gentype(part_counts.begin(), part_counts.begin() + n_parts),
            static_cast<ShapeIdx>(partition_trees_.size()));
        if (is_new)
            partition_trees_.push_back(build_partition_tree_(part_counts.data(), n_parts));
        return it->second;
    }

    namespace
    {
        /// Evaluate log P(final = c_root | all infecteds in c_root).
        ///
        /// Recursion: pi(t|t) = 1 - sum_{w < t} C(t,w) * pi(w|w) * prod_s q_s(w)^(t_s - w_s).
        /// Here q_s(w) is the probability that a member of part s escapes the community and every member in w.
        double evaluate_partition_tree_(Scratch &scratch,
                                        const PartitionTree &dp,
                                        const vector<double> &log_commu_esc,
                                        const matrix &log_p2p_esc)
        {
            const MemberIdx K = dp.sub_vecs.empty()
                                    ? 0
                                    : static_cast<MemberIdx>(dp.sub_vecs.front().size());

            auto &log_pi = scratch.log_pi;
            log_pi.assign(dp.sub_vecs.size(), 0.); // log pi(0|0) = 0

            // log q_s(w) per node w and part s, at w * K + s
            auto &log_q = scratch.log_q;
            log_q.assign(dp.sub_vecs.size() * K, 0.);
            for (size_t wi = 0; wi < dp.sub_vecs.size(); ++wi)
                for (MemberIdx s = 0; s < K; ++s)
                {
                    double log_q_s = log_commu_esc[s];
                    for (MemberIdx u = 0; u < K; ++u)
                        log_q_s += dp.sub_vecs[wi][u] * log_p2p_esc[u][s];
                    log_q[wi * K + s] = log_q_s;
                }

            for (size_t ti = 1; ti < dp.sub_vecs.size(); ++ti)
            {
                const auto &t = dp.sub_vecs[ti];
                double sum = 0.;

                for (const auto &wc : dp.ancestors[ti])
                {
                    const SubVecIdx wi = wc.first;
                    const double log_C = wc.second;
                    const auto &w = dp.sub_vecs[wi];

                    double expo = 0.;
                    for (MemberIdx s = 0; s < K; ++s)
                        if (t[s] != w[s])
                            expo += (t[s] - w[s]) * log_q[wi * K + s];

                    sum += std::exp(log_C + log_pi[wi] + expo);
                }
                log_pi[ti] = (sum < 1.0) ? std::log1p(-sum) : -inf::DOUBLE;
            }

            return log_pi.back();
        }

        //////////////////////////////  Partitioned log-likelihood  //////////////////////////////

        /// Compute the grouped chain-binomial log-likelihood.
        ///
        /// Seropositives sharing a contribution enumerate partition sequences, not permutations.
        /// @note Subset bitmasks do not carry over: integer partitions have no closed-form count.
        /// @param scratch Buffers of the calling thread.
        /// @param tree Partition tree of the part counts.
        /// @param n_parts Number of parts.
        /// @param part_indexes Member representative of each part.
        /// @param commu_infec_rate Community infection rate per member.
        /// @param p2p_infec_rate Person-to-person infection rate per member pair.
        double log_lik_partitioned_(Scratch &scratch,
                                    const PartitionTree &tree,
                                    MemberIdx n_parts,
                                    const std::vector<MemberIdx> &part_indexes,
                                    const std::vector<double> &commu_infec_rate,
                                    const matrix &p2p_infec_rate)
        {
            // log-escape probabilities per part
            auto &log_commu_esc = scratch.log_commu_esc;
            log_commu_esc.assign(n_parts, 0.);
            auto &log_p2p_esc = scratch.log_p2p_esc;
            log_p2p_esc.assign(n_parts, vector<double>(n_parts, 0.));

            for (MemberIdx t = 0; t < n_parts; ++t)
            {
                log_commu_esc[t] = -commu_infec_rate[part_indexes[t]];
                for (MemberIdx s = 0; s < n_parts; ++s)
                    log_p2p_esc[t][s] = -p2p_infec_rate[part_indexes[t]][part_indexes[s]];
            }

            return evaluate_partition_tree_(scratch, tree, log_commu_esc, log_p2p_esc);
        }
    }

    //////////////////////////////  Cost limit  //////////////////////////////

    void assert_tractable(const std::vector<MemberIdx> &part_counts,
                          MemberIdx n_parts,
                          const std::string &where)
    {
        std::uint64_t terms = 1, nodes = 1;
        for (MemberIdx s = 0; s < n_parts && terms - nodes <= MAX_TERMS; ++s)
        {
            const std::uint64_t c = part_counts[s];
            terms *= (c + 1) * (c + 2) / 2;
            nodes *= c + 1;
        }
        if (terms - nodes > MAX_TERMS)
        {
            const size_t n_infec = std::accumulate(part_counts.begin(), part_counts.begin() + n_parts, size_t{0});
            msg::stop(
                "The " + std::to_string(n_infec) + " infected members of " + where + " fall into " +
                std::to_string(n_parts) + " risk partitions, too many to evaluate.\n" +
                "Household members share a partition if their group and contact risks match.\n" +
                msg::alert("Computation time grows exponentially with the number of partitions.") +
                "\nReduce the number of strata.");
        }
    }

    //////////////////////////////  Log-likelihood dispatcher  //////////////////////////////

    double log_likelihood(const LookUpTables &lut,
                          Scratch &scratch,
                          MemberIdx n_infec, MemberIdx n_parts,
                          ShapeIdx shape_idx,
                          const std::vector<MemberIdx> &part_indexes,
                          const std::vector<double> &commu_infec_rate,
                          const matrix &p2p_infec_rate)
    {
        //////  Seronegative contribution  //////

        double ll = 0.;
        for (size_t i = n_infec; i < commu_infec_rate.size(); ++i)
        {
            // community : i escaped from community
            ll -= commu_infec_rate[i];
            // person-to-person : i escaped from k
            for (size_t k = 0; k < n_infec; ++k)
                ll -= p2p_infec_rate[k][i];
        }

        //////  Seropositive contribution  //////

        // no infected
        if (n_infec == 0)
            return ll;
        // single seropositive
        else if (n_infec == 1)
            return ll + math::log1mexp(commu_infec_rate[0]);
        // homogeneous case
        else if (n_parts < 2)
            return ll + log_lik_homogeneous_(scratch, n_infec, commu_infec_rate[0], p2p_infec_rate[0][1]);
        // partitioned case
        else if (n_parts < n_infec)
            return ll + log_lik_partitioned_(
                            scratch, lut.partition_trees_[shape_idx], n_parts, part_indexes,
                            commu_infec_rate, p2p_infec_rate);
        // heterogeneous case
        else
            return ll + log_lik_heterogeneous_(
                            scratch, n_infec, commu_infec_rate, p2p_infec_rate);
    }
}
