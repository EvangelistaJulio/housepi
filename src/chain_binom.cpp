/* housepi - Epidemiological Analysis of Household Outbreaks.
 * Copyright 2024-2026 Institut Pasteur (Paris), France.
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <map>

#include "tools/math.h"
#include "tools/nan.h"
#include "chain_binom.h"

//////////////////////////////  chain_binom  //////////////////////////////

namespace housepi::chain_binom
{
    using std::vector;

    template <typename T>
    using matrix = vector<vector<T>>;

    namespace
    {
        /// Return factorial of `n`.
        inline unsigned long factorial(unsigned long n)
        {
            unsigned long result = 1;
            for (unsigned long i = 2; i <= n; ++i)
                result *= i;
            return result;
        }

        //////////////////////////////  Transmission generations  //////////////////////////////

        /// Return every ranked offspring generation of `n` infected members, `n` above 0.
        /// - `n = 3`: [[0, 0, 0], [0, 0, 1], [0, 1, 1], [0, 1, 2]].
        /// @todo Make it cache friendly.
        vector<gentype> recursive_span_generation_(unsigned n)
        {
            vector<gentype> output(1, gentype(n, 0));

            for (unsigned root = n - 1; root > 0; --root)
                for (gentype g : recursive_span_generation_(n - root))
                {
                    for (unsigned i = 0; i < n - root; ++i)
                        ++g[i];
                    output.emplace_back(root, 0);
                    output.back().insert(output.back().end(), g.begin(), g.end());
                }
            return output;
        }
    }

    void LookUpTables::ensure_homogeneous(MemberIdx n)
    {
        if (n < num_generations_.size())
            return;

        unsigned old_size(num_generations_.size());
        num_generations_.resize(n + 1);
        gen_multiplicity_.resize(n + 1);
        generations_.resize(n + 1);

        if (old_size == 0)
        {
            num_generations_[0] = 0;
            gen_multiplicity_[0].push_back(0);
            generations_[0].push_back({});
            old_size = 1;
        }

        vector<unsigned long> aux_num_vec;
        vector<gentype> aux_gen_vec;

        for (unsigned k = old_size; k <= n; ++k)
        {
            aux_num_vec.clear();
            aux_gen_vec.clear();

            for (const gentype &g : recursive_span_generation_(k))
            {
                std::map<MemberIdx, MemberIdx> freq_map = {};
                gentype counts = {};
                unsigned long multiplicities = factorial(k);

                for (MemberIdx val : g)
                    if (freq_map.count(val))
                        ++freq_map[val];
                    else
                        freq_map.emplace(val, 1);

                for (const auto &[val, c] : freq_map)
                {
                    counts.push_back(c);
                    multiplicities /= factorial(c);
                }

                aux_num_vec.push_back(multiplicities);
                aux_gen_vec.push_back(counts);
            }

            size_t tot_num_gen = aux_gen_vec.size();
            num_generations_[k] = tot_num_gen;
            gen_multiplicity_[k].resize(tot_num_gen);
            generations_[k].resize(tot_num_gen);

            for (size_t c = 0; c < tot_num_gen; ++c)
            {
                gen_multiplicity_[k][c] = aux_num_vec[c];
                generations_[k][c] = aux_gen_vec[c];
            }
        }
    }

    //////////////////////////////  Homogeneous log-likelihood  //////////////////////////////

    namespace
    {
        /// Compute the homogeneous chain-binomial log-likelihood.
        /// Based on @cite longini_household_1982
        /// @param multiplicity Multiplicity of each generation of `n_infec` infected.
        /// @param generations Count-compressed generations of `n_infec` infected.
        double log_lik_homogeneous_(const vector<unsigned long> &multiplicity,
                                    const vector<gentype> &generations,
                                    MemberIdx n_infec,
                                    double commu_infec_rate,
                                    double p2p_infec_rate)
        {
            const double log_commu_esc_proba = -commu_infec_rate;
            const double log_commu_infec_proba = math::log1mexp(commu_infec_rate);
            const double log_p2p_esc_proba = -p2p_infec_rate;

            // log-probability of infection by a generation of n_prev infectors, per generation size
            thread_local vector<double> log_p2p_infec_proba;
            log_p2p_infec_proba.assign(n_infec + 1, 0.);
            for (unsigned n_prev = 1; n_prev <= n_infec; ++n_prev)
                log_p2p_infec_proba[n_prev] = math::log1mexp(n_prev * p2p_infec_rate);

            //////  Seropositive contribution  //////

            // sum over all possible generations
            double sero_positive_term = 0., tmp = 0.;
            for (size_t c = 0; c < generations.size(); ++c)
            {
                // generation counts
                const auto &gen_counts = generations[c];

                // generation (0) : infected by community
                tmp = gen_counts[0] * log_commu_infec_proba;
                // generation (0) : escaped by community
                tmp += (n_infec - gen_counts[0]) * log_commu_esc_proba;

                // generation (i) : within house, i > 0
                unsigned long n_earlier = 0; // members of generations 0 <= g < i - 1
                for (size_t i = 1; i < gen_counts.size(); ++i)
                {
                    // generation (i) : infected by previous generation i - 1
                    tmp += gen_counts[i] * log_p2p_infec_proba[gen_counts[i - 1]];

                    // generation (i) : escaped from other generations g, i. e. 0 <= g < i - 1
                    tmp += gen_counts[i] * n_earlier * log_p2p_esc_proba;

                    n_earlier += gen_counts[i - 1];
                }

                // weighted by the number of equivalent configurations
                sero_positive_term += multiplicity[c] * std::exp(tmp);
            }

            return std::log(sero_positive_term);
        }

        //////////////////////////////  Heterogeneous log-likelihood  //////////////////////////////

        /// Compute the heterogeneous chain-binomial log-likelihood over subset bitmasks.
        /// Based on @cite addy_generalized_1991
        double log_lik_heterogeneous_(MemberIdx n_infec,
                                      const std::vector<double> &commu_infec_rate,
                                      const std::vector<std::vector<double>> &p2p_infec_rate)
        {
            const size_t Nsets = size_t(1) << n_infec; // number of subsets R, 2^n_infec

            // log-escape probabilities q of every subset R, built incrementally
            thread_local vector<double> log_q;
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

            thread_local vector<double> log_pi;
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

            std::map<gentype, std::uint32_t> idx_of;
            for (std::uint32_t i = 0; i < sub_vecs.size(); ++i)
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
        /// Recursion: pi(t|t) = 1 - sum_{w < t} C(t,w) * pi(w|w) * prod_s q_s(w)^(t_s - w_s).
        double evaluate_partition_tree_(const PartitionTree &dp,
                                        const vector<double> &log_commu_esc,
                                        const matrix<double> &log_p2p_esc)
        {
            const MemberIdx K = dp.sub_vecs.empty()
                                    ? 0
                                    : static_cast<MemberIdx>(dp.sub_vecs.front().size());

            thread_local vector<double> log_pi;
            log_pi.assign(dp.sub_vecs.size(), 0.); // log pi(0|0) = 0

            for (size_t ti = 1; ti < dp.sub_vecs.size(); ++ti)
            {
                const auto &t = dp.sub_vecs[ti];
                double sum = 0.;

                for (const auto &wc : dp.ancestors[ti])
                {
                    const std::uint32_t wi = wc.first;
                    const double log_C = wc.second;
                    const auto &w = dp.sub_vecs[wi];

                    double expo = 0.;
                    for (MemberIdx s = 0; s < K; ++s)
                    {
                        if (t[s] == w[s])
                            continue;
                        double log_q_s = log_commu_esc[s];
                        for (MemberIdx u = 0; u < K; ++u)
                            log_q_s += w[u] * log_p2p_esc[u][s];
                        expo += (t[s] - w[s]) * log_q_s;
                    }

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
        /// @param tree Partition tree of the part counts.
        /// @param n_parts Number of parts.
        /// @param part_indexes Member representative of each part.
        /// @param commu_infec_rate Community infection rate per member.
        /// @param p2p_infec_rate Person-to-person infection rate per member pair.
        double log_lik_partitioned_(const PartitionTree &tree,
                                    MemberIdx n_parts,
                                    const std::vector<MemberIdx> &part_indexes,
                                    const std::vector<double> &commu_infec_rate,
                                    const std::vector<std::vector<double>> &p2p_infec_rate)
        {
            // log-escape probabilities per part
            thread_local vector<double> log_commu_esc;
            log_commu_esc.assign(n_parts, 0.);
            thread_local matrix<double> log_p2p_esc;
            log_p2p_esc.assign(n_parts, vector<double>(n_parts, 0.));

            for (MemberIdx t = 0; t < n_parts; ++t)
            {
                log_commu_esc[t] = -commu_infec_rate[part_indexes[t]];
                for (MemberIdx s = 0; s < n_parts; ++s)
                    log_p2p_esc[t][s] = -p2p_infec_rate[part_indexes[t]][part_indexes[s]];
            }

            return evaluate_partition_tree_(tree, log_commu_esc, log_p2p_esc);
        }
    }

    //////////////////////////////  Log-likelihood dispatcher  //////////////////////////////

    double log_likelihood(const LookUpTables &lut,
                          MemberIdx n_infec, MemberIdx n_parts,
                          ShapeIdx shape_idx,
                          const std::vector<MemberIdx> &part_indexes,
                          const std::vector<double> &commu_infec_rate,
                          const std::vector<std::vector<double>> &p2p_infec_rate)
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
            return ll + log_lik_homogeneous_(
                            lut.gen_multiplicity_[n_infec], lut.generations_[n_infec],
                            n_infec, commu_infec_rate[0], p2p_infec_rate[0][1]);
        // partitioned case
        else if (n_parts < n_infec)
            return ll + log_lik_partitioned_(
                            lut.partition_trees_[shape_idx], n_parts, part_indexes,
                            commu_infec_rate, p2p_infec_rate);
        // heterogeneous case
        else
            return ll + log_lik_heterogeneous_(
                            n_infec, commu_infec_rate, p2p_infec_rate);
    }
}
