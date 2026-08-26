/* housepi - Epidemiological Analysis of Household Outbreaks.
 *
 * Copyright © 2024-2026 Institut Pasteur (Paris), CNRS and INSERM.
 * Distributed under the terms of the GNU General Public License (GPLv3).
 * See the LICENSE file for details.
 */

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <numeric>
#include <map>
#include <sstream>
#include <stdexcept>

#include "binomial.h"
#include "tools/messages.h"
#include "tools/nan.h"

/////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////// *** CHAIN-BINOMIAL AUXILIARY NAMESPACE *** //////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////

namespace housepi::chain_binom
{
    using std::vector;
    using gentype = vector<MemberIndex>;

    template <typename T>
    using matrix = vector<vector<T>>;

    /// Standard definition for the factorial.
    inline unsigned long factorial(unsigned long n)
    {
        unsigned long result = 1;
        for (unsigned long i = 2; i <= n; ++i)
            result *= i;
        return result;
    }

    /// Binomial coefficient of N elements in k positions
    inline unsigned long binom_coeff(unsigned N, unsigned k)
    {
        unsigned long result = 1;
        // note: use symmetry
        if (k > N - k)
            k = N - k;
        for (unsigned i = 1; i <= k; ++i)
            result = result * (N - k + i) / i;
        return result;
    }

    /////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////// *** EscapePattern : Class *** /////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////

    class EscapePattern
    {
    public:
        EscapePattern(unsigned max_n_infec)
            : max_n_infec_(max_n_infec) {}
        ~EscapePattern() {}

        /// Add a new entry for the escape matrix.
        void add_entry(unsigned new_entry)
        {
            entries_.push_back(new_entry);
            ++n_entries_;
        }
        /// Add a new pattern index.
        void add_pattern(unsigned long new_pattern_idx)
        {
            pattern_idxs_.push_back(new_pattern_idx);
            ++n_patterns_;
        }

        /// String with the formula of the likelihood.
        std::string str() const;

    private:
        const unsigned max_n_infec_;
        size_t n_entries_ = 0;
        size_t n_patterns_ = 0;

        vector<unsigned> entries_ = {};
        vector<unsigned long> pattern_idxs_ = {};
    };

    /////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////  PartitionTree : Class  ////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////

    /// Multitype partition tree for the respective chain-binomial.
    /// Based on @cite ball_unified_1986
    /// UPDATE: this could be flattened.
    struct PartitionTree
    {
        vector<gentype> sub_vecs;
        vector<vector<std::pair<std::uint32_t, double>>> ancestors;
    };

    /////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////  Look-Up Tables  ///////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////

    struct LookUpTables
    {
        /////////////////////////////////////////////
        /////////  Transmission Generations  ////////
        /////////////////////////////////////////////

        vector<unsigned long> num_generations;         /// n. of generations
        matrix<unsigned long> multiplicity_generation; /// n. of multiplicites per generation
        matrix<gentype> generations;                   /// generation vectors

        ///////////////////////////////////////////
        /////////  Cached Partition Trees  ////////
        ///////////////////////////////////////////

        /// Pre-build partition trees indexed by equivalence class.
        vector<PartitionTree> partition_trees;

        /////////////////////////////////
        /////////  Combinations  ////////
        /////////////////////////////////

        matrix<unsigned long> num_combinations; /// LUT wiht n.e number of combinations
        matrix<matrix<unsigned>> combinations;  /// LUT wiht combinations and complements

        ////////////////////////////////////
        /////////  Escape Patterns  ////////
        ////////////////////////////////////

        vector<EscapePattern> chain_binomial;
        matrix<EscapePattern> escape_patterns;
    } lut;

    /////////////////////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////  Transmission Generations  ///////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////

    /// All possible ranked offspring generations for a given number of infected individuals n > 0.
    /// @example
    ///	n_infec = 3
    ///	return [[0, 0, 0], [0, 0, 1], [0, 1, 1], [0, 1, 2]]
    /// UPDATE: this must be updated to be cache friendly.
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

    void clear_lookuptables()
    {
        lut.num_generations.clear();
        lut.multiplicity_generation.clear();
        lut.generations.clear();
        lut.partition_trees.clear();
    }

    /// Ensure count-compressed generations are available up to n infected.
    void ensure_homogeneous(MemberIndex n)
    {
        if (n < lut.num_generations.size())
            return;

        unsigned old_size(lut.num_generations.size());
        lut.num_generations.resize(n + 1);
        lut.multiplicity_generation.resize(n + 1);
        lut.generations.resize(n + 1);

        if (old_size == 0)
        {
            lut.num_generations[0] = 0;
            lut.multiplicity_generation[0].push_back(0);
            lut.generations[0].push_back({});
            old_size = 1;
        }

        vector<unsigned long> aux_num_vec;
        vector<gentype> aux_gen_vec;

        for (unsigned k = old_size; k <= n; ++k)
        {
            aux_num_vec.clear();
            aux_gen_vec.clear();

            for (gentype g : recursive_span_generation_(k))
            {
                std::map<MemberIndex, MemberIndex> freq_map = {};
                gentype counts = {};
                unsigned long multiplicities = factorial(k);

                for (MemberIndex val : g)
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
            lut.num_generations[k] = tot_num_gen;
            lut.multiplicity_generation[k].resize(tot_num_gen);
            lut.generations[k].resize(tot_num_gen);

            for (size_t c = 0; c < tot_num_gen; ++c)
            {
                lut.multiplicity_generation[k][c] = aux_num_vec[c];
                lut.generations[k][c] = aux_gen_vec[c];
            }
        }
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////  Homogeneous Chain-Binomial LogLikelihood  //////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////

    /// Homogeneous chain binomial log-likelihood.
    /// Based on @cite longini_household_1982
    double log_lik_homogeneous_(MemberIndex n_infec,
                                double commu_infec_rate,
                                double p2p_infec_rate)
    {
        const double log_commu_esc_proba = -commu_infec_rate;
        const double log_commu_infec_proba = std::log1p(-std::exp(log_commu_esc_proba));
        const double log_p2p_esc_proba = -p2p_infec_rate;

        ///////////////////////////////////////////////////
        /////////  (+) Sero-Positive Contribution  ////////
        ///////////////////////////////////////////////////

        // sum over all possible generations
        double sero_positive_term = 0., tmp = 0.;
        for (unsigned long c = 0; c < lut.num_generations[n_infec]; ++c)
        {
            // reference to the LUT
            const auto &gen_counts = lut.generations[n_infec][c];

            // generation (0) : infected by community
            tmp = gen_counts[0] * log_commu_infec_proba;
            // generation (0) : escaped by community
            tmp += (n_infec - gen_counts[0]) * log_commu_esc_proba;

            // generation (i) : within house, i > 0
            for (size_t i = 1; i < gen_counts.size(); ++i)
            {
                // generation (i) : infected by previous generation i - 1
                tmp += gen_counts[i] * std::log1p(-std::exp(gen_counts[i - 1] * log_p2p_esc_proba));

                // generation (i) : escaped from other generations g, i. e. 0 <= g < i - 1
                for (size_t j = 0; j < i - 1; ++j)
                    tmp += gen_counts[i] * gen_counts[j] * log_p2p_esc_proba;
            }

            // times n. of equivalent configurations returning the same kind of generation
            sero_positive_term += lut.multiplicity_generation[n_infec][c] * std::exp(tmp);
        }

        return std::log(sero_positive_term);
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////  COMBINATIONS Recursive Span  /////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////

    /// Algorithm to generate all possible combinations of N elements in k positions.
    /// The N - k complements are appended at the end. Elements indexed from 1 to N.
    /// Based on @cite jimenez_algorithms_1998
    void recursive_span_combinations_(unsigned N,
                                      unsigned k,
                                      vector<unsigned> &current,
                                      matrix<unsigned> &results,
                                      unsigned start,
                                      unsigned loop)
    {
        if (loop + 1 > k)
        {
            // complement of the combination
            // mark elements present in the list
            vector<bool> missing(N, true);
            for (unsigned i = 0; i < k; ++i)
                missing[current[i]] = false;

            // store complement
            for (unsigned i = 0; i < N; ++i)
                if (missing[i])
                    current[loop++] = i;

            // append current combination & complement to results
            results.push_back(current);
        }
        else
            for (unsigned i = start; i < N; ++i)
            {
                current[loop] = i;
                recursive_span_combinations_(N, k, current, results, i + 1, loop + 1);
            }
    }

    /// Store all combinations + complement of N elements in k positions.
    void store_combinations_(unsigned N, unsigned k)
    {
        // generate combinations with the recursive algorithm
        matrix<unsigned> results;
        vector<unsigned> current(N);
        recursive_span_combinations_(N, k, current, results, 0, 0);

        unsigned long num_combin = binom_coeff(N, k);
        lut.num_combinations[N][k] = num_combin;
        lut.combinations[N][k].resize(num_combin);

        // store all results inside the LookUp Table
        unsigned long c = 0;
        for (const vector<unsigned> &res_c : results)
        {
            lut.combinations[N][k][c].resize(N);

            // store combination + complement
            for (unsigned i = 0; i < N; ++i)
                lut.combinations[N][k][c][i] = res_c[i];
            ++c;
        }
    }

    /// Initialize the Look Up Table of combinations and their counts.
    void init_combinations_lookuptable_(unsigned max_N)
    {
        // Assure that LUT are empty before initializing
        lut.num_combinations.clear();
        lut.combinations.clear();

        lut.num_combinations.resize(max_N + 1);
        lut.combinations.resize(max_N + 1);

        for (unsigned N = 0; N <= max_N; ++N)
        {
            lut.num_combinations[N].resize(max_N);
            lut.combinations[N].resize(max_N);

            for (unsigned k = 0; k < N; ++k)
                store_combinations_(N, k);
        }
    }

    /// Clear the Look Up Table of combinations and their counts
    void clear_combinations_lookuptable_()
    {
        lut.num_combinations.clear();
        lut.num_combinations.shrink_to_fit();
        lut.combinations.clear();
        lut.combinations.shrink_to_fit();
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////  Escape Patterns Recursive Functions  /////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////

    /// String with the formula of the likelihood.
    std::string EscapePattern::str() const
    {
        std::ostringstream oss;

        for (size_t c = 0; c < n_entries_; ++c)
        {
            unsigned u = entries_[c];
            if (u < max_n_infec_)
                oss << "B" << u + 1;
            else
                oss << "Q" << u / max_n_infec_ << "_" << u % max_n_infec_ + 1;

            if (c < n_entries_ - 1)
                oss << "*";
        }

        if (n_patterns_ > 0)
        {
            if (n_entries_ > 0)
                oss << "*(";
            oss << "1";
            for (size_t c = 0; c < n_patterns_; ++c)
                oss << " - " << lut.escape_patterns[max_n_infec_][pattern_idxs_[c]].str();
            if (n_entries_ > 0)
                oss << ")";
        }

        return oss.str();
    }

    /// Recursive function to generate all possible escape patterns for a given number of seropositives N_star.
    void recursive_escape_patterns_function_(unsigned N_star,
                                             unsigned N,
                                             vector<unsigned> access_matrix,
                                             EscapePattern &father,
                                             vector<EscapePattern> &results)
    {
        for (unsigned k = 0; k < N; ++k)
            for (unsigned long c = 0; c < lut.num_combinations[N][k]; ++c)
            {
                EscapePattern new_child(N_star);

                // associate entries for this pattern
                vector<unsigned> &indidx = lut.combinations[N][k][c];

                for (unsigned i = k; i < N; ++i)
                {
                    //  community
                    new_child.add_entry(access_matrix[indidx[i]]);

                    // person-to-person
                    for (unsigned j = 0; j < k; ++j)
                        new_child.add_entry(access_matrix[indidx[i] + N * (indidx[j] + 1)]);
                }

                if (k > 0)
                {
                    // find the submatrix of the given infectious
                    // note: always at least 1 infector (community)
                    vector<unsigned> sub_access_matrix((k + 1) * k);

                    // community (inf=0) always at row=0
                    for (unsigned row = 0, inf = 0; row <= k; inf = 1 + indidx[row++])
                        for (unsigned col = 0; col < k; ++col)
                            sub_access_matrix[col + k * row] = access_matrix[indidx[col] + N * inf];

                    recursive_escape_patterns_function_(N_star, k, sub_access_matrix, new_child, results);
                    sub_access_matrix.clear();
                }

                father.add_pattern(results.size()); // attach pattern index
                results.push_back(new_child);       // store pattern
            }
    }

    /// Initialize the Look Up Table of escape patterns and the access object.
    void init_escape_patterns_lookuptable_(unsigned max_n_infec)
    {
        // lookuptable with all combinations + complements
        init_combinations_lookuptable_(max_n_infec);

        // assure that LUT are empty before initializing
        lut.chain_binomial.clear();
        lut.escape_patterns.clear();

        // lookuptable with escape patterns linked to chain-binomial
        // case N_star == 0, for consistency
        lut.chain_binomial.emplace_back(0);
        lut.escape_patterns.push_back(vector<EscapePattern>());

        // case N_star > 0
        for (unsigned n_infec = 1; n_infec <= max_n_infec; ++n_infec)
        {
            vector<unsigned> access_matrix((n_infec + 1) * n_infec);
            std::iota(access_matrix.begin(), access_matrix.end(), 0);

            // store to lookup table
            lut.chain_binomial.emplace_back(n_infec);
            lut.escape_patterns.emplace_back(vector<EscapePattern>());

            recursive_escape_patterns_function_(n_infec, n_infec, access_matrix,
                                                lut.chain_binomial[n_infec],
                                                lut.escape_patterns[n_infec]);
        }
    }

    std::string aux_debugger(unsigned max_n_infec)
    {
        unsigned limit_max_n_infec = 8;
        if (max_n_infec < limit_max_n_infec)
        {
            init_escape_patterns_lookuptable_(max_n_infec);
            std::string output = lut.chain_binomial[max_n_infec].str();
            clear_combinations_lookuptable_();
            return output;
        }
        msg::stop(
            "Maximum number of infected members is fixed to " +
            std::to_string(limit_max_n_infec) + ".");
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////  Heterogeneous Chain-Binomial LogLikelihood  /////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////

    /// Heterogeneous final-size using subset bitmasks.
    /// Based on @cite addy_generalized_1991
    double log_lik_heterogeneous_(MemberIndex n_infec,
                                  const std::vector<double> &commu_infec_rate,
                                  const std::vector<std::vector<double>> &p2p_infec_rate)
    {
        const size_t Nsets = size_t(1) << n_infec; // 2^n_infec = total number of subsets R

        // pre-compute incrementally all values of log-escape probabilities q
        // for all subsets R
        thread_local vector<double> log_q;
        log_q.assign(Nsets * n_infec, 0.);

        // case R = 0 (empty): only community escape [O(n)]
        for (size_t j = 0; j < n_infec; ++j)
            log_q[j] = -commu_infec_rate[j];

        // case R > 0 (set with at least 1 infected) [O(2^n * n)]
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
                    if ((SmR >> j) & size_t(1)) // if j in S \ R?
                        expo += log_q[R * n_infec + j];
                sum += std::exp(log_pi[R] + expo);
                if (R == 0) // last subset
                    break;
            }
            log_pi[S] = (sum < 1.0) ? std::log1p(-sum) : -inf::DOUBLE;
        }

        return log_pi[Nsets - 1];
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////  PartitionTree : build and evaluate  //////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////

    /// Build the partition tree for the infected partition held in c_root[0, K).
    PartitionTree build_partition_tree_(const MemberIndex *c_root,
                                        MemberIndex K)
    {
        size_t total = 1;
        for (MemberIndex s = 0; s < K; ++s)
            total *= (c_root[s] + 1);

        // enumerate all v <= c_root via odometer increment
        vector<gentype> sub_vecs;
        sub_vecs.reserve(total);
        {
            gentype cur(K, 0);
            while (true)
            {
                sub_vecs.push_back(cur);
                MemberIndex s = 0;
                while (s < K && cur[s] == c_root[s])
                    cur[s++] = 0;
                if (s == K)
                    break;
                ++cur[s];
            }
        }

        // linear extension of the componentwise order:
        // |v|_1 ascending, lex as tiebreak
        std::sort(sub_vecs.begin(), sub_vecs.end(),
                  [](const gentype &a, const gentype &b)
                  {
                      unsigned sa = 0, sb = 0;
                      for (MemberIndex x : a)
                          sa += x;
                      for (MemberIndex x : b)
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
                    for (MemberIndex s = 0; s < K; ++s)
                        log_C += std::lgamma(t[s] + 1) - std::lgamma(w[s] + 1) -
                                 std::lgamma(t[s] - w[s] + 1);
                    dp.ancestors[ti].emplace_back(idx_of[w], log_C);
                }
                MemberIndex s = 0;
                while (s < K && w[s] == t[s])
                    w[s++] = 0;
                if (s == K)
                    break;
                ++w[s];
            }
        }

        return dp;
    }

    void ensure_partitioned(HouseIndex idx,
                            const std::vector<MemberIndex> &part_counts,
                            MemberIndex n_parts)
    {
        if (lut.partition_trees.size() <= idx)
            lut.partition_trees.resize(idx + 1);
        lut.partition_trees[idx] = build_partition_tree_(part_counts.data(), n_parts);
    }

    /// Evaluate log P(final = c_root | all infecteds in c_root).
    /// Ball recursion: pi(t|t) = 1 - sum_{w < t} C(t,w) * pi(w|w) * prod_s q_s(w)^(t_s - w_s).
    double evaluate_partition_tree_(const PartitionTree &dp,
                                    const vector<double> &log_commu_esc,
                                    const matrix<double> &log_p2p_esc)
    {
        const MemberIndex K = dp.sub_vecs.empty()
                                  ? 0
                                  : static_cast<MemberIndex>(dp.sub_vecs.front().size());

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
                for (MemberIndex s = 0; s < K; ++s)
                {
                    if (t[s] == w[s])
                        continue;
                    double log_q_s = log_commu_esc[s];
                    for (MemberIndex u = 0; u < K; ++u)
                        log_q_s += w[u] * log_p2p_esc[u][s];
                    expo += (t[s] - w[s]) * log_q_s;
                }

                sum += std::exp(log_C + log_pi[wi] + expo);
            }
            log_pi[ti] = (sum < 1.0) ? std::log1p(-sum) : -inf::DOUBLE;
        }

        return log_pi.back();
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////  Partitioned Chain-Binomial LogLikelihood  ////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////

    /// Grouped chain binomial log-likelihood.
    ///
    /// When multiple seropositives share the same contribution,
    /// use partition sequences per generation instead of all individual permutations.
    /// NOTE: I don't think I can use bitmasks tricks here because of
    /// non analytical solution for number of integer partitions
    ///
    /// @param n_parts number of parts
    /// @param part_indexes indexes of household member representative of part
    /// @param commu_infec_rate per-individual community infection rates
    /// @param p2p_infec_rate per-individual p2p infection rates
    double log_lik_partitioned_(HouseIndex dp_idx,
                                MemberIndex n_parts,
                                const std::vector<MemberIndex> &part_indexes,
                                const std::vector<double> &commu_infec_rate,
                                const std::vector<std::vector<double>> &p2p_infec_rate)
    {
        // pre-compute log-escape probabilities per type
        thread_local vector<double> log_commu_esc;
        log_commu_esc.assign(n_parts, 0.);
        thread_local matrix<double> log_p2p_esc;
        log_p2p_esc.assign(n_parts, vector<double>(n_parts, 0.));

        for (MemberIndex t = 0; t < n_parts; ++t)
        {
            log_commu_esc[t] = -commu_infec_rate[part_indexes[t]];
            for (MemberIndex s = 0; s < n_parts; ++s)
                log_p2p_esc[t][s] = -p2p_infec_rate[part_indexes[t]][part_indexes[s]];
        }

        return evaluate_partition_tree_(lut.partition_trees[dp_idx],
                                        log_commu_esc, log_p2p_esc);
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////////
    //////////////////////////  Chain-Binomial Log-Likelihood (Dispatcher)  ///////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////

    double log_likelihood(HouseIndex dp_idx,
                          MemberIndex n_infec, MemberIndex n_parts,
                          const std::vector<MemberIndex> &part_indexes,
                          const std::vector<double> &commu_infec_rate,
                          const std::vector<std::vector<double>> &p2p_infec_rate)
    {
        ///////////////////////////////////////////////////
        /////////  (-) Sero-Negative Contribution  ////////
        ///////////////////////////////////////////////////

        double ll = 0.;
        for (size_t i = n_infec; i < commu_infec_rate.size(); ++i)
        {
            // community : i escaped from community
            ll -= commu_infec_rate[i];
            // person-to-person : i escaped from k
            for (size_t k = 0; k < n_infec; ++k)
                ll -= p2p_infec_rate[k][i];
        }

        ///////////////////////////////////////////////////
        /////////  (+) Sero-Positive Contribution  ////////
        ///////////////////////////////////////////////////

        // no infected
        if (n_infec == 0)
            return ll;
        // single seropositive
        else if (n_infec == 1)
            return ll + std::log1p(-std::exp(-commu_infec_rate[0]));
        // homogeneous case
        else if (n_parts < 2)
            return ll + log_lik_homogeneous_(
                            n_infec, commu_infec_rate[0], p2p_infec_rate[0][1]);
        // partitioned case
        else if (n_parts < n_infec)
            return ll + log_lik_partitioned_(
                            dp_idx, n_parts, part_indexes,
                            commu_infec_rate, p2p_infec_rate);
        // heterogeneous case
        else
            return ll + log_lik_heterogeneous_(
                            n_infec, commu_infec_rate, p2p_infec_rate);
    }
}
