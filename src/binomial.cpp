/* housepi - Epidemiological Analysis of Household Outbreaks.
 *
 * Copyright © 2026 Institut Pasteur (Paris), CNRS and INSERM.
 * Distributed under the terms of the GNU General Public License (GPLv3).
 * See the LICENSE file for details.
 */

// #include <string>
// #include <sstream>
#include <algorithm>
#include <cmath>
#include <functional> // UPDATE: drop this dependency
#include <map>

#include "tools/messages.h"
#include "tools/operators.h" // UPDATE: drop this dependency
#include "binomial.h"

/////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////// *** BINOMIAL AUXILIARY NAMESPACE *** /////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////

namespace housepi::bino
{
    using std::vector;
    using gentype = vector<unsigned>;

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
    /////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////// *** EscapePattern : Class *** /////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////

    class EscapePattern
    {
    public:
        EscapePattern(unsigned N_star)
            : N_star_(N_star) {}
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
        /// Compute the likelihood starting from the current escape pattern.
        ///	NOTE: it requires initialization of the lookup table.
        double likelihood(const vector<double> &esc_mtx) const;

        /*
        /// String with the formula of the likelihood.
        const std::string EscapePattern::str();
        */

    private:
        const unsigned N_star_;
        size_t n_entries_ = 0;
        size_t n_patterns_ = 0;

        vector<unsigned> entries_ = {};
        vector<unsigned long> pattern_idxs_ = {};
    };

    /////////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////  *** GroupedTree : Class ***  /////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////

    struct GroupedTree
    {
        struct Node
        {
            double log_weight;          // pre-computed lgamma multinomial coefficient
            unsigned gen;               // generation index
            vector<unsigned> comp;      // composition at this node (size K)
            vector<unsigned> prev;      // prev_by_type (size K)
            vector<unsigned> cumul;     // cumul_by_type (size K)
            vector<unsigned> remaining; // remaining before subtracting comp (size K)
            vector<size_t> children;    // indices into flat node array
        };

        vector<Node> nodes;
        vector<size_t> roots; // top-level node indices

        double evaluate_node_(
            size_t idx,
            const vector<double> &log_commu_esc,
            const matrix<double> &log_p2p_esc) const;

        double evaluate(
            const vector<double> &log_commu_esc,
            const matrix<double> &log_p2p_esc) const;
    };

    /////////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////  LookUpTables  ////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////

    struct LookUpTables
    {
        /////////////////////////////////////////////
        /////////  Transmission Generations  ////////
        /////////////////////////////////////////////

        vector<unsigned long> num_generations;         /// n. of generations
        matrix<unsigned long> multiplicity_generation; /// n. of multiplicites per generation
        matrix<gentype> generations;                   /// generation vectors

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

        /////////////////////////////////////////
        /////////  Cached Grouped Trees  ////////
        /////////////////////////////////////////

        /// Grouped chain-binomial cached trees
        /// Key: (n_infec, type_counts); Value: one tree per generation pattern
        std::map<vector<unsigned>, vector<GroupedTree>> grouped_patterns;
    } lut;

    /////////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////  TRANSMISSION GENERATIONS  ///////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////
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
                    g[i] += 1;
                output.push_back(gentype(root, 0) + g);
            }
        return output;
    }

    /// Initialize the Look Up Table of count-compressed generations.
    /// UPDATE: Avoid regenerating if already ok...
    void init_generations_lookuptable(unsigned max_N)
    {
        // Assure that LUT are empty before initializing
        lut.num_generations.clear();
        lut.multiplicity_generation.clear();
        lut.generations.clear();
        lut.grouped_patterns.clear();

        // Allocate LUT memory
        lut.num_generations.resize(max_N + 1);
        lut.multiplicity_generation.resize(max_N + 1);
        lut.generations.resize(max_N + 1);

        // Empty generation n = 0
        // note : just for access clarity
        lut.num_generations[0] = 0;
        lut.multiplicity_generation[0].push_back(0);
        lut.generations[0].push_back({});

        // Count-compression of all possible ranked offpsring generations.
        size_t tot_num_gen = 0;
        vector<unsigned long> aux_num_vec;
        vector<gentype> aux_gen_vec;

        for (unsigned n = 1; n <= max_N; ++n)
        {
            aux_num_vec.clear();
            aux_gen_vec.clear();

            for (gentype g : recursive_span_generation_(n))
            {
                std::map<unsigned, unsigned> freq_map = {};
                gentype counts = {};
                unsigned long multiplicities = factorial(n);

                // I take advantage of the order granted by the map
                for (unsigned val : g)
                    if (freq_map.count(val))
                        ++freq_map[val];
                    else
                        freq_map.emplace(val, 1);

                // assign counts
                for (const auto &[val, c] : freq_map)
                {
                    counts.push_back(c);
                    multiplicities /= factorial(c);
                }

                // add to the list
                aux_num_vec.push_back(multiplicities);
                aux_gen_vec.push_back(counts);
            }

            // allocate lookup table
            // note : how to obtain tot_num_gen a priori
            tot_num_gen = aux_gen_vec.size();
            lut.num_generations[n] = tot_num_gen;
            lut.multiplicity_generation[n].resize(tot_num_gen);
            lut.generations[n].resize(tot_num_gen);

            // copy from auxiliary vectors
            for (size_t c = 0; c < tot_num_gen; ++c)
            {
                lut.multiplicity_generation[n][c] = aux_num_vec[c];
                lut.generations[n][c] = aux_gen_vec[c];
            }
        }
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////  HOMOGENEOUS Chain-Binomisl LogLikelihood  //////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////

    /// Homogeneous chain binomial log-likelihood.
    double chain_binom_homogeneous_log_lik(
        unsigned n_infec, unsigned house_size,
        double commu_escape_rate, double p2p_escape_rate)

    {
        const double log_commu_esc_proba = -commu_escape_rate;

        if (n_infec == 0)
            return house_size * log_commu_esc_proba;

        const double log_commu_infec_proba = std::log1p(-std::exp(log_commu_esc_proba));
        const double log_p2p_esc_proba = -p2p_escape_rate;

        ///////////////////////////////////////////////////
        /////////  (-) Sero-Negative Contribution  ////////
        ///////////////////////////////////////////////////

        const double log_sero_neg_term = (house_size - n_infec) *
                                         (log_commu_esc_proba + n_infec * log_p2p_esc_proba);

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

        return log_sero_neg_term + std::log(sero_positive_term);
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////  COMBINATIONS Recursive Span  /////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////

    /// Algorithm to generate all possible combinations of N elements in k positions, plus the N - k complements.
    /// The elements are indexed from 1 to N.
    ///
    /// Ref: (modified from)
    /// ----
    /// Daniel Angel Jimenez
    /// https://www.cs.utexas.edu/users/djimenez/utsa/cs3343/lecture25.html
    void recursive_span_combinations_(unsigned N, unsigned k, vector<unsigned> &current,
                                      matrix<unsigned> &results, unsigned start, unsigned loop)
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
        for (const gentype &res_c : results)
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
    ///////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////  ESCAPE PATTERNS Recursive Functions  /////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////

    /*
    /// String with the formula of the likelihood.
    const std::string EscapePattern::str()
    {
        std::ostringstream oss;

        for (size_t c = 0; c < n_entries_; ++c)
        {
            unsigned u = entries_[c];
            if (u < N_star_)
                oss << "B" << u + 1;
            else
                oss << "Q" << u / N_star_ << "_" << u % N_star_ + 1;

            if (c < n_entries_ - 1)
                oss << "*";
        }

        if (n_patterns_ > 0)
        {
            if (n_entries_ > 0)
                oss << "*(";
            oss << "1";
            for (size_t c = 0; c < n_patterns_; ++c)
                oss << " - " << lut.EscapePatterns[N_star_][pattern_idxs_[c]].str();
            if (n_entries_ > 0)
                oss << ")";
        }

        return oss.str();
    }
    */

    double EscapePattern::likelihood(const vector<double> &esc_mtx) const
    {
        double output = 1.;

        for (size_t c = 0; c < n_patterns_; ++c)
            output -= lut.escape_patterns[N_star_][pattern_idxs_[c]].likelihood(esc_mtx);

        for (size_t c = 0; c < n_entries_; ++c)
            output *= esc_mtx[entries_[c]];

        return output;
    }

    /// Recursive function to generate all possible escape patterns for a given number of seropositives N_star.
    void recursive_escape_patterns_function_(unsigned N_star, unsigned N, vector<unsigned> access_matrix,
                                             EscapePattern &father, vector<EscapePattern> &results)
    {
        for (unsigned k = 0; k < N; ++k)
            for (unsigned long c = 0; c < lut.num_combinations[N][k]; ++c)
            {
                EscapePattern new_child(N_star);

                // associate entries for this pattern
                gentype &indidx = lut.combinations[N][k][c];

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
    /// UPDATE: Avoid regenerating if already ok...
    void init_escape_patterns_lookuptable(unsigned max_N, bool clean_combinations)
    {
        // lookuptable with all combinations + complements
        init_combinations_lookuptable_(max_N);

        // assure that LUT are empty before initializing
        lut.chain_binomial.clear();
        lut.escape_patterns.clear();

        // lookuptable with escape patterns linked to chain-binomial
        // case N_star == 0, for consistency
        lut.chain_binomial.emplace_back(0);
        lut.escape_patterns.push_back(vector<EscapePattern>());

        // case N_star > 0
        for (unsigned N_star = 1; N_star <= max_N; ++N_star)
        {
            vector<unsigned> access_matrix((N_star + 1) * N_star);
            for (unsigned i = 0; i < (N_star + 1) * N_star; ++i)
                access_matrix[i] = i;

            // store to lookup table
            lut.chain_binomial.emplace_back(N_star);
            lut.escape_patterns.emplace_back(vector<EscapePattern>());

            recursive_escape_patterns_function_(N_star, N_star, access_matrix,
                                                lut.chain_binomial[N_star],
                                                lut.escape_patterns[N_star]);
        }

        // deallocate memory of combinations lookup table
        if (clean_combinations)
            clear_combinations_lookuptable_();
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////  HETEROGENEOUS Chain-Binomisl LogLikelihood  /////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////

    /*
    /// UPDATE: to be moved to a separate file.
    /// Debug mode, to check values of the likelihood and print formulation.
    ///
    /// g++ -Wall -Wextra -std=gnu++11 -g -O2 -o auxiliary_binomial.o auxiliary_binomial.cpp ./auxiliary_binomial.o 7
    ///
    /// Auxiliary function to perform debugs
    void auxiliary_chain_binomial_debugger(unsigned max_num, unsigned limit)
    {
        using namespace bax;
        for (unsigned N_star = 0; N_star <= max_num; ++N_star)
        {
            size_t n_patterns = lut.escape_patterns[N_star].size();
            std::cerr << "\n<<" << N_star << ">>"
                      << " n.terms = " << (n_patterns == 0 ? 0 : 1 + (n_patterns - 1) / 2) << std::endl;

            // elements for the escape matrix //
            vector<double> commu_esc_rate(N_star);
            for (unsigned i = 0; i < N_star; ++i)
            {
                commu_esc_rate[i] = -std::log(0.15 + 0.025 * i);
                if (N_star < limit)
                    std::cerr << "B" << i + 1 << " = " << std::exp(-commu_esc_rate[i]) << std::endl;
            }
            matrix<double> p2p_esc_rate(N_star);
            for (unsigned row = 0; row < N_star; ++row)
            {
                p2p_esc_rate[row].resize(N_star);
                for (unsigned col = 0; col < N_star; ++col)
                {
                    p2p_esc_rate[row][col] = -std::log(0.25 + 0.0125 * col + 0.0125 * row);
                    if (N_star < limit)
                        std::cerr << "Q" << row + 1 << "_" << col + 1 << " = "
                                  << std::exp(-p2p_esc_proba[row][col]) << std::endl;
                }
            }

            if (N_star < limit)
                std::cerr << lut.chain_binomial[N_star].str() << std::endl;

            auto start = std::chrono::high_resolution_clock::now();
            std::cerr << "likelihood: " << std::exp(chain_binom_heterogeneous_log_lik(N_star, commu_esc_rate, p2p_esc_rate));
            auto end = std::chrono::high_resolution_clock::now();
            std::chrono::duration<double> elapsed = end - start;
            std::cerr << " | time: " << elapsed.count() << " seconds" << std::endl;
        }

        std::cerr << std::endl;
    }
    */

    /// Heterogeneous chain binomial log-likelihood for ordered probabilities (seropositives, seronegatives).
    double chain_binom_heterogeneous_log_lik(
        unsigned n_infec,
        const std::vector<double> &commu_escape_rate,
        const std::vector<std::vector<double>> &p2p_escape_rate)
    {
        ///////////////////////////////////////////////////
        /////////  (-) Sero-Negative Contribution  ////////
        ///////////////////////////////////////////////////

        double log_neg_term = 0.;
        for (size_t i = n_infec; i < commu_escape_rate.size(); ++i)
        {
            // community : i escaped from community
            log_neg_term -= commu_escape_rate[i];
            // person-to-person : i escaped from k
            for (unsigned k = 0; k < n_infec; ++k)
                log_neg_term -= p2p_escape_rate[k][i];
        }

        if (n_infec == 0)
            return log_neg_term;

        ///////////////////////////////////////////////////
        /////////  (+) Sero-Positive Contribution  ////////
        ///////////////////////////////////////////////////

        // single seropositive
        if (n_infec == 1)
            return log_neg_term + std::log1p(-std::exp(-commu_escape_rate[0]));

        // flatten sub-matrix with only seropositives
        // UPDATE: can I avoid this re-allocation?
        thread_local vector<double> escape_matrix;
        escape_matrix.assign((n_infec + 1) * n_infec, 1.);

        for (unsigned i = 0; i < n_infec; ++i)
        {
            // community : i escaped from community
            escape_matrix[i] = std::exp(-commu_escape_rate[i]);

            // person-to-person : i escaped from k
            for (unsigned k = 0; k < n_infec; ++k)
                escape_matrix[n_infec * (k + 1) + i] = std::exp(-p2p_escape_rate[k][i]);

            // infector always "escapes" himself (check if actually relevant)
            escape_matrix[n_infec * (i + 1) + i] = 1.;
        }

        return log_neg_term + std::log(lut.chain_binomial[n_infec].likelihood(escape_matrix));
    }

    /////////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////  GroupedTree : implementations  ////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////

    /// Recursive enumeration of partitions of n into K parts, where each
    /// part t contributes at most remaining[tau] individuals.
    /// Calls callback(composition) for each valid composition.
    void enumerate_partitions_(
        unsigned type_idx,
        unsigned remaining_n,
        const vector<unsigned> &remaining,
        vector<unsigned> &current,
        unsigned K,
        const std::function<void(const vector<unsigned> &)> &callback)
    {
        if (type_idx == K - 1)
        {
            if (remaining_n <= remaining[type_idx])
            {
                current[type_idx] = remaining_n;
                callback(current);
            }
            return;
        }
        unsigned max_c = std::min(remaining_n, remaining[type_idx]);
        for (unsigned c = 0; c <= max_c; ++c)
        {
            current[type_idx] = c;
            enumerate_partitions_(type_idx + 1, remaining_n - c,
                                  remaining, current, K, callback);
        }
    }

    double GroupedTree::evaluate_node_(
        size_t idx,
        const vector<double> &log_commu_esc,
        const matrix<double> &log_p2p_esc) const
    {
        const auto &node = nodes[idx];
        const unsigned K = node.comp.size();

        // compute log probability from cached integer coefficients + runtime escape rates
        double log_prob = 0.;

        if (node.gen == 0)
            for (unsigned t = 0; t < K; ++t)
            {
                if (node.comp[t] > 0)
                    log_prob += node.comp[t] * std::log1p(-std::exp(log_commu_esc[t]));
                if (node.remaining[t] > node.comp[t])
                    log_prob += (node.remaining[t] - node.comp[t]) * log_commu_esc[t];
            }
        else
            for (unsigned t = 0; t < K; ++t)
            {
                if (node.comp[t] == 0)
                    continue;

                double log_esc_prev = 0.;
                for (unsigned s = 0; s < K; ++s)
                {
                    log_prob += node.comp[t] * node.cumul[s] * log_p2p_esc[s][t];
                    log_esc_prev += node.prev[s] * log_p2p_esc[s][t];
                }

                log_prob += node.comp[t] * std::log1p(-std::exp(log_esc_prev));
            }

        // evaluate children (next generation)
        double child_result = 0.;
        if (node.children.empty())
            child_result = 1.0;
        else
            for (size_t ci : node.children)
                child_result += evaluate_node_(ci, log_commu_esc, log_p2p_esc);

        return std::exp(node.log_weight + log_prob) * child_result;
    }

    double GroupedTree::evaluate(
        const vector<double> &log_commu_esc,
        const matrix<double> &log_p2p_esc) const
    {
        double sum = 0.;
        for (size_t r : roots)
            sum += evaluate_node_(r, log_commu_esc, log_p2p_esc);
        return sum;
    }

    /// Recursively build the cached tree for one generation pattern.
    void build_grouped_tree_recursive_(
        unsigned g,
        unsigned num_gens,
        const gentype &gen_counts,
        vector<unsigned> &remaining,
        const vector<unsigned> &prev_by_type,
        const vector<unsigned> &cumul_by_type,
        unsigned K,
        GroupedTree &tree,
        vector<size_t> &parent_children)
    {
        if (g == num_gens)
            return;

        unsigned n_g = gen_counts[g];

        auto callback = [&](const vector<unsigned> &comp)
        {
            GroupedTree::Node node;
            node.gen = g;
            node.comp = comp;
            node.prev = prev_by_type;
            node.cumul = cumul_by_type;
            node.remaining.assign(remaining.begin(), remaining.end());

            // Pre-compute log_weight (multinomial coefficient)
            node.log_weight = 0.;
            for (unsigned t = 0; t < K; ++t)
                node.log_weight += std::lgamma(remaining[t] + 1) - std::lgamma(comp[t] + 1) -
                                   std::lgamma(remaining[t] - comp[t] + 1);

            size_t my_idx = tree.nodes.size();
            tree.nodes.push_back(std::move(node));
            parent_children.push_back(my_idx);

            // Recurse into next generation
            // NOTE: collect children into a local vector because tree.nodes
            // may reallocate during recursion, invalidating references.
            vector<unsigned> new_remaining(K), new_cumul(K);
            for (unsigned t = 0; t < K; ++t)
            {
                new_remaining[t] = remaining[t] - comp[t];
                new_cumul[t] = cumul_by_type[t] + (g > 0 ? prev_by_type[t] : 0);
            }

            vector<size_t> my_children;
            build_grouped_tree_recursive_(
                g + 1, num_gens, gen_counts,
                new_remaining, comp, new_cumul, K,
                tree, my_children);
            tree.nodes[my_idx].children = std::move(my_children);
        };

        vector<unsigned> comp(K, 0);
        std::function<void(const vector<unsigned> &)> cb = callback;
        enumerate_partitions_(0, n_g, remaining, comp, K, cb);
    }

    /// Build a GroupedTree for one generation pattern.
    GroupedTree build_grouped_tree_(
        const gentype &gen_counts,
        unsigned K,
        const vector<unsigned> &type_counts)
    {
        GroupedTree tree;
        vector<unsigned> remaining(type_counts.begin(), type_counts.end());
        vector<unsigned> prev(K, 0), cumul(K, 0);

        build_grouped_tree_recursive_(
            0, gen_counts.size(), gen_counts,
            remaining, prev, cumul, K,
            tree, tree.roots);

        return tree;
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////  GROUPED Chain-Binomial LogLikelihood  //////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////

    /// Grouped chain binomial log-likelihood.
    ///
    /// When multiple seropositives share the same contribution,
    /// use partition sequences per generation instead of all individual permutations.
    ///
    /// @param n_infec number of infected members
    /// @param n_parts number of parts
    /// @param part_counts counts of in-house seropositives parts
    /// @param part_indexes indexes of household member representative of part
    /// @param commu_esc_rate per-individual community escape rates
    /// @param p2p_esc_rate per-individual p2p escape rates
    double chain_binom_grouped_log_lik(
        size_t n_infec,
        size_t n_parts,
        const std::vector<size_t> &part_counts,
        const std::vector<size_t> &part_indexes,
        const std::vector<double> &commu_esc_rate,
        const std::vector<std::vector<double>> &p2p_esc_rate)
    {
        ///////////////////////////////////////////////////
        /////////  (-) Sero-Negative Contribution  ////////
        ///////////////////////////////////////////////////

        double log_neg_term = 0.;
        for (size_t i = n_infec; i < commu_esc_rate.size(); ++i)
        {
            log_neg_term -= commu_esc_rate[i];
            for (unsigned k = 0; k < n_infec; ++k)
                log_neg_term -= p2p_esc_rate[k][i];
        }

        if (n_infec == 0)
            return log_neg_term;

        ///////////////////////////////////////////////////
        /////////  (+) Sero-Positive Contribution  ////////
        ///////////////////////////////////////////////////

        // single seropositive
        if (n_infec == 1)
            return log_neg_term + std::log1p(-std::exp(-commu_esc_rate[0]));

        // pre-compute log-escape probabilities per type
        thread_local vector<double> log_commu_esc;
        log_commu_esc.assign(n_parts, 0.);
        thread_local matrix<double> log_p2p_esc;
        log_p2p_esc.assign(n_parts, vector<double>(n_parts, 0.));

        for (unsigned t = 0; t < n_parts; ++t)
        {
            log_commu_esc[t] = -commu_esc_rate[part_indexes[t]];
            for (unsigned s = 0; s < n_parts; ++s)
                log_p2p_esc[t][s] = -p2p_esc_rate[part_indexes[t]][part_indexes[s]];
        }

        // lookup or build cached trees for this (n_infec, type_counts)
        // UPDATE: I am keeping this just to cast size_t to unsigned
        vector<unsigned> tc_key(part_counts.begin(), part_counts.begin() + n_parts);

        auto it = lut.grouped_patterns.find(tc_key);
        if (it == lut.grouped_patterns.end())
        {
            vector<GroupedTree> trees;
            trees.reserve(lut.num_generations[n_infec]);
            for (unsigned long c = 0; c < lut.num_generations[n_infec]; ++c)
                trees.push_back(build_grouped_tree_(
                    lut.generations[n_infec][c], n_parts, tc_key));
            it = lut.grouped_patterns.emplace(tc_key, std::move(trees)).first;
        }

        // sum over generation patterns using cached trees
        double pos_term = 0.;
        for (unsigned long c = 0; c < lut.num_generations[n_infec]; ++c)
            pos_term += it->second[c].evaluate(log_commu_esc, log_p2p_esc);

        return log_neg_term + std::log(pos_term);
    }
}
