/* housepi - Epidemiological Analysis of Household Outbreaks.
 *
 * Copyright © 2026 Institut Pasteur (Paris), CNRS and INSERM.
 * Distributed under the terms of the GNU General Public License (GPLv3).
 * See the LICENSE file for details.
 */

// #include <string>
// #include <sstream>
#include <map>

#include "tools/operators.h" // UPDATE: drop this dependency
#include "binomial.h"

/////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////// *** BINOMIAL AUXILIARY NAMESPACE *** /////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////

namespace bax
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
        if (N == 0 || k == N)
            return 1;
        // note: use symmetry
        if (k > N - k)
            k = N - k;
        unsigned long result = 1;
        for (unsigned i = 1; i <= k; ++i)
            result = result * (N - k + i) / i;
        return result;
    }

    /////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////  LookUp Tables  ////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////

    // Transmission Generations
    vector<unsigned long> Num_Gener_LookupTable;     /// n. of generations
    matrix<unsigned long> Multipl_Gener_LookupTable; /// n. of multiplicites per generation
    matrix<gentype> Generations_LookupTable;         /// generation vectors

    // Combinations
    matrix<unsigned long> Num_Combin_LookupTable;      /// LUT wiht n.e number of combinations
    matrix<matrix<unsigned>> Combinations_LookupTable; /// LUT wiht combinations and complements

    class EscapePattern; // NOTE: forward declaration of EscapePattern

    // EscapePatterns
    vector<EscapePattern> ChainBinomial_LookupTable;
    matrix<EscapePattern> EscapePatterns_LookupTable;

    /////////////////////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////  TRANSMISSION GENERATIONS  ///////////////////////////////////
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
            for (auto g : recursive_span_generation_(n - root))
            {
                for (unsigned i = 0; i < n - root; ++i)
                    g[i] += 1;
                output.push_back(gentype(root, 0) + g);
            }
        return output;
    }

    /// Initialize the Look Up Table of count-compressed generations.
    void init_generations_lookuptable(unsigned max_N)
    {
        // Assure that LUT are empty before initializing
        Num_Gener_LookupTable.clear();
        Multipl_Gener_LookupTable.clear();
        Generations_LookupTable.clear();

        // Allocate LUT memory
        Num_Gener_LookupTable.resize(max_N + 1);
        Multipl_Gener_LookupTable.resize(max_N + 1);
        Generations_LookupTable.resize(max_N + 1);

        // Empty generation n = 0
        // note : just for access clarity
        Num_Gener_LookupTable[0] = 0;
        Multipl_Gener_LookupTable[0].push_back(0);
        Generations_LookupTable[0].push_back({});

        // Count-compression of all possible ranked offpsring generations.
        size_t tot_num_gen = 0;
        vector<unsigned long> aux_num_lut;
        vector<gentype> aux_gen_lut;

        for (unsigned n = 1; n <= max_N; ++n)
        {
            aux_num_lut.clear();
            aux_gen_lut.clear();

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
                aux_num_lut.push_back(multiplicities);
                aux_gen_lut.push_back(counts);
            }

            // allocate lookup table
            // note : how to obtain tot_num_gen a priori
            tot_num_gen = aux_gen_lut.size();
            Num_Gener_LookupTable[n] = tot_num_gen;
            Multipl_Gener_LookupTable[n].resize(tot_num_gen);
            Generations_LookupTable[n].resize(tot_num_gen);

            // copy from auxiliary vectors
            for (size_t c = 0; c < tot_num_gen; ++c)
            {
                Multipl_Gener_LookupTable[n][c] = aux_num_lut[c];
                Generations_LookupTable[n][c] = aux_gen_lut[c];
            }
        }
    }

    unsigned get_num_generations(unsigned n_positive)
    {
        return Num_Gener_LookupTable[n_positive];
    }
    std::vector<unsigned> get_generation_vector(unsigned n_positive, unsigned i)
    {
        return Generations_LookupTable[n_positive][i];
    }
    unsigned get_generation_multipl(unsigned n_positive, unsigned i)
    {
        return Multipl_Gener_LookupTable[n_positive][i];
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////  Recursive Span Combinations  /////////////////////////////////
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
        Num_Combin_LookupTable[N][k] = num_combin;
        Combinations_LookupTable[N][k].resize(num_combin);

        // store all results inside the LookUp Table
        unsigned long c = 0;
        for (const auto &res_c : results)
        {
            Combinations_LookupTable[N][k][c].resize(N);

            // store combination + complement
            for (unsigned i = 0; i < N; ++i)
                Combinations_LookupTable[N][k][c][i] = res_c[i];
            ++c;
        }
    }

    /// Initialize the Look Up Table of combinations and their counts.
    void init_combinations_lookuptable_(unsigned max_N)
    {
        // Assure that LUT are empty before initializing
        Num_Combin_LookupTable.clear();
        Combinations_LookupTable.clear();

        Num_Combin_LookupTable.resize(max_N + 1);
        Combinations_LookupTable.resize(max_N + 1);

        for (unsigned N = 0; N <= max_N; ++N)
        {
            Num_Combin_LookupTable[N].resize(max_N);
            Combinations_LookupTable[N].resize(max_N);

            for (unsigned k = 0; k < N; ++k)
                store_combinations_(N, k);
        }
    }

    /// Clear the Look Up Table of combinations and their counts
    void clear_combinations_lookuptable_()
    {
        Num_Combin_LookupTable.clear();
        Num_Combin_LookupTable.shrink_to_fit();
        Combinations_LookupTable.clear();
        Combinations_LookupTable.shrink_to_fit();
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
        double likelihood(const vector<double> &esc_mtx) const
        {
            double output = 1.;

            for (size_t c = 0; c < n_patterns_; ++c)
                output -= EscapePatterns_LookupTable[N_star_][pattern_idxs_[c]].likelihood(esc_mtx);

            for (size_t c = 0; c < n_entries_; ++c)
                output *= esc_mtx[entries_[c]];

            return output;
        }

        /*
        /// String with the formula of the likelihood.
        const std::string str()
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
                    oss << " - " << EscapePatterns_LookupTable[N_star_][pattern_idxs_[c]].str();
                if (n_entries_ > 0)
                    oss << ")";
            }

            return oss.str();
        }
        */

    private:
        const unsigned N_star_;
        size_t n_entries_ = 0;
        size_t n_patterns_ = 0;

        vector<unsigned> entries_ = {};
        vector<unsigned long> pattern_idxs_ = {};
    };

    ///////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////  Recursive Escape Patterns Functions  /////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////

    /// Recursive function to generate all possible escape patterns for a given number of seropositives N_star.
    void recursive_escape_patterns_function_(unsigned N_star, unsigned N, vector<unsigned> access_matrix,
                                             EscapePattern &father, vector<EscapePattern> &results)
    {
        for (unsigned k = 0; k < N; ++k)
            for (unsigned long c = 0; c < Num_Combin_LookupTable[N][k]; ++c)
            {
                EscapePattern new_child(N_star);

                // associate entries for this pattern
                auto &indidx = Combinations_LookupTable[N][k][c];

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
    void init_escape_patterns_lookuptable(unsigned max_N, bool clean_combinations)
    {
        // lookuptable with all combinations + complements
        init_combinations_lookuptable_(max_N);

        // assure that LUT are empty before initializing
        ChainBinomial_LookupTable.clear();
        EscapePatterns_LookupTable.clear();

        // lookuptable with escape patterns linked to chain-binomial
        // case N_star == 0, for consistency
        ChainBinomial_LookupTable.emplace_back(0);
        EscapePatterns_LookupTable.push_back(vector<EscapePattern>());

        // case N_star > 0
        for (unsigned N_star = 1; N_star <= max_N; ++N_star)
        {
            vector<unsigned> access_matrix((N_star + 1) * N_star);
            for (unsigned i = 0; i < (N_star + 1) * N_star; ++i)
                access_matrix[i] = i;

            // store to lookup table
            ChainBinomial_LookupTable.emplace_back(N_star);
            EscapePatterns_LookupTable.emplace_back(vector<EscapePattern>());

            recursive_escape_patterns_function_(N_star, N_star, access_matrix,
                                                ChainBinomial_LookupTable[N_star],
                                                EscapePatterns_LookupTable[N_star]);
        }

        // deallocate memory of combinations lookup table
        if (clean_combinations)
            clear_combinations_lookuptable_();
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////  Heterogeneous Chain-Binomisl LogLikelihood  /////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////

    /// Heterogeneous chain binomial log-likelihood for ordered probabilities (seropositives, seronegatives).
    double chain_binomial_heter_log_lik(unsigned N_star,
                                        const vector<double> &commu_escape_rate,
                                        const std::vector<std::vector<double>> &p2p_escape_rate)

    {
        unsigned hh_size = commu_escape_rate.size();
        double ll = 0.;

        // (-) seronegatives contribution
        for (unsigned i = N_star; i < hh_size; ++i)
        {
            // community : i escaped from community
            ll -= commu_escape_rate[i];
            // person-to-person : i escaped from k
            for (unsigned k = 0; k < N_star; ++k)
                ll -= p2p_escape_rate[k][i];
        }

        // (+) seropositives contribution
        if (N_star == 1) // single case
            ll += std::log1p(-std::exp(-commu_escape_rate[0]));

        else if (N_star > 1) // multiple cases
        {
            // flatten sub-matrix with only seropositives
            // UPDATE: can I avoid this allocation?
            vector<double> escape_matrix((N_star + 1) * N_star, 1.);

            for (unsigned i = 0; i < N_star; ++i)
            {
                // community : i escaped from community
                escape_matrix[i] = std::exp(-commu_escape_rate[i]);

                // person-to-person : i escaped from k
                for (unsigned k = 0; k < N_star; ++k)
                    escape_matrix[N_star * (k + 1) + i] = std::exp(-p2p_escape_rate[k][i]);

                // infector always "escapes" himself (check if actually relevant)
                escape_matrix[N_star * (i + 1) + i] = 1.;
            }

            ll += std::log(ChainBinomial_LookupTable[N_star].likelihood(escape_matrix));
        }

        return ll;
    }

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
            size_t n_patterns = EscapePatterns_LookupTable[N_star].size();
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
                std::cerr << ChainBinomial_LookupTable[N_star].str() << std::endl;

            auto start = std::chrono::high_resolution_clock::now();
            std::cerr << "likelihood: " << std::exp(chain_binomial_heter_log_lik(N_star, commu_esc_rate, p2p_esc_rate));
            auto end = std::chrono::high_resolution_clock::now();
            std::chrono::duration<double> elapsed = end - start;
            std::cerr << " | time: " << elapsed.count() << " seconds" << std::endl;
        }

        std::cerr << std::endl;
    }
    */
}
