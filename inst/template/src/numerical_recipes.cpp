/*
 *        __
 *       / /_  ____  __  __ ____ __  ____  __
 *      / __ \/ __ \/ / / / ___/ _ \/ __ \/ /
 *     / / / / /_/ / /_/ (__  )  __/ /_/ / /
 *    /_/ /_/\____/\____/____/\___/ ____/_/
 *                               /_/
 *
 * Copyright (C) 2025 Francesco Camaglia, Institut Pasteur
 *
 * numerical_recipes.cpp
 *
 */

#include <string>
#include <sstream>
#include <map>

#include "tools/operators.h"
#include "numerical_recipes.h"
#include "binomial.h"

/////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////// *** NUMERICAL RECIPSE NAMESPACE *** /////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////
namespace numr
{
    using std::exp;
    using std::fabs;
    using std::lgamma;
    using std::log;
    using std::pow;

    /////////////////////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////  COMBINATORICS FUNCTIONS  ////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////

    /// Binomial coefficient of N elements in k positions
    unsigned long binom_coeff(unsigned short N, unsigned short k)
    {
        if (N == 0 || k == N)
            return 1;

        unsigned long result = N;
        unsigned short stop = (k > N - k) ? k : N - k; // use symmetry
        for (unsigned short i = 1; i < stop; ++i)
            result = (result * (N - i)) / (i + 1);
        return result;
    }

    /// Logarithm of binomial coefficient.
    double log_binom_coeff(unsigned short N, unsigned short k)
    {
        if (N == 0 || k == N)
            return 0.;

        double result = log(N);
        unsigned short stop = (k > N - k) ? k : N - k; // use symmetry
        for (unsigned short i = 1; i < stop; ++i)
            result += log(N - i) - log(i + 1);
        return result;
    }

    /// Log-ratio of binomial coefficients.
    double logratio_binom_coeff(unsigned short N, unsigned short ko, unsigned short kn)
    {
        if (N == 0 || ko == kn)
            return 0;

        double result = 0.;
        if (ko > kn)
            std::swap(ko, kn);

        for (unsigned short i = ko + 1; i <= kn; ++i)
            result += log(N - i + 1) - log(i);
        return result;
    }

    /// UPDATE: To better reference this:
    /// Press, Teukolsky, Vetterling, and Flannery
    /// Numerical Recipes 3rd Edition: The Art of Scientific Computing
    namespace gamma_p_impl
    {
        /// Lower regulairzed incomplete gamma function using series expansion.
        double series(double shape, double z)
        {
            double sum = 1. / shape;
            double term = sum;
            for (unsigned m = 1; m <= MAX_ITER; ++m)
            {
                term *= z / (shape + m);
                sum += term;

                if (term < EPSILON * sum)
                    return sum * exp(-z + shape * log(z) - lgamma(shape));
            }

            throw std::overflow_error(
                "Overflow in gamma_p function for x/scale=" + std::to_string(z) +
                ", shape=" + std::to_string(shape) + ".");
        }

        /// Upper regularized incomplete gamma function using continued fraction.
        double continued_fraction(double shape, double z)
        {
            double b = z + 1. - shape;
            double c = 1. / MIN_DBL;
            double d = 1. / b;
            double h = d;
            double delta = 0.;

            for (unsigned m = 1; m <= MAX_ITER; ++m)
            {
                double an = (shape - m) * m;
                b += 2.;
                d = an * d + b;
                if (fabs(d) < MIN_DBL)
                    d = MIN_DBL;
                c = b + an / c;
                if (fabs(c) < MIN_DBL)
                    c = MIN_DBL;
                d = 1. / d;
                delta = c * d;
                h *= delta;

                if (fabs(delta - 1.) < EPSILON)
                    return 1. - exp(-z + shape * log(z) - lgamma(shape)) * h;
            }

            throw std::overflow_error(
                "Overflow in gamma_p function for x/scale=" + std::to_string(z) +
                ", shape=" + std::to_string(shape) + ".");
        }
    }

    /// Regularized incomplete gamma function.
    double gamma_p(double shape, double z)
    {
        if (z < shape + 10.)
            return gamma_p_impl::series(shape, z);
        else
            return gamma_p_impl::continued_fraction(shape, z);
    }

    /////////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////  REGULARISED INCOMPLETE BETA FUNCTION  ////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////

    /// UPDATE: to better reference this
    /// (C) Copr. 1986-92 Numerical Recipes Software Vn2-5.
    namespace ibeta_impl
    {
        double continued_fraction(double a, double b, double x)
        {
            const double qab = a + b,
                         qap = a + 1.,
                         qam = a - 1.;
            double c = 1.,
                   d = 1. - qab * x / qap;

            if (fabs(d) < MIN_DBL)
                d = MIN_DBL;
            d = 1. / d;

            double h = d;
            for (unsigned m = 1; m <= MAX_ITER; ++m)
            {
                unsigned m2 = 2 * m;
                double aa = (b - m) * m * x / ((qam + m2) * (a + m2));

                d = 1. + aa * d;
                if (fabs(d) < MIN_DBL)
                    d = MIN_DBL;

                c = 1. + aa / c;
                if (fabs(c) < MIN_DBL)
                    c = MIN_DBL;

                d = 1. / d;
                h *= d * c;
                aa = -(a + m) * (qab + m) * x / ((a + m2) * (qap + m2));

                d = 1. + aa * d;
                if (fabs(d) < MIN_DBL)
                    d = MIN_DBL;

                c = 1. + aa / c;
                if (fabs(c) < MIN_DBL)
                    c = MIN_DBL;

                d = 1. / d;

                double del = d * c;
                h *= del;

                if (fabs(del - 1.) < EPSILON)
                    return h;
            }

            throw std::overflow_error(
                "Overflow in ibeta function for x=" + std::to_string(x) +
                ", a=" + std::to_string(a) + ", b=" + std::to_string(b) + ".");
        }
    }

    /// Regularized incomplete beta function.
    double ibeta(double a, double b, double x)
    {
        if (x == 0.)
            return 0.;
        else if (x == 1.)
            return 1.;
        else
        {
            const double lbt = a * log(x) + b * log(1.0 - x) - lbeta(a, b);

            if (x < (a + 1.) / (a + b + 2.))
                return exp(lbt + log(ibeta_impl::continued_fraction(a, b, x) / a));
            else
                return 1. - exp(lbt + log(ibeta_impl::continued_fraction(b, a, 1. - x) / b));
        }
    }
} // NUMR NAMESPACE

/////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////// *** BINOMIAL AUXILIARY NAMESPACE *** /////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////
namespace bax
{
    /////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////  LookUp Tables  ////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////

    // Transmission Generations
    std::vector<ulong> Num_Gener_LookupTable; /// n. of generations
    matrix<ulong> Multipl_Gener_LookupTable;  /// n. of multiplicites per generation
    matrix<gentype> Generations_LookupTable;  /// generation std::vectors

    // Combinations
    matrix<ulong> Num_Combin_LookupTable;            /// LUT wiht n.e number of combinations
    matrix<matrix<ushort>> Combinations_LookupTable; /// LUT wiht combinations and complements

    // EscapePatterns
    std::vector<EscapePattern> ChainBinomial_LookupTable;
    matrix<EscapePattern> EscapePatterns_LookupTable;

    /////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////  Escape Pattern : Class   ///////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////

    /// Add a new entry for the escape matrix.
    void EscapePattern::add_entry(ushort new_entry)
    {
        entries_.push_back(new_entry);
        ++n_entries_;
    }

    /// Add a new pattern index.
    void EscapePattern::add_pattern(ulong new_pattern_idx)
    {
        pattern_idxs_.push_back(new_pattern_idx);
        ++n_patterns_;
    }

    /// Compute the loglikelihood starting from the current escape pattern.
    ///	NOTE: it requires initialization of the lookup table.
    double EscapePattern::likelihood(const dbl_vec &esc_mtx)
    {
        double output = 1.;

        for (size_t c = 0; c < n_patterns_; ++c)
            output -= EscapePatterns_LookupTable[N_star_][pattern_idxs_[c]].likelihood(esc_mtx);

        for (size_t c = 0; c < n_entries_; ++c)
            output *= esc_mtx[entries_[c]];

        return output;
    }

    /// String with the formula of the loglikelihood.
    const std::string EscapePattern::str()
    {
        std::ostringstream oss;

        for (size_t c = 0; c < n_entries_; ++c)
        {
            ushort u = entries_[c];
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

    /////////////////////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////  TRANSMISSION GENERATIONS  ///////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////

    /// All possible ranked offspring generations for a given number of infected individuals n > 0.
    /// @example
    ///	n_infec = 3
    ///	return [[0, 0, 0], [0, 0, 1], [0, 1, 1], [0, 1, 2]]
    std::vector<gentype> recursive_span_generation(ushort n)
    {
        std::vector<gentype> output(1, gentype(n, 0));

        for (ushort root = n - 1; root > 0; --root)
            for (gentype g : recursive_span_generation(n - root))
            {
                for (int i = 0; i < n - root; ++i)
                    g[i] += 1;
                output.push_back(gentype(root, 0) + g);
            }

        return output;
    }

    /// Initialize the Look Up Table of count-compressed generations.
    void init_generations_lookuptable(ushort max_N)
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
        std::vector<ulong> aux_num_lut;
        std::vector<gentype> aux_gen_lut;

        for (int n = 1; n <= max_N; ++n)
        {
            aux_num_lut.clear();
            aux_gen_lut.clear();

            for (gentype g : recursive_span_generation(n))
            {
                std::map<ushort, ushort> freq_map = {};
                gentype counts = {};
                ulong multiplicities = numr::factorial(n);

                // I take advantage of the order granted by the map
                for (ushort val : g)
                    if (freq_map.count(val))
                        ++freq_map[val];
                    else
                        freq_map.emplace(val, 1);

                // assign counts
                for (const auto &x : freq_map)
                {
                    counts.push_back(x.second);
                    multiplicities /= numr::factorial(x.second);
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

            // copy from auxiliary std::vectors
            for (size_t c = 0; c < tot_num_gen; ++c)
            {
                Multipl_Gener_LookupTable[n][c] = aux_num_lut[c];
                Generations_LookupTable[n][c] = aux_gen_lut[c];
            }
        }
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
    void recursive_span_combinations(ushort N, ushort k, ushort_vec &current, std::vector<ushort_vec> &results, ushort start, ushort loop)
    {
        if (loop > k - 1)
        {
            // complement of the combination
            // mark elements present in the list
            std::vector<bool> missing(N, true);
            for (ushort i = 0; i < k; ++i)
                missing[current[i]] = false;

            // store complement
            for (ushort i = 0; i < N; ++i)
                if (missing[i])
                    current[loop++] = i;

            // append current combination & complement to results
            results.push_back(current);
        }
        else
            for (ushort i = start; i < N; ++i)
            {
                current[loop] = i;
                recursive_span_combinations(N, k, current, results, i + 1, loop + 1);
            }
    }

    /// Store all combinations + complement of N elements in k positions.
    void store_combinations(ushort N, ushort k)
    {
        // Generate combinations with the recurisve algorithm
        std::vector<ushort_vec> results;
        ushort_vec current(N);
        recursive_span_combinations(N, k, current, results, 0, 0);

        // Store all results inside the Lookup Table
        ulong c = 0;
        for (ushort_vec current : results)
        {
            Combinations_LookupTable[N][k][c].resize(N);

            // store current combination + complement
            for (ushort i = 0; i < N; ++i)
                Combinations_LookupTable[N][k][c][i] = current[i];
            ++c;
        }
    }

    /// Initialize the Look Up Table of combinations and their counts.
    void init_combinations_lookuptable(ushort max_N)
    {
        ulong num_combin;

        // Assure that LUT are empty before initializing
        Num_Combin_LookupTable.clear();
        Combinations_LookupTable.clear();

        Num_Combin_LookupTable.resize(max_N + 1);
        Combinations_LookupTable.resize(max_N + 1);

        for (ushort N = 0; N <= max_N; ++N)
        {
            Num_Combin_LookupTable[N].resize(max_N);
            Combinations_LookupTable[N].resize(max_N);

            for (ushort k = 0; k < N; ++k)
            {
                num_combin = numr::binom_coeff(N, k);
                Num_Combin_LookupTable[N][k] = num_combin;
                Combinations_LookupTable[N][k].resize(num_combin);
                store_combinations(N, k);
            }
        }
    }

    /// Clear the Look Up Table of combinations and their counts
    void clear_combinations_lookuptable()
    {
        Num_Combin_LookupTable.clear();
        Num_Combin_LookupTable.shrink_to_fit();
        Combinations_LookupTable.clear();
        Combinations_LookupTable.shrink_to_fit();
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////  Recursive Escape Patterns Functions  /////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////

    /// Recursive function to generate all possible escape patterns for a given number of seropositives N_star.
    void recursive_escape_patterns_function(ushort N_star, ushort N, ushort_vec access_matrix, EscapePattern &father, std::vector<EscapePattern> &results)
    {
        for (ushort k = 0; k < N; ++k)
            for (ulong c = 0; c < Num_Combin_LookupTable[N][k]; ++c)
            {
                EscapePattern new_child(N_star);

                // associate entries for this pattern
                std::vector<ushort> &indidx = Combinations_LookupTable[N][k][c];
                for (ushort i = k; i < N; ++i)
                {
                    //  community
                    new_child.add_entry(access_matrix[indidx[i]]);

                    // person-to-person
                    for (ushort j = 0; j < k; ++j)
                        new_child.add_entry(access_matrix[indidx[i] + N * (indidx[j] + 1)]);
                }

                if (k > 0)
                {
                    // find the submatrix of the given infectious
                    // note: always at least 1 infector (community)
                    ushort_vec sub_access_matrix((k + 1) * k);

                    // community (inf=0) always at row=0
                    for (ushort row = 0, inf = 0; row <= k; inf = 1 + indidx[row++])
                        for (ushort col = 0; col < k; ++col)
                            sub_access_matrix[col + k * row] = access_matrix[indidx[col] + N * inf];

                    recursive_escape_patterns_function(N_star, k, sub_access_matrix, new_child, results);
                    sub_access_matrix.clear();
                }

                father.add_pattern(results.size()); // attach pattern index
                results.push_back(new_child);       // store pattern
            }
    }

    /// Initialize the Look Up Table of escape patterns and the access object.
    void init_escape_patterns_lookuptable(ushort max_N, bool clean_combinations)
    {
        // lookuptable with all combinations + complements
        init_combinations_lookuptable(max_N);

        // assure that LUT are empty before initializing
        ChainBinomial_LookupTable.clear();
        EscapePatterns_LookupTable.clear();

        // lookuptable with escape patterns linked to chain-binomial
        // case N_star == 0, for consistency
        ChainBinomial_LookupTable.emplace_back(0);
        EscapePatterns_LookupTable.push_back(std::vector<EscapePattern>());

        // case N_star > 0
        for (ushort N_star = 1; N_star <= max_N; ++N_star)
        {
            ushort_vec access_matrix((N_star + 1) * N_star);
            for (ushort i = 0; i < (N_star + 1) * N_star; ++i)
                access_matrix[i] = i;

            // store to lookup table
            ChainBinomial_LookupTable.emplace_back(N_star);
            EscapePatterns_LookupTable.push_back(std::vector<EscapePattern>());

            recursive_escape_patterns_function(N_star, N_star, access_matrix,
                                               ChainBinomial_LookupTable[N_star],
                                               EscapePatterns_LookupTable[N_star]);
        }

        // deallocate memory of combinations lookup table
        if (clean_combinations)
            clear_combinations_lookuptable();
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////  Heterogeneous Chain-Binomisl LogLikelihood  /////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////

    /// Heterogeneous chain binomial loglikelihood for ordered probabilities (seropositives, seronegatives).
    /// @param N_star n. of seropositves
    /// @param commu_esc_proba ordered community escape probabilities
    /// @param p2p_esc_proba ordered person-to-person escape probabilities
    double chain_binomial_heter_loglikelihood(ushort N_star,
                                              const dbl_vec &commu_esc_proba,
                                              const matrix<double> &p2p_esc_proba)

    {
        ushort hh_size = commu_esc_proba.size();
        double ll = 0.;

        // (-) seronegatives contribution
        for (ushort i = N_star; i < hh_size; ++i)
        {
            // community : i escaped from community
            ll += std::log(commu_esc_proba[i]);

            // person-to-person : i escaped from k
            for (ushort k = 0; k < N_star; ++k)
                ll += std::log(p2p_esc_proba[k][i]);
        }

        // (+) seropositives contribution
        if (N_star == 1) // single case
            ll += std::log(1. - commu_esc_proba[0]);

        else if (N_star > 1) // multiple cases
        {
            // flatten sub-matrix with only seropositives
            dbl_vec escape_matrix((N_star + 1) * N_star, 1.);

            for (ushort i = 0; i < N_star; ++i)
            {
                // community : i escaped from community
                escape_matrix[i] = commu_esc_proba[i];

                // person-to-person : i escaped from k
                for (ushort k = 0; k < N_star; ++k)
                    escape_matrix[N_star * (k + 1) + i] = p2p_esc_proba[k][i];
            }

            ll += std::log(ChainBinomial_LookupTable[N_star].likelihood(escape_matrix));
        }

        return ll;
    }
}

/*
/// Debug mode, to check values of the likelihood and print formulation.
///
/// g++ -Wall -Wextra -std=gnu++11 -g -O2 -o auxiliary_binomial.o auxiliary_binomial.cpp ./auxiliary_binomial.o 7
///
/// Auxiliary function to perform debugs
void auxiliary_chain_binomial_debugger(ushort max_num, ushort limit)
{
    using namespace bax;
    for (ushort N_star = 0; N_star <= max_num; ++N_star)
    {
        size_t n_patterns = EscapePatterns_LookupTable[N_star].size();
        std::cerr << "\n<<" << N_star << ">>"
                    << " n.terms = " << ((n_patterns == 0) ? 0 : 1 + (n_patterns - 1) / 2) << std::endl;

        // elements for the escape matrix //
        dbl_vec commu_esc_proba(N_star);
        for (ushort i = 0; i < N_star; ++i)
        {
            commu_esc_proba[i] = 0.15 + 0.025 * i;
            if (N_star < limit)
                std::cerr << "B" << i + 1 << " = " << commu_esc_proba[i] << std::endl;
        }
        matrix<double> p2p_esc_proba(N_star);
        for (ushort row = 0; row < N_star; ++row)
        {
            p2p_esc_proba[row].resize(N_star);
            for (ushort col = 0; col < N_star; ++col)
            {
                p2p_esc_proba[row][col] = 0.25 + 0.0125 * col + 0.0125 * row;
                if (N_star < limit)
                    std::cerr << "Q" << row + 1 << "_" << col + 1 << " = " << p2p_esc_proba[row][col] << std::endl;
            }
        }

        if (N_star < limit)
            std::cerr << ChainBinomial_LookupTable[N_star].str() << std::endl;

        auto start = std::chrono::high_resolution_clock::now();
        std::cerr << "likelihood: " << std::exp(chain_binomial_heter_loglikelihood(N_star, commu_esc_proba, p2p_esc_proba));
        auto end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> elapsed = end - start;
        std::cerr << " | time: " << elapsed.count() << " seconds" << std::endl;
    }

    std::cerr << std::endl;
}

int main(int argc, char *argv[])
{
    using namespace bax;

    int max_num = 4;
    if (argc > 1)
        max_num = std::stoi(argv[1]);

    init_escape_patterns_lookuptable(max_num);
    auxiliary_chain_binomial_debugger(max_num, 4);
}
*/
