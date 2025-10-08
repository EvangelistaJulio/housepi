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
 * binomial.h
 *
 */

#ifndef HOUSEPI__BINOMIAL__H
#define HOUSEPI__BINOMIAL__H

#include "tools/types.h"

/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////  NAMESPACE : bax  /////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

namespace bax
{
    using ushort = unsigned short;
    using ulong = unsigned long;
    using ushort_vec = std::vector<ushort>;
    using gentype = std::vector<ushort>;

    /////////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////// *** EscapePattern : Class *** /////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////

    class EscapePattern
    {
    public:
        EscapePattern(ushort N_star)
            : N_star_(N_star) {}
        ~EscapePattern() {}

        void add_entry(ushort new_entry);
        void add_pattern(ulong new_pattern_idx);
        double likelihood(const dbl_vec &esc_mtx);
        const std::string str();

    private:
        const ushort N_star_;
        size_t n_entries_ = 0;
        size_t n_patterns_ = 0;

        ushort_vec entries_ = {};
        std::vector<ulong> pattern_idxs_ = {};
    };

    /////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////  LookUp Tables  ////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////

    // Transmission Generations
    extern std::vector<ulong> Num_Gener_LookupTable; /// n. of generations
    extern matrix<ulong> Multipl_Gener_LookupTable;  /// n. of multiplicites per generation
    extern matrix<gentype> Generations_LookupTable;  /// generation std::vectors

    // Combinations
    extern matrix<ulong> Num_Combin_LookupTable;            /// LUT wiht n.e number of combinations
    extern matrix<matrix<ushort>> Combinations_LookupTable; /// LUT wiht combinations and complements

    // EscapePatterns
    extern std::vector<EscapePattern> ChainBinomial_LookupTable;
    extern matrix<EscapePattern> EscapePatterns_LookupTable;

    /////////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////  STANDARD FUNCTIONS   /////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////

    // functions for the generations lookup table
    void init_generations_lookuptable(ushort max_N);
    std::vector<gentype> recursive_span_generation(ushort n_infec);

    // functions for the combinations lookup table
    void clear_combinations_lookuptable();
    void init_combinations_lookuptable(ushort max_N);
    void store_combinations(ushort N, ushort k);
    void recursive_span_combinations(ushort N, ushort k,
                                     ushort_vec &current,
                                     std::vector<ushort_vec> &results,
                                     ushort start, ushort loop);

    // recursive likelihood formulation lookuptable
    void init_escape_patterns_lookuptable(ushort max_N,
                                          bool clean_combinations = true);
    void recursive_escape_patterns_function(ushort N_star, ushort N,
                                            ushort_vec access_matrix,
                                            EscapePattern &father,
                                            std::vector<EscapePattern> &results);

    // actual chain binomial likelihood
    double chain_binomial_heter_loglikelihood(ushort N_star,
                                              const dbl_vec &commu_esc_proba,
                                              const matrix<double> &p2p_esc_proba);
}

#endif // HOUSEPI__BINOMIAL__H
