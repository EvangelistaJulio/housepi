/*      __
 *     / /_  ____  __  __ ____ __  ____  __
 *    / __ \/ __ \/ / / / ___/ _ \/ __ \/ /
 *   / / / / /_/ / /_/ (__  )  __/ /_/ / /
 *  /_/ /_/\____/\____/____/\___/ ____/_/
 *                             /_/
 *
 * Copyright (C) 2025 Francesco Camaglia, Institut Pasteur
 *
 * binomial.h
 *
 */

#ifndef HOUSEPI__BINOMIAL__H
#define HOUSEPI__BINOMIAL__H

#include "tools.h"

namespace bax
{
    /////////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////  TYPEDEF  /////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////

    using ushort_vec = std::vector<ushort>;
    using gentype = std::vector<ushort>;

    /////////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////  STANDARD FUNCTIONS   /////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////

    class EscapePattern
    {
    public:
        EscapePattern(ushort N_star)
            : _N_star(N_star) {}
        ~EscapePattern() {}

        void add_entry(ushort new_entry);
        void add_pattern(ulong new_pattern_idx);
        double likelihood(const dbl_vec &esc_mtx);
        const std::string show();

    private:
        const ushort _N_star;
        size_t _n_entries = 0;
        size_t _n_patterns = 0;

        ushort_vec _entries = {};
        std::vector<ulong> _pattern_idxs = {};
    };

    /////////////////////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////  GLOBAL VARIABLES  ///////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////

    //===== Transmission Generations : Lookup Tables =====//
    // numbers of generations
    extern std::vector<ulong> Num_Gener_LookupTable;
    // n. of multiplicites per generation
    extern matrix<ulong> Multipl_Gener_LookupTable;
    // generation std::vectors
    extern matrix<gentype> Generations_LookupTable;

    //===== Combinations : Lookup Tables =====//
    // lookup-table to store the number of combinations
    extern matrix<ulong> Num_Combin_LookupTable;
    // lookup-table to store combinations and complements
    extern std::vector<std::vector<matrix<ushort>>> Combinations_LookupTable;

    //===== EscapePatterns : Lookup Tables =====//
    // lookup-table to store the number of combinations
    extern std::vector<EscapePattern> ChainBinomial_LookupTable;
    // lookup-table to store combinations and complements
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
    void init_escape_patterns_lookuptable(ushort max_N);
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
