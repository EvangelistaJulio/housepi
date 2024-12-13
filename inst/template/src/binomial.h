/*
 * binomial.h
 *
 * Copyright (C) 2024 Francesco Camaglia, Institut Pasteur
 */

#ifndef HOUSEPI__BINOMIAL__H
#define HOUSEPI__BINOMIAL__H

#include <vector>
#include <string>
#include <map>
#include <cmath>
#include <iostream>
#include <sstream>
#include <iomanip>

namespace bax
{
    /////////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////  TYPEDEF  /////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////

    using ushort = unsigned short int;      // alias for unsigned short int;
    using ulong = unsigned long int;        // alias for unsigned long int;
    using vec_ushort = std::vector<ushort>; // alias for a std::vector of ushort
    using gentype = std::vector<ushort>;    // alias for a generation std::vector

    template<typename T>
    using matrix = std::vector<std::vector<T>>;

    /////////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////  STANDARD FUNCTIONS   /////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////

    class EscapePattern
    {
    public:
        EscapePattern(ushort N_star_) : N_star(N_star_) {}
        ~EscapePattern() {}

        void add_entry(ushort new_entry);
        void add_pattern(ulong new_pattern_idx);
        double likelihood(const std::vector<double> &esc_mtx);
        const std::string show();

    private:
        const ushort N_star;
        size_t n_entries = 0;
        size_t n_patterns = 0;

        vec_ushort entries = {};
        std::vector<ulong> pattern_idxs = {};
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

    ulong factorial(ulong n);
    ulong binomial_coefficient(ushort N, ushort k);

    // functions for the generations lookup table
    void init_generations_lookuptable(ushort max_N);
    std::vector<gentype> recursive_span_generation(ushort n_infec);

    // functions for the combinations lookup table
    void clear_combinations_lookuptable();
    void init_combinations_lookuptable(ushort max_N);
    void store_combinations(ushort N, ushort k);
    void recursive_span_combinations(ushort N, ushort k,
                                     vec_ushort &current,
                                     std::vector<vec_ushort> &results,
                                     ushort start, ushort loop);

    // recursive likelihood formulation lookuptable
    void init_escape_patterns_lookuptable(ushort max_N);
    void recursive_escape_patterns_function(ushort N_star, ushort N,
                                            vec_ushort access_matrix,
                                            EscapePattern &father,
                                            std::vector<EscapePattern> &results);

    // actual chain binomial likelihood
    double chain_binomial_heter_loglikelihood(ushort N_star,
                                              const std::vector<double> &commu_esc_proba,
                                              const matrix<double> &p2p_esc_proba);
}

#endif // HOUSEPI__BINOMIAL__H
