/*
 * binomial.h
 *
 * Copyright (C) 2024 Francesco Camaglia, Institut Pasteur
 */

#ifndef BINOMIAL__H
#define BINOMIAL__H

#include <vector>
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

    using std::vector;
    using std::string;
    typedef unsigned short int ushort;      // alias for unsigned short int;
    typedef unsigned long int ulong;        // alias for unsigned long int;
    typedef vector<ushort> vec_ushort; // alias for a vector of ushort
    typedef vector<ushort> gentype;    // alias for a generation vector

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
        double likelihood(vector<double> &esc_mtx);
        string show();

    private:
        const ushort N_star;
        size_t num_entries = 0;
        size_t num_patterns = 0;

        vec_ushort entries = {};
        vector<ulong> pattern_idxs = {};
    };

    /////////////////////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////  GLOBAL VARIABLES  ///////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////

    //===== Transmission Generations : Lookup Tables =====//
    // numbers of generations
    extern vector<ulong> Num_Gener_LookupTable;
    // n. of multiplicites per generation
    extern vector<vector<ulong>> Multipl_Gener_LookupTable;
    // generation vectors
    extern vector<vector<gentype>> Generations_LookupTable;

    //===== Combinations : Lookup Tables =====//
    // lookup-table to store the number of combinations
    extern vector<vector<ulong>> Num_Combin_LookupTable;
    // lookup-table to store combinations and complements
    extern vector<vector<vector<vector<ushort>>>> Combinations_LookupTable;

    //===== EscapePatterns : Lookup Tables =====//
    // lookup-table to store the number of combinations
    extern vector<EscapePattern> ChainBinomial_LookupTable;
    // lookup-table to store combinations and complements
    extern vector<vector<EscapePattern>> EscapePatterns_LookupTable;

    /////////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////  STANDARD FUNCTIONS   /////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////

    ulong factorial(ulong n);
    ulong binomial_coefficient(ushort N, ushort k);

    // functions for the generations lookup table
    vector<gentype> recursive_span_generation(ushort n_infec);
    void init_generations_lookuptable(ushort max_N);

    // functions for the combinations lookup table
    void recursive_span_combinations(ushort N, ushort k, vec_ushort &current, vector<vec_ushort> &results, ushort start, ushort loop);
    void store_combinations(ushort N, ushort k);
    void init_combinations_lookuptable(ushort max_N);
    void clear_combinations_lookuptable();

    // recursive likelihood formulation
    void recursive_escape_patterns_function(ushort N_star, ushort N, vec_ushort access_matrix, EscapePattern &father, vector<EscapePattern> &results);
    void init_escape_patterns_lookuptable(ushort max_N);

    double chain_binomial_heter_loglikelihood(ushort N_star, vector<double> &commu_esc_proba, vector<vector<double>> &p2p_esc_proba);
}

#endif // BINOMIAL__H