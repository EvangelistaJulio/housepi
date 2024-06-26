// auxiliary_binomial.h //
// Copyright (C) 2024 Francesco Camaglia, Institut Pasteur //

#ifndef AUXILIARY_BINOMIAL__H
#define AUXILIARY_BINOMIAL__H

#include <vector>
#include <map>
#include <cmath>
#include <iostream>

// Binomial Auxiliary Namespace
namespace bax
{
    /////////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////  TYPEDEF  /////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////

    // alias for a vector of ushort
    typedef std::vector<ushort> vec_ushort;
    // alias for a generation vector
    typedef std::vector<ushort> gentype;

    /////////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////  STANDARD FUNCTIONS   /////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////

    class EscapePattern
    {
    public:
        const ushort N_star;
        vec_ushort entries = {};
        std::vector<int> pattern_idxs = {};

        EscapePattern(ushort N_star_) : N_star(N_star_) {}

        double esc_prod(std::vector<double> &esc_mtx);
        std::string show(bool latex = false);
    };

    /////////////////////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////  GLOBAL VARIABLES  ///////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////

    //===== Transmission Generations : Lookup Tables =====//
    // numbers of generations
    extern std::vector<unsigned> Num_Gener_LookupTable;
    // n. of multiplicites per generation
    extern std::vector<std::vector<unsigned>> Multipl_Gener_LookupTable;
    // generation vectors
    extern std::vector<std::vector<gentype>> Generations_LookupTable;

    //===== Combinations : Lookup Tables =====//
    // lookup-table to store the number of combinations
    extern std::vector<std::vector<unsigned>> Num_Combin_LookupTable;
    // lookup-table to store combinations and complements
    extern std::vector<std::vector<std::vector<std::vector<ushort>>>> Combinations_LookupTable;

    //===== EscapePatterns : Lookup Tables =====//
    // lookup-table to store the number of combinations
    extern std::vector<EscapePattern> Chain_Binomial_Heter_LookupTable;
    // lookup-table to store combinations and complements
    extern std::vector<std::vector<EscapePattern>> Store_Escape_Patterns_LookupTable;

    /////////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////  STANDARD FUNCTIONS   /////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////

    unsigned factorial(unsigned n);
    unsigned binomial_coefficient(const ushort N, const ushort k);

    // functions for the generations lookup table
    std::vector<gentype> recursive_span_generation(ushort n_infec);
    void init_generations_lookuptable(const ushort max_N);

    // functions for the combinations lookup table
    void recursive_span_combinations(ushort N, ushort k, vec_ushort &current, std::vector<vec_ushort> &results, ushort start, ushort loop);
    void store_combinations(const ushort N, const ushort k);
    void init_combinations_lookuptable(const ushort max_N);
    void clear_combinations_lookuptable();

    // recursive likelihood formulation
    void recursive_escape_patterns_function(ushort N_star, vec_ushort access_matrix, EscapePattern &father, std::vector<EscapePattern> &results);
    void init_escape_patterns_lookuptable(const ushort max_N);
    double chain_binomial_heter_likelihood(ushort N_star, std::vector<double> &commu_esc_proba, std::vector<std::vector<double>> &p2p_esc_proba);
}
#endif // AUXILIARY_BINOMIAL__H