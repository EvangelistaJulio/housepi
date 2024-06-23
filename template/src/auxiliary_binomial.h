// auxiliary_binomial.h //
// Copyright (C) 2024 Francesco Camaglia, Institut Pasteur //

#ifndef AUXILIARY_BINOMIAL__H
#define AUXILIARY_BINOMIAL__H

#include <vector>
#include <map>
#include <iostream>

// Binomial Auxiliary Namespace
namespace bax
{

    /////////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////  TYPEDEF  /////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////

    // types to optimize memory usage
    typedef unsigned short int ushort;
    // alias for a vector of ushort
    typedef std::vector<ushort> vec_ushort;
    // alias for a generation vector
    typedef std::vector<ushort> gentype;

    /////////////////////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////  GLOBAL VARIABLES  ///////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////

    // lookup-table with generations
    extern gentype **Generations_LookupTable;
    // lookup-table with number of generations
    extern unsigned *Num_Gener_LookupTable;
    // lookup-table with n. of multiplicites per generation
    extern unsigned **Multipl_Gener_LookupTable;

    // lookup-table with combinations and complements
    extern std::vector<std::vector<std::vector<std::vector<ushort>>>> Combinations_LookupTable;
    // lookup-table with the number of combinations
    extern std::vector<std::vector<unsigned>> Num_Combinations_LookupTable;

    /////////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////  STANDARD FUNCTIONS   /////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////

    unsigned factorial(unsigned n);
    unsigned binomial_coefficient(const ushort N, const ushort k);

    // functions for the generations lookup table
    std::vector<gentype> generation_span(ushort n_infec);
    void init_generations_lookuptable(const ushort max_N);

    // functions for the combinations lookup table
    void recursive_combinations_span(ushort N, ushort k, vec_ushort &current, ushort start, ushort loop, std::vector<vec_ushort> &results);
    void store_combinations(const ushort N, const ushort k);
    void init_combinations_lookuptable(const ushort max_N);

    // recursive likelihood formulation
    double recursive_function(ushort N, std::vector<double> &scape_matrix);
    double ll_chain_binomial_heterogeneous(ushort N, ushort k, ushort *indidx, std::vector<double> &escape_matrix);

}
#endif // AUXILIARY_BINOMIAL__H