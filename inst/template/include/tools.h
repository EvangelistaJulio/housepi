/*      __
 *     / /_  ____  __  __ ____ __  ____  __
 *    / __ \/ __ \/ / / / ___/ _ \/ __ \/ /
 *   / / / / /_/ / /_/ (__  )  __/ /_/ / /
 *  /_/ /_/\____/\____/____/\___/ ____/_/
 *                             /_/
 *
 * Copyright (C) 2025 Francesco Camaglia, Institut Pasteur
 *
 * tools.h
 *
 */

#ifndef HOUSEPI__TOOLS__H
#define HOUSEPI__TOOLS__H

#include <fstream>
#include <iomanip>
#include <sstream>
#include <iostream>

#include <cmath>
#include <algorithm>
#include <limits>
#include <random>
#include <string>
#include <vector>
#include <map>
#include <unordered_map>
#include <functional>
#include <memory>

/// DEVELOPER: flag used to develop Rcpp chunks
// #define _HOUSEPI_RCPP

#ifdef _HOUSEPI_RCPP
#include <Rcpp.h>
#endif

#ifdef _OPENMP
#include <omp.h>
#endif

/////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////  GLOBAL VARIABLES  ///////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

// attribute aliases
#define HSP_ATTR_ARGUMENT true // default argument
#define HSP_ATTR_RANGE false   // range argument

// variable parts aliases
#define HSP_VAR_FUNCTION false // pass to Variable function
#define HSP_VAR_MODEL true     //  pass to Variable model

// NA value for double
#ifdef _HOUSEPI_RCPP
const double HSP_NA_DOUBLE = NA_REAL;
#else
const double HSP_NA_DOUBLE = std::numeric_limits<double>::quiet_NaN();
#endif

// NA value for string
const std::string HSP_NA_STRING = "<NA>";
// prefic for hidden objects
const std::string HSP_HIDE_PRFX = "..";

// chain column name in draws
const std::string HSP_CHAIN_COL = HSP_HIDE_PRFX + "chain" + HSP_HIDE_PRFX;

// group rates
// alias for reference value
constexpr double REFERENCE_VALUE = 1.0;
// reference choice
const std::string REFERENCE = HSP_HIDE_PRFX + "reference" + HSP_HIDE_PRFX;
// contact rate
const std::string CONTACT_RATE = HSP_HIDE_PRFX + "contact" + HSP_HIDE_PRFX;
// missing rate
const std::string MISSING_RATE = HSP_HIDE_PRFX + "missing_rate" + HSP_HIDE_PRFX;

// infectivity profile
// augmented instantaneous infectivity profile
const std::string HSP_PROFILE_INST = HSP_HIDE_PRFX + "profile_inst" + HSP_HIDE_PRFX;
// augmented cumulative infectivity profile
const std::string HSP_PROFILE_CUMUL = HSP_HIDE_PRFX + "profile_cumul" + HSP_HIDE_PRFX;

// advanced setup constants
constexpr unsigned MAX_STR_PRINT = 5;        // maximum options to print
constexpr unsigned MIN_STR_DIST = 3;         // minimum string distance for find closest
constexpr unsigned DEFAULT_MAX_SIZE = 30;    // default size for lookup table
constexpr unsigned MAX_NUM_ATTEMPTS = 100;   // max n. of attempts for initialization
constexpr double TARGET_ACCEPT_RATE = 0.234; // target acceptation rate for MCMC

// file formatters
constexpr char C_DELIM = ',';      // file delimiter
constexpr char C_LIST_DELIM = '/'; // dictionary delimiter

//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//

/////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////  32bit/64bit  ///////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

/* Use pseudo-random number generator for 64bit architecture. */
#if __x86_64__ || __ppc64__ || __aarch64__ || _WIN64
using rnd_gen = std::mt19937_64;
#else
using rnd_gen = std::mt19937;
#endif

/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////  TYPEDEF  /////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

using int_vec = std::vector<int>;
using dbl_vec = std::vector<double>;
using str_vec = std::vector<std::string>;

// matrix
template <typename T>
using matrix = std::vector<std::vector<T>>;

// ordered dictionaries
using dict_str = std::map<std::string, std::string>;
using dict_str_vec = std::map<std::string, str_vec>;
using dict_dbl_vec = std::map<std::string, dbl_vec>;

// unordered dictionaries
template <typename T>
using u_dict = std::unordered_map<std::string, T>;

//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//

/////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////  OPERATOR OVERLOADS  ///////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

/// Check if container containes given value (% overload)
template <typename Container, typename T>
bool operator%(const Container &c, const T &value)
{
    return std::find(c.begin(), c.end(), value) != c.end();
}

/// Concatenate containers (+ overload)
template <typename T>
std::vector<T> operator+(const std::vector<T> &x, const std::vector<T> &y)
{
    std::vector<T> vec;
    vec.reserve(x.size() + y.size());
    vec.insert(vec.end(), x.begin(), x.end());
    vec.insert(vec.end(), y.begin(), y.end());
    return vec;
}

/// Concatenate containers (+ overload)
template <typename T>
std::vector<T> &operator+=(std::vector<T> &x, const std::vector<T> &y)
{
    x.reserve(x.size() + y.size());
    x.insert(x.end(), y.begin(), y.end());
    return x;
}

/// @brief Overload minus operator.
/// Remove all elements of first containers contained in second.
template <typename T>
std::vector<T> operator-(std::vector<T> v, const std::vector<T> &w)
{
    auto new_end = std::remove_if(v.begin(), v.end(),
                                  [&w](const T &i)
                                  { return w % i; });
    v.erase(new_end, v.end());
    return v;
}

/// @brief Safely access an element in a u_dict by key.
/// If the key does not exist, an exception is thrown.
/// @tparam T
/// @param m
/// @param k
/// @return
template <typename T>
inline T &safe_at(u_dict<T> &m, const char *what, const std::string &k)
{
    auto it = m.find(k);
    if (it != m.end())
        return it->second;
    throw std::out_of_range(std::string("Unknown ") + what + " \"" + k + "\".");
}
template <typename T>
inline const T &safe_at(const u_dict<T> &m, const char *what, const std::string &k)
{
    auto it = m.find(k);
    if (it != m.end())
        return it->second;
    throw std::out_of_range(std::string("Unknown ") + what + " \"" + k + "\".");
}

//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//

/////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////  Basic Tools Functions  //////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

namespace tools
{
    // print warnings to stderr
    void warning(const std::string &txt);
    // print error message to stderr
    void error(const std::string &txt);

    // internal undefined string values
    bool is_na(const std::string &txt);
    // internal undefined double values
    bool is_na(double x);

    // default prng seed for the setup
    std::int32_t default_setup_seed();

    // Welford online algorithm for mean and variance
    std::pair<double, double> compute_stats(const dbl_vec &data);

    namespace str
    {
        // convert to upper\lower case
        std::string tolow(const std::string s);
        std::string toup(const std::string s);
        str_vec tolow(const str_vec &svec);
        str_vec toup(const str_vec &svec);

        // make list composed of distinct entries
        str_vec drop_duplicate(const str_vec &v);
        // drop private strings from list
        str_vec drop_private(str_vec v);

        // double formatter for the case NA
        std::string format(double x);
        // boolean formatter
        std::string format(bool condition);
        // format a list of string options
        std::string format(const str_vec &options, size_t max_print);

        // calculate Levenshtein distance
        unsigned short levenshtein(const std::string &s1, const std::string &s2);
        // find the closest match
        std::string find_closest(const std::string &input,
                                 const str_vec &options,
                                 unsigned short min_dist,
                                 // UPDATE: I find this annoying at the moment...
                                 bool show_avail_options = false);

        // converter functions
        bool tobool(const std::string &tkn);
        int toi(const std::string &input, const std::string &errtxt = "");
        unsigned long toul(const std::string &input, const std::string &errtxt = "");
        std::int32_t toi32(const std::string &input, const std::string &errtxt = "");
        double tod(const std::string &input, const std::string &errtxt = "");
    }

    // argument string parser
    str_vec parse_list_args(const std::string &tkn, const char delim);

    // read csv list
    str_vec headers_from_deflist(const dict_str &deflist);
    dict_str read_csv_list(const std::string &pathtofile);
    // read csv dataframe
    str_vec headers_from_deftable(const dict_str_vec &deftable);
    dict_str_vec read_csv_table(const std::string &pathtofile);

    // check for unused entries from a vector of options and throw exception
    void check_for_missing_entries(const str_vec &entries,
                                   const str_vec &options,
                                   const std::string &where = "column");
    void check_for_unused_entries(const str_vec &entries,
                                  const str_vec &options,
                                  const std::string &where = "column");

#ifdef _HOUSEPI_RCPP
    Rcpp::String stoRString(const std::string &x);
    dict_str rlist_to_map(SEXP input);
    dict_str_vec rdataframe_to_map(SEXP input);
#endif
}

#endif // HOUSEPI__TOOLS__H
