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
/////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////  32bit/64bit  ///////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

// Use pseudo-random number generator for 64bit architecture.
#if __x86_64__ || __ppc64__ || __aarch64__ || _WIN64
using rnd_gen = std::mt19937_64;
#else
using rnd_gen = std::mt19937;
#endif

/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////  TYPEDEF  /////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

using observ_dbl_ptr = const double *;

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

/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////  GLOBAL VARIABLES  ///////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

const double HSP_NA_DOUBLE = std::numeric_limits<double>::quiet_NaN(); /// NA value for double.
#define HSP_NA_STRING "<NA>"                                           /// NA value for string.

// private object names
#define HSP_HIDE_PRFX ".."               /// prefic for hidden objects
#define HSP_CHAIN_COL "..chain.."        /// chain column name in draws
#define HSP_IDENTITY_FUNC "..identity.." /// identity function name
// group rates
constexpr double HSP_REFER_VALUE = 1.0;     /// alias for reference value
#define HSP_REFERENCE "..reference.."       /// reference choice
#define HSP_CONTACT_RATE "..contact.."      /// contact rate
#define HSP_MISSING_RATE "..missing_rate.." /// missing rate
// infectivity profile
#define HSP_PROFILE_INST "..profile_inst.."   /// augmented instantaneous infectivity profile
#define HSP_PROFILE_CUMUL "..profile_cumul.." /// augmented cumulative infectivity profile

// advanced setup constants
constexpr unsigned HSP_MAX_STR_PRINT = 5;       /// maximum options to print
constexpr unsigned HSP_MIN_STR_DIST = 3;        /// minimum string distance for find closest
constexpr unsigned HSP_DEFAULT_MAX_SIZE = 30;   /// default size for lookup table
constexpr unsigned HSP_MAX_NUM_ATTEMPTS = 100;  /// max n. of attempts for initialization
constexpr double HSP_TARGET_ACCEPT_RATE = 0.25; /// target acceptation rate for MCMC
constexpr double HSP_TIME_STEP = 0.05;          /// time step for real time simulations

// file formatters
#define HSP_DELIM ','      /// file delimiter
#define HSP_LIST_DELIM ':' /// dictionary delimiter

// special parameter names
#define HSP_PARAM_LOC "loc" /// location parameter
#define HSP_PARAM_MIN "min" /// minimum value
#define HSP_PARAM_MAX "max" /// maximum value
#define HSP_RANGE_KEYS str_vec({HSP_PARAM_MIN, HSP_PARAM_MAX})

// declaration for typed file inputs
#define HSP_IS_TYPE_PARAM(s) ((s) == "P")  /// == Parameter type
#define HSP_IS_TYPE_DATA(s) ((s) == "D")   /// == Data type
#define HSP_IS_TYPE_SCALAR(s) ((s) == "S") /// == Scalar type
#define HSP_IS_TYPE_AUGM(s) ((s) == "A")   /// == Augmented type

// discretization types
#define HSP_DISCR_NONE 0    /// none
#define HSP_DISCR_LEFT 1    /// left
#define HSP_DISCR_CENTER 2  /// center
#define HSP_DISCR_DEFINED 3 /// defined

// status aliases
#define HSP_INFECTED_STATUS 1 // infected status (>=)
#define HSP_HEALTHY_STATUS 0  // healthy status (==)
#define HSP_MISSING_STATUS -1 // missing status (<=)

#define HSP_IS_INFECTED(s) ((s) >= HSP_INFECTED_STATUS)
#define HSP_IS_HEALTHY(s) ((s) == HSP_HEALTHY_STATUS)
#define HSP_IS_MISSING(s) ((s) <= HSP_MISSING_STATUS)

/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////  OPERATOR OVERLOADS  ///////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

/// Check if container containes given value (% overload)
template <typename Container, typename T>
bool operator%(const Container &c, const T &x)
{
    return std::find(c.begin(), c.end(), x) != c.end();
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
/// If the key does not exist, an exception is thrown with the key name.
/// @tparam T
/// @param m
/// @param k
/// @return
template <typename T>
inline T &safe_at(u_dict<T> &m, const std::string &k, const char *what = "key")
{
    auto it = m.find(k);
    if (it != m.end())
        return it->second;
    throw std::out_of_range(std::string("Unknown ") + what + " \"" + k + "\".");
}
template <typename T>
inline const T &safe_at(const u_dict<T> &m, const std::string &k, const char *what = "key")
{
    auto it = m.find(k);
    if (it != m.end())
        return it->second;
    throw std::out_of_range(std::string("Unknown ") + what + " \"" + k + "\".");
}

/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////  NAMESPACE : tools  ////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

namespace tools
{
    /// Display a single-line warning message.
    void warning(const std::string &txt);
    /// Display a single-line error message.
    void error(const std::string &txt);

    // internal undefined double values
    bool is_na(double x);
    // internal undefined string values
    bool is_na(const std::string &txt);
    // vector of undefined values of type t
    template <typename T>
    bool is_any_na(const std::vector<T> &vec)
    {
        for (const T &x : vec)
            if (is_na(x))
                return true;
        return false;
    }

    namespace str
    {
        // convert to upper\lower case
        std::string tolow(const std::string s);
        std::string toup(const std::string s);
        str_vec tolow(const str_vec &svec);
        str_vec toup(const str_vec &svec);

        // converter functions
        bool tobool(const std::string &tkn);
        int toi(const std::string &input, const std::string &errtxt = "");
        unsigned long toul(const std::string &input, const std::string &errtxt = "");
        std::int32_t toi32(const std::string &input, const std::string &errtxt = "");
        double tod(const std::string &input, const std::string &errtxt = "");

        /// Format a dobule value (it deals with NA).
        std::string format(double x);
        /// Format a boolean value.
        std::string format(bool condition);
        /// Format a list of string options.
        std::string format(const str_vec &options, size_t max_print = HSP_MAX_STR_PRINT);
        /// Format a 2D matrix.
        std::string format(const matrix<double> &vec);

        // make list composed of distinct entries
        str_vec unique(const str_vec &v);
        // drop private strings from list
        str_vec drop_private(str_vec v);

        // find the closest match
        std::string find_closest(const std::string &input,
                                 const str_vec &options,
                                 unsigned short min_dist,
                                 // UPDATE: I find this annoying at the moment...
                                 bool show_avail_options = false);

        // check for unused entries from a vector of options and throw exception
        void assert_missing_entries(const str_vec &entries,
                             const str_vec &options,
                             const std::string &where = "column");
        void assert_unused_entries(const str_vec &entries,
                            const str_vec &options,
                            const std::string &where = "column");
    }

    namespace read
    {
        // read csv list
        str_vec list_headers(const dict_str &deflist);
        dict_str csv_list(const std::string &pathtofile);
        // read csv dataframe
        str_vec table_headers(const dict_str_vec &deftable);
        dict_str_vec csv_table(const std::string &pathtofile);
    }
    
    #ifdef _HOUSEPI_RCPP
            Rcpp::String stoRString(const std::string &x);
            dict_str rlist_to_map(SEXP input);
            dict_str_vec rdataframe_to_map(SEXP input);
    #endif
}

#endif // HOUSEPI__TOOLS__H
