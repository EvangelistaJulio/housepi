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

#include <limits>
#include <memory>
#include <random>
#include <string>
#include <vector>
#include <map>
#include <unordered_map>
#include <functional>
#include <algorithm>

/* DEVELOPER: flag used to predispose to different wrapper
#ifndef _HOUSEPI_RCPP
#define _HOUSEPI_RCPP
#endif
*/

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

// status aliases
#define INFECTED_STATUS 1 // infected status (>=)
#define HEALTHY_STATUS 0  // healthy status (==)
#define MISSING_STATUS -1 // missing status (<=)

// group rates
constexpr double REFERENCE_VALUE = 1.0;                      // alias for reference value
const std::string REFERENCE = "REFERENCE";       // reference choice
const std::string CONTACT_RATE = "Contact";      // contact rate
const std::string MISSING_RATE = "MISSING_RATE"; // missing rate

// auxiliary infectivity profile
const std::string AUX_IP_INST_PROFILE = "__profile_inst__";   // augmented instantaneous infectivity profile
const std::string AUX_IP_CUMUL_PROFILE = "__profile_cumul__"; // augmented cumulative infectivity profile

// constants
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

using ushort = unsigned short int;
using ulong = unsigned long int;

using int_vec = std::vector<int>;
using dbl_vec = std::vector<double>;
using str_vec = std::vector<std::string>;
using dict_str_str = std::map<std::string, std::string>;
using dict_str_vec = std::map<std::string, str_vec>;
using dict_dbl_vec = std::map<std::string, dbl_vec>;

template <typename T>
using matrix = std::vector<std::vector<T>>;

//----- Aliases for improved readibility -----//
using indiv = unsigned short int;

using _access_key_ = std::string;

using par_key = _access_key_;
using augm_key = _access_key_;
using data_key = _access_key_;
using rate_key = _access_key_;
using cache_key = _access_key_;

#include <unordered_map>
template <typename T>
using par_key_map = std::unordered_map<par_key, T>;
template <typename T>
using augm_key_map = std::unordered_map<augm_key, T>;
template <typename T>
using data_key_map = std::unordered_map<data_key, T>;
template <typename T>
using rate_key_map = std::unordered_map<rate_key, T>;
template <typename T>
using cache_key_map = std::unordered_map<cache_key, T>;

template <typename T>
using rate_key_mapvec = rate_key_map<std::vector<T>>;

//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//

/////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////  OPERATOR OVERLOADS  ///////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

/* Concatenate with overloaded + operator. */
template <typename T>
std::vector<T> operator+(const std::vector<T> &x, const std::vector<T> &y)
{
    std::vector<T> vec;
    vec.reserve(x.size() + y.size());
    vec.insert(vec.end(), x.begin(), x.end());
    vec.insert(vec.end(), y.begin(), y.end());
    return vec;
}

/* Append with overloaded += operator. */
template <typename T>
std::vector<T> &operator+=(std::vector<T> &x, const std::vector<T> &y)
{
    x.reserve(x.size() + y.size());
    x.insert(x.end(), y.begin(), y.end());
    return x;
}

// check if a given value is inside a container
template <typename Container, typename T>
bool contains(const Container &c, const T &value)
{
    return std::find(c.begin(), c.end(), value) != c.end();
}

// return all elements of first vector, but those contained in second
template <typename T>
std::vector<T> subtr(std::vector<T> v, const std::vector<T> &w)
{
    auto new_end = std::remove_if(v.begin(), v.end(),
                                  [&w](const T &i)
                                  { return contains(w, i); });
    v.erase(new_end, v.end());
    return v;
}

/////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////  Basic Tools Functions  //////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

namespace tools
{
    // internal undefined string values
    bool is_na(const std::string &txt);
    // internal undefined double values
    bool is_na(double x);

    // output formatter
    const std::string percentage(double fraction, ushort precision = 1);
    const std::string truth(bool condition);

    // stodut print
    void print_stdout(const std::string &txt);
    void print_stderr(const std::string &txt);

    // print warnings to stderr
    void warning(const std::string &txt);

    // print error message to stderr
    void error(const std::string &txt);

    // convert to upper\lower case
    const std::string stolower(const std::string s);
    const std::string stoupper(const std::string s);

    // converter functions
    bool stobool(const std::string &tkn);
    int stoi(const std::string &input, const std::string &errtxt = "");
    unsigned long stoul(const std::string &input, const std::string &errtxt = "");
    std::int32_t stoi32(const std::string &input, const std::string &errtxt = "");
    double stod(const std::string &input, const std::string &errtxt = "");

    // formatter for the case NA
    const std::string na_to_str(double x);

    // default prng seed for the setup
    std::int32_t default_setup_seed();

    // format a list of string options
    const std::string show_list(const str_vec &options,
                                size_t max_print = 5);
    // calculate Levenshtein distance
    ushort levenshtein(const std::string &s1, const std::string &s2);
    // find the closest match
    const std::string find_closest(const std::string &input,
                                   const str_vec &options,
                                   ushort min_dist = 2);

    // argument string parser
    const str_vec parse_list_args(const std::string &tkn,
                                  const char delim = C_LIST_DELIM);

    // read csv files
    void read_csv_list(const std::string &pathtofile,
                       dict_str_str &csv_parsed_list);
    void read_csv_table(const std::string &pathtofile,
                        dict_str_vec &csv_parsed_table);

    const str_vec headers_from_deflist(const dict_str_str &deflist);
    const str_vec headers_from_deftable(const dict_str_vec &deftable);

    // Welford online algorithm for mean and variance
    const std::pair<double, double> compute_stats(const dbl_vec &data);

    // check for unused entries from a vector of options and throw exception
    void check_for_missing_entries(const str_vec &entries, const str_vec &options);
    void check_for_unused_entries(const str_vec &entries, const str_vec &options);
}

//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//

/////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////  Progress Bar Namespace  /////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

/* Class for custom progress bar. */
namespace progress
{
    class Bar
    {
    public:
        //===== CONSTRUCTOR =====//
        Bar();
        ~Bar();

        //===== SETTER =====//
        void set_text(const std::string &text);
        void set_thintime(ulong thintime);
        void open(ulong total);
        void update();
        void crash() const;
        void close();

    private:
        class Impl; // Pointer to Implementation
        std::unique_ptr<Impl> pImpl;
    };
}

#endif // HOUSEPI__TOOLS__H
