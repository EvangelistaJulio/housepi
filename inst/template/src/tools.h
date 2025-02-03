/*
 * tools.h
 *
 * Copyright (C) 2025 Francesco Camaglia, Institut Pasteur
 */

#ifndef HOUSEPI__TOOLS__H
#define HOUSEPI__TOOLS__H

#include <memory>
#include <random>
#include <string>
#include <vector>
#include <map>
#include <cmath>
#include <algorithm>

#ifdef _HOUSEPI_RCPP
#include <Rcpp.h>
#endif

#ifdef _OPENMP
#include <omp.h>
#endif

// UPDATE: Test for pre-compilation
// #define DEBUG_TEMPORARY_MACRO

/////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////  GLOBAL VARIABLES  ///////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

// UPDATE: NaN \ NA have to be revised
const double NULL_DOUBLE = 1e10;        // undefined double
const std::string NULL_STRING = "<NA>"; // undefined string

// constants
const unsigned MAX_NUM_ATTEMPTS = 50;    // max n. of attempts for initialization
const double TARGET_ACCEPT_RATE = 0.234; // target acceptation rate for MCMC

// status aliases
const int INFECTED_STATUS = 1; // infected status (>=)
const int HEALTHY_STATUS = 0;  // healthy status (==)
const int MISSING_STATUS = -1; // missing status (<=)

// group rates
const std::string REFERENCE = "REFERENCE";       // name for reference choice
const std::string CONTACT_RATE = "Contact";      // name for Contact rate
const double REFERENCE_VALUE = 1.0;              // auxiliary reference value for the parameters
const std::string MISSING_RATE = "MISSING_RATE"; // name for Contact rate

// file formatters
const char C_DELIM = ',';      // file delimiter
const char C_LIST_DELIM = '/'; // dictionary delimiter
const char C_DECOR = '-';      // decorator for prints

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
//////////////////////////////////////////  MACROS  /////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

#define UNUSED(x) (void)(x) // MACRO to avoid unused compiler warnings

/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////  TYPEDEF  /////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

using ushort = unsigned short int;
using ulong = unsigned long int;
template <typename T>
using matrix = std::vector<std::vector<T>>;

using dbl_vec = std::vector<double>;
using str_vec = std::vector<std::string>;
using dict_str_vec = std::map<std::string, str_vec>;

// Aliases for improved readibility
using indiv = unsigned short int;

#ifdef DEBUG_TEMPORARY_MACRO
using __access_ID__ = unsigned short int;
const __access_ID__ __undef_ID__ = SIZE_MAX; // value for undefined choice
#else
using __access_ID__ = std::string;
const __access_ID__ __undef_ID__ = NULL_STRING; // name for undefined choice
#endif

using parID = __access_ID__;
using augmID = __access_ID__;
using dataID = __access_ID__;
using rateID = __access_ID__;
using extraID = __access_ID__;

const parID NULL_PARAMETER = __undef_ID__;
const augmID NULL_AUGMENTED = __undef_ID__;
const dataID NULL_DATACOLUMN = __undef_ID__;
const rateID NULL_GROUPRATE = __undef_ID__;
const extraID NULL_EXTRAVARIA = __undef_ID__;

#ifdef DEBUG_TEMPORARY_MACRO
template <typename T>
using parIDmap = std::vector<T>;
template <typename T>
using augmIDmap = std::vector<T>;
template <typename T>
using dataIDmap = std::vector<T>;
template <typename T>
using rateIDmap = std::vector<T>;
template <typename T>
using extraIDmap = std::vector<T>;
#else
#include <unordered_map>
template <typename T>
using parIDmap = std::unordered_map<parID, T>;
template <typename T>
using augmIDmap = std::unordered_map<augmID, T>;
template <typename T>
using dataIDmap = std::unordered_map<dataID, T>;
template <typename T>
using rateIDmap = std::unordered_map<rateID, T>;
template <typename T>
using extraIDmap = std::unordered_map<extraID, T>;
#endif

template <typename T>
using rateIDmapvec = rateIDmap<std::vector<T>>;

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
    void print_decorator(const std::string &title,
                         const std::string &txt = NULL_STRING,
                         char decor = C_DECOR,
                         bool verbose = true);

    // print warnings to stderr
    void warning(const std::string &txt);

    // print error message to stderr
    void error(const std::string &txt);

    // converter functions
    bool stobool(const std::string &tkn);
    int stoi(const std::string &input, const std::string &errtxt = "");
    unsigned long stoul(const std::string &input, const std::string &errtxt = "");
    unsigned long long stoull(const std::string &input, const std::string &errtxt = "");
    double stod(const std::string &input, const std::string &errtxt = "");

    // formatter for the case NA
    double natod(const std::string &input, const std::string &errtxt = "");
    double natod(double x);
    const std::string natostring(double x);

    // argument string parser
    const str_vec parse_list_args(const std::string &tkn,
                                  const char delim = C_LIST_DELIM);

    // read csv files
    void read_csv_list(const std::string &pathtofile,
                       std::map<std::string, std::string> &csv_parsed_list);
    void read_csv_table(const std::string &pathtofile,
                        dict_str_vec &csv_parsed_table,
                        const str_vec &required_headers = {});

    // Welford online algorithm for mean and variance
    void compute_stats(const dbl_vec &data, double &avg, double &var);

    // calculate Levenshtein distance
    ushort levenshtein(const std::string &s1, const std::string &s2);
    // find the closest match
    const std::string find_closest(const std::string &input,
                                   const str_vec &options,
                                   ushort min_dist = 2);
    // check for unused entries from a vector of options and throw exception
    void check_unused_entries(const str_vec &entries, const str_vec &options);
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

namespace progress
{
    /* Class for custom progress bar. */
    class Bar
    {
    public:
        //===== CONSTRUCTOR =====//
        Bar(const std::string &text = "", bool verbose = true);
        ~Bar();

        //===== SETTER =====//
        void set_text(const std::string &text);
        void update(double progress);
        void finalize();

    private:
        // Pointer to Implementation
        class Impl;
        std::unique_ptr<Impl> pImpl;
    };
}

#endif // HOUSEPI__TOOLS__H
