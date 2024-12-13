/*
 * tools.h
 *
 * Copyright (C) 2024 Francesco Camaglia, Institut Pasteur
 */

#ifndef HOUSEPI__TOOLS__H
#define HOUSEPI__TOOLS__H

#include <memory>
#include <cmath>
#include <random>
#include <algorithm>
#include <string>
#include <vector>
#include <map>
#include <unordered_map>
#include <unordered_set>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <chrono>

#ifdef _HOUSEPI_RCPP
#include <Rcpp.h>
#else
#include <iostream>
#endif

#ifdef _OPENMP
#include <omp.h>
#endif

// Test for pre-compilation
//#define DEBUG_TEMPORARY_MACRO

/////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////  GLOBAL VARIABLES  ///////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

// FIXME :NaN and NA to be revised
constexpr double INF_DOUBLE = 1e10;   // numerical infinity value
const std::string UNDEFINED = "<NA>"; // value for undefined choice

// constants
constexpr int MAX_NUMBER_ATTEMPTS = 50;      // max n. of attempts for initialization
constexpr double TARGET_ACCEPT_RATE = 0.234; // target acceptation rate for MCMC

// status aliases
constexpr int INFECTIOUS = 1;      // status of infectious (>=)
constexpr int SUSCEPTIBLE = 0;     // status of susceptible (==)
constexpr int MISSING_STATUS = -1; // missing status (==)

// group rates
const std::string REFERENCE = "REFERENCE";       // name for reference choice
const std::string CONTACT_RATE = "Contact";      // name for Contact rate
constexpr double REFERENCE_VALUE = 1.0;          // auxiliary reference value for the parameters
const std::string MISSING_RATE = "MISSING_RATE"; // name for Contact rate

// file formatters
constexpr char C_DELIM = ',';      // file delimiter
constexpr char C_LIST_DELIM = '/'; // dictionary delimiter
constexpr char C_DECOR = '-';      // decorator for prints

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

using dict_args = std::unordered_map<std::string, double>;
using dict_args_ptr = std::unordered_map<std::string, const double *>;

using str_vec = std::vector<std::string>;
using dict_str_vec = std::map<std::string, str_vec>;

// Aliases for improved readibility
using indiv = unsigned short int;

#ifdef DEBUG_TEMPORARY_MACRO
using __access_ID__ = size_t;
const __access_ID__ __undef_ID__ = SIZE_MAX; // value for undefined choice
#else
using __access_ID__ = std::string;
const __access_ID__ __undef_ID__ = UNDEFINED; // name for undefined choice
#endif

using parID = __access_ID__;
using augmID = __access_ID__;
using dataID = __access_ID__;
using rateID = __access_ID__;
using extraID = __access_ID__;

const parID UNDEF_PARAMETER = __undef_ID__;
const parID UNDEF_AUGMENTED = __undef_ID__;
const parID UNDEF_DATACOLUMN = __undef_ID__;
const parID UNDEF_GROUPRATE = __undef_ID__;
const parID UNDEF_EXTRAVARIA = __undef_ID__;

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

    // output formatter
    const std::string percentage(double fraction, ushort precision = 1);
    const std::string truth(bool condition);

    // stodut print
    void print_stdout(const std::string &txt);
    void print_stderr(const std::string &txt);
    void print_decorator(const std::string &title,
                         const std::string &txt = UNDEFINED,
                         char decor = C_DECOR,
                         bool verbose = true);

    // print warnings to stderr
    void warning(const std::string &txt);

    // print error message to stderr
    void error(const std::string &txt);

    // converter functions
    int stoi(const std::string &input, const std::string &errtxt = "");
    double stod(const std::string &input, const std::string &errtxt = "");
    unsigned long stoul(const std::string &input, const std::string &errtxt = "");
    unsigned long long stoull(const std::string &input, const std::string &errtxt = "");
    bool stobool(const std::string &tkn);
    const std::string strcleaner(std::string tkn);

    void read_csv_table(const std::string &pathtofile,
                        dict_str_vec &csv_parsed_table,
                        const str_vec &required_headers = {});
    const str_vec parse_list_args(const std::string &tkn,
                                  const char delim = C_LIST_DELIM);

    // Calculate Levenshtein distance
    ushort levenshtein(const std::string &s1, const std::string &s2);
    // Find the closest match
    const std::string find_closest(const std::string &input,
                                   const str_vec &options,
                                   ushort min_dist = 2);

}

//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//

/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////  Progress Bar Class  //////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

/* Forward declaration of implementation. */
class CustomProgressBarImpl;

/* Class for custom progress bar. */
class CustomProgressBar
{
public:
    CustomProgressBar(bool verbose = true);
    ~CustomProgressBar();

    void set_text(const std::string &text);
    void initialize(const std::string &text = "");

    void update(double progress);
    void finalize();

private:
    //===== Pointer to Implementation =====//
    std::unique_ptr<CustomProgressBarImpl> pImpl;
};

#endif // HOUSEPI__TOOLS__H
