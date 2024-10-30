/*
 * tools.h
 *
 * Copyright (C) 2024 Francesco Camaglia, Institut Pasteur
 */

#ifndef UTILS__H
#define UTILS__H

#include <memory>
#include <cmath>
#include <random>
#include <algorithm>
#include <string>
#include <vector>
#include <set>
#include <unordered_set>
#include <map>
#include <unordered_map>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <chrono>

#ifdef USE_RCPP
#include <Rcpp.h>
#else
#include <iostream>
#endif

#ifdef _OPENMP
#include <omp.h>
#endif

/////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////  32bit/64bit  ///////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

/* Use pseudo-random number generator for 64bit architecture. */
#if __x86_64__ || __ppc64__ || __aarch64__ || _WIN64
using sys_rnd_gen = std::mt19937_64;
#else
using sys_rnd_gen = std::mt19937;
#endif

/////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////  MACROS  /////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

#define UNUSED(x) (void)(x) // MACRO to avoid unused compiler warnings

/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////  TYPEDEF  /////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

typedef unsigned short int ushort;
typedef unsigned long int ulong;
typedef std::vector<std::vector<double>> squarematrix;
typedef std::unordered_map<std::string, double> dictargs;
typedef std::unordered_map<std::string, const double *> dictargs_ptr;

/////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////  GLOBAL VARIABLES  ///////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

// FIXME :NaN and NA to be revised
constexpr double INF_DOUBLE = 1e10;   // numerical infinity value
const std::string UNDEFINED = "<NA>"; // name for undefined choice

// constants
constexpr int MAX_NUMBER_ATTEMPTS = 50;            // max n. of attempts for initialization
constexpr double TARGET_ACCEPT_RATE = 0.234;       // target acceptation rate for MCMC
const std::string IO_MISMATCH = "__io_mismatch__"; // mismatch in type/key/val dictionary

// group rates
const std::string REFERENCE = "__reference__"; // name for reference choice
const std::string CONTACT_RATE = "Contact";    // name for Contact rate
constexpr double REFERENCE_VALUE = 1.0;        // auxiliary reference value for the parameters
const std::string EMPTY_RATE = "__empty__";    // name for Contact rate

// file formatters
constexpr char C_DELIM = ',';      // file delimiter
constexpr char C_LIST_DELIM = '/'; // dictionary delimiter
constexpr char C_DECOR = '-';      // decorator for prints
constexpr char C_PROG = '=';       // decorator for progress bar

//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//

/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////  Input/Output Functions  //////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

namespace tools
{
    // internal undefined string values
    bool is_na(const std::string &txt);

    // output formatter
    std::string percent(double num, double denom);
    std::string truth(bool condition);

    // stodut print
    void print_stdout(const std::string &txt);
    void print_decorator(const std::string &title,
                         const std::string &txt = UNDEFINED,
                         char decor = C_DECOR,
                         bool verbose = true);

    // print progress bar to stdout
    void progress_bar_init(const std::string &txt, bool verbose = true);
    void progress_bar_run(const std::string &txt, ulong progress, ulong total, bool verbose = true);
    void progress_bar_end(const std::string &txt, bool verbose = true);

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
    std::string strcleaner(std::string tkn);

#ifdef USE_RCPP
    void rdataframe_to_map(Rcpp::DataFrame df, std::map<std::string, std::vector<std::string>> &output);
#endif

    void read_csv_table(const std::string &pathtofile,
                        std::map<std::string, std::vector<std::string>> &csv_parsed_table,
                        const std::vector<std::string> &required_headers = {});
    std::vector<std::string> parse_list_args(const std::string &tkn, const char delim = C_LIST_DELIM);
}

#endif // UTILS__H