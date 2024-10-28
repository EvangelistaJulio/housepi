/*
 * utils.h
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
typedef std::unordered_map<std::string, double> dictargs;
typedef std::unordered_map<std::string, const double *> dictargs_ptr;
typedef std::vector<std::vector<double>> squarematrix;

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
    void rlist_to_map(Rcpp::List list, std::map<std::string, std::string> &output);
    void rdataframe_to_map(Rcpp::DataFrame df, std::map<std::string, std::vector<std::string>> &output);
#endif

    std::vector<std::string> read_column_headers(std::ifstream &input, const std::vector<std::string> &required_headers);
    bool read_row_by_headers(std::ifstream &input, const std::vector<std::string> &headers, std::map<std::string, std::string> &this_row);
    void read_csv_list(const std::string &pathtofile, std::map<std::string, std::string> &csv_parsed_list);
    void read_csv_table(const std::string &pathtofile,
                        std::map<std::string, std::vector<std::string>> &csv_parsed_table,
                        const std::vector<std::string> &required_headers = {});

    std::vector<std::string> parse_list_args(const std::string &tkn, const char delim = C_LIST_DELIM);
    dictargs merge_list_args(const std::string &tkn_keys, const std::string &tkn_vals);
}

//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//

/////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////  Base Function : class //////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

// FIXME :
// for the moment it is not possible to provide default pointer argumets
// a solution could be point to internal, at that point the whole pointer thing could be collapsed

/* Parent class to all function classes where the arguements are stored in arrays. */
class Base_Function
{
public:
    //=====Constructors=====//
    // NOTE : to provide a default order is the only way  I found to enforce an order
    // this is necessary if I want to access arguments by index and not by name (faster...?)
    // FIXME : 1) check if by index faster; 2) collapse overloaded constructors
    Base_Function() {};
    Base_Function(const std::string &name,
                  const dictargs &args_val,
                  const std::vector<std::string> &order,
                  const dictargs &domain = {});
    Base_Function(const std::string &name,
                  const dictargs_ptr &args_ptr,
                  const std::vector<std::string> &order,
                  const dictargs &domain = {});
    virtual ~Base_Function() {};

    //===== Getter (functional) =====//
    std::string get_name() const { return _name_; }
    bool valid(double x) const { return (this->*_valid_range_method_)(x); }
    double get_arg(int i) const { return (this->*_get_arg_method_)(i); }
    std::string get_arg_name(int i) const { return _args_names_[i]; }

    //===== Setter (for internal arguments) =====//
    void set_arg(int i, double new_value) { (this->*_set_arg_method_)(i, new_value); }

    //===== Show =====//
    std::string show_arg_names() const;
    std::string show_arg_values() const;
    std::string show_range_min() const { return (_range_min_ > -INF_DOUBLE) ? std::to_string(_range_min_) : UNDEFINED; }
    std::string show_range_max() const { return (_range_max_ < INF_DOUBLE) ? std::to_string(_range_min_) : UNDEFINED; }

protected:
    //===== Initial non default arguments =====//
    void assign_args(dictargs new_args_val);
    void assign_args(dictargs_ptr new_args_ptr);

    int set_range(const dictargs &new_range);
    double get_range_min() const { return _range_min_; }
    double get_range_max() const { return _range_max_; }

private:
    //===== Arguments =====//
    std::string _name_ = UNDEFINED;
    std::vector<std::string> _args_names_ = {}; // ordered list of argument names
    double *_args_val_ = nullptr;               // internal arguments (fixed at initialization)
    const double **_args_ptr_ = nullptr;        // external arguments (pointer to variable)

    //===== Argument getter : functional methods =====//
    double (Base_Function::*_get_arg_method_)(int) const;
    double _get_arg_val_(int i) const { return _args_val_[i]; }
    double _get_arg_ptr_(int i) const { return *(_args_ptr_[i]); }

    //===== Argument setter : functional methods =====//
    // note : the set_arg is only used for the adaptive
    // for this reason a pointer argument can't be set, but just assigned
    void (Base_Function::*_set_arg_method_)(int, double);
    void _set_arg_val_(int i, double new_value) { _args_val_[i] = new_value; }
    void _set_arg_ptr_(int i, double new_value)
    {
        UNUSED(new_value);
        std::ostringstream oss;
        oss << "Inconsistent change of pointer argument'" << _args_names_[i]
            << "' within '" << get_name() << "'.";
        tools::error(oss.str());
    }

    //===== Valid : functional methods =====//
    double _range_min_ = -INF_DOUBLE; // init to -infty
    double _range_max_ = INF_DOUBLE;  // init to +infty

    // alternative valid range methods
    bool (Base_Function::*_valid_range_method_)(double) const = &Base_Function::_valid_range_always_;
    bool _valid_range_always_(double x) const
    {
        UNUSED(x);
        return true;
    }
    bool _valid_range_min_(double x) const { return (x > _range_min_); }
    bool _valid_range_max_(double x) const { return (x < _range_max_); }
    bool _valid_range_both_(double x) const { return _valid_range_min_(x) && _valid_range_max_(x); }
};

#endif // UTILS__H