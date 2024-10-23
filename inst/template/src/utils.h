/*
 * utils.h
 *
 * Copyright (C) 2024 Francesco Camaglia, Institut Pasteur
 */

#ifndef UTILS__H
#define UTILS__H

#include <random>
#include <algorithm>
#include <string>
#include <vector>
#include <set>
#include <unordered_map>
#include <cmath>
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
constexpr double INF_DOUBLE = 1e10;         // numerical infinity value
const std::string UNDEFINED = "<NA>";       // name for undefined choice
const std::string EMPTY_RATE = "__empty__"; // name for Contact rate

// constants
constexpr int MAX_NUMBER_ATTEMPTS = 50;            // max n. of attempts for initialization
constexpr double TARGET_ACCEPT_RATE = 0.234;       // target acceptation rate for MCMC
const std::string IO_MISMATCH = "__io_mismatch__"; // name for undefined choice
constexpr double REFERENCE_VALUE = 1.0;            // auxiliary reference value for the parameters
const std::string REFERENCE = "__reference__";     // name for reference choice
const std::string CONTACT_RATE = "Contact";        // name for Contact rate

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
    constexpr char C_DELIM = ',';      // file delimiter
    constexpr char C_LIST_DELIM = '/'; // dictionary delimiter
    constexpr char C_DECOR = '-';      // decorator for prints
    constexpr char C_PROG = '=';       // decorator for progress bar

    bool is_na(const std::string &txt);
    void statistics(const std::vector<double> &data, double &avg, double &var);

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
    double get_range_min() const { return (_range_min_ > -INF_DOUBLE) ? _range_min_ : std::numeric_limits<double>::quiet_NaN(); }
    double get_range_max() const { return (_range_max_ < INF_DOUBLE) ? _range_max_ : std::numeric_limits<double>::quiet_NaN(); }
    double get_arg(int i) const { return (this->*_get_arg_method_)(i); }
    std::string get_arg_name(int i) const { return _args_names_[i]; }

    //===== Setter (for internal arguments) =====//
    void set_arg(int i, double new_value) { (this->*_set_arg_method_)(i, new_value); }

    //===== Show =====//
    std::string show_arg_names() const;
    std::string show_arg_values() const;

protected:
    //===== Initial non default arguments =====//
    void assign_args(dictargs new_args_val);
    void assign_args(dictargs_ptr new_args_ptr);

    std::string _name_ = UNDEFINED;
    std::vector<std::string> _args_names_ = {}; // ordered list of argument names
    double *_args_val_ = nullptr;               // internal arguments (fixed at initialization)
    const double **_args_ptr_ = nullptr;        // external arguments (pointer to variable)

    double (Base_Function::*_get_arg_method_)(int) const;
    double _get_arg_val_(int i) const { return _args_val_[i]; }
    double _get_arg_ptr_(int i) const { return *(_args_ptr_[i]); }

    //===== Valid : functional methods =====//
    double _range_min_ = -INF_DOUBLE; // init to -infty
    double _range_max_ = INF_DOUBLE;  // init to +infty
    int set_range(const dictargs &new_range);

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
};

#endif // UTILS__H