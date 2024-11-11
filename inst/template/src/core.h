/*
 * core.cpp
 *
 * Copyright (C) 2024 Francesco Camaglia, Institut Pasteur
 */

#include <string>
#include <chrono>
#include <vector>
#include <algorithm>
#include <map>
#include <memory>

#include "varia.h"

//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//

/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////  Parameter : class  ///////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

class Parameter : public Base_Variable
{
public:
    //===== Constructor =====//
    // NOTE using rule_from = INF_DOUBLE, i. e. default move _propose_from_previous_
    Parameter(const std::string &param_name = UNDEFINED, double value = std::numeric_limits<double>::quiet_NaN(),
              const std::string &prior_name = UNDEFINED, dictargs prior_args = {}, dictargs range = {},
              const std::string &rule_name = UNDEFINED, dictargs rule_args = {}) : Base_Variable(param_name, value)
    {
        assign_distribution(prior_name, prior_args, range);
        assign_proposer(rule_name, rule_args);
    }

    //===== Diagnostics =====//
    ulong dgn_count = 0;
    std::vector<double> dgn_value = {};
    void init_diagnostics(ulong sample_size);
    void update_diagnostics();
    std::string show_diagnostics() const;
    const std::vector<double> &get_dgn_value() const { return dgn_value; }
    double get_dgn_value(size_t i) const { return dgn_value[i]; }
};

//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//

/////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////  Monte_Carlo_Core : class  //////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

/* This is a wrapper class containing aliases for the model parameters. */
class Monte_Carlo_Core
{
public:
    //===== Parameter Attributes =====//
    std::vector<std::string> param_file_headers = {"name", "value", "infer",
                                                   "prior", "args", "argsval",
                                                   "rule", "rule_args", "rule_argsval",
                                                   "min", "max"};
    std::vector<std::string> param_to_be_inferred = {};
    std::vector<std::string> param_to_be_static = {};

    //===== Setup Attributes =====//
    unsigned long long setup_seed = 1; // random generator seed
    bool setup_rnd_init_params = true; // random initialization of parameters

    ulong setup_num_chains = 1;   // n. chains
    ulong setup_num_iter = 10000; // n. iterations
    ulong setup_num_adapt = 0;    // n. (finite) adaptive steps
    ulong setup_thinning = 1;     // n. thinning iterations
    ulong setup_burn_in = 0;      // n. burnin iterations

    double setup_inv_T = 1.0;       // inverse temperature
    unsigned setup_num_threads = 1; // number of OMP threads

    ulong setup_c_sample_size; // single chain sample size
    ulong setup_c_num_sweeps;  // single chain n. of sweeps

    //===== Diagnostics Attributes=====//
    ulong dgn_count = 0;
    std::vector<double> dgn_glob_loglike = {};

    //===== CONSTRUCTOR =====//
    Monte_Carlo_Core() : _uniform_distr_(0., 1.)
    {
#ifdef _OPENMP
        _max_num_available_threads = omp_get_max_threads();
#endif
    }
    ~Monte_Carlo_Core() {};

    //===== READ =====//
    void read_setup(const std::string &pathtofile = UNDEFINED);
    void read_params(const std::string &pathtofile = UNDEFINED);

    //===== SETUP : setter =====//
    void set_num_threads(unsigned n = 1);

    //===== PARAM : getter =====//
    void disable_assert_param();
    void assert_param(const std::string &parID) const;
    double get_param(const std::string &parID) const { return (this->*_func_get_param_)(parID); }
    const double *get_param_pointer(const std::string &parID) const { return (this->*_func_get_param_pointer_)(parID); }
    bool get_infer(const std::string &parID) const { return params.at(parID).get_infer(); }
    double get_accept_rate(const std::string &parID) const { return params.at(parID).get_accept_rate(); }

    //===== PARAM SETTER =====//
    void init_rand_params();
    void reset_accept_rate(const std::string &parID) { params.at(parID).reset_accept_rate(); }
    const double *assign_scalar_pointer(double val);

    //===== Metropolis-Hastings moves =====//
    bool propose_valid_param_move(const std::string &parID);
    void restore_param_value(const std::string &parID) { params.at(parID).restore_value(); }
    void update_param_accept_rate(const std::string &parID, bool is_accepted) { params.at(parID).update_accept_rate(is_accepted); }
    double param_move_logratio(const std::string &parID) const { return params.at(parID).get_move_logratio(); }
    void adapt_param_rule(const std::string &parID, double ratio) { params.at(parID).adapt_rule(ratio); }
    void init_param_diagnostics(ulong sample_size);
    void update_param_diagnostics();
    void init_param_move_counts();

    void dgn_glob_ll_stats(double &avg, double &var) const;

    //===== Pseudo-Random Number Generator =====//
    sys_rnd_gen &get_prng() { return (this->*_func_prng_)(); }
    double uniform() { return _uniform_distr_(get_prng()); }

    //===== Show =====//
    std::string show_param(const std::string &parID) const { return parID + C_DELIM + params.at(parID).show(); }

#ifdef USE_RCPP
    //===== Rcpp IMPORTS =====//
    void import_setup(Rcpp::List setup);
    void import_params(Rcpp::DataFrame params);

    //===== Rcpp EXPORTS =====//
    Rcpp::List export_setup() const;
    Rcpp::DataFrame export_params() const;
    Rcpp::DataFrame export_traces() const;
#else
    void show_setup() const;
    void show_all_params() const;
    void show_results() const;
    void print_infer_state(const std::string &output_filename) const;
#endif

    //===== Build : core functions =====//
    void build_setup();
    void build_params();

private:
    //===== Internal Attributes =====//
    std::unordered_map<std::string, Parameter> params = {};          // parameters dictionary
    std::uniform_real_distribution<double> _uniform_distr_;          // uniform distribution
    std::unordered_map<double, Parameter> _scalar_type_params_ = {}; // scalar types dictionary // FIXME : temporary

    bool _is_undef_params_ = false; // undefined parameters valid only if infer

    //===== get param (functional methods) =====//
    double (Monte_Carlo_Core::*_func_get_param_)(const std::string &) const = &Monte_Carlo_Core::_get_param_assert_;
    double _get_param_default_(const std::string &parID) const { return params.at(parID).get_value(); };
    double _get_param_assert_(const std::string &parID) const;

    //===== get param pointer (functional methods) =====//
    const double *(Monte_Carlo_Core::*_func_get_param_pointer_)(const std::string &) const = &Monte_Carlo_Core::_get_param_pointer_assert_;
    const double *_get_param_pointer_default_(const std::string &parID) const { return params.at(parID).get_pointer(); }
    const double *_get_param_pointer_assert_(const std::string &parID) const;

    //===== Pseudo Random Number Generator (auxiliary methods) =====//
    // FIXME : this should change if Rcpp according to developer guides
    sys_rnd_gen prng;                            // Mersenne Twister 19937 generator (32 or 64 bit)
    std::vector<sys_rnd_gen> prng_parallel = {}; // parallel generators (for multi-threading)
    sys_rnd_gen &(Monte_Carlo_Core::*_func_prng_)();
    sys_rnd_gen &_prng_unique_thread_() { return prng; }
#ifdef _OPENMP
    unsigned _max_num_available_threads = 1;
    sys_rnd_gen &_prng_parallel_threads_() { return prng_parallel[omp_get_thread_num()]; }
#endif
    void set_prng_seed(int seed);

    //----- build : auxiliary attributes -----//
    std::map<std::string, std::string> _setup_deflist_ = {};
    std::map<std::string, std::vector<std::string>> _param_deftable_ = {};

    //===== PARAM SETTER =====//
    void assign_param(const std::string &parID, const Parameter &new_param);
};
