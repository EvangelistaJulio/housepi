/*
 * mcmc_build.h
 *
 * Copyright (C) 2024 Francesco Camaglia, Institut Pasteur
 */

#ifndef HOUSEPI__MCMC__BUILD__H
#define HOUSEPI__MCMC__BUILD__H

#include <string>
#include <chrono>
#include <vector>
#include <algorithm>
#include <map>
#include <memory>

#ifdef _HOUSEPI_RCPP
#include <Rcpp.h>
#endif

#include "variable.h"
#include "household.h"

//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//

#ifdef _HOUSEPI_RCPP
void rdataframe_to_map(const Rcpp::DataFrame &df, dict_str_vec &output);
#endif

/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////  Parameter : class  ///////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

class Parameter : public Base_Variable
{
public:
    //===== Constructor =====//
    // NOTE using rule_from = INF_DOUBLE, i. e. default move _propose_from_previous_
    Parameter(const std::string &name = UNDEFINED,
              double value = std::numeric_limits<double>::quiet_NaN(),
              bool is_discr = false,
              const std::string &prior_name = UNDEFINED,
              dict_args prior_args = {},
              dict_args range = {},
              const std::string &rule_name = UNDEFINED,
              double rule_scale = INF_DOUBLE) : Base_Variable(name, value, is_discr)
    {
        assign_distribution(prior_name, prior_args, range);
        assign_proposer(rule_name, rule_scale);
    }

    //----- Diagnostics -----//
    ulong dgn_count = 0;
    std::vector<double> dgn_value = {};

    void init_diagnostics(ulong sample_size);
    void update_diagnostics();
    std::string show_diagnostics() const;
};

//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//

/* Implementation class. */
class Monte_Carlo_Build
{
public:
    //=========================//
    //===== LOGLIKELIHOOD =====//
    //=========================//

    double buffer = INF_DOUBLE;
    double global_loglike = 0.0;

    double get_global_loglike() const;
    double get_delta_global_loglike() const;

    void set_global_loglike(double new_ll);
    void restore_global_loglike();
    void increment_global_loglike(double delta);
    void dgn_glob_ll_stats(double &avg, double &var) const;

    //=================//
    //===== SETUP =====//
    //=================//

    unsigned long long setup_seed = 1; // random generator seed
    bool setup_rnd_init_params = true; // random initialization of parameters

    ulong setup_n_chains = 1;    // n. chains
    ulong setup_n_steps = 10000; // n. iterations
    ulong setup_n_adapt = 0;     // n. (finite) adaptive steps
    ulong setup_thin = 1;        // n. thinning iterations
    ulong setup_n_burnin = 0;    // n. burnin iterations

    double setup_inv_T = 1.0;     // inverse temperature
    unsigned setup_n_threads = 1; // number of OMP threads

    ulong setup_c_sample_size; // single chain sample size
    ulong setup_c_n_sweeps;    // single chain n. of sweeps

    //=====================//
    //===== PARAMETER =====//
    //=====================//

    //----- Parameter - Attributes -----//
    str_vec param_file_headers = {"name", "value",
                                  "infer", "discrete",
                                  "prior", "args", "argsval",
                                  "rule", "rule_scale",
                                  "min", "max"};
    std::vector<parID> param_to_be_inferred = {};
    std::vector<parID> param_to_be_static = {};

    //----- Parameter - Getter -----//
    bool assert_param(const parID &p) const;
    double get_param(const parID &p) const;
    const double *get_param_pointer(const parID &p) const;
    bool get_infer(const parID &p) const;

    //----- Parameter - Setter -----//
    void assign_param(const parID &p, const Parameter &new_param);
    void init_rand_params();
    const double *assign_scalar_pointer(double val);

    //----- Parameter - Metropolis-Hastings -----//
    bool propose_valid_param_move(const parID &p);
    void restore_param_value(const parID &p);
    void update_param_accept_rate(const parID &p, bool is_accepted);
    double param_move_logratio(const parID &p) const;
    void adapt_param_rule(const parID &p, double ratio);

    //===== HOUSEHOLD DATA : Attributes =====//
    size_t n_individuals = 0;
    size_t n_households = 0;
    std::vector<Household> house_vector = {};

    //----- Data Reports : Attributes -----//
    str_vec data_headers = {};
    ulong _n_contacts_ = 0;
    std::map<std::string, size_t> _report_data_;
    std::map<std::string, std::string> _report_pathtofile_;

    //===== AUGMENTED : Attributes =====//
    std::unordered_set<augmID> augm_list_full = {};
    std::unordered_set<augmID> augm_list_infer = {};
    std::unordered_set<augmID> augm_list_model = {};
    std::map<int, std::unordered_set<augmID>> augm_status_map_infer = {};
    std::map<int, std::unordered_set<augmID>> augm_status_map_model = {};
    bool is_augm_infer = false;
    bool is_augm_model = false;

    //===== CONSTRUCTOR =====//
    Monte_Carlo_Build() : _uniform_distr_(0., 1.)
    {
#ifdef _OPENMP
        _max_n_available_threads = omp_get_max_threads();
#endif
    }
    ~Monte_Carlo_Build() {};

    //===== READ =====//
    void read_setup(const std::string &pathtofile);
    void read_params(const std::string &pathtofile);
    void read_augm(const std::string &pathtofile);

    //===== SHOW =====//
    void show_setup() const;
    std::string show_param(const parID &p) const;
    void show_all_params() const;
    void show_results() const;
    void show_summary() const;
    void show_full_build() const;

    //===== Cpp PRINTS =====//
    void print_infer_state(const std::string &pathtofile) const;
    void print_augm_variables(const std::string &pathtofile) const;

#ifdef _HOUSEPI_RCPP
    //===== Rcpp IMPORTS =====//
    void import_setup(const Rcpp::List &x);
    void import_params(const Rcpp::DataFrame &x);
    void import_augm(const Rcpp::DataFrame &x);

    //===== Rcpp EXPORTS =====//
    Rcpp::List export_setup() const;
    Rcpp::DataFrame export_params() const;
    Rcpp::DataFrame export_data() const;
    Rcpp::DataFrame export_traces() const;
    Rcpp::NumericMatrix export_household_loglike_traces() const;
#endif

    //===== SETUP =====//
    void set_n_threads(unsigned n = 1);

    //===== DIAGNOSTICS =====//
    ulong dgn_count = 0;
    std::vector<double> dgn_glob_loglike = {};

    //===== SETUP =====//
    void assign_house(const std::string &ID);
    void update_report_data();
    size_t report_data(const std::string &key) const;

    //===== AUGMENTED =====//
    void clear_all_augm(Household &house, indiv i);
    void assign_all_augm(Household &house, indiv i);

    void disable_assert();
    void init_move_counts();
    void init_diagnostics(bool is_power_post);
    void update_diagnostics(double glob_loglike, bool is_power_post);

    //----- Metropolis : augmented moves -----//
    double compute_augm_model_loglike(const Household &house) const;

    double (Monte_Carlo_Build::*_func_augm_model_loglike_)(const Household &) const = &Monte_Carlo_Build::_no_augm_model_loglike_;
    double _no_augm_model_loglike_(const Household &house) const;
    double _augm_model_loglike_(const Household &house) const;

    void (Monte_Carlo_Build::*_func_augm_model_by_definition_)(const Household &, indiv, double &) const = &Monte_Carlo_Build::_augm_model_by_status_;
    void _augm_model_by_status_(const Household &house, indiv i, double &delta) const;
    void _augm_model_generic_(const Household &house, indiv i, double &delta) const;

    //===== PARAMETERS =====//
    parIDmap<Parameter> param_map = {};                              // parameters dictionary
    std::unordered_map<double, Parameter> _scalar_type_params_ = {}; // scalar types dictionary // FIXME : temporary

    bool _is_undef_params_ = false; // undefined parameters valid only if infer

    //----- get param (functional methods) -----//
    double (Monte_Carlo_Build::*_func_get_param_)(const parID &) const = &Monte_Carlo_Build::_get_param_assert_;
    double _get_param_default_(const parID &p) const;
    double _get_param_assert_(const parID &p) const;

    //----- get param pointer (functional methods) -----//
    const double *(Monte_Carlo_Build::*_func_get_param_pointer_)(const parID &) const = &Monte_Carlo_Build::_get_param_pointer_assert_;
    const double *_get_param_pointer_default_(const parID &p) const;
    const double *_get_param_pointer_assert_(const parID &p) const;

    //===== RANDOM NUMBER GENERATOR =====//
    rnd_gen prng;                            // Mersenne Twister 19937 generator (32 or 64 bit)
    std::vector<rnd_gen> prng_parallel = {}; // parallel generators (for multi-threading)

    rnd_gen &get_prng();
    void set_prng_seed(int seed);

    rnd_gen &(Monte_Carlo_Build::*_func_prng_)();
    rnd_gen &_prng_unique_thread_() { return prng; }
#ifdef _OPENMP
    unsigned _max_n_available_threads = 1;
    rnd_gen &_prng_parallel_threads_() { return prng_parallel[omp_get_thread_num()]; }
#endif

    //===== UNIFORM SAMPLER =====//
    std::uniform_real_distribution<double> _uniform_distr_;
    double uniform();

    //===== BUILD =====//
    void build_augm();
    void build_setup();
    void build_params();

    //----- build : auxiliary attributes -----//
    std::map<std::string, std::string> _setup_deflist_ = {};
    dict_str_vec _param_deftable_ = {};
    dict_str_vec _augm_deftable_ = {};

    //----- build : type-specified variables -----//
    double type_to_double(const std::string &tkn_val,
                          const std::string &tkn_typ,
                          const Household &house, indiv i);
    const double *type_to_doubleptr(const std::string &tkn_val,
                                    const std::string &tkn_typ,
                                    const Household &house, indiv i);
    dict_args typedict_to_sdargs(const std::string &tkn_key,
                                 const std::string &tkn_val,
                                 const std::string &tkn_typ,
                                 const Household &house, indiv i);
    dict_args_ptr typedict_to_sdargs_ptr(const std::string &tkn_key,
                                         const std::string &tkn_val,
                                         const std::string &tkn_typ,
                                         const Household &house, indiv i);
};

#endif // HOUSEPI__MCMC__BUILD__H
