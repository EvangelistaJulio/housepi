/*
 * mcmc_build.h
 *
 * Copyright (C) 2025 Francesco Camaglia, Institut Pasteur
 */

#ifndef HOUSEPI__MCMC__BUILD__H
#define HOUSEPI__MCMC__BUILD__H

#include <memory>
#include <string>
#include <vector>
#include <map>
#include <unordered_set>
#include <algorithm>

#ifdef _HOUSEPI_RCPP
#include <Rcpp.h>
#endif

#include "tools.h"
#include "distribution.h"
#include "variable.h"
#include "household.h"
#include "mcmc.h"

/////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////  GLOBAL VARIABLES  ///////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

// declaration for typed file inputs
const std::string TYPE_PARAM = "P";  // Parameter type declaration
const std::string TYPE_DATA = "D";   // Data type declaration
const std::string TYPE_SCALAR = "S"; // Scalar type declaration

//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//

/* Monte Carlo Markov Chain implementation class. */
class MCMC::Impl
{
public:
    //=======================//
    //===== CONSTRUCTOR =====//
    //=======================//

    Impl() : _uniform_distr_(0., 1.),
             _infec_profile_(std::make_unique<Distribution>()) {}
    ~Impl() = default;

    //===================================//
    //===== RANDOM NUMBER GENERATOR =====//
    //===================================//

    rnd_gen prng;                            // Mersenne Twister 19937 generator (32 or 64 bit)
    std::vector<rnd_gen> prng_parallel = {}; // parallel generators (for multi-threading)

    rnd_gen &get_prng();
    void set_prng_seed();

    rnd_gen &(MCMC::Impl::*_func_prng_)();
    rnd_gen &_prng_unique_thread_() { return prng; }
#ifdef _OPENMP
    rnd_gen &_prng_parallel_threads_() { return prng_parallel[omp_get_thread_num()]; }
#endif

    //----- uniform sampler -----//
    std::uniform_real_distribution<double> _uniform_distr_;
    double uniform();

    //=========================//
    //===== LOGLIKELIHOOD =====//
    //=========================//

    double buffer = NULL_DOUBLE;
    double global_loglike = 0.0;

    double get_global_loglike() const;

    void set_global_loglike(double new_ll);
    void increment_global_loglike(double delta);
    void dgn_glob_ll_stats(double &avg, double &var) const;

    //=================//
    //===== SETUP =====//
    //=================//

    //----- Setup - Attributes -----//
    const str_vec setup_file_headers = {"sample_size", "seed", "n_chains",
                                        "n_steps", "n_adapt", "thin",
                                        "n_burnin", "rand_init", "n_threads",
                                        "inv_temper", "id_trace", "progress"};

    ulong setup_seed;             // random generator seed
    bool setup_rand_init_params;  // random initialization of parameters
    bool setup_progress = true; // if to show progressbar

    ulong setup_n_chains; // n. chains
    ulong setup_n_steps;  // n. iterations
    ulong setup_n_adapt;  // n. (finite) adaptive steps
    ulong setup_thin;     // n. thinning iterations
    ulong setup_n_burnin; // n. burnin iterations

    double setup_inv_T;       // inverse temperature
    unsigned setup_n_threads; // number of OMP threads

    std::string setup_id_trace; // augmented trace to save

    ulong setup_c_sample_size; // single chain sample size
    ulong setup_c_n_sweeps;    // single chain n. of sweeps

    //----- Setup - Setter -----//
    void set_n_threads();
    void set_inver_temper(double inv_T);

    //=====================//
    //===== PARAMETER =====//
    //=====================//

    //----- Parameter - Attributes -----//
    const str_vec param_file_headers = {"name", "value",
                                        "constant", "discrete",
                                        "prior", "arg_name", "arg_value",
                                        "rule", "rule_scale",
                                        "min", "max"};
    std::vector<parID> param_const_list = {};
    std::vector<parID> param_standard_infer_list = {};
    std::vector<parID> param_extra_infer_list = {};

    parIDmap<Variable> param_map = {}; // parameters dictionary
    bool _is_undef_params_ = false;    // undefined parameters valid only if infer

    //----- Parameter - Getter -----//
    double get_param(const parID &p) const;
    const double *get_param_pointer(const parID &p) const;

    //----- Parameter - Setter -----//
    void init_parameters();

    //----- Parameter - Metropolis-Hastings -----//
    bool propose_valid_param_move(const parID &p);
    bool accept_valid_param_move(const parID &p);

    void update_param_accept_rate(const parID &p, bool is_accepted);
    void adapt_param_rule(const parID &p, double ratio);

    void restore_standard_param_value(const parID &p);
    void restore_extra_param_value(const parID &p);

    //==========================//
    //===== HOUSEHOLD DATA =====//
    //==========================//

    //----- Household Data - Attributes -----//
    str_vec data_headers = {};
    size_t n_individuals = 0;
    size_t n_households = 0;
    std::vector<Household> house_vector = {};

    //----- Data - Reports -----//
    std::map<std::string, ulong> _report_indiv_ = {};
    std::map<std::string, ulong> _report_house_ = {};
    std::map<std::string, std::pair<ushort, ulong>> _report_hh_size_ = {};

    void fill_report_data();
    ulong get_report_indiv(const std::string &key) const;
    ulong get_report_house(const std::string &key) const;
    std::pair<ushort, ulong> get_report_hh_size(const std::string &key) const;

    //=====================//
    //===== AUGMENTED =====//
    //=====================//

    //----- Augmented - Attributes -----//
    const str_vec augm_file_headers = {"name", "status",
                                       "value", "value_at",
                                       "constant", "discrete", "init_rand",
                                       "model", "arg_name", "arg_value", "arg_at",
                                       "rule", "rule_scale", "rule_scale_at",
                                       "rule_from", "rule_from_at",
                                       "min", "min_at",
                                       "max", "max_at"};
    std::unordered_set<augmID> augm_list_full = {};
    std::unordered_set<augmID> augm_list_model = {};
    std::unordered_set<augmID> augm_infer_list = {};
    std::map<int, std::unordered_set<augmID>> augm_status_map_model = {};
    std::map<int, std::unordered_set<augmID>> augm_status_map_infer = {};

    //----- Augmented - Setter -----//
    void clear_all_augm(Household &house, indiv i);
    void assign_all_augm(Household &house, indiv i);

    //----- Augmented - State-Checkers -----//
    bool is_augm_any() const;
    bool is_augm_sample() const;
    bool is_augm_model() const;
    bool is_augm_generic() const;

    //===============================//
    //===== INFECTIVITY PROFILE =====//
    //===============================//

    //----- Infectivity Profile - Attributes -----//
    const str_vec infec_profile_file_headers = {"name", "discrete",
                                                "arg_name", "arg_value", "arg_at",
                                                "min", "min_at",
                                                "max", "max_at"};

    std::shared_ptr<Distribution> _infec_profile_ = nullptr;

    //----- Infectivity Profile - Setter -----//
    void switchboard_infec_profile(const std::string &name);
    void assign_infec_profile_args(const str_vec &keys,
                                   const str_vec &vals,
                                   const str_vec &typs);
    void assign_infec_profile_min(const std::string &min_value,
                                  const std::string &min_at);
    void assign_infec_profile_max(const std::string &max_value,
                                  const std::string &max_at);
    void build_infec_profile();

    void fix_param_to_be_inferred_list(const std::string &par_name);

    //===========================//
    //===== IMPORTS/EXPORTS =====//
    //===========================//

    std::map<std::string, std::string> _report_pathtofile_ = {};

#ifdef _HOUSEPI_RCPP

    //----- Rcpp IMPORTS -----//
    void import_setup(const Rcpp::List &input);
    void import_infectivity_profile(const Rcpp::List &input);
    void import_parameter(const Rcpp::DataFrame &input);
    void import_augmented(const Rcpp::DataFrame &input);
    void import_dataset(const Rcpp::DataFrame &input);

    //----- Rcpp EXPORTS -----//

    size_t _aux_row = 0;
    str_vec _aux_data_headers = {};
    Rcpp::CharacterVector _aux_df_index = {};
    Rcpp::CharacterVector _aux_df_house_index = {};
    Rcpp::IntegerVector _aux_df_status = {};
    std::vector<Rcpp::NumericVector> _aux_df_dbl = {};

    void _aux_export_dataset_init_();
    void _aux_export_dataset_update_(const Household &house);
    Rcpp::DataFrame _aux_export_dataset_finalize_();
    Rcpp::DataFrame export_dataset();

    Rcpp::DataFrame export_data_report() const;
    Rcpp::DataFrame export_house_size_report() const;

    Rcpp::List export_setup() const;
    Rcpp::DataFrame export_parameter() const;
    Rcpp::NumericMatrix export_traces() const;
    Rcpp::NumericMatrix export_household_loglike_matrix() const;
    Rcpp::DataFrame export_infectivity_profile(const Rcpp::List &opts,
                                               const Rcpp::List &xlim,
                                               int n_bins) const;

#else
    //----- READ -----//
    void read_setup(const std::string &pathtofile);
    void read_infec_profile(const std::string &pathtofile);
    void read_params(const std::string &pathtofile);
    void read_augm(const std::string &pathtofile);
    void read_data(const std::string &pathtofile);

    //----- SHOW -----//
    void show_summary() const;
    void show_setup() const;
    void show_all_params() const;
    void show_full_build() const;
    void show_results() const;

    //----- Cpp Prints -----//
    std::ofstream aux_data_ofstream;
    void _aux_print_dataset_init_(const std::string &pathtofile);
    void _aux_print_dataset_update_(const Household &house);
    void _aux_print_dataset_finalize_();
    
    void print_infer_state(const std::string &pathtofile) const;
#endif

    //=================//
    //===== BUILD =====//
    //=================//

    void build_setup();
    void build_params();
    void build_augm();
    void build_data(const dict_str_vec &datatable);

    //----- build : auxiliary attributes -----//
    std::map<std::string, std::string> _setup_deflist_ = {};
    std::map<std::string, std::string> _infec_profile_deflist_ = {};
    dict_str_vec _param_deftable_ = {};
    dict_str_vec _augm_deftable_ = {};

    //----- build : type-specified variables -----//
    double type_to_double(const std::string &tkn_val,
                          const std::string &tkn_typ,
                          const Household &house, indiv i);
    void assign_augm_type_args(Household &house, indiv i, augmID a,
                               const str_vec &keys,
                               const str_vec &vals,
                               const str_vec &typs);
    void assign_augm_type_range(Household &house, indiv i, augmID a,
                                const std::string &tkn_min_val,
                                const std::string &tkn_min_typ,
                                const std::string &tkn_max_val,
                                const std::string &tkn_max_typ);

    //----- build : auxiliary methods -----//
    const dict_str_vec _build_augm_name_status_check_();
    void _build_augm_allocate_univocal_defs_(const dict_str_vec &check);
    void _build_augm_clear_houses_and_set_id_trace_();

    //======================//
    //===== METROPOLIS =====//
    //======================//

    ulong _chain_ = 1;   // current chain
    ulong _rung_ = 0;    // current thermodynamic rung
    ulong _n_rungs_ = 0; // number of rungs

    bool is_power_post() const;
    const std::string fmt_progtxt();

    void assert_proper_dataset();
    void assert_proper_mcmc_setup();

    //----- Containers for sample values (aka diagnostics) -----//
    ulong dgn_count = 0;
    dbl_vec dgn_glob_loglike = {};
    size_t dgn_H_star = 0; // index of house to follow up
    indiv dgn_i_star = 0;  // index of individual to follow up
    parIDmap<dbl_vec> dgn_param_value = {};
    augmIDmap<dbl_vec> dgn_augm_value = {};

    void init_move_counts();
    void init_diagnostics();
    void update_diagnostics();

    //----- Metropolis - augmented moves -----//
    double compute_augm_model_loglike(const Household &house) const;

    void (Impl::*_func_augm_model_by_definition_)(const Household &, indiv, double &) const = &Impl::_augm_model_generic_;
    void _augm_model_by_status_(const Household &house, indiv i, double &delta) const;
    void _augm_model_generic_(const Household &house, indiv i, double &delta) const;
};

#endif // HOUSEPI__MCMC__BUILD__H
