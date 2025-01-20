/*
 * mcmc_build.h
 *
 * Copyright (C) 2025 Francesco Camaglia, Institut Pasteur
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

#include "mcmc.h"
#include "variable.h"
#include "household.h"

//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//

/* Implementation class. */
class MCMC::Impl
{
public:
    //=======================//
    //===== CONSTRUCTOR =====//
    //=======================//

    Impl();
    ~Impl();

    //=========================//
    //===== LOGLIKELIHOOD =====//
    //=========================//

    double buffer = NULL_DOUBLE;
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

    unsigned long long setup_seed = 1;  // random generator seed
    bool setup_rand_init_params = true; // random initialization of parameters

    ulong setup_n_chains = 1;    // n. chains
    ulong setup_n_steps = 10000; // n. iterations
    ulong setup_n_adapt = 0;     // n. (finite) adaptive steps
    ulong setup_thin = 1;        // n. thinning iterations
    ulong setup_n_burnin = 0;    // n. burnin iterations

    double setup_inv_T = 1.0;     // inverse temperature
    unsigned setup_n_threads = 1; // number of OMP threads

    std::string setup_id_trace = NULL_STRING; // augmented trace to save

    ulong setup_c_sample_size; // single chain sample size
    ulong setup_c_n_sweeps;    // single chain n. of sweeps

    void set_n_threads(unsigned n = 1);
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
    std::vector<parID> param_to_be_inferred = {};
    std::vector<parID> param_to_be_constant = {};

    //----- Parameter - Getter -----//
    double get_param(const parID &p) const;
    const double *get_param_pointer(const parID &p) const;

    //----- Parameter - Setter -----//
    void init_parameters();

    //----- Parameter - Metropolis-Hastings -----//
    bool propose_valid_param_move(const parID &p);
    void restore_param_value(const parID &p);
    void update_param_accept_rate(const parID &p, bool is_accepted);
    double param_move_logratio(const parID &p) const;
    void adapt_param_rule(const parID &p, double ratio);

    //==========================//
    //===== HOUSEHOLD DATA =====//
    //==========================//

    //----- Household Data - Attributes -----//
    size_t n_individuals = 0;
    size_t n_households = 0;
    std::vector<Household> house_vector = {};
    str_vec data_headers = {};
    ulong _n_contacts_ = 0;
    std::map<std::string, size_t> _report_data_;
    std::map<std::string, std::string> _report_pathtofile_;

    void update_report_data();
    size_t report_data(const std::string &key) const;

    //=====================//
    //===== AUGMENTED =====//
    //=====================//

    //----- Augmented - Attributes -----//
    const str_vec augm_file_headers = {"name", "status",
                                       "value", "value_type",
                                       "discrete", "init_rand",
                                       "model", "arg_name", "arg_value", "arg_type",
                                       "rule", "rule_scale", "rule_scale_type",
                                       "rule_from", "rule_from_type",
                                       "min", "min_type",
                                       "max", "max_type"};
    std::unordered_set<augmID> augm_list_full = {};
    std::unordered_set<augmID> augm_list_model = {};
    std::unordered_set<augmID> augm_list_to_be_inferred = {};
    std::map<int, std::unordered_set<augmID>> augm_status_map_model = {};
    std::map<int, std::unordered_set<augmID>> augm_status_map_infer = {};

    //----- Augmented - Setter -----//
    void clear_all_augm(Household &house, indiv i);
    void assign_all_augm(Household &house, indiv i);

    //----- Augmented - State-Checkers -----//
    bool is_augm_infer() const;
    bool is_augm_model() const;

    //================//
    //===== READ =====//
    //================//

    void read_setup(const std::string &pathtofile);
    void read_params(const std::string &pathtofile);
    void read_augm(const std::string &pathtofile);

    //================//
    //===== SHOW =====//
    //================//

    void show_setup() const;
    void show_all_params() const;
    void show_results() const;
    void show_summary() const;
    void show_full_build() const;

    //======================//
    //===== Cpp PRINTS =====//
    //======================//

    void print_infer_state(const std::string &pathtofile) const;
    void print_augm_variables(const std::string &pathtofile) const;

#ifdef _HOUSEPI_RCPP
    //========================//
    //===== Rcpp IMPORTS =====//
    //========================//
    void import_setup(const Rcpp::List &x);
    void import_params(const Rcpp::DataFrame &x);
    void import_augm(const Rcpp::DataFrame &x);

    //========================//
    //===== Rcpp EXPORTS =====//
    //========================//
    Rcpp::List export_setup() const;
    Rcpp::DataFrame export_params() const;
    Rcpp::DataFrame export_data() const;
    Rcpp::DataFrame export_traces() const;
    Rcpp::DataFrame export_augm_traces() const;
    Rcpp::NumericMatrix export_household_loglike_traces() const;
#endif

    //===================================//
    //===== RANDOM NUMBER GENERATOR =====//
    //===================================//

    rnd_gen prng;                            // Mersenne Twister 19937 generator (32 or 64 bit)
    std::vector<rnd_gen> prng_parallel = {}; // parallel generators (for multi-threading)

    rnd_gen &get_prng();
    void set_prng_seed(int seed);

    rnd_gen &(MCMC::Impl::*_func_prng_)();
    rnd_gen &_prng_unique_thread_() { return prng; }
#ifdef _OPENMP
    unsigned _max_n_available_threads = 1;
    rnd_gen &_prng_parallel_threads_() { return prng_parallel[omp_get_thread_num()]; }
#endif

    //===========================//
    //===== UNIFORM SAMPLER =====//
    //===========================//

    std::uniform_real_distribution<double> _uniform_distr_;
    double uniform();

    //=================//
    //===== BUILD =====//
    //=================//

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
    void assign_augm_type_args(Household &house, indiv i, augmID a,
                               const std::string &tkn_key,
                               const std::string &tkn_val,
                               const std::string &tkn_typ);
    void assign_augm_type_range(Household &house, indiv i, augmID a,
                                const std::string &tkn_min_val,
                                const std::string &tkn_min_typ,
                                const std::string &tkn_max_val,
                                const std::string &tkn_max_typ);

    //======================//
    //===== PARAMETERS =====//
    //======================//

    parIDmap<Base_Variable> param_map = {}; // parameters dictionary
    bool _is_undef_params_ = false;         // undefined parameters valid only if infer

    //===============================//
	//===== INFECTIVITY PROFILE =====//
    //===============================//

	class InfecProfile;
    class InfecProfile_Ferretti20;
	std::unique_ptr<InfecProfile> _infec_profile_;
	void assign_infectivity_profile(const std::string &name);

    //======================//
    //===== METROPOLIS =====//
    //======================//

	const std::string fmt_progtxt(const std::string &txt, bool is_power_post = false);

	//===== INFER =====//
	ulong _chain_ = 1;			 // current chain
	ulong _tot_thin_sweeps_ = 0; // total number of sweeps
	ulong _rung_ = 0;			 // current thermodynamic rung
	ulong _n_rungs_ = 0;		 // number of rungs

    //----- Containers for sample values (aka diagnostics) -----//
    ulong dgn_count = 0;
    dbl_vec dgn_glob_loglike = {};
    size_t dgn_H_star = 0; // index of house to follow up
    indiv dgn_i_star = 0;  // index of individual to follow up
    parIDmap<dbl_vec> dgn_param_value = {};
    augmIDmap<dbl_vec> dgn_augm_value = {};

    void init_move_counts();
    void init_diagnostics(bool is_power_post);
    void update_diagnostics(bool is_power_post);

    //----- Metropolis - State-Checker -----//
    bool is_augm_generic() const;

    //----- Metropolis - augmented moves -----//
    double compute_augm_model_loglike(const Household &house) const;

    void (Impl::*_func_augm_model_by_definition_)(const Household &, indiv, double &) const = &Impl::_augm_model_generic_;
    void _augm_model_by_status_(const Household &house, indiv i, double &delta) const;
    void _augm_model_generic_(const Household &house, indiv i, double &delta) const;
};

#endif // HOUSEPI__MCMC__BUILD__H
