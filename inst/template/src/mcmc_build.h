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

    Impl() : _uniform_distr_(0., 1.), _infec_(*this) {}
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
    const std::pair<double, double> dgn_glob_ll_stats() const;

    //=================//
    //===== SETUP =====//
    //=================//

    //----- Setup - Attributes -----//
    const str_vec setup_file_headers = {"sample_size", "seed", "n_chains",
                                        "n_steps", "n_adapt", "thin",
                                        "n_burnin", "rand_init", "n_threads",
                                        "inv_temper", "id_trace", "progress"};

    ulong setup_seed;            // random generator seed
    bool setup_rand_init_params; // random initialization of parameters
    bool setup_progress = true;  // if to show progressbar

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

    void set_n_threads();

    //=====================//
    //===== PARAMETER =====//
    //=====================//

    //----- Parameter - Attributes -----//
    const str_vec param_file_headers = {"name", "value",
                                        "constant", "discrete",
                                        "prior", "arg_name", "arg_value",
                                        "min", "max",
                                        "rule", "scale"};
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

    void restore_parameter(const parID &p, bool extra_depen);

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
                                       "min", "min_at",
                                       "max", "max_at",
                                       "rule", "scale", "scale_at",
                                       "from", "from_at"};
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

    void fix_param_to_be_inferred_list(const std::string &par_name);

    class InfecProfile
    {
    public:
        std::map<std::string, std::string> deflist = {};
        std::shared_ptr<Distribution> profile = nullptr;

        //===== CONSTRUCTOR =====//
        InfecProfile(MCMC::Impl &ptr) : profile(std::make_unique<Distribution>()),
                                        pImpl(ptr) {}
        ~InfecProfile() = default;

        //===== SETTER =====//
        void build();
        void init_lookup();
        void restore_lookup();
        void update_lookup();
        void assign_template(std::map<std::string, std::string> &x);

        //===== GETTER =====//
        bool is_lookup() const;
        double inst(double t, const dbl_vec &opts) const;
        double cumul(double t, const dbl_vec &opts) const;

        void assert_valid_time(double t) const;

    private:
        MCMC::Impl &pImpl; // pointer to external class
        const str_vec headers = {"name", "discrete",
                                 "arg_name", "arg_value", "arg_at",
                                 "min", "min_at",
                                 "max", "max_at"};
        dbl_vec lookup_inst = {};
        dbl_vec lookup_cumul = {};
        dbl_vec lookup_buffer_inst = {};
        dbl_vec lookup_buffer_cumul = {};

        void assign_profile(const std::string &name);
        void assign_args(const str_vec &keys,
                         const str_vec &vals,
                         const str_vec &typs);
        void assign_min(bool is_discrete,
                        const std::string &min_value,
                        const std::string &min_at);
        void assign_max(const std::string &max_value,
                        const std::string &max_at);

	//===== Functional methods for getter methods =====//
    double (InfecProfile::*_func_inst_)(double, const dbl_vec &) const = &InfecProfile::_default_inst_;
    double (InfecProfile::*_func_cumul_)(double, const dbl_vec &) const = &InfecProfile::_default_cumul_;

    double _default_inst_(double t, const dbl_vec &opts) const;
    double _default_cumul_(double t, const dbl_vec &opts) const;
    double _lookup_inst_(double t, const dbl_vec &opts) const;
    double _lookup_cumul_(double t, const dbl_vec &opts) const;

    };
    InfecProfile _infec_;

    void update_infec_profile();
    void restore_infec_profile();

    //===========================//
    //===== IMPORTS/EXPORTS =====//
    //===========================//

    std::map<std::string, std::string> _report_pathtofile_ = {};

#ifdef _HOUSEPI_RCPP

    //----- Rcpp IMPORTS -----//
    void import_setup(const Rcpp::List &input);
    void import_infec_profile(const Rcpp::List &input);
    void import_parameter(const Rcpp::DataFrame &input);
    void import_augmented(const Rcpp::DataFrame &input);
    void import_dataset(const Rcpp::DataFrame &input);

    //----- Rcpp EXPORTS -----//
    Rcpp::DataFrame export_dataset();
    Rcpp::DataFrame export_data_report() const;
    Rcpp::DataFrame export_house_size_report() const;

    Rcpp::List export_setup() const;
    Rcpp::DataFrame export_parameter() const;
    Rcpp::NumericMatrix export_traces() const;
    Rcpp::NumericMatrix export_house_loglike_mtx() const;
    Rcpp::List export_infec_profile(const Rcpp::List &opts,
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

    ulong _chain_ = 1;       // current chain
    ulong _rung_ = 0;        // current thermodynamic rung
    ulong setup_n_rungs = 0; // number of rungs

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
