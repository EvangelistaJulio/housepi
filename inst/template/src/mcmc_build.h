/*      __
 *     / /_  ____  __  __ ____ __  ____  __
 *    / __ \/ __ \/ / / / ___/ _ \/ __ \/ /
 *   / / / / /_/ / /_/ (__  )  __/ /_/ / /
 *  /_/ /_/\____/\____/____/\___/ ____/_/
 *                             /_/
 *
 * Copyright (C) 2025 Francesco Camaglia, Institut Pasteur
 *
 * mcmc_build.h
 *
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
const std::string TYPE_AUGM = "A";   // Augmented type declaration

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

    Impl();
    ~Impl();

    //===================================//
    //===== RANDOM NUMBER GENERATOR =====//
    //===================================//

    rnd_gen &get_prng();
    void set_prng_seed();
    double uniform();

    //=================//
    //===== SETUP =====//
    //=================//

    //----- Setup - Attributes -----//
    const str_vec setup_file_headers = {"sample_size", "seed", "n_chains",
                                        "n_steps", "n_adapt", "thin",
                                        "n_burnin", "rand_init", "n_threads",
                                        "inv_temper", "id_trace", "progress"};

    std::int32_t setup_seed;            // random generator seed
    bool setup_rand_init_params = true; // random initialization of parameters
    ulong setup_progress = 1;           // thinning seconds for progress bar

    ulong setup_n_chains = 1; // n. chains
    ulong setup_n_steps = 1;  // n. iterations
    ulong setup_n_adapt = 0;  // n. (finite) adaptive steps
    ulong setup_thin = 1;     // n. thinning iterations
    ulong setup_n_burnin = 0; // n. burnin iterations
    ulong setup_n_rungs = 0;  // n. of rungs

    double setup_inv_T = 1.;      // inverse temperature
    unsigned setup_n_threads = 1; // number of OMP threads

    std::string setup_id_trace = HOUSEPI_NA_STRING; // augmented trace to save

    ulong setup_c_sample_size = 0; // single chain sample size
    ulong setup_c_n_sweeps = 0;    // single chain n. of sweeps

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
    void make_param_extra_dependent(const std::string &par_name);

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
    std::vector<size_t> restrict_hh_indexes = {};
    std::unordered_map<std::string, size_t> _hhid_to_index = {};

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

    private:
        MCMC::Impl &pImpl; // pointer to external class
        const str_vec _headers = {"name", "discrete",
                                  "arg_name", "arg_value", "arg_at",
                                  "min", "min_at",
                                  "max", "max_at"};
        // lookup table
        size_t _size_lut = 0;
        dbl_vec _lookup_inst = {};
        dbl_vec _lookup_cumul = {};
        dbl_vec _lookup_buffer_inst = {};
        dbl_vec _lookup_buffer_cumul = {};

        void _assign_profile(const std::string &name);
        void _assign_args(const str_vec &keys,
                          const str_vec &vals,
                          const str_vec &typs);
        void _assign_range(const std::string &range_name,
                           const std::string &range_at,
                           const std::string &range_value);

        //===== Functional methods for getter methods =====//
        double (InfecProfile::*_do_inst)(double, const dbl_vec &) const = &InfecProfile::__default_inst;
        double (InfecProfile::*_do_cumul)(double, const dbl_vec &) const = &InfecProfile::__default_cumul;

        double __default_inst(double t, const dbl_vec &opts) const;
        double __default_cumul(double t, const dbl_vec &opts) const;
        double __lookup_inst(double t, const dbl_vec &opts) const;
        double __lookup_cumul(double t, const dbl_vec &opts) const;
    };
    InfecProfile _infec_profile;

    void update_infec_profile();
    void restore_infec_profile();

    //===========================//
    //===== IMPORTS/EXPORTS =====//
    //===========================//

#ifdef _HOUSEPI_RCPP

    //----- Rcpp IMPORTS -----//
    void import_setup(const Rcpp::List &input);
    void import_infec_profile(const Rcpp::List &input);
    void import_param(const Rcpp::DataFrame &input);
    void import_augm(const Rcpp::DataFrame &input);
    void import_dataset(const Rcpp::DataFrame &input);
    void import_multiplicity(const Rcpp::DataFrame &input);

    //----- Rcpp EXPORTS -----//
    Rcpp::DataFrame export_dataset();
    Rcpp::DataFrame export_data_report() const;
    Rcpp::DataFrame export_house_size_report() const;

    Rcpp::List export_setup() const;
    Rcpp::DataFrame export_param() const;
    Rcpp::NumericMatrix export_traces() const;
    Rcpp::NumericMatrix export_house_loglike_mtx() const;
    Rcpp::List showcase_infec_profile(const Rcpp::List &opts,
                                      int n_bins) const;

#else
    //----- READ -----//
    void read_dataset(const std::string &pathtofile);
    void read_setup(const std::string &pathtofile);
    void read_param(const std::string &pathtofile);
    void read_infec_profile(const std::string &pathtofile);
    void read_augm(const std::string &pathtofile);
    void read_multiplicity(const std::string &pathtofile);

    //----- SHOW -----//
    void show_summary() const;
    void show_setup() const;
    void show_all_params() const;
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
    void build_param();
    void build_augm();
    void build_dataset(const dict_str_vec &input);
    void build_multiplicity(const dict_str_vec &input);

    //=========================//
    //===== LOGLIKELIHOOD =====//
    //=========================//

    double full_loglike = 0.; // current full loglikelihood
    double full_buffer = 0.;  // buffer loglikelihood

    double get_full_loglike() const;

    void set_full_loglike(double new_ll);
    void increment_full_loglike(double delta);
    const std::pair<double, double> smpl_full_loglike_stats() const;

    //===============================//
    //===== METROPOLIS-HASTINGS =====//
    //===============================//

    ulong chain = 1; // current chain
    ulong rung = 0;  // current thermodynamic rung

    bool is_power_post() const;
    const std::string fmt_progtxt();

    void assert_proper_dataset();
    void assert_proper_mcmc_setup();

    //----- Containers for posterior sample values -----//
    ulong smpl_count = 0;
    dbl_vec smpl_full_loglike = {};
    size_t smpl_H_star = 0; // index of house to follow up
    indiv smpl_i_star = 0;  // index of individual to follow up
    parIDmap<dbl_vec> smpl_param_value = {};
    augmIDmap<dbl_vec> smpl_augm_value = {};

    void init_move_counts();
    void init_diagnostics();
    void update_diagnostics();

    //----- Metropolis - augmented moves -----//
    double compute_augm_model_loglike(const Household &house) const;

private:
    //===================================//
    //===== RANDOM NUMBER GENERATOR =====//
    //===================================//

    rnd_gen _prng; // Mersenne Twister 19937 generator (32 or 64 bit)
    rnd_gen &(MCMC::Impl::*_do_get_prng)();
    rnd_gen &_prng_unique_thread() { return _prng; }

#ifdef _OPENMP
    std::vector<rnd_gen> _prng_parallel = {}; // parallel generators (for multi-threading)
    rnd_gen &_prng_parallel_threads() { return _prng_parallel[omp_get_thread_num()]; }
#endif

    //----- uniform sampler -----//
    std::uniform_real_distribution<double> _uniform_distr;

    //=================//
    //===== BUILD =====//
    //=================//

    dict_str_str _setup_deflist = {};
    dict_str_vec _param_deftable = {};
    dict_str_vec _augm_deftable = {};

    //----- build augmented : auxiliary methods -----//
    const dict_str_vec _build_augm_name_status_check();
    void _build_augm_allocate_univocal_defs(const dict_str_vec &check);
    void _build_augm_clear_houses_and_set_id_trace();

    double _type_to_double(const std::string &tkn_val,
                           const std::string &tkn_typ,
                           const Household &house, indiv i);
    void _assign_augm_model(Household &house, indiv i, augmID a,
                            const str_vec &keys,
                            const str_vec &vals,
                            const str_vec &typs,
                            const std::vector<bool> &is_arg);

    //======================//
    //===== METROPOLIS =====//
    //======================//

    void (Impl::*_do_augm_model_by_definition)(const Household &, indiv, double &) const = &Impl::__augm_model_generic;
    void __augm_model_by_status(const Household &house, indiv i, double &delta) const;
    void __augm_model_generic(const Household &house, indiv i, double &delta) const;
};

#endif // HOUSEPI__MCMC__BUILD__H
