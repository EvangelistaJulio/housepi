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

    Impl(std::function<double(Household &)> callback_house_loglike,
         std::function<void(Household &, indiv)> callback_update_cache)
        : _infec_profile(*this),
          _uniform_distr(0., 1.),
          _callback_house_loglike(std::move(callback_house_loglike)),
          _callback_update_cache(std::move(callback_update_cache))
    {
        setup_seed = tools::default_setup_seed();
    }
    ~Impl() = default;

    //===================================//
    //===== RANDOM NUMBER GENERATOR =====//
    //===================================//

    // reference to the pseudo random number generator.
    rnd_gen &get_prng() { return (this->*_do_get_prng)(); }
    // draw a random variable uniformly distributed
    double uniform() { return _uniform_distr(get_prng()); }
    void set_prng_seed();

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
                                        "const", "discr",
                                        "prior", "arg_key", "arg_val",
                                        "min", "max",
                                        "rule", "scale"};
    std::vector<parID> param_const_list = {};
    std::vector<parID> param_standard_infer_list = {};
    std::vector<parID> param_cache_infer_list = {};

    parIDmap<Variable> param_map = {}; // parameters dictionary
    bool _is_undef_params_ = false;    // undefined parameters valid only if infer

    // get value of the parameter
    double get_param(const parID &p) const { return param_map.at(p).get(); }
    // get pointer to the parameter
    const double *get_param_pointer(const parID &p) const { return param_map.at(p).get_pointer(); }

    //----- Parameter - Setter -----//
    void init_parameters();
    void make_param_cache_dependent(const std::string &par_name);

    //----- Parameter - Metropolis-Hastings -----//
    bool propose_valid_param_move(const parID &p);
    bool accept_valid_param_move(const parID &p);

    void update_param_accept_rate(const parID &p, bool is_accepted);
    void adapt_param_rule(const parID &p, double ratio);

    void recover_parameter(const parID &p, bool cache_depen);

    //==========================//
    //===== HOUSEHOLD DATA =====//
    //==========================//

    //----- Household Data - Attributes -----//
    str_vec data_headers = {};
    size_t n_individuals = 0;
    size_t n_households = 0;
    std::vector<dataIDmap<dbl_vec>> indiv_datalist = {};
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
                                       "func", "farg_key", "farg_val", "farg_at",
                                       "value", "value_at", "init_rand",
                                       "const", "discr",
                                       "model", "arg_key", "arg_val", "arg_at",
                                       "min", "min_at",
                                       "max", "max_at",
                                       "rule", "scale", "scale_at",
                                       "from", "from_at"};
    std::vector<augmID> augm_list_full = {};
    std::vector<augmID> augm_list_model = {};
    std::vector<augmID> augm_list_infer = {};
    std::vector<augmID> augm_list_func = {};
    std::map<int, std::vector<augmID>> augm_status_map_model = {};
    std::map<int, std::vector<augmID>> augm_status_map_infer = {};

    //----- Augmented - Setter -----//
    void init_augmented_variables();
    void clear_all_augm(Individual &member);
    void assign_all_augm(Individual &member);

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
        InfecProfile(MCMC::Impl &ptr)
            : profile(std::make_unique<Distribution>()),
              pImpl(ptr) {}
        ~InfecProfile() = default;

        //===== SETTER =====//
        void build();
        void init_lookup();
        void recover_lookup();
        void update_lookup();
        void assign_template(std::map<std::string, std::string> &x);
        void append_augm_definition();

        //===== GETTER =====//
        bool is_lookup() const { return !_lookup_inst.empty(); }
        bool is_augmented() const { return !_augm_infec_profile_deflist.empty(); }
        
        double inst(Individual &member, double t) { return (this->*_do_inst)(member, t); }
        double cumul(Individual &member, double t) { return (this->*_do_cumul)(member, t); }

    private:
        MCMC::Impl &pImpl; // pointer to external class
        const str_vec _headers = {"name", "discr",
                                  "arg_key", "arg_val", "arg_at",
                                  "min", "min_at", "max", "max_at"};
        // lookup table
        size_t _size_lut = 0;
        dbl_vec _lookup_inst = {};
        dbl_vec _lookup_cumul = {};
        dbl_vec _lookup_buffer_inst = {};
        dbl_vec _lookup_buffer_cumul = {};
        dict_str_str _augm_infec_profile_deflist = {};

        void _build_generic();
        void _build_augmented();

        void _assign_single_attribute(const std::string &key,
                                      const std::string &val,
                                      const std::string &at,
                                      bool is_arg);
        void _assign_arguments(const str_vec &keys,
                               const str_vec &vals,
                               const str_vec &typs);

        //===== Functional methods for getter methods =====//
        double (InfecProfile::*_do_inst)(Individual &, double) = &InfecProfile::__default_inst;
        double (InfecProfile::*_do_cumul)(Individual &, double) = &InfecProfile::__default_cumul;

        double __default_inst(Individual &member, double t);
        double __default_cumul(Individual &member, double t);
        double __lookup_inst(Individual &member, double t);
        double __lookup_cumul(Individual &member, double t);
        double __augm_inst(Individual &member, double t);
        double __augm_cumul(Individual &member, double t);
        double __error(Individual &member, double t);
    } _infec_profile;

    void update_infec_profile();
    void recover_infec_profile();
    double infec_profile_inst(Individual &member, double t) { return _infec_profile.inst(member, t);}
    double infec_profile_cumul(Individual &member, double t) { return _infec_profile.cumul(member, t);}

    //===========================//
    //===== IMPORTS/EXPORTS =====//
    //===========================//

#ifdef _HOUSEPI_RCPP

    //----- Rcpp IMPORTS -----//
    void import_setup(const Rcpp::List &input);
    void import_dataset(const Rcpp::DataFrame &input);
    void import_datalist(const Rcpp::DataFrame &input);
    void import_multiplicity(const Rcpp::DataFrame &input);
    void import_param(const Rcpp::DataFrame &input);
    void import_infec_profile(const Rcpp::List &input);
    void import_augm(const Rcpp::DataFrame &input);

    //----- Rcpp EXPORTS -----//
    Rcpp::DataFrame export_dataset();
    Rcpp::DataFrame export_data_report() const;
    Rcpp::DataFrame export_house_size_report() const;

    Rcpp::List export_setup() const;
    Rcpp::DataFrame export_param() const;
    Rcpp::NumericMatrix export_traces() const;
    Rcpp::NumericMatrix export_house_loglike_mtx() const;
    Rcpp::List showcase_infec_profile(int n_bins) const;

#else
    //----- READ -----//
    void read_setup(const std::string &pathtofile);
    void read_dataset(const std::string &pathtofile);
    void read_datalist(const std::string &pathtofile);
    void read_multiplicity(const std::string &pathtofile);
    void read_param(const std::string &pathtofile);
    void read_augm(const std::string &pathtofile);
    void read_infec_profile(const std::string &pathtofile);

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
    void build_dataset(const dict_str_vec &input);
    void build_datalist(const dict_str_vec &input, const std::string &list_name = HOUSEPI_NA_STRING);
    void build_multiplicity(const dict_str_vec &input);
    void build_param();
    void build_augm();

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

    // Flag ot tell if power posteriors
    bool is_power_post() const { return setup_n_rungs > 0; }
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
    void update_full_loglike(const bool cache_depen);
    double update_house_loglike(size_t R, const bool is_augm_ll, const bool cache_depen);
    void run_monte_carlo(progress::Bar &progbar);
    void move_all();
    void move_all_and_adapt_parameters(ulong step);

    //----- Metropolis - parameter moves -----//
    bool move_single_parameter(const parID &p, bool cache_depen);

    //----- Metropolis - augmented moves -----//
    void move_all_augm();

private:
    //===================================//
    //===== RANDOM NUMBER GENERATOR =====//
    //===================================//

    rnd_gen _prng; // Mersenne Twister 19937 generator (32 or 64 bit)
    rnd_gen &(MCMC::Impl::*_do_get_prng)();
    rnd_gen &_prng_unique_thread() { return _prng; }

    //----- parallelization ----//
    std::vector<std::string> _error_message = {};

#ifdef _OPENMP
    std::vector<rnd_gen> _prng_parallel = {}; // parallel generators (for multi-threading)
    rnd_gen &_prng_parallel_threads() { return _prng_parallel[omp_get_thread_num()]; }
#endif

    //----- uniform sampler -----//
    std::uniform_real_distribution<double> _uniform_distr;

    //=====================//
    //===== CALLBACKS =====//
    //=====================//

    std::function<double(Household &)> _callback_house_loglike;
    std::function<void(Household &, indiv)> _callback_update_cache;

    double _compute_house_obs_loglike(Household &h) { return _callback_house_loglike(h); }
    void _update_cache_dependecies(Household &h, indiv i) { _callback_update_cache(h, i); }

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

    double _type_to_double(const Individual &member,
                           const std::string &tkn_val,
                           const std::string &tkn_typ);
    void _assign_augm_attributes(Individual &member,
                                 const augmID &a,
                                 const str_vec &keys,
                                 const str_vec &vals,
                                 const str_vec &typs,
                                 const std::vector<bool> &is_arg,
                                 bool is_var_model);

    void _assign_augm_func_args(Individual &member,
                                const augmID &a,
                                const str_vec &keys,
                                const str_vec &vals,
                                const str_vec &typs);
    void _assign_augm_func_range(Individual &member,
                                 const augmID &a,
                                 const str_vec &vals,
                                 const str_vec &typs);

    void _assign_augm_model_args(Individual &member,
                                 const augmID &a,
                                 const str_vec &keys,
                                 const str_vec &vals,
                                 const str_vec &typs);
    void _assign_augm_model_range(Individual &member,
                                  const augmID &a,
                                  const str_vec &vals,
                                  const str_vec &typs);

    //======================//
    //===== METROPOLIS =====//
    //======================//

    double _compute_augm_model_loglike(Household &house);

    //----- Metropolis : augmented moves -----//
    double move_augm_single_house(size_t H);
    double _metropolis_augm_single_move(Household &house, indiv i, const augmID &a);

    double (MCMC::Impl::*_do_move_augm_indiv)(Household &, indiv) = &MCMC::Impl::__move_augm_indiv_generic;
    double __move_augm_indiv_generic(Household &house, indiv i);
    double __move_augm_indiv_by_status(Household &house, indiv i);
};

#endif // HOUSEPI__MCMC__BUILD__H
