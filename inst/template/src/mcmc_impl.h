/*      __
 *     / /_  ____  __  __ ____ __  ____  __
 *    / __ \/ __ \/ / / / ___/ _ \/ __ \/ /
 *   / / / / /_/ / /_/ (__  )  __/ /_/ / /
 *  /_/ /_/\____/\____/____/\___/ ____/_/
 *                             /_/
 *
 * Copyright (C) 2025 Francesco Camaglia, Institut Pasteur
 *
 * mcmc_impl.h
 *
 */
#ifndef HOUSEPI__MCMC__IMPL__H
#define HOUSEPI__MCMC__IMPL__H

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

    Impl(std::function<void()> func_init_inference,
         std::function<void()> func_init_simulation,
         std::function<double(Household &)> func_household_loglike,
         std::function<void(Household &)> func_simulate_house,
         std::function<void(Household &, Individual &)> func_update_cache)
        : setup(*this), param(*this), dataset(*this), datalist(*this),
          group_rates(*this), augm(*this), infec_profile(*this),
          callback_init_inference_(std::move(func_init_inference)),
          callback_init_simulation_(std::move(func_init_simulation)),
          callback_house_loglike_(std::move(func_household_loglike)),
          callback_simulate_house_(std::move(func_simulate_house)),
          callback_update_cache_(std::move(func_update_cache)) {}
    ~Impl() = default;

    //=================//
    //===== SETUP =====//
    //=================//

    /* Implementation of Setup. */
    struct ImplSetup
    {
        ImplSetup(MCMC::Impl &ptr)
            : pImpl_(ptr), unif_(0., 1.) { seed = tools::default_setup_seed(); }
        ~ImplSetup() = default;

        const str_vec file_headers = {"sample_size", "seed", "n_chains",
                                      "n_steps", "n_adapt", "thin",
                                      "n_burnin", "rand_init", "n_threads",
                                      "id_trace", "progress"};

        std::int32_t seed;            // random generator seed
        bool rand_init_params = true; // random initialization of parameters
        ulong progress = 1;           // thinning seconds for progress bar
        ulong n_chains = 1;           // n. chains
        ulong n_steps = 1;            // n. iterations
        ulong n_adapt = 0;            // n. (finite) adaptive steps
        ulong thin = 1;               // n. thinning iterations
        ulong n_burnin = 0;           // n. burnin iterations
        ulong n_rungs = 0;            // n. of rungs
        ulong c_sample_size = 0;      // single chain sample size
        ulong c_n_sweeps = 0;         // single chain n. of sweeps
        double inv_T = 1.;            // inverse temperature
        unsigned n_threads = 1;       // number of OMP threads
        str_vec id_trace = {};        // augmented trace to save

        // state variables for MCMC sampling
        ulong chain = 1; // current chain
        ulong rung = 0;  // current thermodynamic rung

        bool is_power_post() const { return n_rungs > 0; }
        const std::string fmt_progtxt();

        void build(const dict_str_str &deflist);
        void assert_mcmc_setup();
        void set_n_threads();

        //----- Progress bar -----//
        void prog_set_text(const std::string &text = HSP_NA_STRING)
        {
            progressbar_.set_text((text == HSP_NA_STRING) ? fmt_progtxt() : text);
        }
        void prog_declare(const std::string &text = HSP_NA_STRING)
        {
            prog_set_text(text);
            progressbar_.set_thintime(progress);
        }
        void prog_open(ulong total) { progressbar_.open(total); }
        void prog_update() { progressbar_.update(); }
        void prog_close() { progressbar_.close(); }
        void prog_crash() { progressbar_.crash(); }

        //----- Pseudo Random Number Geenerator -----//
        void set_prng_seed();
        rnd_gen &prng() { return (this->*get_prng_)(); }
        double unif() { return unif_(prng()); }

#ifdef _HOUSEPI_RCPP
        void rcpp_import(SEXP input);
        Rcpp::List rcpp_export() const;

        Rcpp::IntegerVector traces_chain_tags() const;
        Rcpp::CharacterVector traces_iter_names() const;
#else
        void read(const std::string &pathtofile);
        std::string show() const;
#endif

    private:
        const std::string section_name_ = "setup";
        MCMC::Impl &pImpl_; // pointer to external class

        progress::Bar progressbar_; // progress bar

        rnd_gen prng_; // Mersenne Twister 19937 generator (32 or 64 bit)
        rnd_gen &(MCMC::Impl::ImplSetup::*get_prng_)();
        rnd_gen &do_prng_unique_thread_() { return prng_; }
        std::uniform_real_distribution<double> unif_; // uniform sampler

#ifdef _OPENMP
        std::vector<rnd_gen> prng_parallel_ = {}; // parallel generators (for multi-threading)
        rnd_gen &do_prng_parallel_threads_() { return prng_parallel_[omp_get_thread_num()]; }
#endif
    } setup;

    std::string error_message = {}; // error message for the last exception

    //=====================//
    //===== PARAMETER =====//
    //=====================//

    /* Implementation of Parameters. */
    struct ImplParameter
    {
        const str_vec file_headers = {"name", "value",
                                      "const", "discr",
                                      "prior", "arg_key", "arg_val",
                                      "min", "max",
                                      "rule", "scale"};
        std::vector<par_key> list_const = {};
        std::vector<par_key> list_sample_standard = {};
        std::vector<par_key> list_sample_cache = {};
        std::vector<par_key> list_sample = {};

        par_key_map<Variable> dictionary = {}; // parameters dictionary
        par_key_map<dbl_vec> smpl_values = {}; // values of the sample

        //===== CONSTRUCTOR =====//
        ImplParameter(MCMC::Impl &ptr)
            : pImpl_(ptr) {}
        ~ImplParameter() = default;

        //===== GETTER =====//
        // get value of the parameter
        double get(const par_key &p) const { return dictionary.at(p).get(); }
        // get pointer to the parameter
        const double *get_ptr(const par_key &p) const { return dictionary.at(p).get_ptr(); }

        // flag to check if there are parameters to sample
        bool is_sample() const { return !list_sample.empty(); }

        //===== SETTER =====//
        void build(const dict_str_vec &input);
        void make_cache_dependent(const par_key &p);
        void init();

        void init_traces(ulong sample_size);
        void update_traces(ulong smpl_count);

#ifdef _HOUSEPI_RCPP
        void rcpp_import(SEXP input);
        Rcpp::NumericVector rcpp_snapshot() const;
        Rcpp::DataFrame rcpp_export() const;
#else
        void read(const std::string &pathtofile);
        void show() const;
        void show_results() const;
#endif
    private:
        const std::string section_name_ = "parameter";
        MCMC::Impl &pImpl_; // pointer to external class
    } param;

    //=====================//
    //=====  DATASET  =====//
    //=====================//

    /* Implementation of Dataset. */
    struct ImplDataSet
    {
        //===== CONSTRUCTOR =====//
        ImplDataSet(MCMC::Impl &ptr)
            : pImpl_(ptr) {}
        ~ImplDataSet() = default;

        Household &operator()(size_t H) { return house_vector_[H]; }
        const Household &operator()(const std::string &house_id) const
        {
            if (!hhid_to_index_.count(house_id))
                throw std::invalid_argument(
                    "Index \"" + house_id + "\" does not match any household in the dataset.");
            return house_vector_.at(hhid_to_index_.at(house_id));
        }

        //===== GETTER =====//
        size_t n_houses() const { return house_vector_.size(); }
        size_t n_restricted_houses() const { return restrict_indexes_.size(); }
        size_t n_individuals() const { return n_individuals_; }
        void fill_report_data();
        ulong get_report_indiv(const std::string &key) const;
        ulong get_report_house(const std::string &key) const;
        std::pair<ushort, ulong> get_report_hh_size(const std::string &key) const;
        bool is_empty() const { return house_vector_.empty(); }
        const str_vec &headers() const { return headers_; }
        Household copy(size_t H) const { return house_vector_[H]; }

        //===== SETTER =====//
        void build(const dict_str_vec &input);
        void build_multiplicity(const dict_str_vec &input);
        void assign_cache_integer_vector(const cache_key &e);
        void assign_cache_double_vector(const cache_key &e);
        void assign_cache_double_matrix(const cache_key &e);
        void init_traces(ulong sample_size);
        void update_traces(ulong smpl_count);

        //===== ITERATORS =====//
        auto begin() { return house_vector_.begin(); }
        auto begin() const { return house_vector_.begin(); }
        auto end() { return house_vector_.end(); }
        auto end() const { return house_vector_.end(); }

        // T can be Individual or const Individual
        template <typename T>
        class ViewIterator
        {
            // conditional to select the vector type if constant or not
            using VectorType = std::conditional_t<std::is_const<T>::value,
                                                  const std::vector<Household>,
                                                  std::vector<Household>>;

            VectorType &pVector;
            const std::vector<size_t> &indexes;
            size_t pos;

        public:
            ViewIterator(VectorType &m, const std::vector<size_t> &i, size_t p)
                : pVector(m), indexes(i), pos(p) {}

            T &operator*() const { return pVector[indexes[pos]]; }
            T *operator->() const { return &pVector[indexes[pos]]; }

            ViewIterator &operator++()
            {
                ++pos;
                return *this;
            }
            ViewIterator operator++(int)
            {
                ViewIterator tmp = *this;
                ++(*this);
                return tmp;
            }

            bool operator==(const ViewIterator &other) const { return pos == other.pos; }
            bool operator!=(const ViewIterator &other) const { return !(*this == other); }
        };

        template <typename T>
        class View
        {
            // conditional to select the vector type if constant or not
            using VectorType = std::conditional_t<std::is_const<T>::value, const std::vector<Household>, std::vector<Household>>;
            using iterator = ViewIterator<T>;

        public:
            VectorType &vec;
            const std::vector<size_t> &indexes;

            View(VectorType &m, const std::vector<size_t> &i)
                : vec(m), indexes(i) {}

            iterator begin() { return iterator(vec, indexes, 0); }
            iterator end() { return iterator(vec, indexes, indexes.size()); }
        };

        View<Household> restricted() { return View<Household>(house_vector_, restrict_indexes_); }
        View<const Household> restricted() const { return View<const Household>(house_vector_, restrict_indexes_); }
        Household &restricted(size_t R) { return house_vector_[restrict_indexes_[R]]; }

        //===== IMPORT/EXPORT =====//
#ifdef _HOUSEPI_RCPP
        void rcpp_import(SEXP input);
        void rcpp_import_multiplicity(SEXP input);

        Rcpp::DataFrame rcpp_snapshot(const std::string &house_id) const;
        Rcpp::DataFrame rcpp_export() const;
        Rcpp::DataFrame rcpp_export_report() const;
        Rcpp::DataFrame rcpp_export_size_report() const;
        Rcpp::NumericMatrix rcpp_export_traces() const;
#else
        void read(const std::string &pathtofile);
        void read_multiplicity(const std::string &pathtofile);
        void show_summary() const;

        void print_open(const std::string &pathtofile);
        void print_update(Household &house);
        void print_close();
#endif
    private:
        const std::string section_name_ = "dataset";
        MCMC::Impl &pImpl_; // pointer to external class

        str_vec headers_ = {};

        size_t n_individuals_ = 0;
        std::vector<Household> house_vector_ = {};
        std::unordered_map<std::string, size_t> hhid_to_index_ = {};
        std::vector<size_t> restrict_indexes_ = {};

        std::vector<dbl_vec> smpl_values_ = {};

        std::map<std::string, ulong> report_indiv_ = {};
        std::map<std::string, ulong> report_house_ = {};
        std::map<std::string, std::pair<ushort, ulong>> report_size_ = {};

        std::ofstream aux_dataset_ofstream_;
    } dataset;

    //======================//
    //=====  DATALIST  =====//
    //======================//

    /* Implementation of Individual Datalist. */
    struct ImplDataList
    {
        //===== CONSTRUCTOR =====//
        ImplDataList(MCMC::Impl &ptr)
            : pImpl_(ptr) {}
        ~ImplDataList() = default;

        //===== GETTER =====//
        bool is_empty() const { return indiv_datalist_.empty(); }
        const str_vec &headers() const { return headers_; }

        //===== SETTER =====//
        void build(const dict_str_vec &input);

        //===== IMPORT/EXPORT =====//
#ifdef _HOUSEPI_RCPP
        void rcpp_import(SEXP input);

        Rcpp::DataFrame rcpp_snapshot(const std::string &house_id) const;
        Rcpp::DataFrame rcpp_export() const;
#else
        void read(const std::string &pathtofile);
#endif
    private:
        const std::string section_name_ = "datalist";
        MCMC::Impl &pImpl_; // pointer to external class
        str_vec headers_ = {};
        std::vector<data_key_map<dbl_vec>> indiv_datalist_ = {};
    } datalist;

    //=======================//
    //===== GROUP RATES =====//
    //=======================//

    /* Implementation of Group Rates. */
    struct ImplGroupRates
    {
        //===== CONSTRUCTOR =====//
        ImplGroupRates(MCMC::Impl &ptr)
            : pImpl_(ptr) {}
        ~ImplGroupRates() = default;

        //===== GETTER =====//
        bool is_empty() const { return headers_.empty(); }

        //===== SETTER =====//
        void build(const dict_str_vec &input);

        //===== IMPORT/EXPORT =====//
#ifdef _HOUSEPI_RCPP
        void rcpp_import(SEXP input);
        Rcpp::DataFrame rcpp_snapshot(const std::string &house_id) const;
#else
        void read(const std::string &pathtofile);
#endif
    private:
        const std::string section_name_ = "group rates";
        MCMC::Impl &pImpl_; // pointer to external class

        const str_vec file_headers_ = {"group_rate", "param_name"};
        str_vec headers_ = {};

        void assign_rate_as_ptrs_(Individual &member,
                                  const rate_key_mapvec<par_key> &these_rates);
    } group_rates;

    //=====================//
    //===== AUGMENTED =====//
    //=====================//

    /* Implementation of Augmented Variables. */
    struct ImplAugmented
    {
        dict_str_vec deftable = {};
        const str_vec file_headers = {"name", "status",
                                      "func", "farg_key", "farg_val", "farg_at",
                                      "value", "value_at", "init_rand",
                                      "const", "discr",
                                      "model", "arg_key", "arg_val", "arg_at",
                                      "min", "min_at",
                                      "max", "max_at",
                                      "rule", "scale", "scale_at",
                                      "from", "from_at"};

        std::vector<augm_key> list_full = {};
        std::vector<augm_key> list_model = {};
        std::vector<augm_key> list_sample = {};
        std::vector<augm_key> list_func = {};
        std::map<int, std::vector<augm_key>> status_map_model = {};
        std::map<int, std::vector<augm_key>> status_map_infer = {};

        // optional augmented posterior sample
        std::vector<const Individual *> smpl_individuals = {};
        std::vector<augm_key_map<dbl_vec>> smpl_values = {};

        //===== CONSTRUCTOR =====//
        ImplAugmented(MCMC::Impl &ptr)
            : pImpl_(ptr) {}
        ~ImplAugmented() = default;

        //===== SETTER =====//
        void build();
        void init();

        void clear_all(Individual &member);
        void assign_all(Individual &member);

        void init_traces(ulong sample_size);
        void update_traces(ulong smpl_count);

        //===== GETTER =====//
        const std::vector<augm_key> list_export() const
        {
            return list_sample + list_func - str_vec({AUX_IP_INST_PROFILE, AUX_IP_CUMUL_PROFILE});
        }

        /* Flag if there are ant augmented variable definitions. */
        bool is_any() const { return !list_full.empty(); }
        /* Flag if there are non-constant augmented variable definitions. */
        bool is_sample() const { return !list_sample.empty(); }
        /* Flag if there are augmented variable definitions with a model. */
        bool is_model() const { return !list_model.empty(); }
        /* Flag if augmented variable definition are generic or not (i. e. individual status specific). */
        bool is_generic() const { return status_map_infer.empty(); }
        /* Flag if augmented traces are stored. */
        bool is_traces() const { return !smpl_individuals.empty(); }

#ifdef _HOUSEPI_RCPP
        void rcpp_import(SEXP input);
        Rcpp::DataFrame rcpp_snapshot(const std::string &house_id) const;
#else
        void read(const std::string &pathtofile);
#endif

    private:
        const std::string section_name_ = "augmented variables";
        MCMC::Impl &pImpl_; // pointer to external class

        const dict_str_vec check_name_status_(const dict_str_vec &deftable);
        void allocate_univocal_defs_(const dict_str_vec &check);
        void set_id_trace_(str_vec &id_trace);

        double type_to_double_(const Individual &member,
                               const std::string &tkn_val, const std::string &tkn_typ);
        void assign_attributes_(Individual &member, const augm_key &a,
                                const str_vec &keys, const str_vec &vals, const str_vec &typs,
                                const std::vector<bool> &is_arg, bool is_var_model);

        void assign_func_args_(Individual &member, const augm_key &a,
                               const str_vec &keys, const str_vec &vals, const str_vec &typs);
        void assign_func_range_(Individual &member, const augm_key &a,
                                const str_vec &vals, const str_vec &typs);

        void assign_model_args_(Individual &member, const augm_key &a,
                                const str_vec &keys, const str_vec &vals, const str_vec &typs);
        void assign_model_range_(Individual &member, const augm_key &a,
                                 const str_vec &vals, const str_vec &typs);
    } augm;

    //===============================//
    //===== INFECTIVITY PROFILE =====//
    //===============================//

    /* Implementation of Infectivity Profile. */
    struct ImplInfecProfile
    {
    public:
        //===== CONSTRUCTOR =====//
        ImplInfecProfile(MCMC::Impl &ptr)
            : pImpl_(ptr),
              profile_(std::make_unique<Distribution>()) {}
        ~ImplInfecProfile() = default;

        //===== SETTER =====//
        void build(const dict_str_str &deflist);
        void assign_template(const dict_str_str &deflist);
        void append_augm_definition();
        void init();
        void update(bool is_cache_depen);
        void recover(bool is_cache_depen);

        //===== GETTER =====//
        bool is_lookup() const { return !lookup_inst_.empty(); }
        bool is_augmented() const { return !augm_infec_profile_deflist_.empty(); }

        double inst(const Individual &member, double t) { return (this->*inst_)(member, t); }
        double cumul(const Individual &member, double t) { return (this->*cumul_)(member, t); }

#ifdef _HOUSEPI_RCPP
        void rcpp_import(SEXP input);

        Rcpp::List rcpp_snapshot(const std::string &house_id) const;
        Rcpp::List rcpp_showcase(int n_bins) const;
#else
        void read(const std::string &pathtofile);
#endif

    private:
        const std::string section_name_ = "infectivity profile";
        MCMC::Impl &pImpl_; // pointer to external class
        std::shared_ptr<Distribution> profile_ = nullptr;
        dict_str_str deflist_template_ = {};

        const str_vec headers_ = {"name", "discr",
                                  "arg_key", "arg_val", "arg_at",
                                  "max", "max_at"};

        // lookup table
        size_t size_lut_ = 0;
        dbl_vec lookup_inst_ = {};
        dbl_vec lookup_cumul_ = {};
        dbl_vec lookup_buffer_inst_ = {};
        dbl_vec lookup_buffer_cumul_ = {};
        // augmented function
        dict_str_vec augm_infec_profile_deflist_ = {};

        void init_lookup_();

        void build_generic_profile_(const dict_str_str &deflist);
        void build_augmented_profile_(const dict_str_str &deflist);

        void assign_single_attribute_(const std::string &key, const std::string &val, const std::string &at,
                                      bool is_arg);
        void assign_arguments_(const str_vec &keys, const str_vec &vals, const str_vec &typs);

        //===== Functional methods for getter methods =====//
        double (ImplInfecProfile::*inst_)(const Individual &, double) = &ImplInfecProfile::do_inst_default_;
        double do_inst_default_(const Individual &member, double t);
        double do_inst_lookup_(const Individual &member, double t);
        double do_inst_augm_(const Individual &member, double t);

        double (ImplInfecProfile::*cumul_)(const Individual &, double) = &ImplInfecProfile::do_cumul_default_;
        double do_cumul_default_(const Individual &member, double t);
        double do_cumul_lookup_(const Individual &member, double t);
        double do_cumul_augm_(const Individual &member, double t);
    } infec_profile;

    //=================================//
    //=====  METROPOLIS-HASTINGS  =====//
    //=================================//

    double evaluate();
    void infer();

#ifdef _HOUSEPI_RCPP
    //----- C++ Print -----//
    Rcpp::DataFrame simulate();
    Rcpp::List evidence(const std::string &scheme, ulong n_rungs);

    //----- Rcpp EXPORTS -----//
    Rcpp::NumericMatrix rcpp_export_traces() const;
    Rcpp::List rcpp_export_house_snapshot(const std::string &house_id) const;
#else
    //----- C++ Print -----//
    void print_infer_state(const std::string &output_filename = HSP_NA_STRING) const;
    void print_simulated_dataset(const std::string &output_filename = HSP_NA_STRING);
    void print_evidence(const std::string &scheme, ulong n_rungs);
#endif

private:
    //=======================================//
    //===== CALLBACKS (VIRTUAL METHODS) =====//
    //=======================================//

    std::function<void()> callback_init_inference_;
    std::function<void()> callback_init_simulation_;
    std::function<double(Household &)> callback_house_loglike_;
    std::function<void(Household &)> callback_simulate_house_;
    std::function<void(Household &, Individual &)> callback_update_cache_;

    void init_model_for_inferencee_() { callback_init_inference_(); }
    void init_model_for_simulation_() { callback_init_simulation_(); }
    void simulate_single_household_(Household &h) { callback_simulate_house_(h); }
    void update_cache_dependecies_(Household &h, Individual &i) { callback_update_cache_(h, i); }

    double compute_and_set_household_obs_loglike_(Household &h)
    {
        const double output = callback_house_loglike_(h);
        h.set_obs_loglike(output);
        return output;
    }

    //=================================//
    //=====  METROPOLIS-HASTINGS  =====//
    //=================================//

    double full_loglike_ = 0.; // current full loglikelihood
    double full_buffer_ = 0.;  // buffer loglikelihood

    ulong smpl_count_ = 0;
    dbl_vec smpl_full_loglike_ = {};

    void assert_proper_dataset_();
    void init_monte_carlo_();

    /* Set new full loglikelihood value and buffer old value. */
    void set_full_loglike_(double new_ll)
    {
        full_buffer_ = full_loglike_;
        full_loglike_ = new_ll;
    }

    /* Increment full loglikelohood by a delta, buffer old value. */
    void increment_full_loglike_(double delta) { set_full_loglike_(full_loglike_ + delta); }

    void init_move_counts_();
    void init_infer_diagnostics_();
    void update_posterior_sample_();
    void update_full_loglike_(bool cache_depen);
    double update_house_loglike_(Household &house, bool is_augm_ll, bool is_cache_depen);
    void run_monte_carlo_();
    void move_all_();
    void move_all_and_adapt_params_(ulong step);

    //----- Metropolis - parameter moves -----//
    bool move_single_parameter_(const par_key &p, bool cache_depen);
    bool propose_valid_param_move_(const par_key &p);
    bool accept_valid_param_move_(const par_key &p);
    void recover_parameter_(const par_key &p, bool cache_depen);

    //----- Metropolis : augmented moves -----//
    void move_all_augm_();
    double move_augm_single_house_(Household &house);
    double metropolis_augm_single_move_(Household &house, Individual &member, const augm_key &a);

    double (MCMC::Impl::*move_augm_indiv_)(Household &, Individual &) = &MCMC::Impl::do_move_augm_indiv_generic_;
    double do_move_augm_indiv_generic_(Household &house, Individual &member);
    double do_move_augm_indiv_by_status_(Household &house, Individual &member);

    /* Compute statistics from sample of full loglikelihood. */
    const std::pair<double, double> smpl_full_loglike_stats_() const { return tools::compute_stats(smpl_full_loglike_); }

    //----- Metropolis : power posteriors -----//
    const std::vector<dbl_vec> power_posteriors_(const std::string &scheme, ulong N);

    //====================//
    //===== SIMULATE =====//
    //====================//

    void aux_init_simulate_();
};

#endif // HOUSEPI__MCMC__IMPL__H
