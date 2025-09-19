/*
 *        __
 *       / /_  ____  __  __ ____ __  ____  __
 *      / __ \/ __ \/ / / / ___/ _ \/ __ \/ /
 *     / / / / /_/ / /_/ (__  )  __/ /_/ / /
 *    /_/ /_/\____/\____/____/\___/ ____/_/
 *                               /_/
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
#include "progress.h"
#include "distribution.h"
#include "variable.h"
#include "household.h"
#include "mcmc.h"

/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////  MCMC : Impl : class  ////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

/// Monte Carlo Markov Chain implementation class.
class MCMC::Impl
{
public:
    /////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////  MCMC : Impl : constructor  ///////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////

    Impl(std::function<void()> func_init_inference,
         std::function<void()> func_init_simulation,
         std::function<double(Household &)> func_household_loglike,
         std::function<void(Household &)> func_simulate_house,
         std::function<void(Household &, Individual &)> func_update_cache)
        : setup(*this), param(*this), dataset(*this), datalong(*this),
          group_rates(*this), augm(*this), infec_profile(*this),
          callback_init_inference_(std::move(func_init_inference)),
          callback_init_simulation_(std::move(func_init_simulation)),
          callback_house_loglike_(std::move(func_household_loglike)),
          callback_simulate_house_(std::move(func_simulate_house)),
          callback_update_cache_(std::move(func_update_cache)) {}
    ~Impl() = default;

    /////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////  MCMC : Impl : Setup : struct  /////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////

    /// Implementation of Setup.
    struct Setup
    {
        Setup(MCMC::Impl &ptr)
            : pImpl_(ptr), unif_(0., 1.) {}
        ~Setup() = default;

        const str_vec file_headers = {"sample_size", "seed", "n_chains",
                                      "n_steps", "n_adapt", "thin",
                                      "n_burnin", "rand_init", "n_threads",
                                      "id_draws", "progress"};

        std::int32_t seed;               /// random generator seed
        bool rand_init_params = true;    /// random initialization of parameters
        unsigned long progress = 1;      /// thinning seconds for progress bar
        unsigned long n_chains = 1;      /// n. chains
        unsigned long n_steps = 1;       /// n. iterations
        unsigned long n_adapt = 0;       /// n. (finite) adaptive steps
        unsigned long thin = 1;          /// n. thinning iterations
        unsigned long n_burnin = 0;      /// n. burnin iterations
        unsigned long n_rungs = 0;       /// n. of rungs
        unsigned long c_sample_size = 0; /// single chain sample size
        unsigned long c_n_sweeps = 0;    /// single chain n. of sweeps
        double inv_T = 1.;               /// inverse temperature
        unsigned n_threads = 1;          /// n. of OMP threads
        str_vec id_draws = {};           /// augmented trace to save

        unsigned long chain = 1; /// current MCMC chain
        unsigned long rung = 0;  /// current MCMC thermodynamic rung

        bool is_power_post() const { return n_rungs > 0; }
        std::string fmt_progtxt();

        void build(const dict_str &deflist);
        void assert_mcmc_setup();
        void set_n_threads();

        // progress bar
        void prog_set_text(const std::string &text = HSP_NA_STRING)
        {
            progressbar_.set_text((text == HSP_NA_STRING) ? fmt_progtxt() : text);
        }
        void prog_declare(const std::string &text = HSP_NA_STRING)
        {
            prog_set_text(text);
            progressbar_.set_thintime(progress);
        }
        void prog_open(unsigned long total) { progressbar_.open(total); }
        void prog_update() { progressbar_.update(); }
        void prog_close() { progressbar_.close(); }
        void prog_crash() { progressbar_.crash(); }

        // Pseudo Random Number Geenerator
        void set_prng_seed();
        rnd_gen &prng() { return (this->*get_prng_)(); }
        double unif() { return unif_(prng()); }

#ifdef _HOUSEPI_RCPP
        void rcpp_import(SEXP input);
        Rcpp::List rcpp_export() const;

        Rcpp::IntegerVector draws_chain_tags() const;
        Rcpp::CharacterVector draws_iter_names() const;
#else
        void read(const std::string &pathtofile);
        void show() const;
#endif

    private:
        const std::string section_name_ = "setup";
        MCMC::Impl &pImpl_; // pointer to external class

        progress::Bar progressbar_; // progress bar

        rnd_gen prng_; // Mersenne Twister 19937 generator (32 or 64 bit)
        rnd_gen &(MCMC::Impl::Setup::*get_prng_)();
        rnd_gen &do_prng_unique_thread_() { return prng_; }
        std::uniform_real_distribution<double> unif_; // uniform sampler

#ifdef _OPENMP
        std::vector<rnd_gen> prng_parallel_ = {}; // parallel generators (for multi-threading)
        rnd_gen &do_prng_parallel_threads_() { return prng_parallel_[omp_get_thread_num()]; }
#endif
    } setup;

    std::string error_message = {}; // error message for the last exception

    /////////////////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////  MCMC : Impl : Parameter : struct  ///////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////

    /// Implementation of Parameters.
    struct Parameter
    {
        const str_vec file_headers = {"name", "value",
                                      "const", "discr",
                                      "prior", "arg_key", "arg_val",
                                      "min", "max",
                                      "rule", "scale"};
        str_vec list_const = {};
        str_vec list_sample_standard = {};
        str_vec list_sample_cache = {};
        str_vec list_sample = {};

        // values of the sample
        u_dict<dbl_vec> smpl_values = {};

        /////////////////////////////////////////////////////////////////////////////////////////
        /////////////////////////  MCMC : Impl : Parameter : constructor  ///////////////////////
        /////////////////////////////////////////////////////////////////////////////////////////

        Parameter(MCMC::Impl &ptr)
            : pImpl_(ptr) {}
        ~Parameter() = default;

        /////////////////////////////////////////////////////////////////////////////////////////
        ///////////////////////////  MCMC : Impl : Parameter : getter  //////////////////////////
        /////////////////////////////////////////////////////////////////////////////////////////

        Variable &at(const std::string &p) { return safe_at(dictionary_, "parameter", p); }
        const Variable &at(const std::string &p) const { return safe_at(dictionary_, "parameter", p); }
        // get value of the parameter
        double get(const std::string &p) const { return this->at(p).get(); }
        // get pointer to the parameter
        observ_dbl_ptr get_ptr(const std::string &p) const { return this->at(p).get_ptr(); }

        // flag to check if there are parameters to sample
        bool is_any_sample() const { return !list_sample.empty(); }

        /////////////////////////////////////////////////////////////////////////////////////////
        ///////////////////////////  MCMC : Impl : Parameter : setter  //////////////////////////
        /////////////////////////////////////////////////////////////////////////////////////////

        void build(const dict_str_vec &input);
        void make_cache_dependent(const std::string &p);
        void init();

        void init_draws(unsigned long sample_size);
        void update_draws(unsigned long smpl_count);
        void clear_draws();

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

        u_dict<Variable> dictionary_ = {}; // parameters dictionary
    } param;

    /////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////  MCMC : Impl : DataWide : struct  ////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////

    /// Implementation of Dataset.
    struct DataWide
    {
        /////////////////////////////////////////////////////////////////////////////////////////
        /////////////////////////  MCMC : Impl : DataSet : constructor  /////////////////////////
        /////////////////////////////////////////////////////////////////////////////////////////

        DataWide(MCMC::Impl &ptr)
            : pImpl_(ptr) {}
        ~DataWide() = default;

        Household &operator()(size_t H) { return house_vector_[H]; }
        const Household &operator()(const std::string &house_id) const
        {
            if (!hhid_to_index_.count(house_id))
                throw std::invalid_argument(
                    "Index \"" + house_id + "\" does not match any household in the dataset.");
            return house_vector_.at(hhid_to_index_.at(house_id));
        }

        /////////////////////////////////////////////////////////////////////////////////////////
        ////////////////////////////  MCMC : Impl : DataWide : getter  ///////////////////////////
        /////////////////////////////////////////////////////////////////////////////////////////

        size_t n_houses() const { return house_vector_.size(); }
        size_t n_restricted_houses() const { return restrict_indexes_.size(); }
        size_t n_individuals() const { return n_individuals_; }
        void fill_report_data();
        unsigned long get_report_indiv(const std::string &key) const;
        unsigned long get_report_house(const std::string &key) const;
        std::pair<unsigned short, unsigned long> get_report_hh_size(const std::string &key) const;
        bool is_empty() const { return house_vector_.empty(); }
        const str_vec &headers() const { return headers_; }
        Household copy(size_t H) const { return house_vector_[H]; }

        /////////////////////////////////////////////////////////////////////////////////////////
        ////////////////////////////  MCMC : Impl : DataWide : setter  ///////////////////////////
        /////////////////////////////////////////////////////////////////////////////////////////

        void build(const dict_str_vec &input);
        void build_multiplicity(const dict_str_vec &input);
        void init_cache_integer_vector(const std::string &e);
        void init_cache_size_t_vector(const std::string &e);
        void init_cache_double_vector(const std::string &e);
        void init_cache_double_matrix(const std::string &e);

        void init_draws(unsigned long sample_size);
        void update_draws(unsigned long smpl_count);
        void clear_draws();

        /////////////////////////////////////////////////////////////////////////////////////////
        ///////////////////////////  MCMC : Impl : DataWide : iterators  /////////////////////////
        /////////////////////////////////////////////////////////////////////////////////////////

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

        /////////////////////////////////////////////////////////////////////////////////////////
        /////////////////////////  MCMC : Impl : DataWide : import/export  ///////////////////////
        /////////////////////////////////////////////////////////////////////////////////////////

#ifdef _HOUSEPI_RCPP
        void rcpp_import(SEXP input);
        void rcpp_import_multiplicity(SEXP input);

        Rcpp::DataFrame rcpp_snapshot(const std::string &house_id) const;
        Rcpp::DataFrame rcpp_export(const str_vec &columns, bool is_status) const;
        Rcpp::DataFrame rcpp_summarise_cases() const;
        Rcpp::DataFrame rcpp_summarise_house_sizes() const;
        Rcpp::NumericMatrix rcpp_export_draws() const;

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

        std::map<std::string, unsigned long> report_indiv_ = {};
        std::map<std::string, unsigned long> report_house_ = {};
        std::map<std::string, std::pair<unsigned short, unsigned long>> report_size_ = {};

        std::ofstream aux_ofstream_;
    } dataset;

    /////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////  MCMC : Impl : DataLong : struct  ///////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////

    /// Implementation of Individual Datalist.
    struct DataLong
    {
        /////////////////////////////////////////////////////////////////////////////////////////
        //////////////////////////  MCMC : Impl : DataLong : constructor  ///////////////////////
        /////////////////////////////////////////////////////////////////////////////////////////

        DataLong(MCMC::Impl &ptr)
            : pImpl_(ptr) {}
        ~DataLong() = default;

        /////////////////////////////////////////////////////////////////////////////////////////
        ////////////////////////////  MCMC : Impl : DataLong : getter  //////////////////////////
        /////////////////////////////////////////////////////////////////////////////////////////

        bool is_empty() const { return indiv_datalist_.empty(); }
        const str_vec &headers() const { return headers_; }

        /////////////////////////////////////////////////////////////////////////////////////////
        ////////////////////////////  MCMC : Impl : DataLong : setter  //////////////////////////
        /////////////////////////////////////////////////////////////////////////////////////////

        void build(const dict_str_vec &input);

        /////////////////////////////////////////////////////////////////////////////////////////
        /////////////////////////  MCMC : Impl : DataLong : import/export  //////////////////////
        /////////////////////////////////////////////////////////////////////////////////////////

#ifdef _HOUSEPI_RCPP
        void rcpp_import(SEXP input);

        Rcpp::DataFrame rcpp_snapshot(const std::string &house_id) const;
        Rcpp::DataFrame rcpp_export() const;
#else
        void read(const std::string &pathtofile);
#endif
    private:
        const std::string section_name_ = "datalong";
        MCMC::Impl &pImpl_; // pointer to external class
        str_vec headers_ = {};
        std::vector<u_dict<dbl_vec>> indiv_datalist_ = {};
    } datalong;

    /////////////////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////  MCMC : Impl : GroupRates : class  ///////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////

    /// Implementation of Group Rates.
    struct GroupRates
    {
        /////////////////////////////////////////////////////////////////////////////////////////
        //////////////////////////  MCMC : Impl : GroupRates : constructor  /////////////////////
        /////////////////////////////////////////////////////////////////////////////////////////

        GroupRates(MCMC::Impl &ptr)
            : pImpl_(ptr) {}
        ~GroupRates() = default;

        /////////////////////////////////////////////////////////////////////////////////////////
        ////////////////////////////  MCMC : Impl : GroupRates : getter  ////////////////////////
        /////////////////////////////////////////////////////////////////////////////////////////

        bool is_empty() const { return headers_.empty(); }

        /////////////////////////////////////////////////////////////////////////////////////////
        ////////////////////////////  MCMC : Impl : GroupRates : setter  ////////////////////////
        /////////////////////////////////////////////////////////////////////////////////////////

        void build(const dict_str_vec &input);

        /////////////////////////////////////////////////////////////////////////////////////////
        ////////////////////////  MCMC : Impl : GroupRates : import/export  /////////////////////
        /////////////////////////////////////////////////////////////////////////////////////////

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

        void init_rate_as_ptrs_(Individual &member,
                                const u_dict<str_vec> &these_rates);
    } group_rates;

    /////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////  MCMC : Impl : Augmented : class  ///////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////

    /// Implementation of Augmented Variables.
    struct Augmented
    {
        dict_str_vec deftable = {};
        const str_vec file_headers = {"name", "status",
                                      "func", "fcumul",
                                      "farg_key", "farg_val", "farg_at",
                                      "value", "value_at", "init_rand",
                                      "const", "discr",
                                      "model", "arg_key", "arg_val", "arg_at",
                                      "min", "min_at",
                                      "max", "max_at",
                                      "rule", "scale"};

        str_vec list_model = {};
        str_vec list_const = {};
        str_vec list_sample = {};
        str_vec list_func = {};
        str_vec list_draws = {};
        str_vec list_full = {};
        str_vec list_export = {};
        std::map<int, str_vec> status_map_model = {};
        std::map<int, str_vec> status_map_infer = {};

        // optional augmented posterior sample
        std::vector<const Individual *> smpl_individuals = {};
        std::vector<u_dict<dbl_vec>> smpl_values = {};

        /////////////////////////////////////////////////////////////////////////////////////////
        //////////////////////////  MCMC : Impl : Augmented : constructor  //////////////////////
        /////////////////////////////////////////////////////////////////////////////////////////

        Augmented(MCMC::Impl &ptr)
            : pImpl_(ptr) {}
        ~Augmented() = default;

        /////////////////////////////////////////////////////////////////////////////////////////
        ////////////////////////////  MCMC : Impl : Augmented : setter  /////////////////////////
        /////////////////////////////////////////////////////////////////////////////////////////

        void build();
        void init();

        void clear_all(Individual &member);
        void init_all(Individual &member);

        void init_draws(unsigned long sample_size);
        void update_draws(unsigned long smpl_count);
        void clear_draws();

        /////////////////////////////////////////////////////////////////////////////////////////
        ////////////////////////////  MCMC : Impl : Augmented : getter  /////////////////////////
        /////////////////////////////////////////////////////////////////////////////////////////

        /// Flag if there are ant augmented variable definitions.
        bool is_any() const { return !(list_sample.empty() && list_func.empty() && list_const.empty()); }
        /// Flag if there are non-constant augmented variable definitions.
        bool is_any_sample() const { return !list_sample.empty(); }
        /// Flag if there are augmented variable definitions with a model.
        bool is_any_model() const { return !list_model.empty(); }
        /// Flag if augmented variable definitions are generic.
        bool is_generic() const { return status_map_infer.empty(); }
        /// Flag if augmented draws are stored.
        bool is_any_draw() const { return !smpl_individuals.empty(); }

#ifdef _HOUSEPI_RCPP
        void rcpp_import(SEXP input);
        Rcpp::DataFrame rcpp_export_data_wide() const;
        Rcpp::DataFrame rcpp_snapshot(const std::string &house_id) const;
#else
        void read(const std::string &pathtofile);
#endif

    private:
        const std::string section_name_ = "augmented variables";
        MCMC::Impl &pImpl_; // pointer to external class

        dict_str_vec check_name_status_(const dict_str_vec &deftable);
        void allocate_univocal_defs_(const dict_str_vec &check);
        void set_id_draws_(str_vec &id_draws);

        double type_to_double_(const Individual &member,
                               const std::string &tkn_value,
                               const std::string &tkn_type);
        void init_attribute_list_(bool is_var_model,
                                  Individual &member,
                                  const std::string &a,
                                  const str_vec &all_keys,
                                  const str_vec &all_values,
                                  const str_vec &all_types);
        Attribute read_single_attribute_(Individual &member,
                                         const std::string &a,
                                         const std::string &value,
                                         const std::string &type);

        void build_function_(size_t k);
        void build_function_identity_helper_(size_t k);
        void build_var_to_sample_(size_t k);
        void build_var_with_model_(size_t k);

        void init_single_var_dependencies_(Individual &member, size_t k);
    } augm;

    /////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////  MCMC : Impl : InfecProfile : class  //////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////

    /// Implementation of Infectivity Profile.
    struct InfecProfile
    {
    public:
        /////////////////////////////////////////////////////////////////////////////////////////
        ////////////////////////  MCMC : Impl : InfecProfile : constructor  /////////////////////
        /////////////////////////////////////////////////////////////////////////////////////////

        InfecProfile(MCMC::Impl &ptr)
            : pImpl_(ptr),
              profile_(std::make_unique<Distribution>()) {}
        ~InfecProfile() = default;

        /////////////////////////////////////////////////////////////////////////////////////////
        ///////////////////////////  MCMC : Impl : InfecProfile : setter  ///////////////////////
        /////////////////////////////////////////////////////////////////////////////////////////

        void build(const dict_str &deflist);
        void init_template(dict_str deflist);
        void append_augm_definition();
        void update(bool is_cache_depen);
        void recover(bool is_cache_depen);

        /////////////////////////////////////////////////////////////////////////////////////////
        ///////////////////////////  MCMC : Impl : InfecProfile : getter  ///////////////////////
        /////////////////////////////////////////////////////////////////////////////////////////

        bool is_lookup() const { return inst_ == &InfecProfile::do_inst_lookup_; }
        bool is_augmented() const { return inst_ == &InfecProfile::do_inst_augm_; }

        double inst(const Individual *member, double t) { return (this->*inst_)(member, t); }
        double cumul(const Individual *member, double t) { return (this->*cumul_)(member, t); }

#ifdef _HOUSEPI_RCPP
        void rcpp_import(SEXP input);
        Rcpp::List rcpp_export() const;
        Rcpp::List rcpp_snapshot(const std::string &house_id) const;
#else
        void read(const std::string &pathtofile);
#endif

    private:
        const std::string section_name_ = "infectivity profile";
        MCMC::Impl &pImpl_; // pointer to external class
        std::shared_ptr<Distribution> profile_;

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
        void clear_lookup_();

        void build_generic_profile_(const dict_str &deflist);
        void build_augmented_profile_(const dict_str &deflist);

        void init_single_attribute_(const std::string &key, const std::string &val, const std::string &at);
        void init_arguments_(const str_vec &keys, const str_vec &vals, const str_vec &typs);

        //===== Functional methods for getter methods =====//
        double (InfecProfile::*inst_)(const Individual *, double) = &InfecProfile::do_inst_default_;
        double do_inst_default_(const Individual *, double t) { return profile_->pdf(t); }
        double do_inst_augm_(const Individual *i, double t) { return i->get_augm_variable(HSP_PROFILE_INST)(t); }
        double do_inst_lookup_(const Individual *, double t);

        double (InfecProfile::*cumul_)(const Individual *, double) = &InfecProfile::do_cumul_default_;
        double do_cumul_default_(const Individual *, double t) { return profile_->cdf(t); }
        double do_cumul_augm_(const Individual *i, double t) { return i->get_augm_variable(HSP_PROFILE_CUMUL)(t); }
        double do_cumul_lookup_(const Individual *, double t);
    } infec_profile;

    /////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////  MCMC : Impl : Metropolis-Hastings  ///////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////

    double evaluate();
    void infer();

    /// Release memory.
    void clear_memory() { clear_posterior_sample_(); }

#ifdef _HOUSEPI_RCPP
    //----- C++ Print -----//
    Rcpp::List simulate(size_t n_rep);
    Rcpp::DataFrame simulate_repetition();

    Rcpp::List evidence(const std::string &scheme, unsigned long n_rungs);

    //----- Rcpp EXPORTS -----//
    Rcpp::NumericMatrix rcpp_export_draws() const;
    Rcpp::List rcpp_export_house_snapshot(const std::string &house_id) const;
#else
    //----- C++ Print -----//
    void print_infer_state(const std::string &output_filename = HSP_NA_STRING) const;
    void print_simulated_data_wide(const std::string &output_filename = HSP_NA_STRING);
    void print_evidence(const std::string &scheme, unsigned long n_rungs);
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

    void init_model_for_inference_() { callback_init_inference_(); }
    void init_model_for_simulation_() { callback_init_simulation_(); }
    void simulate_single_household_(Household &h) { callback_simulate_house_(h); }
    void update_cache_dependecies_(Household &h, Individual &i) { callback_update_cache_(h, i); }

    double compute_and_set_household_obs_loglike_(Household &h)
    {
        const double output = callback_house_loglike_(h);
        h.set_obs_loglike(output);
        return output;
    }

    /////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////  MCMC : Impl : Metropolis-Hastings [private]  //////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////

    double full_loglike_ = 0.; /// current full loglikelihood
    double full_buffer_ = 0.;  /// buffer loglikelihood

    unsigned long smpl_count_ = 0;
    dbl_vec smpl_full_loglike_ = {};

    void assert_proper_dataset_();
    void init_monte_carlo_();

    /// Set new full loglikelihood value and buffer old value.
    void set_full_loglike_(double new_ll)
    {
        full_buffer_ = full_loglike_;
        full_loglike_ = new_ll;
    }

    /// Increment full loglikelohood by a delta, buffer old value.
    void increment_full_loglike_(double delta) { set_full_loglike_(full_loglike_ + delta); }

    void init_move_counts_();

    void update_full_loglike_(bool cache_depen);
    double update_house_loglike_(Household &house, bool is_augm_ll, bool is_cache_depen);

    void run_monte_carlo_();
    void move_all_();
    void move_all_and_adapt_params_(unsigned long step);

    void init_posterior_sample_();
    void update_posterior_sample_();
    void clear_posterior_sample_();

    //----- Metropolis - parameter moves -----//
    bool move_single_parameter_(const std::string &p, bool cache_depen);
    bool propose_valid_param_move_(const std::string &p);
    bool accept_valid_param_move_(const std::string &p);
    void recover_parameter_(const std::string &p, bool cache_depen);

    //----- Metropolis : augmented moves -----//
    void move_all_augm_();
    double move_augm_single_house_(Household &house);
    double metropolis_augm_single_move_(Household &house, Individual &member, const std::string &a);

    /// Compute statistics from sample of full loglikelihood.
    std::pair<double, double> smpl_full_loglike_stats_() const { return tools::compute_stats(smpl_full_loglike_); }

    //----- Metropolis : power posteriors -----//
    std::vector<dbl_vec> power_posteriors_(const std::string &scheme, unsigned long N);

    //====================//
    //===== SIMULATE =====//
    //====================//

    void aux_init_simulate_();
};

#endif // HOUSEPI__MCMC__IMPL__H
