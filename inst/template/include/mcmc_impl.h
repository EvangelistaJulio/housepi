/*
 *        __
 *       / /_  ____  __  __ ____ __  ____  __
 *      / __ \/ __ \/ / / / ___/ _ \/ __ \/ /
 *     / / / / /_/ / /_/ (__  )  __/ /_/ / /
 *    /_/ /_/\____/\____/____/\___/ ____/_/
 *                               /_/
 *
 * housepi - Epidemiological Analysis of Household Outbreaks.
 *
 * Author: Francesco Camaglia
 * Copyright © 2025 Institut Pasteur (Paris), CNRS and INSERM.
 *
 * Distributed under the terms of the GNU General Public License (GPLv3).
 * See the LICENSE file for details.
 *
 * mcmc_impl.h
 *
 */

#ifndef HOUSEPI__MCMC__IMPL__H
#define HOUSEPI__MCMC__IMPL__H

#include <memory>
#include <functional>
#include <atomic>

#include "tools/types.h"
#include "tools/na_utils.h"
#include "tools/messages.h"
#include "tools/progress.h"
#include "tools/infprof.h"
#include "function.h"
#include "distribution.h"
#include "variable.h"
#include "household.h"
#include "mcmc.h"

#ifdef _OPENMP
#include <omp.h>
#endif

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
         std::function<double(Household &)> func_household_log_lik,
         std::function<void(Household &)> func_simulate_house,
         std::function<void(Household &, Individual &)> func_update_cache)
        : setup(), param(*this),
          data_wide(*this), multiplicity(*this), data_long(*this),
          group_rates(*this), infec_profile(*this), augm(*this),
          callback_init_inference_(std::move(func_init_inference)),
          callback_init_simulation_(std::move(func_init_simulation)),
          callback_house_log_lik_(std::move(func_household_log_lik)),
          callback_simulate_house_(std::move(func_simulate_house)),
          callback_update_cache_(std::move(func_update_cache)) {}
    ~Impl() = default;

    /////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////  MCMC : Impl : Setup : struct  /////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////

    /// Implementation of Setup.
    struct Setup
    {
        Setup() : unif_(0., 1.) {}
        ~Setup() = default;

        const str_vec file_headers = {"steps", "burnin", "thin",
                                      "adapt", "chains", "threads",
                                      "id_draws", "seed", "progress"};

        unsigned long n_steps = 1;  /// n. iterations
        unsigned long n_burnin = 0; /// n. burnin iterations
        unsigned long n_thin = 1;   /// n. thinning iterations
        unsigned long n_adapt = 0;  /// n. (finite) adaptive steps
        unsigned long n_chains = 1; /// n. chains
        unsigned n_threads = 1;     /// n. of OMP threads
        std::int32_t seed;          /// random generator seed
        unsigned long progress = 1; /// thinning seconds for progress bar
        unsigned long n_rungs = 0;  /// n. of rungs

        unsigned long c_sample_size = 0; /// single chain sample size
        unsigned long c_n_sweeps = 0;    /// single chain n. of sweeps
        double inv_T = 1.;               /// inverse temperature
        str_vec id_draws = {};           /// augmented trace to save

        unsigned long chain = 1; /// current MCMC chain
        unsigned long rung = 0;  /// current MCMC thermodynamic rung

        bool is_power_post() const { return n_rungs > 0; }
        std::string fmt_progtxt();

        void set_n_threads();
        void build(const u_dict<str_vec> &deftable);
        void assert_mcmc_setup();
        void drop_mcmc_adapt();

        // progress bar
        void prog_set_text(const std::string &text = hsp::na::STRING)
        {
            progressbar_.set_text(text == hsp::na::STRING ? fmt_progtxt() : text);
        }
        void prog_declare(const std::string &text = hsp::na::STRING)
        {
            prog_set_text(text);
            progressbar_.set_thintime(progress);
        }
        void prog_open_for_inference() { progressbar_.open(c_n_sweeps * n_chains); }
        void prog_open_for_power_post() { progressbar_.open(c_n_sweeps * n_chains * (n_rungs + 1)); }
        void prog_update() { progressbar_.update(); }
        void prog_close() { progressbar_.close(); }
        void prog_crash() { progressbar_.crash(); }

        // Pseudo Random Number Geenerator
        void set_prng_seed();
        rnd_gen &prng() { return (this->*get_prng_)(); }
        double unif() { return unif_(prng()); }

    private:
        hsp::progress::Bar progressbar_; // progress bar

        rnd_gen prng_; // Mersenne Twister 19937 generator (32 or 64 bit)
        rnd_gen &(Setup::*get_prng_)();
        rnd_gen &do_prng_unique_thread_() { return prng_; }
        std::uniform_real_distribution<double> unif_; // uniform sampler

#ifdef _OPENMP
        std::vector<rnd_gen> prng_parallel_ = {}; // parallel generators (for multi-threading)
        rnd_gen &do_prng_parallel_threads_() { return prng_parallel_[omp_get_thread_num()]; }
#endif
    } setup;

    /////////////////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////  MCMC : Impl : Parameter : struct  ///////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////

    /// Implementation of Parameters.
    struct Parameter
    {
        /// NOTE: list_const + list_sample = all parameters in dictionary_
        /// list of constant parameters
        str_vec list_const = {};
        /// list of parameters to sample
        str_vec list_sample = {};

        /// NOTE: list_sample_standard + list_sample_cached = list_sample
        str_vec list_sample_standard = {};
        str_vec list_sample_cached = {};

        /// NOTE: list_adaptable_sample_standard + list_non_adapt_sample_standard = list_sample_standard
        str_vec list_adaptable_sample_standard = {};
        str_vec list_non_adapt_sample_standard = {};
        /// NOTE: list_adaptable_sample_cached + list_non_adapt_sample_cached = list_sample_cached
        str_vec list_adaptable_sample_cached = {};
        str_vec list_non_adapt_sample_cached = {};

        u_dict<dbl_vec> smpl_values = {}; /// values of the sample
        u_dict<dbl_vec> smpl_priors = {}; /// prior values of the sample

        /////////////////////////////////////////////////////////////////////////////////////////
        /////////////////////////  MCMC : Impl : Parameter : constructor  ///////////////////////
        /////////////////////////////////////////////////////////////////////////////////////////

        Parameter(MCMC::Impl &ptr)
            : pImpl_(ptr) {}
        ~Parameter() = default;

        /////////////////////////////////////////////////////////////////////////////////////////
        ///////////////////////////  MCMC : Impl : Parameter : getter  //////////////////////////
        /////////////////////////////////////////////////////////////////////////////////////////

        Variable &at(const std::string &p) { return safe_at(dictionary_, p, "parameter"); }
        const Variable &at(const std::string &p) const { return safe_at(dictionary_, p, "parameter"); }
        // get value of the parameter
        double get(const std::string &p) const { return at(p).get(); }

        // flag to check if there are parameters to sample
        bool is_any_sample() const { return !list_sample.empty(); }

        /////////////////////////////////////////////////////////////////////////////////////////
        ///////////////////////////  MCMC : Impl : Parameter : setter  //////////////////////////
        /////////////////////////////////////////////////////////////////////////////////////////

        void build(const u_dict<str_vec> &deftable);
        void build_for_expec(const u_dict<dbl_vec> &deftable);
        void make_cache_dependent(const std::string &p);

        /// Randomly initialize values of the parameters to be inferred.
        void user_reset()
        {
            for (const std::string &p : list_sample)
                at(p).user_reset(pImpl_.setup.prng());
        }
        /// Load parameters values of the n-th posterrior sample.
        void load_sample_values(size_t n)
        {
            for (const auto &pair : smpl_values)
                at(pair.first).set(pair.second[n]);
        }

        void init_draws(unsigned long sample_size);
        void store_draws(unsigned long smpl_count);

    private:
        MCMC::Impl &pImpl_; // pointer to external class

        u_dict<Variable> dictionary_ = {}; // parameters dictionary

        void clear_lists_();
        void build_constant_parameter_(const u_dict<str_vec> &deftable, size_t k);
        void build_parameter_to_sample_(const u_dict<str_vec> &deftable, size_t k);
    } param;

    /////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////  MCMC : Impl : DataWide : struct  ////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////

    /// Implementation of individidual data wide format.
    struct DataWide
    {
        str_vec headers = {};

        // list of all infected statuses in the dataset
        int_vec all_infected_status_list = {};
        std::vector<dbl_vec> smpl_values = {};

        ReportIndividual report_indiv;
        ReportHousehold report_house;
        ReportHouseSizes report_hh_size;

        /////////////////////////////////////////////////////////////////////////////////////////
        /////////////////////////  MCMC : Impl : DataSet : constructor  /////////////////////////
        /////////////////////////////////////////////////////////////////////////////////////////

        DataWide(MCMC::Impl &ptr)
            : pImpl_(ptr) {}
        ~DataWide() = default;

        Household &operator()(size_t H) { return house_vector_[H]; }
        const Household &operator()(size_t H) const { return house_vector_.at(H); }
        Household &operator()(const std::string &house_id) { return (*this)(house_id_to_index_(house_id)); }
        const Household &operator()(const std::string &house_id) const { return (*this)(house_id_to_index_(house_id)); }

        /////////////////////////////////////////////////////////////////////////////////////////
        ////////////////////////////  MCMC : Impl : DataWide : getter  ///////////////////////////
        /////////////////////////////////////////////////////////////////////////////////////////

        bool is_empty() const { return house_vector_.empty(); }
        Household copy(size_t H) const { return house_vector_[H]; }

        /////////////////////////////////////////////////////////////////////////////////////////
        ////////////////////////////  MCMC : Impl : DataWide : setter  ///////////////////////////
        /////////////////////////////////////////////////////////////////////////////////////////

        void build(const u_dict<str_vec> &deftable);

        void init_draws(unsigned long sample_size);
        void store_draws(unsigned long smpl_count);

        /////////////////////////////////////////////////////////////////////////////////////////
        ///////////////////////////  MCMC : Impl : DataWide : iterators  /////////////////////////
        /////////////////////////////////////////////////////////////////////////////////////////

        auto begin() { return house_vector_.begin(); }
        auto begin() const { return house_vector_.begin(); }
        auto end() { return house_vector_.end(); }
        auto end() const { return house_vector_.end(); }

        template <typename T>
        class View
        {
            // conditional to select the vector type if constant or not
            using VectorRef = std::conditional_t<std::is_const<T>::value,
                                                 const std::vector<Household> &,
                                                 std::vector<Household> &>;

        public:
            VectorRef vec;
            const std::vector<size_t> &indexes; /// NOTE: no need to copy here

            View(VectorRef m, const std::vector<size_t> &i)
                : vec(m), indexes(i) {}

            class iterator
            {
            public:
                iterator(View *v, size_t p) : view_(v), pos_(p) {}

                T &operator*() const { return view_->vec[view_->indexes[pos_]]; }
                T *operator->() const { return &view_->vec[view_->indexes[pos_]]; }

                iterator &operator++()
                {
                    ++pos_;
                    return *this;
                }
                iterator operator++(int)
                {
                    iterator tmp = *this;
                    ++(*this);
                    return tmp;
                }

                bool operator==(const iterator &other) const { return pos_ == other.pos_; }
                bool operator!=(const iterator &other) const { return !(*this == other); }

            private:
                View *view_;
                size_t pos_;
            };

            iterator begin() { return iterator(this, 0); }
            iterator end() { return iterator(this, indexes.size()); }
        };

        View<Household> restricted() { return View<Household>(house_vector_, pImpl_.multiplicity.head_class_house_indexes); }
        View<const Household> restricted() const { return View<const Household>(house_vector_, pImpl_.multiplicity.head_class_house_indexes); }
        Household &restricted(size_t R) { return house_vector_[pImpl_.multiplicity.head_class_house_indexes[R]]; }

    private:
        MCMC::Impl &pImpl_; // pointer to external class

        std::vector<Household> house_vector_ = {};
        std::unordered_map<std::string, size_t> hhid_to_index_ = {};

        size_t house_id_to_index_(const std::string &house_id) const
        {
            if (hhid_to_index_.count(house_id))
                return hhid_to_index_.at(house_id);
            throw std::invalid_argument(
                "Index '" + house_id + "' did not match any household in the data.");
        }

        /// Clear data report.
        void build_helper_clear_data_();
        /// Use the loaded data to fill-in the respective report.
        void build_helper_fill_report_data_();
    } data_wide;

    /////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////  MCMC : Impl : Multiplicity : struct  /////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////

    /// Implementation of multiplicity optimization.
    struct Multiplicity
    {
        std::vector<size_t> head_class_house_indexes = {};

        /////////////////////////////////////////////////////////////////////////////////////////
        ///////////////////////  MCMC : Impl : Multiplicity : constructor  //////////////////////
        /////////////////////////////////////////////////////////////////////////////////////////

        Multiplicity(MCMC::Impl &ptr)
            : pImpl_(ptr) {}
        ~Multiplicity() = default;

        void build(const u_dict<str_vec> &deftable);

    private:
        MCMC::Impl &pImpl_; // pointer to external class

    } multiplicity;

    /////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////  MCMC : Impl : DataLong : struct  ///////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////

    /// Implementation of additional individual data long format.
    struct DataLong
    {
        str_vec headers = {};
        std::vector<u_dict<dbl_vec>> indiv_datalist = {};

        /////////////////////////////////////////////////////////////////////////////////////////
        //////////////////////////  MCMC : Impl : DataLong : constructor  ///////////////////////
        /////////////////////////////////////////////////////////////////////////////////////////

        DataLong(MCMC::Impl &ptr)
            : pImpl_(ptr) {}
        ~DataLong() = default;

        /////////////////////////////////////////////////////////////////////////////////////////
        ////////////////////////////  MCMC : Impl : DataLong : setter  //////////////////////////
        /////////////////////////////////////////////////////////////////////////////////////////

        void build(const u_dict<str_vec> &deftable);

    private:
        MCMC::Impl &pImpl_; // pointer to external class
    } data_long;

    /////////////////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////  MCMC : Impl : GroupRates : class  ///////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////

    /// Implementation of Group Rates.
    struct GroupRates
    {
        str_vec headers;

        /////////////////////////////////////////////////////////////////////////////////////////
        //////////////////////////  MCMC : Impl : GroupRates : constructor  /////////////////////
        /////////////////////////////////////////////////////////////////////////////////////////

        GroupRates(MCMC::Impl &ptr)
            : pImpl_(ptr) {}
        ~GroupRates() = default;

        /////////////////////////////////////////////////////////////////////////////////////////
        ////////////////////////////  MCMC : Impl : GroupRates : setter  ////////////////////////
        /////////////////////////////////////////////////////////////////////////////////////////

        void build(const u_dict<str_vec> &deftable);

    private:
        MCMC::Impl &pImpl_; // pointer to external class

        const str_vec file_headers_ = {"group_rate", "param_name"};

        str_vec helper_assert_definitions_(const u_dict<str_vec> &deftable);
        void helper_init_all_rates_(Individual &member, const u_dict<str_vec> &these_rates);
    } group_rates;

    /////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////  MCMC : Impl : InfecProfile : class  //////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////

    /// Implementation of Infectivity Profile.
    struct InfecProfile
    {
    public:
        std::unique_ptr<Distribution> generic_profile;

        /////////////////////////////////////////////////////////////////////////////////////////
        ////////////////////////  MCMC : Impl : InfecProfile : constructor  /////////////////////
        /////////////////////////////////////////////////////////////////////////////////////////

        InfecProfile(MCMC::Impl &ptr)
            : generic_profile(std::make_unique<Distribution>()),
              pImpl_(ptr) {}
        ~InfecProfile() = default;

        /////////////////////////////////////////////////////////////////////////////////////////
        ///////////////////////////  MCMC : Impl : InfecProfile : setter  ///////////////////////
        /////////////////////////////////////////////////////////////////////////////////////////

        void build(const u_dict<str_vec> &deftable);
        void init_template(u_dict<str_vec> deftable);
        void update();
        void recover();
        void request_discrete() { is_discrete_requested_ = true; };
        void transfer_augm_definitions(u_dict<str_vec> &target) const;

        /////////////////////////////////////////////////////////////////////////////////////////
        ///////////////////////////  MCMC : Impl : InfecProfile : getter  ///////////////////////
        /////////////////////////////////////////////////////////////////////////////////////////

        bool is_lookup() const { return inst_ == &InfecProfile::do_inst_lookup_; }
        bool is_augmented() const { return inst_ == &InfecProfile::do_inst_augm_; }

        double inst(const Individual *member, double t) { return (this->*inst_)(member, t); }
        double cumul(const Individual *member, double t) { return (this->*cumul_)(member, t); }

    private:
        MCMC::Impl &pImpl_; // pointer to external class

        // lookup table
        bool is_discrete_requested_ = false;
        size_t size_lut_ = 0;
        dbl_vec lookup_inst_ = {};
        dbl_vec lookup_cumul_ = {};
        dbl_vec lookup_buffer_inst_ = {};
        dbl_vec lookup_buffer_cumul_ = {};
        // augmented function
        u_dict<str_vec> augm_profile_definitions_ = {};

        void init_lookup_();
        void clear_lookup_();

        void build_generic_profile_(const u_dict<str_vec> &deftable);
        void build_augmented_profile_(const u_dict<str_vec> &deftable);

        //===== Functional methods for getter methods =====//
        double (InfecProfile::*inst_)(const Individual *, double) = &InfecProfile::do_inst_default_;
        double do_inst_default_(const Individual *, double t) { return generic_profile->pdf(t); }
        double do_inst_augm_(const Individual *i, double t) { return i->get_augm_variable(hsp::infprof::INSTANTANEOUS)(t); }
        double do_inst_lookup_(const Individual *, double t);

        double (InfecProfile::*cumul_)(const Individual *, double) = &InfecProfile::do_cumul_default_;
        double do_cumul_default_(const Individual *, double t) { return generic_profile->cdf(t); }
        double do_cumul_augm_(const Individual *i, double t) { return i->get_augm_variable(hsp::infprof::CUMULATIVE)(t); }
        double do_cumul_lookup_(const Individual *, double t);
    } infec_profile;

    /////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////  MCMC : Impl : Augmented : class  ///////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////

    /// Implementation of Augmented Variables.
    struct Augmented
    {
        str_vec list_draws = {};
        str_vec list_full = {};
        str_vec list_export = {};
        str_vec user_list_discrete = {};

        u_dict<std::vector<size_t>> status_to_def_row = {}; // definition rows by status
        std::unordered_map<int, u_dict<str_vec>> status_map_model = {};
        std::unordered_map<int, str_vec> status_map_sample = {};

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

        void build(const u_dict<str_vec> &deftable);
        void init();
        void user_reset();

        void clear_all(Individual &member);
        void init_all(Individual &member);
        void init_missing(Individual &member);
        void init_model_logratios(Individual &member);
        void user_reset(Individual &member);

        void init_draws(unsigned long sample_size);
        void store_draws(unsigned long smpl_count);

        /////////////////////////////////////////////////////////////////////////////////////////
        ////////////////////////////  MCMC : Impl : Augmented : getter  /////////////////////////
        /////////////////////////////////////////////////////////////////////////////////////////

        /// Flag if there are any augmented definitions.
        bool is_any() const { return !list_full.empty(); }
        /// Flag if there are non-constant augmented definitions.
        bool is_any_sample() const { return !status_map_sample.empty(); }
        /// Flag if there are augmented definitions with a model.
        bool is_any_model() const { return !status_map_model.empty(); }
        /// Flag if augmented draws are stored.
        bool is_any_draw() const { return !smpl_individuals.empty(); }

    private:
        MCMC::Impl &pImpl_;             // pointer to external class
        u_dict<str_vec> deftable_ = {}; // store definition table

        void helper_assert_no_reserved_();
        void helper_build_definition_by_status_();
        void helper_build_definition_order_();
        void set_definition_dependency_order_for_(const std::string &s);
        std::vector<size_t> get_unsorted_definition_rows_(const std::string &s, bool is_pre_process) const;
        void set_sample_map_for_(const std::string &s, size_t k);
        void set_function_or_identity_for_(size_t k);
        void helper_build_model_status_maps_();

        void helper_init_attributes_(Individual &member, size_t k);
        std::vector<size_t> init_get_definition_order_(Individual &member);
        void init_model_attributes_(Individual &member, size_t k);
        void init_func_attributes_(Individual &member, size_t k);
        void helper_init_id_draws_();
    } augm;

    /////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////  MCMC : Impl : Errors  /////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////

    struct Errors
    {
    public:
        /// Reset error of the critical scope.
        void reset()
        {
            has_error_.exchange(false);
            message_.clear();
        }
        // Store first error encountered (thread-safe)
        void set_if_first(const std::string &message)
        {
            if (!has_error_.exchange(true, std::memory_order_relaxed))
            {
#pragma omp critical
                message_ = message;
            }
        }
        void set_if_first(const std::exception &e) { set_if_first(e.what()); }
        // Wheter it should be skipped or not.
        bool is_pending()
        {
            return has_error_.load(std::memory_order_relaxed);
        }
        /// If an error was caught in critical scope, throw an exception.
        void throw_exception() const
        {
            if (has_error_)
                throw std::runtime_error(message_);
        }

    private:
        std::atomic<bool> has_error_{false}; /// Wheter error was raised during multi-threading
        std::string message_ = {};           /// Error message for last exception during multi-threading
    } atomic_error;

    /////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////  MCMC : Impl : Metropolis-Hastings  ///////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////

    unsigned long smpl_count = 0;    /// count of store samples
    dbl_vec smpl_total_log_lik = {}; /// vector of total log-likelihood samples

    void init_monte_carlo();
    void init_evaluate();
    void init_model_for_inference() { callback_init_inference_(); }
    void run_monte_carlo();
    std::vector<dbl_vec> power_posteriors(const std::string &scheme, unsigned long N);
    double evaluate_single_house_log_lik(const std::string &house_id);
    double update_total_log_lik(bool is_cache_depen = true);

    /////////////////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////  MCMC : Impl : Simulation  ///////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////

    /// Auxiliary method to initialize simulations.
    // UPDATE: once I allow for multi-treading and a flag other than n_chains = 0, I can drop this.
    void helper_init_simulation()
    {
        assert_non_empty_data_wide_(); // check if data is empty
        setup.n_chains = 0;            // use 0 chains as a state variable
        setup.n_threads = 1;           // impose no parallelization
        setup.set_n_threads();         // change accordingly
        init_model_for_simulation_();  // initialize model for simulation
    }
    void simulate_single_household_in_place(Household &h) { callback_simulate_house_(h); }

private:
    /////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////  MCMC : Impl : definition helpers [private]  //////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////

    mutable Attribute storage_; /// Auxiliary storage attribute used for reading.
    const Attribute &read_numeric_(const hsp::parse::Definition::Entry &node, const Individual &member = Individual()) const;
    const Attribute &read_numeric_(const std::string &expr, const Individual &member = Individual()) const;
    const Attribute &read_attribute_(const hsp::parse::Definition::Entry &node, const Individual &member = Individual()) const;
    const Attribute &read_attribute_(const std::string &expr, const Individual &member = Individual()) const;
    void make_param_cache_depen_(const hsp::parse::Definition::Entry &node);

    /////////////////////////////////////////////////////////////////////////////////////////////
    //////////////////////////  MCMC : Impl : Callbacks Virtual Methods  ////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////

    std::function<void()> callback_init_inference_;
    std::function<void()> callback_init_simulation_;
    std::function<double(Household &)> callback_house_log_lik_;
    std::function<void(Household &)> callback_simulate_house_;
    std::function<void(Household &, Individual &)> callback_update_cache_;

    void init_model_for_simulation_() { callback_init_simulation_(); }
    void update_cache_dependecies_(Household &h, Individual &i) { callback_update_cache_(h, i); }

    double materialize_house_obs_log_lik_(Household &h)
    {
        const double output = callback_house_log_lik_(h);
        h.set_obs_log_lik(output);
        return output;
    }

    /////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////  MCMC : Impl : Metropolis-Hastings [private]  //////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////

    double total_log_lik_ = 0.; /// current total log-likelihood
    double total_buffer_ = 0.;  /// buffer total log-likelihood

    /// Assert data_wide is not empty.
    void assert_non_empty_data_wide_()
    {
        if (data_wide.report_house.total == 0)
            throw std::length_error("Unable to proceed: no data was found.");
    }

    /// Warn if all individuals have missing status.
    void warn_if_all_missing_status_()
    {
        if (data_wide.report_indiv.total == data_wide.report_indiv.missing)
            hsp::msg::warning("All individuals in the data had missing status.");
    }

    /// Set new total log-likelihood value and buffer old value.
    void set_total_log_lik_(double new_ll)
    {
        total_buffer_ = total_log_lik_;
        total_log_lik_ = new_ll;
    }

    /// Increment total log-likelohood by a delta, buffer old value.
    void increment_total_log_lik_(double delta) { set_total_log_lik_(total_log_lik_ + delta); }

    void init_move_counts_();

    double atomic_materialize_house_log_lik_(Household &house, bool is_augm_ll, bool is_cache_depen = true);

    void move_all_();
    void adapt_proposal_of_sampled_params_(unsigned long step);

    void init_posterior_sample_();
    void store_posterior_sample_();
    void clear_posterior_sample_();

    //----- Metropolis - parameter moves -----//
    bool move_single_parameter_(const std::string &p, bool is_cache_depen);
    bool propose_valid_param_move_(const std::string &p);
    bool accept_valid_param_move_(const std::string &p);
    void recover_parameter_(const std::string &p, bool is_cache_depen);

    //----- Metropolis : augmented moves -----//
    void move_all_augm_();
    double atomic_move_augm_single_house_(Household &house);
    double metropolis_augm_single_move_(Household &house, Individual &member, const std::string &a);

    //----- Metropolis : power posteriors -----//
    void power_post_zero_(dbl_vec &itm, dbl_vec &avg, dbl_vec &var);
    void power_post_fractions_(dbl_vec &itm, dbl_vec &avg, dbl_vec &var);
    void power_post_adaptive_(dbl_vec &itm, dbl_vec &avg, dbl_vec &var);
};

#endif // HOUSEPI__MCMC__IMPL__H
