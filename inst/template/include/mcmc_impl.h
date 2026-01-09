/* housepi - Epidemiological Analysis of Household Outbreaks.
 *
 * Copyright © 2026 Institut Pasteur (Paris), CNRS and INSERM.
 * Distributed under the terms of the GNU General Public License (GPLv3).
 * See the LICENSE file for details.
 */

#ifndef HOUSEPI__MCMC__IMPL__H
#define HOUSEPI__MCMC__IMPL__H

#include <memory>
#include <functional>

#include "tools/types.h"
#include "tools/na_utils.h"
#include "tools/messages.h"
#include "tools/atomic.h"

#include "attribute.h"
#include "variable.h"
#include "mcmc.h"

/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////  *** MCMC : Impl : class ***  ////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

namespace hsp
{
    namespace atomic
    {
        struct Error;
    }
}

class Individual;
class Household;
struct Setup;
struct DataWide;
struct DataLong;
class Distribution;

class MCMC::Impl
{
public:
    std::unique_ptr<hsp::atomic::Error> atomic_error;
    std::shared_ptr<Setup> setup;
    std::shared_ptr<DataWide> data_wide;
    std::shared_ptr<DataLong> data_long;

    /////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////  MCMC : Impl : constructor  ///////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////

    Impl(std::function<void()> func_init_inference,
         std::function<void()> func_init_simulation,
         std::function<double(Household &)> func_household_log_lik,
         std::function<void(Household &)> func_simulate_house,
         std::function<void(Household &, Individual &)> func_update_cache);
    ~Impl();

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
        void user_reset();
        /// Load parameters values of the n-th posterrior sample.
        void load_sample_values(size_t n);

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

        InfecProfile(MCMC::Impl &ptr);
        ~InfecProfile();

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

        void build_generic_profile_(const u_dict<str_vec> &input);
        void build_augmented_profile_(const u_dict<str_vec> &input);
        std::unique_ptr<Distribution> helper_build_profile_(const u_dict<std::string> &input);

        //===== Functional methods for getter methods =====//
        double (InfecProfile::*inst_)(const Individual *, double) = &InfecProfile::do_inst_default_;
        double do_inst_default_(const Individual *, double t);
        double do_inst_lookup_(const Individual *, double t);
        double do_inst_augm_(const Individual *i, double t);

        double (InfecProfile::*cumul_)(const Individual *, double) = &InfecProfile::do_cumul_default_;
        double do_cumul_default_(const Individual *, double t);
        double do_cumul_lookup_(const Individual *, double t);
        double do_cumul_augm_(const Individual *i, double t);
    } infec_profile;

    /////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////  MCMC : Impl : Augmented : class  ///////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////

    /// Implementation of Augmented variables and functions.
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
    /// UPDATE: once I allow for multi-treading and a flag other than n_chains = 0, I can drop this.
    void helper_init_simulation();
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
    double materialize_house_obs_log_lik_(Household &h);

    /////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////  MCMC : Impl : Metropolis-Hastings [private]  //////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////

    double total_log_lik_ = 0.; /// current total log-likelihood
    double total_buffer_ = 0.;  /// buffer total log-likelihood

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

    //----- Metropolis - parameter moves -----//
    bool move_single_parameter_(const std::string &p, bool is_cache_depen);
    bool propose_valid_param_move_(const std::string &p);
    bool accept_valid_param_move_(const std::string &p);
    void recover_parameter_(const std::string &p, bool is_cache_depen);

    //----- Metropolis : augmented moves -----//
    void move_all_augm_();
    double atomic_move_augm_single_house_(size_t H);
    double metropolis_augm_single_move_(Household &house, Individual &member, const std::string &a);

    //----- Metropolis : power posteriors -----//
    void power_post_zero_(dbl_vec &itm, dbl_vec &avg, dbl_vec &var);
    void power_post_fractions_(dbl_vec &itm, dbl_vec &avg, dbl_vec &var);
    void power_post_adaptive_(dbl_vec &itm, dbl_vec &avg, dbl_vec &var);
};

#endif // HOUSEPI__MCMC__IMPL__H
