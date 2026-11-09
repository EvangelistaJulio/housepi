/* housepi - Epidemiological Analysis of Household Outbreaks.
 * Copyright 2024-2026 Institut Pasteur (Paris), France.
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#pragma once

#include <string>
#include <utility>
#include <vector>

#include "tools/messages.h"
#include "tools/vector.h"
#include "data/individual.h"
#include "mcmc/base.h"

//////////////////////////////  ModelBase  //////////////////////////////

namespace housepi
{
    class Individual;
    class Household;

#define MAKE_HOUSEPI_CACHE_FORWARDER(NAME)                                    \
    auto &cache_##NAME(const Household &house, size_t e) noexcept             \
    {                                                                         \
        return mcmc_->data.state.cache_##NAME(house, e);                      \
    }                                                                         \
    const auto &cache_##NAME(const Household &house, size_t e) const noexcept \
    {                                                                         \
        return mcmc_->data.state.cache_##NAME(house, e);                      \
    }

    /// Named quantities of one household that a model exports, indexed by member.
    struct HouseQuantities
    {
        /// Collect quantities of a household of `n` members.
        explicit HouseQuantities(size_t n) noexcept : n_(n) {}

        /// Add `values` as `name`, one per member.
        void vector(std::string name, std::vector<double> values)
        {
            if (values.size() != n_)
                msg::stop("Exported vector " + msg::val(name) + " is not sized by household members.");
            vectors.emplace_back(std::move(name), std::move(values));
        }
        /// Add `values` as `name`, one row and one column per member.
        void matrix(std::string name, housepi::matrix values)
        {
            bool sized = values.size() == n_;
            for (const std::vector<double> &row : values)
                sized = sized && row.size() == n_;
            if (!sized)
                msg::stop("Exported matrix " + msg::val(name) + " is not sized by household members.");
            matrices.emplace_back(std::move(name), std::move(values));
        }

        /// Vectors, by name.
        std::vector<std::pair<std::string, std::vector<double>>> vectors = {};
        /// Matrices, by name.
        std::vector<std::pair<std::string, housepi::matrix>> matrices = {};

    private:
        /// Household size.
        size_t n_;
    };

    /// Transmission model: per-household log-likelihood and simulation hooks.
    ///
    /// * built-in and custom models derive from it through the `HOUSEPI_*` macros;
    /// * the accessors below read the engine the model is bound to.
    class ModelBase
    {
    public:
        //////  Constructor  //////

        virtual ~ModelBase() = default;

        /// Bind the engine `mcmc` the model reads.
        void bind_solver(MCMCBase *mcmc) noexcept { mcmc_ = mcmc; }

        //////  Model-declared keys  //////

        /// Names of the keys per slot, overridden by `HOUSEPI_DECLARE_<SLOT>_KEYS`, empty otherwise.
        virtual std::vector<std::string> data_keys() const { return {}; }
        virtual std::vector<std::string> house_member_keys() const { return {}; }
        virtual std::vector<std::string> house_member_vec_keys() const { return {}; }
        virtual std::vector<std::string> house_double_vec_keys() const { return {}; }
        virtual std::vector<std::string> house_double_mtx_keys() const { return {}; }
        virtual std::vector<std::string> long_keys() const { return {}; }
        virtual std::vector<std::string> param_keys() const { return {}; }
        virtual std::vector<std::string> group_rate_keys() const { return {}; }
        virtual std::vector<std::string> latent_keys() const { return {}; }

        //////  Model  //////

        /// Whether the model uses contact rates.
        virtual bool has_contact_rate() const { return false; }
        /// Whether the model requires an infectivity profile definition.
        virtual bool has_infec_profile() const { return false; }
        /// Whether households enter the study through a case.
        virtual bool has_case_ascertainment() const { return false; }

        /// Recompute the cached quantities of `house` read by the transmission model.
        virtual void update_house_transmission_cache(const Household & /*house*/) {}
        /// Recompute the cached quantities of member `indiv` read by the transmission model.
        virtual void update_member_transmission_cache(const Household & /*house*/,
                                                      const Individual & /*indiv*/) {}
        /// Compute the transmission log-likelihood of `house`.
        virtual double compute_house_transm_log_lik(const Household & /*house*/) = 0;
        /// Export named quantities of `house` at the current state into `out`.
        /// * called on an equivalence-class representative, once every class is evaluated;
        /// * default -> nothing.
        /// @note A cache written by compute_house_transm_log_lik() may hold a rejected proposal.
        virtual void export_house_quantities(const Household & /*house*/, HouseQuantities & /*out*/) {}

        /// Initialize the model for evaluation (cache variables, lookup tables).
        virtual void init_model_for_evaluation() {}
        /// Initialize household `house` for evaluation.
        /// * runs right after init_model_for_evaluation().
        virtual void init_house_for_evaluation(const Household & /*house*/) {}

        /// Initialize the model for simulation (cache variables, lookup tables).
        virtual void init_model_for_simulation() {}
        /// Simulate household `house`.
        virtual void simulate_single_house(Household & /*house*/) = 0;

        //////  Setter: model requests  //////

        /// Declare the model definition requirements.
        virtual void init_model_requirements() {}

        /// Request that baseline parameter `p` and group rates `rates` share one reference.
        /// @param contacts Whether contact rates scale the same term as `rates`.
        void request_rate_identifiability(size_t p, std::vector<size_t> rates, bool contacts = false)
        {
            mcmc_->rates.request_identifiability(p, std::move(rates), contacts);
        }

        /// Declare the group rates `vec` a rule may stratify by infected status.
        /// @note Leave out rates read before a member's status is drawn in simulation.
        void request_rate_status_stratified(std::vector<size_t> vec)
        {
            mcmc_->rates.request_status_stratified(std::move(vec));
        }

        /// Request discrete parameters `vec`.
        void request_param_discrete(std::vector<size_t> vec) { mcmc_->param.request_discrete(vec); }
        /// Request positive parameters `vec`.
        void request_param_positive(std::vector<size_t> vec) { mcmc_->param.request_positive(vec); }
        /// Request parameters `vec` in one block random walk, joining any block holding one of them.
        /// * sampled with the default proposal -> moves in the block;
        /// * pinned, or a declared proposal -> left out.
        void request_param_rw_block(std::vector<size_t> vec) { mcmc_->param.request_rw_block(std::move(vec)); }

        /// Request discrete latents `vec`.
        void request_latent_discrete(std::vector<size_t> vec) { mcmc_->augm.request_discrete(vec); }
        /// Request positive latents `vec`.
        void request_latent_positive(std::vector<size_t> vec) { mcmc_->augm.request_positive(vec); }
        /// Declare latents `vec` as generated and read by transmission, cache included.
        /// A latent listed here takes no model: transmission gives its density.
        void request_latent_transmission(std::vector<size_t> vec) { mcmc_->augm.request_transmission(vec); }

        /// Request a discrete infectivity profile.
        void request_infec_profile_discrete() { mcmc_->infec_profile.request_discrete(); }

        /// Request discrete data columns `vec`.
        void request_data_discrete(std::vector<size_t> vec) { mcmc_->data.dataset().request_discrete(vec); }
        /// Request positive data columns `vec`.
        void request_data_positive(std::vector<size_t> vec) { mcmc_->data.dataset().request_positive(vec); }
        /// Request data columns `vec` with a value for every member, escaped ones included.
        void request_data_every_member(std::vector<size_t> vec) { mcmc_->data.dataset().request_everyone(vec); }

        /// Throw if some member has missing status.
        void assert_no_missing_status() const { mcmc_->data.dataset().assert_no_missing_status(); }

        //////  Setter: model calls  //////

        /// Set member `indiv` infected, at the status the status process draws.
        /// * no status process -> the lowest infected status.
        /// @note Call init_missing_latent() after to initialize individual status-specific variables.
        void set_infected(Household &house, Individual &indiv)
        {
            house.set_new_status_(indiv, mcmc_->status_process.draw(indiv, mcmc_->settings.prng()));
        }
        /// Set member `indiv` escaped.
        void set_escaped(Household &house, Individual &indiv)
        {
            house.set_new_status_(indiv, status::ESCAPED);
        }

        /// Set latent `a` of member `m` to `x`; throw if not allowed.
        void set_latent_value(size_t a, const Individual &m, double x) { mcmc_->augm.at(a, m).set(x); }

        /// Introduce infection in `house` through latent `a`.
        /// Set every other member to escaped.
        /// * member: drawn uniformly;
        /// * time: `t0`.
        /// @todo Generalize to user-choices other than uniform.
        /// @note This procedure is validated by SBC on case-ascertained.
        /// @return The member introducing the infection.
        MemberIdx set_introduction(Household &house, size_t a, double t0)
        {
            const MemberIdx i0 = uniform() * house.size();
            Individual &indiv_0 = house.at(i0);
            set_infected(house, indiv_0);
            set_latent_value(a, indiv_0, t0);
            init_missing_latent(indiv_0);

            for (Individual &indiv : house)
                if (indiv != indiv_0)
                    set_escaped(house, indiv);
            return i0;
        }

        /// Initialize the latents and observations of member `indiv`.
        /// * latent set by the user -> skipped.
        /// @todo Allow order-dependent assignments in simulation.
        void init_missing_latent(Individual &indiv)
        {
            mcmc_->augm.init_definitions_for_simulation(indiv);
            mcmc_->augm.draw_sample_for_simulation(indiv, mcmc_->settings.prng());
            mcmc_->augm.sample_obs(indiv, mcmc_->settings.prng());
        }

        //////  Getter  //////

        /// Draw a uniform variate on [0, 1].
        double uniform() const { return mcmc_->settings.uniform(); }
        /// Draw a unit-rate exponential variate.
        double exponential() const { return mcmc_->settings.exponential(); }

        /// User-configured time step of the simulation.
        double real_time_step() const { return mcmc_->settings.time_step; }
        /// Number of thread slots, the size of a per-thread pool.
        size_t n_threads() const { return mcmc_->settings.n_threads; }
        /// Slot of the calling thread in a per-thread pool.
        size_t thread_slot() const { return mcmc_->settings.thread_slot(); }

        /// Value of parameter `p`.
        double par(size_t p) const { return mcmc_->param.at(p).value(); }
        /// Latent `a` of member `m`.
        double aug(size_t a, const Individual &m) const { return mcmc_->aug(a, m); }
        /// Datum `d` of member `m`.
        double data(size_t d, const Individual &m) const { return mcmc_->data.at(d, m).value(); }
        /// Group rate `g` of member `m`.
        double group_rate(size_t g, const Individual &m) const { return mcmc_->rates.group().value(g, m); }

        /// Name of member `m`, for error messages.
        std::string str(const Individual &m) const { return mcmc_->data.dataset().str(m); }
        /// Name of household `h`, for error messages.
        std::string str(const Household &h) const { return mcmc_->data.dataset().str(h); }

        /// Contact rate from member `from` to member `to` of `house`.
        double contact_rate(const Household &house, MemberIdx from, MemberIdx to) const
        {
            return mcmc_->rates.contact().value(house, from, to);
        }
        /// Whether members `i` and `j` of `house` are equivalent according to rates.
        bool is_same_rate_partition(const Household &house, MemberIdx i, MemberIdx j) const
        {
            return mcmc_->rates.same_partition(house, i, j);
        }
        /// Contact rate within the partition of member `i` (self-contact).
        double intra_partition_contact_rate(const Household &house, MemberIdx i) const
        {
            return mcmc_->rates.contact().intra_value(house, i);
        }

        /// Instantaneous infectivity profile of member `m` at time `t`.
        double infec_profile_inst(const Individual &m, double t) const
        {
            return mcmc_->infec_profile.inst(&m, t);
        }
        /// Cumulative infectivity profile of member `m` at time `t`.
        double infec_profile_cumul(const Individual &m, double t) const
        {
            return mcmc_->infec_profile.cumul(&m, t);
        }

        //////  Cache  //////

        /// @name Model cache of the equivalence class of `house`, by cache key.
        /// * recovered on a rejected latent move;
        /// * recovered on a rejected move of a parameter feeding the infectivity profile.
        /// @note Values computed from other parameters may hold a rejected proposal: recompute them.
        ///@{

        /// Cached member value.
        MAKE_HOUSEPI_CACHE_FORWARDER(member)
        /// Cached member vector.
        MAKE_HOUSEPI_CACHE_FORWARDER(member_vector)
        /// Cached double vector.
        MAKE_HOUSEPI_CACHE_FORWARDER(double_vector)
        /// Cached double matrix.
        MAKE_HOUSEPI_CACHE_FORWARDER(double_matrix)
        ///@}

    private:
        /// Engine the model reads.
        MCMCBase *mcmc_;
    };
}

#undef MAKE_HOUSEPI_CACHE_FORWARDER
