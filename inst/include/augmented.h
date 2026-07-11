/* housepi - Epidemiological Analysis of Household Outbreaks.
 *
 * Copyright © 2024-2026 Institut Pasteur (Paris), CNRS and INSERM.
 * Distributed under the terms of the GNU General Public License (GPLv3).
 * See the LICENSE file for details.
 */

#pragma once

#include <cstdint>
#include <string>
#include <unordered_set>
#include <unordered_map>
#include <utility>
#include <vector>

#include "tools/keys.h"
#include "tools/logproba.h"
#include "tools/nan.h"
#include "tools/random.h"
#include "tools/types.h"

#include "individual.h"
#include "variable.h"

/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////  *** Augmented : struct ***  /////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

namespace housepi
{
    struct InfecProfile;
    struct DataWide;
    class Household;
    class MCMCBase;

    namespace parse
    {
        struct Entry;
    }

    /// Max resamples per member (see Augmented::restart_obs_aware).
    static constexpr size_t N_INIT_ATTEMPTS_ = 100;

    struct Augmented : public KeyModule
    {
        /// Discriminates definition rows: augmented latent vs observation flavours.
        enum class Kind : std::uint8_t
        {
            Augm,
            Obs,
            House,
            Status
        };

        std::vector<size_t> list_export = {};       /// augmented to export
        std::vector<size_t> list_obs_export = {};   /// observation to export
        std::vector<size_t> list_house_export = {}; /// household-observation to export

        /// Auxiliary block for individual augmented objects.
        struct IndivBlock
        {
            size_t idx = 0;                 /// index of status-specific definition
            std::vector<Variable> vec = {}; /// associated variables
            LogProba log_prior;             /// augmented log-prior
            LogProba log_obs;               /// observation log-likelihood
        };
        std::vector<IndivBlock> varia;

        /// Household-level observations, per house per house-export slot. The single
        /// source of truth for a derived column: standing-in data cells point at it.
        /// NOTE: built for both inference (data-bound) and simulation (live).
        std::vector<std::vector<Variable>> house_varia;

        ///////////////////////////
        //////  CONSTRUCTOR  //////
        ///////////////////////////

        Augmented(MCMCBase *mcmc);
        ~Augmented();

        //////////////////////
        //////  SETTER  //////
        //////////////////////

        /// Assign augmented variables to all individuals.
        void build(InfecProfile &infec_profile, const str_table &input);

        /// Set all augmented definitions to empty for all individuals.
        void clear_all();
        /// Set individual augmented definitions to empty. Used in simulations.
        void clear_definitions(Individual &m);
        /// Initialize augmented variables only if missing.
        void init_definitions_for_simulation(Individual &m);
        /// Initialize augmented dependecies and priors.
        /// NOTE: evaluation-only, never run for simulations.
        void init_for_evaluation();

        /// Request discrete augmented definition.
        void request_discrete(size_t a) { user_discrete_set_.emplace(a); }

        /// Sample observation values from their distributions.
        /// NOTE: used for simulations.
        void sample_obs(Individual &member, rnd_gen &prng);

        /// Sample the status observation, if any.
        /// NOTE: used for simulations.
        void simulate_infected_status(Household &house, Individual &member, rnd_gen &prng);

        /// Initialize household-level observations for simulations.
        void init_house_for_simulation();
        /// Clear every member's observation Variables of a household before re-sampling.
        /// FIXME: I don't think I need this anymore
        void reset_obs_data_for_simulation(Household &house);

        /// Build the household observations that stand in for required data columns and
        /// point each member's read-surface at its Variable before evaluation.
        /// NOTE: only pure-data (constant) household observations are eligible.
        void materialize_data_bound_house_obs(DataWide &data_wide);

        /// Reset augmented initial values for given individual.
        void restart(Individual &member, rnd_gen &prng)
        {
            for (size_t a : list_sample(member))
                varia_of_(member).vec[a].restart(prng);
        }
        /// Initialize augmented initial values for all individuals.
        /// NOTE: skips automatically if no augmented to sample.
        void restart(rnd_gen &prng);
        /// Optimize all variables (see Variable::optimize).
        void optimize();

        /// Resolve individual ids for augmented draws and return tracked pointers.
        std::vector<const Individual *> init_individual_sample(const std::vector<std::string> &id_draws) const;

        //////////////////////
        //////  GETTER  //////
        //////////////////////

        /// Access the augmented variable of given key for a member.
        Variable &at(size_t a, const Individual &m) noexcept { return varia_of_(m).vec[a]; }
        const Variable &at(size_t a, const Individual &m) const noexcept { return varia_of_(m).vec[a]; }

        const std::vector<size_t> &list_sample(const Individual &m) const { return block_of_(m).sample; }
        const std::vector<size_t> &list_adaptive(const Individual &m) const { return block_of_(m).adaptive; }
        const std::vector<size_t> &list_draws(const Individual &m) const { return block_of_(m).draws; }
        const std::vector<size_t> &list_priors(const Individual &m) const { return block_of_(m).priors; }
        const std::vector<size_t> &list_obs(const Individual &m) const { return block_of_(m).obs; }

        /// Export column name of the household-level observation.
        const std::string &house_obs_name(size_t b) const { return deftable_.c_name[list_house_export[b]]; }

        /// Compute and store log-prior associated to augmented variables.
        /// NOTE: to be called only if any prior is defined.
        double commit_augm_log_prior(Individual &member);
        /// Compute and store log-ratio of prior of single augmented variable.
        double commit_augm_log_prior_ratio(Individual &member, size_t a);

        /// Compute and store observation log-likelihood of observation nodes.
        /// NOTE: to be called only if any observation node is defined.
        double commit_obs_log_lik(Individual &member);
        /// Compute and store log-ratio of observation node.
        double commit_obs_log_lik_ratio(Individual &member, size_t a);

        /// Read-only individual observation log-likelihood.
        double eval_log_obs(const Individual &member) const
        {
            double ll = 0.0;
            for (size_t a : list_obs(member))
                ll += varia_of_(member).vec[a].log_obs();
            return ll;
        }

        /// Throw initialization error for associated to observation model.
        [[noreturn]] void throw_obs_init_error(const Individual &member) const;

        /////////////////////
        //////  FLAGS  //////
        /////////////////////

        /// If there are non-constant augmented definitions.
        bool is_any_sample() const noexcept { return is_any_sample_; }
        /// If there are augmented definitions with an adaptive proposal.
        bool is_any_adaptive() const noexcept { return is_any_adaptive_; }
        /// If there are augmented definitions with a prior.
        bool is_any_prior() const noexcept { return is_any_prior_; }
        /// If moving given parameter may change the augmented log-prior.
        bool param_affects_augm_prior(size_t p) const noexcept { return prior_param_depen_[p]; }
        /// If there are definitions for observations.
        bool is_any_obs() const noexcept { return is_any_obs_; }
        /// If a status observation is defined.
        bool is_any_status_obs() const noexcept { return !na::is(status_obs_def_); }
        /// Whether the given name is an observation node (kind "obs").
        bool is_obs_node(const std::string &name) const noexcept;
        /// If moving given parameter may change the observation log-likelihood.
        bool param_affects_obs_log_lik(size_t p) const noexcept { return obs_param_depen_[p]; }
        /// If moving augmented 'a' for this member parameterizes any obs node.
        bool affects_obs(const Individual &m, size_t a) const { return !block_of_(m).obs_depen[a].empty(); }

    private:
        MCMCBase *mcmc_;               /// pointer to owning MCMCBase // UPDATE drop
        bool is_any_sample_ = false;   /// at least one augmented is to sample
        bool is_any_adaptive_ = false; /// at least one augmented has adaptive proposal
        bool is_any_prior_ = false;    /// at least one augmented prior is provided
        bool is_any_obs_ = false;      /// at least one observation node is provided

        size_t status_obs_def_ = na::INDEX; /// deftable row of status observation
        size_t status_obs_key_ = na::INDEX; /// key-holder slot of status observation

        /// User defined discrete parameters
        std::unordered_set<size_t> user_discrete_set_ = {};
        /// Status to definition row.
        std::unordered_map<int, std::vector<size_t>> map_status_to_def_row_ = {};
        /// Status to internal index.
        std::unordered_map<int, size_t> map_status_to_block_ = {};

        /// Parameters whose moves can change the augmented log-prior.
        std::vector<bool> prior_param_depen_ = {};
        /// Parameters whose moves can change the observation log-likelihood.
        std::vector<bool> obs_param_depen_ = {};

        /// Auxiliary definition table.
        struct DefTable_
        {
            using col = std::vector<std::string>;
            col c_name = {}, c_status = {};
            col c_func = {}, c_fcumul = {};
            col c_value = {}, c_const = {}, c_discr = {};
            col c_prior = {}, c_min = {}, c_max = {};
            col c_proposal = {};
            std::vector<Kind> c_kind = {}; /// augm | obs | house | status

            /// Clear all definitions.
            void clear();
            /// Append definitions to current table.
            void append(const str_table &input);
            /// Columns presenting definitions possibly containing sources.
            std::vector<const col *> source_cols() const
            {
                return {&c_value, &c_prior, &c_func, &c_min, &c_max};
            }
        } deftable_;

        /// Auxiliary block for status-specific priors.
        struct StatusBlock_
        {
            StatusBlock_(int status_) { status = status_; }

            int status; /// associated status

            std::vector<size_t> order = {}; /// definition order as in the deftable

            std::vector<size_t> priors = {};   /// keys for priors
            std::vector<size_t> sample = {};   /// keys for sample
            std::vector<size_t> adaptive = {}; /// keys for adaptive (subset of sample)
            std::vector<size_t> draws = {};    /// keys for draws
            /// key -> key and its transitive consumers
            std::vector<std::vector<size_t>> consumers = {};

            std::vector<size_t> obs = {};                    /// keys for observation nodes
            std::vector<std::vector<size_t>> obs_depen = {}; /// latent key -> obs nodes conditioned on it
        };
        std::vector<StatusBlock_> status_blocks_;

        /// Access IndivBlock for a member.
        IndivBlock &varia_of_(const Individual &m) noexcept { return varia[Index::of(m)]; }
        const IndivBlock &varia_of_(const Individual &m) const noexcept { return varia[Index::of(m)]; }
        /// Status-specific block governing augmented definitions of member.
        StatusBlock_ &block_at_(int status) { return status_blocks_[map_status_to_block_.at(status)]; }
        const StatusBlock_ &block_of_(const Individual &m) const { return status_blocks_[varia_of_(m).idx]; }

        /// Check augmented names and assert no reserved names are used.
        void helper_build_assert_input_(const str_table &input);

        /// Collect parameters feeding any augmented prior. Infectivity profile excluded.
        /// NOTE: call once parameters are registered, i.e. from init_for_evaluation.
        /// NOTE: prior parameters may be registered after the augmented definitions.
        void set_prior_param_depen_();

        /// Eventually append infectivity profile definitions.
        void helper_build_append_(InfecProfile &infec_profile);

        /// Assert augmented status/generic definitions are non degenerate.
        /// Assert status-specific definitions are only for infected status.
        /// NOTE: generic definitions are associated to MISSING status.
        void helper_build_definition_by_status_();

        /// Deal with the definition 'k' of augmented to be sampled (const = FALSE).
        /// NOTE: this assume definitions are univocal and status correct.
        void helper_build_set_allocation_order_();

        /// Build the augmented producer graph (name -> the @A names it reads across
        /// all source columns) for the given definition rows.
        str_table build_augm_graph_(const std::vector<size_t> &rows) const;

        /// Set dependecies ordering for ininitializiation fo given status.
        /// At the same time assert that definition are complet and acyclic.
        void set_definition_dependency_order_for_(int status, int err_status);
        /// Extract status-specific definitions row indexes.
        std::vector<size_t> get_unsorted_definition_rows_(int status) const;
        /// Set sample and draws list for given status.
        void set_sample_list_(int status);
        /// Set sample and draws list for given definition.
        void set_sample_list_for_single_definition_(int status, size_t k);

        /// Inspect if constant definition.
        std::pair<bool, bool> set_sample_list_for_function_(int status, size_t k);

        /// Deal with the definition of augmented variable with a prior.
        /// NOTE: process priors only for posterior sampling, not simulations.
        void helper_init_prior_list_();
        /// Return dictionary of ordered priors.
        str_table get_ordered_priors_(const StatusBlock_ &block) const;

        /// Deal with observation definitions.
        /// Build per-status obs lists and augmented inverse-dependency map.
        void helper_init_observation_list_();

        /// Assert status observation as correct dependency and output.
        void assert_status_obs_consistent_();

        /// Collect household-level deftable rows into list_house_export.
        void helper_init_house_obs_list_();
        /// Deal with household-level observations.
        void helper_init_house_attributes_(Household &house, size_t b);
        /// Build the (reduction) variable for household observation 'k' over 'house'.
        void build_house_obs_variable_(Variable &hvar, Household &house, size_t k);
        /// Point each standing-in data cell of a house at its household-obs Variable.
        void wire_house_obs_cells_(Household &house);
        /// Whether household observation 'k' references only data and numeric sources.
        bool house_obs_is_pure_data_(size_t k) const;
        /// Whether a parsed definition subtree references any group rate.
        static bool entry_has_group_rate_(const parse::Entry &node);

        /// Whether household observation 'k' references a generated observation column.
        bool house_obs_is_live_(size_t k) const;
        /// Reject an augmented variable that reads a live data-bound household observation.
        void assert_no_live_house_obs_cycle_() const;
        /// Reject a definition whose prior reads an observation column.
        void assert_no_obs_entry_collision_() const;
        /// Reject a group/contact rate keyed on an observation column: it feeds the
        /// transmission model the observation is generated from. Simulation only.
        void assert_no_obs_rate_collision_() const;

        /// Get the definition order for current member according to its status.
        /// NOTE: if status hasn't been observed in data, build on demand.
        std::vector<size_t> helper_materialize_definition_order_(Individual &member);
        /// Assign prior arguments and range. Same for function, otherwise assign proposer.
        /// NOTE: a Variable must already be allocated.
        /// NOTE: observations are generated by simulations.
        void helper_init_attributes_(Individual &member, size_t k, bool simulation = false);
        /// Auxiliary method to assign attributes to augmented prior.
        void helper_init_prior_attributes_(Individual &member, size_t k);
        /// Auxiliary method to assign attributes to augmented function.
        void helper_init_func_attributes_(Individual &member, size_t k);
    };
}
