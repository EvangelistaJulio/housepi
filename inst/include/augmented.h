/* housepi - Epidemiological Analysis of Household Outbreaks.
 *
 * Copyright © 2024-2026 Institut Pasteur (Paris), CNRS and INSERM.
 * Distributed under the terms of the GNU General Public License (GPLv3).
 * See the LICENSE file for details.
 */

#pragma once

#include <string>
#include <unordered_set>
#include <unordered_map>
#include <utility>
#include <vector>

#include "tools/keys.h"
#include "tools/logdensity.h"
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
    class Individual;
    class MCMCBase;

    struct Augmented : public KeyModule
    {
        std::vector<size_t> list_export = {};     /// augmented keys to export
        std::vector<size_t> list_obs_export = {}; /// observation keys to export

        /// Auxiliary block for individual augmented objects.
        /// UPDATE: this should become private.
        struct IndivBlock
        {
            size_t idx = 0;                 /// index of status-specific definition
            std::vector<Variable> vec = {}; /// associated variables
            LogDensity log_prior;           /// log-prior
            LogDensity log_obs;             /// observation log-likelihood
        };
        std::vector<IndivBlock> varia;

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
        void sample_obs(rnd_gen &prng, Individual &member);

        /// Reset augmented initial values for given individual.
        void restart(rnd_gen &prng, Individual &m);
        /// Initialize augmented initial values for all individuals.
        void restart(rnd_gen &prng);
        /// Optimize all variables (see Variable::optimize).
        void optimize();

        /// Resolve individual ids for augmented draws and return tracked pointers.
        std::vector<const Individual *> init_individual_sample(const std::vector<std::string> &id_draws) const;

        //////////////////////
        //////  GETTER  //////
        //////////////////////

        /// Access the augmented variable of given key for a member.
        Variable &at(size_t a, const Individual &member) noexcept { return varia[member.index_].vec[a]; }
        const Variable &at(size_t a, const Individual &member) const noexcept { return varia[member.index_].vec[a]; }

        const std::vector<size_t> &list_sample(const Individual &member) const;
        const std::vector<size_t> &list_adaptive(const Individual &member) const;
        const std::vector<size_t> &list_draws(const Individual &member) const;
        const std::vector<size_t> &list_priors(const Individual &member) const;
        const std::vector<size_t> &list_obs(const Individual &member) const;

        /// Compute and store log-prior associated to augmented variables.
        /// NOTE: to be called only if any prior is defined.
        double commit_augm_log_prior(Individual &member);
        /// Compute and store log-ratio of prior of single augmented variable.
        double commit_augm_log_prior_ratio(Individual &member, size_t a);

        /// Compute and store observation log-likelihood of observation nodes.
        /// NOTE: to be called only if any observation node is defined.
        double commit_log_lik_obs(Individual &member);
        /// Compute and store log-ratio of observation node.
        double commit_log_lik_obs_ratio(Individual &member, size_t a);

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
        /// If moving given parameter may change the observation log-likelihood.
        bool param_affects_obs_log_lik(size_t p) const noexcept { return obs_param_depen_[p]; }
        /// If moving augmented 'a' for this member parameterizes any obs node.
        bool affects_obs(const Individual &member, size_t a) const
        {
            return !status_blocks_[varia[member.index_].idx].obs_depen[a].empty();
        }

    private:
        MCMCBase *mcmc_;               /// pointer to owning MCMCBase // UPDATE drop
        bool is_any_sample_ = false;   /// at least one augmented is to sample
        bool is_any_adaptive_ = false; /// at least one augmented has adaptive proposal
        bool is_any_prior_ = false;    /// at least one augmented prior is provided
        bool is_any_obs_ = false;      /// at least one observation node is provided

        /// Parameters whose moves can change the augmented log-prior.
        std::vector<bool> prior_param_depen_ = {};
        /// Parameters whose moves can change the observation log-likelihood.
        std::vector<bool> obs_param_depen_ = {};

        /// User defined discrete parameters
        std::unordered_set<size_t> user_discrete_set_ = {};
        /// Status to definition row.
        std::unordered_map<int, std::vector<size_t>> map_status_to_def_row_ = {};
        /// Status to internal index.
        std::unordered_map<int, size_t> map_status_to_block_ = {};

        /// Auxiliary definition table.
        struct DefTable_
        {
            using col = std::vector<std::string>;
            col c_name = {}, c_status = {};
            col c_func = {}, c_fcumul = {};
            col c_value = {}, c_const = {}, c_discr = {};
            col c_prior = {}, c_min = {}, c_max = {};
            col c_proposal = {}, c_scale = {};
            col c_kind = {}; /// "prior" (augmented) or "obs" (observation)

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

            std::vector<size_t> priors = {};             /// keys for priors
            std::vector<size_t> sample = {};             /// keys for sample
            std::vector<size_t> adaptive = {};           /// keys for adaptive (subset of sample)
            std::vector<size_t> draws = {};              /// keys for draws
            std::vector<std::vector<size_t>> depen = {}; /// keys for respective dependecies

            std::vector<size_t> obs = {};                    /// keys for observation nodes
            std::vector<std::vector<size_t>> obs_depen = {}; /// latent key -> obs nodes conditioned on it
        };
        std::vector<StatusBlock_> status_blocks_;

        /// Check augmented names and assert no reserved names are used.
        void helper_build_assert_augm_input_(const str_table &input);

        /// Collect parameters feeding any augmented prior. Infectivity profile excluded.
        /// NOTE: call once parameters are registered, i.e. from init_for_evaluation.
        /// NOTE: prior parameters may be registered after the augmented definitions.
        void set_prior_param_depen_();

        /// Eventually append infectivity profile definitions.
        void helper_build_append_infec_profile_(InfecProfile &infec_profile);

        /// Assert augmented status/generic definitions are non degenerate.
        /// Assert status-specific definitions are only for infected status.
        /// NOTE: generic definitions are associated to MISSING status.
        void helper_build_definition_by_status_();

        /// Deal with the definition 'k' of augmented to be sampled (const = FALSE).
        /// NOTE: this assume definitions are univocal and status correct.
        void helper_build_set_allocation_order_();

        /// Set dependecies ordering for ininitializiation fo given status.
        /// At the same time assert that definition are complet and acyclic.
        void set_definition_dependency_order_for_(int status, int err_status);
        /// Extract status-specific definitions row indexes.
        std::vector<size_t> get_unsorted_definition_rows_(int status) const;
        /// Set sample and draws list for given status.
        void set_sample_list_(int status);
        /// Set sample and draws list for given definition.
        void set_sample_list_for_single_defininition_(int status, size_t k);

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
