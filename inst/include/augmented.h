/* housepi - Epidemiological Analysis of Household Outbreaks.
 *
 * Copyright © 2026 Institut Pasteur (Paris), CNRS and INSERM.
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
#include "tools/random.h"
#include "tools/types.h"

#include "indiv_storage.h"

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
        std::vector<size_t> list_export = {}; /// list of keys to export

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
        void clear_definitions(Individual &member);
        /// Initialize augmented variables only if missing.
        void init_missing_definitions(Individual &member);
        /// Initialize augmented dependecies for evaluations.
        /// NOTE: not to be called before simulation which builds on demand.
        void init();

        /// Request discrete augmented definition.
        void request_discrete(size_t a) { user_discrete_set_.emplace(a); }

        /// Reset augmented initial values for given individual.
        void user_reset(rnd_gen &prng, Individual &member);
        /// Initialize augmented initial values for all individuals.
        void user_reset(rnd_gen &prng);

        /// Resolve individual ids for augmented draws and return tracked pointers.
        std::vector<const Individual *> init_individual_sample(const std::vector<std::string> &id_draws) const;

        //////////////////////
        //////  GETTER  //////
        //////////////////////

        const std::vector<size_t> &list_sample(const Individual &member) const;
        const std::vector<size_t> &list_draws(const Individual &member) const;
        const std::vector<size_t> &list_models(const Individual &member) const;

        /// Compute and store log-likelihood associated to augmented variable models.
        /// NOTE: to be called only if any model is defined.
        double commit_augm_model_log_lik(Individual &member);
        /// Change in log-likelihood associated to given augmented.
        double commit_augm_model_log_ratio(Individual &member, size_t a);

        /////////////////////
        //////  FLAGS  //////
        /////////////////////

        /// Whether there are non-constant augmented definitions.
        bool is_any_sample() const noexcept { return is_any_sample_; }
        /// Whether there are augmented definitions with a model.
        bool is_any_model() const noexcept { return is_any_model_; }

    private:
        MCMCBase *mcmc_;                      /// pointer to owning MCMCBase // UPDATE drop
        std::vector<Indiv::Augm> indiv_augm_; /// External storage for individual augmented data
        str_table deftable_ = {};             /// store definition table (no names)
        bool is_any_sample_ = false;          /// whether at least one augmented model is provided
        bool is_any_model_ = false;           /// whether at least one augmented model is provided

        /// User defined discrete parameters
        std::unordered_set<size_t> user_discrete_set_ = {};
        /// Status to definition row.
        std::unordered_map<int, std::vector<size_t>> map_status_to_def_row_ = {};
        /// Status to internal index.
        std::unordered_map<int, size_t> map_status_to_block_ = {};

        /// Auxiliary block for status-specific models.
        struct StatusBlock_
        {
            StatusBlock_(int status_) { status = status_; }

            int status; /// associated status

            std::vector<size_t> order = {}; /// definition order as in the deftable

            std::vector<size_t> models = {};             /// keys for models
            std::vector<size_t> sample = {};             /// keys for sample
            std::vector<size_t> draws = {};              /// keys for draws
            std::vector<std::vector<size_t>> depen = {}; /// keys for respective dependecies
        };
        std::vector<StatusBlock_> status_blocks_;

        /// Allocate external augmented storage and set paugm_ pointers on all individuals.
        void helper_build_allocate_storage_();

        /// Check augmented names and assert no reserved names are used.
        void helper_assert_augm_input_(const str_table &input);

        /// Eventually append infectivity profile definitions.
        void helper_build_append_infec_profile_(InfecProfile &infec_profile);

        /// Assert that augmented status/generic definitions are non degenerate.
        /// Also, assert that status-specific definitions are for infected status only.
        /// By convention, generic defintion are associated to missing status.
        void helper_build_definition_by_status_();

        /// Deal with the definition 'k' of augmented to be sampled (const = FALSE).
        /// NOTE: this assume definitions are univocal and status correct.
        void helper_build_allocation_order_();

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

        /// Deal with the definition of augmented variable with a model.
        /// NOTE: process models only for posterior sampling, not simulations.
        void helper_init_model_list_();
        /// Return dictionary of ordered models.
        str_table get_ordered_models_(const StatusBlock_ &block) const;

        /// Get the definition order for current member according to its status.
        /// NOTE: if status hasn't been observed in data, build on demand.
        std::vector<size_t> helper_materialize_definition_order_(Individual &member);
        /// Assign model arguments and range. Same for function, otherwise assign proposer.
        /// NOTE: a Variable must already be allocated.
        void helper_init_attributes_(Individual &member, size_t k);
        /// Auxiliary method to assign attributes to augmented model.
        void helper_init_model_attributes_(Individual &member, size_t k);
        /// Auxiliary method to assign attributes to augmented function.
        void helper_init_func_attributes_(Individual &member, size_t k);
    };
}
