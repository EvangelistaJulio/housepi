/* housepi - Epidemiological Analysis of Household Outbreaks.
 *
 * Copyright © 2024-2026 Institut Pasteur (Paris), CNRS and INSERM.
 * Distributed under the terms of the GNU General Public License (GPLv3).
 * See the LICENSE file for details.
 */

#pragma once

#include <string>
#include <unordered_set>
#include <vector>

#include "tools/keys.h"
#include "tools/parser.h" // UPDATE: this shouldn't be public
#include "tools/random.h"
#include "tools/types.h"

#include "function.h"
#include "variable.h"

/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////  *** Parameter : struct ***  ////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

namespace housepi
{
    class Individual;

    /// Implementation of Parameters.
    struct Parameter : public KeyModule
    {
        /// NOTE:
        /// list_const + list_sample = all parameters in param_vec_
        /// list_func < list_const (function-based parameters are deterministic)
        /// list_sample_standard + list_sample_cached = list_sample
        /// list_adaptable_sample_standard + list_non_adapt_sample_standard = list_sample_standard
        /// list_adaptable_sample_cached + list_non_adapt_sample_cached = list_sample_cached

        std::vector<size_t> list_const = {};  /// list of constant parameters
        std::vector<size_t> list_func = {};   /// list of function-based parameters
        std::vector<size_t> list_sample = {}; /// list of parameters to sample

        std::vector<size_t> list_sample_standard = {};
        std::vector<size_t> list_sample_cached = {};

        std::vector<size_t> list_adaptable_sample_standard = {};
        std::vector<size_t> list_non_adapt_sample_standard = {};
        std::vector<size_t> list_adaptable_sample_cached = {};
        std::vector<size_t> list_non_adapt_sample_cached = {};

        std::vector<std::vector<double>> simulation_values = {}; /// simulation values

        //////////////////////
        //////  SETTER  //////
        //////////////////////

        /// Assign parameter definitions to dictionary.
        /// Expect unique definitions.
        void build(const str_table &deftable);
        /// Build parameter dictionaries for expected measure on posterior sample.
        /// Expect a non-empty table.
        void build_for_predic(const dbl_table &deftable);
        /// Adjust lists of parameters to be inferred making parameter as dependent on cache.
        void make_cache_dependent(size_t p);
        /// Adjust lists of parameters to be inferred making parameter as dependent on cache.
        void make_cache_dependent(const parse::Entry &node);

        /// Request discrete parameter definition.
        void request_discrete(size_t p) { user_discrete_set_.emplace(p); }

        /// Randomly initialize values of the parameters to be inferred.
        void user_reset(rnd_gen &prng);
        /// Load parameters values of the n-th posterior sample.
        void load_sample_values(size_t n);

        /// Propose a valid move (return true), otherwise update refusal (false).
        bool propose_valid_move(size_t p, rnd_gen &prng);
        /// Adpat the proposal of model parameters.
        void adapt_proposal_of_sampled(unsigned long step, unsigned long tot);

        //////////////////////
        //////  GETTER  //////
        //////////////////////

        /// Get reference to parameter Variable class.
        Variable &at(size_t p) { return param_vec_[p]; }
        /// Get reference to parameter constant Variable class.
        const Variable &at(size_t p) const { return param_vec_[p]; }
        /// Get reference to parameter constant Variable class.

        /// Flag to check if there are parameters to sample.
        bool is_any_sample() const noexcept { return !list_sample.empty(); }

        ///////////////////////
        //////  HELPERS  //////
        ///////////////////////
        /// UPDATE: these are just temporary fixes.

        /// Links used just by read_attribute.
        void bind_to_data_wide(const KeyHolder *khold) { pDataWideKeys_ = khold; }
        void bind_to_group_rate(const KeyHolder *khold) { pGroupRateKeys_ = khold; }
        void bind_to_augm(const KeyHolder *khold) { pAugmKeys_ = khold; }

        /// Auxiliary method to read a definition as a single attribute.
        const Attribute &read_attribute(const parse::Entry &node, const Individual *member_ptr = nullptr) const;
        /// Auxiliary method to convert a typed string into a numeric attribute.
        const Attribute &read_numeric(const parse::Entry &node, const Individual *member_ptr = nullptr) const;

        /// Auxiliary method to read a definition as a single attribute (overload for string).
        const Attribute &read_attribute(const std::string &expr, const Individual *member_ptr = nullptr) const
        {
            return read_attribute(parse::Definition(expr).root, member_ptr);
        }
        /// Auxiliary method to read a definition as a single attribute (overload for string).
        const Attribute &read_numeric(const std::string &expr, const Individual *member_ptr = nullptr) const
        {
            return read_numeric(parse::Definition(expr).root, member_ptr);
        }

    private:
        std::vector<Variable> param_vec_ = {};              /// parameters vector
        std::unordered_set<size_t> user_discrete_set_ = {}; /// user defined discrete parameters

        /// Initialize parameter vector.
        void init_param_vec_();
        /// Drop current lists of parameters.
        void clear_lists_();
        /// Helper to build constant parameter.
        void helper_build_constant_(const str_table &deftable, size_t k);
        /// Helper to build parameter to sample.
        void helper_build_to_sample_(const str_table &deftable, size_t k);
        /// Helper to build function-based parameter.
        void helper_build_function_(const str_table &deftable, size_t k);

        /// UPDATE: waiting for an updated read_attribute framework...
        mutable Attribute storage_;       /// auxiliary storage attribute used for reading.
        const KeyHolder *pDataWideKeys_;  /// pointer to external registry
        const KeyHolder *pGroupRateKeys_; /// pointer to external registry
        const KeyHolder *pAugmKeys_;      /// pointer to external registry

        void helper_store_read_function_(const parse::Entry &node, const Individual *member_ptr) const;
    };
}
