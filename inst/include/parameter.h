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
    struct Translator;

    /// Implementation of Parameters.
    struct Parameter : public KeyModule
    {
        /// NOTE:
        /// list_const + list_sample = all parameters in param_vec_
        /// list_func < list_const (function-based parameters are deterministic)
        /// list_sample = list_sample_standard + list_sample_cached

        std::vector<size_t> list_const = {};    /// constant parameters
        std::vector<size_t> list_func = {};     /// function-based parameters
        std::vector<size_t> list_sample = {};   /// parameters to sample
        std::vector<size_t> list_adaptive = {}; /// parameters to sample with adaptive proposal

        std::vector<size_t> list_sample_standard = {};
        std::vector<size_t> list_sample_cached = {};

        std::vector<std::vector<double>> simulation_values = {}; /// simulation values

        //////////////////////
        //////  SETTER  //////
        //////////////////////

        /// Assign parameter definitions to dictionary.
        /// Expect unique definitions.
        void build(const str_table &deftable, const Translator &translator);
        /// Build parameter dictionaries for expected measure on posterior sample.
        /// Expect a non-empty table.
        void build_for_predic(const dbl_table &deftable);
        /// Adjust lists of parameters to be inferred making parameter as dependent on cache.
        void make_cache_dependent(size_t p);

        /// Request discrete parameter definition.
        void request_discrete(size_t p) { user_discrete_set_.emplace(p); }

        /// Randomly initialize values of the parameters to be inferred.
        void restart(rnd_gen &prng);
        /// Optimize all variables (see Variable::optimize).
        void optimize();
        /// Load parameters values of the n-th posterior sample.
        void load_sample_values(size_t n);

        //////////////////////
        //////  GETTER  //////
        //////////////////////

        /// Get reference to parameter Variable class.
        Variable &at(size_t p) noexcept { return param_vec_[p]; }
        /// Get reference to parameter constant Variable class.
        const Variable &at(size_t p) const noexcept { return param_vec_[p]; }
        /// Get reference to parameter constant Variable class.

        /// Flag to check if there are parameters to sample.
        bool is_any_sample() const noexcept { return !list_sample.empty(); }

    private:
        std::vector<Variable> param_vec_ = {};              /// parameters vector
        std::unordered_set<size_t> user_discrete_set_ = {}; /// user defined discrete parameters

        /// Initialize parameter vector.
        void init_param_vec_();
        /// Drop current lists of parameters.
        void clear_lists_();
        /// Helper to build constant parameter.
        void helper_build_constant_(const str_table &deftable, size_t k, const Translator &translator);
        /// Helper to build parameter to sample.
        void helper_build_to_sample_(const str_table &deftable, size_t k, const Translator &translator);
        /// Helper to build function-based parameter.
        void helper_build_function_(const str_table &deftable, size_t k, const Translator &translator);
    };
}
