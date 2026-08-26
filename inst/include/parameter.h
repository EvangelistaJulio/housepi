/* housepi - Epidemiological Analysis of Household Outbreaks.
 *
 * Copyright © 2024-2026 Institut Pasteur (Paris), CNRS and INSERM.
 * Distributed under the terms of the GNU General Public License (GPLv3).
 * See the LICENSE file for details.
 */

#pragma once

#include <string>
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
    struct ContactRisk;
    struct GroupRisk;
    struct Translator;

    /// Implementation of Parameters.
    struct Parameter : public KeyModuleWithRequests
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

        /// Request that a baseline parameter and the relative risks scaling it share one reference.
        void request_risk_identifiability(size_t p, std::vector<size_t> risks, bool contacts);

        /// Assign parameter definitions to dictionary.
        /// Expect unique definitions.
        void build(const str_table &deftable, const Translator &translator);
        /// Build parameter dictionaries for expected measure on posterior sample.
        /// Expect a non-empty table.
        void build_for_predic(const dbl_table &deftable);
        /// Mark parameter as a producer of transmission cache.
        void mark_produces_transm_cache(size_t p);

        /// Assert that at most one scale is free among a baseline and the risks scaling it.
        /// Pin the baseline to the reference value when it is the only way to resolve.
        void assert_risk_identifiability(const GroupRisk &group_risk, const ContactRisk &contact_risk);

        /// Randomly initialize values of the parameters to be inferred.
        void restart(rnd_gen &prng);
        /// Finalize all variables (see Variable::finalize).
        void finalize();
        /// Load parameters values of the n-th posterior sample.
        void load_sample_values(size_t n);

        /// Baseline parameter and the relative risks scaling it, as declared by the model.
        struct RiskAnchor
        {
            size_t baseline;           /// parameter scaled by the risks
            std::vector<size_t> risks; /// (relative) risk keys
            bool contacts;             /// whether contact risks scale the same term
        };

        //////////////////////
        //////  GETTER  //////
        //////////////////////

        /// Baselines declared by the model, with the risks anchored to them.
        const std::vector<RiskAnchor> &risk_anchors() const noexcept { return risk_anchors_; }

        /// Get reference to parameter Variable class.
        Variable &at(size_t p) noexcept { return param_vec_[p]; }
        /// Get reference to parameter constant Variable class.
        const Variable &at(size_t p) const noexcept { return param_vec_[p]; }

        /// Whether parameter storage was allocated.
        bool is_allocated() const noexcept { return !param_vec_.empty(); }

        /// Flag to check if there are parameters to sample.
        bool is_any_sample() const noexcept { return !list_sample.empty(); }

        /// Total log-prior over sampled parameters.
        double total_log_prior() const
        {
            double out = 0.;
            for (size_t p : list_sample)
                out += at(p).log_pdf();
            return out;
        }

    private:
        std::vector<Variable> param_vec_ = {};      /// parameters vector
        std::vector<RiskAnchor> risk_anchors_ = {}; /// model requested identifiability

        /// Initialize parameter vector.
        void init_param_vec_();
        /// Turn a parameter into the constant reference value, dropping it from the sample.
        void helper_pin_to_reference_(size_t p);
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
