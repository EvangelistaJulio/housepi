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

#include "variable.h"

//////////////////////////////  Parameter  //////////////////////////////

namespace housepi
{
    struct Translator;

    /// Model parameters of a run.
    ///
    /// - `list_const` + `list_sample`: every parameter;
    /// - `list_func`: subset of `list_const`, deterministic function parameters;
    /// - `list_sample_standard` + `list_sample_cached`: `list_sample`.
    struct Parameter : public KeyModuleWithRequests
    {
        std::vector<size_t> list_const = {};    /// Constant parameters.
        std::vector<size_t> list_func = {};     /// Function parameters.
        std::vector<size_t> list_sample = {};   /// Parameters to sample.
        std::vector<size_t> list_adaptive = {}; /// Parameters to sample with an adaptive proposal.

        std::vector<size_t> list_sample_standard = {}; /// Sampled parameters feeding no transmission cache.
        std::vector<size_t> list_sample_cached = {};   /// Sampled parameters feeding the transmission cache.

        std::vector<std::vector<double>> simulation_values = {}; /// Simulation values, per parameter.

        //////  Setter  //////

        /// Build parameters from definition table `deftable`.
        /// @note Expects one definition per parameter.
        void build(const str_table &deftable, const Translator &translator);
        /// Pin the parameters `deftable` holds a column for, for measures on a posterior sample.
        /// @note Expects a non-empty table.
        void build_for_predic(const dbl_table &deftable);
        /// Mark parameter `p` as a producer of the transmission cache.
        void mark_produces_transm_cache(size_t p);

        /// Initialize parameters to sample, for inference.
        void restart(rnd_gen &prng);
        /// Draw parameters to sample from their prior, for simulation.
        /// Unlike restart(), it uses no proposal.
        void draw_sample(rnd_gen &prng)
        {
            for (size_t p : list_sample)
                at(p).draw_value(prng);
        }
        /// Finalize every parameter (see Variable::finalize).
        void finalize();
        /// Load parameter values of posterior sample `n`.
        void load_sample_values(size_t n);

        //////  Getter  //////

        /// Parameter `p`.
        Variable &at(size_t p) noexcept { return param_vec_[p]; }
        const Variable &at(size_t p) const noexcept { return param_vec_[p]; }

        /// Whether parameter storage is allocated.
        bool is_allocated() const noexcept { return !param_vec_.empty(); }

        /// Whether parameter `p` is pinned and cannot move.
        /// - constant by definition;
        /// - held by the simulation table of the current run.
        bool is_pinned(size_t p) const noexcept
        {
            return at(p).is_constant() || !simulation_values[p].empty();
        }

        /// Whether some parameter is sampled.
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
        std::vector<Variable> param_vec_ = {}; /// Parameters, by key.

        /// Allocate one named parameter per key.
        void init_param_vec_();
        /// Drop parameter `p` from every sampling list and register it as constant.
        void switch_to_const_(size_t p);
        /// Drop every parameter list.
        void clear_lists_();
        /// Build constant parameter of row `k`.
        void helper_build_constant_(const str_table &deftable, size_t k, const Translator &translator);
        /// Build parameter to sample of row `k`.
        void helper_build_to_sample_(const str_table &deftable, size_t k, const Translator &translator);
        /// Build function parameter of row `k`.
        void helper_build_function_(const str_table &deftable, size_t k, const Translator &translator);
    };
}
