/* housepi - Epidemiological Analysis of Household Outbreaks.
 * Copyright 2024-2026 Institut Pasteur (Paris), France.
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#pragma once

#include <string>
#include <string_view>
#include <vector>

#include "tools/dictionary.h"
#include "tools/keys.h"
#include "tools/random.h"
#include "parameter/definition.h"
#include "parameter/scheme.h"
#include "parameter/state.h"
#include "variable.h"

//////////////////////////////  Parameter  //////////////////////////////

namespace housepi
{
    struct Translator;

    /// Model parameters of a run.
    ///
    /// * definitions, stored;
    /// * sampling scheme, built from definitions;
    /// * parameter variables, bound to definitions.
    struct Parameter : public KeyHolderWithRequests
    {
        using KeyHolderWithRequests::KeyHolderWithRequests;

        param::State state; /// Variables moved by inference and simulation.

        //////  Setter  //////

        /// Build definitions from `input`, then their sampling scheme and variables.
        /// @note Expects one definition per parameter.
        void build(const str_table &input, const Translator &translator);
        /// Pin the parameters `input` holds a column for, for measures on a posterior sample.
        /// @note Expects a non-empty table.
        void build_for_predic(const dbl_table &input);
        /// Mark parameter `p` as a producer of the transmission cache.
        void mark_produces_transm_cache(size_t p) { scheme_.mark_produces_transm_cache(p); }
        /// Request parameters `vec` in one block random walk (see Rates::list_param_rw_blocks).
        void request_rw_block(std::vector<size_t> vec) { rw_blocks_.push_back(std::move(vec)); }

        /// Initialize parameters to sample, for inference.
        void restart(RndGen &prng) { param::restart(scheme_, state, prng); }
        /// Draw parameters to sample from their prior, for simulation.
        /// Unlike restart(), it uses no proposal.
        void draw_sample(RndGen &prng) { param::draw_sample(scheme_, state, prng); }
        /// Finalize every parameter (see Variable::finalize).
        void finalize() { param::finalize(scheme_, state); }
        /// Load pinned values of simulation repetition `r`.
        void load_sample_values(size_t r) { param::load_sample_values(scheme_, state, r); }

        //////  Getter  //////

        /// Sampling scheme of parameters.
        const param::Scheme &scheme() const noexcept { return scheme_; }
        /// Parameter sets requested in one block random walk each.
        const std::vector<std::vector<size_t>> &rw_blocks() const noexcept { return rw_blocks_; }

        /// Parameter `p`.
        Variable &at(size_t p) noexcept { return state.vec[p]; }
        const Variable &at(size_t p) const noexcept { return state.vec[p]; }

        /// Whether parameter `pstr` is constant.
        bool is_constant(std::string_view pstr) const { return at(consult(pstr)).is_constant(); }

        /// Whether parameter storage is allocated.
        bool is_allocated() const noexcept { return !state.vec.empty(); }

        /// Whether parameter `p` is pinned and cannot move.
        /// * constant by definition;
        /// * held by the simulation table of the current run.
        bool is_pinned(size_t p) const noexcept
        {
            return at(p).is_constant() || !scheme_.simulation_values[p].empty();
        }

        /// Whether the prior or the bounds of parameter `p` read a parameter of `others`, other than `p`.
        bool reads_any(size_t p, const std::vector<size_t> &others) const;

        /// Whether parameter `p` may be read by rates alone.
        /// * model-declared -> no, the model reads it;
        /// * read by the function or bounds of another parameter -> no;
        /// * feeding the transmission cache -> no.
        bool is_rate_only(size_t p) const;

        /// Whether some parameter is sampled.
        bool is_any_sample() const noexcept { return !scheme_.sample.empty(); }

    private:
        param::DefTable def_;                             /// Parameter definitions.
        param::Scheme scheme_;                            /// Sampling scheme built from definitions.
        std::vector<std::vector<size_t>> rw_blocks_ = {}; /// Parameter sets the model requests in one block each.

        /// Bind Variables to definitions, producer before consumer.
        /// * each Variable checked against the requests on its key.
        void bind_(const Translator &translator);
        /// Bind constant `var` to definition `row`.
        static void bind_constant_(Variable &var, const param::DefRow &row, const Translator &translator);
        /// Bind sampled `var` to definition `row`: initial value, prior, proposal.
        static void bind_to_sample_(Variable &var, const param::DefRow &row, const Translator &translator);
        /// Bind function `var` to definition `row`.
        /// * terminal function -> alias of another parameter.
        static void bind_function_(Variable &var, const param::DefRow &row, const Translator &translator);
    };
}
