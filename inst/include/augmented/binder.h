/* housepi - Epidemiological Analysis of Household Outbreaks.
 * Copyright 2024-2026 Institut Pasteur (Paris), France.
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#pragma once

#include <string>

#include "tools/keys.h"
#include "augmented/definition.h"
#include "augmented/scheme.h"
#include "augmented/state.h"

//////////////////////////////  augm::Binder  //////////////////////////////

namespace housepi
{
    struct Data;
    struct Translator;
    class Household;
}

namespace housepi::augm
{
    /// Build Variables according to definitions over members and households.
    ///
    /// Translator resolves the sources a definition reads.
    struct Binder
    {
        Binder(const DefTable &def,
               const KeyHolderWithRequests &keys,
               const Translator &translator,
               Data &data)
            : def_(def),
              keys_(keys),
              translator_(translator),
              data_(data) {}

        /// Bind definition row `k` to `anchor`.
        ///
        /// Where the operands come from depends on the row's level:
        /// * household row -> per-member cells, folded over the members of `house`;
        /// * any other row -> `anchor` alone.
        /// @param house Household of a household row, nullptr for any other row.
        void bind_variable(Variable &var, size_t k,
                           const Individual &anchor,
                           const Household *house) const;

        /// Initialize Variable `var` of definition row `k` for member `indiv`.
        ///
        /// Binds, in order:
        /// * value or function;
        /// * model;
        /// * proposal, in inference only.
        ///
        /// In simulation:
        /// * a generated observation is a placeholder;
        /// * the value of a sampled latent pins the Variable to it.
        /// @note In simulation, never call it on a bare latent with no value: transmission sets it.
        void bind_member(Variable &var, size_t k,
                         const Individual &indiv,
                         bool simulation) const;

        /// Build Variables of deterministic observation `dobs` over `house`.
        ///
        /// Each Variable is wired to the data cells it stands in for:
        /// * household level -> one Variable, reduced over members;
        /// * member level -> one Variable per member.
        ///
        /// With `for_evaluation`, each Variable is checked against column requests, then fixed.
        /// @note Every deterministic observation must be registered as a data key.
        void bind_deterministic_obs(State &state,
                                    const DeterministicObs &dobs,
                                    const Household &house,
                                    bool for_evaluation) const;

    private:
        const DefTable &def_;               /// Definitions, one row per latent or observation.
        const KeyHolderWithRequests &keys_; /// Latent keys and their model requests.
        const Translator &translator_;      /// Resolves definition entries into Attributes.
        Data &data_;                        /// Data columns, whose read-cells deterministic observations fill.

        /// Assign model, arguments and domain to definition row `k`.
        void bind_distr_(Variable &var, size_t k,
                         const Individual &indiv) const;

        /// Assert that `var` meets model requests of data column `d`.
        /// * NA `var` -> not checked.
        /// @param context Entity named in the error message.
        void assert_datum_requests_(size_t d,
                                    const Variable &var,
                                    const std::string &context) const;
    };
}
