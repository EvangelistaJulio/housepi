/* housepi - Epidemiological Analysis of Household Outbreaks.
 *
 * Copyright © 2024-2026 Institut Pasteur (Paris), CNRS and INSERM.
 * Distributed under the terms of the GNU General Public License (GPLv3).
 * See the LICENSE file for details.
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
    struct DataWide;
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
               const KeyModuleWithRequests &keys,
               const Translator &translator,
               DataWide &data_wide)
            : def_(def),
              keys_(keys),
              translator_(translator),
              data_wide_(data_wide) {}

        /// Bind definition row `k` to `anchor`.
        ///
        /// Where the operands come from depends on the row's level:
        /// - household row: per-member cells, folded over the members of `house`;
        /// - any other row: `anchor` alone.
        /// @param house Household of a household row, nullptr for any other row.
        void bind_variable(Variable &var, size_t k,
                           const Individual &anchor,
                           const Household *house) const;

        /// Initialize Variable `var` of definition row `k` for member `indiv`.
        ///
        /// Binds, in order:
        /// - value or function;
        /// - model;
        /// - proposal, in inference only.
        ///
        /// In simulation:
        /// - a generated observation is a placeholder;
        /// - a value paired with a model pins the Variable to that value.
        void bind_member(Variable &var, size_t k,
                         const Individual &indiv,
                         bool simulation) const;

        /// Build Variables of surrogate column `surr` over `house`.
        ///
        /// Each Variable is wired to the data cells it stands in for:
        /// - household level: one Variable, reduced over members;
        /// - member level: one Variable per member.
        ///
        /// With `for_evaluation`, each Variable is checked against column requests, then fixed.
        /// @note Every surrogate column must be registered as a data key.
        void bind_surrogate(State &state,
                            const Surrogate &surr,
                            const Household &house,
                            bool for_evaluation) const;

    private:
        const DefTable &def_;               /// Definitions, one row per latent or observation.
        const KeyModuleWithRequests &keys_; /// Latent keys and their model requests.
        const Translator &translator_;      /// Resolves definition entries into Attributes.
        DataWide &data_wide_;               /// Data columns, whose cells surrogates fill.

        /// Assign model, arguments and domain to definition row `k`.
        void bind_distr_(Variable &var, size_t k,
                         const Individual &indiv) const;

        /// Assert that `var` meets model requests of data column `d`.
        /// - NA `var`: not checked.
        /// @param context Entity named in the error message.
        void assert_datum_requests_(size_t d,
                                    const Variable &var,
                                    const std::string &context) const;
    };
}
