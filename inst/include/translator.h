/* housepi - Epidemiological Analysis of Household Outbreaks.
 *
 * Copyright © 2024-2026 Institut Pasteur (Paris), CNRS and INSERM.
 * Distributed under the terms of the GNU General Public License (GPLv3).
 * See the LICENSE file for details.
 */

#pragma once

#include <string>

#include "function.h"

/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////  *** Translator : struct ***  ////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

namespace housepi
{
    namespace parse
    {
        struct Entry;
    }
    struct KeyHolder;
    struct Augmented;
    struct Parameter;
    class Individual;

    /// Cross-module translator of parsed definition into Attribute.
    struct Translator
    {
        Translator(Parameter &param, Augmented &augm) : param_(param), augm_(augm) {}

        /// Bind key holders for member-specific source resolution.
        void bind_data_wide(const KeyHolder *khold) noexcept { pDataWideKeys_ = khold; }
        void bind_group_rate(const KeyHolder *khold) noexcept { pGroupRateKeys_ = khold; }

        /// Translate a definition node or a string to a generic attribute.
        const Attribute &to_attribute(const parse::Entry &node, const Individual *member_ptr = nullptr) const;
        const Attribute &to_attribute(const std::string &expr, const Individual *member_ptr = nullptr) const;

        /// Translate a definition node or a string to a numeric constant attribute.
        const Attribute &to_numeric(const parse::Entry &node, const Individual *member_ptr = nullptr) const;
        const Attribute &to_numeric(const std::string &expr, const Individual *member_ptr = nullptr) const;

        /// Mark a referenced non-constant parameter as cache-dependent.
        void make_cache_dependent(const parse::Entry &node) const;

    private:
        Parameter &param_;                /// parameter key-module
        Augmented &augm_;                 /// augmented key-module
        const KeyHolder *pDataWideKeys_;  /// data-wide key-holder
        const KeyHolder *pGroupRateKeys_; /// group-rates key-holder
        mutable Attribute storage_;       /// auxiliary storage attribute

        void helper_store_read_function_(const parse::Entry &node, const Individual *member_ptr) const;
    };
}
