/* housepi - Epidemiological Analysis of Household Outbreaks.
 * Copyright 2024-2026 Institut Pasteur (Paris), France.
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#pragma once

#include <memory>
#include <string>
#include <vector>

#include "tools/dictionary.h"
#include "data/long.h"
#include "data/state.h"
#include "data/wide.h"

//////////////////////////////  Data  //////////////////////////////

namespace housepi
{
    class ModelBase;

    /// Data of a run.
    ///
    /// - imported tables, shared copy-on-write across engines;
    /// - run state, owned by one engine.
    struct Data
    {
        data::State state; /// Values written by inference and simulation.

        /// Register the data keys declared by `model`.
        explicit Data(const ModelBase &model)
            : store_(std::make_shared<Store_>(model)) {}

        //////  Setter  //////

        /// Share the imported tables of `other`, then reset the read-cells over them.
        void share(const Data &other)
        {
            store_ = other.store_;
            state.reset_readcells(wide());
        }

        /// Build the wide-format table from `input`, then reset the read-cells over it.
        void build_wide(const str_table &input)
        {
            wide().build(input);
            state.reset_readcells(wide());
        }
        /// Build the long-format table from `input`, matched to the wide-format one.
        void build_long(const str_table &input) { long_format().build(wide(), input); }

        /// Point every read-cell back at its imported value.
        /// - `produced_cols`: data columns produced by observations, declared alongside.
        /// @note Call before rebuilding Variables: wired cells would point at stale ones.
        void reset_readcells(const std::vector<std::string> &produced_cols = {})
        {
            wide().add_produced_cols(produced_cols);
            state.reset_readcells(wide());
        }
        /// Finalize every Variable wired into a read-cell (see data::State::finalize_readcells).
        void finalize() { state.finalize_readcells(); }

        //////  Getter  //////

        /// Read-cell `d` of member `m`: imported value or wired Variable.
        const Attribute &at(size_t d, const Individual &m) const noexcept { return state.readcell(d, m); }

        /// Wide-format table.
        data::Wide &wide() noexcept { return store_->wide; }
        const data::Wide &wide() const noexcept { return store_->wide; }
        /// Long-format table.
        data::Long &long_format() noexcept { return store_->long_format; }
        const data::Long &long_format() const noexcept { return store_->long_format; }

        /// Whether another engine shares the imported tables.
        bool is_shared() const noexcept { return store_.use_count() > 1; }

    private:
        /// Imported tables.
        struct Store_
        {
            /// Register the data keys declared by `model`.
            explicit Store_(const ModelBase &model);

            data::Wide wide;        /// Wide-format table.
            data::Long long_format; /// Long-format table.
        };
        std::shared_ptr<Store_> store_; /// Imported tables, shared copy-on-write.
    };
}
