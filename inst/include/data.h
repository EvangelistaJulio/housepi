/* housepi - Epidemiological Analysis of Household Outbreaks.
 * Copyright 2024-2026 Institut Pasteur (Paris), France.
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#pragma once

#include <array>
#include <cstdint>
#include <memory>
#include <string>
#include <vector>

#include "tools/dictionary.h"
#include "data/dataset.h"
#include "data/long_dataset.h"
#include "data/state.h"

//////////////////////////////  Data  //////////////////////////////

namespace housepi
{
    class ModelBase;

    /// Data of a run.
    ///
    /// * imported tables, shared copy-on-write across engines;
    /// * run state and equivalence classes, owned by one engine.
    struct Data
    {
        data::State state; /// Values written by inference and simulation.

        /// Definitions reading data columns that tell equivalence classes apart.
        enum class Source : std::uint8_t
        {
            Model,
            GroupRate,
            ContactRate,
            Latent
        };

        /// Register the data keys declared by `model`, also equivalence-class keys.
        explicit Data(const ModelBase &model);

        //////  Setter  //////

        /// Share the imported tables of `other`, then reset the read-cells over them.
        void share(const Data &other)
        {
            store_ = other.store_;
            state.reset_readcells(dataset());
        }

        /// Build the dataset from `input`, then reset the read-cells over it.
        /// * equivalence-class keys of definitions cleared.
        void build_dataset(const str_table &input);
        /// Build the long dataset from `input`, matched to the dataset.
        void build_long_dataset(const str_table &input) { long_dataset().build(dataset(), input); }

        /// Point every read-cell back at its imported value.
        /// * `produced_cols`: data columns produced by observations, declared alongside.
        /// @note Call before rebuilding Variables: wired cells would point at stale ones.
        void reset_readcells(const std::vector<std::string> &produced_cols = {})
        {
            dataset().add_produced_cols(produced_cols);
            state.reset_readcells(dataset());
        }
        /// Finalize every Variable wired into a read-cell (see data::State::finalize_readcells).
        void finalize() { state.finalize_readcells(); }

        //////  Equivalence classes  //////

        /// Set the equivalence-class keys of `source` to data columns `names`.
        void set_eq_class_keys(Source source, std::vector<std::string> names)
        {
            eq_class_keys_[static_cast<size_t>(source)] = std::move(names);
        }
        /// Partition households into equivalence classes.
        /// * equivalent households: same status and key cells per member, in sorted member order;
        /// * equivalent households share one likelihood, computed on the class representative;
        /// * `force_singleton` -> one class per household;
        /// * keys read through the read-cells of `state`.
        void init_eq_classes(bool force_singleton);

        /// Number of equivalence classes.
        size_t n_eq_classes() const noexcept { return state.eq_classes_.size(); }
        /// Number of households in equivalence class `R`.
        size_t eq_class_size(size_t R) const noexcept { return state.eq_classes_[R].size_; }
        /// Representative household of equivalence class `R`.
        Household &eq_class_rep(size_t R) noexcept { return dataset().at(state.eq_classes_[R].rep_); }
        const Household &eq_class_rep(size_t R) const noexcept { return dataset().at(state.eq_classes_[R].rep_); }
        /// Member of `house` paired with each member of its class representative, by representative index.
        /// * equivalent members pair up in sorted order, so ties pair arbitrarily.
        std::vector<MemberIdx> map_rep_members(const Household &house) const;

        //////  Getter  //////

        /// Read-cell `d` of member `m`: imported value or wired Variable.
        const Attribute &at(size_t d, const Individual &m) const noexcept { return state.readcell(d, m); }

        /// Dataset: one row per member.
        data::Dataset &dataset() noexcept { return store_->dataset; }
        const data::Dataset &dataset() const noexcept { return store_->dataset; }
        /// Long dataset: any number of rows per member.
        data::LongDataset &long_dataset() noexcept { return store_->long_dataset; }
        const data::LongDataset &long_dataset() const noexcept { return store_->long_dataset; }

        /// Whether another engine shares the imported tables.
        bool is_shared() const noexcept { return store_.use_count() > 1; }

    private:
        /// Imported tables.
        struct Store_
        {
            /// Register the data keys declared by `model`.
            explicit Store_(const ModelBase &model);

            data::Dataset dataset;          /// Dataset: one row per member.
            data::LongDataset long_dataset; /// Long dataset: any number of rows per member.
        };
        std::shared_ptr<Store_> store_; /// Imported tables, shared copy-on-write.

        /// Data columns telling equivalence classes apart, by source.
        std::array<std::vector<std::string>, 4> eq_class_keys_ = {};
    };
}
