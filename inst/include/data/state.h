/* housepi - Epidemiological Analysis of Household Outbreaks.
 * Copyright 2024-2026 Institut Pasteur (Paris), France.
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#pragma once

#include <vector>

#include "tools/logproba.h"
#include "tools/vector.h"
#include "data/household.h"
#include "data/individual.h"
#include "function.h"

//////////////////////////////  data::State  //////////////////////////////

namespace housepi::data
{
#define MAKE_HOUSEPI_CACHE_ACCESSOR(TYPE, NAME)                                                                  \
    TYPE &cache_##NAME(const Household &house, size_t e) noexcept { return current_[house.eq_class()].NAME[e]; } \
    const TYPE &cache_##NAME(const Household &house, size_t e) const noexcept { return current_[house.eq_class()].NAME[e]; }

    struct Wide;

    /// Mutable per-run storage of data, written by inference and simulation.
    ///
    /// * read-cells, per member and column;
    /// * transmission log-likelihood, per equivalence class;
    /// * model cache, current and last saved, per equivalence class;
    /// * whether the last two hold an evaluation, per equivalence class.
    struct State
    {
        //////  Setter  //////

        /// Point every read-cell at the imported value of `wide`.
        /// @note Call before rebuilding Variables: wired cells would point at stale ones.
        void reset_readcells(const Wide &wide);
        /// Point read-cell `d` of member `indiv` at the Variable of a surrogate column.
        void wire_readcell(size_t d, const Individual &indiv, const Attribute *p) noexcept
        {
            readcell_[Index::of(indiv)][d].init_pointer(p);
        }
        /// Finalize every Variable wired into a read-cell.
        /// * evaluation only -> simulation re-wires the cells live.
        void finalize_readcells();

        /// Allocate transmission log-likelihoods, one per equivalence class of `wide`.
        /// * clears the evaluated flags.
        void allocate_transm_log_lik(const Wide &wide);
        /// Allocate model caches, per equivalence class of `wide`, sized by household and cache keys.
        /// * clears the evaluated flags.
        void allocate_cache(const Wide &wide);
        /// Flag whether transmission log-likelihoods and model caches hold an evaluation of every class.
        void set_evaluated(bool evaluated) noexcept { evaluated_.assign(evaluated_.size(), evaluated); }
        /// Flag the equivalence class of `house` as evaluated.
        void set_evaluated(const Household &house) noexcept { evaluated_[house.eq_class()] = true; }

        /// Save every model cache.
        void save_cache() { previous_ = current_; }
        /// Recover every model cache from the last save.
        void recover_cache() noexcept { current_.swap(previous_); }
        /// Save the model cache of member `indiv`.
        void save_member_cache(const Household &house, const Individual &indiv);
        /// Recover the model cache of member `indiv` from the last save.
        void recover_member_cache(const Household &house, const Individual &indiv);

        //////  Getter  //////

        /// Read-cell `d` of member `indiv`: imported value or wired Variable, resolved live by a data read.
        const Attribute &readcell(size_t d, const Individual &indiv) const noexcept
        {
            return readcell_[Index::of(indiv)][d];
        }

        /// Transmission log-likelihood of the equivalence class of `house`.
        LogProba &transm_log_lik(const Household &house) noexcept { return transm_log_lik_[house.eq_class()]; }
        /// Whether transmission log-likelihood and model cache hold an evaluation of the class of `house`.
        /// * class not allocated in this state -> false.
        bool is_evaluated(const Household &house) const noexcept
        {
            return house.eq_class() < evaluated_.size() && evaluated_[house.eq_class()];
        }

        /// Cached member value.
        MAKE_HOUSEPI_CACHE_ACCESSOR(MemberIdx, member)
        /// Cached member vector.
        MAKE_HOUSEPI_CACHE_ACCESSOR(std::vector<MemberIdx>, member_vector)
        /// Cached double vector.
        MAKE_HOUSEPI_CACHE_ACCESSOR(std::vector<double>, double_vector)
        /// Cached double matrix.
        MAKE_HOUSEPI_CACHE_ACCESSOR(matrix, double_matrix)

    private:
        /// Model cache of one equivalence class, by cache key.
        struct Cache_
        {
            std::vector<MemberIdx> member = {};                     /// Member values.
            std::vector<std::vector<MemberIdx>> member_vector = {}; /// Member vectors.
            std::vector<std::vector<double>> double_vector = {};    /// Double vectors.
            std::vector<matrix> double_matrix = {};                 /// Double matrices.
        };

        std::vector<std::vector<Attribute>> readcell_ = {}; /// Read-cells, by member index and column.
        std::vector<LogProba> transm_log_lik_ = {};         /// Transmission log-likelihood, per equivalence class.
        std::vector<Cache_> current_ = {};                  /// Model cache, per equivalence class.
        std::vector<Cache_> previous_ = {};                 /// Model cache at the last save, per equivalence class.
        std::vector<bool> evaluated_ = {};                  /// Whether evaluated, per equivalence class.
    };
}

#undef MAKE_HOUSEPI_CACHE_ACCESSOR
