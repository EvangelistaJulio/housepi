/* housepi - Epidemiological Analysis of Household Outbreaks.
 * Copyright 2024-2026 Institut Pasteur (Paris), France.
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#pragma once

#include <array>
#include <cstdint>
#include <optional>
#include <string>
#include <string_view>
#include <vector>

#include "tools/dictionary.h"
#include "tools/messages.h"
#include "tools/parser.h"
#include "data/status.h"

//////////////////////////////  augm::DefTable  //////////////////////////////

namespace housepi::augm
{
    /// Level of a definition row: latent, or one of the observation levels.
    enum class Kind : std::uint8_t
    {
        Latent,
        Member,
        House
    };

    /// Map a row level to its error-message family.
    inline ErrType kind_to_err_type(Kind kind) noexcept
    {
        switch (kind)
        {
        case Kind::Member:
            return ErrType::Observation;
        case Kind::House:
            return ErrType::HouseholdObservation;
        default:
            return ErrType::Latent;
        }
    }

    /// Definition parsed once: every source column as an AST, flags typed.
    struct DefRow
    {
        std::string name;                   /// Latent or observation name.
        Kind kind = Kind::Latent;           /// Row level.
        StatusLvl status = status::MISSING; /// Declared status, MISSING for a generic definition.

        parse::Definition value, func, model, min, max, proposal;

        bool has_value = false;
        bool has_func = false;
        bool has_model = false;
        bool has_proposal = false;
        bool is_const = false;
        bool is_discr = false;
        bool is_fcumul = false;

        /// Parse row `r` of definition table `input`.
        DefRow(const str_table &input, size_t r);

        /// Whether inference samples the row.
        /// * with a model or a proposal;
        /// * never when constant or a function.
        bool is_sampled() const noexcept { return !is_const && !has_func && (has_model || has_proposal); }
        /// Whether the row is sampled flat within its bounds, with no model.
        /// Transmission generates it in simulation, else its value pins it.
        bool is_flat() const noexcept { return is_sampled() && !has_model; }
        /// Whether a value pins the sampled row in simulation.
        bool is_sim_override() const noexcept { return has_value && is_sampled(); }

        /// List the columns that may read sources.
        std::array<const parse::Definition *, 5> source_cols() const noexcept
        {
            return {&value, &model, &func, &min, &max};
        }
        /// List the columns a run reads sources from.
        /// * inference -> every column, except the value of a sampled row;
        /// * simulation -> value and function of an override or a flat row, every column otherwise.
        std::vector<const parse::Definition *> source_cols(bool simulation) const;
        /// Select the column the definition resolves through in simulation.
        /// * with a value, an override included -> the value;
        /// * otherwise -> the function.
        const parse::Definition &sim_source() const noexcept { return has_value ? value : func; }
    };

    /// Latent and observation definitions, immutable once built.
    struct DefTable
    {
        std::vector<DefRow> rows = {};

        const DefRow &operator[](size_t k) const noexcept { return rows[k]; }
        size_t size() const noexcept { return rows.size(); }
        bool empty() const noexcept { return rows.empty(); }
        void clear() noexcept { rows.clear(); }
        /// Append definitions from `input`.
        void append(const str_table &input);

        /// Whether a row declares `name`, of level `kind` if given.
        /// @todo The level is imprecise when statuses declare `name` apart.
        bool has_row(std::string_view name, std::optional<Kind> kind = std::nullopt) const noexcept;
        /// Whether any declaration of `name` carrying a model, or sampled, is continuous.
        bool has_continuous_declaration(std::string_view name) const noexcept;
        /// Whether definition row `k` reads only data and numeric sources.
        bool is_pure_data(size_t k) const;
    };
}

namespace housepi::str
{
    inline std::string format(augm::Kind kind) { return format(augm::kind_to_err_type(kind)); }
}
