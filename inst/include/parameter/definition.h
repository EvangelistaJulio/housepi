/* housepi - Epidemiological Analysis of Household Outbreaks.
 * Copyright 2024-2026 Institut Pasteur (Paris), France.
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#pragma once

#include <string>
#include <vector>

#include "tools/dictionary.h"
#include "tools/parser.h"

//////////////////////////////  param::DefTable  //////////////////////////////

namespace housepi::param
{
    /// Parameter definition parsed once: every source column as an AST, flags typed.
    struct DefRow
    {
        std::string name; /// Parameter name.

        parse::Definition value, prior, func, min, max, proposal;

        bool is_discr = false;

        /// Parse row `r` of definition table `input`.
        /// * function reading a latent or a data column -> error.
        DefRow(const str_table &input, size_t r);

        /// Whether the row holds a function.
        bool has_func() const noexcept { return !func.is_na(); }
        /// Whether the row holds a constant.
        bool is_const() const noexcept { return !value.is_na(); }
    };

    /// Parameter definitions, immutable once built.
    struct DefTable
    {
        std::vector<DefRow> rows = {};

        const DefRow &operator[](size_t k) const noexcept { return rows[k]; }
        /// Replace definitions with those of `input`.
        void build(const str_table &input);

        /// List rows producer before consumer, through `func`, `min` and `max`.
        /// * no function parameter -> row order.
        std::vector<size_t> order() const;
    };
}
