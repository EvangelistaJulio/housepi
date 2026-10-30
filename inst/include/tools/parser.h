/* housepi - Epidemiological Analysis of Household Outbreaks.
 * Copyright 2024-2026 Institut Pasteur (Paris), France.
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#pragma once

#include <string>
#include <string_view>
#include <vector>

#include "tools/nan.h"

namespace housepi::parse
{
    /// Single-letter code tagging a definition's terminal type.
    enum class Type : char
    {
        Latent = 'A',
        Dataset = 'D',
        GroupRate = 'G',
        Numeric = 'N',
        Parameter = 'P',
    };

    /// Whether type code `s` equals tag `t`.
    inline bool is_(std::string_view s, Type t) { return s.size() == 1 && s.front() == static_cast<char>(t); }

    /// Whether type code `s` is latent.
    inline bool is_latent(std::string_view s) { return is_(s, Type::Latent); }
    /// Whether type code `s` is a dataset column.
    inline bool is_dataset(std::string_view s) { return is_(s, Type::Dataset); }
    /// Whether type code `s` is group rate.
    inline bool is_group_rate(std::string_view s) { return is_(s, Type::GroupRate); }
    /// Whether type code `s` is numeric.
    inline bool is_numeric(std::string_view s) { return is_(s, Type::Numeric); }
    /// Whether type code `s` is parameter.
    inline bool is_parameter(std::string_view s) { return is_(s, Type::Parameter); }

    /// Definition entry as an AST node.
    struct Entry
    {
        std::string func;            /// Function name, empty for a terminal.
        std::string key;             /// Argument name in the parent.
        std::string type;            /// Terminal type.
        std::string value;           /// Terminal value.
        std::vector<Entry> children; /// Arguments of a function node.

        Entry() = default;
        explicit Entry(const std::string &value, const std::string &type);
        ~Entry() = default;

        /// Rebuild the expression string, recursively.
        std::string str() const;
        /// Human-readable string, for errors.
        std::string err_str() const;

        /// Whether the node is a terminal, not a function.
        bool is_terminal() const noexcept { return func.empty(); }
        /// Whether the node is a function with default arguments.
        bool is_func_default() const noexcept { return !is_terminal() && children.empty(); }
        /// Whether the node is a terminal holding NA.
        bool is_empty() const noexcept { return is_terminal() && na::is(value); }
        /// Whether the node is callable implicitly on argument `point`.
        bool is_callable_implicitly(std::string_view point) const noexcept;

        /// Append producer names of a given type, consumed by this node and its arguments.
        void append_latent_producers(std::vector<std::string> &stack) const;
        void append_dataset_producers(std::vector<std::string> &stack) const;
        void append_parameter_producers(std::vector<std::string> &stack) const;
    };

    /// Abstract syntax tree of a definition.
    struct Definition
    {
        Entry root; /// Root node.

        /// Parse expression `expr` into a flattened tree.
        Definition(std::string_view expr);
        ~Definition() = default;

        /// Whether the root is a terminal, not a function.
        bool is_terminal() const noexcept { return root.is_terminal(); }
        /// Whether the terminal value is NA.
        bool is_na() const noexcept { return is_terminal() && na::is(root.value); }
        /// Whether the root is callable implicitly on argument `point`.
        bool is_callable_implicitly(std::string_view point) const noexcept { return root.is_callable_implicitly(point); }

        /// Names of the given type this definition consumes, i.e. its producers.
        std::vector<std::string> list_latent_producers() const;
        std::vector<std::string> list_dataset_producers() const;
        std::vector<std::string> list_parameter_producers() const;

        /// Whether this definition consumes any latent, data or parameter name.
        bool consumes_latent() const { return !list_latent_producers().empty(); }
        bool consumes_dataset() const { return !list_dataset_producers().empty(); }
        bool consumes_parameter() const { return !list_parameter_producers().empty(); }
    };
}
