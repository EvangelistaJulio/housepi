/* housepi - Epidemiological Analysis of Household Outbreaks.
 *
 * Copyright © 2024-2026 Institut Pasteur (Paris), CNRS and INSERM.
 * Distributed under the terms of the GNU General Public License (GPLv3).
 * See the LICENSE file for details.
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
        DataWide = 'D',
        GroupRisk = 'G',
        Numeric = 'N',
        Parameter = 'P',
    };

    /// Wheter type code equals the given tag.
    inline bool is_(std::string_view s, Type t) { return s.size() == 1 && s.front() == static_cast<char>(t); }

    /// Wheter definitions is 'latent'.
    inline bool is_latent(std::string_view s) { return is_(s, Type::Latent); }
    /// Wheter definitions is 'data'.
    inline bool is_data_wide(std::string_view s) { return is_(s, Type::DataWide); }
    /// Wheter definitions is 'group risk'.
    inline bool is_group_risk(std::string_view s) { return is_(s, Type::GroupRisk); }
    /// Wheter definitions is 'numeric'.
    inline bool is_numeric(std::string_view s) { return is_(s, Type::Numeric); }
    /// Wheter definitions is 'parameter'.
    inline bool is_parameter(std::string_view s) { return is_(s, Type::Parameter); }

    /// Definition entry as an AST node.
    struct Entry
    {
        std::string func;            /// function node name if non-empty
        std::string key;             /// name of argument from parent
        std::string type;            /// terminal type
        std::string value;           /// terminal value
        bool is_bracketed = false;   /// wheter node was wrapped in []
        std::vector<Entry> children; /// arguments of a function node.

        Entry() = default;
        explicit Entry(const std::string &value_, const std::string &type_);
        ~Entry() = default;

        /// Rebuild string from node with recursion.
        std::string str() const;
        /// Human readable string used for errors.
        std::string err_str() const;

        /// Wheter node is a terminal node. Otherwise it is a function node.
        bool is_terminal() const noexcept { return func.empty(); }
        /// Wheter node is a function with default arguments.
        bool is_func_default() const noexcept { return !is_terminal() && children.empty(); }
        /// Wheter node is a terminal node. Otherwise it is a function node.
        bool is_empty() const noexcept { return is_terminal() && na::is(value); }
        /// Wheter is implicitely callable on argument 'point'.
        bool is_callable_implicitly(std::string_view point) const noexcept;

        /// Append producer names of a given type, consumed by this node and its arguments.
        void append_latent_producers(std::vector<std::string> &stack) const;
        void append_data_wide_producers(std::vector<std::string> &stack) const;
        void append_parameter_producers(std::vector<std::string> &stack) const;
    };

    /// Abstract Syntax Tree (AST) for parsing definitions of functions.
    struct Definition
    {
        Entry root; // parsed Entry root

        /// Construct and parse expression into AST (and flattens it).
        Definition(std::string_view expr);
        ~Definition() = default;

        /// Wheter root node is a terminal node. Otherwise it is a function node.
        bool is_terminal() const noexcept { return root.is_terminal(); }
        /// Wheter terminal value is NA.
        bool is_na() const noexcept { return is_terminal() && na::is(root.value); }
        /// Wheter is implicitely callable.
        bool is_callable_implicitly(std::string_view point) const noexcept { return root.is_callable_implicitly(point); }

        /// Names of the given type this definition consumes, i.e. its producers.
        std::vector<std::string> list_latent_producers() const;
        std::vector<std::string> list_data_wide_producers() const;
        std::vector<std::string> list_parameter_producers() const;

        /// Wheter this definition consumes any latent, data, or parameter name.
        bool consumes_latent() const { return !list_latent_producers().empty(); }
        bool consumes_data_wide() const { return !list_data_wide_producers().empty(); }
        bool consumes_parameter() const { return !list_parameter_producers().empty(); }
    };
}
