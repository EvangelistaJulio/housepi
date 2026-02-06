/* housepi - Epidemiological Analysis of Household Outbreaks.
 *
 * Copyright © 2026 Institut Pasteur (Paris), CNRS and INSERM.
 * Distributed under the terms of the GNU General Public License (GPLv3).
 * See the LICENSE file for details.
 */

#pragma once

#include <string>
#include <string_view>
#include <vector>

#include "tools/nan.h"

namespace hsp::parse
{
    /// Wheter definitions is 'augmented'.
    inline bool is_augmented(std::string_view s) { return s == "A"; }
    /// Wheter definitions is 'data'.
    inline bool is_data_wide(std::string_view s) { return s == "D"; }
    /// Wheter definitions is 'group rate'.
    inline bool is_group_rate(std::string_view s) { return s == "G"; }
    /// Wheter definitions is 'numeric'.
    inline bool is_numeric(std::string_view s) { return s == "N"; }
    /// Wheter definitions is 'parameter'.
    inline bool is_parameter(std::string_view s) { return s == "P"; }

    /// Definition entry as an AST node.
    struct Entry
    {
        std::string func;            /// function node name if non-empty
        std::string key;             /// name of argument from parent
        std::string type;            /// terminal type
        std::string value;           /// terminal value
        bool is_braced = false;      /// wheter node was wrapped in {}
        std::vector<Entry> children; /// children if function node

        Entry() = default;
        explicit Entry(const std::string &value_, const std::string &type_);
        ~Entry() = default;

        /// Rebuild string from node with recursion.
        std::string str() const;

        /// Wheter node is a terminal node. Otherwise it is a function node.
        bool is_terminal() const noexcept { return func.empty(); }
        /// Wheter node is a function with default arguments.
        bool is_func_default() const noexcept { return !is_terminal() && children.empty(); }
        /// Wheter node is a terminal node. Otherwise it is a function node.
        bool is_empty() const noexcept { return is_terminal() && na::is(value); }
        /// Wheter is implicitely callable on argument 'point'.
        bool is_callable_implicitly(std::string_view point) const noexcept;

        /// List of dependencies per type contained in current and children definitions.
        void append_augmented_depen(std::vector<std::string> &stack) const;
        void append_data_wide_depen(std::vector<std::string> &stack) const;
        void append_parameter_depen(std::vector<std::string> &stack) const;
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

        /// List of dependencies per type contained in currenet definition.
        std::vector<std::string> list_augmented_depen() const;
        std::vector<std::string> list_data_wide_depen() const;
        std::vector<std::string> list_parameter_depen() const;

        /// Wheter there is any dependency on augmented, data, or parameters.
        bool is_augmented_depen() const { return !list_augmented_depen().empty(); }
        bool is_data_wide_depen() const { return !list_data_wide_depen().empty(); }
        bool is_parameter_depen() const { return !list_parameter_depen().empty(); }
    };
}
