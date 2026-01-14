/* housepi - Epidemiological Analysis of Household Outbreaks.
 *
 * Copyright © 2026 Institut Pasteur (Paris), CNRS and INSERM.
 * Distributed under the terms of the GNU General Public License (GPLv3).
 * See the LICENSE file for details.
 */

#ifndef HOUSEPI__TOOLS__PARSER__H
#define HOUSEPI__TOOLS__PARSER__H

#include <string>
#include <vector>

namespace hsp::parse
{
    /// Definition entry as an AST node.
    struct Entry
    {
        std::string func;            /// function node name if non-empty
        std::string key;             /// name of argument from parent
        std::string type;            /// terminal type
        std::string value;           /// terminal value
        bool is_braced = false;      /// wheter node was wrapped in {}
        std::vector<Entry> children; /// children if function node

        Entry() {}
        explicit Entry(const std::string &value_, const std::string &type_);
        ~Entry() = default;

        /// Rebuild string from node with recursion.
        std::string str() const;

        /// Wheter node is a terminal node. Otherwise it is a function node.
        bool is_terminal() const { return func.empty(); }
        /// Wheter node is a function with default arguments.
        bool is_func_default() const { return !is_terminal() && children.empty(); }
        /// Wheter node is a terminal node. Otherwise it is a function node.
        bool is_empty() const;
        /// Wheter is implicitely callable.
        bool is_callable_implicitly() const;

        /// List of dependencies per type contained in current and children definitions.
        void append_augmented_depen(std::vector<std::string> &stack) const;
        void append_data_depen(std::vector<std::string> &stack) const;
        void append_parameter_depen(std::vector<std::string> &stack) const;
    };

    /// Abstract Syntax Tree (AST) for parsing definitions of functions.
    struct Definition
    {
        Entry root; // parsed Entry root

        /// Construct and parse expression into AST (and flattens it).
        Definition(const std::string &expr);
        Definition(const char *expr) : Definition(std::string(expr)) {}
        ~Definition() {}

        /// Wheter root node is a terminal node. Otherwise it is a function node.
        bool is_terminal() const { return root.is_terminal(); }
        /// Wheter terminal value is NA.
        bool is_na() const;
        /// Wheter is implicitely callable.
        bool is_callable_implicitly() const { return root.is_callable_implicitly(); }

        /// List of dependencies per type contained in currenet definition.
        std::vector<std::string> list_augmented_depen() const;
        std::vector<std::string> list_data_depen() const;
        std::vector<std::string> list_parameter_depen() const;

        /// Wheter there is any dependency on augmented, data, or parameters.
        bool is_augmented_depen() const { return !list_augmented_depen().empty(); }
        bool is_data_depen() const { return !list_data_depen().empty(); }
        bool is_parameter_depen() const { return !list_parameter_depen().empty(); }
    };
}

#endif // HOUSEPI__TOOLS__PARSER__H
