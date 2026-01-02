/*
 *        __
 *       / /_  ____  __  __ ____ __  ____  __
 *      / __ \/ __ \/ / / / ___/ _ \/ __ \/ /
 *     / / / / /_/ / /_/ (__  )  __/ /_/ / /
 *    /_/ /_/\____/\____/____/\___/ ____/_/
 *                               /_/
 *
 * housepi - Epidemiological Analysis of Household Outbreaks.
 *
 * Author: Francesco Camaglia
 * Copyright © 2026 Institut Pasteur (Paris), CNRS and INSERM.
 *
 * Distributed under the terms of the GNU General Public License (GPLv3).
 * See the LICENSE file for details.
 *
 * tools/parser.h
 *
 */

#ifndef HOUSEPI__TOOLS__PARSER__H
#define HOUSEPI__TOOLS__PARSER__H

#include "tools/types.h"
#include "tools/na_utils.h"

namespace hsp
{
    namespace parse
    {
        /// Abstract Syntax Tree (AST) for parsing definitions of functions.
        class Definition
        {
        public:
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
                Entry(const std::string &value_, const std::string &type_);

                /// Rebuild string from node with recursion.
                std::string str() const;

                /// Wheter node is a terminal node. Otherwise it is a function node.
                bool is_terminal() const { return func.empty(); }
                /// Wheter node is a function with default arguments.
                bool is_func_default() const { return !is_terminal() && children.empty(); }
                /// Wheter node is a terminal node. Otherwise it is a function node.
                bool is_empty() const { return is_terminal() && hsp::na::is(value); }
                /// Wheter is implicitely callable.
                bool is_callable_implicitly() const;

                /// List of dependencies per type contained in current and children definitions.
                void append_augmented_depen(str_vec &stack) const;
                void append_data_depen(str_vec &stack) const;
                void append_parameter_depen(str_vec &stack) const;
            } root; // parsed Entry root

            /////////////////////////////////////////////////////////////////////////////////////////////
            //////////////////////////////  Definition : constructor  /////////////////////////////
            /////////////////////////////////////////////////////////////////////////////////////////////

            /// Construct and parse expression into AST (and flattens it).
            Definition(const std::string &expr);
            Definition(const char *expr) : Definition(std::string(expr)) {}

            /////////////////////////////////////////////////////////////////////////////////////////////
            ////////////////////////////////  Definition : getter  ////////////////////////////////
            /////////////////////////////////////////////////////////////////////////////////////////////

            /// Reconstruct original string.
            std::string str() const { return root.str(); }

            /// Wheter root node is a terminal node. Otherwise it is a function node.
            bool is_terminal() const { return root.is_terminal(); }
            /// Wheter terminal value is NA.
            bool is_na() const { return is_terminal() && hsp::na::is(root.value); }
            /// Wheter is implicitely callable.
            bool is_callable_implicitly() const { return root.is_callable_implicitly(); }

            /// List of dependencies per type contained in currenet definition.
            str_vec list_augmented_depen() const;
            str_vec list_data_depen() const;
            str_vec list_parameter_depen() const;

            /// Wheter there is any dependency on augmented, data, or parameters.
            bool is_augmented_depen() const { return !list_augmented_depen().empty(); }
            bool is_data_depen() const { return !list_data_depen().empty(); }
            bool is_parameter_depen() const { return !list_parameter_depen().empty(); }
        };
    }
}
#endif // HOUSEPI__TOOLS__PARSER__H
