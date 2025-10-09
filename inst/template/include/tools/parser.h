/*
 *        __
 *       / /_  ____  __  __ ____ __  ____  __
 *      / __ \/ __ \/ / / / ___/ _ \/ __ \/ /
 *     / / / / /_/ / /_/ (__  )  __/ /_/ / /
 *    /_/ /_/\____/\____/____/\___/ ____/_/
 *                               /_/
 *
 * Copyright (C) 2025 Francesco Camaglia, Institut Pasteur
 *
 * tools/parser.h
 *
 */

#ifndef HOUSEPI__TOOLS__PARSER__H
#define HOUSEPI__TOOLS__PARSER__H

#include "tools/types.h"
#include "tools/na_utils.h"

/// Abstract Syntax Tree (AST) for parsing definitions of functions.
class DefinitionParser
{
public:
    /// AST node.
    struct ASTNode
    {
        std::string func;              /// function node name if non-empty
        std::string key;               /// name of argument from parent
        std::string type;              /// terminal type
        std::string value;             /// terminal value
        bool is_braced = false;        /// wheter node was wrapped in {}
        std::vector<ASTNode> children; /// children if function node

        ASTNode() {}
        ASTNode(const std::string &value_, const std::string &type_);

        /// Rebuild string from node using recursion.
        std::string str() const;

        /// Wheter node is a terminal node. Otherwise it is a function node.
        bool is_terminal() const { return func.empty(); }
        /// Wheter node is a function with default arguments.
        bool is_func_default() const { return !is_terminal() && children.empty(); }
        /// Wheter node is a terminal node. Otherwise it is a function node.
        bool is_empty() const { return is_terminal() && hsp::na::is(value); }
        /// Wheter is implicitely callable.
        bool is_implicit_callable() const;

        /// List of dependencies per type contained in current and children definitions.
        void append_augmented_depen(str_vec &stack) const;
        void append_data_depen(str_vec &stack) const;
        void append_parameter_depen(str_vec &stack) const;
    } root; // parsed AST root

    /////////////////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////  DefinitionParser : constructor  /////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////

    /// Construct and parse expression into AST (and flattens it).
    DefinitionParser(const std::string &expr);
    DefinitionParser(const char *expr) : DefinitionParser(std::string(expr)) {}

    /////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////  DefinitionParser : getter  ////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////

    /// Reconstruct original string.
    std::string str() const { return root.str(); }

    /// Wheter root node is a terminal node. Otherwise it is a function node.
    bool is_terminal() const { return root.is_terminal(); }
    /// Wheter terminal value is NA.
    bool is_na() const { return is_terminal() && hsp::na::is(root.value); }
    /// Wheter is implicitely callable.
    bool is_implicit_callable() const { return root.is_implicit_callable(); }

    /// List of dependencies per type contained in currenet definition.
    str_vec list_augmented_depen() const;
    str_vec list_data_depen() const;
    str_vec list_parameter_depen() const;

    /// Wheter there is any dependency on augmented, data, or parameters.
    bool is_augmented_depen() const { return !list_augmented_depen().empty(); }
    bool is_data_depen() const { return !list_data_depen().empty(); }
    bool is_parameter_depen() const { return !list_parameter_depen().empty(); }
};

#endif // HOUSEPI__TOOLS__PARSER__H
