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
 * parser.cpp
 *
 */

#include "parser.h"

/// Rebuild string from node using recursion.
std::string DefinitionParser::ASTNode::str() const
{
    if (is_empty())
        return "";
    else if (is_func_default())
        return func;
    else if (is_terminal())
        return value + "|" + type;

    std::ostringstream ss_args;
    std::ostringstream ss_children;
    for (size_t i = 0; i < children.size(); ++i)
    {
        const ASTNode &child = children[i];
        if (i > 0)
        {
            ss_args << HSP_LIST_DELIM;
            ss_children << HSP_LIST_DELIM;
        }
        ss_args << child.key;

        const std::string cs = child.str();
        if (child.is_braced)
            ss_children << "{" << cs << "}";
        else
            ss_children << cs;
    }

    return func + "|" + ss_args.str() + "|" + ss_children.str();
}

/// Split on sep at top level (ignore separators inside braces).
str_vec split_top_(const std::string &s)
{
    str_vec entries;
    std::string current;
    int level = 0;

    for (char c : s)
    {
        // check the level accroding to braces
        if (c == '{')
            ++level;
        else if (c == '}')
            --level;
        // parse only top level
        if (c == HSP_LIST_DELIM && level == 0)
        {
            entries.push_back(current);
            current.clear();
        }
        else
            current.push_back(c);
    }
    if (!current.empty())
        entries.push_back(current);

    return entries;
}

/// Parse a single expression into an AST node.
DefinitionParser::ASTNode parse_expression_(const std::string &expr)
{
    // empty node
    if (tools::is_na(expr))
        return DefinitionParser::ASTNode();

    // pure function node
    const size_t bar_1 = expr.find('|');
    if (bar_1 == std::string::npos)
    {
        DefinitionParser::ASTNode node;
        node.func = expr;
        return node;
    }

    const std::string expr_head = expr.substr(0, bar_1);
    const std::string expr_tail = expr.substr(bar_1 + 1);

    // if not found, then it's a terminal node {value|type}
    const size_t bar_2 = expr_tail.find('|');
    if (bar_2 == std::string::npos)
        return DefinitionParser::ASTNode(expr_head, expr_tail);

    // else it is a non-default function node
    // e. g. "func|arg_1:...:arg_n|child_1:...child_n"
    const str_vec args = split_top_(expr_tail.substr(0, bar_2));
    const str_vec children = split_top_(expr_tail.substr(bar_2 + 1));

    if (args.size() != children.size())
        throw std::runtime_error(
            "Mismatch in the number of keys and arguments of definition \"" + expr + "\".");

    // populate the parent node with children nodes
    DefinitionParser::ASTNode node;
    node.func = expr_head;
    for (size_t i = 0; i < children.size(); ++i)
    {
        const std::string c = children[i];

        if (c.size() >= 2 && c.front() == '{' && c.back() == '}')
        {
            // Function Child Node
            // parse inner expression and keep braced flag
            DefinitionParser::ASTNode child = parse_expression_(c.substr(1, c.size() - 2));
            child.is_braced = true;
            child.key = args[i]; // FIXME check not NA
            node.children.push_back(std::move(child));
        }
        else
        {
            // Terminal Child Node
            size_t barc = c.find('|');
            if (barc == std::string::npos)
                throw std::runtime_error(
                    "Missing terminal \'|\' in definition \"" + expr + "\".");
            DefinitionParser::ASTNode child(c.substr(0, barc), c.substr(barc + 1));
            child.key = args[i]; // FIXME check not NA
            node.children.push_back(std::move(child));
        }
    }

    return node;
}

/// Construct and parse expression into AST (and flattens it).
DefinitionParser::DefinitionParser(const std::string &expr) { root = parse_expression_(expr); }

#define MAKE_PARSER_LIST_TYPE_DEPEN(type_check, type_flag)                \
    str_vec DefinitionParser::list_##type_check##_depen() const           \
    {                                                                     \
        str_vec output = {};                                              \
        for (const auto &node : *this)                                    \
            if (node.is_terminal() && HSP_IS_TYPE_##type_flag(node.type)) \
                output.push_back(node.value);                             \
        return tools::str::unique(output);                                \
    }

MAKE_PARSER_LIST_TYPE_DEPEN(augm, AUGM)
MAKE_PARSER_LIST_TYPE_DEPEN(data, DATA)
MAKE_PARSER_LIST_TYPE_DEPEN(param, PARAM)
