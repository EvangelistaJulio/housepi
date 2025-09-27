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
 * parser.h
 *
 */

#ifndef HOUSEPI__PARSER__H
#define HOUSEPI__PARSER__H

#include "tools.h"

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
        ASTNode(const std::string &value_, const std::string &type_)
        {
            value = value_;
            type = tools::str::toup(type_);
        }

        /// Rebuild string from node using recursion.
        std::string str() const;
        
        /// Wheter node is a terminal node. Otherwise it is a function node.
        bool is_terminal() const { return func.empty(); }
        /// Wheter node is a function with default arguments.
        bool is_func_default() const { return !is_terminal() && children.empty(); }
        /// Wheter node is a terminal node. Otherwise it is a function node.
        bool is_empty() const { return is_terminal() && tools::is_na(value); }
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
    bool is_na() const { return is_terminal() && tools::is_na(root.value); }

    /// List of dependencies per type contained in currenet definition.
    str_vec list_augm_depen() const;
    str_vec list_data_depen() const; 
    str_vec list_param_depen() const; 

    /// Wheter there is any dependency on augmented, data, or parameters.
    bool is_augm_depen() const { return !list_augm_depen().empty(); }
    bool is_data_depen() const { return !list_data_depen().empty(); }
    bool is_param_depen() const { return !list_param_depen().empty(); }

    /// Flat iterator interface to access all nodes in the tree.
    /// NOTE: just used to explore dependecies...
    /// UPDATE: drop and provide a recursive approach through nodes.
    class ViewIterator
    {
    public:
        ViewIterator() = default;
        explicit ViewIterator(const ASTNode *root)
        {
            if (root)
                stack_.push_back({root, 0});
        }

        const ASTNode &operator*() const { return *stack_.back().first; }
        const ASTNode *operator->() const { return &(*stack_.back().first); }

        ViewIterator &operator++()
        {
            advance_();
            return *this;
        }
        ViewIterator operator++(int)
        {
            auto tmp = *this;
            ++(*this);
            return tmp;
        }

        bool operator==(const ViewIterator &other) const { return stack_ == other.stack_; }
        bool operator!=(const ViewIterator &other) const { return !(*this == other); }

    private:
        // stack holds (node pointer, child index)
        std::vector<std::pair<const ASTNode *, size_t>> stack_;

        void advance_()
        {
            if (stack_.empty())
                return;

            auto &node = stack_.back().first;
            auto &idx = stack_.back().second;

            if (idx < node->children.size())
            {
                // descend into next child
                const ASTNode *child = &node->children[idx++];
                stack_.push_back({child, 0});
            }
            else
            {
                // backtrack
                stack_.pop_back();
                if (!stack_.empty())
                    advance_();
            }
        }
    };

    ViewIterator begin() const { return ViewIterator(&root); }
    ViewIterator end() const { return ViewIterator(); }
};

#endif // HOUSEPI__PARSER__H
