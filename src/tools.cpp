/* housepi - Epidemiological Analysis of Household Outbreaks.
 * Copyright 2024-2026 Institut Pasteur (Paris), France.
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#include <sstream>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <string>
#include <unordered_set>
#include <vector>

#include "tools/dag.h"
#include "tools/keys.h"
#include "tools/math.h"
#include "tools/messages.h"
#include "tools/nan.h"
#include "tools/parser.h"
#include "tools/string.h"
#include "tools/vector.h"
#include "function.h"
#include "variable.h"

//////////////////////////////  KeyHolder  //////////////////////////////

namespace housepi
{
    void KeyHolder::clear_optional_() { keys_.resize(n_declared_); }

    void KeyHolder::add_(std::string_view name)
    {
        if (!count(name))
            keys_.emplace_back(name);
    }

    std::string KeyHolder::consult(size_t index) const
    {
        if (index >= keys_.size())
            msg::stop("Registration failed: key index out of bounds.");
        return keys_[index];
    }
    std::vector<std::string> KeyHolder::consult(const std::vector<size_t> &indexes) const
    {
        std::vector<std::string> output;
        output.reserve(indexes.size());
        for (size_t s : indexes)
            output.push_back(consult(s));
        return output;
    }

    size_t KeyHolder::consult(std::string_view name) const
    {
        const size_t k = find(name);
        if (na::is(k))
            msg::stop("Registration failed: key name " + msg::val(std::string(name)) + " was not found.");
        return k;
    }
    std::vector<size_t> KeyHolder::consult(const std::vector<std::string> &names) const
    {
        std::vector<size_t> output;
        output.reserve(names.size());
        for (const std::string &s : names)
            output.push_back(consult(s));
        return output;
    }

    size_t KeyHolder::find(std::string_view name) const noexcept
    {
        const auto it = std::find(keys_.begin(), keys_.end(), name);
        return it == keys_.end() ? na::INDEX : static_cast<size_t>(it - keys_.begin());
    }
    std::vector<size_t> KeyHolder::find(const std::vector<std::string> &names) const
    {
        std::vector<size_t> output;
        output.reserve(names.size());
        for (const std::string &s : names)
            output.push_back(find(s));
        return output;
    }

    void KeyHolderWithRequests::request_discrete(std::vector<size_t> vec)
    {
        vect::append_unique(discr_, vec);
    }
    void KeyHolderWithRequests::request_positive(std::vector<size_t> vec)
    {
        vect::append_unique(posit_, vec);
    }
    bool KeyHolderWithRequests::is_required_discrete(size_t p) const noexcept
    {
        return vect::contains(discr_, p);
    }
    bool KeyHolderWithRequests::is_required_positive(size_t p) const noexcept
    {
        return vect::contains(posit_, p);
    }

    template <class Def>
    void KeyHolderWithRequests::assert_requests(size_t p, const Def &def) const
    {
        if (is_required_discrete(p) && !def.is_discrete())
            msg::stop(def.str() +
                      " was not discrete-valued, as requested by the model.");
        if (is_required_positive(p) && !def.is_positive())
            msg::stop(def.str() +
                      " was not strictly positive, as requested by the model.");
    }
    template void KeyHolderWithRequests::assert_requests(size_t, const Attribute &) const;
    template void KeyHolderWithRequests::assert_requests(size_t, const Variable &) const;
}

///////////////////////////  messages  ///////////////////////////

namespace housepi::msg
{
    namespace
    {
        /// Innermost warning buffer of the calling thread.
        thread_local WarningBuffer *innermost_ = nullptr;
    }

    WarningBuffer::WarningBuffer() noexcept
        : outer_(innermost_) { innermost_ = this; }
    WarningBuffer::~WarningBuffer()
    {
        if (is_open_)
            innermost_ = outer_;
    }

    bool WarningBuffer::hold(const std::string &txt)
    {
        if (!innermost_)
            return false;
        if (!vect::contains(innermost_->warnings_, txt))
            innermost_->warnings_.push_back(txt);
        return true;
    }

    void WarningBuffer::flush()
    {
        if (is_open_)
            innermost_ = outer_;
        is_open_ = false;
        for (const std::string &txt : warnings_)
            warn(txt);
        warnings_.clear();
    }
}

///////////////////////////  na values  ///////////////////////////

namespace housepi::na
{
    bool is(std::string_view txt) noexcept
    {
        return txt.empty() ||
               txt == "<NA>" ||
               txt == "NA" ||
               txt == "N/A" ||
               txt == "NaN" ||
               txt == "NAN" ||
               txt == "nan";
    }
}

////////////////////////////  strings  ////////////////////////////

namespace housepi::str
{
    //////  Converters  //////

    std::string toup(std::string s)
    {
        std::transform(s.begin(), s.end(), s.begin(),
                       [](unsigned char c)
                       { return std::toupper(c); });
        return s;
    }

    bool tobool(const std::string &tkn)
    {
        if (tkn == "FALSE" || tkn == "false" || tkn == "False" || tkn == "0")
            return false;
        if (tkn == "TRUE" || tkn == "true" || tkn == "True" || tkn == "1")
            return true;
        return false;
    }

    namespace
    {
        /// Call conversion `func` on `input`, raising an invalid argument as an error naming `errtxt`.
        template <typename Func, typename... Args>
        auto helper_converter_(Func func, const std::string &input, const std::string &errtxt, Args &&...args) -> decltype(func(input, std::forward<Args>(args)...))
        {
            try
            {
                return func(input, std::forward<Args>(args)...);
            }
            catch (const std::invalid_argument &)
            {
                msg::stop(
                    (errtxt.empty() ? "" : "Unable to import " + errtxt + ". ") +
                    "Value " + msg::val(input) + " is not suitable to conversion.");
            }
            catch (const std::out_of_range &)
            {
                msg::stop(
                    (errtxt.empty() ? "" : "Unable to import " + errtxt + ". ") +
                    "Value " + msg::val(input) + " is out of range.");
            }
        }
    }

    double tod(const std::string &input, const std::string &errtxt)
    {
        if (na::is(input))
            return na::DOUBLE;
        // @todo Replace the literal match on "pi" with a lookup of named constants.
        if (input == "pi")
            return M_PI;
        return helper_converter_(static_cast<double (*)(const std::string &, size_t *)>(std::stod), input, errtxt, nullptr);
    }
    int toi(const std::string &input, const std::string &errtxt)
    {
        return helper_converter_(static_cast<int (*)(const std::string &, size_t *, int)>(std::stoi), input, errtxt, nullptr, 10);
    }
    std::uint64_t toull(const std::string &input, const std::string &errtxt)
    {
        return helper_converter_(static_cast<unsigned long long (*)(const std::string &, size_t *, int)>(std::stoull), input, errtxt, nullptr, 10);
    }

    //////  Formatters  //////

    std::string format(bool condition) { return condition ? "TRUE" : "FALSE"; }

    std::string format(double x)
    {
        if (na::is(x))
            return "NA";
        char buf[32];
        std::snprintf(buf, sizeof buf, "%g", x);
        return buf;
    }

    std::string format(const std::vector<std::string> &options,
                       bool add_full_stop,
                       size_t max_print)
    {
        if (options.empty())
            return "";

        max_print = std::min(options.size(), max_print);
        std::string output = msg::val(options[0]);
        for (size_t i = 1; i < max_print; ++i)
            output += ", " + msg::val(options[i]);
        if (options.size() > max_print)
            output += ", ...";
        else if (add_full_stop)
            output += ".";
        return output;
    }

    std::string format(ErrType type) noexcept
    {
        switch (type)
        {
        case ErrType::Numeric:
            return "numeric";
        case ErrType::Pointer:
            return "pointer";
        case ErrType::Function:
            return "function";
        case ErrType::Operator:
            return "operator";
        case ErrType::Parameter:
            return "parameter";
        case ErrType::Latent:
            return "latent";
        case ErrType::Observation:
            return "observation";
        case ErrType::HouseholdObservation:
            return "household observation";
        case ErrType::GroupRate:
            return "group rate";
        case ErrType::Data:
            return "data";
        case ErrType::Distribution:
            return "distribution";
        case ErrType::InfecProfile:
            return "infectivity profile";
        case ErrType::Proposal:
            return "proposal";
        default: // case ErrType::Empty:
            return "empty";
        }
    }

    //////  String vectors  //////

    bool starts_with(std::string_view s, std::string_view pattern) noexcept
    {
        return s.size() >= pattern.size() && s.compare(0, pattern.size(), pattern) == 0;
    }

    std::vector<std::string> unique(const std::vector<std::string> &vec)
    {
        std::vector<std::string> out;
        out.reserve(vec.size());
        std::unordered_set<std::string> seen;
        seen.reserve(vec.size());
        for (const std::string &s : vec)
            if (seen.insert(s).second)
                out.push_back(s);
        return out;
    }

    namespace
    {
        /// Return Levenshtein distance between `seq_A` and `seq_B`.
        /// * one empty sequence -> length of the other.
        /// Based on @cite wikipedia_levenshtein_2026
        size_t helper_levenshtein_(std::string_view seq_A, std::string_view seq_B)
        {
            size_t len_A = seq_A.size();
            size_t len_B = seq_B.size();
            std::vector<size_t> prev(len_B + 1), curr(len_B + 1);

            for (size_t b = 0; b <= len_B; ++b)
                prev[b] = b;

            for (size_t a = 1; a <= len_A; ++a)
            {
                curr[0] = a;
                for (size_t b = 1; b <= len_B; ++b)
                {
                    size_t cost = (seq_A[a - 1] == seq_B[b - 1]) ? 0 : 1;
                    curr[b] = std::min({prev[b] + 1, curr[b - 1] + 1, prev[b - 1] + cost});
                }
                std::swap(prev, curr);
            }

            return prev[len_B];
        }
    }

    std::string find_closest(std::string_view input,
                             const std::vector<std::string> &options,
                             std::string_view otherwise,
                             size_t min_dist)
    {
        std::string closest = "";
        // shorter inputs tolerate fewer edits
        size_t dist_closest = std::min(min_dist, input.size() / 2);

        if (dist_closest > 0)
            for (const std::string &x : options)
            {
                const size_t d = helper_levenshtein_(input, x);
                if (d <= dist_closest)
                {
                    dist_closest = d;
                    closest = x;
                }
            }

        if (!closest.empty())
            return "\nDid you mean " + msg::val(closest) + "?";
        return "\n" + std::string(otherwise);
    }

    //////  List parsers  //////

    std::string clean(std::string tkn)
    {
        // forbidden characters
        for (char r : {'"', '\r', '\n'})
            tkn.erase(remove(tkn.begin(), tkn.end(), r), tkn.end());
        return tkn;
    }

    std::vector<std::string> split_list(const std::string &tkn)
    {
        std::vector<std::string> output = {};
        std::istringstream args_ss(clean(tkn));
        std::string tmp;

        // read the entries of the list
        while (std::getline(args_ss, tmp, LIST_DELIMITER))
            output.push_back(tmp);
        return output;
    }

    std::vector<std::string> parse_names(std::string input)
    {
        std::vector<std::string> result;
        size_t start = 0;

        while (start < input.size())
        {
            // skip whitespace
            while (start < input.size() &&
                   std::isspace(static_cast<unsigned char>(input[start])))
                ++start;

            if (start >= input.size())
                break;

            // find end of identifier (alphanumeric + underscore)
            size_t end = start;
            while (end < input.size() &&
                   (std::isalnum(static_cast<unsigned char>(input[end])) || input[end] == '_'))
                ++end;

            // add the identifier if found
            if (end > start)
                result.push_back(input.substr(start, end - start));

            // skip to next comma or assignment
            while (end < input.size() && input[end] != ',')
                ++end;

            start = end + 1;
        }

        return result;
    }
}

////////////////////////  directed acyclic graphs  ////////////////////////

namespace housepi::dag
{
    void walk_graph(const std::string &entry,
                    const prod_graph &graph,
                    std::vector<std::string> &visited,
                    std::vector<std::string> &recursion_stack,
                    std::vector<std::string> &ordered,
                    ErrType err)
    {
        if (!graph.count(entry))
            msg::stop(
                "Could not find the " + str::format(err) +
                " definition " + msg::val(entry) + ".");

        if (vect::contains(recursion_stack, entry))
            msg::stop(
                "Found inadmissible cyclic " + str::format(err) +
                " definition " + msg::val(entry) + ".");

        if (vect::contains(visited, entry))
            return;

        recursion_stack.push_back(entry);
        for (const std::string &next : graph.at(entry))
            walk_graph(next, graph, visited, recursion_stack, ordered, err);

        // the cycle guard keeps 'entry' unique on the stack, nested calls balance out
        // 'entry' is the top on return
        recursion_stack.pop_back();
        visited.push_back(entry);
        ordered.push_back(entry);
    }

    std::vector<std::string> order_producers_first(const prod_graph &graph, ErrType err)
    {
        std::vector<std::string> visited, recursion_stack, ordered;
        for (const auto &[a, vec] : graph)
            if (!vect::contains(visited, a))
                walk_graph(a, graph, visited, recursion_stack, ordered, err);
        return ordered;
    }

    std::vector<std::string> transitive_producers(const std::string &start,
                                                  const prod_graph &graph,
                                                  ErrType err)
    {
        // post-order: producers first, 'start' last
        std::vector<std::string> visited, recursion_stack, ordered;
        walk_graph(start, graph, visited, recursion_stack, ordered, err);
        return ordered;
    }

    std::vector<std::string> transitive_consumers(const std::string &start,
                                                  const prod_graph &graph,
                                                  ErrType err)
    {
        // consumer graph: each node maps to the nodes reading it
        prod_graph consumer_graph;
        for (const auto &[node, producers] : graph)
        {
            for (const std::string &prod : producers)
                consumer_graph[prod].insert(node);

            // ensure isolated nodes are present
            consumer_graph.try_emplace(node);
        }

        // walk downstream of start
        std::vector<std::string> visited, recursion_stack, ordered;
        for (const std::string &cons : consumer_graph.at(start))
            walk_graph(cons, consumer_graph, visited, recursion_stack, ordered, err);

        std::reverse(ordered.begin(), ordered.end());
        return vect::concatenate(std::vector<std::string>({start}), ordered);
    }
}

//////////////////////////////  parse  //////////////////////////////

/// Define `Entry::append_<type>_producers`, appending producers of the type to `stack`.
#define MAKE_NODE_LIST_PRODUCERS(type_flag)                                           \
    void Entry::append_##type_flag##_producers(std::vector<std::string> &stack) const \
    {                                                                                 \
        if (is_terminal() && is_##type_flag(type))                                    \
            stack.push_back(value);                                                   \
        else                                                                          \
            for (const auto &child : children)                                        \
                child.append_##type_flag##_producers(stack);                          \
    }

/// Define `Definition::list_<type>_producers`, listing distinct producers of the type.
#define MAKE_PARSER_LIST_PRODUCERS(type_flag)                                 \
    std::vector<std::string> Definition::list_##type_flag##_producers() const \
    {                                                                         \
        std::vector<std::string> stack = {};                                  \
        root.append_##type_flag##_producers(stack);                           \
        return str::unique(stack);                                            \
    }

namespace housepi::parse
{
    namespace
    {
        /// Split `s` at top-level delimiters, skipping those inside brackets.
        ///
        /// * n delimiters give n + 1 entries, empty ones included: an unnamed argument keeps its key slot.
        std::vector<std::string> split_top_(const std::string &s)
        {
            std::vector<std::string> entries;
            std::string current;
            int level = 0;

            for (char c : s)
            {
                // bracket depth
                if (c == '[')
                    ++level;
                else if (c == ']')
                    --level;
                // parse only top level
                if (c == str::LIST_DELIMITER && level == 0)
                {
                    entries.push_back(current);
                    current.clear();
                }
                else
                    current.push_back(c);
            }
            entries.push_back(current);

            return entries;
        }

        /// Parse a single expression into an AST node.
        Entry parse_expression_(std::string_view input)
        {
            // empty node
            if (na::is(input))
                return Entry();

            // pure function node
            const std::string expr(input);
            size_t bar_1 = expr.find('|');
            if (bar_1 == std::string::npos)
            {
                Entry node;
                node.func = expr;
                return node;
            }

            const std::string expr_head = expr.substr(0, bar_1);
            const std::string expr_tail = expr.substr(bar_1 + 1);

            // no second bar: terminal node "value|type"
            size_t bar_2 = expr_tail.find('|');
            if (bar_2 == std::string::npos)
                return Entry(expr_head, expr_tail);

            // non-default function node, e.g. "func|arg_1:...:arg_n|child_1:...child_n"
            const std::vector<std::string> args = split_top_(expr_tail.substr(0, bar_2));
            const std::vector<std::string> children = split_top_(expr_tail.substr(bar_2 + 1));

            // children nodes
            Entry node;
            node.func = expr_head;
            node.children.reserve(children.size());
            // every child is wrapped in []
            for (size_t i = 0; i < children.size(); ++i)
            {
                const std::string &c = children[i];
                Entry child = parse_expression_(c.substr(1, c.size() - 2));
                child.key = args[i];
                node.children.emplace_back(std::move(child));
            }

            // "x" argument goes last
            auto it = std::find_if(node.children.begin(), node.children.end(),
                                   [&](const Entry &e)
                                   { return e.key == "x"; });

            if (it != node.children.end())
                std::rotate(it, it + 1, node.children.end());

            return node;
        }
    }

    Entry::Entry(const std::string &value, const std::string &type)
        : type(str::toup(type)), value(value) {}

    std::string Entry::str() const
    {
        if (is_empty())
            return "";
        if (is_func_default())
            return func;
        if (is_terminal())
            return value + "|" + type;

        std::string str_args, str_children;
        for (size_t i = 0; i < children.size(); ++i)
        {
            const auto &child = children[i];
            if (i > 0)
            {
                str_args += str::LIST_DELIMITER;
                str_children += str::LIST_DELIMITER;
            }
            str_args += child.key;

            str_children += "[" + child.str() + "]";
        }
        return func + "|" + str_args + "|" + str_children;
    }

    std::string Entry::err_str() const
    {
        if (is_empty())
            return "empty expression";
        if (is_terminal())
        {
            std::string str_type;
            if (is_latent(type))
                str_type = "latent";
            else if (is_data_wide(type))
                str_type = "data";
            else if (is_group_rate(type))
                str_type = "group_rate";
            else if (is_parameter(type))
                str_type = "parameter";
            else // is_numeric
                str_type = "numeric";
            return str_type + " " + msg::val(value);
        }

        return "expression " + msg::code(str());
    }

    bool Entry::is_callable_implicitly(std::string_view point) const noexcept
    {
        if (is_terminal())
            return true;

        bool is_x_arg = false;
        for (const auto &child : children)
        {
            if (!child.is_callable_implicitly(point))
                return false;
            if (child.key == point)
                is_x_arg = true;
        }
        return is_x_arg;
    }

    MAKE_NODE_LIST_PRODUCERS(latent)
    MAKE_NODE_LIST_PRODUCERS(data_wide)
    MAKE_NODE_LIST_PRODUCERS(parameter)

    Definition::Definition(std::string_view expr) { root = parse_expression_(expr); }

    MAKE_PARSER_LIST_PRODUCERS(latent)
    MAKE_PARSER_LIST_PRODUCERS(data_wide)
    MAKE_PARSER_LIST_PRODUCERS(parameter)
}
