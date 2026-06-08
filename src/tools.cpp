/* housepi - Epidemiological Analysis of Household Outbreaks.
 *
 * Copyright © 2024-2026 Institut Pasteur (Paris), CNRS and INSERM.
 * Distributed under the terms of the GNU General Public License (GPLv3).
 * See the LICENSE file for details.
 */

#include <fstream>
#include <sstream>

#include <algorithm>
#include <cmath>
#include <limits>
#include <string>
#include <unordered_set>
#include <vector>

#include "tools/keys.h"
#include "tools/messages.h"
#include "tools/nan.h"
#include "tools/operators.h"
#include "tools/parser.h"
#include "tools/read.h"
#include "tools/string.h"
#include "tools/types.h"

namespace housepi::msg
{
    /// Set of already issued warnings.
    std::unordered_set<std::string> issued_warnings_set_;

    bool is_new_warning(const std::string &txt) noexcept { return !issued_warnings_set_.count(txt); }
    void add_new_warning(const std::string &txt) noexcept { issued_warnings_set_.insert(txt); }
}

/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////  Housepi Model Keys  ///////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

namespace housepi
{
    void KeyHolder::init(const std::vector<std::string> &keys)
    {
        declared_ = keys;
        optional_.clear();
        size_ = declared_.size();
    }

    void KeyHolder::clear_optional()
    {
        optional_.clear();
        size_ = declared_.size();
    }

    void KeyHolder::add(std::string_view name)
    {
        if (count(name))
            return;
        optional_.emplace_back(name);
        ++size_;
    }

    std::vector<std::string> KeyHolder::all() const
    {
        std::vector<std::string> vec;
        vec.reserve(declared_.size() + optional_.size());
        for (const std::string &s : declared_)
            vec.push_back(s);
        for (const std::string &s : optional_)
            vec.push_back(s);
        return vec;
    }

    std::string KeyHolder::as_string(size_t index) const
    {
        if (index < declared_.size())
            return declared_[index];
        index -= declared_.size();
        if (index < optional_.size())
            return optional_[index];
        throw std::out_of_range("Registration failed: key index out of range.");
    }
    std::vector<std::string> KeyHolder::as_string(std::vector<size_t> indexes) const
    {
        std::vector<std::string> output;
        output.reserve(indexes.size());
        for (size_t s : indexes)
            output.push_back(as_string(s));
        return output;
    }

    size_t KeyHolder::from_string(std::string_view name) const
    {
        auto it = std::find(declared_.begin(), declared_.end(), name);
        if (it != declared_.end())
            return static_cast<size_t>(std::distance(declared_.begin(), it));
        auto it2 = std::find(optional_.begin(), optional_.end(), name);
        if (it2 != optional_.end())
            return declared_.size() + static_cast<size_t>(std::distance(optional_.begin(), it2));
        throw std::out_of_range("Registration failed: key name '" + std::string(name) + "' was not found.");
    }
    std::vector<size_t> KeyHolder::from_string(const std::vector<std::string> &names) const
    {
        std::vector<size_t> output;
        output.reserve(names.size());
        for (const std::string &s : names)
            output.push_back(from_string(s));
        return output;
    }

    bool KeyHolder::count(std::string_view name) const
    {
        return std::find(declared_.begin(), declared_.end(), name) != declared_.end() ||
               std::find(optional_.begin(), optional_.end(), name) != optional_.end();
    }
}

/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////  Not Available Recognition Functions  //////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

namespace housepi::na
{
    bool is(double x) noexcept { return std::isnan(x); }
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

/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////  String Manipulation Functions  /////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

#ifndef M_PI
#define M_PI 3.14159265358979323846264338327950288 /// pi
#endif

namespace housepi::str
{
    /////////////////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////  Converter Functions  ////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////

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

#ifndef NDEBUG
        throw std::invalid_argument(
            "Unrecognized boolean value '" + tkn + "'.");
#endif
        return false;
    }

    /// Template function to wrap conversion functions with try-catch.
    /// Using function pointer casting.
    template <typename Func, typename... Args>
    auto helper_converter_(Func func, const std::string &input, const std::string &errtxt, Args &&...args) -> decltype(func(input, std::forward<Args>(args)...))
    {
        try
        {
            return func(input, std::forward<Args>(args)...);
        }
        catch (const std::invalid_argument &e)
        {
            throw std::invalid_argument(
                (errtxt.empty() ? "" : "Unable to import " + errtxt + ". ") +
                "Value '" + input + "' is not suitable to conversion.");
        }
    }

    double tod(const std::string &input, const std::string &errtxt)
    {
        if (na::is(input))
            return na::DOUBLE;
        // UPDATE: this is not a robust approach...
        if (input == "pi")
            return M_PI;
        return helper_converter_(static_cast<double (*)(const std::string &, size_t *)>(std::stod), input, errtxt, nullptr);
    }
    int toi(const std::string &input, const std::string &errtxt)
    {
        return helper_converter_(static_cast<int (*)(const std::string &, size_t *, int)>(std::stoi), input, errtxt, nullptr, 10);
    }
    unsigned long toul(const std::string &input, const std::string &errtxt)
    {
        return helper_converter_(static_cast<unsigned long (*)(const std::string &, size_t *, int)>(std::stoul), input, errtxt, nullptr, 10);
    }

    /////////////////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////  Formatter Functions  ////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////

    std::string format(bool condition) { return condition ? "TRUE" : "FALSE"; }

    std::string format(double x) { return na::is(x) ? "NA" : std::to_string(x); }

    std::string format(const std::vector<std::string> &options,
                       size_t max_print, bool add_dot)
    {
        if (options.empty())
            return "";

        max_print = std::min(options.size(), max_print);
        std::string output = "'" + options[0] + "'";
        for (size_t i = 1; i < max_print; ++i)
            output += ", '" + options[i] + "'";
        if (options.size() > max_print)
            output += add_dot ? ", .." : ", ...";
        else if (add_dot)
            output += ".";
        return output;
    }

    /////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////  String Vector Functions  //////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////

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

    std::vector<std::string> drop_private(std::vector<std::string> vec)
    {
        auto new_end = std::remove_if(vec.begin(), vec.end(), is_private);
        vec.erase(new_end, vec.end());
        return vec;
    }

    /// Calcuate the Levenshtein distance.
    /// Distance is equal to seq_B length if seq_A is empty, and viceversa.
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

    std::string find_closest(std::string_view input,
                             const std::vector<std::string> &options,
                             size_t min_dist,
                             bool show_avail_options)
    {
        std::string closest = "";
        size_t dist_closest = min_dist;

        for (const std::string &x : options)
        {
            size_t d = helper_levenshtein_(input, x);
            if (d <= dist_closest)
            {
                dist_closest = d;
                closest = x;
            }
        }

        // prepare a message
        if (options.empty())
            return "";
        if (dist_closest == 0)
            return "";
        if (!closest.empty())
            return " Did you mean '" + closest + "'?";
        if (!show_avail_options)
            return "";
        if (options.size() == 1)
            return " Available option is " + format(options, 1);
        return "\nAvailable options are: " + format(options);
    }
}

/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////  CSV Reader Functions  /////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

namespace housepi::read
{
    std::string strcleaner(std::string tkn)
    {
        // forbidden characters
        for (char r : {'"', '\r', '\n'})
            tkn.erase(remove(tkn.begin(), tkn.end(), r), tkn.end());
        return tkn;
    }

    std::vector<std::string> parse_list(const std::string &tkn)
    {
        std::vector<std::string> output = {};
        std::istringstream args_ss(strcleaner(tkn));
        std::string tmp;

        // rea the entries of the dictionary
        while (std::getline(args_ss, tmp, LIST_DELIMITER))
            output.push_back(tmp);
        return output;
    }

    // UPDATE: can't I merge this and the previous 'parse_list'?
    std::vector<std::string> helper_parse_names(std::string input)
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

/////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////  Dependencies Checking ///////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////

namespace housepi::dependency
{
    void walk_tree(const std::string &entry,
                   const str_table &tree,
                   std::vector<std::string> &visited,
                   std::vector<std::string> &recursion_stack,
                   std::vector<std::string> &ordered)
    {
        if (!tree.count(entry))
            throw std::invalid_argument(
                "Could not find a definition for '" + entry + "'");

        if (std::find(recursion_stack.begin(), recursion_stack.end(), entry) != recursion_stack.end())
            throw std::invalid_argument(
                "Found inadmissible cyclic definition for '" + entry + "'");

        if (std::find(visited.begin(), visited.end(), entry) != visited.end())
            return;

        recursion_stack.push_back(entry);
        for (const std::string &dep : tree.at(entry))
            walk_tree(dep, tree, visited, recursion_stack, ordered);

        recursion_stack = recursion_stack - std::vector<std::string>({entry});
        visited.push_back(entry);
        ordered.push_back(entry);
    }

    std::vector<std::string> order_dependencies(const str_table &graph)
    {
        std::vector<std::string> visited, recursion_stack, ordered;
        for (const auto &[a, vec] : graph)
            if (std::find(visited.begin(), visited.end(), a) == visited.end())
                walk_tree(a, graph, visited, recursion_stack, ordered);
        return ordered;
    }

    std::vector<std::string> order_conditioned_on(const std::string &start, const str_table &graph)
    {
        // construct a reversed graph for dependents
        str_table reverse_graph;
        for (const auto &[node, dependents] : graph)
        {
            for (const std::string &dep : dependents)
                reverse_graph[dep].push_back(node);

            // ensure isolated nodes are present
            if (!reverse_graph.count(node))
                reverse_graph[node] = {};
        }

        // walk dependents of start
        std::vector<std::string> visited, recursion_stack, ordered;
        for (const std::string &dep : reverse_graph.at(start))
            walk_tree(dep, reverse_graph, visited, recursion_stack, ordered);

        std::reverse(ordered.begin(), ordered.end());
        return std::vector<std::string>({start}) + ordered;
    }
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////  Definition Parser /////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////

/// Append dependencies of chosen type to stack.
#define MAKE_NODE_LIST_TYPE_DEPEN(type_flag)                                      \
    void Entry::append_##type_flag##_depen(std::vector<std::string> &stack) const \
    {                                                                             \
        if (is_terminal() && is_##type_flag(type))                                \
            stack.push_back(value);                                               \
        else                                                                      \
            for (const auto &child : children)                                    \
                child.append_##type_flag##_depen(stack);                          \
    }

/// Return unique list of dependencies of chosen type.
#define MAKE_PARSER_LIST_TYPE_DEPEN(type_flag)                            \
    std::vector<std::string> Definition::list_##type_flag##_depen() const \
    {                                                                     \
        std::vector<std::string> stack = {};                              \
        root.append_##type_flag##_depen(stack);                           \
        return str::unique(stack);                                        \
    }

namespace housepi::parse
{
    /// Split on sep at top level (ignore separators inside braces).
    std::vector<std::string> split_top_(const std::string &s)
    {
        std::vector<std::string> entries;
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
            if (c == read::LIST_DELIMITER && level == 0)
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

        // if not found, then it's a terminal node {value|type}
        size_t bar_2 = expr_tail.find('|');
        if (bar_2 == std::string::npos)
            return Entry(expr_head, expr_tail);

        // otherwise it is a non-default function node
        // e. g. "func|arg_1:...:arg_n|child_1:...child_n"
        const std::vector<std::string> args = split_top_(expr_tail.substr(0, bar_2));
        const std::vector<std::string> children = split_top_(expr_tail.substr(bar_2 + 1));

        if (args.size() != children.size())
            throw std::runtime_error(
                "Mismatch in the number of keys and arguments of definition '" + expr + "'.");

        // populate the parent node with children nodes
        Entry node;
        node.func = expr_head;
        node.children.reserve(children.size());
        for (size_t i = 0; i < children.size(); ++i)
        {
            const std::string &c = children[i];

            if (c.size() >= 2 && c.front() == '{' && c.back() == '}')
            {
                // function child node
                // parse inner expression and keep braced flag
                auto child = parse_expression_(c.substr(1, c.size() - 2));
                child.is_braced = true;
                child.key = args[i];
                node.children.emplace_back(child);
            }
            else
            {
                // terminal child node
                size_t barc = c.find('|');
                if (barc == std::string::npos)
                    throw std::runtime_error(
                        "Missing terminal '|' in definition '" + expr + "'.");
                Entry child(c.substr(0, barc), c.substr(barc + 1));
                child.key = args[i];
                node.children.emplace_back(child);
            }
        }

        // ensure "x" arguments is defined at the end
        auto it = std::find_if(node.children.begin(), node.children.end(),
                               [&](const Entry &e)
                               { return e.key == "x"; });

        if (it != node.children.end())
            std::rotate(it, it + 1, node.children.end());

        return node;
    }

    Entry::Entry(const std::string &value_, const std::string &type_)
    {
        value = value_;
        type = str::toup(type_);
    }

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
                str_args += read::DELIMITER;
                str_children += read::DELIMITER;
            }
            str_args += child.key;

            const std::string &cs = child.str();
            if (child.is_braced)
                str_children += "{" + cs + "}";
            else
                str_children += cs;
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
            if (is_augmented(type))
                str_type = "augmented";
            else if (is_data_wide(type))
                str_type = "data";
            else if (is_group_rate(type))
                str_type = "group_rate";
            else if (is_parameter(type))
                str_type = "parameter";
            else // is_numeric
                str_type = "numeric";
            return str_type + " '" + value + "'";
        }

        return "expression '" + str() + "'";
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

    MAKE_NODE_LIST_TYPE_DEPEN(augmented)
    MAKE_NODE_LIST_TYPE_DEPEN(data_wide)
    MAKE_NODE_LIST_TYPE_DEPEN(parameter)

    Definition::Definition(std::string_view expr) { root = parse_expression_(expr); }

    MAKE_PARSER_LIST_TYPE_DEPEN(augmented)
    MAKE_PARSER_LIST_TYPE_DEPEN(data_wide)
    MAKE_PARSER_LIST_TYPE_DEPEN(parameter)
}
