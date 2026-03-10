/* housepi - Epidemiological Analysis of Household Outbreaks.
 *
 * Copyright © 2026 Institut Pasteur (Paris), CNRS and INSERM.
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

namespace hsp::msg
{
    /// Set of already issued warnings.
    std::unordered_set<std::string> housepi_issued_warnings_set;

    /// Wheter text represent a previously unseen warning.
    bool is_new_housepi_warning(const std::string &txt) noexcept { return !housepi_issued_warnings_set.count(txt); }
    /// Mark given text as an already seen warning.
    void add_new_housepi_warning(const std::string &txt) noexcept { housepi_issued_warnings_set.insert(txt); }
}

/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////  Housepi Model Keys  ///////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

namespace hsp::key
{
    std::vector<std::string> Registry::get_all() const
    {
        std::vector<std::string> vec;
        vec.reserve(declared_.size() + optional_.size());
        for (const std::string &s : declared_)
            vec.push_back(s);
        for (const std::string &s : optional_)
            vec.push_back(s);
        return vec;
    }

    std::string Registry::as_string(size_t index) const
    {
        if (index < declared_.size())
            return declared_[index];
        index -= declared_.size();
        if (index < optional_.size())
            return optional_[index];
        throw std::out_of_range("Registration failed: key index out of range.");
    }
    std::vector<std::string> Registry::as_string(std::vector<size_t> indexes)
    {
        std::vector<std::string> output;
        output.reserve(indexes.size());
        for (size_t s : indexes)
            output.push_back(as_string(s));
        return output;
    }

    size_t Registry::from_string(std::string_view name) const
    {
        auto it = std::find(declared_.begin(), declared_.end(), name);
        if (it != declared_.end())
            return static_cast<size_t>(std::distance(declared_.begin(), it));
        auto it2 = std::find(optional_.begin(), optional_.end(), name);
        if (it2 != optional_.end())
            return declared_.size() + static_cast<size_t>(std::distance(optional_.begin(), it2));
        throw std::out_of_range("Registration failed: key name '" + std::string(name) + "' was not found.");
    }
    std::vector<size_t> Registry::from_string(const std::vector<std::string> &names) const
    {
        std::vector<size_t> output;
        output.reserve(names.size());
        for (const std::string &s : names)
            output.push_back(from_string(s));
        return output;
    }

    bool Registry::count(std::string_view name) const
    {
        return std::find(declared_.begin(), declared_.end(), name) != declared_.end() ||
               std::find(optional_.begin(), optional_.end(), name) != optional_.end();
    }

    /// Parse comma-separated enum names from stringified VA_ARGS.
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

/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////  Not Available Recognition Functions  //////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

namespace hsp::na
{
    /// Wheter double represents NA.
    bool is(double x) noexcept { return std::isnan(x); }
    /// Wheter string represents NA.
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

namespace hsp::str
{
    /////////////////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////  Converter Functions  ////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////

    /// Convert to lower case.
    std::string tolow(std::string s)
    {
        std::transform(s.begin(), s.end(), s.begin(),
                       [](unsigned char c)
                       { return std::tolower(c); });
        return s;
    }

    /// Convert to upper case.
    std::string toup(std::string s)
    {
        std::transform(s.begin(), s.end(), s.begin(),
                       [](unsigned char c)
                       { return std::toupper(c); });
        return s;
    }

    /// Convert string to boolean variable. NA is assumed as FALSE.
    bool tobool(const std::string &tkn)
    {
        if (na::is(tkn))
            return false;
        if (tkn == "FALSE" || tkn == "false" || tkn == "False" || tkn == "0")
            return false;
        if (tkn == "TRUE" || tkn == "true" || tkn == "True" || tkn == "1")
            return true;

        throw std::invalid_argument(
            "Unrecognized boolean value '" + tkn + "'.");
    }

    /// Template function to wrap conversion functions with try-catch.
    template <typename Func, typename... Args>
    auto converter_(Func func, const std::string &input, const std::string &errtxt, Args &&...args) -> decltype(func(input, std::forward<Args>(args)...))
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

    /// Converter to double using stod with function pointer casting.
    double tod(const std::string &input, const std::string &errtxt)
    {
        if (na::is(input))
            return na::DOUBLE;
        return converter_(static_cast<double (*)(const std::string &, size_t *)>(std::stod), input, errtxt, nullptr);
    }

    /// Converter to int using stoi with function pointer casting.
    int toi(const std::string &input, const std::string &errtxt)
    {
        return converter_(static_cast<int (*)(const std::string &, size_t *, int)>(std::stoi), input, errtxt, nullptr, 10);
    }

    /// Converter to unsigned long using stoul with function pointer casting.
    unsigned long toul(const std::string &input, const std::string &errtxt)
    {
        return converter_(static_cast<unsigned long (*)(const std::string &, size_t *, int)>(std::stoul), input, errtxt, nullptr, 10);
    }

    /// Converter to unsigned 32bit integer (for R compatibility).
    std::int32_t toi32(const std::string &input, const std::string &errtxt)
    {
        return static_cast<std::int32_t>(toul(input, errtxt));
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

    /// Wheter a name has been tagged as private.
    bool is_private(std::string_view s) noexcept
    {
        return s.size() >= HIDDEN_PREFIX.size() &&
               s.compare(0, HIDDEN_PREFIX.size(), HIDDEN_PREFIX) == 0;
    }

    /// Make list composed of distinct entries.
    std::vector<std::string> unique(const std::vector<std::string> &vec)
    {
        std::vector<std::string> out;
        out.reserve(vec.size());
        for (const std::string &s : vec)
            if (std::find(out.begin(), out.end(), s) == out.end())
                out.push_back(s);
        return out;
    }

    /// Drop private entries from list.
    std::vector<std::string> drop_private(std::vector<std::string> vec)
    {
        auto new_end = std::remove_if(vec.begin(), vec.end(), is_private);
        vec.erase(new_end, vec.end());
        return vec;
    }

    /// Calcuate the Levenshtein distance.
    /// Distance is equal to seq_B length if seq_A is empty, and viceversa.
    /// (modified from) Ref.
    /// -----
    /// Hjelmqvist, Sten (26 March 2012), Fast, memory efficient Levenshtein algorithm.
    size_t levenshtein(std::string_view seq_A, std::string_view seq_B)
    {
        size_t len_A = seq_A.size();
        size_t len_B = seq_B.size();

        if (len_A == 0)
            return len_B;
        if (len_B == 0)
            return len_A;

        std::vector<size_t> prev(len_B + 1), curr(len_B + 1);

        for (size_t b = 0; b <= len_B; ++b)
            prev[b] = b;

        for (size_t a = 1; a <= len_A; ++a)
        {
            curr[0] = a;
            for (size_t b = 1; b <= len_B; ++b)
            {
                size_t cost = (seq_A[a - 1] == seq_B[b - 1]) ? 0 : 1;
                curr[b] = std::min(prev[b] + 1,
                                   std::min(curr[b - 1] + 1,
                                            prev[b - 1] + cost));
            }
            std::swap(prev, curr);
        }

        return prev[len_B];
    }

    /// Function to find the closest match.
    std::string find_closest(std::string_view input,
                             const std::vector<std::string> &options,
                             size_t min_dist,
                             bool show_avail_options)
    {
        std::string closest = "";
        size_t dist_closest = min_dist;

        for (const std::string &x : options)
        {
            size_t d = levenshtein(input, x);
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

    /// Check for missing entries given a vector of options.
    void assert_no_missing_entries(const std::vector<std::string> &entries,
                                   const std::vector<std::string> &options,
                                   const std::string &where)
    {
        if (options.empty())
            return;

        std::vector<std::string> missing = {};

        for (const std::string &x : options)
            if (std::find(entries.begin(), entries.end(), x) == entries.end())
                missing.push_back(x);

        if (missing.size() == 1)
            throw std::invalid_argument(
                "Missing " + where + " " + str::format(missing, 1));

        if (!missing.empty())
            throw std::invalid_argument(
                "Missing " + where + "s: " + str::format(missing));
    }

    /// Check for unused entries given a vector of options.
    void assert_no_unused_entries(const std::vector<std::string> &entries,
                                  const std::vector<std::string> &options,
                                  const std::string &where)
    {
        std::vector<std::string> unused = {};
        for (const std::string &x : entries)
            if (std::find(options.begin(), options.end(), x) == options.end())
                unused.push_back(x);

        if (unused.size() == 1)
            throw std::invalid_argument(
                "Unknown " + where + " " + str::format(unused, 1) +
                str::find_closest(unused[0], options));

        if (!unused.empty())
            throw std::invalid_argument(
                "Unknown " + where + "s: " + str::format(unused) +
                str::find_closest("", options, 0));
    }
}

/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////  CSV Reader Functions  /////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

namespace hsp::read
{
    /// Clean string from forbidden characters and convert to UNDEFINED.
    const std::string strcleaner_(std::string tkn)
    {
        // forbidden characters
        for (char r : {'"', '\r', '\n'})
            tkn.erase(remove(tkn.begin(), tkn.end(), r), tkn.end());
        // check for undefined
        if (na::is(tkn))
            tkn = "";
        return tkn;
    }

    /// Read column headers and compare to required headers.
    std::vector<std::string> read_column_headers_(std::ifstream &input)
    {
        // clean headers
        std::vector<std::string> headers = {};

        // read first line
        std::string line, tkn;
        std::getline(input, line);
        std::istringstream iss(line);

        // get headers from first line
        while (std::getline(iss, tkn, DELIMITER))
            headers.push_back(strcleaner_(tkn));

        return headers;
    }

    /// Load files from the row with headers.
    bool read_row_by_headers_(std::ifstream &input,
                              const std::vector<std::string> &headers,
                              dict<std::string> &this_row)
    {
        std::string line, tkn;

        if (!std::getline(input, line))
            return false;

        std::istringstream iss(line);
        size_t n_cols = 0;

        while (std::getline(iss, tkn, DELIMITER))
            this_row.emplace(headers[n_cols++], strcleaner_(tkn));

        // last column is empty
        for (size_t i = n_cols; i < headers.size(); i++)
            this_row.emplace(headers[i], "");

        return true;
    }

    /// Parse the list from the string.
    std::vector<std::string> parse_list(const std::string &tkn)
    {
        std::vector<std::string> output = {};
        std::istringstream args_ss(strcleaner_(tkn));
        std::string tmp;

        // rea the entries of the dictionary
        while (std::getline(args_ss, tmp, LIST_DELIMITER))
            output.push_back(tmp);
        return output;
    }

    /// Read csv table in pathtofile and store it parsed in csv_parsed_table.
    str_table csv_table(const std::string &pathtofile)
    {
        if (na::is(pathtofile))
            return {};

        std::ifstream input(pathtofile.c_str());
        if (!input.is_open())
            throw std::runtime_error(
                "Unable to open file '" + pathtofile + "'.");

        // read and check for valid headers
        const std::vector<std::string> headers = read_column_headers_(input);

        str_table output = {};
        dict<std::string> this_row = {};
        while (read_row_by_headers_(input, headers, this_row))
        {
            for (const auto &[key, val] : this_row)
                output[key].push_back(val);
            this_row.clear();
        }

        // ensure all columns have the same number of rows
        size_t n_rows = csv_n_rows(output);
        for (const auto &[col, vec] : output)
            if (vec.size() != n_rows)
                throw std::runtime_error(
                    "Inconsistent number of rows found in column '" + col +
                    "' while reading file '" + pathtofile + "'.");

        return output;
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

namespace hsp::parse
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

    /// Definition Entry constructor.
    Entry::Entry(const std::string &value_, const std::string &type_)
    {
        value = value_;
        type = str::toup(type_);
    }

    /// Rebuild string from node using recursion.
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

    /// Wheter is implicitely callable.
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

    /// Construct and parse expression into AST (and flattens it).
    Definition::Definition(std::string_view expr) { root = parse_expression_(expr); }

    MAKE_PARSER_LIST_TYPE_DEPEN(augmented)
    MAKE_PARSER_LIST_TYPE_DEPEN(data_wide)
    MAKE_PARSER_LIST_TYPE_DEPEN(parameter)
}
