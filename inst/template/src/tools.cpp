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
 * tools.cpp
 *
 */

#include <fstream>
#include <iomanip>
#include <sstream>
#include <algorithm>
#include <unordered_set>
#include <chrono>

#include "tools/types.h"
#include "tools/error.h"
#include "tools/na_utils.h"
#include "tools/operators.h"
#include "tools/string.h"
#include "tools/read.h"

#include "tools/progress.h"

#include "tools/params.h"
#include "tools/sources.h"
#include "tools/parser.h"

#ifdef _HOUSEPI_RCPP
#include <Rcpp.h>
#else
#include <iostream>
#endif

/////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////// *** TOOLS NAMESPACE *** ///////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////
namespace hsp
{
    /////////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////  Stdout Message Function  ////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////

    // TOOLS::MSG NAMESPACE
    namespace msg
    {
        /// Set of already issued warnings.
        std::unordered_set<std::string> housepi_issued_warnings_set;

        /// Throw an error message.
        /// Ref:
        /// -----
        /// Suppressing Call Stack Info in Rcpp-Generated Errors and Warnings
        /// Michael Weylandt — written Mar 27, 2018
        /// https://gallery.rcpp.org/articles/quiet-stop-and-warning/
        void error(const std::string &txt)
        {
#ifdef _HOUSEPI_RCPP
            throw Rcpp::exception(txt.c_str(), false);
#else
            std::cerr << "Error: " + txt + "\nExit..."
                      << std::endl;
            std::exit(EXIT_FAILURE);
#endif
        }

        /// Show a message to stdout.
        void message(const std::string &txt)
        {
#ifdef _HOUSEPI_RCPP
            Rcpp::Rcout << txt;
            Rcpp::Rcout.flush();
#else
            std::cout << txt;
            std::cout.flush();
#endif
        }

        /// Display a warning message if not already displayed before.
        /// Ref:
        /// -----
        /// Suppressing Call Stack Info in Rcpp-Generated Errors and Warnings
        /// Michael Weylandt — written Mar 27, 2018
        /// https://gallery.rcpp.org/articles/quiet-stop-and-warning/
        void warning(const std::string &txt)
        {
            if (!housepi_issued_warnings_set.count(txt))
            {
#ifdef _HOUSEPI_RCPP
                Rf_warningcall(R_NilValue, "%s", txt.c_str());
#else
                std::cerr << "Warning: " << txt << std::endl;
#endif
                housepi_issued_warnings_set.insert(txt);
            }
        }
    }

    /////////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////  String Manipulation Functions  /////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////

    // TOOLS::STR NAMESPACE
    namespace str
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

        /// Convert std::string to boolean variable. NA is assumed as FALSE.
        bool tobool(const std::string &tkn)
        {
            if (na::is(tkn))
                return false;
            else if (tkn == "FALSE" || tkn == "false" || tkn == "False" || tkn == "0")
                return false;
            else if (tkn == "TRUE" || tkn == "true" || tkn == "True" || tkn == "1")
                return true;
            else
                throw std::invalid_argument(
                    "Unrecognized boolean value \"" + tkn + "\".");
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
                std::string sms = (errtxt.empty()) ? "" : "Unable to import " + errtxt + ". ";
                throw std::invalid_argument(sms + "Value \"" + input + "\" is not suitable to conversion.");
            }
        }

        /// Converter to double using std::stod with function pointer casting.
        double tod(const std::string &input, const std::string &errtxt)
        {
            if (na::is(input))
                return na::DOUBLE;
            return converter_(static_cast<double (*)(const std::string &, std::size_t *)>(std::stod), input, errtxt, nullptr);
        }

        /// Converter to int using std::stoi with function pointer casting.
        int toi(const std::string &input, const std::string &errtxt)
        {
            return converter_(static_cast<int (*)(const std::string &, std::size_t *, int)>(std::stoi), input, errtxt, nullptr, 10);
        }

        /// Converter to unsigned long using std::stoul with function pointer casting.
        unsigned long toul(const std::string &input, const std::string &errtxt)
        {
            return converter_(static_cast<unsigned long (*)(const std::string &, std::size_t *, int)>(std::stoul), input, errtxt, nullptr, 10);
        }

        /// Converter to unsigned 32bit integer (for R compatibility).
        std::int32_t toi32(const std::string &input, const std::string &errtxt)
        {
            return static_cast<std::int32_t>(toul(input, errtxt));
        }

        /////////////////////////////////////////////////////////////////////////////////////////////
        //////////////////////////////////  Formatter Functions  ////////////////////////////////////
        /////////////////////////////////////////////////////////////////////////////////////////////

        /// Format a list of string options.
        std::string format(const str_vec &options, size_t max_print)
        {
            if (options.empty())
                return "";
            else
            {
                max_print = std::min(options.size(), max_print);

                std::string output = "\"" + options[0] + "\"";
                for (size_t i = 1; i < max_print; ++i)
                    output += ", \"" + options[i] + "\"";
                if (options.size() > max_print)
                    output += ", ..";
                return output + ".";
            }
        }

        /// Format a 2D matrix.
        std::string format(const matrix<double> &vec)
        {
            if (vec.empty())
                return "";
            else
            {
                std::string output = "[";
                for (size_t i = 0; i < vec.size(); ++i)
                {
                    output += "[";
                    for (size_t j = 0; j < vec[i].size(); ++j)
                    {
                        output += std::to_string(vec[i][j]);
                        if (j < vec[i].size() - 1)
                            output += ", ";
                    }
                    output += "]";
                    if (i < vec.size() - 1)
                        output += ",\n";
                }
                return output + "]\n";
            }
        }

        /////////////////////////////////////////////////////////////////////////////////////////////
        ////////////////////////////////  String Vector Functions  //////////////////////////////////
        /////////////////////////////////////////////////////////////////////////////////////////////

        /// make list composed of distinct entries
        str_vec unique(const str_vec &vec)
        {
            str_vec out;
            out.reserve(vec.size());
            for (const std::string &s : vec)
                if (!(out % s))
                    out.push_back(s);
            return out;
        }

        /// Remove private entries from list.
        str_vec drop_private(str_vec vec)
        {
            auto new_end = std::remove_if(
                vec.begin(), vec.end(),
                [](const std::string &a)
                { return a.size() >= std::string(HIDDEN_PREFIX).size() &&
                         a.compare(0, std::string(HIDDEN_PREFIX).size(), HIDDEN_PREFIX) == 0; });
            vec.erase(new_end, vec.end());
            return vec;
        }

        // Function to calculate Levenshtein distance
        /// Levenshtein distance calculated the dynamic programming approach.
        /// Ref.
        /// -----
        /// https://github.com/guilhermeagostinelli/levenshtein/blob/master/levenshtein.cpp
        unsigned short levenshtein(const std::string &seq_A, const std::string &seq_B)
        {
            unsigned short len_A = seq_A.size();
            unsigned short len_B = seq_B.size();
            // Verification matrix i.e. 2D array which will store the calculated distance
            matrix<unsigned short> verif(len_A + 1, std::vector<unsigned short>(len_B + 1, 0));

            // If one of the words has zero length, the distance is equal to the size of the other word.
            if (len_A == 0)
                return len_B;
            if (len_B == 0)
                return len_A;

            // Sets the first row and the first column of the verification matrix
            for (unsigned short a = 0; a <= len_A; ++a)
                verif[a][0] = a;
            for (unsigned short b = 0; b <= len_B; ++b)
                verif[0][b] = b;

            // Verification step / matrix filling.
            for (unsigned short a = 1; a <= len_A; ++a)
                for (unsigned short b = 1; b <= len_B; ++b)
                    if (seq_A[a - 1] == seq_B[b - 1])
                        verif[a][b] = verif[a - 1][b - 1];
                    else
                        verif[a][b] = std::min({verif[a - 1][b - 1], verif[a - 1][b], verif[a][b - 1]}) + 1;

            // Last position of the matrix contains the distance
            return verif[len_A][len_B];
        }

        /// Function to find the closest match.
        std::string find_closest(const std::string &input,
                                 const str_vec &options,
                                 unsigned short min_dist,
                                 bool show_avail_options)
        {
            std::string closest = "";
            unsigned short dist_closest = min_dist;

            for (const std::string &x : options)
            {
                unsigned short distance = levenshtein(input, x);
                if (distance <= dist_closest)
                {
                    dist_closest = distance;
                    closest = x;
                }
            }

            // prepare a message
            if (options.empty())
                return "";
            else if (dist_closest == 0)
                return "";
            else if (closest.empty())
            {
                if (show_avail_options)
                {
                    if (options.size() == 1)
                        return " Available option is " + format(options, 1);
                    else
                        return "\nAvailable options are: " + format(options);
                }
                return "";
            }
            else
                return " Did you mean \"" + closest + "\"?";
        }

        /// Check for missing entries given a vector of options.
        void assert_no_missing_entries(const str_vec &entries,
                                       const str_vec &options,
                                       const std::string &where)
        {
            if (!options.empty())
            {
                str_vec missing = {};

                for (const std::string &x : options)
                    if (!(entries % x))
                        missing.push_back(x);

                if (!missing.empty())
                {
                    if (missing.size() == 1)
                        throw std::invalid_argument(
                            "Missing " + where + " " + str::format(missing, 1));
                    else
                        throw std::invalid_argument(
                            "Missing " + where + "s: " + str::format(missing));
                }
            }
        }

        /// Check for unused entries given a vector of options.
        void assert_no_unused_entries(const str_vec &entries,
                                      const str_vec &options,
                                      const std::string &where)
        {
            str_vec unused = {};
            for (const std::string &x : entries)
                if (!(options % x))
                    unused.push_back(x);

            if (unused.size() == 1)
                throw std::invalid_argument(
                    "Unknown " + where + " " + str::format(unused, 1) +
                    str::find_closest(unused[0], options));

            else if (unused.size() > 1)
                throw std::invalid_argument(
                    "Unknown " + where + "s: " + str::format(unused) +
                    str::find_closest(hsp::na::STRING, options, 0));
        }
    } // TOOLS::STR NAMESPACE

    /////////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////  CSV Reader Functions  /////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////

    // TOOLS::READ NAMESPACE
    namespace read
    {
        /// Clean std::string from forbidden characters and convert to UNDEFINED.
        const std::string strcleaner_(std::string tkn)
        {
            // forbidden characters
            for (char r : {'"', '\r', '\n'})
                tkn.erase(remove(tkn.begin(), tkn.end(), r), tkn.end());
            // check for undefined
            if (na::is(tkn))
                tkn = hsp::na::STRING;
            return tkn;
        }

        /// Read column headers and compare to required headers.
        str_vec read_column_headers_(std::ifstream &input)
        {
            // clean headers
            str_vec headers = {};

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
        bool read_row_by_headers_(std::ifstream &input, const str_vec &headers, u_dict<std::string> &this_row)
        {
            std::string line, tkn;

            if (std::getline(input, line))
            {
                std::istringstream iss(line);
                size_t n_cols = 0;

                while (std::getline(iss, tkn, DELIMITER))
                    this_row.emplace(headers[n_cols++], strcleaner_(tkn));

                // last column is empty
                for (size_t i = n_cols; i < headers.size(); i++)
                    this_row.emplace(headers[i], hsp::na::STRING);

                return true;
            }

            return false;
        }

        /// Parse the list from the string.
        str_vec parse_list(const std::string &tkn)
        {
            str_vec output = {};
            std::istringstream args_ss(strcleaner_(tkn));
            std::string tmp;

            // rea the entries of the dictionary
            while (std::getline(args_ss, tmp, LIST_DELIMITER))
                output.push_back(tmp);
            return output;
        }

        /// Extract headers from deflist.
        str_vec list_headers(const u_dict<std::string> &deflist)
        {
            str_vec entries;
            entries.reserve(deflist.size());
            for (const auto &pair : deflist)
                entries.push_back(pair.first);
            return entries;
        }

        /// Read csv list store in pathtofile and store it parsed in csv_parsed_list.
        u_dict<std::string> csv_list(const std::string &pathtofile)
        {
            u_dict<std::string> output = {};

            if (!na::is(pathtofile))
            {
                std::ifstream input(pathtofile.c_str());

                if (input.is_open())
                {
                    // Create a std::stringstream of the current line
                    for (std::string line; std::getline(input, line);)
                    {
                        std::string name, choice;
                        std::istringstream iss(line);

                        // name
                        if (!std::getline(iss, name, DELIMITER))
                            throw std::runtime_error(
                                "Unable to extract list name. Check file line \"" + line + "\".");
                        // choice
                        if (!std::getline(iss, choice))
                            throw std::runtime_error(
                                "Unable to extract list value. Check file line \"" + line + "\".");

                        // add the new setup
                        output.emplace(strcleaner_(name), strcleaner_(choice));
                    }
                }
                else
                    throw std::runtime_error(
                        "Unable to open file \"" + pathtofile + "\".");
            }

            return output;
        }

        /// Extract headers from deftable.
        str_vec table_headers(const u_dict<str_vec> &deftable)
        {
            str_vec entries;
            entries.reserve(deftable.size());
            for (const auto &pair : deftable)
                entries.push_back(pair.first);
            return entries;
        }

        /// Read csv table in pathtofile and store it parsed in csv_parsed_table.
        u_dict<str_vec> csv_table(const std::string &pathtofile)
        {
            u_dict<str_vec> output = {};

            if (!na::is(pathtofile))
            {
                std::ifstream input(pathtofile.c_str());

                if (input.is_open())
                {
                    // read and check for valid headers
                    const str_vec headers = read_column_headers_(input);

                    u_dict<std::string> this_row = {};
                    while (read_row_by_headers_(input, headers, this_row))
                    {
                        for (const auto &x : this_row)
                            output[x.first].push_back(x.second);
                        this_row.clear();
                    }
                }
                else
                    throw std::runtime_error(
                        "Unable to open file \"" + pathtofile + "\".");
            }

            return output;
        }
    } // TOOLS::READ NAMESPACE
} // TOOLS NAMESPACE

/////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////// *** PROGRESS NAMESPACE *** //////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////
namespace progress
{
    using namespace std::chrono;

    /////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////  Progress Bar Functions  ///////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////

    /// Format elapsed time in seconds.
    std::string format(std::chrono::seconds scs)
    {
        std::ostringstream oss;

        auto dys = duration_cast<duration<int, std::ratio<86400>>>(scs);
        scs -= duration_cast<seconds>(dys);

        auto hrs = duration_cast<hours>(scs);
        scs -= duration_cast<seconds>(hrs);

        auto mns = duration_cast<minutes>(scs);
        scs -= duration_cast<seconds>(mns);

        // show large time scales only if relevant
        if (dys.count() > 0)
            oss << dys.count() << "d"
                << std::setfill('0') << std::setw(2) << hrs.count() << "h";

        else if (hrs.count() > 0)
            oss << hrs.count() << "h"
                << std::setfill('0') << std::setw(2) << mns.count() << "m";

        else if (mns.count() > 0)
            oss << mns.count() << "m"
                << std::setfill('0') << std::setw(2) << scs.count() << "s";

        else
            oss << scs.count() << "s";

        return oss.str();
    }
    /////////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////// *** Bar : Impl : class *** ///////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////

    /// Implementation class for progress bar.
    ///
    /// -----
    /// Ref :
    /// OTHMAN EL HAMMOUCHI, April 15, 2024
    /// "Fancy thread safe progress bar in R for compiled code"
    /// https://othmanelhammouchi.com/posts/r-parallel-progress/
    class Bar::Impl
    {
    private:
        const std::string symbol_ = u8"\u25A0"; // bar decorator
        unsigned barwidth_ = 35;                // bar total size
        bool verbose_ = true;                   // wheter to show bar
        bool is_bar_print_ = false;             // wheter bar is printing

        unsigned long state_ = 0; // current iteration
        unsigned long total_ = 0; // total number of iterations
        std::string text_ = "";   // current text to show

        steady_clock::time_point start_;
        steady_clock::time_point lastprint_;
        unsigned long thintime_ = 0;

        /// Function to display the progress bar.
        void bar_(const std::string &sms) const { (this->*which_bar_)(sms); }

        void (Impl::*which_bar_)(const std::string &) const = &Impl::do_bar_carriage_return_;
        void do_bar_carriage_return_(const std::string &sms) const { hsp::msg::message("\r" + sms); }
        void do_bar_new_line_(const std::string &sms) const { hsp::msg::message(sms + "\n"); }

    public:
        /////////////////////////////////////////////////////////////////////////////////////////////
        /////////////////////////////////  Bar : Impl : constructor /////////////////////////////////
        /////////////////////////////////////////////////////////////////////////////////////////////

        Impl()
        {
#ifdef _HOUSEPI_RCPP
            Rcpp::checkUserInterrupt();
            // change barwidth according to console width
            barwidth_ = Rcpp::as<unsigned>(Rcpp::Function("getOption")("width")) / 2;

            // change if R interactive environment
            bool is_interactive = Rcpp::as<bool>(Rcpp::Function("interactive")());
            if (is_interactive)
                which_bar_ = &Impl::do_bar_carriage_return_;
            else
                which_bar_ = &Impl::do_bar_new_line_;
#endif
        };
        ~Impl() {};

        /////////////////////////////////////////////////////////////////////////////////////////////
        ////////////////////////////////////  Bar : Impl : setter ///////////////////////////////////
        /////////////////////////////////////////////////////////////////////////////////////////////

        /// Set description text of progress bar.
        void set_text(const std::string &text) { text_ = text; }

        /// Set description text of progress bar.
        void set_thintime(unsigned long thintime)
        {
            verbose_ = thintime > 0;
            thintime_ = thintime - 1;
        }

        /// Start progressbar clock and display firs message.
        void open(unsigned long total)
        {
#ifdef _HOUSEPI_RCPP
            Rcpp::checkUserInterrupt();
#endif
            if (total == 0)
                throw std::invalid_argument(
                    "Progressbar total should be greater than 0.");

            // initialize internal fields
            total_ = total;
            start_ = steady_clock::now();
            lastprint_ = start_;
            state_ = 0;
            is_bar_print_ = false;
        }

        /// Update progress bar given the progress state.
        void update()
        {
#ifdef _HOUSEPI_RCPP
            Rcpp::checkUserInterrupt();
#endif
            // update internal state
            ++state_;

            if (verbose_)
            {
                const double progress = static_cast<double>(state_) / total_;
                steady_clock::time_point now = steady_clock::now();
                seconds elapsed = duration_cast<seconds>(now - start_);
                seconds missing = duration_cast<seconds>(elapsed * (1. / progress - 1.));

                // print update message
                if (static_cast<unsigned long>(duration_cast<seconds>(now - lastprint_).count()) > thintime_)
                {
                    const unsigned where = progress * barwidth_;
                    std::ostringstream oss;
                    oss << text_ << " |";
                    for (unsigned u = 0; u < where; ++u)
                        oss << symbol_;
                    oss << std::string(barwidth_ - where, ' ') << "| "
                        << std::fixed << std::setprecision(0) << 100. * progress << "% "
                        << format(elapsed) << "<" << format(missing)
                        << std::string(10, ' ');
                    bar_(oss.str());

                    lastprint_ = steady_clock::now();
                    is_bar_print_ = true;
                }
                else if (!is_bar_print_) // eventually print initial message
                    if (static_cast<unsigned long>(missing.count()) > thintime_)
                    {
                        std::ostringstream oss;
                        oss << text_ << " |"
                            << std::string(barwidth_, ' ') << "| 0% ";
                        if (which_bar_ == &Impl::do_bar_new_line_)
                            oss << "\n";
                        hsp::msg::message(oss.str());

                        is_bar_print_ = true;
                    }
            }
        }

        /// End the progress bar because of a crash.
        void crash() const
        {
#ifdef _HOUSEPI_RCPP
            Rcpp::checkUserInterrupt();
#endif
            if (is_bar_print_)
                hsp::msg::message("\n");
        }

        /// Close the progress bar.
        void close()
        {
#ifdef _HOUSEPI_RCPP
            Rcpp::checkUserInterrupt();
#endif
            if (verbose_)
            {
                seconds elapsed = duration_cast<seconds>(steady_clock::now() - start_);

                if (static_cast<unsigned long>(elapsed.count()) > thintime_)
                {
                    std::ostringstream oss;
                    oss << text_ << " |";
                    for (unsigned u = 0; u < barwidth_; u++)
                        oss << symbol_;
                    oss << "| " << format(elapsed)
                        << std::string(25, ' ');
                    if (which_bar_ == &Impl::do_bar_carriage_return_)
                        oss << "\n";
                    bar_(oss.str());
                }
            }

            is_bar_print_ = false;
            state_ = total_;
        }
    };

    ///////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////  Bar : Constructor  ///////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////

    Bar::Bar()
        : pImpl_(std::make_unique<Impl>()) {}
    Bar::~Bar() = default;

    ///////////////////////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////  Bar : Setter  /////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////

    void Bar::set_text(const std::string &text) { pImpl_->set_text(text); }
    void Bar::set_thintime(unsigned long thintime) { pImpl_->set_thintime(thintime); }
    void Bar::open(unsigned long total) { pImpl_->open(total); }
    void Bar::update() { pImpl_->update(); }
    void Bar::crash() const { pImpl_->crash(); }
    void Bar::close() { pImpl_->close(); }
} // PROGRESS NAMESPACE

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
        if (c == hsp::read::LIST_DELIMITER && level == 0)
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
    if (hsp::na::is(expr))
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
            child.key = args[i];
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
            child.key = args[i];
            node.children.push_back(std::move(child));
        }
    }

    return node;
}

/// Construct and parse expression into AST (and flattens it).
DefinitionParser::DefinitionParser(const std::string &expr) { root = parse_expression_(expr); }

DefinitionParser::ASTNode::ASTNode(const std::string &value_, const std::string &type_)
{
    value = value_;
    type = hsp::str::toup(type_);
}

/// Rebuild string from node using recursion.
std::string DefinitionParser::ASTNode::str() const
{
    if (is_empty())
        return "";
    else if (is_func_default())
        return func;
    else if (is_terminal())
        return value + "|" + type;

    std::string str_args, str_children;
    for (size_t i = 0; i < children.size(); ++i)
    {
        const ASTNode &child = children[i];
        if (i > 0)
        {
            str_args += hsp::read::DELIMITER;
            str_children += hsp::read::DELIMITER;
        }
        str_args += child.key;

        const std::string cs = child.str();
        if (child.is_braced)
            str_children += "{" + cs + "}";
        else
            str_children + cs;
    }
    return func + "|" + str_args + "|" + str_children;
}

/// Wheter is implicitely callable.
bool DefinitionParser::ASTNode::is_implicit_callable() const
{
    if (!is_terminal())
    {
        bool is_x_arg = false;
        for (const ASTNode &child : children)
        {
            if (!child.is_implicit_callable())
                return false;
            else if (child.key == hsp::params::X)
                is_x_arg = true;
        }
        if (!is_x_arg)
            return false;
    }
    return true;
}

/// Append dependencies of chosen type to stack.
#define MAKE_NODE_LIST_TYPE_DEPEN(type_flag)                                         \
    void DefinitionParser::ASTNode::append_##type_flag##_depen(str_vec &stack) const \
    {                                                                                \
        if (is_terminal() && hsp::sources::is_##type_flag(type))                     \
            stack.push_back(value);                                                  \
        else                                                                         \
            for (const ASTNode &child : children)                                    \
                child.append_##type_flag##_depen(stack);                             \
    }

MAKE_NODE_LIST_TYPE_DEPEN(augmented)
MAKE_NODE_LIST_TYPE_DEPEN(data)
MAKE_NODE_LIST_TYPE_DEPEN(parameter)

/// Return unique list of dependencies of chosen type.
#define MAKE_PARSER_LIST_TYPE_DEPEN(type_flag)                 \
    str_vec DefinitionParser::list_##type_flag##_depen() const \
    {                                                          \
        str_vec stack = {};                                    \
        root.append_##type_flag##_depen(stack);                \
        return hsp::str::unique(stack);                        \
    }

MAKE_PARSER_LIST_TYPE_DEPEN(augmented)
MAKE_PARSER_LIST_TYPE_DEPEN(data)
MAKE_PARSER_LIST_TYPE_DEPEN(parameter)
