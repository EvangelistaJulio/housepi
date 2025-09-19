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

#include <algorithm>
#include <unordered_set>
#include <chrono>

#include "tools.h"

// TOOLS NAMESPACE
namespace tools
{
    std::unordered_set<std::string> _issued_warnings_;

    /// Display a single-line warning message.
    void warning(const std::string &txt)
    {
        if (_issued_warnings_.find(txt) == _issued_warnings_.end())
        {
#ifdef _HOUSEPI_RCPP
            // Rcpp::warning(txt);
            Rf_warningcall(R_NilValue, "%s", txt.c_str());
#else
            std::cerr << "Warning: " + txt << std::endl;
#endif
            _issued_warnings_.insert(txt);
        }
    }

    /// Display a single-line error message.
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

    /////////////////////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////  Auxiliary Functions  ////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////

    /// Internal check for std::string representing NA.
    bool is_na(const std::string &txt)
    {
        if (txt.size() == 0)
            return true;
        else if (txt == HSP_NA_STRING)
            return true;
        else
            return (str_vec({"<NA>", "NA", "NAN"}) % str::toup(txt));
    }

    /// Internal check for double representing NA.
    bool is_na(double x)
    {
        if (std::isnan(x))
            return true;
        return false;
    }

    /// Default seed for pseudo random number generator.
    std::int32_t default_setup_seed()
    {
        auto time_raw = std::chrono::steady_clock::now().time_since_epoch().count();
        // ensure time_now is 64bit
        std::uint64_t time_now = static_cast<std::uint64_t>(time_raw);
        // xor between upper and lower part of time stamp
        std::uint64_t mixed = time_now ^ (time_now >> 32);
        // note: constrain to 32bit for compatibility with R (mask sign bit)
        return static_cast<std::int32_t>(mixed & 0x7FFFFFFF);
    }

    /// Welford online's algorithm for mean and variance.
    /// Ref :
    /// ----
    /// https://en.wikipedia.org/wiki/Algorithms_for_calculating_variance
    std::pair<double, double> compute_stats(const dbl_vec &data)
    {
        size_t i = 0;
        size_t N = data.size();
        double x, delta, delta2;
        double sqr_diff = 0.;

        double avg = 0.;
        while (i < N)
        {
            x = data[i++];
            delta = x - avg;
            avg += delta / i;
            delta2 = x - avg;
            sqr_diff += delta * delta2;
        }

        return {avg, sqr_diff / (N - 1.)};
    }

    //-----------------------------------------------------------------------------------------------//
    //-----------------------------------------------------------------------------------------------//
    //-----------------------------------------------------------------------------------------------//
    //-----------------------------------------------------------------------------------------------//
    //-----------------------------------------------------------------------------------------------//
    //-----------------------------------------------------------------------------------------------//

    /////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////  String Manipulation Functions  /////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////

    // TOOLS::STR NAMESPACE
    namespace str
    {
        /// Convert to lower case.
        std::string tolow(const std::string s)
        {
            std::string out = s;
            std::transform(out.begin(), out.end(), out.begin(), ::tolower);
            return out;
        }

        /// Convert to upper case.
        std::string toup(const std::string s)
        {
            std::string out = s;
            std::transform(out.begin(), out.end(), out.begin(), ::toupper);
            return out;
        }

        /// Convert vector of strings to lower case.
        str_vec tolow(const str_vec &svec)
        {
            str_vec out;
            out.reserve(svec.size());
            for (const auto &s : svec)
                out.push_back(tolow(s));
            return out;
        }

        /// Convert vector of strings to upper case.
        str_vec toup(const str_vec &svec)
        {
            str_vec out;
            out.reserve(svec.size());
            for (const auto &s : svec)
                out.push_back(toup(s));
            return out;
        }

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
                { return a.size() >= std::string(HSP_HIDE_PRFX).size() &&
                         a.compare(0, std::string(HSP_HIDE_PRFX).size(), HSP_HIDE_PRFX) == 0; });
            vec.erase(new_end, vec.end());
            return vec;
        }

        /// Format a dobule value (it deals with NA).
        std::string format(double x)
        {
            return is_na(x) ? HSP_NA_STRING : std::to_string(x);
        }
        /// Format a boolean value.
        std::string format(bool condition)
        {
#ifdef _HOUSEPI_RCPP
            return (condition ? "true" : "false");
#else
            return (condition ? "TRUE" : "FALSE");
#endif
        }
        /// Format a list of string options.
        std::string format(const str_vec &options, size_t max_print)
        {
            if (options.empty())
                return "";
            else
            {
                std::ostringstream oss;
                max_print = std::min(options.size(), max_print);

                oss << "\"" << options[0] << "\"";
                for (size_t i = 1; i < max_print; ++i)
                    oss << ", \"" << options[i] << "\"";
                if (options.size() > max_print)
                    oss << ", ..";
                oss << ".";
                return oss.str();
            }
        }

        /// Format a 2D matrix.
        std::string format(const matrix<double> &vec)
        {
            if (vec.empty())
                return "";
            else
            {
                std::ostringstream oss;
                oss << "[";
                for (size_t i = 0; i < vec.size(); ++i)
                {
                    oss << "[";
                    for (size_t j = 0; j < vec[i].size(); ++j)
                    {
                        oss << std::fixed << std::setprecision(4) << vec[i][j];
                        if (j < vec[i].size() - 1)
                            oss << ", ";
                    }
                    oss << "]";
                    if (i < vec.size() - 1)
                        oss << ",\n";
                }
                oss << "]\n";
                return oss.str();
            }
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
                        return "\nAvailable option is " + format(options, 1);
                    else
                        return "\nAvailable options are: " + format(options, HSP_MAX_STR_PRINT);
                }
                return "";
            }
            else
                return "\nDid you mean \"" + closest + "\"?";
        }

        /////////////////////////////////////////////////////////////////////////////////////////////////
        ////////////////////////////////////  Converter Functions ///////////////////////////////////////
        /////////////////////////////////////////////////////////////////////////////////////////////////

        /// Convert std::string to boolean variable. NA is assumed as FALSE.
        bool tobool(const std::string &tkn)
        {
            if (is_na(tkn))
                return false;
            else if (toup(tkn) == "FALSE" || tkn == "0")
                return false;
            else if (toup(tkn) == "TRUE" || tkn == "1")
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

        /// Converter to double using std::stod with function pointer casting.
        double tod(const std::string &input, const std::string &errtxt)
        {
            if (tools::is_na(input))
                return HSP_NA_DOUBLE;
            return converter_(static_cast<double (*)(const std::string &, std::size_t *)>(std::stod), input, errtxt, nullptr);
        }
    } // TOOLS::STR NAMESPACE

    //-----------------------------------------------------------------------------------------------//
    //-----------------------------------------------------------------------------------------------//
    //-----------------------------------------------------------------------------------------------//
    //-----------------------------------------------------------------------------------------------//
    //-----------------------------------------------------------------------------------------------//
    //-----------------------------------------------------------------------------------------------//

    /// Clean std::string from forbidden characters and convert to UNDEFINED.
    const std::string strcleaner_(std::string tkn)
    {
        // forbidden characters
        for (char r : {'"', '\r', '\n'})
            tkn.erase(remove(tkn.begin(), tkn.end(), r), tkn.end());
        // check for undefined
        if (is_na(tkn))
            tkn = HSP_NA_STRING;
        return tkn;
    }

    /// Read column headers and compare to required headers.
    str_vec read_column_headers_(std::ifstream &input)
    {
        // clean headers
        str_vec headers = {};

        // read first line
        std::string line, tkn;
        getline(input, line);
        std::istringstream iss(line);

        // get headers from first line
        while (getline(iss, tkn, HSP_DELIM))
            headers.push_back(strcleaner_(tkn));

        return headers;
    }

    /// Load files from the row with headers.
    bool read_row_by_headers_(std::ifstream &input, const str_vec &headers, dict_str &this_row)
    {
        std::string line, tkn;

        if (getline(input, line))
        {
            std::istringstream iss(line);
            size_t n_cols = 0;

            while (getline(iss, tkn, HSP_DELIM))
                this_row.emplace(headers[n_cols++], strcleaner_(tkn));

            // last column is empty
            for (size_t i = n_cols; i < headers.size(); i++)
                this_row.emplace(headers[i], HSP_NA_STRING);

            return true;
        }

        return false;
    }

    //-----------------------------------------------------------------------------------------------//
    //-----------------------------------------------------------------------------------------------//
    //-----------------------------------------------------------------------------------------------//
    //-----------------------------------------------------------------------------------------------//
    //-----------------------------------------------------------------------------------------------//
    //-----------------------------------------------------------------------------------------------//

    /// Parse the list from a given token according to delimiter.
    str_vec parse_arg_list(const std::string &tkn, const char delim)
    {
        std::istringstream args_ss(tkn);
        str_vec output = {};
        std::string tmp;

        // rea the entries of the dictionary
        while (getline(args_ss, tmp, delim))
            output.push_back(tmp);

        return output;
    }

    /// Assert the parsed list can be used a list of keys.
    str_vec parse_keys(const std::string &tkn, const char delim)
    {
        str_vec output = parse_arg_list(tkn, delim);
        if (output.empty())
            return {};
        if (output.size() == 1 && tools::is_na(output[0]))
            return {};
        else if (is_any_na(output))
            throw std::invalid_argument("Argument keys cannot be NA.");
        else if (tools::str::unique(output).size() != output.size())
            throw std::invalid_argument("Argument keys must be unique.");
        else
            return output;
    }

    /// Extract headers from deflist.
    str_vec headers_from_deflist(const dict_str &deflist)
    {
        str_vec entries = {};
        for (const auto &pair : deflist)
            entries.push_back(pair.first);
        return entries;
    }

    /// Read csv list store in pathtofile and store it parsed in csv_parsed_list.
    dict_str read_csv_list(const std::string &pathtofile)
    {
        dict_str output = {};

        if (!is_na(pathtofile))
        {
            std::ifstream input(pathtofile.c_str());

            if (input.is_open())
            {
                // Create a std::stringstream of the current line
                for (std::string line; getline(input, line);)
                {
                    std::string name, choice;
                    std::istringstream iss(line);

                    // name
                    if (!getline(iss, name, HSP_DELIM))
                        throw std::runtime_error(
                            "Unable to extract list name. Check file line \"" + line + "\".");
                    // choice
                    if (!getline(iss, choice))
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
    str_vec headers_from_deftable(const dict_str_vec &deftable)
    {
        str_vec entries = {};
        for (const auto &pair : deftable)
            entries.push_back(pair.first);
        return entries;
    }

    /// Read csv table in pathtofile and store it parsed in csv_parsed_table.
    dict_str_vec read_csv_table(const std::string &pathtofile)
    {
        dict_str_vec output = {};

        if (!is_na(pathtofile))
        {
            std::ifstream input(pathtofile.c_str());

            if (input.is_open())
            {
                // read and check for valid headers
                const str_vec headers = read_column_headers_(input);

                dict_str this_row = {};
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

    //-----------------------------------------------------------------------------------------------//
    //-----------------------------------------------------------------------------------------------//
    //-----------------------------------------------------------------------------------------------//
    //-----------------------------------------------------------------------------------------------//
    //-----------------------------------------------------------------------------------------------//
    //-----------------------------------------------------------------------------------------------//

    /// Check for missing entries given a vector of options.
    void check_for_missing_entries(const str_vec &entries,
                                   const str_vec &options,
                                   const std::string &where)
    {
        // proceed only if there's a given vector of options
        if (!options.empty())
        {
            //----- check for missing entries -----//
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
                        "Missing " + where + "s: " + str::format(missing, HSP_MAX_STR_PRINT));
            }
        }
    }

    /// Check for unused entries given a vector of options.
    void check_for_unused_entries(const str_vec &entries,
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
                str::find_closest(unused[0], options, HSP_MIN_STR_DIST));

        else if (unused.size() > 1)
            throw std::invalid_argument(
                "Unknown " + where + "s: " + str::format(unused, HSP_MAX_STR_PRINT) +
                str::find_closest(HSP_NA_STRING, options, 0));
    }

#ifdef _HOUSEPI_RCPP
    /// Format C++ string into Rcpp string.
    Rcpp::String stoRString(const std::string &x)
    {
        if (tools::is_na(x))
            return HSP_NA_STRING;
        return Rcpp::String(x);
    }

    /// Convert R List into map<string, string>.
    dict_str rlist_to_map(SEXP input)
    {
        using namespace Rcpp;
        dict_str output = {};

        if (!Rf_isNull(input))
        {
            List list(input);
            CharacterVector rownames = list.names();

            for (R_xlen_t i = 0; i < list.size(); ++i)
            {
                CharacterVector element = list[i];
                std::string value = as<std::string>(wrap(element[0]));
                output[as<std::string>(rownames[i])] = value;
            }
        }

        return output;
    }

    /// Convert R dataframe into map<string, vector<string>>.
    dict_str_vec rdataframe_to_map(SEXP input)
    {
        using namespace Rcpp;
        dict_str_vec output = {};

        if (!Rf_isNull(input))
        {
            Rcpp::DataFrame df(input);
            CharacterVector colnames = df.names();

            for (R_xlen_t i = 0; i < df.size(); ++i)
            {
                CharacterVector col_str = as<CharacterVector>(df[i]);
                str_vec col_data(col_str.begin(), col_str.end());
                output[as<std::string>(colnames[i])] = std::move(col_data);
            }
        }

        return output;
    }
#endif
} // TOOLS NAMESPACE
