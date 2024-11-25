/*
 * tools.cpp
 *
 * Copyright (C) 2024 Francesco Camaglia, Institut Pasteur
 */

#include "tools.h"

// auxiliary namespace for I/O custom functions
namespace tools
{
    /////////////////////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////  Auxiliary Variables  ////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////

    std::unordered_set<std::string> _issued_warnings_;
    constexpr ulong DECOR_LEFT = 3;   // left size of decorator
    constexpr ulong DECOR_WIDTH = 65; // total size of decorator

    /////////////////////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////  Auxiliary Functions  ////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////

    /* An internal check for std::string representing NA. */
    bool is_na(const std::string &txt)
    {
        if (txt.size() == 0)
        {
            return true;
        }
        else
        {
            std::vector<std::string> forbidden = {UNDEFINED, "NA", "nan", "NaN", "NAN"};
            for (const std::string &txtnan : forbidden)
            {
                if (txt == txtnan)
                {
                    return true;
                }
            }
        }

        return false;
    }

    /* A formatted std::string for percentage of num/denom */
    const std::string percentage(double fraction, ushort precision)
    {
        std::ostringstream oss;
        oss << std::fixed << std::setprecision(precision)
            << (100.0 * fraction) << "%";
        return oss.str();
    }

    /* A formatted std::string for TRUE/FALSE values. */
    const std::string truth(bool condition)
    {
#ifdef _HOUSEPI_RCPP
        return (condition ? "true" : "false");
#else
        return (condition ? "TRUE" : "FALSE");
#endif
    }

    /////////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////  RCPP Instructions  ///////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////

    /* Function to change stdout in case of RCPP. */
    void print_stdout(const std::string &txt)
    {
#ifdef _HOUSEPI_RCPP
        Rcpp::Rcout << txt;
        Rcpp::Rcout.flush();
#else
        std::cout << txt;
        std::cout.flush();
#endif
    }

    /* Internal helper function to change stderr in case of RCPP. */
    void print_stderr(const std::string &txt)
    {
#ifdef _HOUSEPI_RCPP
        Rcpp::Rcerr << txt;
        Rcpp::Rcerr.flush();
#else
        std::cerr << txt;
        std::cerr.flush();
#endif
    }

    /////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////  Stdout Print Functions  ///////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////

    /* Display a decorator for a section (if VERBOSE) */
    void print_decorator(const std::string &title, const std::string &txt, char decor, bool verbose)
    {
        using std::string;

        if (verbose)
        {
            string sms = string(DECOR_LEFT, decor) + " " + title + " " + string(DECOR_WIDTH - title.length() - DECOR_LEFT, decor) + "\n";
            if (!is_na(txt))
            {
                sms += txt + "\n";
            }
            print_stdout(sms);
        }
    }

    //-----------------------------------------------------------------------------------------------//
    //-----------------------------------------------------------------------------------------------//
    //-----------------------------------------------------------------------------------------------//
    //-----------------------------------------------------------------------------------------------//
    //-----------------------------------------------------------------------------------------------//
    //-----------------------------------------------------------------------------------------------//

    /////////////////////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////  Warning Functions  //////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////

    /* Raise a warning. */
    // FIXME : Throw and check exceptions instead depending on Rcpp
    void warning(const std::string &txt)
    {
        if (_issued_warnings_.find(txt) == _issued_warnings_.end())
        {
#ifdef _HOUSEPI_RCPP
            Rcpp::warning(txt);
#else
            std::cerr << "\nWarning: " + txt << std::endl;
#endif
            _issued_warnings_.insert(txt);
        }
    }

    /////////////////////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////  Error Functions  ////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////

    /* Display an error message. */
    // FIXME : omp critical is wrong! Throw and check exceptions instead depending on Rcpp
    void error(const std::string &txt)
    {
#ifdef _HOUSEPI_RCPP
        Rcpp::stop(txt);
#else
        std::cerr << "\nError: " + txt + "\nExit...\n"
                  << std::endl;
        std::exit(EXIT_FAILURE);
#endif
    }

    /////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////  Converter Functions ///////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////

    /* Template function to wrap conversion functions with try-catch */
    template <typename Func, typename... Args>
    auto _safe_converter_(Func func, const std::string &input, const std::string &errtxt, Args &&...args) -> decltype(func(input, std::forward<Args>(args)...))
    {
        try
        {
            return func(input, std::forward<Args>(args)...);
        }
        catch (const std::invalid_argument &e)
        {
            throw std::invalid_argument("Impossible to import " + errtxt + ": invalid value '" + input + "'.");
        }
    }

    /* Converter to int using std::stoi with function pointer casting */
    int stoi(const std::string &input, const std::string &errtxt)
    {
        return _safe_converter_(static_cast<int (*)(const std::string &, std::size_t *, int)>(std::stoi), input, errtxt, nullptr, 10);
    }

    /* Converter to double using std::stod with function pointer casting */
    double stod(const std::string &input, const std::string &errtxt)
    {
        return _safe_converter_(static_cast<double (*)(const std::string &, std::size_t *)>(std::stod), input, errtxt, nullptr);
    }

    /* Converter to unsigned long using std::stoul with function pointer casting */
    unsigned long stoul(const std::string &input, const std::string &errtxt)
    {
        return _safe_converter_(static_cast<ulong (*)(const std::string &, std::size_t *, int)>(std::stoul), input, errtxt, nullptr, 10);
    }

    /* Converter to unsigned long long using std::stoull with function pointer casting */
    unsigned long long stoull(const std::string &input, const std::string &errtxt)
    {
        return _safe_converter_(static_cast<unsigned long long (*)(const std::string &, std::size_t *, int)>(std::stoull), input, errtxt, nullptr, 10);
    }

    /* Convert std::string to boolean variable. */
    bool stobool(const std::string &tkn)
    {
        if (tkn == "TRUE" || tkn == "true" || tkn == "True" || tkn == "1")
        {
            return true;
        }
        else if (tkn == "FALSE" || tkn == "false" || tkn == "False" || tkn == "0")
        {
            return false;
        }
        else
        {
            throw std::invalid_argument("Unrecognized boolean value '" + tkn + "'.");
        }
    }

    //-----------------------------------------------------------------------------------------------//
    //-----------------------------------------------------------------------------------------------//
    //-----------------------------------------------------------------------------------------------//
    //-----------------------------------------------------------------------------------------------//
    //-----------------------------------------------------------------------------------------------//
    //-----------------------------------------------------------------------------------------------//

#ifdef _HOUSEPI_RCPP
    void rdataframe_to_map(const Rcpp::DataFrame &df, std::map<std::string, std::vector<std::string>> &output)
    {
        using namespace Rcpp;

        CharacterVector keys = df.names();

        output.clear();
        for (int i = 0; i < df.size(); ++i)
        {
            CharacterVector col_str = as<CharacterVector>(df[i]);
            std::vector<std::string> col_data(col_str.begin(), col_str.end());
            output[as<std::string>(keys[i])] = std::move(col_data);
        }
    }
#endif

    //-----------------------------------------------------------------------------------------------//
    //-----------------------------------------------------------------------------------------------//
    //-----------------------------------------------------------------------------------------------//
    //-----------------------------------------------------------------------------------------------//
    //-----------------------------------------------------------------------------------------------//
    //-----------------------------------------------------------------------------------------------//

    /* Clean std::string from forbidden characters and convert to UNDEFINED. */
    const std::string strcleaner(std::string tkn)
    {
        // forbidden characters
        for (char r : {'"', '\r', '\n'})
        {
            tkn.erase(remove(tkn.begin(), tkn.end(), r), tkn.end());
        }

        // check for undefined
        if (is_na(tkn))
        {
            tkn = UNDEFINED;
        }

        return tkn;
    }

    /* Read column headers and compare to required headers. */
    std::vector<std::string> read_column_headers(std::ifstream &input, const std::vector<std::string> &required_headers)
    {
        // clean headers
        std::vector<std::string> headers = {};

        // read first line
        std::string line, tkn;
        getline(input, line);
        std::istringstream iss(line);

        // get headers from first line
        while (getline(iss, tkn, C_DELIM))
        {
            headers.push_back(strcleaner(tkn));
        }

        // check that all required headers are present
        if (!required_headers.empty())
        {
            bool is_missing = false;
            std::string sms;

            // check the validity of the headers
            for (const std::string &req : required_headers)
            {
                if (std::find(headers.begin(), headers.end(), req) == headers.end())
                {
                    if (!is_missing)
                    {
                        // point to the faulty file
                        sms += "Missing column(s): '" + req + "'";
                        is_missing = true;
                    }
                    else
                    {
                        sms += ", '" + req + "'";
                    }
                }
            }

            // raise the error if something was missing
            if (is_missing)
            {
                throw std::invalid_argument(sms);
            }
        }

        return headers;
    }

    /* Load files from the row with headers */
    bool read_row_by_headers(std::ifstream &input, const std::vector<std::string> &headers, std::map<std::string, std::string> &this_row)
    {
        std::string line, tkn;

        if (getline(input, line))
        {
            std::istringstream iss(line);
            size_t n_cols = 0;

            while (getline(iss, tkn, C_DELIM))
            {
                this_row.emplace(headers[n_cols++], strcleaner(tkn));
            }

            // last column is empty
            for (size_t i = n_cols; i < headers.size(); i++)
            {
                this_row.emplace(headers[i], UNDEFINED);
            }

            return true;
        }
        else
        {
            return false;
        }
    }

    //-----------------------------------------------------------------------------------------------//
    //-----------------------------------------------------------------------------------------------//
    //-----------------------------------------------------------------------------------------------//
    //-----------------------------------------------------------------------------------------------//
    //-----------------------------------------------------------------------------------------------//
    //-----------------------------------------------------------------------------------------------//

    /* Parse the list from the std::string. */
    const std::vector<std::string> parse_list_args(const std::string &tkn, const char delim)
    {
        std::istringstream args_ss(tkn);
        std::vector<std::string> output = {};
        std::string tmp;

        // rea the entries of the dictionary
        while (getline(args_ss, tmp, delim))
        {
            output.push_back(tmp);
        }

        return output;
    }

    /* Read csv table store in pathtofile and store it parsed in csv_parsed_table. */
    void read_csv_table(const std::string &pathtofile,
                        std::map<std::string, std::vector<std::string>> &csv_parsed_table,
                        const std::vector<std::string> &required_headers)
    {
        try
        {
            if (!is_na(pathtofile))
            {
                std::ifstream input(pathtofile.c_str());

                if (input.is_open())
                {
                    std::vector<std::string> headers = read_column_headers(input, required_headers);
                    std::map<std::string, std::string> this_row = {};

                    csv_parsed_table.clear(); // FIXME: is this correct?
                    while (read_row_by_headers(input, headers, this_row))
                    {
                        for (const auto &x : this_row)
                        {
                            csv_parsed_table[x.first].push_back(x.second);
                        }
                        this_row.clear();
                    }
                }
                else
                {
                    throw std::runtime_error(
                        "Impossible to open file '" + pathtofile + "'.");
                }
            }
        }
        catch (const std::exception &e)
        {
            error(e.what());
        }
    }
} // TOOLS NAMESPACE
