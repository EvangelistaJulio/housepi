/*
 * tools.cpp
 *
 * Copyright (C) 2025 Francesco Camaglia, Institut Pasteur
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
        str_vec forbidden = {UNDEFINED, "NA", "nan", "NaN", "NAN"};

        if (txt.size() == 0)
            return true;
        else
            for (const std::string &txtnan : forbidden)
                if (txt == txtnan)
                    return true;

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
                sms += txt + "\n";

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
    void error(const std::string &txt)
    {
#ifdef _HOUSEPI_RCPP
        print_stderr("\n");
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
            throw std::invalid_argument(
                "Impossible to import " + errtxt + ". Value '" + input + "' is not suitabe to conversion.");
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

    /* Converter checking for NA values formatter for NA */
    double natod(const std::string &input)
    {
        if (is_na(input))
            return std::numeric_limits<double>::quiet_NaN();
        else
            return stod(input, "");
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
            return true;

        else if (tkn == "FALSE" || tkn == "false" || tkn == "False" || tkn == "0")
            return false;

        else
            throw std::invalid_argument(
                "Unrecognized boolean value '" + tkn + "'.");
    }    

    //-----------------------------------------------------------------------------------------------//
    //-----------------------------------------------------------------------------------------------//
    //-----------------------------------------------------------------------------------------------//
    //-----------------------------------------------------------------------------------------------//
    //-----------------------------------------------------------------------------------------------//
    //-----------------------------------------------------------------------------------------------//

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
            tkn.erase(remove(tkn.begin(), tkn.end(), r), tkn.end());

        // check for undefined
        if (is_na(tkn))
            tkn = UNDEFINED;

        return tkn;
    }

    /* Read column headers and compare to required headers. */
    str_vec read_column_headers(std::ifstream &input, const str_vec &required_headers)
    {
        // clean headers
        str_vec headers = {};

        // read first line
        std::string line, tkn;
        getline(input, line);
        std::istringstream iss(line);

        // get headers from first line
        while (getline(iss, tkn, C_DELIM))
            headers.push_back(strcleaner(tkn));

        // check that all required headers are present
        if (!required_headers.empty())
        {
            std::string sms = "";

            // check the validity of the headers
            for (const std::string &req : required_headers)
                if (std::find(headers.begin(), headers.end(), req) == headers.end())
                    sms += ((sms.empty()) ? "Missing column(s)" : ",") + (" '" + req + "'");

            // raise the error if something was missing
            if (!sms.empty())
                throw std::invalid_argument(sms + ".");
        }

        return headers;
    }

    /* Load files from the row with headers */
    bool read_row_by_headers(std::ifstream &input, const str_vec &headers, std::map<std::string, std::string> &this_row)
    {
        std::string line, tkn;

        if (getline(input, line))
        {
            std::istringstream iss(line);
            size_t n_cols = 0;

            while (getline(iss, tkn, C_DELIM))
                this_row.emplace(headers[n_cols++], strcleaner(tkn));

            // last column is empty
            for (size_t i = n_cols; i < headers.size(); i++)
                this_row.emplace(headers[i], UNDEFINED);

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
    const str_vec parse_list_args(const std::string &tkn, const char delim)
    {
        std::istringstream args_ss(tkn);
        str_vec output = {};
        std::string tmp;

        // rea the entries of the dictionary
        while (getline(args_ss, tmp, delim))
            output.push_back(tmp);

        return output;
    }

    /* Read csv list store in pathtofile and store it parsed in csv_parsed_list. */
    void read_csv_list(const std::string &pathtofile, std::map<std::string, std::string> &csv_parsed_list)
    {
        try
        {
            if (!is_na(pathtofile))
            {
                std::ifstream input(pathtofile.c_str());

                if (input.is_open())
                {
                    csv_parsed_list.clear();

                    // Create a std::stringstream of the current line
                    for (std::string line; getline(input, line);)
                    {
                        std::string name, choice;
                        std::istringstream iss(line);

                        //===== NAME =====//
                        if (!getline(iss, name, C_DELIM))
                            throw std::runtime_error(
                                "Unable to extract list name. Check file line '" + line + "'.");

                        //===== CHOICE =====//
                        if (!getline(iss, choice))
                            throw std::runtime_error(
                                "Unable to extract list value. Check file line '" + line + "'.");

                        // add the new setup
                        csv_parsed_list.emplace(strcleaner(name), strcleaner(choice));
                    }
                }
                else
                {
                    throw std::runtime_error(
                        "Impossible to open file '" + pathtofile + "'.");
                }
            }
        }
        catch (std::exception &e)
        {
            error(e.what());
        }
    }

    /* Read csv table store in pathtofile and store it parsed in csv_parsed_table. */
    void read_csv_table(const std::string &pathtofile,
                        dict_str_vec &csv_parsed_table,
                        const str_vec &required_headers)
    {
        if (!is_na(pathtofile))
        {
            std::ifstream input(pathtofile.c_str());

            if (input.is_open())
            {
                str_vec headers = read_column_headers(input, required_headers);
                std::map<std::string, std::string> this_row = {};

                csv_parsed_table.clear(); // FIXME: is this correct?
                while (read_row_by_headers(input, headers, this_row))
                {
                    for (const auto &x : this_row)
                        csv_parsed_table[x.first].push_back(x.second);
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

    //-----------------------------------------------------------------------------------------------//
    //-----------------------------------------------------------------------------------------------//
    //-----------------------------------------------------------------------------------------------//
    //-----------------------------------------------------------------------------------------------//
    //-----------------------------------------------------------------------------------------------//
    //-----------------------------------------------------------------------------------------------//

    /* Welford online algorithm for mean and variance. */
    void compute_stats(const std::vector<double> &data, double &avg, double &var)
    {
        size_t N = data.size();
        double x, delta, delta2;
        double sqr_diff = 0.0;

        avg = 0.0;
        for (size_t i = 0; i < N; ++i)
        {
            x = data[i];
            delta = x - avg;
            avg += delta / (i + 1);
            delta2 = x - avg;
            sqr_diff += delta * delta2;
        }
        var = sqr_diff / (N - 1);
    }

    //-----------------------------------------------------------------------------------------------//
    //-----------------------------------------------------------------------------------------------//
    //-----------------------------------------------------------------------------------------------//
    //-----------------------------------------------------------------------------------------------//
    //-----------------------------------------------------------------------------------------------//
    //-----------------------------------------------------------------------------------------------//

    // Function to calculate Levenshtein distance
    /* Levenshtein distance calculated the dynamic programming approach.
        Ref.
        -----
        https://github.com/guilhermeagostinelli/levenshtein/blob/master/levenshtein.cpp
    */
    ushort levenshtein(const std::string &seq_A, const std::string &seq_B)
    {
        ushort len_A = seq_A.size();
        ushort len_B = seq_B.size();
        // Verification matrix i.e. 2D array which will store the calculated distance
        matrix<ushort> verif(len_A + 1, std::vector<ushort>(len_B + 1, 0));

        // If one of the words has zero length, the distance is equal to the size of the other word.
        if (len_A == 0)
            return len_B;
        if (len_B == 0)
            return len_A;

        // Sets the first row and the first column of the verification matrix
        for (ushort a = 0; a <= len_A; ++a)
            verif[a][0] = a;
        for (ushort b = 0; b <= len_B; ++b)
            verif[0][b] = b;

        // Verification step / matrix filling.
        for (ushort a = 1; a <= len_A; ++a)
            for (ushort b = 1; b <= len_B; ++b)
                if (seq_A[a - 1] == seq_B[b - 1])
                    verif[a][b] = verif[a - 1][b - 1];
                else
                    verif[a][b] = std::min({verif[a - 1][b - 1], verif[a - 1][b], verif[a][b - 1]}) + 1;

        // Last position of the matrix contains the distance
        return verif[len_A][len_B];
    }

    /* Function to find the closest match. */
    const std::string find_closest(const std::string &input,
                                   const str_vec &options,
                                   ushort min_dist)
    {
        std::string closest = "";
        ushort dist_closest = min_dist;

        for (const std::string &x : options)
        {
            ushort distance = levenshtein(input, x);
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
            std::ostringstream oss;
            oss << " Available option(s): '" << options[0] << "'";
            for (size_t i = 1; i < options.size(); ++i)
                oss << ", '" << options[i] << "'";
            oss << ".";
            return oss.str();
        }
        else
            return " Did you mean '" + closest + "'?";
    }
} // TOOLS NAMESPACE
