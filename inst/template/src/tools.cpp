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
    constexpr ulong DECOR_LEFT = 3;     // left size of decorator
    constexpr ulong DECOR_WIDTH = 65;   // total size of decorator
    constexpr ulong PROGBAR_WIDTH = 30; // total size of progress bar
    ulong _ProgBar_Thin_;               // thinning of progress bar
    std::chrono::steady_clock::time_point _ProgBar_Start_;
    std::chrono::steady_clock::time_point _ProgBar_LastTime_;

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

    /* A formatted std::string for percent of num/denom */
    std::string percent(double num, double den)
    {
        std::ostringstream oss;
        oss << std::fixed << std::setprecision(1)
            << "(" << (100.0 * num / den) << "%)";
        return oss.str();
    }

    /* A formatted std::string for TRUE/FALSE values. */
    std::string truth(bool condition)
    {
#ifdef USE_RCPP
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
#ifdef USE_RCPP
        Rcpp::Rcout << txt;
        Rcpp::Rcout.flush();
#else
        std::cout << txt;
        std::cout.flush();
#endif
    }

    /* Internal helper function to change stderr in case of RCPP. */
    void _helper_stderr_(const std::string &txt)
    {
#ifdef USE_RCPP
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
        if (verbose)
        {
            std::ostringstream oss;
            oss << std::string(DECOR_LEFT, decor) << " " << title << " "
                << std::string(DECOR_WIDTH - title.length() - DECOR_LEFT, decor) << "\n";
            if (!is_na(txt))
            {
                oss << txt << "\n";
            }
            print_stdout(oss.str());
        }
    }

    //-----------------------------------------------------------------------------------------------//
    //-----------------------------------------------------------------------------------------------//
    //-----------------------------------------------------------------------------------------------//
    //-----------------------------------------------------------------------------------------------//
    //-----------------------------------------------------------------------------------------------//
    //-----------------------------------------------------------------------------------------------//

    /////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////  Progress Bar Functions  ///////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////

    /* Format the  */
    std::string format_millisecs(long long elapsed)
    {
        using namespace std::chrono;

        auto dys = duration_cast<duration<int, std::ratio<86400>>>(milliseconds(elapsed));
        elapsed -= duration_cast<milliseconds>(dys).count();

        auto hrs = duration_cast<hours>(milliseconds(elapsed));
        elapsed -= duration_cast<milliseconds>(hrs).count();

        auto mns = duration_cast<minutes>(milliseconds(elapsed));
        elapsed -= duration_cast<milliseconds>(mns).count();

        auto scs = duration_cast<seconds>(milliseconds(elapsed));
        elapsed -= duration_cast<milliseconds>(scs).count();

        std::ostringstream oss;

        // show large time scales only if the case
        if (dys.count() > 1)
        {
            oss << dys.count() << "d";
        }
        if (hrs.count() > 1)
        {
            oss << std::setfill('0') << std::setw(2) << hrs.count() << "h";
        }

        oss << std::setfill('0') << std::setw(2) << mns.count() << "m"
            << std::setfill('0') << std::setw(2) << scs.count() << "s";

        return oss.str();
    }

    /* Initialize progress bar. */
    void progress_bar_init(const std::string &txt, bool verbose)
    {
        using namespace std::chrono;

        if (verbose)
        {
            std::ostringstream oss;
            oss << "\r"
                << txt << ": [" + std::string(PROGBAR_WIDTH, '.') << "] "
                << percent(0., 1.) << " 00m00s<??m??s"
                << std::string(20, ' ');
            _helper_stderr_(oss.str());

            _ProgBar_Thin_ = 5;
            _ProgBar_Start_ = steady_clock::now();
            _ProgBar_LastTime_ = steady_clock::now();
        }

#ifdef USE_RCPP
        Rcpp::checkUserInterrupt();
#endif
    }

    /* Percentual progress bar (progress/total). */
    void progress_bar_run(const std::string &txt, ulong progress, ulong total, bool verbose)
    {
        using namespace std::chrono;

        if (verbose)
        {
            if (progress % _ProgBar_Thin_ == 0)
            {
                ulong where = static_cast<double>(progress) * PROGBAR_WIDTH / total;
                std::chrono::steady_clock::time_point now = steady_clock::now();
                long long time_elapsed = duration_cast<milliseconds>(now - _ProgBar_Start_).count();
                long long time_lastchecked = duration_cast<milliseconds>(now - _ProgBar_LastTime_).count();
                long long time_missing = time_elapsed * (static_cast<double>(total) / progress - 1.);
                _ProgBar_LastTime_ = steady_clock::now();

                std::ostringstream oss;
                oss << "\r"
                    << txt << ": [" << std::string(where, C_PROG) << std::string(PROGBAR_WIDTH - where, '.') << "] "
                    << percent(progress, total)
                    << " " << format_millisecs(time_elapsed) << "<"
                    << ((time_elapsed > 4000) ? format_millisecs(time_missing) : "??m??s")
                    << std::string(20, ' ');
                _helper_stderr_(oss.str());

                // adapt progressbar thinning
                if (time_lastchecked < 900)
                {
                    _ProgBar_Thin_ *= 2;
                }
                else if (time_lastchecked > 1100)
                {
                    _ProgBar_Thin_ /= 3;
                }
            }
        }

#ifdef USE_RCPP
        Rcpp::checkUserInterrupt();
#endif
    }

    /*  End the progress bar. */
    void progress_bar_end(const std::string &txt, bool verbose)
    {
        using namespace std::chrono;

        if (verbose)
        {
            auto time_elapsed = duration_cast<milliseconds>(steady_clock::now() - _ProgBar_Start_).count();
            std::ostringstream oss;
            oss << "\r"
                << txt << ": [" << std::string(PROGBAR_WIDTH, C_PROG) << "] "
                << format_millisecs(time_elapsed)
                << std::string(25, ' ') << "\n";
            _helper_stderr_(oss.str());
        }

#ifdef USE_RCPP
        Rcpp::checkUserInterrupt();
#endif
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
#ifdef USE_RCPP
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
#ifdef USE_RCPP
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
            throw std::invalid_argument("Invalid value '" + input + "'. Impossible to import " + errtxt + ".");
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

#ifdef USE_RCPP
    void rdataframe_to_map(Rcpp::DataFrame df, std::map<std::string, std::vector<std::string>> &output)
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
    std::string strcleaner(std::string tkn)
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
            std::ostringstream oss;

            // check the validity of the headers
            for (const std::string &req : required_headers)
            {
                if (std::find(headers.begin(), headers.end(), req) == headers.end())
                {
                    if (!is_missing)
                    {
                        // point to the faulty file
                        oss << "Missing column(s): '" + req + "'";
                        is_missing = true;
                    }
                    else
                    {
                        oss << ", '" << req << "'";
                    }
                }
            }

            // raise the error if something was missing
            if (is_missing)
            {
                throw std::invalid_argument(oss.str());
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
            size_t num_cols = 0;

            while (getline(iss, tkn, C_DELIM))
            {
                this_row.emplace(headers[num_cols++], strcleaner(tkn));
            }

            // last column is empty
            for (size_t i = num_cols; i < headers.size(); i++)
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
    std::vector<std::string> parse_list_args(const std::string &tkn, const char delim)
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
                throw std::runtime_error("Impossible to open file '" + pathtofile + "'.");
            }
        }
    }
} // TOOLS NAMESPACE
