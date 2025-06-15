/*      __
 *     / /_  ____  __  __ ____ __  ____  __
 *    / __ \/ __ \/ / / / ___/ _ \/ __ \/ /
 *   / / / / /_/ / /_/ (__  )  __/ /_/ / /
 *  /_/ /_/\____/\____/____/\___/ ____/_/
 *                             /_/
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

// namespace for I/O custom functions
namespace tools
{
    /////////////////////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////  Auxiliary Variables  ////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////

    std::unordered_set<std::string> _issued_warnings_;

    /////////////////////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////  Auxiliary Functions  ////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////

    /* Internal check for std::string representing NA. */
    bool is_na(const std::string &txt)
    {
        const str_vec forbidden = {"<NA>", "NA", "NAN"};

        if (txt.size() == 0)
            return true;
        else if (txt == HOUSEPI_NA_STRING)
            return true;
        else
            for (const std::string &txtnan : forbidden)
                if (stoupper(txt) == txtnan)
                    return true;

        return false;
    }

    /* Internal check for double representing NA. */
    bool is_na(double x)
    {
#ifdef _HOUSEPI_RCPP
        if (Rcpp::traits::is_na<REALSXP>(x))
            return true;
#else
        if (std::isnan(x))
            return true;
#endif
        return false;
    }

    /* A formatted std::string for percentage of num/denom */
    const std::string percentage(double fraction, ushort precision)
    {
        std::ostringstream oss;
        oss << std::fixed << std::setprecision(precision)
            << (100. * fraction) << "%";
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

    //-----------------------------------------------------------------------------------------------//
    //-----------------------------------------------------------------------------------------------//
    //-----------------------------------------------------------------------------------------------//
    //-----------------------------------------------------------------------------------------------//
    //-----------------------------------------------------------------------------------------------//
    //-----------------------------------------------------------------------------------------------//

    /* Reference:
     * ---------
     * Suppressing Call Stack Info in Rcpp-Generated Errors and Warnings
     * Michael Weylandt — written Mar 27, 2018
     * https://gallery.rcpp.org/articles/quiet-stop-and-warning/
     */

    /////////////////////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////  Warning Functions  //////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////

    /* Display a single-line warning message. */
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

    /////////////////////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////  Error Functions  ////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////

    /* Display a single-line error message. */
    void error(const std::string &txt)
    {
#ifdef _HOUSEPI_RCPP
        throw Rcpp::exception(txt.c_str(), false);
        // Rcpp::stop(txt);
#else
        std::cerr << "Error: " + txt + "\nExit..."
                  << std::endl;
        std::exit(EXIT_FAILURE);
#endif
    }

    /////////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////  String Converter Functions ///////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////

    /* Convert to lower case. */
    const std::string stolower(const std::string s)
    {
        std::string out = s;
        std::transform(out.begin(), out.end(), out.begin(), ::tolower);
        return out;
    }

    /* Convert to upper case. */
    const std::string stoupper(const std::string s)
    {
        std::string out = s;
        std::transform(out.begin(), out.end(), out.begin(), ::toupper);
        return out;
    }

    /////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////  Converter Functions ///////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////

    /* Template function to wrap conversion functions with try-catch */
    template <typename Func, typename... Args>
    auto _safe_converter(Func func, const std::string &input, const std::string &errtxt, Args &&...args) -> decltype(func(input, std::forward<Args>(args)...))
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

    /* Converter to int using std::stoi with function pointer casting */
    int stoi(const std::string &input, const std::string &errtxt)
    {
        return _safe_converter(static_cast<int (*)(const std::string &, std::size_t *, int)>(std::stoi), input, errtxt, nullptr, 10);
    }

    /* Converter to unsigned long using std::stoul with function pointer casting */
    unsigned long stoul(const std::string &input, const std::string &errtxt)
    {
        return _safe_converter(static_cast<ulong (*)(const std::string &, std::size_t *, int)>(std::stoul), input, errtxt, nullptr, 10);
    }

    /* Converter to unsigned 32bit integer (for R compatibility). */
    std::int32_t stoi32(const std::string &input, const std::string &errtxt)
    {
        return static_cast<std::int32_t>(stoul(input, errtxt));
    }

    /* Converter to double using std::stod with function pointer casting */
    double stod(const std::string &input, const std::string &errtxt)
    {
        return _safe_converter(static_cast<double (*)(const std::string &, std::size_t *)>(std::stod), input, errtxt, nullptr);
    }

    /* Converter to string checking for NA values. */
    const std::string na_to_str(double x)
    {
        return is_na(x) ? HOUSEPI_NA_STRING : std::to_string(x);
    }

    /* Convert std::string to boolean variable. */
    bool stobool(const std::string &tkn)
    {
        if (stoupper(tkn) == "TRUE" || tkn == "1")
            return true;
        else if (stoupper(tkn) == "FALSE" || tkn == "0")
            return false;
        throw std::invalid_argument(
            "Unrecognized boolean value \"" + tkn + "\".");
    }

    //-----------------------------------------------------------------------------------------------//
    //-----------------------------------------------------------------------------------------------//
    //-----------------------------------------------------------------------------------------------//
    //-----------------------------------------------------------------------------------------------//
    //-----------------------------------------------------------------------------------------------//
    //-----------------------------------------------------------------------------------------------//

    /* Default seed for pseudo random number generator.  */
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
            tkn = HOUSEPI_NA_STRING;

        return tkn;
    }

    /* Read column headers and compare to required headers. */
    str_vec read_column_headers(std::ifstream &input)
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

        return headers;
    }

    /* Load files from the row with headers */
    bool read_row_by_headers(std::ifstream &input, const str_vec &headers, dict_str_str &this_row)
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
                this_row.emplace(headers[i], HOUSEPI_NA_STRING);

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

    /* Extract headers from deflist. */
    const str_vec headers_from_deflist(const dict_str_str &deflist)
    {
        str_vec entries = {};
        for (const auto &pair : deflist)
            entries.push_back(pair.first);
        return entries;
    }

    /* Read csv list store in pathtofile and store it parsed in csv_parsed_list. */
    void read_csv_list(const std::string &pathtofile,
                       dict_str_str &csv_parsed_list)
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
                            "Unable to extract list name. Check file line \"" + line + "\".");

                    //===== CHOICE =====//
                    if (!getline(iss, choice))
                        throw std::runtime_error(
                            "Unable to extract list value. Check file line \"" + line + "\".");

                    // add the new setup
                    csv_parsed_list.emplace(strcleaner(name), strcleaner(choice));
                }
            }
            else
                throw std::runtime_error(
                    "Unable to open file \"" + pathtofile + "\".");
        }
    }

    /* Extract headers from deftable. */
    const str_vec headers_from_deftable(const dict_str_vec &deftable)
    {
        str_vec entries = {};
        for (const auto &pair : deftable)
            entries.push_back(pair.first);
        return entries;
    }

    /* Read csv table in pathtofile and store it parsed in csv_parsed_table. */
    void read_csv_table(const std::string &pathtofile,
                        dict_str_vec &csv_parsed_table)
    {
        if (!is_na(pathtofile))
        {
            std::ifstream input(pathtofile.c_str());

            if (input.is_open())
            {
                // read and check for valid headers
                const str_vec headers = read_column_headers(input);

                // clean all content of previous table
                csv_parsed_table.clear();

                dict_str_str this_row = {};
                while (read_row_by_headers(input, headers, this_row))
                {
                    for (const auto &x : this_row)
                        csv_parsed_table[x.first].push_back(x.second);
                    this_row.clear();
                }
            }
            else
                throw std::runtime_error(
                    "Unable to open file \"" + pathtofile + "\".");
        }
    }

    //-----------------------------------------------------------------------------------------------//
    //-----------------------------------------------------------------------------------------------//
    //-----------------------------------------------------------------------------------------------//
    //-----------------------------------------------------------------------------------------------//
    //-----------------------------------------------------------------------------------------------//
    //-----------------------------------------------------------------------------------------------//

    /* Welford online's algorithm for mean and variance.
        Ref :
        ----
        https://en.wikipedia.org/wiki/Algorithms_for_calculating_variance
    */
    const std::pair<double, double> compute_stats(const dbl_vec &data)
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

    /* Format a list of string options. */
    const std::string show_list(const str_vec &options, size_t max_print)
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
            if (options.size() == 1)
                return " Available option is " + show_list(options);
            else
                return " Available options are: " + show_list(options);
        else
            return " Did you mean \"" + closest + "\"?";
    }

    /* Check for missing entries given a vector of options. */
    void check_for_missing_entries(const str_vec &entries, const str_vec &options)
    {
        // proceed only if there's a given vector of options
        if (!options.empty())
        {
            //----- check for missing entries -----//
            str_vec missing = {};

            for (const std::string &x : options)
                if (!contains(entries, x))
                    missing.push_back(x);

            if (!missing.empty())
            {
                if (missing.size() == 1)
                    throw std::invalid_argument(
                        "Missing column " + show_list(missing));
                else
                    throw std::invalid_argument(
                        "Missing columns: " + show_list(missing));
            }
        }
    }

    /* Check for unused entries given a vector of options. */
    void check_for_unused_entries(const str_vec &entries, const str_vec &options)
    {
        //----- check for unused entries -----//
        str_vec unused = {};
        for (const std::string &x : entries)
            if (!contains(options, x))
                unused.push_back(x);

        if (unused.size() == 1)
            throw std::invalid_argument(
                "Unknown option " + show_list(unused) + find_closest(unused[0], options));

        else if (unused.size() > 1)
            throw std::invalid_argument(
                "Unknown options: " + show_list(unused) + find_closest(HOUSEPI_NA_STRING, options, 0));
    }

} // TOOLS NAMESPACE

//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//

namespace progress
{
    using namespace std::chrono;

    /////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////  Progress Bar Functions  ///////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////

    /* Format elapsed time in seconds. */
    const std::string fmt_time(std::chrono::seconds scs)
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

    /* Implementation class for progress bar.
        Ref :
        OTHMAN EL HAMMOUCHI, April 15, 2024
        "Fancy thread safe progress bar in R for compiled code"
        -----
        https://othmanelhammouchi.com/posts/r-parallel-progress/
    */
    class Bar::Impl
    {
    public:
        //=======================//
        //===== CONSTRUCTOR =====//
        //=======================//

        Impl(const std::string &text, ulong thintime)
        {
#ifdef _HOUSEPI_RCPP
            Rcpp::checkUserInterrupt();
            // change barwidth according to console width
            _barwidth = Rcpp::as<unsigned>(Rcpp::Function("getOption")("width")) / 2;
#endif

            // load options
            _text = text;
            _verbose = thintime > 0;
            _thintime = thintime - 1;
        };
        ~Impl() {};

        //==================//
        //===== SETTER =====//
        //==================//

        /* Set description text of progress bar. */
        void set_text(const std::string &text) { _text = text; }

        /* Start progressbar clock and display firs message. */
        void open(ulong total)
        {
#ifdef _HOUSEPI_RCPP
            Rcpp::checkUserInterrupt();
#endif

            // initialize internal fields
            _total = total;
            _start = steady_clock::now();
            _lastprint = _start;
            _state = 0;
        }

        /* Update progress bar given the progress state. */
        void update()
        {
#ifdef _HOUSEPI_RCPP
            Rcpp::checkUserInterrupt();
#endif
            if (_verbose)
            {
                const double progress = static_cast<double>(_state) / _total;
                steady_clock::time_point now = steady_clock::now();
                seconds elapsed = duration_cast<seconds>(now - _start);
                seconds missing = duration_cast<seconds>(elapsed * (1. / progress - 1.));

                // print update message
                if (static_cast<ulong>(duration_cast<seconds>(now - _lastprint).count()) > _thintime)
                {
                    const unsigned where = progress * _barwidth;
                    std::ostringstream oss;
                    oss << "\r" << _text << " |";
                    for (unsigned u = 0; u < where; ++u)
                        oss << _symbol;
                    oss << std::string(_barwidth - where, ' ') << "| "
                        << tools::percentage(progress, 0) << " "
                        << fmt_time(elapsed) << "<" << fmt_time(missing)
                        << std::string(10, ' ');
                    _display(oss.str());

                    _lastprint = steady_clock::now();
                }
                else if (!_is_opened())
                    // eventually print initial message
                    if (static_cast<ulong>(missing.count()) > _thintime)
                    {
                        std::ostringstream oss;
                        oss << _text << " |" << std::string(_barwidth, ' ') << "| "
                            << tools::percentage(0., 0);
                        _display(oss.str());
                    }
            }

            ++_state;
        }

        /* End the progress bar because of a crash. */
        void crash() const
        {
#ifdef _HOUSEPI_RCPP
            Rcpp::checkUserInterrupt();
#endif
            if ((_verbose) && (!_is_closed()))
                _display("\n");
        }

        /*  Close the progress bar. */
        void close()
        {
#ifdef _HOUSEPI_RCPP
            Rcpp::checkUserInterrupt();
#endif
            if (_verbose)
            {
                seconds elapsed = duration_cast<seconds>(steady_clock::now() - _start);

                if (static_cast<ulong>(elapsed.count()) > _thintime)
                {
                    std::ostringstream oss;
                    oss << "\r" << _text << " |";
                    for (unsigned u = 0; u < _barwidth; u++)
                        oss << _symbol;
                    oss << "| " << fmt_time(elapsed)
                        << std::string(25, ' ') << "\n";
                    _display(oss.str());
                }

                _state = _total;
            }
        }

    private:
        const std::string _symbol = u8"\u25A0"; // bar decorator
        unsigned _barwidth = 40;                // bar total size
        bool _verbose = true;                   // wheter to show bar

        ulong _state = 0;       // current iteration
        ulong _total = 0;       // total number of iterations
        std::string _text = ""; // current text to show

        steady_clock::time_point _start;
        steady_clock::time_point _lastprint;
        ulong _thintime = 0;

        /* Function to display the progress bar. */
        void _display(const std::string &sms) const { tools::print_stdout(sms); }

        /* Tell if progressbar has been closed. */
        bool _is_opened() const { return (_total > 0) && (_state > 0); }

        /* Tell if progressbar has been closed. */
        bool _is_closed() const { return (_total > 0) && (_state == _total); }
    };

    //-----------------------------------------------------------------------------------------------//
    //-----------------------------------------------------------------------------------------------//
    //-----------------------------------------------------------------------------------------------//
    //-----------------------------------------------------------------------------------------------//
    //-----------------------------------------------------------------------------------------------//
    //-----------------------------------------------------------------------------------------------//

    ///////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////  Bar : Constructor  ///////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////

    Bar::Bar(const std::string &text, ulong thintime)
        : pImpl(std::make_unique<Impl>(text, thintime)) {}
    Bar::~Bar() = default;

    ///////////////////////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////  Bar : Setter  /////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////

    void Bar::set_text(const std::string &text) { pImpl->set_text(text); }
    void Bar::open(ulong total) { pImpl->open(total); }
    void Bar::update() { pImpl->update(); }
    void Bar::crash() const { pImpl->crash(); }
    void Bar::close() { pImpl->close(); }

} // PROGRESS NAMESPACE
