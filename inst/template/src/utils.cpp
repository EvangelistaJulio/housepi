/*
 * utils.cpp
 *
 * Copyright (C) 2024 Francesco Camaglia, Institut Pasteur
 */

#include "utils.h"

// auxiliary namespace for I/O custom functions
namespace tools
{
    /////////////////////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////  Auxiliary Variables  ////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////

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

    /* Welford online algorithm for mean and variance. */
    void statistics(const std::vector<double> &data, double &avg, double &var)
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

            _ProgBar_Thin_ = 1;
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
                if (time_lastchecked < 950)
                {
                    ++_ProgBar_Thin_;
                }
                else if (time_lastchecked > 1050)
                {
                    --_ProgBar_Thin_;
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

#ifdef USE_RCPP
        Rcpp::warning(txt);
#else
        std::cerr << "\nWarning: " + txt << std::endl;
#endif
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

    // Template function to wrap conversion functions with try-catch
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
        // throw std::runtime_error("Unhandled error in safe_converter");
    }

    // Converter to int using std::stoi with function pointer casting
    int stoi(const std::string &input, const std::string &errtxt)
    {
        return _safe_converter_(static_cast<int (*)(const std::string &, std::size_t *, int)>(std::stoi), input, errtxt, nullptr, 10);
    }

    // Converter to double using std::stod with function pointer casting
    double stod(const std::string &input, const std::string &errtxt)
    {
        return _safe_converter_(static_cast<double (*)(const std::string &, std::size_t *)>(std::stod), input, errtxt, nullptr);
    }

    // Converter to unsigned long using std::stoul with function pointer casting
    unsigned long stoul(const std::string &input, const std::string &errtxt)
    {
        return _safe_converter_(static_cast<ulong (*)(const std::string &, std::size_t *, int)>(std::stoul), input, errtxt, nullptr, 10);
    }

    // Converter to unsigned long long using std::stoull with function pointer casting
    unsigned long long stoull(const std::string &input, const std::string &errtxt)
    {
        return _safe_converter_(static_cast<unsigned long long (*)(const std::string &, std::size_t *, int)>(std::stoull), input, errtxt, nullptr, 10);
    }
}

//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//

/////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////  Base Function : methods  //////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

/* Overload constructor for argument VALUES. */
Base_Function::Base_Function(const std::string &name,
                             const dictargs &args_val,
                             const std::vector<std::string> &order,
                             const dictargs &domain)
{
    _name_ = name;
    _args_names_ = order;
    _args_val_ = new double[_args_names_.size()];
    int i = 0;
    for (const std::string &argID : order)
    {
        _args_val_[i++] = args_val.at(argID);
    }
    // assign getter
    _get_arg_method_ = &Base_Function::_get_arg_val_;
    // assign setter
    _set_arg_method_ = &Base_Function::_set_arg_val_;
    // assign domain
    set_range(domain);
}

/* Overload constructor for argument POINTERS. */
Base_Function::Base_Function(const std::string &name,
                             const dictargs_ptr &args_ptr,
                             const std::vector<std::string> &order,
                             const dictargs &domain)
{
    _name_ = name;
    _args_names_ = order;
    _args_ptr_ = new const double *[args_ptr.size()];
    int i = 0;
    for (const std::string &argID : order)
    {
        _args_ptr_[i++] = args_ptr.at(argID);
    }
    // assign getter
    _get_arg_method_ = &Base_Function::_get_arg_ptr_;
    // assign setter
    _set_arg_method_ = &Base_Function::_set_arg_ptr_;
    // assign domain
    set_range(domain);
}

/* Overload assignment of non default POINTER arguments. */
void Base_Function::assign_args(dictargs args_val)
{
    for (const auto &pair : args_val)
    {
        auto it = std::find(_args_names_.begin(), _args_names_.end(), pair.first);
        if (it == _args_names_.end())
        {
            std::ostringstream oss;
            oss << "Unknown argument '" << pair.first << "' for '" << get_name() << "'.";
            throw std::invalid_argument(oss.str());
        }
        else
        {
            _args_val_[std::distance(_args_names_.begin(), it)] = pair.second;
        }
    }
}

/* Overload assignment of non default POINTER arguments. */
void Base_Function::assign_args(dictargs_ptr args_ptr)
{
    for (const auto &pair : args_ptr)
    {
        auto it = std::find(_args_names_.begin(), _args_names_.end(), pair.first);
        if (it == _args_names_.end())
        {
            std::ostringstream oss;
            oss << "Unknown argument '" << pair.first << "' for '" << get_name() << "'.";
            throw std::invalid_argument(oss.str());
        }
        else
        {
            _args_ptr_[std::distance(_args_names_.begin(), it)] = pair.second;
        }
    }
}

//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//

/////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////  Base Function : set range  ////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

/* Set the range of validity of the value, and the corresponding validation function. */
int Base_Function::set_range(const dictargs &new_range)
{
    double min_val = _range_min_;
    double max_val = _range_max_;
    bool is_new_min = false;
    bool is_new_max = false;
    int output = 0;

    //===== Read the new range =====//
    for (auto const &x : new_range)
    {
        if (x.first == "min")
        {
            if (is_new_min)
            {
                throw std::invalid_argument("Multiple definition for 'min' in '" + get_name() + "'.");
            }
            else
            {
                min_val = x.second;
                is_new_min = true;
            }
        }
        else if (x.first == "max")
        {
            if (is_new_max)
            {
                throw std::invalid_argument("Multiple definition for 'max' in '" + get_name() + "'.");
            }
            else
            {
                max_val = x.second;
                is_new_max = true;
            }
        }
        else
        {
            std::ostringstream oss;
            oss << "Invalid range keyword '" << x.first << "'. Provide only 'min' and/or 'max'.";
            throw std::invalid_argument(oss.str());
        }
    }

    //===== Adjust new min value and domain =====//
    if (std::find(_args_names_.begin(), _args_names_.end(), "min") != _args_names_.end())
    {
        //----- Overwrite min if it exists as argument (e .g. uniform distribution) -----//
        min_val = get_arg(std::find(_args_names_.begin(), _args_names_.end(), "min") - _args_names_.begin());
        if (is_new_min)
        {
            std::ostringstream oss;
            oss << "Overwriting range min with argument min " << min_val << " in '" << get_name() << "'.";
            tools::warning(oss.str());
        }
        is_new_min = true;
    }
    else if (is_new_min)
    {
        //----- Check compatibility with min from domain (e .g. uniform distribution) -----//
        if (min_val < _range_min_)
        {
            std::ostringstream oss;
            oss << "Invalid min value " << min_val << " for distribution '" << get_name() << "'.";
            throw std::invalid_argument(oss.str());
        }
    }

    //===== Adjust new max value and domain =====//
    if (std::find(_args_names_.begin(), _args_names_.end(), "max") != _args_names_.end())
    {
        //----- Overwrite max if it exists as argument -----//
        max_val = get_arg(std::find(_args_names_.begin(), _args_names_.end(), "max") - _args_names_.begin());
        if (is_new_max)
        {
            std::ostringstream oss;
            oss << "Overwriting range max with argument max " << max_val << " in '" << get_name() << "'.";
            tools::warning(oss.str());
        }
        is_new_max = true;
    }
    else if (is_new_max)
    {
        //----- Check compatibility with max from domain -----//
        if (max_val > _range_max_)
        {
            std::ostringstream oss;
            oss << "Invalid max value" << max_val << " for distribution '" << get_name() << "'.";
            throw std::invalid_argument(oss.str());
        }
    }

    //===== Assign range and respective validation method =====//

    // at least one from previous range (domain) should be finite
    // FIXME : to be moved to NA...
    is_new_min = (is_new_min != -INF_DOUBLE) || (_range_min_ != -INF_DOUBLE);
    is_new_max = (is_new_max != INF_DOUBLE) || (_range_max_ != INF_DOUBLE);

    if ((!is_new_min) && (!is_new_max))
    {
        output = 0;
        _valid_range_method_ = &Base_Function::_valid_range_always_;
    }
    else if ((is_new_min) && (!is_new_max))
    {
        output = 1;
        _range_min_ = min_val;
        _valid_range_method_ = &Base_Function::_valid_range_min_;
    }
    else if ((!is_new_min) && (is_new_max))
    {
        output = 2;
        _range_max_ = max_val;
        _valid_range_method_ = &Base_Function::_valid_range_max_;
    }
    else // ((is_new_min) && (is_new_max))
    {
        output = 3;
        if (_range_min_ >= _range_max_)
        {
            std::ostringstream oss;
            oss << "Provided an inconsistent range {min:" << _range_min_ << ", max:" << _range_max_ << "}.";
            throw std::invalid_argument(oss.str());
        }
        else
        {
            _range_min_ = min_val;
            _range_max_ = max_val;
            _valid_range_method_ = &Base_Function::_valid_range_both_;
        }
    }

    return output;
}

//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//

/////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////  Base Function : show args  ////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

/* Return formatted argument names of the distribution. */
std::string Base_Function::show_arg_names() const
{
    if (_args_names_.empty())
    {
        return UNDEFINED;
    }
    else
    {
        std::ostringstream oss;
        oss << _args_names_[0];
        for (size_t i = 1; i < _args_names_.size(); ++i)
        {
            oss << tools::C_LIST_DELIM << _args_names_[i];
        }
        return oss.str();
    }
}

/* Return formatted argument names of the distribution. */
std::string Base_Function::show_arg_values() const
{
    if (_args_names_.empty())
    {
        return UNDEFINED;
    }
    else
    {
        std::ostringstream oss;
        oss << get_arg(0);
        for (size_t i = 1; i < _args_names_.size(); ++i)
        {
            oss << tools::C_LIST_DELIM << get_arg(i);
        }
        return oss.str();
    }
}