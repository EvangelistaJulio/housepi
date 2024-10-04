/*
 * utils.cpp
 *
 * Copyright (C) 2024 Francesco Camaglia, Institut Pasteur
 */

#include "utils.h"

// auxiliary namespace for I/O custom functions
namespace basics
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

    /* An internal check for string representing NA. */
    bool is_na(const string &txt)
    {
        if ((txt == UNDEFINED) || (txt == "NA") || (txt.size() == 0))
        {
            return true;
        }
        else
        {
            return false;
        }
    }

    /* Welford online algorithm for mean and variance. */
    void statistics(const vector<double> &data, double &avg, double &var)
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

    /* A formatted string for percent of num/denom */
    string percent(double num, double den)
    {
        std::ostringstream oss;
        oss << std::fixed << std::setprecision(1)
            << "(" << (100.0 * num / den) << "%)";
        return oss.str();
    }

    /* A formatted string for TRUE/FALSE values. */
    string truth(bool condition)
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
    void print_stdout(const string &txt)
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
    void _helper_stderr_(const string &txt)
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
    void print_decorator(const string &title, const string &txt, char decor, bool verbose)
    {
        if (verbose)
        {
            std::ostringstream oss;
            oss << string(DECOR_LEFT, decor) << " " << title << " "
                << string(DECOR_WIDTH - title.length() - DECOR_LEFT, decor) << "\n";
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
    string format_millisecs(long long elapsed)
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
    void progress_bar_init(const string &txt, bool verbose)
    {
        using namespace std::chrono;

        if (verbose)
        {
            std::ostringstream oss;
            oss << "\r"
                << txt << ": [" + string(PROGBAR_WIDTH, '.') << "] "
                << percent(0., 1.) << " 00m00s<??m??s"
                << string(20, ' ');
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
    void progress_bar_run(const string &txt, ulong progress, ulong total, bool verbose)
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

                std::ostringstream oss;
                oss << "\r"
                    << txt << ": [" << string(where, C_PROG) << string(PROGBAR_WIDTH - where, '.') << "] "
                    << percent(progress, total)
                    << " " << format_millisecs(time_elapsed) << "<"
                    << ((time_elapsed > 5000) ? format_millisecs(time_missing) : "??m??s")
                    << string(20, ' ');
                _helper_stderr_(oss.str());

                // adapt progressbar thinning
                if (time_lastchecked < 1925)
                {
                    ++_ProgBar_Thin_;
                }
                else if (time_lastchecked > 2075)
                {
                    --_ProgBar_Thin_;
                }
                _ProgBar_LastTime_ = steady_clock::now();
            }
        }

#ifdef USE_RCPP
        Rcpp::checkUserInterrupt();
#endif
    }

    /*  End the progress bar. */
    void progress_bar_end(const string &txt, bool verbose)
    {
        using namespace std::chrono;

        if (verbose)
        {
            auto time_elapsed = duration_cast<milliseconds>(steady_clock::now() - _ProgBar_Start_).count();
            std::ostringstream oss;
            oss << "\r"
                << txt << ": [" << string(PROGBAR_WIDTH, C_PROG) << "] "
                << format_millisecs(time_elapsed)
                << string(25, ' ') << "\n";
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

    /* Raise a warning if it wasn't raised before. */
    // FIXME : omp critical is so wrong! Throw and check exceptions instead depending on Rcpp
    void raise_warning(const string &txt, bool &is_raised_once)
    {
        if (is_raised_once == false)
        {
#pragma omp critical
            {
#ifdef USE_RCPP
                Rcpp::warning(txt);
#else
                std::cerr << "\nWarning: " + txt << std::endl;
#endif
                is_raised_once = true;
            }
        }
    }

    /////////////////////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////  Error Functions  ////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////

    /* Display an error message. */
    // FIXME : omp critical is so wrong! Throw and check exceptions instead depending on Rcpp
    void raise_error(const string &txt)
    {
#pragma omp critical
        {
#ifdef USE_RCPP
            Rcpp::stop(txt);
#else
            std::cerr << "\nError: " + txt + "\nExit...\n"
                      << std::endl;
            std::exit(EXIT_FAILURE);
#endif
        }
    }
}

using namespace basics;

/////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////  WARNINGS  /////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

bool Unused_Range_Warning = false; // warning for invalid proposer range declaration

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
Base_Function::Base_Function(const string &name,
                             const dictargs &args_val,
                             const vector<string> &order,
                             const dictargs &domain)
{
    _name_ = name;
    _args_names_ = order;
    _args_val_ = new double[_args_names_.size()];
    int i = 0;
    for (const string &argID : order)
    {
        _args_val_[i++] = args_val.at(argID);
    }
    // assign getter
    _get_arg_method_ = &Base_Function::_get_arg_val_;
    // assign setter
    _set_arg_method_ = &Base_Function::_set_arg_val_;
    // assign domain
    set_range(domain, true);
}

/* Overload constructor for argument POINTERS. */
Base_Function::Base_Function(const string &name,
                             const dictargs_ptr &args_ptr,
                             const vector<string> &order,
                             const dictargs &domain)
{
    _name_ = name;
    _args_names_ = order;
    _args_ptr_ = new const double *[args_ptr.size()];
    int i = 0;
    for (const string &argID : order)
    {
        _args_ptr_[i++] = args_ptr.at(argID);
    }
    // assign getter
    _get_arg_method_ = &Base_Function::_get_arg_ptr_;
    // assign setter
    _set_arg_method_ = &Base_Function::_set_arg_ptr_;
    // assign domain
    set_range(domain, true);
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
            raise_error(oss.str());
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
            raise_error(oss.str());
        }
        else
        {
            _args_ptr_[std::distance(_args_names_.begin(), it)] = pair.second;
        }
    }
}

/* Set the range of validity of the value, and the corresponding validation function. */
void Base_Function::set_range(dictargs new_range, bool is_domain)
{
    bool is_min = false;
    bool is_max = false;
    double new_min = -INF_DOUBLE;
    double new_max = INF_DOUBLE;

    //----- Read the range -----//
    for (auto const &x : new_range)
    {
        if (x.first == "min")
        {
            if (is_min)
            {
                raise_error("Multiple definition for 'min'.");
            }
            else
            {
                new_min = x.second;
                is_min = true;
            }
        }
        else if (x.first == "max")
        {
            if (is_max)
            {
                raise_error("Multiple definition for 'max'.");
            }
            else
            {
                new_max = x.second;
                is_max = true;
            }
        }
        else
        {
            std::ostringstream oss;
            oss << "Invalid range keyword '" << x.first << "'. Provide only 'min' and/or 'max'.";
            raise_error(oss.str());
        }
    }

    if (!is_domain)
    {
        //----- Overwrite min if it exists as argument -----//
        if (std::find(_args_names_.begin(), _args_names_.end(), "min") != _args_names_.end())
        {
            new_min = get_arg(std::find(_args_names_.begin(), _args_names_.end(), "min") - _args_names_.begin());
            if (is_min)
            {
                std::ostringstream oss;
                oss << "Overwriting range min with argument min " << new_min << " in '" << get_name() << "'.";
                raise_warning(oss.str(), Unused_Range_Warning);
            }
            is_min = true;
        }

        //----- Overwrite max if it exists as argument -----//
        if (std::find(_args_names_.begin(), _args_names_.end(), "max") != _args_names_.end())
        {
            new_max = get_arg(std::find(_args_names_.begin(), _args_names_.end(), "max") - _args_names_.begin());
            if (is_max)
            {
                std::ostringstream oss;
                oss << "Overwriting range max with argument max " << new_max << " in '" << get_name() << "'.";
                raise_warning(oss.str(), Unused_Range_Warning);
            }
            is_max = true;
        }

        //----- Check compatibility with domain -----//
        if ((is_min) && (new_min < _range_min_))
        {
            std::ostringstream oss;
            oss << "Invalid min value " << new_min << " for distribution '" << get_name() << "'.";
            raise_error(oss.str());
        }
        if ((is_max) && (new_max > _range_max_))
        {
            std::ostringstream oss;
            oss << "Invalid max value" << new_max << " for distribution '" << get_name() << "'.";
            raise_error(oss.str());
        }
    }

    //----- Assign valid range functional method -----//
    if ((!is_min) && (!is_max))
    {
        _valid_range_method_ = &Base_Function::_valid_range_always_;
    }
    else if ((is_min) && (!is_max))
    {
        _range_min_ = new_min;
        _valid_range_method_ = &Base_Function::_valid_range_min_;
    }
    else if ((!is_min) && (is_max))
    {
        _range_max_ = new_max;
        _valid_range_method_ = &Base_Function::_valid_range_max_;
    }
    else // ((is_min) && (is_max))
    {
        if (_range_min_ >= _range_max_)
        {
            std::ostringstream oss;
            oss << "Provided an inconsistent range {min:" << _range_min_ << ", max:" << _range_max_ << "}.";
            raise_error(oss.str());
        }
        else
        {
            _range_min_ = new_min;
            _range_max_ = new_max;
            _valid_range_method_ = &Base_Function::_valid_range_both_;
        }
    }
}

/* Return formatted argument names of the distribution. */
string Base_Function::show_arg_names() const
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
            oss << C_LIST_DELIM << _args_names_[i];
        }
        return oss.str();
    }
}

/* Return formatted argument names of the distribution. */
string Base_Function::show_arg_values() const
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
            oss << C_LIST_DELIM << get_arg(i);
        }
        return oss.str();
    }
}