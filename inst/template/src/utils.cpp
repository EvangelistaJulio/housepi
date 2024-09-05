/*
 * utils.cpp
 *
 * Copyright (C) 2024 Francesco Camaglia, Institut Pasteur
 */

#include "utils.h"

// auxiliary namespace for I/O custom functions
namespace screen
{
    /////////////////////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////  Auxiliary Variables  ////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////

    std::chrono::steady_clock::time_point _Prog_Bar_Start_;

    /////////////////////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////  Auxiliary Functions  ////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////

    /* A formatted string for percent of num/denom */
    std::string percent(int num, int den)
    {
        std::ostringstream oss;
        oss << std::fixed << std::setprecision(1)
            << "(" << (100.0 * num / den) << "%)";
        return oss.str();
    }

    /* A formatted string for TRUE/FALSE values. */
    std::string truth(bool condition)
    {
        return (condition ? "TRUE" : "FALSE");
    }

    /////////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////  RCPP Instructions  ///////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////

    /* Function to change stdout in case of RCPP. */
    void print_stdout(const std::string txt)
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
    void _helper_stderr_(const std::string txt)
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
    void show_decorator(const std::string title, char decor, const std::string txt, bool verbose)
    {
        if (verbose)
        {
            int sx = 3;
            int dx = DECOR_WIDTH - title.length() - sx; // right half bar
            std::ostringstream oss;
            oss << std::string(sx, decor) << " " << title << " "
                << std::string(dx, decor) + "\n";
            if (txt != UNDEFINED)
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
                << std::fixed << std::setprecision(2) << 0.0 << " %"
                << " 00m00s<??m??s"
                << std::string(20, ' ');
            _helper_stderr_(oss.str());

            _Prog_Bar_Start_ = steady_clock::now();
        }

#ifdef USE_RCPP
        Rcpp::checkUserInterrupt();
#endif
    }

    /* Percentual progress bar (progress/total). */
    void progress_bar_run(const std::string &txt, int progress, int total, bool verbose)
    {
        using namespace std::chrono;

        if (verbose)
        {
            if ((total < PROGBAR_THIN) || (progress % PROGBAR_THIN == 0))
            {
                int where = static_cast<double>(progress) * PROGBAR_WIDTH / total;
                auto time_elapsed = duration_cast<milliseconds>(steady_clock::now() - _Prog_Bar_Start_).count();
                double time_missing = static_cast<double>(time_elapsed) * (static_cast<double>(total) / progress - 1.);

                std::ostringstream oss;
                oss << "\r"
                    << txt << ": [" << std::string(where, C_DECOR) << std::string(PROGBAR_WIDTH - where, '.') << "] "
                    << percent(progress, total)
                    << " " << format_millisecs(time_elapsed) << "<" << format_millisecs(time_missing)
                    << std::string(20, ' ');
                _helper_stderr_(oss.str());
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
            auto time_elapsed = duration_cast<milliseconds>(steady_clock::now() - _Prog_Bar_Start_).count();
            std::ostringstream oss;
            oss << "\r"
                << txt << ": [" << std::string(PROGBAR_WIDTH, C_DECOR) << "] "
                << percent(100, 100)
                << " " << format_millisecs(time_elapsed)
                << std::string(20, ' ') << "\n";
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
    void raise_warning(const std::string &txt, bool &is_raised_once)
    {
        if (is_raised_once == false)
        {
#ifdef USE_RCPP
            Rcpp::warning(txt);
#else
            std::cerr << "\nWarning: " + txt << std::endl;
#endif
            is_raised_once = true;
        }
    }

    /////////////////////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////  Error Functions  ////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////

    /* Display an error message. */
    void raise_error(const std::string &txt)
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

//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//

/////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////  Base Function : methods  //////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

/* Overload constructor for INTERNAL arguments. */
Base_Function::Base_Function(sdargs new_args_int, std::vector<std::string> order)
{
    args_name = order;                          // developer assigment of the order of arguments
    args_int = new double[new_args_int.size()]; // dynamic allocation
    for (size_t i = 0; i < args_name.size(); i++)
    {
        args_int[i] = new_args_int.at(args_name[i]);
    }
    _get_arg_func_ = &Base_Function::_get_arg_int_; // assign getter
    _set_arg_func_ = &Base_Function::_set_arg_int_; // assign setter
}

/* Overload constructor for EXTERNAL arguments. */
Base_Function::Base_Function(sdargs_ext new_args_ext, std::vector<std::string> order)
{
    args_name = order;                                  // developer assigment of the order of arguments
    args_ext = new const double *[new_args_ext.size()]; // dynamic allocation
    for (size_t i = 0; i < args_name.size(); i++)
    {
        args_ext[i] = new_args_ext.at(args_name[i]);
    }
    _get_arg_func_ = &Base_Function::_get_arg_ext_; // assign getter
    _set_arg_func_ = &Base_Function::_set_arg_ext_; // assign setter
}

/* Return formatted argument names of the distribution. */
std::string Base_Function::show_arg_names() const
{
    if (args_name.empty())
    {
        return UNDEFINED;
    }
    else
    {
        std::ostringstream oss;
        oss << args_name[0];
        for (size_t x = 1; x < args_name.size(); x++)
        {
            oss << screen::C_LIST_DELIM << args_name[x];
        }
        return oss.str();
    }
}

/* Return formatted argument names of the distribution. */
std::string Base_Function::show_arg_values() const
{
    if (args_name.empty())
    {
        return UNDEFINED;
    }
    else
    {
        std::ostringstream oss;
        oss << get_arg(0);
        for (size_t x = 1; x < args_name.size(); x++)
        {
            oss << screen::C_LIST_DELIM << get_arg(x);
        }
        return oss.str();
    }
}

//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//

/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////  Uniform Distribution : methods  //////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

/* uniform probability density function */
double Uniform_Distr::pdf(double x)
{
    double low = get_arg(0);
    double up = get_arg(1);
    (void)x;
    return pow(up - low, -1);
}

/* uniform probability density function */
double Uniform_Distr::log_pdf(double x)
{
    double low = get_arg(0);
    double up = get_arg(1);
    (void)x;
    return -log(up - low);
}

/* uniform probability density function */
double Uniform_Distr::logratio(double x, double y)
{
    (void)x;
    (void)y;
    return 0.0;
}

/* uniform random generator */
double Uniform_Distr::rand(sys_rnd_gen &rng)
{
    double low = get_arg(0);
    double up = get_arg(1);
    std::uniform_real_distribution<double> distr_gen(low, up);
    return distr_gen(rng);
}

/////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////  Normal Distribution : methods  /////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

/* normal probability density function */
double Normal_Distr::pdf(double x)
{
    double mean = get_arg(0);
    double sd = get_arg(1);
    return exp(-0.5 * pow((x - mean) / sd, 2)) * pow(sd, -1) * 0.398942280;
}

/* logarithm of the normal probability density function */
double Normal_Distr::log_pdf(double x)
{
    double mean = get_arg(0);
    double sd = get_arg(1);
    return -0.5 * pow((x - mean) / sd, 2) - log(sd) - 0.918938533;
}

/* log-ratio of the normal probability density function */
double Normal_Distr::logratio(double x, double y)
{
    double mean = get_arg(0);
    double sd = get_arg(1);
    return 0.5 * (pow((x - mean) / sd, 2) - pow((y - mean) / sd, 2));
}

/* normal random generator */
double Normal_Distr::rand(sys_rnd_gen &rng)
{
    double mean = get_arg(0);
    double sd = get_arg(1);
    std::normal_distribution<double> distr_gen(mean, sd);
    return distr_gen(rng);
}

/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////  Normal Distribution (discrete) : methods  ////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

double Normal_Distr_int::pdf(double x)
{
    double mean = get_arg(0);
    double sd = get_arg(1);
    int r = round(x);
    return 0.5 * (erf(0.707106781 * (r + 0.5 - mean) / sd) - erf(0.707106781 * (r - 0.5 - mean) / sd));
}

/////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////  Lognormal Distribution : methods  ///////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

/* lognormal probability density function */
double Lognorm_Distr::pdf(double x)
{
    double meanlog = get_arg(0);
    double sdlog = get_arg(1);
    return exp(-0.5 * pow((log(x) - meanlog) / sdlog, 2)) * pow(x * sdlog, -1) * 0.398942280;
}

/* logarithm of the lognormal probability density function */
double Lognorm_Distr::log_pdf(double x)
{
    double meanlog = get_arg(0);
    double sdlog = get_arg(1);
    return -0.5 * pow((log(x) - meanlog) / sdlog, 2) - log(x) - log(sdlog) - 0.918938533;
}

/* log-ratio of the lognormal probability density function */
double Lognorm_Distr::logratio(double x, double y)
{
    double meanlog = get_arg(0);
    double sdlog = get_arg(1);
    return 0.5 * (pow((log(x) - meanlog) / sdlog, 2) - pow((log(y) - meanlog) / sdlog, 2)) + log(x) - log(y);
}

/* lognormal random generator */
double Lognorm_Distr::rand(sys_rnd_gen &rng)
{
    double meanlog = get_arg(0);
    double sdlog = get_arg(1);
    std::lognormal_distribution<double> distr_gen(meanlog, sdlog);
    return distr_gen(rng);
}

/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////  Lognormal Distribution (discrete) : methods  /////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

double Lognorm_Distr_int::pdf(double x)
{
    double meanlog = get_arg(0);
    double sdlog = get_arg(1);
    double r = floor(x);
    double aux = (x < 1) ? -1 : erf(0.707106781 * (log(r) - meanlog) / sdlog);
    return 0.5 * (erf(0.707106781 * (log(r + 1) - meanlog) / sdlog) - aux);
}

/////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////  Gamma Distribution : methods  /////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

/* gamma probability density function */
double Gamma_Distr::pdf(double x)
{
    double shape = get_arg(0);
    double scale = get_arg(1);
    return exp(-x / scale) * pow(x, shape - 1.) / (tgamma(shape) * pow(scale, shape));
}

/* logarithm of the gamma probability density function */
double Gamma_Distr::log_pdf(double x)
{
    double shape = get_arg(0);
    double scale = get_arg(1);
    return (shape - 1) * log(x) - lgamma(shape) - shape * log(scale) - x / scale;
}

/* log-ratio of the gamma probability density function */
double Gamma_Distr::logratio(double x, double y)
{
    double shape = get_arg(0);
    double scale = get_arg(1);
    return (shape - 1) * (log(y / x)) + (x - y) / scale;
}

/* gamma random generator */
double Gamma_Distr::rand(sys_rnd_gen &rng)
{
    double shape = get_arg(0);
    double scale = get_arg(1);
    std::gamma_distribution<double> distr_gen(shape, scale);
    return distr_gen(rng);
}

/////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////  Weibull Distribution : methods  ////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

/* weibull probability density function */
double Weibull_Distr::pdf(double x)
{
    double shape = get_arg(0);
    double scale = get_arg(1);
    return (shape / scale) * pow((x / scale), shape - 1) * exp(-pow(x / scale, shape));
}

/* logarithm of the weibull probability density function */
double Weibull_Distr::log_pdf(double x)
{
    double shape = get_arg(0);
    double scale = get_arg(1);
    return log(shape / scale) + (shape - 1) * log(x / scale) - pow(x / scale, shape);
}

/* log-ratio of the weibull probability density function */
double Weibull_Distr::logratio(double x, double y)
{
    double shape = get_arg(0);
    double scale = get_arg(1);
    return (shape - 1) * (log(y / x)) + pow(x / scale, shape) - pow(y / scale, shape);
}

/* weibull random generator */
double Weibull_Distr::rand(sys_rnd_gen &rng)
{
    double shape = get_arg(0);
    double scale = get_arg(1);
    std::weibull_distribution<double> distr_gen(shape, scale);
    return distr_gen(rng);
}

/////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////  Weibull Distribution (discrete) : methods  ////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

double Weibull_Distr_int::pdf(double x)
{
    double shape = get_arg(0);
    double scale = get_arg(1);
    double r = floor(x);
    return exp(-pow(r / scale, shape)) - exp(-pow((r + 1) / scale, shape));
}

/////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////  Rayleigh Distribution : methods  ///////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

/* rayleigh probability density function */
double Rayleigh_Distr::pdf(double x)
{
    double scale = get_arg(0);
    return x * pow(scale, -2) * exp(-0.5 * pow(x / scale, 2));
}

/* logarithm of the rayleigh probability density function */
double Rayleigh_Distr::log_pdf(double x)
{
    double scale = get_arg(0);
    return log(x) - 2 * log(scale) - 0.5 * pow(x / scale, 2);
}

/* log-ratio of the rayleigh probability density function */
double Rayleigh_Distr::logratio(double x, double y)
{
    double scale = get_arg(0);
    return log(y / x) + 0.5 * pow(scale, -2) * (pow(x, 2) - pow(y, 2));
}

/* rayleigh random generator */
double Rayleigh_Distr::rand(sys_rnd_gen &rng)
{
    double scale = get_arg(0);
    std::uniform_real_distribution<double> unif_gen(0, 1);
    return scale * sqrt(-2 * unif_gen(rng));
}

/////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////  Rayleigh Distribution (discrete) : methods  ///////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

double Rayleigh_Distr_int::pdf(double x)
{
    double scale = get_arg(0);
    double r = floor(x);
    return exp(-0.5 * pow(r / scale, 2)) - exp(-0.5 * pow((r + 1) / scale, 2));
}

/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////  Exponential Distribution : methods  //////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

/* exponential probability density function */
double Exponential_Distr::pdf(double x)
{
    double rate = get_arg(0);
    return rate * exp(-rate * x);
}

/* logarithm of the exponential probability density function */
double Exponential_Distr::log_pdf(double x)
{
    double rate = get_arg(0);
    return log(rate) - rate * x;
}

/* log-ratio of the exponential probability density function */
double Exponential_Distr::logratio(double x, double y)
{
    double rate = get_arg(0);
    return rate * (x - y);
}

/* exponential random generator */
double Exponential_Distr::rand(sys_rnd_gen &rng)
{
    double rate = get_arg(0);
    std::exponential_distribution<double> distr_gen(rate);
    return distr_gen(rng);
}

/////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////  Exponential Distribution (discrete) : methods  ////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

double Exponential_Distr_int::pdf(double x)
{
    double rate = get_arg(0);
    return exp(-rate * floor(x)) * (1 - exp(-rate));
}

double Exponential_Distr_int::log_pdf(double x)
{
    double rate = get_arg(0);
    return -rate * floor(x) + log(1 - exp(-rate));
}

double Exponential_Distr_int::logratio(double x, double y)
{
    double rate = get_arg(0);
    return rate * (floor(x) - floor(y));
}

//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//