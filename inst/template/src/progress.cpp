/*
 * tools.cpp
 *
 * Copyright (C) 2024 Francesco Camaglia, Institut Pasteur
 */

#include "tools.h"

using namespace tools;
using namespace std::chrono;
using std::string;

/////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////  Progress Bar Functions  ///////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

/* Format the  */
string format_millisecs(long long elapsed)
{
    std::ostringstream oss;

    auto dys = duration_cast<duration<int, std::ratio<86400>>>(milliseconds(elapsed));
    elapsed -= duration_cast<milliseconds>(dys).count();

    auto hrs = duration_cast<hours>(milliseconds(elapsed));
    elapsed -= duration_cast<milliseconds>(hrs).count();

    auto mns = duration_cast<minutes>(milliseconds(elapsed));
    elapsed -= duration_cast<milliseconds>(mns).count();

    auto scs = duration_cast<seconds>(milliseconds(elapsed));
    elapsed -= duration_cast<milliseconds>(scs).count();

    // show large time scales only if the case
    if (dys.count() > 1)
    {
        oss << dys.count() << "d"
            << std::setfill('0') << std::setw(2) << hrs.count() << "h"
            << std::setfill('0') << std::setw(2) << mns.count() << "m"
            << std::setfill('0') << std::setw(2) << scs.count() << "s";
    }
    else if (hrs.count() > 1)
    {
        oss << hrs.count() << "h"
            << std::setfill('0') << std::setw(2) << mns.count() << "m"
            << std::setfill('0') << std::setw(2) << scs.count() << "s";
    }
    else if (mns.count() > 1)
    {
        oss << mns.count() << "m"
            << std::setfill('0') << std::setw(2) << scs.count() << "s";
    }
    else
    {
        oss << scs.count() << "s";
    }

    return oss.str();
}

class CustomProgressBarImpl
{
private:
    bool verbose = true;
    std::string text = "";
    steady_clock::time_point start;
    steady_clock::time_point lastprint;
    const long long thintime = 999;
    const std::string symbol = u8"\u25A0"; // decorator for progress bar
    unsigned barwidth = 30;                // total size of progress bar

public:
    CustomProgressBarImpl(bool verbose = true)
    {
        this->verbose = verbose;
#ifdef _HOUSEPI_RCPP
        Rcpp::Environment cli = Rcpp::Environment::namespace_env("cli");
        Rcpp::Function console_width = cli["console_width"];
        barwidth = Rcpp::as<unsigned>(console_width()) / 2;
#endif
    };
    ~CustomProgressBarImpl() {};

    /* Change text displayed */
    void set_text(const string &text)
    {
        this->text = text;
    }

    /* Initialize progress bar. */
    void initialize()
    {
        if (verbose)
        {
            string sms = "init |" + string(barwidth, ' ') + "| " + percent(0., 1.) + " 00m00s<??m??s";
            print_stderr(sms);

            start = steady_clock::now();
            lastprint = start;
        }

#ifdef _HOUSEPI_RCPP
        Rcpp::checkUserInterrupt();
#endif
    }

    /* Percentual progress bar (progress/total). */
    void update(double progress)
    {
        if (verbose)
        {
            steady_clock::time_point now = steady_clock::now();
            long long time_lastchecked = duration_cast<milliseconds>(now - lastprint).count();

            if (time_lastchecked > thintime)
            {
                lastprint = steady_clock::now();

                long long time_elapsed = duration_cast<milliseconds>(now - start).count();
                long long time_missing = time_elapsed * (1. / progress - 1.);

                unsigned where = static_cast<double>(progress) * barwidth;
                string sms = text + " |";
                for (unsigned u = 0; u < where; ++u)
                {
                    sms += symbol;
                }
                sms += string(barwidth - where, ' ') + "| ";
                sms += percent(progress) + " " + format_millisecs(time_elapsed);
                sms += "<" + format_millisecs(time_missing);
                print_stderr("\r" + sms + string(10, ' '));
            }
        }

#ifdef _HOUSEPI_RCPP
        Rcpp::checkUserInterrupt();
#endif
    }

    /*  End the progress bar. */
    void finalize()
    {
        if (verbose)
        {
            auto time_elapsed = duration_cast<milliseconds>(steady_clock::now() - start).count();
            string sms = "end |";
            for (unsigned u = 0; u < barwidth; u++)
            {
                sms += symbol;
            }
            sms += "| " + format_millisecs(time_elapsed);
            print_stderr("\r" + sms + string(25, ' ') + "\n");
        }

#ifdef _HOUSEPI_RCPP
        Rcpp::checkUserInterrupt();
#endif
    }
};

//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//

CustomProgressBar::CustomProgressBar(bool verbose) : pImpl(std::make_unique<CustomProgressBarImpl>(verbose)) {}
CustomProgressBar::~CustomProgressBar() = default;

void CustomProgressBar::set_text(const std::string &text)
{
    pImpl->set_text(text);
}
void CustomProgressBar::initialize()
{
    pImpl->initialize();
}
void CustomProgressBar::update(double progress)
{
    pImpl->update(progress);
}
void CustomProgressBar::finalize()
{
    pImpl->finalize();
}
