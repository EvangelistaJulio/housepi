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

/* Format elapsed time in seconds. */
const string fmt_time(std::chrono::seconds scs)
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
    {
        oss << dys.count() << "d"
            << std::setfill('0') << std::setw(2) << hrs.count() << "h"
            << std::setfill('0') << std::setw(2) << mns.count() << "m"
            << std::setfill('0') << std::setw(2) << scs.count() << "s";
    }
    else if (hrs.count() > 0)
    {
        oss << hrs.count() << "h"
            << std::setfill('0') << std::setw(2) << mns.count() << "m"
            << std::setfill('0') << std::setw(2) << scs.count() << "s";
    }
    else if (mns.count() > 0)
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
    string text = "";
    const string symbol = u8"\u25A0"; // decorator for progress bar
    unsigned barwidth = 30;           // total size of progress bar

    steady_clock::time_point start;
    steady_clock::time_point lastprint;
    const int64_t thintime = 0;

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

    void set_text(const std::string &text)
    {
        this->text = text;
    }

    /* Initialize progress bar. */
    void initialize(const string &text)
    {
        set_text(text);

#ifdef _HOUSEPI_RCPP
        Rcpp::checkUserInterrupt();
#endif
        if (verbose)
        {
            string sms = text + " |" + string(barwidth, ' ') + "| " + percentage(0.0, 0);
            sms += " 0s<?";
            display(sms);
            start = steady_clock::now();
            lastprint = start;
        }
    }

    /* Update progress bar given the progress state (step / total). */
    void update(double progress)
    {
#ifdef _HOUSEPI_RCPP
        Rcpp::checkUserInterrupt();
#endif
        if (verbose)
        {
            steady_clock::time_point now = steady_clock::now();

            if (duration_cast<seconds>(now - lastprint).count() > thintime)
            {
                unsigned where = progress * barwidth;
                string sms = text + " |";
                for (unsigned u = 0; u < where; ++u)
                {
                    sms += symbol;
                }
                sms += string(barwidth - where, ' ') + "| " + percentage(progress, 0);

                seconds elapsed = duration_cast<seconds>(now - start);
                seconds missing = duration_cast<seconds>(elapsed * (1. / progress - 1.));
                sms += " " + fmt_time(elapsed) + "<" + fmt_time(missing);
                display("\r" + sms + string(10, ' '));

                lastprint = steady_clock::now();
            }
        }
    }

    /*  End the progress bar. */
    void finalize()
    {
#ifdef _HOUSEPI_RCPP
        Rcpp::checkUserInterrupt();
#endif
        if (verbose)
        {
            seconds elapsed = duration_cast<seconds>(steady_clock::now() - start);
            string sms = text + " |";
            for (unsigned u = 0; u < barwidth; u++)
            {
                sms += symbol;
            }
            sms += "| " + fmt_time(elapsed);
            display("\r" + sms + string(25, ' ') + "\n");
        }
    }

private:
    void display(const std::string &sms)
    {
        print_stdout(sms);
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

void CustomProgressBar::set_text(const std::string &text) { pImpl->set_text(text); }
void CustomProgressBar::initialize(const std::string &text) { pImpl->initialize(text); }
void CustomProgressBar::update(double progress) { pImpl->update(progress); }
void CustomProgressBar::finalize() { pImpl->finalize(); }
