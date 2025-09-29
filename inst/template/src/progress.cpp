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
 * progress.cpp
 *
 */

#include <algorithm>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <iostream>
#include <chrono>

#include "tools.h"
#include "progress.h"

namespace progress
{
    using namespace std::chrono;

    /////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////  Progress Bar Functions  ///////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////

    /// Format elapsed time in seconds.
    std::string format(std::chrono::seconds scs)
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
    /////////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////// *** Bar : Impl : class *** ///////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////

    /// Implementation class for progress bar.
    ///
    /// -----
    /// Ref :
    /// OTHMAN EL HAMMOUCHI, April 15, 2024
    /// "Fancy thread safe progress bar in R for compiled code"
    /// https://othmanelhammouchi.com/posts/r-parallel-progress/
    class Bar::Impl
    {
    private:
        const std::string symbol_ = u8"\u25A0"; // bar decorator
        unsigned barwidth_ = 35;                // bar total size
        bool verbose_ = true;                   // wheter to show bar
        bool is_bar_print_ = false;             // wheter bar is printing

        unsigned long state_ = 0; // current iteration
        unsigned long total_ = 0; // total number of iterations
        std::string text_ = "";   // current text to show

        steady_clock::time_point start_;
        steady_clock::time_point lastprint_;
        unsigned long thintime_ = 0;

        void print_(const std::string &sms) const
        {
#ifdef _HOUSEPI_RCPP
            Rcpp::Rcout << sms;
            Rcpp::Rcout.flush();
#else
            std::cout << sms;
            std::cout.flush();
#endif
        }

        /// Function to display the progress bar.
        void bar_(const std::string &sms) const { (this->*which_bar_)(sms); }

        void (Impl::*which_bar_)(const std::string &) const = &Impl::do_bar_carriage_return_;
        void do_bar_carriage_return_(const std::string &sms) const { print_("\r" + sms); }
        void do_bar_new_line_(const std::string &sms) const { print_(sms + "\n"); }

    public:
        /////////////////////////////////////////////////////////////////////////////////////////////
        /////////////////////////////////  Bar : Impl : constructor /////////////////////////////////
        /////////////////////////////////////////////////////////////////////////////////////////////

        Impl()
        {
#ifdef _HOUSEPI_RCPP
            Rcpp::checkUserInterrupt();
            // change barwidth according to console width
            barwidth_ = Rcpp::as<unsigned>(Rcpp::Function("getOption")("width")) / 2;

            // change if R interactive environment
            bool is_interactive = Rcpp::as<bool>(Rcpp::Function("interactive")());
            if (is_interactive)
                which_bar_ = &Impl::do_bar_carriage_return_;
            else
                which_bar_ = &Impl::do_bar_new_line_;
#endif
        };
        ~Impl() {};

        /////////////////////////////////////////////////////////////////////////////////////////////
        ////////////////////////////////////  Bar : Impl : setter ///////////////////////////////////
        /////////////////////////////////////////////////////////////////////////////////////////////

        /// Set description text of progress bar.
        void set_text(const std::string &text) { text_ = text; }

        /// Set description text of progress bar.
        void set_thintime(unsigned long thintime)
        {
            verbose_ = thintime > 0;
            thintime_ = thintime - 1;
        }

        /// Start progressbar clock and display firs message.
        void open(unsigned long total)
        {
#ifdef _HOUSEPI_RCPP
            Rcpp::checkUserInterrupt();
#endif
            if (total == 0)
                throw std::invalid_argument(
                    "Progressbar total should be greater than 0.");

            // initialize internal fields
            total_ = total;
            start_ = steady_clock::now();
            lastprint_ = start_;
            state_ = 0;
            is_bar_print_ = false;
        }

        /// Update progress bar given the progress state.
        void update()
        {
#ifdef _HOUSEPI_RCPP
            Rcpp::checkUserInterrupt();
#endif
            // update internal state
            ++state_;

            if (verbose_)
            {
                const double progress = static_cast<double>(state_) / total_;
                steady_clock::time_point now = steady_clock::now();
                seconds elapsed = duration_cast<seconds>(now - start_);
                seconds missing = duration_cast<seconds>(elapsed * (1. / progress - 1.));

                // print update message
                if (static_cast<unsigned long>(duration_cast<seconds>(now - lastprint_).count()) > thintime_)
                {
                    const unsigned where = progress * barwidth_;
                    std::ostringstream oss;
                    oss << text_ << " |";
                    for (unsigned u = 0; u < where; ++u)
                        oss << symbol_;
                    oss << std::string(barwidth_ - where, ' ') << "| "
                        << std::fixed << std::setprecision(0) << 100. * progress << "% "
                        << format(elapsed) << "<" << format(missing)
                        << std::string(10, ' ');
                    bar_(oss.str());

                    lastprint_ = steady_clock::now();
                    is_bar_print_ = true;
                }
                else if (!is_bar_print_) // eventually print initial message
                    if (static_cast<unsigned long>(missing.count()) > thintime_)
                    {
                        std::ostringstream oss;
                        oss << text_ << " |"
                            << std::string(barwidth_, ' ') << "| 0% ";
                        if (which_bar_ == &Impl::do_bar_new_line_)
                            oss << "\n";
                        print_(oss.str());

                        is_bar_print_ = true;
                    }
            }
        }

        /// End the progress bar because of a crash.
        void crash() const
        {
#ifdef _HOUSEPI_RCPP
            Rcpp::checkUserInterrupt();
#endif
            if (is_bar_print_)
                print_("\n");
        }

        /// Close the progress bar.
        void close()
        {
#ifdef _HOUSEPI_RCPP
            Rcpp::checkUserInterrupt();
#endif
            if (verbose_)
            {
                seconds elapsed = duration_cast<seconds>(steady_clock::now() - start_);

                if (static_cast<unsigned long>(elapsed.count()) > thintime_)
                {
                    std::ostringstream oss;
                    oss << text_ << " |";
                    for (unsigned u = 0; u < barwidth_; u++)
                        oss << symbol_;
                    oss << "| " << format(elapsed)
                        << std::string(25, ' ');
                    if (which_bar_ == &Impl::do_bar_carriage_return_)
                        oss << "\n";
                    bar_(oss.str());
                }
            }

            is_bar_print_ = false;
            state_ = total_;
        }
    };

    ///////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////  Bar : Constructor  ///////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////

    Bar::Bar()
        : pImpl_(std::make_unique<Impl>()) {}
    Bar::~Bar() = default;

    ///////////////////////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////  Bar : Setter  /////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////

    void Bar::set_text(const std::string &text) { pImpl_->set_text(text); }
    void Bar::set_thintime(unsigned long thintime) { pImpl_->set_thintime(thintime); }
    void Bar::open(unsigned long total) { pImpl_->open(total); }
    void Bar::update() { pImpl_->update(); }
    void Bar::crash() const { pImpl_->crash(); }
    void Bar::close() { pImpl_->close(); }
} // PROGRESS NAMESPACE
