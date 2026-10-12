/* housepi - Epidemiological Analysis of Household Outbreaks.
 * Copyright 2024-2026 Institut Pasteur (Paris), France.
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#pragma once

#include <algorithm>
#include <chrono>
#include <cmath>
#include <iomanip>
#include <sstream>
#include <string>
#include <utility>

#include "tools/messages.h"

//////////////////////////////  ProgressBar  //////////////////////////////

namespace housepi
{
    /// Format `secs` in whole seconds, showing large time scales only when relevant.
    /// * e.g. `12s`, `3m05s`, `1h02m`, `2d03h`.
    inline std::string format_secs(double secs)
    {
        const long t = static_cast<long>(secs);
        const long dys = t / 86400, hrs = t / 3600 % 24, mns = t / 60 % 60, scs = t % 60;
        std::ostringstream oss;
        oss << std::setfill('0');
        if (dys > 0)
            oss << dys << "d" << std::setw(2) << hrs << "h";
        else if (hrs > 0)
            oss << hrs << "h" << std::setw(2) << mns << "m";
        else if (mns > 0)
            oss << mns << "m" << std::setw(2) << scs << "s";
        else
            oss << scs << "s";
        return oss.str();
    }

    /// Display the progress of a sequence of runs, with elapsed and remaining time.
    /// Based on @cite hammouchi_fancy_2024
    ///
    /// * one thread prints, the threads running engines only count their iterations;
    /// * each run fills one unit per engine, after the units of past runs;
    /// * open line -> `text |###   | 50% 3s<3s`, closing line -> `text |######| 6s`;
    /// * closed with no line printed -> prints nothing.
    class ProgressBar
    {
        using Clock = std::chrono::steady_clock;

    public:
        //////  Constructor  //////

        /// Open a bar named `text` over `total` units, printed every `thintime` seconds.
        /// @param barwidth Symbols in a full bar.
        /// @param symbol Character filling the bar.
        ProgressBar(std::string text,
                    double total,
                    double thintime,
                    unsigned barwidth,
                    bool is_interactive,
                    std::string symbol)
            : interactive_(is_interactive),
              symbol_(std::move(symbol)),
              barwidth_(barwidth),
              text_(std::move(text)),
              total_(total),
              thintime_(thintime), next_(thintime) {}

        //////  Setter  //////

        /// Add `units` filled by a completed run.
        void advance(double units) { offset_ += units; }

        /// Print `done` units of the current run, every `thintime` seconds from opening.
        /// * prints missed by a busy thread -> skipped.
        void update(double done)
        {
            const double secs = secs_();
            if (secs < next_)
                return;
            next_ = thintime_ * (std::floor(secs / thintime_) + 1.);

            const double prog = std::min(1., (offset_ + done) / total_);
            std::ostringstream oss;
            line_(oss, prog);
            oss << std::fixed << std::setprecision(0) << 100. * prog << "% "
                << format_secs(secs) << "<"
                << (secs >= silencer_ && prog > 0. ? format_secs(secs * (1. - prog) / prog) : "?s")
                << std::setw(interactive_ ? 10 : 0) << ""; // wipe a longer previous line
            print_(oss, false);
            has_bar_print_ = true;
        }

        /// End the line of a bar left open by an error, once.
        /// * non-interactive -> no-op, every print ends its line.
        void crash()
        {
            if (has_bar_print_ && interactive_)
                msg::info("\n"); // # nocov
            has_bar_print_ = false;
        }

        /// Close the progress bar.
        void close()
        {
            // close an opened bar only
            if (!has_bar_print_)
                return;
            has_bar_print_ = false;

            std::ostringstream oss;
            line_(oss, 1.);
            oss << format_secs(secs_()) << std::setw(interactive_ ? 15 : 0) << "";
            print_(oss, true);
        }

    private:
        static constexpr const char *ansi_reset_ = "\033[39m"; /// Default foreground colour.
        static constexpr const char *ansi_green_ = "\033[32m"; /// Green foreground colour.

        const bool interactive_;                 /// Whether prints overwrite the current line.
        const std::string symbol_;               /// Character filling the bar.
        unsigned barwidth_ = 33;                 /// Bar total size.
        double silencer_ = 3.;                   /// Minimal time to wait before estimating.
        bool has_bar_print_ = false;             /// Whether the bar has printed once.
        std::string text_ = "";                  /// Text heading the bar.
        double total_ = 1.;                      /// Units of a full bar.
        double offset_ = 0.;                     /// Units filled by past runs.
        double thintime_ = 1.;                   /// Seconds between two prints.
        double next_ = 1.;                       /// Seconds from opening of the next print.
        Clock::time_point start_ = Clock::now(); /// Time of opening.

        /// Seconds since opening.
        double secs_() const { return std::chrono::duration<double>(Clock::now() - start_).count(); }

        /// Write into `os` the text and the bar filled at progress `prog`.
        /// * interactive -> led by a carriage return, bar in green.
        void line_(std::ostream &os, double prog) const
        {
            const unsigned where = prog * barwidth_;
            os << (interactive_ ? "\r" : "") << text_ << " |" << (interactive_ ? ansi_green_ : "");
            for (unsigned u = 0; u < where; ++u)
                os << symbol_;
            os << (interactive_ ? ansi_reset_ : "") << std::setw(barwidth_ - where) << "" << "| ";
        }

        /// Print line `oss`, ended on `end` or when non-interactive.
        void print_(std::ostringstream &oss, bool end) const
        {
            if (end || !interactive_)
                oss << "\n";
            msg::info(oss.str());
        }
    };
}
