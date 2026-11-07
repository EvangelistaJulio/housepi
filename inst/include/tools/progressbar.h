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
    /// Format `secs` rounded up to the smallest time scale shown.
    /// * large time scales shown only when relevant, e.g. `12s`, `3m05s`, `1h02m`, `2d03h`;
    /// * not finite -> `?s`.
    inline std::string format_secs(double secs)
    {
        if (!std::isfinite(secs))
            return "?s";
        long t = static_cast<long>(std::ceil(secs));
        // seconds hidden, round up to minutes
        if (t >= 3600)
            t = (t + 59) / 60 * 60;
        // minutes hidden, round up to hours
        if (t >= 86400)
            t = (t + 3599) / 3600 * 3600;
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
    /// * one thread prints, the threads running engines only count their steps;
    /// * each run fills one unit per engine, after the units of past runs;
    /// * open line -> `text |###   | 50% 3s<3s`, closing line -> `text |######| 6s`;
    /// * document -> closing line alone, never rewritten;
    /// * closed before its first print time -> prints nothing.
    class ProgressBar
    {
        using Clock = std::chrono::steady_clock;

    public:
        //////  Constructor  //////

        /// Open a bar named `text` over `total` units, printed every `thintime` seconds.
        /// @param barwidth Symbols in a full bar.
        /// @param is_document Whether the output is a knitted document, which overrides `is_interactive`.
        /// @param is_coloured Whether the output renders ANSI colours.
        /// @param symbol Character filling the bar.
        ProgressBar(std::string text,
                    double total,
                    double thintime,
                    unsigned barwidth,
                    bool is_interactive,
                    bool is_document,
                    bool is_coloured,
                    std::string symbol)
            : interactive_(is_interactive && !is_document),
              document_(is_document),
              coloured_(is_coloured),
              symbol_(std::move(symbol)),
              text_(std::move(text)),
              barwidth_(barwidth),
              total_(total),
              thintime_(thintime), next_(thintime) {}

        //////  Setter  //////

        /// Add `units` filled by a completed run.
        void advance(double units) { offset_ += units; }

        /// Print `done` units of the current run, every `thintime` seconds from opening.
        /// * prints missed by a busy thread -> skipped;
        /// * document -> no print, the closing line held due.
        void update(double done)
        {
            const double secs = secs_();
            if (secs < next_)
                return;
            next_ = thintime_ * (std::floor(secs / thintime_) + 1.);
            has_line_ = true;
            if (document_)
                return;

            const double prog = std::min(1., (offset_ + done) / total_);
            std::ostringstream oss;
            line_(oss, prog);
            oss << std::fixed << std::setprecision(0) << 100. * prog << "% "
                << format_secs(secs) << "<"
                << (secs >= silencer_ && prog > 0. ? format_secs(secs * (1. - prog) / prog) : "?s")
                // wipe a longer previous line
                << std::setw(interactive_ ? 10 : 0) << "";
            print_(oss, false);
        }

        /// End the line of a bar left open by an error, once.
        /// * non-interactive -> no-op, every print ends its line.
        void crash()
        {
            if (has_line_ && interactive_)
                msg::info("\n"); // # nocov
            has_line_ = false;
        }

        /// Close the progress bar.
        void close()
        {
            // close a bar past its first print time only
            if (!has_line_)
                return;
            has_line_ = false;

            std::ostringstream oss;
            line_(oss, 1.);
            oss << format_secs(secs_()) << std::setw(interactive_ ? 15 : 0) << "";
            print_(oss, true);
        }

    private:
        /// Default foreground colour.
        static constexpr const char *ANSI_RESET = "\033[39m";
        /// Green foreground colour.
        static constexpr const char *ANSI_GREEN = "\033[32m";

        /// Whether prints overwrite the current line.
        const bool interactive_;
        /// Whether only the closing line prints.
        const bool document_;
        /// Whether the bar is filled in green.
        const bool coloured_;
        /// Character filling the bar.
        const std::string symbol_;
        /// Text heading the bar.
        const std::string text_;
        /// Bar total size.
        const unsigned barwidth_ = 33;
        /// Minimal time to wait before estimating.
        const double silencer_ = 3.;
        /// Units of a full bar.
        const double total_ = 1.;
        /// Seconds between two prints.
        const double thintime_ = 1.;

        /// Units filled by past runs.
        double offset_ = 0.;
        /// Seconds from opening of the next print.
        double next_ = 1.;
        /// Whether the bar is past its first print time.
        bool has_line_ = false;
        /// Time of opening.
        Clock::time_point start_ = Clock::now();

        /// Seconds since opening.
        double secs_() const { return std::chrono::duration<double>(Clock::now() - start_).count(); }

        /// Write into `os` the text and the bar filled at progress `prog`.
        /// * interactive -> led by a carriage return;
        /// * coloured -> bar in green.
        void line_(std::ostream &os, double prog) const
        {
            const unsigned where = prog * barwidth_;
            os << (interactive_ ? "\r" : "") << text_ << " |" << (coloured_ ? ANSI_GREEN : "");
            for (unsigned u = 0; u < where; ++u)
                os << symbol_;
            os << (coloured_ ? ANSI_RESET : "") << std::setw(barwidth_ - where) << "" << "| ";
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
