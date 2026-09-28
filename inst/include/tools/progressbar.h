/* housepi - Epidemiological Analysis of Household Outbreaks.
 *
 * Copyright © 2024-2026 Institut Pasteur (Paris), CNRS and INSERM.
 * Distributed under the terms of the GNU General Public License (GPLv3).
 * See the LICENSE file for details.
 */

#pragma once

#include <chrono>
#include <cstdint>
#include <functional>
#include <iomanip>
#include <sstream>
#include <string>

#include "tools/messages.h"

//////////////////////////////  ProgressBar  //////////////////////////////

namespace housepi
{
    /// Display the progress of an iteration count, with elapsed and remaining time.
    /// Based on @cite hammouchi_fancy_2024
    class ProgressBar
    {
        enum class Mode : std::uint8_t
        {
            Interactive,
            Static
        } mode_ = Mode::Interactive;

        using Clock = std::chrono::steady_clock;
        using Seconds = std::chrono::seconds;

    public:
        //////  Constructor  //////

        ProgressBar(
            unsigned barwidth = 33,
            bool is_interactive = true,
            std::function<void()> interrupt_func = [] {})
            : barwidth_(barwidth), interrupt_func_(std::move(interrupt_func))
        {
            if (!is_interactive)
                mode_ = Mode::Static;
        };
        ~ProgressBar() = default;

        //////  Setter  //////

        /// Set description text of progress bar.
        void set_text(const std::string &text) { text_ = text; }

        /// Start progress bar clock and display first message.
        void open(unsigned long total, unsigned long thintime)
        {
            check_user_interrupt_();

            // initialize internal fields
            total_ = total;
            start_ = Clock::now();
            lastprint_ = start_;
            state_ = 0;
            has_bar_print_ = false;
            verbose_ = thintime > 0;
            thintime_ = thintime - 1;
        }

        /// Update progress bar given the progress state.
        void update()
        {
            using std::chrono::duration_cast;

            check_user_interrupt_();
            ++state_; // update internal state
            if (!verbose_)
                return;

            const double progress = static_cast<double>(state_) / static_cast<double>(total_);
            const auto now = Clock::now();
            const auto elapsed = duration_cast<Seconds>(now - start_);
            const auto missing = duration_cast<Seconds>(elapsed * (1. - progress) / progress);

            // skip if expected to finish soon (less than 1 second)
            if (missing.count() < 1)
                return;

            // print update message if time since last print greater than thin-time
            if (static_cast<unsigned long>(duration_cast<Seconds>(now - lastprint_).count()) > thintime_)
            {
                std::ostringstream oss;
                oss << text_
                    << bar_(progress)
                    << std::fixed << std::setprecision(0) << 100. * progress << "% "
                    << format_time_(elapsed) << "<"
                    << ((elapsed.count() > 2) ? format_time_(missing) : "?s");
                oss << std::string(10, ' ');
                print_(oss.str());

                lastprint_ = Clock::now();
                has_bar_print_ = true;
            }
            else if (!has_bar_print_ && static_cast<unsigned long>(missing.count()) > thintime_)
            {
                std::ostringstream oss;
                oss << text_
                    << bar_(0.)
                    << std::fixed << std::setprecision(0) << 0. << "% "
                    << "0s<?s"
                    << std::string(10, ' ');

                // initial message prints no '\r', a '\n' in static mode
                if (mode_ == Mode::Static)
                    oss << "\n";
                msg::info(oss.str());

                has_bar_print_ = true;
            }
        }

        /// End the progress bar because of a crash.
        void crash() const
        {
            if (has_bar_print_)
                msg::info("\n");
        }

        /// Close the progress bar, returning the seconds since `open()`.
        double close()
        {
            const auto now = Clock::now();
            const double spent = std::chrono::duration<double>(now - start_).count();

            // close an opened bar only
            if (!has_bar_print_)
                return spent;
            has_bar_print_ = false;

            const auto elapsed = std::chrono::duration_cast<Seconds>(now - start_);
            if (static_cast<unsigned long>(elapsed.count()) > thintime_)
            {
                std::ostringstream oss;
                oss << text_
                    << bar_(1)
                    << format_time_(elapsed)
                    << std::string(25, ' ');
                // # nocov start
                if (mode_ == Mode::Interactive)
                    oss << "\n";
                // # nocov end
                print_(oss.str());
            }

            check_user_interrupt_();
            return spent;
        }

    private:
        static constexpr const char *ansi_reset_ = "\033[39m"; /// Default foreground colour.
        static constexpr const char *ansi_green_ = "\033[32m"; /// Green foreground colour.

        const std::string symbol_ = "■";       /// Bar decorator.
        unsigned barwidth_ = 33;               /// Bar total size.
        std::function<void()> interrupt_func_; /// User-interrupt callback.
        bool verbose_ = true;                  /// Whether to show the bar.
        bool has_bar_print_ = false;           /// Whether the bar has printed once.
        unsigned long state_ = 0;              /// Current iteration.
        unsigned long total_ = 0;              /// Total number of iterations.
        std::string text_ = "";                /// Current text to show.
        Clock::time_point start_;              /// Time of `open()`.
        Clock::time_point lastprint_;          /// Time of the last print.
        unsigned long thintime_ = 0;           /// Seconds between two prints.

        /// Format elapsed time `scs`, showing large time scales only when relevant.
        static std::string format_time_(Seconds scs)
        {
            using namespace std::chrono;
            std::ostringstream oss;

            const auto dys = duration_cast<duration<int, std::ratio<86400>>>(scs);
            scs -= duration_cast<seconds>(dys);

            const auto hrs = duration_cast<hours>(scs);
            scs -= duration_cast<seconds>(hrs);

            const auto mns = duration_cast<minutes>(scs);
            scs -= duration_cast<seconds>(mns);

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

        /// Colored progress bar.
        std::string bar_(double progress) const
        {
            const unsigned where = progress * barwidth_;
            std::ostringstream oss;
            oss << " |";
            // # nocov start
            if (mode_ == Mode::Interactive)
                oss << ansi_green_;
            // # nocov end
            for (unsigned u = 0; u < where; ++u)
                oss << symbol_;
            // # nocov start
            if (mode_ == Mode::Interactive)
                oss << ansi_reset_;
            // # nocov end
            oss << std::string(barwidth_ - where, ' ') << "| ";
            return oss.str();
        }

        /// Check for user interrupt.
        void check_user_interrupt_() const { interrupt_func_(); }

        /// Print message `sms`, overwriting the current line in interactive mode.
        void print_(const std::string &sms) const
        {
            switch (mode_)
            {
            case Mode::Static:
                msg::info(sms + "\n");
                return;
            default:
                // # nocov start
                msg::info("\r" + sms);
                return;
                // # nocov end
            }
        }
    };
}
