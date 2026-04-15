/* housepi - Epidemiological Analysis of Household Outbreaks.
 *
 * Copyright © 2026 Institut Pasteur (Paris), CNRS and INSERM.
 * Distributed under the terms of the GNU General Public License (GPLv3).
 * See the LICENSE file for details.
 */

#include <fstream>
#include <functional>
#include <iomanip>
#include <sstream>

#include <cstdint>
#include <chrono>
#include <memory>
#include <random>
#include <stdexcept>
#include <string>
#include <vector>

#include "tools/messages.h"
#include "tools/nan.h"
#include "tools/openmp.h"
#include "tools/read.h"
#include "tools/string.h"
#include "tools/types.h"

#include "setup.h"

/////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////  Progress Bar ////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////

namespace housepi
{
    using namespace std::chrono;

    /// Format elapsed time in seconds.
    std::string format_time_(seconds scs)
    {
        std::ostringstream oss;

        const auto dys = duration_cast<duration<int, std::ratio<86400>>>(scs);
        scs -= duration_cast<seconds>(dys);

        const auto hrs = duration_cast<hours>(scs);
        scs -= duration_cast<seconds>(hrs);

        const auto mns = duration_cast<minutes>(scs);
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
    ////////////////////////////////// *** ProgressBar : class *** //////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////

    namespace ansi
    {
        constexpr const char *reset = "\033[39m";
        constexpr const char *green = "\033[32m";
    }

    /// Implementation class for progress bar.
    ///
    /// -----
    /// Adapted from :
    /// OTHMAN EL HAMMOUCHI, April 15, 2024
    /// "Fancy thread safe progress bar in R for compiled code"
    /// https://othmanelhammouchi.com/posts/r-parallel-progress/
    class ProgressBar
    {
        enum class Mode : std::uint8_t
        {
            Interactive,
            Static
        } mode_ = Mode::Interactive;

    public:
        ///////////////////////////
        //////  CONSTRUCTOR  //////
        ///////////////////////////

        ProgressBar(unsigned barwidth = 33, bool is_interactive = true,
                    std::function<void()> interrupt_func = {})
            : barwidth_(barwidth), interrupt_func_(std::move(interrupt_func))
        {
            if (!is_interactive)
                mode_ = Mode::Static;
        };
        ~ProgressBar() = default;

        //////////////////////
        //////  SETTER  //////
        //////////////////////

        /// Set description text of progress bar.
        void set_text(const std::string &text) { text_ = text; }

        /// Start progressbar clock and display firs message.
        void open(unsigned long total, unsigned long thintime)
        {
            check_user_interrupt_();
            if (total == 0)
                throw std::invalid_argument(
                    "Progressbar total should be greater than 0.");

            // initialize internal fields
            total_ = total;
            start_ = steady_clock::now();
            lastprint_ = start_;
            state_ = 0;
            has_bar_print_ = false;
            verbose_ = thintime > 0;
            thintime_ = thintime - 1;
        }

        /// Update progress bar given the progress state.
        void update()
        {
            check_user_interrupt_();
            ++state_; // update internal state
            if (!verbose_)
                return;

            const double progress = static_cast<double>(state_) / static_cast<double>(total_);
            const auto now = steady_clock::now();
            const auto elapsed = duration_cast<seconds>(now - start_);
            const auto missing = duration_cast<seconds>(elapsed * (1. - progress) / progress);

            // skip if expected to finish soon (less than 1 second)
            if (missing.count() < 1)
                return;

            // print update message if time since last print greater than thin-time
            if (static_cast<unsigned long>(duration_cast<seconds>(now - lastprint_).count()) > thintime_)
            {
                std::ostringstream oss;
                oss << text_
                    << bar_(progress)
                    << std::fixed << std::setprecision(0) << 100. * progress << "% "
                    << format_time_(elapsed) << "<";
                if (elapsed.count() > 2)
                    oss << format_time_(missing);
                else 
                    oss << "?s";
                oss << std::string(10, ' ');
                print_(oss.str());

                lastprint_ = steady_clock::now();
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

                // NOTE: initial message must not print '\r', eventually '\n'
                if (mode_ == Mode::Static)
                    oss << "\n";
                msg::info(oss.str());

                has_bar_print_ = true;
            }
        }

        /// End the progress bar because of a crash.
        void crash() const
        {
            check_user_interrupt_();
            if (has_bar_print_)
                msg::info("\n");
        }

        /// Close the progress bar.
        void close()
        {
            check_user_interrupt_();

            // close only if it was opened
            if (!has_bar_print_)
                return;
            has_bar_print_ = false;

            if (!verbose_)
                return;

            const auto elapsed = duration_cast<seconds>(steady_clock::now() - start_);
            if (static_cast<unsigned long>(elapsed.count()) > thintime_)
            {
                std::ostringstream oss;
                oss << text_
                    << bar_(1)
                    << format_time_(elapsed)
                    << std::string(25, ' ');
                if (mode_ == Mode::Interactive)
                    oss << "\n";
                print_(oss.str());
            }
        }

    private:
        const std::string symbol_ = "\u25A0";  /// bar decorator
        unsigned barwidth_ = 33;               /// bar total size
        std::function<void()> interrupt_func_; /// user-interrupt callback
        bool verbose_ = true;                  /// wheter to show bar
        bool has_bar_print_ = false;           /// wheter bar has printed once
        unsigned long state_ = 0;              /// current iteration
        unsigned long total_ = 0;              /// total number of iterations
        std::string text_ = "";                /// current text to show

        steady_clock::time_point start_;
        steady_clock::time_point lastprint_;
        unsigned long thintime_ = 0;

        /// Colored progress bar.
        std::string bar_(double progress) const
        {
            const unsigned where = progress * barwidth_;
            std::ostringstream oss;
            oss << " |";
            if (mode_ == Mode::Interactive)
                oss << ansi::green;
            for (unsigned u = 0; u < where; ++u)
                oss << symbol_;
            if (mode_ == Mode::Interactive)
                oss << ansi::reset;
            oss << std::string(barwidth_ - where, ' ') << "| ";
            return oss.str();
        }

        /// Function to display the progress bar.
        void print_(const std::string &sms) const
        {
            switch (mode_)
            {
            case Mode::Interactive:
                msg::info("\r" + sms);
                return;
            case Mode::Static:
                msg::info(sms + "\n");
                return;
            default:
                __builtin_unreachable();
            }
        }

        void check_user_interrupt_() const
        {
            if (interrupt_func_)
                interrupt_func_();
        }
    };

    /////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////  *** Setup : struct ***  ////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////

    Setup::Setup()
        : uniform_(0., 1.),
          progressbar_(std::make_unique<ProgressBar>()) {}
    Setup::~Setup() = default;

    void Setup::init_R_context(
        unsigned barwidth,
        bool is_interactive,
        std::function<void()> interrupt_func)
    {
        progressbar_ = std::make_unique<ProgressBar>(
            barwidth,
            is_interactive,
            std::move(interrupt_func));
    }

    void Setup::build(const str_table &input)
    {
        const auto time_raw = std::chrono::steady_clock::now().time_since_epoch().count();
        // ensure time_now is 64bit
        const auto time_now = static_cast<std::uint64_t>(time_raw);
        // xor between upper and lower part of time stamp
        const std::uint64_t mixed = time_now ^ (time_now >> 32);
        // default seed for pseudo random number generator
        // note: constrain to 32bit for compatibility with R (mask sign bit)
        seed = static_cast<std::int32_t>(mixed & 0x7FFFFFFF);

        chain = 1;                        // current chain
        rung = 0;                         // current thermodynamic rung
        id_draws.clear();                 // drop previous id
        helper_non_default_build_(input); // read user input
        set_n_threads();                  // define number of threads
        set_prng_seed();                  // pseudo-random number generator seeding
    }

    void Setup::helper_non_default_build_(const housepi::str_table &input)
    {
        if (input.empty())
            return;

        // check no multiple entries if internal error
        if (read::csv_n_rows(input) != 1)
            throw std::runtime_error(
                "Internal error: multiple entries found in MCMC setup definition.");

        // check entries
        const std::vector<std::string> &headers = read::extract_keys(input);
        str::assert_no_unused_entries(headers, file_headers);

        // Random Initialization of Parameters
        if (input.count("progress") && !na::is(input.at("progress")[0]))
            progress = str::toul(input.at("progress")[0]);

        // N. of Threads
        if (input.count("threads") && !na::is(input.at("threads")[0]))
            n_threads = str::toi(input.at("threads")[0], "threads");

        // PRNG Seed
        if (input.count("seed") && !na::is(input.at("seed")[0]))
            seed = str::toi32(input.at("seed")[0], "seed");

        // N. of Chains
        if (input.count("chains") && !na::is(input.at("chains")[0]))
            n_chains = str::toul(input.at("chains")[0], "chains");

        // N. of Steps
        if (input.count("steps") && !na::is(input.at("steps")[0]))
            n_steps = str::toul(input.at("steps")[0], "steps");

        // N. of Adaptive Steps
        if (input.count("adapt") && !na::is(input.at("adapt")[0]))
            n_adapt = str::toul(input.at("adapt")[0], "adpat");

        // N. of Thinning Steps
        if (input.count("thin") && !na::is(input.at("thin")[0]))
            n_thin = str::toul(input.at("thin")[0], "thin");

        // N. of Steps to "Burn-In"
        if (input.count("burn") && !na::is(input.at("burn")[0]))
            n_burn = str::toul(input.at("burn")[0], "burn");

        // Individual ID for Augmented Trace
        if (input.count("id_draws") && !na::is(input.at("id_draws")[0]))
            id_draws = read::parse_list(input.at("id_draws")[0]);
    }

    void Setup::set_n_threads()
    {
        if (n_threads == 0)
            throw std::invalid_argument(
                "N. of threads was not greater than 0.");

#ifdef _OPENMP
        unsigned int max_threads = omp_get_num_procs();

        if (n_threads > max_threads)
        {
            msg::warning(
                "N. of threads reduced to " + std::to_string(max_threads) + " available threads.");
            n_threads = max_threads;
        }

        // declare the global number of omp threads
        omp_set_num_threads(n_threads);

#else
        if (n_threads != 1)
            msg::warning(
                "Multi-threading was disabled. OpenMP is not supported on current machine.");

        n_threads = 1;
#endif
    }

    /// Assign seed for random generators.
    /// Concomitantely create different generators for different parallel threads.
    void Setup::set_prng_seed()
    {
        prng_parallel_.resize(n_threads);
        prng_parallel_[0].seed(seed);
        for (size_t t = 1; t < n_threads; ++t)
        {
            prng_parallel_[t] = prng_parallel_[t - 1];
            prng_parallel_[t].jump();
        }
    }

    /// Return pseudo random number generator.
    rnd_gen &Setup::prng() const
    {
#ifdef _OPENMP
        return prng_parallel_[omp_get_thread_num()];
#else
        return prng_parallel_[0];
#endif
    }

    void Setup::assert_mcmc_setup()
    {
        if (inv_T < 0. || inv_T > 1.)
            throw std::out_of_range(
                "Inverse temperature " + std::to_string(inv_T) + " was out of range [0, 1].");

        if (n_chains < 1)
            throw std::out_of_range(
                "Number of chains was not greater than 0.");

        // Interpret thinning equal to 0 has been interpreted as 1.
        if (n_thin == 0)
            n_thin = 1;

        if (n_steps < 1)
            throw std::out_of_range(
                "N. of steps was not greater than 0.");
        if (n_burn >= n_steps)
            throw std::out_of_range(
                "Burn-in must was not smaller than n. of steps.");
        if (n_thin > n_steps)
            throw std::out_of_range(
                "Thinning was not be smaller than (or equal to) n. of steps. ");

        c_sample_size = (n_steps - n_burn) / n_thin;
        // adjust steps and sample size
        c_sample_size += (n_steps - n_burn) % n_thin > 0;
        n_steps = n_burn + c_sample_size * n_thin;

        if (c_sample_size == 0)
            throw std::out_of_range(
                "Choice of steps, burn-in and thin generated an empty sample.");

        // assign number of sweeps
        // NOTE: because of rounding I can't put all in the same line
        c_n_sweeps = n_adapt / n_thin;
        c_n_sweeps += n_burn / n_thin;
        c_n_sweeps += c_sample_size;
    }

    // Drop adaptation steps.
    void Setup::drop_mcmc_adapt()
    {
        if (n_adapt > 0)
            msg::warning(
                "Dropped adaptation steps: no parametric MCMC proposals were found.");
        n_adapt = 0;

        // re-assign number of sweeps
        c_n_sweeps = n_burn / n_thin;
        c_n_sweeps += c_sample_size;
    }

    void Setup::prog_set_text(std::string text)
    {
        if (!na::is(text))
            ;
        else if (is_power_post())
            text = "power post " + std::to_string(rung) + "/" + std::to_string(n_rungs);
        else if (n_chains > 1)
            text = "sample " + std::to_string(chain) + "/" + std::to_string(n_chains);
        else
            text = "sample";
        progressbar_->set_text(text);
    }
    void Setup::prog_open_for_mcmc_sample()
    {
        prog_set_text();
        progressbar_->open(c_n_sweeps * n_chains, progress);
    }
    void Setup::prog_open_for_power_post()
    {
        prog_set_text();
        progressbar_->open(c_n_sweeps * n_chains * (n_rungs + 1), progress);
    }
    void Setup::prog_open_for_simulation()
    {
        prog_set_text("simulate");
        progressbar_->open(n_steps, progress);
    }

    void Setup::prog_update() { progressbar_->update(); }
    void Setup::prog_close() { progressbar_->close(); }
    void Setup::crash(const std::exception &e)
    {
        progressbar_->crash();
        msg::error(e.what());
    }
}
