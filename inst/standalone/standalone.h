/* housepi - Epidemiological Analysis of Household Outbreaks.
 *
 * Copyright © 2024-2026 Institut Pasteur (Paris), CNRS and INSERM.
 * Distributed under the terms of the GNU General Public License (GPLv3).
 * See the LICENSE file for details.
 */

#pragma once

#include <iosfwd>
#include <string>
#include <vector>

namespace housepi
{
    class Household;
    class MCMCBase;
    class MCMC;

    /// Command-line options of the standalone build.
    struct Options
    {
        bool is_log_lik = false;  /// Evaluate the log-likelihood.
        bool is_infer = false;    /// Run inference.
        bool is_simulate = false; /// Run a simulation.

        /// Paths of the input tables and of the output file.
        struct PathToFiles
        {
            std::string latent = "";
            std::string data_wide = "";
            std::string data_long = "";
            std::string group_risk = "";
            std::string contact_risk = "";
            std::string infec_profile = "";
            std::string output = "";
            std::string param = "";
            std::string setup = "";
        } file;

        /// Parse the command line.
        Options(int argc, char *argv[]);

        /// Import each table into `engine`, when its path is given.
        void read_data_wide(MCMCBase &engine);
        void read_data_long(MCMCBase &engine);
        void read_setup(MCMCBase &engine);
        void read_param(MCMCBase &engine);
        void read_infec_profile(MCMCBase &engine);
        void read_group_risk(MCMCBase &engine);
        void read_contact_risk(MCMCBase &engine);
        void read_latent(MCMCBase &engine);

        /// Assert at most one mode is selected, log-likelihood by default.
        void assert();
    };

    /// Run inference and write its trace to `f`.
    void print_infer_state(MCMC &engine, const std::string &f);

    /// Open `f` and write its header: data columns, then latent columns.
    std::ofstream open_ofstream_datawide(const std::string &f,
                                         const std::vector<std::string> &data_headers,
                                         const std::vector<std::string> &latent_headers);

    /// Simulate and write the data to `pathtofile`.
    void print_simulated_data_wide(MCMC &engine, const std::string &pathtofile);

    /// Entry point of the standalone build, for model `model_name`.
    int standalone(std::string_view model_name, int argc, char *argv[]);
}
