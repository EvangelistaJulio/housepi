/* housepi - Epidemiological Analysis of Household Outbreaks.
 *
 * Copyright © 2024-2026 Institut Pasteur (Paris), CNRS and INSERM.
 * Distributed under the terms of the GNU General Public License (GPLv3).
 * See the LICENSE file for details.
 */

#pragma once

#include <string>
#include <vector>

namespace housepi
{
    class MCMCBase;
    class MCMC;

    struct Options
    {
        bool is_log_lik = false;     // loglikelihood option
        bool is_infer = false;       // inference option
        bool is_simulate = false;    // simulation option

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

        Options(int argc, char *argv[]);

        void read_data_wide(MCMCBase &engine);
        void read_data_long(MCMCBase &engine);
        void read_setup(MCMCBase &engine);
        void read_param(MCMCBase &engine);
        void read_infec_profile(MCMCBase &engine);
        void read_group_risk(MCMCBase &engine);
        void read_contact_risk(MCMCBase &engine);
        void read_latent(MCMCBase &engine);

        void assert();
    };

    /// Initialize the output file with inferred parameters with column names and initial values.
    void print_infer_state(MCMC &engine, const std::string &f);

    /// Print all latent variables of each individual.
    std::ofstream open_ofstream_datawide(const std::string &f,
                                         const std::vector<std::string> &data_headers,
                                         const std::vector<std::string> &latent_headers);

    /// Print all latent variables of each individual.
    void update_ofstream_datawide(std::ofstream &ofs,
                                  const Household &house,
                                  size_t data_wide_size,
                                  const std::vector<size_t> &latent_export);

    /// Print the simulated data to the provided file.
    void print_simulated_data_wide(MCMC &engine, const std::string &pathtofile);

    /// Main function for the standalone version of 'housepi'.
    int standalone(std::string_view model_name, int argc, char *argv[]);
}
