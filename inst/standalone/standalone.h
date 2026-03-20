/* housepi - Epidemiological Analysis of Household Outbreaks.
 *
 * Copyright © 2026 Institut Pasteur (Paris), CNRS and INSERM.
 * Distributed under the terms of the GNU General Public License (GPLv3).
 * See the LICENSE file for details.
 */

#pragma once

#include <string>
#include <vector>

namespace hsp
{
    class MCMCBase;
    class MCMC;

    struct Options
    {
        bool is_log_lik = false;     // loglikelihood option
        bool is_infer = false;       // inference option
        bool is_simulate = false;    // simulation option
        unsigned long num_rungs = 0; // n. rungs power posteriors

        struct PathToFiles
        {
            std::string augm = "";
            std::string data_wide = "";
            std::string data_long = "";
            std::string group_rate = "";
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
        void read_group_rate(MCMCBase &engine);
        void read_augm(MCMCBase &engine);

        bool is_log_evidence() const noexcept { return num_rungs > 0; }
        void assert();
    };

    /// Initialize the output file with inferred parameters with column names and initial values.
    void print_infer_state(MCMC &engine, const std::string &f);

    /// Print all augmented variables of each individual.
    std::ofstream open_ofstream_datawide(const std::string &f,
                                         const std::vector<std::string> &data_headers,
                                         const std::vector<std::string> &augm_headers);

    /// Print all augmented variables of each individual.
    void update_ofstream_datawide(std::ofstream &ofs,
                                  const Household &house,
                                  size_t data_wide_size,
                                  const std::vector<size_t> &augm_export);

    /// Print the simulated data to the provided file.
    void print_simulated_data_wide(MCMC &engine, const std::string &pathtofile);

    /* --FUTURE-LOG-EVIDENCE--
    /// Estimate log_evidence according to a given scheme and number of rungs.
    /// Uses trapezoidal rule to estimate logarithm of marginal likelihood (evidence).
    ///
    /// Ref.
    /// -----
    /// Ref: Friel, N., Hurn, M. & Wyse, J.
    /// Improving power posterior estimation of statistical evidence.
    /// Stat Comput 24, 709–723 (2014).
    double log_evidence(MCMC &engine, std::string_view scheme, unsigned long n_rungs);
    */

    /// Main function for the standalone version of 'housepi'.
    int standalone(MCMC &engine, int argc, char *argv[]);
}
