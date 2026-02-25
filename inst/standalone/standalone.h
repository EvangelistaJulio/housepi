/* housepi - Epidemiological Analysis of Household Outbreaks.
 *
 * Copyright © 2026 Institut Pasteur (Paris), CNRS and INSERM.
 * Distributed under the terms of the GNU General Public License (GPLv3).
 * See the LICENSE file for details.
 */

#pragma once

#include <string>
#include <vector>
#include <fstream>
#include <iostream>
#include <iomanip>

#include "tools/messages.h"

#include "mcmc.h"
#include "data_wide.h"
#include "augmented.h"

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
        std::string multiplicity = "";
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
    void read_multiplicity(MCMCBase &engine);

    bool is_log_evidence() const noexcept { return num_rungs > 0; }
    void assert();
};

/// Initialize the output file with inferred parameters with column names and initial values.
void print_infer_state(const MCMCBase &engine, const std::string &f);

/// Print all augmented variables of each individual.
std::ofstream open_ofstream_datawide(const std::string &f,
                                     const std::vector<std::string> &data_headers,
                                     const std::vector<std::string> &augm_headers);

/// Print all augmented variables of each individual.
void update_ofstream_datawide(std::ofstream &ofs,
                              Household &house,
                              size_t data_wide_size,
                              const std::vector<size_t> &augm_export);

/// Print the simulated data to the provided file.
template <class Model>
void print_simulated_data_wide(MCMC<Model> &engine, const std::string &pathtofile)
{
    const DataWide &data_wide = engine.data_wide;
    const Augmented &augm = engine.augm;

    // initialization
    std::ofstream ofs = open_ofstream_datawide(pathtofile,
                                               data_wide.pRegistry->get_all(),
                                               augm.pRegistry->as_string(augm.list_export));
    engine.helper_init_simulation();
    engine.model.init_model_for_simulation();

    for (Household house : data_wide) // NOTE: deep copy of household
        try
        {
            engine.model.simulate_single_house(house);
            // print the simulated data
            update_ofstream_datawide(ofs, house,
                                     data_wide.pRegistry->size(),
                                     augm.list_export);
        }
        catch (const std::exception &e)
        {
            throw std::runtime_error(
                "While simulating " + house.err_name() + ":\n" + e.what());
        }

    // close ofstream
    if (!ofs.is_open())
        throw std::runtime_error(
            "Output file '" + pathtofile + "' closed unexpectedly.");
    ofs.close();
}

/*
/// Estimate log_evidence according to a given scheme and number of rungs.
/// Uses trapezoidal rule to estimate logarithm of marginal likelihood (evidence).
///
/// Ref.
/// -----
/// Ref: Friel, N., Hurn, M. & Wyse, J.
/// Improving power posterior estimation of statistical evidence.
/// Stat Comput 24, 709–723 (2014).
template <class Model>
double log_evidence(MCMC<Model> &engine, std::string_view scheme, unsigned long n_rungs)
{
    double output = 0;
    const auto pp = engine.power_posteriors(scheme, n_rungs);
    const std::vector<double> &inv_tem = pp[0];
    const std::vector<double> &avg_vec = pp[1];
    const std::vector<double> &var_vec = pp[2];

    for (unsigned long k = 0; k < engine.setup.n_rungs; ++k)
    {
        double d_inv_tem = inv_tem[k + 1] - inv_tem[k];
        output += (avg_vec[k + 1] + avg_vec[k]) * d_inv_tem / 2.;
        output += (var_vec[k] - var_vec[k + 1]) * d_inv_tem * d_inv_tem / 12.;
    }

    return output;
}
*/

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////  * *** * MAIN FUNCTION * *** *  /////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/// Main function for the standalone version of 'housepi'.
template <class Model>
int standalone(int argc, char *argv[])
{
    ///////////////////////
    //////  OPTIONS  //////
    ///////////////////////

    MCMC<Model> engine = MCMC<Model>();

    Options opts(argc, argv);
    opts.read_data_wide(engine);     // 1) Load Data Wide
    opts.read_data_long(engine);     // ... and Data Long
    opts.read_setup(engine);         // 2) Load Setup
    opts.read_param(engine);         // 3) Load Parameters
    opts.read_infec_profile(engine); // 4) Load Infectivity Profile
    opts.read_group_rate(engine);    // 5) Load Group Rates
    engine.contact_rate.build();     // ... and Contact Rates UPDATE: read instead!
    opts.read_augm(engine);          // 6) Load Augmented Variables
    opts.read_multiplicity(engine);  // 7) Load Multiplicities

    //////////////////////////////
    ////////  EXECUTION  /////////
    //////////////////////////////

    try
    {
        if (opts.is_log_lik) // evaluate loglikelihood
            std::cout << engine.evaluate_log_lik() << std::endl;

        if (opts.is_simulate) // simulation an outbreak
            print_simulated_data_wide(engine, opts.file.output);

        if (opts.is_infer) // perform inference
        {
            engine.run_inference();
            print_infer_state(engine, opts.file.output);
        }

        /* --FUTURE-LOG-EVIDENCE--
        if (opts.is_log_evidence())
            std::cout << log_evidence(engine, "adaptive", num_rungs) << std::endl;
        */
    }
    catch (const std::exception &e)
    {
        engine.setup.prog_crash();
        hsp::msg::error(e.what());
    }

    return EXIT_SUCCESS;
}
