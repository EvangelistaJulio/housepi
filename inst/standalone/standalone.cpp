/* housepi - Epidemiological Analysis of Household Outbreaks.
 *
 * Copyright © 2026 Institut Pasteur (Paris), CNRS and INSERM.
 * Distributed under the terms of the GNU General Public License (GPLv3).
 * See the LICENSE file for details.
 */

#include <string>
#include <vector>
#include <unordered_set>
#include <fstream>
#include <iostream>
#include <iomanip>

#include "tools/keys.h"
#include "tools/messages.h"
#include "tools/nan.h"
#include "tools/read.h"
#include "tools/string.h"

#include "mcmc.h"
#include "setup.h"
#include "data_wide.h"
#include "data_long.h"
#include "parameter.h"
#include "group_rate.h"
#include "contact_rate.h"
#include "infec_profile.h"
#include "augmented.h"

/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////  Stdout Message Function  ////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

namespace hsp::msg
{
    /// Throw an error message.
    /// Ref:
    /// -----
    /// Suppressing Call Stack Info in Rcpp-Generated Errors and Warnings
    /// Michael Weylandt — written Mar 27, 2018
    /// https://gallery.rcpp.org/articles/quiet-stop-and-warning/
    [[noreturn]] void error(const std::string &txt)
    {
        std::cerr << "Error: " + txt + "\nExit..." << std::endl;
        std::exit(EXIT_FAILURE);
    }

    /// Print an information message to stdout.
    void info(const std::string &txt) noexcept
    {
        std::cout << txt;
        std::cout.flush();
    }

    /// Display a warning message if not already displayed before.
    /// Ref:
    /// -----
    /// Suppressing Call Stack Info in Rcpp-Generated Errors and Warnings
    /// Michael Weylandt — written Mar 27, 2018
    /// https://gallery.rcpp.org/articles/quiet-stop-and-warning/
    void warning(const std::string &txt)
    {
        if (is_new_housepi_warning(txt))
        {
            std::cerr << "Warning: " << txt << std::endl;
            add_new_housepi_warning(txt);
        }
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////  MCMC : Impl : read  ////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

#define MAKE_READ_METHOD(module, section_name)                            \
    void read_##module(MCMC &model, const std::string &f)                 \
    {                                                                     \
        try                                                               \
        {                                                                 \
            model.module.build(hsp::read::csv_table(f));                  \
        }                                                                 \
        catch (const std::exception &e)                                   \
        {                                                                 \
            hsp::msg::error(std::string("While reading ") +               \
                            section_name + " definitions:\n" + e.what()); \
        }                                                                 \
    }

MAKE_READ_METHOD(setup, "setup")
MAKE_READ_METHOD(data_wide, "data-wide")
MAKE_READ_METHOD(data_long, "data-long")
MAKE_READ_METHOD(param, "parameters")
MAKE_READ_METHOD(group_rate, "group rates")
MAKE_READ_METHOD(infec_profile, "infectivity profile")
MAKE_READ_METHOD(augm, "augmented")

// UPDATE: this will be removed once automatized in C++ reflection
void read_multiplicity(MCMC &model, const std::string &f)
{
    try
    {
        model.data_wide.multiplicity_build(hsp::read::csv_table(f),
                                           model.augm.pRegistry->size() > 0);
    }
    catch (const std::exception &e)
    {
        hsp::msg::error(std::string("While reading ") + "multiplicity" + " definitions:\n" + e.what());
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////  MCMC : print infer state  /////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/// Initialize the output file with inferred parameters with column names and initial values.
void print_infer_state(const MCMC &model, const std::string &output_filename)
{
    if (hsp::na::is(output_filename))
        return;

    const Setup &setup = model.setup;
    const Parameter &param = model.param;
    const Augmented &augm = model.augm;

    std::ofstream ofs;
    ofs.open(output_filename);

    if (!ofs.is_open())
        throw std::runtime_error(
            "Unable to open file: " + output_filename);

    // Headers
    // header log-likelihood
    ofs << "iteration" << hsp::read::DELIMITER << "log_lik";

    // header parameters
    for (const size_t p : param.list_sample)
        ofs << hsp::read::DELIMITER << param.pRegistry->as_string(p);

    // header augmented variable of id trace
    if (augm.is_any_draw())
        for (const AugmContainer *member_ptr : augm.smpl_individuals)
            for (const size_t a : augm.list_draws(*member_ptr))
                ofs << hsp::read::DELIMITER << member_ptr->id()
                    << "_" << augm.pRegistry->as_string(a);

    // header chain
    ofs << hsp::read::DELIMITER << "_chain_" << std::endl;

    // Values
    for (size_t c = 0; c < setup.n_chains; ++c)
        for (size_t k = 0; k < setup.c_sample_size; k++)
        {
            size_t i = k + c * setup.c_sample_size;

            // print iteration
            unsigned long sweep = setup.n_burnin + (k + 1) * setup.n_thin;
            ofs << sweep;

            // print log-likelihood
            ofs << hsp::read::DELIMITER << std::setprecision(2)
                << std::fixed << model.tot_log_lik.smpl_values[i];

            // print parameter values
            for (size_t u = 0; u < param.list_sample.size(); ++u)
                ofs << hsp::read::DELIMITER << std::setprecision(5)
                    << std::fixed << param.smpl_values[u][i];

            // print augmented variable of id trace
            if (!augm.smpl_values.empty())
                for (size_t j = 0; j < augm.smpl_individuals.size(); ++j)
                    for (size_t u = 0; u < augm.list_draws(*augm.smpl_individuals[j]).size(); ++u)
                        ofs << hsp::read::DELIMITER << std::setprecision(5)
                            << std::fixed << augm.smpl_values[j][u][i];

            // print which chain
            ofs << hsp::read::DELIMITER << (c + 1) << std::endl;
        }

    ofs.close();
}

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////  MCMC::Impl : print simulation  /////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/// Print all augmented variables of each individual.
std::ofstream open_ofstream_datawide(const std::string &f,
                                     const std::vector<std::string> &data_headers,
                                     const std::vector<std::string> &augm_headers)
{
    std::ofstream ofs;
    ofs.open(f);

    if (!ofs.is_open())
        throw std::runtime_error("Unable to open file: '" + f + "'");

    ofs << "id"
        << hsp::read::DELIMITER << "house_id"
        << hsp::read::DELIMITER << "status";
    for (const std::string &datastr : data_headers)
        ofs << hsp::read::DELIMITER << datastr;
    for (const std::string &augmstr : augm_headers)
        ofs << hsp::read::DELIMITER << augmstr;
    ofs << std::endl;

    return ofs;
}

/// Print all augmented variables of each individual.
void update_ofstream_datawide(std::ofstream &ofs,
                              Household &house,
                              size_t data_wide_size,
                              const std::vector<size_t> &augm_export)
{
    if (!ofs.is_open())
        throw std::runtime_error("Output file closed unexpectedly.");

    for (const Individual &member : house)
    {
        ofs << member.id()
            << hsp::read::DELIMITER << house.house_id()
            << hsp::read::DELIMITER << member.status();
        for (size_t d = 0; d < data_wide_size; ++d)
            ofs << hsp::read::DELIMITER << member.get_data_wide_vec()[d];
        for (size_t a : augm_export)
        {
            const Variable &xaugm = member.get_augm_variable(a);
            double x = xaugm.is_empty() ? hsp::na::DOUBLE : xaugm.value();
            ofs << hsp::read::DELIMITER << x;
        }
        ofs << std::endl;
    }
}

/// Print the simulated data to the provided file.
void print_simulated_data_wide(MCMC &model, const std::string &pathtofile)
{
    const DataWide &data_wide = model.data_wide;
    const Augmented &augm = model.augm;

    try
    {
        // initialization
        std::ofstream ofs = open_ofstream_datawide(pathtofile,
                                                   data_wide.pRegistry->get_all(),
                                                   augm.pRegistry->as_string(augm.list_export));
        model.helper_init_simulation();
        model.init_model_for_simulation();

        for (Household house : data_wide) // NOTE: deep copy of household
            try
            {
                model.simulate_single_house(house);
                // print the simulated data
                update_ofstream_datawide(ofs, house, data_wide.pRegistry->size(),
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
    catch (const std::exception &e)
    {
        hsp::msg::error(e.what());
    }
}

/// Estimate log_evidence according to a given scheme and number of rungs.
/// Uses trapezoidal rule to estimate logarithm of marginal likelihood (evidence).
///
/// Ref.
/// -----
/// Ref: Friel, N., Hurn, M. & Wyse, J.
/// Improving power posterior estimation of statistical evidence.
/// Stat Comput 24, 709–723 (2014).
double log_evidence(MCMC &model, std::string_view scheme, unsigned long n_rungs)
{
    try
    {
        double output = 0;
        const auto pp = model.power_posteriors(scheme, n_rungs);
        const std::vector<double> &inv_tem = pp[0];
        const std::vector<double> &avg_vec = pp[1];
        const std::vector<double> &var_vec = pp[2];

        for (unsigned long k = 0; k < model.setup.n_rungs; ++k)
        {
            double d_inv_tem = inv_tem[k + 1] - inv_tem[k];
            output += (avg_vec[k + 1] + avg_vec[k]) * d_inv_tem / 2.;
            output += (var_vec[k] - var_vec[k + 1]) * d_inv_tem * d_inv_tem / 12.;
        }

        return output;
    }
    catch (const std::exception &e)
    {
        hsp::msg::error("While estimating log_evidence with power posteriors:\n" +
                        std::string(e.what()));
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////  * *** * MAIN FUNCTION * *** *  /////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/// Main function for the standalone version of 'housepi'.
int standalone(MCMC &user_model, int argc, char *argv[])
{
    using hsp::str::tolow;
    using hsp::str::toup;

    int carg_i = 1;   // command line iterator
    std::string carg; // command line argument

    bool is_loglikelihood = false; // loglikelihood option
    bool is_infer = false;         // inference option
    bool is_simulate = false;      // simulation option
    unsigned long num_rungs = 0;   // n. rungs power posteriors

    // paths/to/files for import
    std::string augm_file = "";
    std::string data_wide_file = "";
    std::string data_long_file = "";
    std::string group_rates_file = "";
    std::string infec_profile_file = "";
    std::string multiplicity_file = "";
    std::string output_file = "";
    std::string param_file = "";
    std::string setup_file = "";

    // check any command was supplied
    if (argc < 2)
        hsp::msg::error("No argument or command found. Check the usage through:\n" +
                        std::string(argv[0]) + " -h");

    // parse commands and options
    while (carg_i < argc)
    {
        carg = std::string(argv[carg_i]); // update command line argument

        // Command line argument asking for help
        // UPDATE : should I produce a manual for the standalone?
        if (toup(carg) == "-H" || tolow(carg) == "--help")
            hsp::msg::error(
                "No Software Manual is available, refer to the Rcpp tutorials.\nApologies for the inconvenience.");

        else if (tolow(carg) == "--log-lik")
            is_loglikelihood = true;

        else if (tolow(carg) == "--infer")
            is_infer = true;

        else if (tolow(carg) == "--simulate")
            is_simulate = true;

        // UPDATE: under development
        // else if (tolow(carg) == "--log-evidence")
        //    num_rungs = hsp::str::toul(argv[++carg_i], "number of rungs");

        else if (toup(carg) == "-A" || tolow(carg) == "--augmented")
            augm_file = argv[++carg_i];

        else if (toup(carg) == "-D" || tolow(carg) == "--dataset")
            data_wide_file = argv[++carg_i];

        else if (toup(carg) == "-L" || tolow(carg) == "--datalong")
            data_long_file = argv[++carg_i];

        else if (toup(carg) == "-G" || tolow(carg) == "--group-rates")
            group_rates_file = argv[++carg_i];

        else if (toup(carg) == "-I" || tolow(carg) == "--infec-profile")
            infec_profile_file = argv[++carg_i];

        else if (toup(carg) == "-M" || tolow(carg) == "--multiplicity")
            multiplicity_file = argv[++carg_i];

        else if (toup(carg) == "-O" || tolow(carg) == "--output")
            output_file = argv[++carg_i];

        else if (toup(carg) == "-P" || tolow(carg) == "--parameter")
            param_file = argv[++carg_i];

        else if (toup(carg) == "-S" || tolow(carg) == "--setup")
            setup_file = argv[++carg_i];

        else
            hsp::msg::error("Unrecognized option `" + carg + "`.");

        carg_i++; // move to the next command line argument
    }

    // check for single option
    bool is_logevidence = num_rungs > 0; // power posterior option
    int aux = is_loglikelihood + is_infer + is_simulate + is_logevidence;
    if (aux > 1)
        hsp::msg::error(
            "Incompatible commands: choose one between --log-lik, --infer, --simulate or --log-evidence.");
    if (aux == 0)
        is_loglikelihood = true; // default

    //////////////////////////////
    //////  INITIALIZATION  //////
    //////////////////////////////

    read_data_wide(user_model, data_wide_file);         // 1) Load Data Wide
    read_data_long(user_model, data_long_file);         // ... and Data Long
    read_setup(user_model, setup_file);                 // 2) Load Setup
    read_param(user_model, param_file);                 // 3) Load Parameters
    read_infec_profile(user_model, infec_profile_file); // 4) Load Infectivity Profile
    read_group_rate(user_model, group_rates_file);      // 5) Load Group Rates
    user_model.contact_rate.build();                    // ... and Contact Rates UPDATE: read instead!
    read_augm(user_model, augm_file);                   // 6) Load Augmented Variables
    read_multiplicity(user_model, multiplicity_file);   // 7) Load Multiplicities

    //////////////////////////////
    ////////  EXECUTION  /////////
    //////////////////////////////

    if (is_loglikelihood) // evaluate loglikelihood
    {
        std::cout << user_model.evaluate_log_lik() << std::endl;
    }
    else if (is_infer) // perform inference
    {
        if (hsp::na::is(output_file))
            hsp::msg::error("Missing output file for inference.");
        user_model.run_inference();
        print_infer_state(user_model, output_file);
    }
    else if (is_simulate) // simulation an outbreak
    {
        if (hsp::na::is(output_file))
            hsp::msg::error("Missing output file for simulation.");
        print_simulated_data_wide(user_model, output_file);
    }
    else if (is_logevidence) // estimate log_evidence
    {
        std::cout << log_evidence(user_model, "adaptive", num_rungs) << std::endl;
    }

    return EXIT_SUCCESS;
}
