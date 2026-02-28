/* housepi - Epidemiological Analysis of Household Outbreaks.
 *
 * Copyright © 2026 Institut Pasteur (Paris), CNRS and INSERM.
 * Distributed under the terms of the GNU General Public License (GPLv3).
 * See the LICENSE file for details.
 */

#include <string>
#include <vector>
#include <fstream>
#include <iostream>
#include <iomanip>

#include "tools/keys.h"
#include "tools/messages.h"
#include "tools/nan.h"
#include "tools/read.h"
#include "tools/string.h"

#include "mcmc.h"
#include "standalone.h"
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

#define STANDALONE_MAKE_READ_METHOD_(module, section_name)                \
    void Options::read_##module(MCMCBase &model)                          \
    {                                                                     \
        try                                                               \
        {                                                                 \
            model.module.build(hsp::read::csv_table(file.module));        \
        }                                                                 \
        catch (const std::exception &e)                                   \
        {                                                                 \
            hsp::msg::error(std::string("While reading ") +               \
                            section_name + " definitions:\n" + e.what()); \
        }                                                                 \
    }

STANDALONE_MAKE_READ_METHOD_(setup, "setup")
STANDALONE_MAKE_READ_METHOD_(data_wide, "data-wide")
STANDALONE_MAKE_READ_METHOD_(data_long, "data-long")
STANDALONE_MAKE_READ_METHOD_(param, "parameters")
STANDALONE_MAKE_READ_METHOD_(group_rate, "group rate")
STANDALONE_MAKE_READ_METHOD_(infec_profile, "infectivity profile")
STANDALONE_MAKE_READ_METHOD_(augm, "augmented")

// UPDATE: this will be removed once automatized in C++ reflection
void Options::read_multiplicity(MCMCBase &model)
{
    try
    {
        model.data_wide.multiplicity_build(
            hsp::read::csv_table(file.multiplicity),
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
void print_infer_state(const MCMCBase &engine, const std::string &out_file)
{
    if (hsp::na::is(out_file))
        return;

    const Setup &setup = engine.setup;
    const Parameter &param = engine.param;
    const Augmented &augm = engine.augm;

    std::ofstream ofs;
    ofs.open(out_file);

    if (!ofs.is_open())
        throw std::runtime_error(
            "Unable to open file: " + out_file);

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
                << std::fixed << engine.tot_log_lik.smpl_values[i];

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

void print_simulated_data_wide(MCMC &engine, const std::string &pathtofile)
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

Options::Options(int argc, char *argv[])
{
    using hsp::str::tolow;
    using hsp::str::toup;

    int carg_i = 1;   // command line iterator
    std::string carg; // command line argument

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
            is_log_lik = true;

        else if (tolow(carg) == "--infer")
            is_infer = true;

        else if (tolow(carg) == "--simulate")
            is_simulate = true;

        /*  --FUTURE-LOG-EVIDENCE--
        else if (tolow(carg) == "--log-evidence")
            num_rungs = hsp::str::toul(argv[++carg_i], "number of rungs");
        */

        else if (toup(carg) == "-A" || tolow(carg) == "--augmented")
            file.augm = argv[++carg_i];

        else if (toup(carg) == "-D" || tolow(carg) == "--dataset")
            file.data_wide = argv[++carg_i];

        else if (toup(carg) == "-L" || tolow(carg) == "--datalong")
            file.data_long = argv[++carg_i];

        else if (toup(carg) == "-G" || tolow(carg) == "--group-rates")
            file.group_rate = argv[++carg_i];

        else if (toup(carg) == "-I" || tolow(carg) == "--infec-profile")
            file.infec_profile = argv[++carg_i];

        else if (toup(carg) == "-M" || tolow(carg) == "--multiplicity")
            file.multiplicity = argv[++carg_i];

        else if (toup(carg) == "-O" || tolow(carg) == "--output")
            file.output = argv[++carg_i];

        else if (toup(carg) == "-P" || tolow(carg) == "--parameter")
            file.param = argv[++carg_i];

        else if (toup(carg) == "-S" || tolow(carg) == "--setup")
            file.setup = argv[++carg_i];

        else
            hsp::msg::error("Unrecognized option `" + carg + "`.");

        carg_i++; // move to the next command line argument
    }
    assert();
}

void Options::assert()
{
    // check for single option
    int aux = is_log_lik + is_infer + is_simulate + is_log_evidence();

    if (aux == 0)
        is_log_lik = true; // default
    if (aux > 1)
        hsp::msg::error(
            "Incompatible commands: choose one between --log-lik, --infer, --simulate or --log-evidence.");

    if (!hsp::na::is(file.output))
        return;
    if (is_infer)
        hsp::msg::error("Missing output file for inference.");
    if (is_simulate)
        hsp::msg::error("Missing output file for simulation.");
}

/* --FUTURE-LOG-EVIDENCE--
double log_evidence(MCMC &engine, std::string_view scheme, unsigned long n_rungs)
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

/// Main function for the standalone version of 'housepi'.
int standalone(MCMC &engine, int argc, char *argv[])
{
    ///////////////////////
    //////  OPTIONS  //////
    ///////////////////////

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
            std::cout << engine.evaluate_log_lik("") << std::endl;

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
