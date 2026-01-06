/* housepi - Epidemiological Analysis of Household Outbreaks.
 *
 * Copyright © 2026 Institut Pasteur (Paris), CNRS and INSERM.
 * Distributed under the terms of the GNU General Public License (GPLv3).
 * See the LICENSE file for details.
 */

#include <string>
#include <fstream>
#include <iostream>
#include <iomanip>

#include "tools/na_utils.h"
#include "tools/messages.h"
#include "tools/string.h"
#include "tools/read.h"
#include "mcmc_impl.h"
#include "mcmc.h"

// Header to the user defined model provided by the Make command
#include "model.h"

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////  standalone print decorator  //////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/// Display a decorator for a section (if VERBOSE).
void print_decorator(const std::string &title, const std::string &txt = hsp::na::STRING)
{
    const unsigned decor_left = 3;   // left size of decorator
    const unsigned decor_width = 65; // total size of decorator
    const char decor = '-';          // decorator for prints

    std::cout << std::string(decor_left, decor)
              << " " << title << " "
              << std::string(decor_width - title.size() - decor_left, decor)
              << (hsp::na::is(txt) ? "" : "\n" + txt)
              << std::endl;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////  MCMC : Impl : read  ////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

#define MAKE_READ_METHOD(module, section_name)                            \
    void MCMC::read_##module(const std::string &f)                        \
    {                                                                     \
        try                                                               \
        {                                                                 \
            pImpl_->module.build(hsp::read::csv_table(f));                \
        }                                                                 \
        catch (const std::exception &e)                                   \
        {                                                                 \
            hsp::msg::error(std::string("While reading ") +               \
                            section_name + " definitions:\n" + e.what()); \
        }                                                                 \
    }

MAKE_READ_METHOD(setup, "setup")
MAKE_READ_METHOD(infec_profile, "infectivity profile")
MAKE_READ_METHOD(param, "parameters")
MAKE_READ_METHOD(data_wide, "data-wide")
MAKE_READ_METHOD(data_long, "data-long")
MAKE_READ_METHOD(group_rates, "group rates")
MAKE_READ_METHOD(augm, "augmented")

// UPDATE: this will be removed once automatized in C++ reflection
void MCMC::read_multiplicity(const std::string &f)
{
    try
    {
        pImpl_->data_wide.multiplicity_build(hsp::read::csv_table(f), pImpl_->augm.is_any());
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
void MCMC::print_infer_state(const std::string &output_filename) const
{
    if (hsp::na::is(output_filename))
        return;

    const Setup &setup = pImpl_->setup;
    const Impl::Parameter &param = pImpl_->param;
    const Impl::Augmented &augm = pImpl_->augm;

    std::ofstream ofs;
    ofs.open(output_filename.c_str());

    if (!ofs.is_open())
        throw std::runtime_error(
            "Unable to open file: " + output_filename);

    // Headers
    // header log-likelihood
    ofs << "iteration" << hsp::read::DELIMITER << "log_lik";

    // header parameters
    for (const std::string &p : param.list_sample)
        ofs << hsp::read::DELIMITER << p;

    // header augmented variable of id trace
    if (augm.is_any_draw())
        for (size_t k = 0; k < augm.smpl_individuals.size(); ++k)
            for (const std::string &a : augm.list_draws)
                ofs << hsp::read::DELIMITER << augm.smpl_individuals[k]->id() + "_" + a;

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
                << std::fixed << pImpl_->smpl_total_log_lik[i];

            // print parameter values
            for (const std::string &p : param.list_sample)
                ofs << hsp::read::DELIMITER << std::setprecision(5)
                    << std::fixed << param.smpl_values.find(p)->second[i];

            // print augmented variable of id trace
            if (!augm.smpl_values.empty())
                for (size_t k = 0; k < augm.smpl_individuals.size(); ++k)
                    for (const std::string &a : augm.list_draws)
                        ofs << hsp::read::DELIMITER << std::setprecision(5)
                            << std::fixed << augm.smpl_values[k].find(a)->second[i];

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
                                     const str_vec &data_headers,
                                     const str_vec &augm_headers)
{
    std::ofstream ofs;
    ofs.open(f.c_str());

    if (!ofs.is_open())
        throw std::runtime_error("Unable to open file: '" + f + "'");

    // print default data headers
    ofs << "id" << hsp::read::DELIMITER << "house_id" << hsp::read::DELIMITER << "status";

    // print other data headers
    for (const std::string &d : data_headers)
        ofs << hsp::read::DELIMITER << d;

    // print augmented variables headers
    for (const std::string &a : augm_headers)
        ofs << hsp::read::DELIMITER << a;

    ofs << std::endl;

    return ofs;
}

/// Print all augmented variables of each individual.
void update_ofstream_datawide(std::ofstream &ofs,
                              Household &house,
                              const str_vec &data_headers,
                              const str_vec &augm_headers)
{
    if (!ofs.is_open())
        throw std::runtime_error("Output file closed unexpectedly.");

    for (const Individual &member : house)
    {
        ofs << member.id()
            << hsp::read::DELIMITER << house.house_id()
            << hsp::read::DELIMITER << member.status();
        for (const std::string &d : data_headers)
            ofs << hsp::read::DELIMITER << member.data(d);

        for (const std::string &a : augm_headers)
            ofs << hsp::read::DELIMITER << member.aug(a);

        ofs << std::endl;
    }
}

/// Print the simulated data to the provided file.
void MCMC::print_simulated_data_wide(const std::string &pathtofile)
{
    const DataWide &data_wide = pImpl_->data_wide;
    const Impl::Augmented &augm = pImpl_->augm;

    try
    {
        // initialization
        std::ofstream ofs = open_ofstream_datawide(pathtofile,
                                                   data_wide.headers,
                                                   augm.list_export);
        pImpl_->helper_init_simulation();

        for (Household house : data_wide) // NOTE: deep copy of household
            try
            {
                // inplace simulation
                pImpl_->simulate_single_household_in_place(house);
                // print the simulated data
                update_ofstream_datawide(ofs, house,
                                         data_wide.headers,
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
double MCMC::log_evidence(const std::string &scheme, unsigned long n_rungs)
{
    try
    {
        double output = 0;
        const std::vector<dbl_vec> pp = pImpl_->power_posteriors(scheme, n_rungs);
        const dbl_vec &inv_tem = pp[0];
        const dbl_vec &avg_vec = pp[1];
        const dbl_vec &var_vec = pp[2];

        for (unsigned long k = 0; k < pImpl_->setup.n_rungs; ++k)
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
int main(int argc, char *argv[])
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
    std::string augm_file = hsp::na::STRING;
    std::string data_wide_file = hsp::na::STRING;
    std::string data_long_file = hsp::na::STRING;
    std::string group_rates_file = hsp::na::STRING;
    std::string infec_profile_file = hsp::na::STRING;
    std::string multiplicity_file = hsp::na::STRING;
    std::string output_file = hsp::na::STRING;
    std::string param_file = hsp::na::STRING;
    std::string setup_file = hsp::na::STRING;

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

    // initialize the user model object
    std::unique_ptr<Model> user_model = std::make_unique<Model>();
    user_model->init_advanced_setup();

    user_model->read_data_wide(data_wide_file);         // 1) Load Data Wide
    user_model->read_data_long(data_long_file);         // ... and Data Long
    user_model->read_setup(setup_file);                 // 2) Load Setup
    user_model->read_param(param_file);                 // 3) Load Parameters
    user_model->read_infec_profile(infec_profile_file); // 4) Load Infectivity Profile
    user_model->read_group_rates(group_rates_file);     // 5) Load Group Rates
    user_model->build_all_rates();                      // ... and Contact Rates (temporary)
    user_model->read_augm(augm_file);                   // 6) Load Augmented Variables
    user_model->read_multiplicity(multiplicity_file);   // 7) Load Multiplicities

    //////////////////////////////
    ////////  EXECUTION  /////////
    //////////////////////////////

    if (is_loglikelihood) // evaluate loglikelihood
    {
        std::cout << user_model->evaluate_log_lik() << std::endl;
    }
    else if (is_infer) // perform inference
    {
        if (hsp::na::is(output_file))
            hsp::msg::error("Missing output file for inference.");
        user_model->infer();
        user_model->print_infer_state(output_file);
    }
    else if (is_simulate) // simulation an outbreak
    {
        if (hsp::na::is(output_file))
            hsp::msg::error("Missing output file for simulation.");
        user_model->print_simulated_data_wide(output_file);
    }
    else if (is_logevidence) // estimate log_evidence
    {
        std::cout << user_model->log_evidence("adaptive", num_rungs) << std::endl;
    }

    return EXIT_SUCCESS;
}
