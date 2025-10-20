/*
 *      __
 *     / /_  ____  __  __ ____ __  ____  __
 *    / __ \/ __ \/ / / / ___/ _ \/ __ \/ /
 *   / / / / /_/ / /_/ (__  )  __/ /_/ / /
 *  /_/ /_/\____/\____/____/\___/ ____/_/
 *                             /_/
 *
 * Copyright (C) 2025 Francesco Camaglia, Institut Pasteur
 *
 * main.cpp
 *
 */

#include <string>
#include <iostream>

#include "tools/na_utils.h"
#include "tools/error.h"
#include "tools/string.h"
#include "mcmc.h"

// Header to the user defined model provided by the Make command
#include "model.h"

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////  standalone print decorator  //////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/// Display a decorator for a section (if VERBOSE).
void print_decorator(const std::string &title, bool verbose = true, const std::string &txt = hsp::na::STRING)
{
    const unsigned decor_left = 3;   // left size of decorator
    const unsigned decor_width = 65; // total size of decorator
    const char decor = '-';          // decorator for prints

    if (verbose)
        std::cout << std::string(decor_left, decor)
                  << " " << title << " "
                  << std::string(decor_width - title.size() - decor_left, decor)
                  << (hsp::na::is(txt) ? "" : "\n" + txt)
                  << std::endl;
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

    bool verbose = false;          // verbose option
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
        if (toup(carg) == "-H" || tolow(carg) == "--help")
            // Show the manual using man and installed manpages
            // UPDATE : produce a manual
            hsp::msg::error(
                "No Software Manual is available, refer to the Rcpp tutorials.\nApologies for the inconvenience.");

        else if (tolow(carg) == "--logevidence")
            num_rungs = hsp::str::toul(argv[++carg_i], "number of rungs");

        else if (tolow(carg) == "--loglikelihood")
            is_loglikelihood = true;

        else if (tolow(carg) == "--infer")
            is_infer = true;

        else if (tolow(carg) == "--simulate")
            is_simulate = true;

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

        else if (toup(carg) == "-V" || tolow(carg) == "--verbose")
            verbose = true;

        else
            hsp::msg::error("Unrecognized option `" + carg + "`.");

        carg_i++; // move to the next command line argument
    }

    // check for command
    bool is_logevidence = (num_rungs > 0); // power posterior option
    int aux = is_loglikelihood + is_infer + is_simulate + is_logevidence;
    if (aux == 0)
        is_loglikelihood = true; // default
    else if (aux > 1)
        hsp::msg::error(
            "Incompatible commands: choose one between --loglikelihood, --infer, --simulate or --logevidence.");

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
        print_decorator("Loglikelihood", verbose);
        std::cout << user_model->evaluate() << std::endl;
    }
    else if (is_infer) // perform inference
    {
        if (!hsp::na::is(output_file))
        {
            print_decorator("Inference", verbose);
            user_model->infer();
            print_decorator("Inference::infer state", verbose, "# save: " + output_file);
            user_model->print_infer_state(output_file);
        }
        else
            hsp::msg::error("Missing output file for inference.");
    }
    else if (is_simulate) // simulation an outbreak
    {
        if (!hsp::na::is(output_file))
        {
            print_decorator("Simulation", verbose, "# save: " + output_file);
            user_model->print_simulated_data_wide(output_file);
        }
        else
            hsp::msg::error("Missing output file for simulation.");
    }
    else if (is_logevidence) // estimate logevidence
    {
        print_decorator("Evidence", verbose);
        std::cout << user_model->logevidence("adaptive", num_rungs) << std::endl;
    }

    return EXIT_SUCCESS;
}
