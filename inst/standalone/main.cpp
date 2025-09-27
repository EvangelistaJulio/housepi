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
#include <map>
#include <cmath>
#include <random>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <time.h>

#include "mcmc.h"

// Header to the user defined model provided by the Make command
#include "model.h"

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////  standalone print decorator  //////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/// Display a decorator for a section (if VERBOSE).
void print_decorator(const std::string &title, bool verbose = true, const std::string &txt = HSP_NA_STRING)
{
    const unsigned decor_left = 3;   // left size of decorator
    const unsigned decor_width = 65; // total size of decorator
    const char decor = '-';          // decorator for prints

    if (verbose)
    {
        std::cout << std::string(decor_left, decor);
        std::cout << " " << title << " ";
        std::cout << std::string(decor_width - title.length() - decor_left, decor);
        if (!tools::is_na(txt))
            std::cout << "\n"
                      << txt;
        std::cout << std::endl;
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
    using tools::str::tolow;
    using tools::str::toup;

    int carg_i = 1;   // command line iterator
    std::string carg; // command line argument

    bool is_loglikelihood = false; // loglikelihood option
    bool is_infer = false;         // inference option
    bool is_simulate = false;      // simulation option
    unsigned long num_rungs = 0;   // n. rungs power posteriors

    // paths/to/files for import
    std::string augm_file = HSP_NA_STRING;
    std::string data_wide_file = HSP_NA_STRING;
    std::string data_long_file = HSP_NA_STRING;
    std::string group_rates_file = HSP_NA_STRING;
    std::string infec_profile_file = HSP_NA_STRING;
    std::string multiplicity_file = HSP_NA_STRING;
    std::string output_file = HSP_NA_STRING;
    std::string param_file = HSP_NA_STRING;
    std::string setup_file = HSP_NA_STRING;

    bool verbose = false; // verbose option

    // check any command was supplied
    if (argc < 2)
    {
        std::ostringstream oss;
        oss << "No argument or command found. Check the usage through:\n"
            << argv[0] << " -h";
        tools::error(oss.str());
    }

    // parse commands and options
    while (carg_i < argc)
    {
        carg = std::string(argv[carg_i]); // update command line argument

        // Command line argument asking for help
        if (toup(carg) == "-H" || tolow(carg) == "--help")
        {
            // Show the manual using man and installed manpages // UPDATE : produce a manual
            std::ostringstream oss;
            oss << "No Software Manual is available, refer to the Rcpp tutorials.\n"
                << "Apologies for the inconvenience.";
            tools::error(oss.str());
        }

        else if (tolow(carg) == "--evidence")
            num_rungs = tools::str::toul(argv[++carg_i], "number of rungs");

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
            tools::error("Unrecognized option `" + carg + "`.");

        carg_i++; // move to the next command line argument
    }

    // check for command
    bool is_evidence = (num_rungs > 0); // power posterior option
    int aux = is_loglikelihood + is_infer + is_simulate + is_evidence;
    if (aux == 0)
        is_loglikelihood = true; // default
    else if (aux > 1)
        tools::error(
            "Incompatible commands: choose one between --loglikelihood, --infer, --simulate or --evidence.");

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

    if (verbose)
    {
        print_decorator("Setup");
        user_model->show_setup();
        print_decorator("Parameters");
        user_model->show_all_params();
        print_decorator("Dataset");
        user_model->show_summary();
        // UPDATE: missing show_augm
    }

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
        if (!tools::is_na(output_file))
        {
            print_decorator("Inference", verbose);
            user_model->infer();
            print_decorator("Inference::infer state", verbose, "# save: " + output_file);
            user_model->print_infer_state(output_file);

            if (verbose)
            {
                print_decorator("Inference::results");
                user_model->show_results();
            }
        }
        else
            tools::error("Missing output file for inference.");
    }
    else if (is_simulate) // simulation an outbreak
    {
        if (!tools::is_na(output_file))
        {
            print_decorator("Simulation", verbose, "# save: " + output_file);
            user_model->print_simulated_data_wide(output_file);
        }
        else
            tools::error("Missing output file for simulation.");
    }
    else if (is_evidence) // estimate evidence
    {
        print_decorator("Evidence", verbose);
        user_model->print_evidence("adaptive", num_rungs);
    }

    return EXIT_SUCCESS;
}
