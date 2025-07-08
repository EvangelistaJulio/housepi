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

//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//

/* Display a decorator for a section (if VERBOSE) */
void print_decorator(const std::string &title, bool verbose = true, const std::string &txt = HSP_NA_STRING)
{
    const ulong decor_left = 3;   // left size of decorator
    const ulong decor_width = 65; // total size of decorator
    const char decor = '-';       // decorator for prints

    if (verbose)
    {
        std::string sms = std::string(decor_left, decor);
        sms += " " + title + " ";
        sms += std::string(decor_width - title.length() - decor_left, decor) + "\n";
        if (!tools::is_na(txt))
            sms += txt + "\n";

        tools::print_stdout(sms);
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////  * *** * MAIN FUNCTION * *** *  /////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/* Main function for the standalone version of 'housepi'. */
int main(int argc, char *argv[])
{
    using tools::stoupper;
    using tools::stolower;

    int carg_i = 1;   // command line iterator
    std::string carg; // command line argument

    bool is_loglikelihood = false; // loglikelihood option
    bool is_infer = false;         // inference option
    bool is_simulate = false;      // simulation option
    ulong num_rungs = 0;           // n. rungs power posteriors

    // paths/to/files for import
    std::string augm_file = HSP_NA_STRING;
    std::string dataset_file = HSP_NA_STRING;
    std::string datalist_file = HSP_NA_STRING;
    std::string group_rates_file = HSP_NA_STRING;
    std::string infec_profile_file = HSP_NA_STRING;
    std::string multiplicity_file = HSP_NA_STRING;
    std::string output_file = HSP_NA_STRING;
    std::string param_file = HSP_NA_STRING;
    std::string setup_file = HSP_NA_STRING;

    bool verbose = false; // verbose option

    //=====Test if some commands were supplied=====//
    if (argc < 2)
    {
        std::ostringstream oss;
        oss << "No argument or command found. Check the usage through:\n"
            << argv[0] << " -h";
        tools::error(oss.str());
    }

    //=====Parse commands and options=====//
    while (carg_i < argc)
    {
        carg = std::string(argv[carg_i]); // update command line argument

        // Command line argument asking for help
        if (stoupper(carg) == "-H" || stolower(carg) == "--help")
        {
            // Show the manual using man and installed manpages
            // UPDATE : produce a manual
            std::ostringstream oss;
            oss << "No Software Manual is available, refer to the Rcpp tutorials.\n"
                << "Apologies for the inconvenience.";
            tools::error(oss.str());
        }

        else if (stolower(carg) == "--evidence")
            num_rungs = tools::stoul(argv[++carg_i], "number of rungs");

        else if (stolower(carg) == "--loglikelihood")
            is_loglikelihood = true;

        else if (stolower(carg) == "--infer")
            is_infer = true;

        else if (stolower(carg) == "--simulate")
            is_simulate = true;

        else if (stoupper(carg) == "-A" || stolower(carg) == "--augmented")
            augm_file = argv[++carg_i];

        else if (stoupper(carg) == "-D" || stolower(carg) == "--dataset")
            dataset_file = argv[++carg_i];

        else if (stoupper(carg) == "-L" || stolower(carg) == "--datalist")
            dataset_file = argv[++carg_i];

        else if (stoupper(carg) == "-G" || stolower(carg) == "--group-rates")
            group_rates_file = argv[++carg_i];

        else if (stoupper(carg) == "-I" || stolower(carg) == "--infec-profile")
            infec_profile_file = argv[++carg_i];

        else if (stoupper(carg) == "-M" || stolower(carg) == "--multiplicity")
            multiplicity_file = argv[++carg_i];

        else if (stoupper(carg) == "-O" || stolower(carg) == "--output")
            output_file = argv[++carg_i];

        else if (stoupper(carg) == "-P" || stolower(carg) == "--parameter")
            param_file = argv[++carg_i];

        else if (stoupper(carg) == "-S" || stolower(carg) == "--setup")
            setup_file = argv[++carg_i];

        else if (stoupper(carg) == "-V" || stolower(carg) == "--verbose")
            verbose = true;

        else
            tools::error("Unrecognized option `" + carg + "`.");

        // Move to the next command line argument
        carg_i++;
    }

    //==========================//
    //===== INITIALIZATION =====//
    //==========================//

    //----- Check for command -----//
    bool is_evidence = (num_rungs > 0); // power posterior option
    int aux = is_loglikelihood + is_infer + is_simulate + is_evidence;
    if (aux == 0)
        is_loglikelihood = true; // default
    else if (aux > 1)
        tools::error(
            "Incompatible commands: choose one between --loglikelihood, --infer, --simulate or --evidence.");

    //----- Initialize the custom Model object -----//
    std::unique_ptr<Model> user_model = std::make_unique<Model>();
    user_model->init_user_advanced_setup();

    // 1) Load Data
    user_model->read_dataset(dataset_file);
    user_model->read_datalist(datalist_file);

    // 2) Load Setup
    user_model->read_setup(setup_file);

    // 3) Load Parameters
    user_model->read_param(param_file);

    // 4) Load Setup
    user_model->read_infec_profile(infec_profile_file);

    // 5) Load Group & Contact Rates
    user_model->read_group_rates(group_rates_file);
    user_model->build_all_rates(); // temporary for contact rates

    // 6) Load Augmented Variables
    user_model->read_augm(augm_file);

    // 7) Load Multiplicities
    user_model->read_multiplicity(multiplicity_file);

    if (verbose)
    {
        //----- SETUP -----//
        print_decorator("Setup");
        user_model->show_setup();

        //----- PARAMETERS -----//
        print_decorator("Parameters");
        user_model->show_all_params();

        //----- DATABASE -----//
        print_decorator("Dataset");
        user_model->show_summary();

        // UPDATE: missing show_augm
    }

    //=====================//
    //===== EXECUTION =====//
    //=====================//

    if (is_loglikelihood)
    { //===== Evaluate Loglikelihood =====//
        print_decorator("Loglikelihood", verbose);
        tools::print_stdout(std::to_string(user_model->evaluate()) + "\n");
    }

    else if (is_infer)
    { //===== Inference =====//
        print_decorator("Inference", verbose);
        user_model->infer();

        if (verbose)
        {
            print_decorator("Inference::results");
            user_model->show_results();
        }

        if (output_file != HSP_NA_STRING)
        {
            print_decorator("Inference::infer state", verbose,
                            "# save: " + output_file);
            user_model->print_infer_state(output_file);
        }
    }

    else if (is_simulate)
    { //===== Simulation Study =====//
        if (output_file != HSP_NA_STRING)
        {
            print_decorator("Simulation", verbose,
                            "# save: " + output_file);
            user_model->print_simulated_dataset(output_file);
        }
        else
            tools::error("Missing output file for simulation.");
    }

    else if (is_evidence)
    { //===== Evidence =====//
        print_decorator("Evidence", verbose);
        user_model->print_evidence("adaptive", num_rungs);
    }

    return EXIT_SUCCESS;
}
