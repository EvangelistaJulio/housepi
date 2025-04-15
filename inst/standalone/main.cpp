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

///////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////  * *** * MAIN FUNCTION * *** *  /////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/* Main function for the standalone version of 'housepi'. */
int main(int argc, char *argv[])
{
    int carg_i = 1;   // command line iterator
    std::string carg; // command line argument

    bool is_loglikelihood = false; // loglikelihood option
    bool is_infer = false;         // inference option
    bool is_simulate = false;      // simulation option
    ulong num_rungs = 0;           // n. rungs power posteriors

    // paths/to/files for import
    std::string augm_file = NULL_STRING;
    std::string dataset_file = NULL_STRING;
    std::string infec_profile_file = NULL_STRING;
    std::string multiplicity_file = NULL_STRING;
    std::string output_file = NULL_STRING;
    std::string param_file = NULL_STRING;
    std::string setup_file = NULL_STRING;
    
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
        if (carg == std::string("-h") or carg == std::string("--help"))
        {
            // Show the manual using man and installed manpages
            // UPDATE : produce a manual
            std::ostringstream oss;
            oss << "No Software Manual is available, refer to the Rcpp tutorials.\n"
                << "Apologies for the inconvenience.";
            tools::error(oss.str());
        }

        else if (carg == std::string("--evidence"))
            num_rungs = tools::stoul(argv[++carg_i], "number of rungs");

        else if (carg == std::string("--loglikelihood"))
            is_loglikelihood = true;

        else if (carg == std::string("--infer"))
            is_infer = true;

        else if (carg == std::string("--simulate"))
            is_simulate = true;

        else if (carg == std::string("-a") || carg == std::string("--augmented"))
            augm_file = argv[++carg_i];

        else if (carg == std::string("-d") || carg == std::string("--dataset"))
            dataset_file = argv[++carg_i];

        else if (carg == std::string("-i") || carg == std::string("--infec-profile"))
            infec_profile_file = argv[++carg_i];

        else if (carg == std::string("-m") || carg == std::string("--multiplicity"))
            multiplicity_file = argv[++carg_i];

        else if (carg == std::string("-o") || carg == std::string("--output"))
            output_file = argv[++carg_i];

        else if (carg == std::string("-p") || carg == std::string("--parameter"))
            param_file = argv[++carg_i];

        else if (carg == std::string("-s") || carg == std::string("--setup"))
            setup_file = argv[++carg_i];

        else if ((carg == std::string("-v")) || carg == std::string("--verbose"))
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

    // 1) Load Data
    user_model->read_dataset(dataset_file);

    // 2) Load Setup
    user_model->read_setup(setup_file);

    // 3) Load Parameters
    user_model->read_param(param_file);

    // 4) Load Setup
    user_model->read_infec_profile(infec_profile_file);

    // 5) Load Group & Contact Rates
    user_model->build_all_rates();

    // 6) Load Augmented Variables
    user_model->read_augm(augm_file);

    // 7) Load Multiplicities
    user_model->read_multiplicity(multiplicity_file);

    if (verbose)
        user_model->show_full_build();

    //=====================//
    //===== EXECUTION =====//
    //=====================//

    if (is_loglikelihood)
    { //===== Evaluate Loglikelihood =====//
        tools::print_decorator("Loglikelihood", NULL_STRING, C_DECOR, verbose);
        tools::print_stdout(std::to_string(user_model->evaluate()) + "\n");
    }

    else if (is_infer)
    { //===== Inference =====//
        tools::print_decorator("Inference", NULL_STRING, C_DECOR, verbose);
        user_model->infer();

        if (verbose)
        {
            tools::print_decorator("Inference::results", NULL_STRING);
            user_model->show_results();
        }

        if (output_file != NULL_STRING)
        {
            tools::print_decorator("Inference::infer state",
                                   "# save: " + output_file,
                                   C_DECOR, verbose);
            user_model->print_infer_state(output_file);
        }
    }

    else if (is_simulate)
    { //===== Simulation Study =====//
        if (output_file != NULL_STRING)
        {
            tools::print_decorator("Simulation",
                                   "# save: " + output_file,
                                   C_DECOR, verbose);
            user_model->print_simulated_dataset(output_file);
        }
        else
            tools::error("Missing output file for simulation.");
    }

    else if (is_evidence)
    { //===== Evidence =====//
        tools::print_decorator("Evidence", NULL_STRING, C_DECOR, verbose);
        user_model->print_evidence("adaptive", num_rungs);
    }

    return EXIT_SUCCESS;
}
