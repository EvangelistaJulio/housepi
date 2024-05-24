// main.cpp //
// Copyright (C) 2024 Francesco Camaglia, Institut Pasteur //

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
#include MODEL_HEADER

bool Verbose = false; // assign externally defined global variable

//------------------------------------------------------------------------------------------------------------//
//------------------------------------------------------------------------------------------------------------//
//------------------------------------------------------------------------------------------------------------//
//------------------------------------------------------------------------------------------------------------//
//------------------------------------------------------------------------------------------------------------//
//------------------------------------------------------------------------------------------------------------//

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////  * ***** * MAIN FUNCTION * ***** *  ///////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/* The Main Function. */
int main(int argc, char *argv[])
{
    int carg_i = 1;          // Command line argument iterator
    double inver_temp = 1.0; // Inverse temperature
    bool rand_init_params = false;
    int num_threads = 1;
    bool is_infer = false;
    bool is_simulate = false;

    // paths/to/files
    std::string setup_pathtofile = UNDEFINED;
    std::string parameters_pathtofile = UNDEFINED;
    std::string extra_pathtofile = UNDEFINED;
    std::string database_pathtofile = UNDEFINED;
    std::string output_pathtofile = UNDEFINED;
    std::string augmented_pathtofile = UNDEFINED;
    std::string save_augm_pathtofile = UNDEFINED;

    // ===== Test if some commands were supplied ===== //
    if (argc < 2)
    {
        std::ostringstream oss;
        oss << "No argument or command found. Check the usage through:\n"
            << argv[0] << " -h";
        display_error_message(oss);
        std::exit(EXIT_FAILURE);
    }

    //=====Load commands and options=====//
    while (carg_i < argc)
    {
        // Command line argument asking for help
        if (std::string(argv[carg_i]) == std::string("-h") or std::string(argv[carg_i]) == std::string("--help"))
        {
            // Show the manual using man and installed manpages // UPDATE : produce a manual
            std::ostringstream oss;
            oss << "We apologize for the inconvenience, but the Software Manual is still missing.\n"
                << "Please, contact the developer.";
            display_error_message(oss);
            std::exit(EXIT_FAILURE);
        }

        else if (std::string(argv[carg_i]) == std::string("-infer"))
        {
            // Command line asking for infer
            is_infer = true;
        }

        else if (std::string(argv[carg_i]) == std::string("-simulate"))
        {
            // Command line asking for simulate
            is_simulate = true;
        }

        else if (std::string(argv[carg_i]) == std::string("-v") || std::string(argv[carg_i]) == std::string("--verbose"))
        {
            Verbose = true;
        }

        else if (std::string(argv[carg_i]) == std::string("--init_rand"))
        {
            rand_init_params = true;
        }

        else if (std::string(argv[carg_i]) == std::string("-s") || std::string(argv[carg_i]) == std::string("--setup"))
        {
            ++carg_i;
            setup_pathtofile = argv[carg_i];
        }

        else if (std::string(argv[carg_i]) == std::string("-p") || std::string(argv[carg_i]) == std::string("--params"))
        {
            ++carg_i;
            parameters_pathtofile = argv[carg_i];
        }

        else if (std::string(argv[carg_i]) == std::string("-e") || std::string(argv[carg_i]) == std::string("--extra"))
        {
            ++carg_i;
            extra_pathtofile = argv[carg_i];
        }

        else if (std::string(argv[carg_i]) == std::string("-d") || std::string(argv[carg_i]) == std::string("--database"))
        {
            ++carg_i;
            database_pathtofile = argv[carg_i];
        }

        else if (std::string(argv[carg_i]) == std::string("-a") || std::string(argv[carg_i]) == std::string("--augmented"))
        {
            ++carg_i;
            augmented_pathtofile = argv[carg_i];
        }

        else if (std::string(argv[carg_i]) == std::string("-o") || std::string(argv[carg_i]) == std::string("--output"))
        {
            ++carg_i;
            output_pathtofile = argv[carg_i];
        }

        else if (std::string(argv[carg_i]) == std::string("-i") || std::string(argv[carg_i]) == std::string("--inver_temp"))
        {
            ++carg_i;
            inver_temp = std::stod(argv[carg_i]);
        }

        else if (std::string(argv[carg_i]) == std::string("--save_augm"))
        {
            ++carg_i;
            save_augm_pathtofile = argv[carg_i];
        }

        else if (std::string(argv[carg_i]) == std::string("--num_threads"))
        {
            ++carg_i;
            num_threads = std::stoi(argv[carg_i]);
        }

        else
        {
            display_error_message("Unrecognized option `" + std::string(argv[carg_i]) + "`.");
            std::exit(EXIT_FAILURE);
        }

        ++carg_i; // Move to the next command line argument
    }
    
    //=====OPTION VALIDITY=====//

    // UPDATE : complete check for options validity

    if (is_infer && is_simulate)
    {
        display_error_message("Incompatible commands: choose one between -infer and -simulate.");
        std::exit(EXIT_FAILURE);
    }

    //=====INITIALIZATION=====//

    display_decorator("Household Survey EPI Modeling");
    display_decorator("Camaglia F., Institut Pasteur (2024)", 0, '-');

    // Initialize Model object with user defined derived class
    Model user_model;

    // 0) Load Parameters and Check (required)
    if (parameters_pathtofile == UNDEFINED)
    {
        display_error_message("Required argument: parameter path/to/file not found.");
        std::exit(EXIT_FAILURE);
    }
    else
    {
        user_model.build_parameters(parameters_pathtofile);
    }

    // 1) Load Data and Check (required)
    if (database_pathtofile == UNDEFINED)
    {
        display_error_message("Required argument: database path/to/file not found.");
        std::exit(EXIT_FAILURE);
    }
    else
    {
        user_model.build_data(database_pathtofile);
    }

    // 2) Load Augmented Variables (optional)
    user_model.build_augmented(augmented_pathtofile);

    // 3) Set number of threads (optional)
    user_model.set_num_threads(num_threads);
    
    // 4) Load Setup (optional)
    user_model.build_setup(setup_pathtofile, is_infer);

    // 5) Load Extra File (optional, user defined)
    // UPDATE: there could be plenty extra files to be store somewhere
    user_model.build_extra(extra_pathtofile);

    //=====EXECUTION=====//

    if (is_infer)
    { // Parameter Inference
        user_model.inference(output_pathtofile, rand_init_params, inver_temp);
        user_model.print_augm_variables(save_augm_pathtofile);
    }
    else
    { // Simulation Study
        user_model.simulation();
        user_model.print_augm_variables(output_pathtofile);
    }

    display_decorator("END");
    display_stdout("\n");

    return EXIT_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------//
//------------------------------------------------------------------------------------------------------------//
//------------------------------------------------------------------------------------------------------------//
//------------------------------------------------------------------------------------------------------------//
//------------------------------------------------------------------------------------------------------------//
//------------------------------------------------------------------------------------------------------------//
