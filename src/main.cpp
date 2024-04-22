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
#include MODEL_HEADER

//=====Define the global variables=====//

bool Verbose = false;

// /////////////////////////////////////////////////////////////////////////// //
// ////////////////////////////// MAIN FUNCTION ////////////////////////////// //
// /////////////////////////////////////////////////////////////////////////// //

int main(int argc, char *argv[])
{
    int carg_i = 1;          // Command line argument iterator
    double inver_temp = 1.0; // Inverse temperature
    bool infer = false, simulate = false;
    bool rand_init_params = false;
    int num_threads =1;
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
        std::string err_text = "No argument or command found. Check the usage through:\n" + (std::string)(argv[0]) + " -h";
        raise_error_message(err_text);
        std::exit(EXIT_FAILURE);
    }

    //=====Load commands and options=====//
    while (carg_i < argc)
    {
        // Command line argument asking for help
        if (std::string(argv[carg_i]) == std::string("-h") or std::string(argv[carg_i]) == std::string("--help"))
        {
            // Show the manual using man and installed manpages // UPDATE : produce a manual
            std::string err_text = "We apologize for the inconvenience, but the Software Manual is still missing.\nPlease, contact the developer.";
            raise_error_message(err_text);
            std::exit(EXIT_FAILURE);
        }

        else if (std::string(argv[carg_i]) == std::string("-infer"))
        {
            // Command line asking for infer
            infer = true;
        }

        else if (std::string(argv[carg_i]) == std::string("-simulate"))
        {
            // Command line asking for simulate
            simulate = true;
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
            std::string aux = argv[carg_i];
            std::string err_text = "Unrecognized option `" + aux + "`.";
            raise_error_message(err_text);
            std::exit(EXIT_FAILURE);
        }

        ++carg_i; // Move to the next command line argument
    }

    //=====INITIALIZATION=====//

    display_decorator("Household Survey EPI Modeling");
    display_decorator("Camaglia F., Institut Pasteur (2024)", '-', 0);

    // Initialize Model object with user defined derived class
    Model user_model;

    // 0) set number of threads
    user_model.set_num_threads(num_threads);

    // 1) Load Setup (optional)
    user_model.build_setup(setup_pathtofile);

    // 2) Load Parameters and Check (required)
    if (parameters_pathtofile == UNDEFINED)
    {
        std::string err_text = "Required argument: parameter path/to/file not found.";
        raise_error_message(err_text);
        std::exit(EXIT_FAILURE);
    }
    else
    {
        user_model.build_parameters(parameters_pathtofile);
    }

    // 3) Load Data and Check (required)
    if (database_pathtofile == UNDEFINED)
    {
        std::string err_text = "Required argument: database path/to/file not found.";
        raise_error_message(err_text);
        std::exit(EXIT_FAILURE);
    }
    else
    {
        user_model.build_data(database_pathtofile);
    }

    // Load Augmented Variables
    user_model.build_augmented(augmented_pathtofile);

    // UPDATE: there could be plenty extra files to be store somewhere
    user_model.build_extra(extra_pathtofile);

    //=====EXECUTION=====//

    //=====check for options validity=====//

    // UPDATE : complete check for options validity

    if (infer && simulate)
    {
        std::string err_text = "Incompatible commands: choose one between -infer and -simulate.";
        raise_error_message(err_text);
        std::exit(EXIT_FAILURE);
    }
    else if (infer)
    { // Parameter Inference
        user_model.inference(output_pathtofile, rand_init_params, inver_temp);
        user_model.print_augm_variables(save_augm_pathtofile);
    }
    else if (simulate)
    { // Simulation Study
        user_model.simulation(output_pathtofile);
    }

    display_decorator("END");
    display_stdout("\n");

    return EXIT_SUCCESS;
}
