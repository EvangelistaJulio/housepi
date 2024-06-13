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
    int carg_i = 1;   // command line argument iterator
    std::string carg; // command line argument

    bool is_infer = true;    // inference option
    bool is_simulate = false; // simulation option

    // paths/to/files
    std::string params_file = UNDEFINED;
    std::string data_file = UNDEFINED;
    std::string augm_file = UNDEFINED;
    std::string setup_file = UNDEFINED;
    std::string extra_file = UNDEFINED;
    std::string output_pathtofile = UNDEFINED;
    std::string save_augm_pathtofile = UNDEFINED;

    //=====Test if some commands were supplied=====//
    if (argc < 2)
    {
        std::ostringstream oss;
        oss << "No argument or command found. Check the usage through:\n"
            << argv[0] << " -h";
        raise_error_message(oss);
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
            oss << "We apologize for the inconvenience, but the Software Manual is still missing.\n"
                << "Please, contact the developer.";
            raise_error_message(oss);
        }

        else if (carg == std::string("-infer"))
        {
            is_infer = true;
        }

        else if (carg == std::string("-simulate"))
        {
            is_simulate = true;
        }

        else if (carg == std::string("-s") || carg == std::string("--setup"))
        {
            setup_file = argv[++carg_i];
        }

        else if (carg == std::string("-p") || carg == std::string("--params"))
        {
            params_file = argv[++carg_i];
        }

        else if (carg == std::string("-e") || carg == std::string("--extra"))
        {
            extra_file = argv[++carg_i];
        }

        else if (carg == std::string("-d") || carg == std::string("--database"))
        {
            data_file = argv[++carg_i];
        }

        else if (carg == std::string("-a") || carg == std::string("--augmented"))
        {
            augm_file = argv[++carg_i];
        }

        else if (carg == std::string("-o") || carg == std::string("--output"))
        {
            output_pathtofile = argv[++carg_i];
        }

        else if (carg == std::string("--save_augm"))
        {
            save_augm_pathtofile = argv[++carg_i];
        }

        else
        {
            raise_error_message("Unrecognized option `" + carg + "`.");
        }

        carg_i++; // Move to the next command line argument
    }

    //=====INITIALIZATION=====//

    // Initialize Model object with user defined derived class
    Model user_model;

    user_model.build(
        is_infer, is_simulate,
        params_file, data_file,
        augm_file, setup_file,
        extra_file);

    //=====EXECUTION=====//

    if (is_infer)
    { // Parameter Inference
        user_model.inference(output_pathtofile);
        user_model.print_augm_variables(save_augm_pathtofile);
    }
    else
    { // Simulation Study
        if (output_pathtofile == UNDEFINED)
        {
            raise_error_message("Undefined output file.");
        }
        else
        {
            user_model.simulation();
            user_model.print_augm_variables(output_pathtofile);
        }
    }

    return EXIT_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------//
//------------------------------------------------------------------------------------------------------------//
//------------------------------------------------------------------------------------------------------------//
//------------------------------------------------------------------------------------------------------------//
//------------------------------------------------------------------------------------------------------------//
//------------------------------------------------------------------------------------------------------------//
