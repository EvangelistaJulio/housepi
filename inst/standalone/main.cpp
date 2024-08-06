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
#include "model.h"

//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////  * ***** * MAIN FUNCTION * ***** *  ///////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/* The Main Function. */
int main(int argc, char *argv[])
{
    int carg_i = 1;   // command line argument iterator
    std::string carg; // command line argument

    bool verbose = false;     // print relavent information to stdout
    bool is_infer = false;    // inference option
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
        screen::raise_error(oss.str());
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
            screen::raise_error(oss.str());
        }

        else if ((carg == std::string("-v")) || carg == std::string("--verbose"))
        {
            verbose = true;
        }

        else if (carg == std::string("--infer"))
        {
            is_infer = true;
        }

        else if (carg == std::string("--simulate"))
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
            screen::raise_error("Unrecognized option `" + carg + "`.");
        }

        carg_i++; // Move to the next command line argument
    }

    //=====INITIALIZATION=====//

    // Initialize Model object with user defined derived class
    Model user_model;

    if (is_infer && is_simulate)
    {
        screen::raise_error("Incompatible commands: choose one between -infer and -simulate.");
    }

    user_model.build(setup_file, params_file, 
                     data_file, augm_file, 
                     extra_file, verbose);

    //=====================//
    //===== EXECUTION =====//
    //=====================//

    if (is_infer)
    { //===== Inference =====//
        screen::display_decorator("Inference", screen::C_DECOR, UNDEFINED, verbose);
        user_model.infer(output_pathtofile);

        if (verbose)
        {
            screen::display_decorator("Inference::results", '-', UNDEFINED);
            user_model.show_param_diagnostics();
            screen::display_decorator("Inference::diagnostics", '-', UNDEFINED);
            user_model.show_results();
        }

        if (save_augm_pathtofile != UNDEFINED)
        {
            screen::display_decorator("Inference::output", '-', "save: " + save_augm_pathtofile, verbose);
            user_model.print_augm_variables(save_augm_pathtofile);
        }
    }

    else if (is_simulate)
    { //===== Simulation Study =====//
        if (output_pathtofile == UNDEFINED)
        {
            screen::raise_error("Missing output file for simulation.");
        }
        else
        {
            screen::display_decorator("Simulation", screen::C_DECOR, "save: " + output_pathtofile, verbose);
            user_model.simulate();
            user_model.print_augm_variables(output_pathtofile);
            if (verbose)
            {
                screen::display_decorator("Simulation::summary", '-');
               user_model.show_summary();
            }
        }
    }

    else
    { //===== Compute Loglikelihood =====//
        screen::display_decorator("Loglikelihood", screen::C_DECOR, UNDEFINED, verbose);
        screen::print_stdout(std::to_string(user_model.evaluate()) + "\n");
    }

    return EXIT_SUCCESS;
}
