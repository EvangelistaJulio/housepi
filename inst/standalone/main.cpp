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
/////////////////////////////////  * *** * MAIN FUNCTION * *** *  /////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/* The Main Function. */
int main(int argc, char *argv[])
{
    using namespace screen;

    int carg_i = 1;   // command line iterator
    std::string carg; // command line argument

    bool verbose = false;     // verbose option
    bool is_infer = false;    // inference option
    bool is_simulate = false; // simulation option
    bool is_evaluate = false; // evaluate option
    bool is_marginal = false; // power posterior option
    ulong num_rungs = 0;

    // paths/to/files
    std::string params_file = UNDEFINED;
    std::string data_file = UNDEFINED;
    std::string augm_file = UNDEFINED;
    std::string setup_file = UNDEFINED;
    std::string extra_file = UNDEFINED;
    std::string output_file = UNDEFINED;
    std::string save_augm_file = UNDEFINED;

    //=====Test if some commands were supplied=====//
    if (argc < 2)
    {
        std::ostringstream oss;
        oss << "No argument or command found. Check the usage through:\n"
            << argv[0] << " -h";
        raise_error(oss.str());
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
            raise_error(oss.str());
        }

        else if ((carg == std::string("-v")) || carg == std::string("--verbose"))
        {
            verbose = true;
        }

        else if (carg == std::string("--evaluate"))
        {
            is_evaluate = true;
        }

        else if (carg == std::string("--infer"))
        {
            is_infer = true;
        }

        else if (carg == std::string("--simulate"))
        {
            is_simulate = true;
        }

        else if (carg == std::string("--marginal"))
        {
            is_marginal = true;
            try
            {
                num_rungs = std::stoul(argv[++carg_i]);
            }
            catch (...)
            {
                raise_error("Unrecognized number of rungs. Declare it as: --marginal <num_rungs>");
            }
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
            output_file = argv[++carg_i];
        }

        else if (carg == std::string("--save_augm"))
        {
            save_augm_file = argv[++carg_i];
        }

        else
        {
            raise_error("Unrecognized option `" + carg + "`.");
        }

        carg_i++; // Move to the next command line argument
    }

    //==========================//
    //===== INITIALIZATION =====//
    //==========================//

    //----- Check for command -----//
    int aux = is_evaluate + is_infer + is_simulate + is_marginal;
    if (aux > 1)
    {
        raise_error(
            "Incompatible commands: choose one between --evaluate, --infer, --simulate or --marginal.");
    }
    else if (aux == 0)
    {
        is_evaluate = true; // default
    }

    //----- Initialize the Model object -----//
    Model user_model;
    user_model.build(setup_file, params_file,
                     data_file, augm_file,
                     extra_file, verbose);

    //=====================//
    //===== EXECUTION =====//
    //=====================//

    if (is_infer)
    { //===== Inference =====//
        print_decorator("Inference", C_DECOR, "save: " + output_file, verbose);
        user_model.infer(output_file, true);

        if (verbose)
        {
            print_decorator("Inference::results", '-', UNDEFINED);
            user_model.show_results();
            print_decorator("Inference::diagnostics", '-', UNDEFINED);
            user_model.show_diagnostics();
        }

        if (save_augm_file != UNDEFINED)
        {
            print_decorator("Inference::output", '-',
                           "save: " + save_augm_file,
                           verbose);
            user_model.print_augm_variables(save_augm_file);
        }
    }

    else if (is_simulate)
    { //===== Simulation Study =====//
        if (output_file == UNDEFINED)
        {
            raise_error("Missing output file for simulation.");
        }
        else
        {
            print_decorator("Simulation", C_DECOR,
                           "save: " + output_file,
                           verbose);
            user_model.simulate();
            user_model.print_augm_variables(output_file);

            if (verbose)
            {
                print_decorator("Simulation::summary", '-');
                user_model.show_summary();
            }
        }
    }

    else if (is_evaluate)
    { //===== Evaluate Loglikelihood =====//
        print_decorator("Loglikelihood", C_DECOR, UNDEFINED, verbose);
        print_stdout(std::to_string(user_model.evaluate()) + "\n");
    }

    else if (is_marginal)
    { //===== Power Posteriors =====//
        print_decorator("Marginal Loglikelihood", C_DECOR, UNDEFINED, verbose);
        print_stdout(std::to_string(user_model.evaluate_marginal(num_rungs)) + "\n");
    }

    return EXIT_SUCCESS;
}
