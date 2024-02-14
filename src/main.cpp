// main.cpp //

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

//=====Define the global variables=====//

bool Verbose = false;

// /////////////////////////////////////////////////////////////////////////// //
// ////////////////////////////// MAIN FUNCTION ////////////////////////////// //
// /////////////////////////////////////////////////////////////////////////// //

int main(int argc, char *argv[])
{
    int carg_i = 1;               // Command line argument iterator
    double inver_temp = INFINITY; // Inverse temperature // UPDATE : temporary
    bool infer = false;
    bool simulate = false;
    bool rand_init_params = false;
    Setup_MCMC_Struct setup = {}; // default setup
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
        std::string aux = argv[0];
        std::string err_text = "The user did not supply any command.\nCheck the list of commands through:\n" + aux + " -h";
        exit_with_failure_error(err_text);
    }

    //=====Load commands and options=====//
    while (carg_i < argc)
    {
        // Command line argument asking for help
        if (std::string(argv[carg_i]) == std::string("-h") or std::string(argv[carg_i]) == std::string("--help"))
        {
            // Show the manual using man and installed manpages // UPDATE : produce a manual
            std::string err_text = "We apologize for the inconvenience, but the Software Manual is still missing.\nPlease, contact the developer.";
            exit_with_failure_error(err_text);
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

        else
        {
            std::string aux = argv[carg_i];
            std::string err_text = "Unrecognized option `" + aux + "`.";
            exit_with_failure_error(err_text);
        }

        ++carg_i; // Move to the next command line argument
    }

    // UPDATE : check for options validity

    if (infer && simulate)
    {
        std::string err_text = "Incompatible options: choose one between -infer and -simulate.";
        exit_with_failure_error(err_text);
    }

    //=====INITIALIZATION=====//

    // Load Setup
    setup = load_setup(setup_pathtofile);

    // Initialize MCMC object
    MCMC mcmc(setup);

    // Load Parameters
    mcmc.build_parameters(parameters_pathtofile);

    // Load Augmented Variables
    mcmc.build_augmented(augmented_pathtofile);

    // Load Data
    mcmc.build_data(database_pathtofile);

    // UPDATE: there could be plenty extra files to be store somewhere
    if (extra_pathtofile != UNDEFINED)
    {
        build_extra(extra_pathtofile);
    }

    //=====EXECUTION=====//

    if (infer) // Parameter Inference
    {
        mcmc.inference(output_pathtofile, rand_init_params, inver_temp);
        if (save_augm_pathtofile != UNDEFINED)
        {
            mcmc.print_augm_variables(save_augm_pathtofile);
        }
    }
    else if (simulate) // Simulation Study
    {
        mcmc.simulation(output_pathtofile);
    }

    if (Verbose)
    {
        std::cout << "\n=================END==================" << std::endl;
    }

    return EXIT_SUCCESS;
}
