// main.cpp //

#include <string>
#include <vector>
#include <map>
#include <cmath>
#include <random>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <time.h>

#include "mcmc.h"

// UPDATE : understand namespace usage and clean up
using namespace std;

//=====Define the global variables=====//

bool Verbose = false;

// ////////////////////////////////////////////////////////////////////////// //
// ////////////////////////////// FILE LOADING ////////////////////////////// //
// ////////////////////////////////////////////////////////////////////////// //

/* Read Setup file */
Setup_MCMC_Struct load_setup(string data_pathtofile)
{
    char delim = ',';
    Setup_MCMC_Struct loader;
    ifstream input(data_pathtofile.c_str());
    map<string, string> safe_load = {};

    if (input.is_open())
    {
        if (Verbose)
        {
            std::cout << "\n================SETUP=================" << std::endl;
            std::cout << "Setup file: " << data_pathtofile << "\n"
                      << std::endl;
        }

        input.clear();            // Clear any errors
        input.seekg(0, ios::beg); // Seek to the beginning

        // Create a stringstream of the current line
        for (string line; getline(input, line);)
        {
            string name, choice;
            istringstream iss(line);

            //=====NAME=====//
            if (!getline(iss, name, delim))
            {
                std::string err_text = "Unable to extract setup name.\nCheck file line: " + line;
                exit_with_failure_error(err_text);
            }

            //=====CHOICE=====//
            if (!getline(iss, choice))
            {
                std::string err_text = "Unable to extract setup choice.\nCheck file line: " + line;
                exit_with_failure_error(err_text);
            }

            // add the new setup
            safe_load.emplace(string_cleaner(name), string_cleaner(choice));
        }

        std::string seed_str = safe_load.at("seed");
        loader.seed = (seed_str == UNDEFINED) ? time(NULL) : stoi(seed_str);
        loader.num_iter = stoi(safe_load.at("num_iter"));
        loader.num_adapt_iter = stoi(safe_load.at("num_adapt_iter"));
        loader.thinning = stoi(safe_load.at("thinning"));
        loader.delta_t = stod(safe_load.at("delta_t"));

        //=====Display=====//
        if (Verbose)
        {
            std::cout << "seed : " << loader.seed << endl;
            std::cout << "num_iter : " << loader.num_iter << endl;
            std::cout << "num_adapt_iter : " << loader.num_adapt_iter << endl;
            std::cout << "thinning : " << loader.thinning << endl;
            std::cout << "delta_t : " << loader.delta_t << endl;
        }
    }
    else
    {
        std::string err_text = "Impossible to read file: " + data_pathtofile + ".";
        exit_with_failure_error(err_text);
    }

    return loader;
}

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
    string setup_pathtofile = UNDEFINED;
    string parameters_pathtofile = UNDEFINED;
    string extra_pathtofile = UNDEFINED;
    string database_pathtofile = UNDEFINED;
    string output_pathtofile = UNDEFINED;
    string augmented_pathtofile = UNDEFINED;

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
        if (string(argv[carg_i]) == string("-h") or string(argv[carg_i]) == string("--help"))
        {
            // Show the manual using man and installed manpages // UPDATE : produce a manual
            std::string err_text = "We apologize for the inconvenience, but the Software Manual is still missing.\nPlease, contact the developer.";
            exit_with_failure_error(err_text);
        }

        else if (string(argv[carg_i]) == string("-infer"))
        {
            // Command line asking for infer
            infer = true;
        }

        else if (string(argv[carg_i]) == string("-simulate"))
        {
            // Command line asking for simulate
            simulate = true;
        }

        else if (string(argv[carg_i]) == string("-v") || string(argv[carg_i]) == string("--verbose"))
        {
            Verbose = true;
        }

        else if (string(argv[carg_i]) == string("--init_rand"))
        {
            rand_init_params = true;
        }

        else if (string(argv[carg_i]) == string("-s") || string(argv[carg_i]) == string("--setup"))
        {
            ++carg_i;
            setup_pathtofile = argv[carg_i];
        }

        else if (string(argv[carg_i]) == string("-p") || string(argv[carg_i]) == string("--params"))
        {
            ++carg_i;
            parameters_pathtofile = argv[carg_i];
        }

        else if (string(argv[carg_i]) == string("-e") || string(argv[carg_i]) == string("--extra"))
        {
            ++carg_i;
            extra_pathtofile = argv[carg_i];
        }

        else if (string(argv[carg_i]) == string("-d") || string(argv[carg_i]) == string("--database"))
        {
            ++carg_i;
            database_pathtofile = argv[carg_i];
        }

        else if (string(argv[carg_i]) == string("-o") || string(argv[carg_i]) == string("--output"))
        {
            ++carg_i;
            output_pathtofile = argv[carg_i];
        }

        else if (string(argv[carg_i]) == string("-i") || string(argv[carg_i]) == string("--inver_temp"))
        {
            ++carg_i;
            inver_temp = stod(argv[carg_i]);
        }

        else if (string(argv[carg_i]) == string("--save_augm"))
        {
            ++carg_i;
            augmented_pathtofile = argv[carg_i];
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
        if (augmented_pathtofile != UNDEFINED)
        {
            mcmc.print_augm_variables(augmented_pathtofile);
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
