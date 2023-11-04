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

// FIXME: this depends on the model
std::vector<double> Cumul_Comm_Infec = {};
std::vector<double> Inst_Comm_Infec = {};

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
                cerr << "\nERROR: Unable to extract setup name." << endl;
                cerr << "Check line: " << line << "\n"
                     << endl;
                cerr << "Exit...\n"
                     << endl;
                exit(EXIT_FAILURE);
            }

            //=====CHOICE=====//
            if (!getline(iss, choice))
            {
                cerr << "\nERROR: Unable to extract setup choice." << endl;
                cerr << "Check line: " << line << "\n"
                     << endl;
                cerr << "Exit...\n"
                     << endl;
                exit(EXIT_FAILURE);
            }

            // add the new setup
            safe_load.emplace(string_cleaner(name), string_cleaner(choice));
        }

        std::string seed_str = safe_load.at("seed");
        if (seed_str == UNDEFINED)
        {
            loader.seed = time(NULL);
        }
        else
        {
            loader.seed = stoi(seed_str);
        }

        loader.num_iter = stoi(safe_load.at("num_iter"));
        loader.num_adapt_iter = stoi(safe_load.at("num_adapt_iter"));
        loader.thinning = stoi(safe_load.at("thinning"));
        loader.delim = load_delim(safe_load.at("delim"));
        loader.dict_delim = load_delim(safe_load.at("dict_delim"));
        loader.delta_t = stod(safe_load.at("delta_t"));

        //=====Check for values validity=====//
        if (loader.num_adapt_iter >= loader.num_iter)
        {
            std::cerr << "\nERROR: N. of adapative iteration >= N. of iterations." << endl;
            std::cerr << "Exit...\n"
                      << std::endl;
            std::exit(EXIT_FAILURE);
        }

        //=====Display=====//
        if (Verbose)
        {
            for (const auto &pair : safe_load)
            {
                std::cout << pair.first << " : " << pair.second << endl;
            }
        }

        // FLAG
        std::exit(1);
    }
    else
    {
        std::cerr << "\nERROR: Impossible to read: " << data_pathtofile << endl;
        std::cerr << "Exit...\n"
                  << std::endl;
        std::exit(EXIT_FAILURE);
    }

    return loader;
}

/* Read Extra File */
void build_extra(string data_pathtofile)
{
    double date, inst_val, cumul_val;

    ifstream input(data_pathtofile.c_str());
    if (input.is_open())
    {
        if (Verbose)
        {
            std::cout << "\n==============COMMUNITY===============" << std::endl;
            std::cout << "Data file: " << data_pathtofile << std::endl;
        }

        // Create a stringstream of the current line
        for (string line; getline(input, line);)
        {
            istringstream iss(line);

            // Store information in variables
            iss >> date >> inst_val >> cumul_val;

            Inst_Comm_Infec.push_back(inst_val);
            Cumul_Comm_Infec.push_back(cumul_val);
        }
    }
    else
    {
        cerr << "\nERROR: Impossible to read: " << data_pathtofile << endl;
        cerr << "Exit...\n"
             << endl;
        exit(EXIT_FAILURE);
    }
}

// /////////////////////////////////////////////////////////////////////////// //
// ////////////////////////////// MAIN FUNCTION ////////////////////////////// //
// /////////////////////////////////////////////////////////////////////////// //

int main(int argc, char *argv[])
{
    int carg_i = 1;             // Command line argument iterator
    double inver_temp = 1.0;    // Inverse temperature // UPDATE : temporary
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
        cerr << "The user did not supply any command.\n";
        cerr << "Check the list of commands through:\n"
             << argv[0] << " -h";
        cerr << "Exit...\n"
             << endl;
        exit(EXIT_FAILURE);
    }

    //=====Load commands and options=====//
    while (carg_i < argc)
    {
        // Command line argument asking for help
        if (string(argv[carg_i]) == string("-h") or string(argv[carg_i]) == string("--help"))
        {
            // Show the manual using man and installed manpages
            cerr << "Software Manual missing. Contact the developer." << endl;
            // End the program without error
            exit(EXIT_SUCCESS);
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

        else if (string(argv[carg_i]) == string("-c") || string(argv[carg_i]) == string("--community"))
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
            cerr << "Unrecognized option :" << argv[carg_i] << endl;
            cerr << "Exit...\n"
                 << endl;
            exit(EXIT_FAILURE);
        }

        ++carg_i; // Move to the next command line argument
    }

    // UPDATE : check for options validity

    if (infer && simulate)
    {
        cerr << "Incompatible options: choose one between -infer and -simulate." << endl;
        cerr << "Exit...\n"
             << endl;
        exit(EXIT_FAILURE);
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
    mcmc.pre_process_data();

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
