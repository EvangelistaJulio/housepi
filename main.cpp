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
            cout << "\n=============SETUP==============" << endl;
            cout << "Setup file: " << data_pathtofile << endl;
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
            string_cleaner(&name);
            string_cleaner(&choice);

            // add the new setup
            safe_load.emplace(name, choice);
        }

        loader.seed = stoi(safe_load.at("seed"));
        loader.num_iter = stoi(safe_load.at("num_iter"));
        loader.num_adapt_iter = stoi(safe_load.at("num_adapt_iter"));
        loader.thinning = stoi(safe_load.at("thinning"));
        loader.delim = load_delim(safe_load.at("delim"));
        loader.dict_delim = load_delim(safe_load.at("dict_delim"));
        loader.delta_t = stod(safe_load.at("delta_t"));
    }
    else
    {
        cerr << "\nERROR: Impossible to read: " << data_pathtofile << endl;
        cerr << "Exit...\n"
             << endl;
        exit(EXIT_FAILURE);
    }

    //=====Check for values validity=====//
    if (loader.num_adapt_iter >= loader.num_iter)
    {
        cerr << "\nERROR: N. of adapative iteration >= N. of iterations." << endl;
        cerr << "Exit...\n"
             << endl;
        exit(EXIT_FAILURE);
    }

    //=====Display=====//
    if (Verbose)
    {
        for (const auto &pair : safe_load)
        {
            cout << pair.first << " : " << pair.second << endl;
        }
    }

    return loader;
}

/* Read Community file */
void build_community(string data_pathtofile)
{
    double date, inst_val, cumul_val;

    ifstream input(data_pathtofile.c_str());
    if (input.is_open())
    {
        if (Verbose)
        {
            cout << "\n=============COMMUNITY==============" << endl;
            cout << "Data file: " << data_pathtofile << endl;
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
    int carg_i = 1;    // Command line argument iterator
    bool infer = false;
    bool simulate = false;
    bool rand_init_params = false;
    Setup_MCMC_Struct setup = {}; // default setup
    string parameters_pathtofile = UNDEFINED;
    string community_pathtofile = UNDEFINED;
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
            cerr << "Manual missing. Contact the developer." << endl;
            // End the program without error
            exit(EXIT_SUCCESS);
        }

        if (string(argv[carg_i]) == string("-v"))
        {
            Verbose = true;
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

        else if (string(argv[carg_i]) == string("-s") || string(argv[carg_i]) == string("--setup"))
        {
            ++carg_i;
            setup = load_setup(argv[carg_i]);
        }

        else if (string(argv[carg_i]) == string("-p") || string(argv[carg_i]) == string("--params"))
        {
            ++carg_i;
            parameters_pathtofile = argv[carg_i];
        }

        else if (string(argv[carg_i]) == string("-c") || string(argv[carg_i]) == string("--community"))
        {
            ++carg_i;            
            community_pathtofile = argv[carg_i];
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

    // ===== Initialize MCMC object ===== //
    MCMC mcmc(setup);

    //=====Load Data=====//
    build_community(community_pathtofile); // UPDATE: there could be plenty extra files to be store somewhere
    mcmc.build_data(database_pathtofile);

    //=====Load Parameters=====//
    mcmc.build_parameters(parameters_pathtofile);
    
    //=====Assign classes from the model=====//
    mcmc.pre_process_data();
    
    //=====EXECUTION=====//

    //=====Parameter Inference=====//
    if (infer)
    {
        mcmc.inference(output_pathtofile, rand_init_params);
        if (augmented_pathtofile != UNDEFINED)
        {
            mcmc.print_augm_variables(augmented_pathtofile);
        }
    }

    //=====Simulation Study=====//
    if (simulate)
    {
        mcmc.simulation(output_pathtofile);
    }

    if (Verbose)
    {
        cout << endl;
    }

    return EXIT_SUCCESS;
}
