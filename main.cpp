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
vector<double> Cumul_Comm_Infec = {};
vector<double> Inst_Comm_Infec = {};

// ////////////////////////////////////////////////////////////////////////// //
// ///////////////////////// AUXILIARY FUNCTIONS //////////////////////////// //
// ////////////////////////////////////////////////////////////////////////// //

/* Cleaner string from forbidden characters*/
void string_cleaner(string *token)
{
    for (char r : {'\r', '\n'})
    {
        token->erase(remove(token->begin(), token->end(), r), token->end());
    }
}

/* Return the dictionary from the string */
kwargs parse_dict_string(string args_str, char delim)
{
    istringstream args_ss(args_str);
    string tmp;
    vector<string> args_stock;
    kwargs output = {};

    while (getline(args_ss, tmp, delim))
    {
        args_stock.push_back(tmp);
    }
    if (args_stock.size() == 0)
    {
        output = {}; // empty dictionary
    }
    else if (args_stock.size() % 2 == 1)
    {
        cerr << "\nERROR: Unrecognized dictionary." << endl;
        cerr << "Exit...\n"
             << endl;
        exit(EXIT_FAILURE);
    }
    else
    {
        for (size_t i = 0; i < args_stock.size(); i += 2)
        {
            try
            {
                output.emplace(args_stock[i], stod(args_stock[i + 1]));
            }
            catch (const exception &)
            {
                cerr << "\nERROR: Invalid entry of a dictionary." << endl;
                cerr << "Exit...\n"
                     << endl;
                exit(EXIT_FAILURE);
            }
        }
    }

    return output;
}

/* function to load boolean variables */
bool stobool(string text)
{
    bool output = false;

    if (text == "TRUE" || text == "true" || text == "True" || text == "1")
    {
        output = true;
    }
    else if (text == "FALSE" || text == "false" || text == "False" || text == "0")
    {
        output = false;
    }
    else
    {
        cerr << "\nERROR: Unrecognized boolean value for: " << text << endl;
        cerr << "Exit...\n"
             << endl;
        exit(EXIT_FAILURE);
    }

    return output;
}

/* function to load delimiter */
char load_delim(string text)
{
    char delim = '_';

    if (text.front() == '"' && text.back() == '"')
    {
        // Remove the double quotes and store it as a string
        delim = text.substr(1, text.length() - 2)[0];
    }
    else
    {
        cerr << "\nERROR: Unable to extract delimiter." << endl;
        cerr << "Check: " << text << "\n"
             << endl;
        cerr << "Exit...\n"
             << endl;
        exit(EXIT_FAILURE);
    }

    return delim;
}

/* auxiliary function to call for error in build parameters */
void build_column_error(string var_name, string &line)
{
    cerr << "\nERROR: Unable to extract `" << var_name << "` column." << endl;
    cerr << "Check line: " << line << "\n"
         << endl;
    cerr << "Exit...\n"
         << endl;
    exit(EXIT_FAILURE);
}

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

/* Read Data file */
void build_data(MCMC *mcmc, string data_pathtofile)
{
    ifstream input(data_pathtofile.c_str());
    if (input.is_open())
    {
        if (Verbose)
        {
            cout << "\n=============DATA===============" << endl;
            cout << "Data file: " << data_pathtofile << endl;
        }

        //=====Read column names for variables in the model=====//
        string line;
        getline(input, line);
        istringstream iss(line);
        string token;
        unordered_map<string, double> data = {};

        // Tokenize the line using the specified delimiter
        while (getline(iss, token, mcmc->setup.delim))
        {
            string_cleaner(&token);
            mcmc->data_headers.push_back(token);
            data.emplace(token, 0);
        }

        //=====Load individuals from the raw with values=====//
        for (string line; getline(input, line);)
        {
            istringstream iss(line);

            // Store information in variables
            size_t n_col = 0;
            while (getline(iss, token, mcmc->setup.delim))
            {
                data.at(mcmc->data_headers[n_col]) = stod(token);
                n_col++;
            }

            // Create an onbject Individual and add it to the model
            Individual this_indiv(data.at("status"), data.at("index"), data);
            mcmc->add_individual(this_indiv);
        }

        //=====assign structure of the community=====//
        mcmc->assign_community_structure();

        // ===== Display ===== //
        if (Verbose)
        {
            cout << "N. of individuals : " << mcmc->num_indiv << endl;
            cout << "N. of households : " << mcmc->num_househ << endl;
            cout << "N. of cases : " << mcmc->num_symptomatic << endl;
        }
    }
    else
    {
        cerr << "\nERROR: Impossible to open: " << data_pathtofile << endl;
        cerr << "Exit...\n"
             << endl;
        exit(EXIT_FAILURE);
    }
}

/* Read Model Parameters file */
void build_parameters(MCMC *mcmc, string data_pathtofile)
{
    char delim = mcmc->setup.delim;
    char dict_delim = mcmc->setup.dict_delim;
    ifstream input(data_pathtofile.c_str());

    if (input.is_open())
    {
        if (Verbose)
        {
            cout << "\n=============PARAMETERS==============" << endl;
            cout << "Model file: " << data_pathtofile << endl;
        }

        // Create a stringstream of the current line
        for (string line; getline(input, line);)
        {
            string name, value_str, infer, strategy, rate_str, range_str, prior_name, prior_args_str;
            double value;
            kwargs prior_args = {}, rate = {}, range = {};
            istringstream iss(line);

            //=====NAME=====//
            if (!getline(iss, name, delim))
            {
                build_column_error("name", line);
            }

            //=====VALUE=====//
            if (!getline(iss, value_str, delim))
            {
                build_column_error("value", line);
            }
            value = stod(value_str);

            //=====INFER=====//
            if (!getline(iss, infer, delim))
            {
                build_column_error("infer", line);
            }

            if (stobool(infer) == true)
            {
                //=====PRIOR=====//
                if (!getline(iss, prior_name, delim))
                {
                    build_column_error("prior", line);
                }

                //=====ARGS=====//
                if (!getline(iss, prior_args_str, delim))
                {
                    build_column_error("args", line);
                }
                prior_args = parse_dict_string(prior_args_str, dict_delim);

                //=====STRATEGY=====//
                if (!getline(iss, strategy, delim))
                {
                    build_column_error("strategy", line);
                }

                //=====STRATEGY_ARGS=====//
                if (!getline(iss, rate_str, delim))
                {
                    build_column_error("strat_args", line);
                }
                rate = parse_dict_string(rate_str, dict_delim);

                //=====RANGE=====//
                if (!getline(iss, range_str, delim))
                {
                    build_column_error("range", line);
                }
                range = parse_dict_string(range_str, dict_delim);

                // add entry to shared parameter dictionary
                mcmc->add_param(name, value, prior_name, prior_args, strategy, rate, range);
            }
            else
            {
                // add entry to shared parameter dictionary
                mcmc->add_param(name, value);
            }
        }

        if (Verbose)
        {
            mcmc->display_parameters();
        }
    }
    else
    {
        cerr << "\nERROR: Impossible to read `" << data_pathtofile << "`" << endl;
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
    bool infer = false;
    bool simulate = false;
    Setup_MCMC_Struct setup = {};
    string parameters_pathtofile = UNDEFINED;
    string community_pathtofile = UNDEFINED;
    string database_pathtofile = UNDEFINED;
    string output_pathtofile = UNDEFINED;
    string augmented_pathtofile = UNDEFINED;

    // Command line argument iterator
    int carg_i = 1;

    //=====Test if some commands were supplied=====//
    if (argc < 2)
    {
        cerr << "The user did not supply any command.\n";
        cerr << "Check the list of commands through:"
             << argv[0] << " -h\n";
        cerr << "Exit...\n"
             << endl;
        exit(EXIT_FAILURE);
    }

    //=====Load commands and options=====//
    while (carg_i < argc)
    {
        // Command line argument asking for help
        if (string(argv[carg_i]) == string("-h") or string(argv[carg_i]) == string("-help"))
        {
            // FIXME: Show the manual using man and installed manpages
            cerr << "Manual missing. Contact the developer." << endl;
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

        else if (string(argv[carg_i]) == string("-p") or string(argv[carg_i]) == string("--params"))
        {
            ++carg_i;
            parameters_pathtofile = argv[carg_i];
        }

        else if (string(argv[carg_i]) == string("-c") or string(argv[carg_i]) == string("--community"))
        {
            ++carg_i;
            community_pathtofile = argv[carg_i];
        }

        else if (string(argv[carg_i]) == string("-d") or string(argv[carg_i]) == string("--database"))
        {
            ++carg_i;
            database_pathtofile = argv[carg_i];
        }

        else if (string(argv[carg_i]) == string("-o") or string(argv[carg_i]) == string("--output"))
        {
            ++carg_i;
            output_pathtofile = argv[carg_i];
        }

        else if (string(argv[carg_i]) == string("-s") or string(argv[carg_i]) == string("--setup"))
        {
            ++carg_i;
            setup = load_setup(argv[carg_i]);
        }

        else if (string(argv[carg_i]) == string("-a") or string(argv[carg_i]) == string("--save_augmented"))
        {
            ++carg_i;
            augmented_pathtofile = argv[carg_i];
        }

        else if (string(argv[carg_i]) == string("-v") or string(argv[carg_i]) == string("--verbose"))
        {
            Verbose = true;
        }

        else
        {
            cerr << "Unrecognized option :" << argv[carg_i] << endl;
            cerr << "Exit...\n"
                 << endl;
            exit(EXIT_FAILURE);
        }

        // Read the next command line argument
        ++carg_i;
    }

    // UPDATE : check for options validity

    if (infer && simulate)
    {
        cerr << "Incompatible options: choose one between -infer and -simulate." << endl;
        cerr << "Exit...\n"
             << endl;
        exit(EXIT_FAILURE);
    }

    // UPDATE : check for path validity

    //=====INITIALIZATION=====//

    //=====Initialize MCMC=====//
    MCMC mcmc(setup);

    //=====Load Data=====//
    build_community(community_pathtofile);
    build_data(&mcmc, database_pathtofile);

    //=====Load Parameters=====//
    build_parameters(&mcmc, parameters_pathtofile);

    //=====Assign Classes from the Model=====//
    mcmc.pre_process_data();

    //=====Parameter Inference=====//
    if (infer)
    {
        mcmc.run(output_pathtofile);
        if (augmented_pathtofile != UNDEFINED)
        {
            mcmc.print_augm_variables(augmented_pathtofile);
        }
    }

    //=====Simulation Study=====//
    if (simulate)
    {
        mcmc.simulate(output_pathtofile);
    }

    if (Verbose)
    {
        cout << endl;
    }

    return EXIT_SUCCESS;
}
