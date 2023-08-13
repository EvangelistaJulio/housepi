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

#include "mcmc.hpp"

using namespace std;

// --- Define the global variables here ---

bool Verbose = true;
Parameter_Dictionary glob_params = {};
std::vector<double> Cumul_Community_Infec = {};
std::vector<double> Inst_Community_Infec = {};

// ////////////////////////////////////////////////////////////////////////// //
// ///////////////////////// AUXILIARY FUNCTIONS //////////////////////////// //
// ////////////////////////////////////////////////////////////////////////// //

kwargs parse_dict_string(std::string args_str, char delimiter = ',')
{
    std::istringstream args_ss(args_str);
    std::string tmp;
    std::vector<std::string> args_stock;
    kwargs output;

    while (std::getline(args_ss, tmp, delimiter))
    {
        args_stock.push_back(tmp);
    }
    if (args_stock.size() == 0 || args_stock.size() % 2 == 1)
    {
        std::cerr << "\nERROR: Unrecognized dictionary." << std::endl;
        std::cerr << "Exit...\n"
                  << std::endl;
        exit(EXIT_FAILURE);
    }
    else
    {
        for (size_t i = 0; i < args_stock.size(); i += 2)
        {
            try
            {
                output.emplace(args_stock[i], std::stod(args_stock[i + 1]));
            }
            catch (const std::exception &)
            {
                std::cerr << "\nERROR: Invalid entry of a dictionary." << std::endl;
                std::cerr << "Exit...\n"
                          << std::endl;
                exit(EXIT_FAILURE);
            }
        }
    }

    return output;
}

bool stobool(std::string text)
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
        std::cerr << "\nERROR: Unrecognized boolean value for: " << text << std::endl;
        std::cerr << "Exit...\n"
                  << std::endl;
        exit(EXIT_FAILURE);
    }

    return output;
}

// ////////////////////////////////////////////////////////////////////////// //
// ////////////////////////////// FILE LOADING ////////////////////////////// //
// ////////////////////////////////////////////////////////////////////////// //

/* Read Community file */
void build_community(std::string data_path_to_file)
{
    double date, inst_val, cumul_val;

    std::ifstream input_file(data_path_to_file.c_str());
    if (input_file)
    {
        std::cout << "\n=============COMMUNITY==============" << std::endl;
        std::cout << "Data file: " << data_path_to_file << std::endl;

        // Create a stringstream of the current line
        for (std::string line; std::getline(input_file, line);)
        {
            std::istringstream iss(line);

            // Store information in variables
            iss >> date >> inst_val >> cumul_val;

            Inst_Community_Infec.push_back(inst_val);
            Cumul_Community_Infec.push_back(cumul_val);
        }
    }
    else
    {
        std::cerr << "\nERROR: Impossible to read: " << data_path_to_file << std::endl;
        std::cerr << "Exit...\n"
                  << std::endl;
        exit(EXIT_FAILURE);
    }
}

/*
    Read Households file
    ---------
    The Households file should be a space-separated table with household ids
    sorted in ascending order and numbered from 0 to n-1 households

    datafile with columns :
    -----------------------
    1: indid (int), 2: this_HH_id (int), 3: HH_size (int), 4: onset_time (double), 5: study_period (double), 6: status (0 to 3),
    7: vacc_status (0 to 3), 8:  age (int), 9: isolation (0 or 1), 10: sex (0 or 1), 11: vacc_delay (int), 12: recovered (0 or 1)

    status :
    ---------
    0: not a case
    1: symptomatic case with known symptom onset
    2: asymptomatic case
    3: symptomatic case with unknown symptom onset
*/
std::vector<Household> build_households(std::string data_path_to_file)
{
    int this_HH_id, indid, age, sex;
    int HH_size, status, vacc_status;
    double onset_time, study_period, vacc_delay;
    bool isolation = false, recovered = false;
    int num_cases = 0, num_subjects = 0;
    double max_study_period = 0.0;
    std::unordered_map<int, Household> safe_load = {};
    std::vector<Household> output = {};

    std::ifstream input_file(data_path_to_file.c_str());
    if (input_file)
    {
        std::cout << "\n=============DATA===============" << std::endl;
        std::cout << "Data file: " << data_path_to_file << std::endl;

        // Create a stringstream of the current line
        for (std::string line; std::getline(input_file, line);)
        {
            std::istringstream iss(line);

            // Store information in variables
            iss >> indid >> this_HH_id >> HH_size >> onset_time >> study_period >> status >> vacc_status >> age >> isolation >> sex >> vacc_delay >> recovered;

            // if this_HH_id doesn't exist, append new household
            if (safe_load.find(this_HH_id) == safe_load.end())
            {
                safe_load.emplace(this_HH_id, Household());
            }

            if ((status != HEALTHY) && (status != SYMPTOMATIC))
            {
                std::cerr << "ERROR: status not in the model: " << status << endl;
                std::cerr << "Exit...\n"
                          << std::endl;
                std::exit(EXIT_FAILURE);
            }

            // Update household object and information parameters
            Individual this_indiv(
                indid,
                status,
                age,
                sex,
                isolation,
                recovered,
                vacc_status,
                study_period,
                vacc_delay,
                onset_time);

            safe_load[this_HH_id].add_individual(this_indiv);
            num_subjects++;
            num_cases += (status == HEALTHY) ? 0 : 1;
            max_study_period = max(max_study_period, study_period);
        }
    }
    else
    {
        std::cerr << "\nERROR: Impossible to read: " << data_path_to_file << std::endl;
        std::cerr << "Exit...\n"
                  << std::endl;
        exit(EXIT_FAILURE);
    }

    // extract only households and forget about their id
    for (const auto &entry : safe_load)
    {
        output.push_back(entry.second);
    }

    // General information
    if (Verbose == true)
    {
        std::cout << "N. of households :" << TAB << output.size() << std::endl;
        std::cout << "N. of individuals :" << TAB << num_subjects << std::endl;
        std::cout << "N. of cases :" << TAB << num_cases << std::endl;
        std::cout << "N. of days :" << TAB << max_study_period << std::endl;
    }

    return (output);
}

/* Read Model Parameters file */
void build_parameters(std::string data_path_to_file, char delimiter = '\t')
{
    std::ifstream input_file(data_path_to_file.c_str());

    if (input_file)
    {
        std::cout << "\n=============PARAMETERS==============" << std::endl;
        std::cout << "Model file: " << data_path_to_file << std::endl;

        // Create a stringstream of the current line
        for (std::string line; std::getline(input_file, line);)
        {
            std::string name, value_str, infer, strategy, rate_str, range_str, prior_name, prior_args_str;
            double value;
            kwargs prior_args = {}, rate = {}, range = {};
            std::istringstream iss(line);

            // --- NAME ---
            if (!getline(iss, name, delimiter))
            {
                std::cerr << "\nERROR: Unable to extract `name` column." << std::endl;
                std::cerr << "Check line: " << line << "\n"
                          << std::endl;
                std::cerr << "Exit...\n"
                          << std::endl;
                exit(EXIT_FAILURE);
            }

            // --- VALUE ---
            if (!getline(iss, value_str, delimiter))
            {
                std::cerr << "\nERROR: Unable to extract `value` column." << std::endl;
                std::cerr << "Check line: " << line << "\n"
                          << std::endl;
                std::cerr << "Exit...\n"
                          << std::endl;
                exit(EXIT_FAILURE);
            }
            value = stod(value_str);

            // --- INFER ---
            if (!getline(iss, infer, delimiter))
            {
                std::cerr << "\nERROR: Unable to extract `infer` column." << std::endl;
                std::cerr << "Check line: " << line << "\n"
                          << std::endl;
                std::cerr << "Exit...\n"
                          << std::endl;
                exit(EXIT_FAILURE);
            }

            if (stobool(infer) == true)
            {
                // --- PRIOR ---
                if (!getline(iss, prior_name, delimiter))
                {
                    std::cerr << "\nERROR: Unable to extract `prior` column." << std::endl;
                    std::cerr << "Check line: " << line << "\n"
                              << std::endl;
                    std::cerr << "Exit...\n"
                              << std::endl;
                    exit(EXIT_FAILURE);
                }

                // --- ARGS ---
                if (!getline(iss, prior_args_str, delimiter))
                {
                    std::cerr << "\nERROR: Unable to extract `args` column." << std::endl;
                    std::cerr << "Check line: " << line << "\n"
                              << std::endl;
                    std::cerr << "Exit...\n"
                              << std::endl;
                    exit(EXIT_FAILURE);
                }
                prior_args = parse_dict_string(prior_args_str);

                // --- STRATEGY ---
                if (!getline(iss, strategy, delimiter))
                {
                    std::cerr << "\nERROR: Unable to extract `strategy` column." << std::endl;
                    std::cerr << "Check line: " << line << "\n"
                              << std::endl;
                    std::cerr << "Exit...\n"
                              << std::endl;
                    exit(EXIT_FAILURE);
                }

                // --- RATE ---
                if (!getline(iss, rate_str, delimiter))
                {
                    std::cerr << "\nERROR: Unable to extract `strategy` column." << std::endl;
                    std::cerr << "Check line: " << line << "\n"
                              << std::endl;
                    std::cerr << "Exit...\n"
                              << std::endl;
                    exit(EXIT_FAILURE);
                }
                rate = parse_dict_string(rate_str);

                if (!getline(iss, range_str, delimiter))
                {
                    std::cerr << "\nERROR: Unable to extract `range` column." << std::endl;
                    std::cerr << "Check line: " << line << "\n"
                              << std::endl;
                    std::cerr << "Exit...\n"
                              << std::endl;
                    exit(EXIT_FAILURE);
                }
                range = parse_dict_string(range_str);

                // add entry to shared paramter dictionary
                Parameter this_param(value, true, strategy, rate, range, prior_name, prior_args);
                glob_params.add_param(name, this_param);
            }
            else
            {
                Parameter this_param(value, false);
                // add entry to shared paramter dictionary}
                glob_params.add_param(name, this_param);
            }
        }
    }
    else
    {
        std::cerr << "\nERROR: Impossible to read: " << data_path_to_file << std::endl;
        std::cerr << "Exit...\n"
                  << std::endl;
        exit(EXIT_FAILURE);
    }

    // --- Display --- //
    if (Verbose == true)
    {
        glob_params.print();
    }
    // --- End-Of-Display --- //
}

// /////////////////////////////////////////////////////////////////////////// //
// ////////////////////////////// MAIN FUNCTION ////////////////////////////// //
// /////////////////////////////////////////////////////////////////////////// //

int main(int argc, char **argv)
{

    // UPDATE : these should be user optional variables
    unsigned seed = DEFAULT_SEED;
    int num_iter = DEFAULT_NUM_ITER;
    int num_param_sweeps = DEFAULT_NUM_PARAM_SWEEEPS;
    int num_augm_sweeps = DEFAULT_NUM_AUGM_SWEEPS;
    int num_adapt_steps = DEFAULT_NUM_ADAPT_STEPS;

    //========== USER parameters ==========
    if (argc < 5)
    {
        // print usage here and return, since that's what you need.
        std::cerr << "Wrong usage" << std::endl;
        std::cerr << argv[0] << " <households:path/to/file>"
                  << " <community:path/to/file>"
                  << " <parameters:path/to/file>"
                  << " <output:path/to/file>"
                  << " <n.iterations:int>"
                  << " <verbose:bool>" << std::endl;
        std::cerr << "Exit...\n"
                  << std::endl;
        std::exit(EXIT_FAILURE);
    }

    // mandatory parameters
    std::string households_path_to_file = argv[1]; // Path/to/file (households)
    std::string community_path_to_file = argv[2];
    std::string parameters_path_to_file = argv[3];
    std::string output_path_to_file = argv[4]; // Name of the output file

    // optional parameters
    if (argc > 5)
    { // N. of iterations in the MCMC
        num_iter = std::stoi(argv[5]);
    }
    if (argc > 6)
    { // Verbose
        Verbose = stobool(argv[6]);
    }

    // UPDATE : check for space and path validity

    // Load data and Paramters
    std::vector<Household> HH_data = build_households(households_path_to_file);
    build_community(community_path_to_file);
    build_parameters(parameters_path_to_file);

    // Initialize MCMC object
    MCMC_Class mcmc_object(seed, HH_data);

    // run MCMC inference
    mcmc_object.infer(
        output_path_to_file,
        num_iter,
        num_param_sweeps,
        num_augm_sweeps,
        num_adapt_steps);

    return EXIT_SUCCESS;
}
