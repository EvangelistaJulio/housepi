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

//=====Define the global variables=====//
// UPDATE: this depends on the model
std::vector<double> Cumul_Community_Infec = {};
std::vector<double> Inst_Community_Infec = {};

// ////////////////////////////////////////////////////////////////////////// //
// ///////////////////////// AUXILIARY FUNCTIONS //////////////////////////// //
// ////////////////////////////////////////////////////////////////////////// //

kwargs parse_dict_string(std::string args_str, char delim)
{
    std::istringstream args_ss(args_str);
    std::string tmp;
    std::vector<std::string> args_stock;
    kwargs output = {};

    while (std::getline(args_ss, tmp, delim))
    {
        args_stock.push_back(tmp);
    }
    if (args_stock.size() == 0)
    {
        output = {}; // empty dictionary
    }
    else if (args_stock.size() % 2 == 1)
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

/* function to load boolean variables */
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

/* function to load delimiter */
char load_delim(std::string text)
{
    char delim = '_';

    if (text.front() == '"' && text.back() == '"')
    {
        // Remove the double quotes and store it as a string
        delim = text.substr(1, text.length() - 2)[0];
    }
    else
    {
        std::cerr << "\nERROR: Unable to extract delimiter." << std::endl;
        std::cerr << "Check: " << text << "\n"
                  << std::endl;
        std::cerr << "Exit...\n"
                  << std::endl;
        exit(EXIT_FAILURE);
    }

    return delim;
}

/* auxiliary function to call for error in build parameters */
void build_column_error(std::string var_name, std::string &line)
{
    std::cerr << "\nERROR: Unable to extract `" << var_name << "` column." << std::endl;
    std::cerr << "Check line: " << line << "\n"
              << std::endl;
    std::cerr << "Exit...\n"
              << std::endl;
    exit(EXIT_FAILURE);
}

// ////////////////////////////////////////////////////////////////////////// //
// ////////////////////////////// FILE LOADING ////////////////////////////// //
// ////////////////////////////////////////////////////////////////////////// //

/* Read Setup file */
Setup_MCMC_Struct load_setup(std::string data_path_to_file)
{
    char delim = ',';
    Setup_MCMC_Struct loader;
    std::ifstream input_file(data_path_to_file.c_str());
    std::map<std::string, std::string> safe_load = {};

    if (input_file)
    {
        std::cout << "\n=============SETUP==============" << std::endl;
        std::cout << "Setup file: " << data_path_to_file << std::endl;

        input_file.clear();                 // Clear any errors
        input_file.seekg(0, std::ios::beg); // Seek to the beginning

        // Create a stringstream of the current line
        for (std::string line; std::getline(input_file, line);)
        {
            std::string name, choice;
            std::istringstream iss(line);

            //=====NAME=====//
            if (!getline(iss, name, delim))
            {
                std::cerr << "\nERROR: Unable to extract setup name." << std::endl;
                std::cerr << "Check line: " << line << "\n"
                          << std::endl;
                std::cerr << "Exit...\n"
                          << std::endl;
                exit(EXIT_FAILURE);
            }

            //=====CHOICE=====//
            if (!getline(iss, choice, '\r')) // FIXME: this might be just a problem of my editor
            {
                std::cerr << "\nERROR: Unable to extract setup choice." << std::endl;
                std::cerr << "Check line: " << line << "\n"
                          << std::endl;
                std::cerr << "Exit...\n"
                          << std::endl;
                exit(EXIT_FAILURE);
            }

            // Trim leading and trailing spaces from key and value
            name.erase(0, name.find_first_not_of(" "));
            name.erase(name.find_last_not_of(" ") + 1);
            choice.erase(0, choice.find_first_not_of(" "));
            choice.erase(choice.find_last_not_of(" ") + 1);

            // add the new setup
            safe_load.emplace(name, choice);
        }

        loader.seed = std::stoi(safe_load.at("seed"));
        loader.num_iter = std::stoi(safe_load.at("num_iter"));
        loader.num_adapt_iter = std::stoi(safe_load.at("num_adapt_iter"));
        loader.thinning = std::stoi(safe_load.at("thinning"));
        loader.is_data_augm = stobool(safe_load.at("is_data_augm"));
        loader.verbose = stobool(safe_load.at("verbose"));
        loader.delim = load_delim(safe_load.at("delim"));
        loader.dict_delim = load_delim(safe_load.at("dict_delim"));
    }
    else
    {
        std::cerr << "\nERROR: Impossible to read: " << data_path_to_file << std::endl;
        std::cerr << "Exit...\n"
                  << std::endl;
        exit(EXIT_FAILURE);
    }

    //=====Check for values validity=====//
    if (loader.num_adapt_iter >= loader.num_iter)
    {
        std::cerr << "\nERROR: N. of adapative iteration >= N. of iterations." << std::endl;
        std::cerr << "Exit...\n"
                  << std::endl;
        exit(EXIT_FAILURE);
    }

    //=====Display=====//
    if (loader.verbose)
    {
        for (const auto &pair : safe_load)
        {
            std::cout << pair.first << " : " << pair.second << std::endl;
        }
    }

    return loader;
}

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
std::vector<Household> build_households(std::string data_path_to_file, bool verbose)
{
    int this_HH_id, indid, age, sex;
    int HH_size, status, vacc_status;
    timetype onset_time, study_period, vacc_delay;
    bool isolation = false, recovered = false;
    int num_cases = 0, num_subjects = 0;
    timetype max_study_period = 0.0;
    std::unordered_map<int, Household> safe_load_HH = {};
    std::vector<Household> data_HH = {};

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
            if (safe_load_HH.find(this_HH_id) == safe_load_HH.end())
            {
                safe_load_HH.emplace(this_HH_id, Household());
                safe_load_HH.at(this_HH_id).set_index(this_HH_id);
            }

            // UPDATE : to be generalized starting from the model
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
                age,
                sex,
                isolation,
                recovered,
                vacc_status,
                study_period,
                vacc_delay,
                onset_time);
            safe_load_HH[this_HH_id].add_individual(this_indiv, status);

            // dataset composition
            if (verbose)
            {
                num_subjects++;
                num_cases += (status == HEALTHY) ? 0 : 1;
                max_study_period = max(max_study_period, study_period);
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

    // extract only households and forget about their id
    for (const auto &entry : safe_load_HH)
    {
        data_HH.push_back(entry.second);
    }

    //=====Display=====//
    if (verbose)
    {
        std::cout << "N. of households : " << data_HH.size() << std::endl;
        std::cout << "N. of individuals : " << num_subjects << std::endl;
        std::cout << "N. of cases : " << num_cases << std::endl;
        std::cout << "N. of days : " << max_study_period << std::endl;
    }

    return (data_HH);
}

/* Read Model Parameters file */
void build_parameters(MCMC *mcmc, std::string data_path_to_file, char delim, char dict_delim, bool verbose)
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
    }
    else
    {
        std::cerr << "\nERROR: Impossible to read `" << data_path_to_file << "`" << std::endl;
        std::cerr << "Exit...\n"
                  << std::endl;
        exit(EXIT_FAILURE);
    }

    //=====Display=====//
    if (verbose)
    {
        mcmc->print();
    }
}

// /////////////////////////////////////////////////////////////////////////// //
// ////////////////////////////// MAIN FUNCTION ////////////////////////////// //
// /////////////////////////////////////////////////////////////////////////// //

int main(int argc, char **argv)
{
    //========== USER parameters ==========//
    if (argc < 6)
    {
        // print usage here and return, since that's what you need.
        std::cerr << "Wrong usage" << std::endl;
        std::cerr << argv[0] << " <households:path/to/file>"
                  << " <community:path/to/file>"
                  << " <parameters:path/to/file>"
                  << " <output:path/to/file>"
                  << " <augmented:path/to/file>";
        std::cerr << std::endl
                  << "Exit...\n"
                  << std::endl;
        std::exit(EXIT_FAILURE);
    }

    // mandatory parameters
    std::string households_path_to_file = argv[1]; // Path/to/file (households)
    std::string community_path_to_file = argv[2];  //
    std::string parameters_path_to_file = argv[3]; //
    std::string output_path_to_file = argv[4];     // Name of the output file
    std::string augmented_path_to_file = argv[5];  // Name of the output file for augmented data

    std::string setup_path_to_file = "default_setup.csv"; // UPDATE : these should be user optional variables

    // UPDATE : check for path validity

    //=====Initialization=====//

    // Load Setup
    Setup_MCMC_Struct setup = load_setup(setup_path_to_file);

    // Load Data
    build_community(community_path_to_file);
    std::vector<Household> HH_data = build_households(households_path_to_file, setup.verbose);

    // Initialize MCMC object
    MCMC mcmc(setup, HH_data);

    // Load Parameters
    build_parameters(&mcmc, parameters_path_to_file, setup.delim, setup.dict_delim, setup.verbose);

    // Initialize infectivity profile
    mcmc.assign_infec_profile();

    //=====Run Monte Carlo Markov Chain=====//
    mcmc.run(output_path_to_file);

    // print last augmented data configuration // FIXME : temporary
    mcmc.open_augm_file(augmented_path_to_file);
    mcmc.print_augm_variables();
    mcmc.close_augm_file();

    return EXIT_SUCCESS;
}
