#include <string>
#include <vector>
#include <unordered_map>
#include <cmath>
#include <random>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>

#include "mcmc_class.hpp"

using namespace std;

// Define the global variables here
bool Verbose = false;
std::unordered_map<std::string, Parameter> Parameter_Dictionary = {};
std::vector<double> Cumul_Community_Infec = {};
std::vector<double> Inst_Community_Infec = {};


// ////////////////////////////////////////////////////////////////////////// //
// ////////////////////////////// FILE LOADING ////////////////////////////// //
// ////////////////////////////////////////////////////////////////////////// //

/*
    Read Community file
    ---------
*/
void build_community (std::string data_file) {
    double date, inst_val, cumul_val;

    std::ifstream infile(data_file.c_str());
    if (infile) {
        std::cout << "=============COMMUNITY==============" << std::endl;
        std::cout << "Read data file: " << data_file << std::endl;

        // Create a stringstream of the current line
        for (std::string line; std::getline(infile, line);) {
            std::istringstream iss(line);
            
            // Store information in variables 
            iss >> date >> inst_val >> cumul_val;

            Inst_Community_Infec.push_back(inst_val);
            Cumul_Community_Infec.push_back(cumul_val);
        }

    } else {
        std::cout << "ERROR: Impossible to read: " << data_file << std::endl;
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
    1: indid (int), 2: this_HH_id (int), 3: HH_size (int), 4: onset_time (double), 5: study_period (double), 6: is_case (0 to 3), 
    7: vacc_status (0 to 3), 8:  age (int), 9: isolation (0 or 1), 10: sex (0 or 1), 11: vacc_delay (int), 12: recovered (0 or 1)

    is_case :
    ---------
    0: not a case
    1: symptomatic case
    2: asymptomatic case
    3: symptomatic case with unknown symptom onset
*/
std::vector<Household> build_households (std::string data_file) {
    int this_HH_id, indid, age, sex;
    int HH_size, is_case, vacc_status;
    double onset_time, study_period, vacc_delay;
    bool isolation = false, recovered = false;
    int num_cases = 0, num_subjects = 0;
    double max_study_period = 0.0;
    std::unordered_map<int, Household> safe_load = {};
    std::vector<Household> output = {};

    std::ifstream infile(data_file.c_str());
    if (infile) {
        std::cout << "=============DATA===============" << std::endl;
        std::cout << "Read data file: " << data_file << std::endl;

        // Create a stringstream of the current line
        for (std::string line; std::getline(infile, line);) {
            std::istringstream iss(line);
            
            // Store information in variables 
            iss >> indid >> this_HH_id >> HH_size >> onset_time >> study_period >> is_case >> vacc_status >> age >> isolation >> sex >> vacc_delay >> recovered;

            // if this_HH_id doesn't exist, append new household
            if (safe_load.find(this_HH_id) == safe_load.end()) {
                safe_load.emplace(this_HH_id, Household());
            } 

            // Update household object and information parameters
            Individual this_indiv (
                indid,
                onset_time,
                is_case,
                age,
                sex, 
                study_period,
                vacc_status,
                isolation,
                vacc_delay,
                recovered
                );

            safe_load[this_HH_id].add_individual(this_indiv);
            num_subjects++;
            if (is_case != 0){
                num_cases++;
            }
            max_study_period = max(max_study_period, study_period);
        }

    } else {
        std::cout << "ERROR: Impossible to read: " << data_file << std::endl;
        exit(EXIT_FAILURE);
    }

    // extract only households and forget about their id
    for (const auto& entry : safe_load) {
        output.push_back(entry.second);
    }

    // General information
    if (Verbose == true){
        std::cout << "Number of households: " << output.size() << std::endl;
        std::cout << "Number of individuals: " << num_subjects << std::endl;
        std::cout << "Number of cases: " << num_cases << std::endl;
        std::cout << "Number of days: " << max_study_period << "\n" << std::endl;
    }

    return(output);
}

/*
    Read Model Parameters file
    ---------

void build_parameters (std::string data_file) {   
    std::string name, infer, prior, args, strategy, range;
    double value, rate;

    std::ifstream infile(data_file.c_str());
    if (infile) {
        std::cout << "=============PARAMETERS==============" << std::endl;
        std::cout << "Read data file: " << data_file << std::endl;

        // Create a stringstream of the current line
        for (std::string line; std::getline(infile, line);) {
            std::istringstream iss(line);

            // Read mandatory variables
            if (!(iss >> name >> value >> infer)) {
                std::cout << "Unable to extract required columns." << std::endl;
                exit(EXIT_FAILURE);
            }

            // Convert "infer" to boolean
            if (infer == "TRUE" || infer == "true" || infer == "True"){
                if (!(iss >> args >> strategy >> rate >> range)) {
                    std::cout << "Unable to extract required columns." << std::endl;
                    exit(EXIT_FAILURE);
                }


                Parameter this_param (name, value, true);
                // add entry to shared paramter dictionary
                Parameter_Dictionary.emplace(name, this_param);

            } else if (infer == "FALSE" || infer == "false" || infer == "False"){
                Parameter this_param (name, value, false);
                // add entry to shared paramter dictionary
                Parameter_Dictionary.emplace(name, this_param);

            } else {
                std::cout << "Unrecognized format value for infer." << std::endl;
                exit(EXIT_FAILURE);
            }
        }   
    }
}
*/

// /////////////////////////////////////////////////////////////////////////// //
// ////////////////////////////// MAIN FUNCTION ////////////////////////////// //
// /////////////////////////////////////////////////////////////////////////// //


int main (int argc, char **argv) {
           
    int num_sweeps = 20; 
    int num_adapt_iter = 0;
    int num_augm_iter = 1;
    unsigned seed = DEFAULT_SEED; // FIXME : this should be a user optional variable  

    //========== USER parameters ==========
    if (argc < 6){
        // print usage here and return, since that's what you need.
        std::cout << "Wrong usage" << std::endl;
        std::exit(EXIT_FAILURE);
    }

    std::string households_file = argv[1];      // Path/to/file (households)
    std::string community_file = argv[2];       // Path/to/file (community)
    std::string output_file = argv[3];          // Name of the output file
    int num_iter = std::stoi(argv[4]);          // N. of iterations in the MCMC

    // --- Temporary --- FIXME: bool from user 
    std::string tmp = argv[5];
    if ((tmp == "yes") || (tmp == "true") || (tmp == "True") || (tmp == "TRUE")) {
        Verbose = true;
    } else {
        Verbose = false;
    }
    // --- End-of-Temporary ---

    // FIXME : check for space and path validity

    //========== Model parameters ==========

    // --- Temporary --- FIXME : load parameters from file  // need of iteration scopes
    {Parameter this_param ("alpha", 0.001, true, "unif", {0.0, 1.0}, "exp", 0.7, {0.0, 1.0});
    Parameter_Dictionary.emplace(this_param.p_name, this_param);    }
    {Parameter this_param ("beta", 0.6, true, "unif", {0.0, 5.0}, "exp", 0.1,  {0.0, 5.0}); 
    Parameter_Dictionary.emplace(this_param.p_name, this_param);}
    {Parameter this_param("delta", 1.0, true, "unif", {-3.0, 3.0}, "norm", 0.2, {-3.0, 3.0});
    Parameter_Dictionary.emplace(this_param.p_name, this_param);}
    {Parameter this_param("rU", 1.0, true, "lognorm", {0.0, 1.0}, "exp", 0.3, {0.0, 5.0});
    Parameter_Dictionary.emplace(this_param.p_name, this_param);}
    {Parameter this_param("r1V", 1.0, true, "lognorm", {0.0, 1.0}, "exp", 1.0, {0.0, 5.0});
    Parameter_Dictionary.emplace(this_param.p_name, this_param);}
    {Parameter this_param("r1V90", 1.0, true, "lognorm", {0.0, 1.0}, "exp", 1.0, {0.0, 0.1});
    Parameter_Dictionary.emplace(this_param.p_name, this_param);}
    {Parameter this_param("r2V", 1.0, true, "lognorm", {0.0, 1.0}, "exp", 1.0, {0.0, 0.1});
    Parameter_Dictionary.emplace(this_param.p_name, this_param);}
    {Parameter this_param("r3V", 1.0, true, "lognorm", {0.0, 1.0}, "exp", 0.9, {0.0, 5.0});
    Parameter_Dictionary.emplace(this_param.p_name, this_param);}
    {Parameter this_param("rC", 1.0, true, "lognorm", {0.0, 1.0}, "exp", 0.1, {0.0, 5.0});
    Parameter_Dictionary.emplace(this_param.p_name, this_param);}
    {Parameter this_param("rQ", 1.0, true, "lognorm", {0.0, 1.0}, "exp", 0.2, {0.0, 5.0});
    Parameter_Dictionary.emplace(this_param.p_name, this_param);}
    {Parameter this_param("rCC", 1.0, false, "lognorm", {0.0, 1.0}, "exp", 0.5, {0.0, 5.0});
    Parameter_Dictionary.emplace(this_param.p_name, this_param);}
    {Parameter this_param("rMC", 1.0, true, "lognorm", {0.0, 1.0}, "exp", 0.5, {0.0, 5.0});
    Parameter_Dictionary.emplace(this_param.p_name, this_param);}
    {Parameter this_param("rFC", 1.0, true, "lognorm", {0.0, 1.0}, "exp", 0.5, {0.0, 5.0});
    Parameter_Dictionary.emplace(this_param.p_name, this_param);}
    {Parameter this_param("rR", 1.0, true, "lognorm", {0.0, 1.0}, "exp", 1.0, {0.0, 5.0});
    Parameter_Dictionary.emplace(this_param.p_name, this_param);}

    {Parameter this_param("tost_tau", 5.783, false);
    Parameter_Dictionary.emplace(this_param.p_name, this_param);}
    {Parameter this_param("tost_mu", -4.0, false);
    Parameter_Dictionary.emplace(this_param.p_name, this_param);}
    {Parameter this_param("tost_sigma", 1.85, false);
    Parameter_Dictionary.emplace(this_param.p_name, this_param);}
    {Parameter this_param("tost_alpha", 5.85, false);
    Parameter_Dictionary.emplace(this_param.p_name, this_param);}

    {Parameter this_param("incub_shape", 3.44898, false);
    Parameter_Dictionary.emplace(this_param.p_name, this_param);}
    {Parameter this_param("incub_scale", 1.13077, false);
    Parameter_Dictionary.emplace(this_param.p_name, this_param);}
    {Parameter this_param("typical_size", 5.0, false);
    Parameter_Dictionary.emplace(this_param.p_name, this_param);}
    // --- End-of-Temporary ---

    //==========Build data==========
    // Load data
    build_community(community_file);
    std::vector<Household> HH_data = build_households(households_file);

    //==========MCMC==========

    // Initialize MCMC object
    MCMC_Class mcmc_object (
        seed, 
        num_iter, 
        num_adapt_iter,
        HH_data,
        num_augm_iter
        );

    mcmc_object.run(output_file, num_sweeps);

    return EXIT_SUCCESS;
}
