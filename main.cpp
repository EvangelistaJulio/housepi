#include <string>
#include <vector>
#include <unordered_map>
#include <cmath>
#include <random>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>

#include "mcmc_object.hpp"

#define DEFAULT_SEED 20210415

using namespace std;

// Define the global variable here
std::unordered_map<std::string, Parameter> Parameter_Dictionary = {};

//----------------------------------------------------------------------
// Load data
// datafile with columns :  indid (int) currhhid (int) hhsize (int) onsetTime (double) study_period (double) is_case (0, 1, 2, 3) vacc_status (0 to 3)  age (int) isolation (0 or 1) sex (0 or 1) vacc_delay (int) recovered (0 or 1)
/* 
    is_case
    0: not a case
    1: symptomatic case
    2: asymptomatic case
    3: symptomatic case with unknown symptom onset
*/
//----------------------------------------------------------------------
std::vector<Household> build_data(std::string data_file, std::string community_infec_file, bool contact_pattern) {
    Individual this_indiv; // FIXME : why this is not std::Individual ?
    std::vector<double> instinfcom;
    std::vector<double> cuminfcom;

    instinfcom.resize(0);
    cuminfcom.resize(0);

    std::ifstream infileinf(community_infec_file.c_str());
    if (infileinf) {
        cout << "=============INFILE==============" << endl;
        cout << "Read data file: " << data_file << endl;

        int date(0);
        double instval = 0.0;
        double cumval = 0.0;

        // while(std::getline(infile, line))
        for (std::string line; std::getline(infileinf, line);) {
            // Create a stringstream of the current line
            std::istringstream in(line);
            
            // Store information in variables 
            in >> date >> instval >> cumval;

            //cout << currhhid << endl;

            instinfcom.push_back(instval);
            cuminfcom.push_back(cumval);
        }
    }

    // Variable with all the data
    std::vector<Household> output;
    output.resize(0);

    // Read file
    /*
    The data file should be a space-separated table with household ids
    sorted in ascending order and numbered from 0 to n-1 households
    */
    std::ifstream infile(data_file.c_str());
    if (infile) {
        cout << "=============DATA===============" << endl;
        cout << "Read data file: " << data_file << endl;

        int indid(0), lasthhid(-1), currhhid(0), hhsize(0), is_case(0), start_follow_up(0), study_period(0), age(0),sex(0),vacc_status(0),isolation(0),vacc_delay(0),recovered(0);
        double onsetTime = 0.0;
        int num_cases = 0, num_households = 0, num_subjects = 0, num_days =0;
        // datafile with columns : 1: indid (int), 2: currhhid (int), 3: hhsize (int), 4: onsetTime (double), 5: study_period (double), 6: is_case (0 or 1), 7: vacc_status (0 to 3), 8:  age (int), 9: isolation (0 or 1), 10: sex (0 or 1), 11: vacc_delay (int), 12: recovered (0 or 1)
        //int nSymptomatic(0);

        // initialize current household
        Household currHH(contact_pattern, instinfcom, cuminfcom);
        for (std::string line; std::getline(infile, line);) {
            // Create a stringstream of the current line
            std::istringstream in(line);
            
            // Store information in variables 
            in >> indid >> currhhid >> hhsize >> onsetTime >> study_period >> is_case >> vacc_status >>  age >> isolation >>sex>>vacc_delay>>recovered;
            cout << recovered << endl;

            // Append previous household to the list of households
            if (currhhid != lasthhid) {
                if (num_households > 0) {
                    num_households++;
                    output.push_back(currHH);
                    // FIXME : check me out, am I providing an empty object of the class?
                    currHH = Household(contact_pattern, instinfcom, cuminfcom); 
                    lasthhid = currhhid;
                } else {
                    num_households++;
                    lasthhid = currhhid;
                }
            }

            // Update household object and information parameters
            this_indiv = Individual(
                indid,
                onsetTime,
                is_case,
                age,
                sex, 
                start_follow_up,
                study_period,
                vacc_status,
                isolation,
                vacc_delay,
                recovered
                );
            currHH.add_individual(this_indiv);
            num_subjects++;
            if (is_case != 0){
                num_cases++;
            }
            num_days = max(num_days, study_period);
        }

        // Add last household
        output.push_back(currHH);

        // General information
        cout << "Number of households: " << num_households << endl;
        cout << "Number of individuals: " << num_subjects << endl;
        cout << "Number of cases: " << num_cases << endl;
        cout << "Number of days: " << num_days << "\n\n";

        return(output);

    }else{
        cout << "ERROR: Cannot open the file." << endl;
        return(output);
    }
}


// ///////////////////////////////////////////// //
// /////////////// MAIN FUNCTION /////////////// //
// ///////////////////////////////////////////// //

int main (int argc, char **argv) {
           
    size_t num_sweeps = 20; 
    size_t num_adaptive_steps = 0;
    size_t num_iter_infec_time = 1;
    unsigned seed = DEFAULT_SEED; // FIXME : this should be a user optional variable  

    //========== USER parameters ==========
    if (argc < 8){
        // print usage here and return, since that's what you need.
        cout << "Wrong usage" << endl;
        std::exit(EXIT_FAILURE);
    }
    
    // Arguments num_sweepssed to main
    std::string tmp = argv[1];  // If we estimate contact pattern parameters
    // --- Temporary --- FIXME: bool from user 
    bool infer_contact = false;
    if ((tmp == "yes") || (tmp == "true")) {
        infer_contact = true;
    }
    // --- End-of-Temporary ---
    size_t num_iter = std::stoi(argv[2]);  // Number of iterations in the MCMC
    std::string name_output = argv[3];  // Name of the output file
    std::string name_data = argv[4];  // Name of the date file
    std::string path_data = argv[5];  // Path data file
    std::string path_to_save = argv[6];  // Path otput file

    tmp = argv[7];
    // --- Temporary --- FIXME: bool from user 
    bool display = false; 
    if ((tmp == "yes") || (tmp == "true")) {
        display = true;
    }
    // --- End-of-Temporary ---

    //==========Output files==========
    // FIXME : check for space and path validity

    //Paths
    std::string condition;
    std::string data_file, output_file, communityfile;

    //Data file
    // File structure :         0: indid, 1: hhid, 2: hhsize, 3: dds, 4: case, 5: study_period, 6: adult
    data_file = path_data + name_data + ".txt";         //Input file 
    output_file = path_to_save + "mcmc_" + name_output + ".txt"; //Output file
    communityfile = path_data + "community_infection_israel.txt";// communityfile

    cout << "Input file: " << data_file << endl;
    cout << "Output file: " << output_file << "\n\n";

    //========== Model parameters ==========

    // FIXME load parameters from file 
    Parameter this_param;
    this_param.init("alpha", 0.001, true, "unif", {0.0, 1.0}, "exp", 0.7, {0.0, 1.0});
    Parameter_Dictionary.emplace(this_param.p_name, this_param);
    this_param.init("beta", 0.6, true, "unif", {0.0, 5.0}, "exp", 0.1,  {0.0, 5.0}); 
    Parameter_Dictionary.emplace(this_param.p_name, this_param);
    this_param.init("delta", 1.0, true, "unif", {-3.0, 3.0}, "norm", 0.2, {-3.0, 3.0});
    Parameter_Dictionary.emplace(this_param.p_name, this_param);
    this_param.init("rU", 1.0, true, "lognorm", {0.0, 1.0}, "exp", 0.3, {0.0, 5.0});
    Parameter_Dictionary.emplace(this_param.p_name, this_param);
    this_param.init("r1V", 1.0, true, "lognorm", {0.0, 1.0}, "exp", 1.0, {0.0, 5.0});
    Parameter_Dictionary.emplace(this_param.p_name, this_param);
    this_param.init("r1V90", 1.0, true, "lognorm", {0.0, 1.0}, "exp", 1.0, {0.0, 0.1});
    Parameter_Dictionary.emplace(this_param.p_name, this_param);
    this_param.init("r2V", 1.0, true, "lognorm", {0.0, 1.0}, "exp", 1.0, {0.0, 0.1});
    Parameter_Dictionary.emplace(this_param.p_name, this_param);
    this_param.init("r3V", 1.0, true, "lognorm", {0.0, 1.0}, "exp", 0.9, {0.0, 5.0});
    Parameter_Dictionary.emplace(this_param.p_name, this_param);
    this_param.init("rC", 1.0, true, "lognorm", {0.0, 1.0}, "exp", 0.1, {0.0, 5.0});
    Parameter_Dictionary.emplace(this_param.p_name, this_param);
    this_param.init("rQ", 1.0, true, "lognorm", {0.0, 1.0}, "exp", 0.2, {0.0, 5.0});
    Parameter_Dictionary.emplace(this_param.p_name, this_param);
    this_param.init("rCC", 1.0, false, "lognorm", {0.0, 1.0}, "exp", 0.5, {0.0, 5.0});
    Parameter_Dictionary.emplace(this_param.p_name, this_param);
    this_param.init("rMC", 1.0, infer_contact, "lognorm", {0.0, 1.0}, "exp", 0.5, {0.0, 5.0});
    Parameter_Dictionary.emplace(this_param.p_name, this_param);
    this_param.init("rFC", 1.0, infer_contact, "lognorm", {0.0, 1.0}, "exp", 0.5, {0.0, 5.0});
    Parameter_Dictionary.emplace(this_param.p_name, this_param);
    this_param.init("rR", 1.0, true, "lognorm", {0.0, 1.0}, "exp", 1.0, {0.0, 5.0});
    Parameter_Dictionary.emplace(this_param.p_name, this_param);

    this_param.init("tost_tau", 5.783, false);
    Parameter_Dictionary.emplace(this_param.p_name, this_param);
    this_param.init("tost_mu", -4.0, false);
    Parameter_Dictionary.emplace(this_param.p_name, this_param);
    this_param.init("tost_sigma", 1.85, false);
    Parameter_Dictionary.emplace(this_param.p_name, this_param);
    this_param.init("tost_alpha", 5.85, false);
    Parameter_Dictionary.emplace(this_param.p_name, this_param);

    this_param.init("incub_shape", 3.44898, false);
    Parameter_Dictionary.emplace(this_param.p_name, this_param);
    this_param.init("incub_scale", 1.13077, false);
    Parameter_Dictionary.emplace(this_param.p_name, this_param);
    this_param.init("typical_size", 5.0, false);
    Parameter_Dictionary.emplace(this_param.p_name, this_param);

    //==========Build data==========
    // Load data
    std::vector<Household> hhData = build_data(data_file, communityfile, infer_contact);

    // Initialize MCMC object
    McmcObject mcmc(
        seed, 
        num_iter, 
        num_adaptive_steps,
        hhData,
        num_iter_infec_time
        );

    //==========MCMC==========
    mcmc.run(
        output_file,
        num_sweeps,
        display
        );

    return 0;

    std::exit(0);
}
