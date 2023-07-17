#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <string>
#include <random>
#include <chrono>
#include <vector>
#include <cmath>
#include <time.h>
#include <stdlib.h>
#include "household.hpp"
#include "mcmc_object.hpp"

using namespace std;

//----------------------------------------------------------------------
// Load data
// datafile with columns :  indid (int) currhhid (int) hhsize (int) onsetTime (double) study_period (double) is_case (0 or 1) vacc_status (0 to 3)  age (int) isolation (0 or 1) sex (0 or 1) vacc_delay (int) recovered (0 or 1)
//----------------------------------------------------------------------
std::vector<Household> build_data(std::string data_file, std::string community_infec_file, bool contact_pattern) {
    Individual this_indiv; // FIXME : why this is not std::Individual ?
    std::vector<double> instinfcom;
    std::vector<double> cuminfcom;

    instinfcom.resize(0);
    cuminfcom.resize(0);

    std::ifstream infileinf(community_infec_file.c_str());
    if(infileinf)
    {
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
    if(infile) {
        cout << "=============DATA===============" << endl;
        cout << "Read data file: " << data_file << endl;

        int indid(0), lasthhid(-1), currhhid(0), hhsize(0), is_case(0), start_follow_up(0), study_period(0), age(0),sex(0),vacc_status(0),isolation(0),vacc_delay(0),recovered(0);
        double onsetTime = 0.0;
        double ddi = 0.0;
        int num_cases =0, num_households = 0, num_subjects = 0, num_days =0;
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
        cout << "Number of cases: " << num_cases<< endl;
        cout << "Number of days: " << num_days << "\n\n";

        return(output);

    }else{
        cout << "ERROR: Cannot open the file." << endl;
        return(output);
    }
}


//----------------------------------------------------------------------
// Run mcmc
//----------------------------------------------------------------------
void runMCMC(
    McmcObject mcmc,
    std::string output_file,
    std::string augm_data_file,
    int pas,
    std::vector<int> id_selec_param,
    std::vector<std::string> paramNames,
    bool display
    ) {

	ofstream output(output_file.c_str());

	int iteration, iter, parID;
	int numberOfIteration = int(mcmc.iteration() / pas);
    
	int n_iter_time_infec = mcmc.get_num_iter_time_infec();

	// Column names of the output file
	std::string colNames="iteration logLik ";
	for (size_t p = 0; p < paramNames.size(); p++){
        colNames += paramNames[p] + " " + paramNames[p] + "_p " + paramNames[p] + "_a ";
    }
	colNames += "data_p data_a";
	output << colNames << endl;

	// MCMC chain
	cout << "=============MCMC===============" << endl;

    // Initial state
    mcmc.init_param_values();
	mcmc.init_augm_data(); // Initialize infection time and symptom onset for all infected individuals
	mcmc.init_log_lik(); 
	cout << "Initial log likelihood: " << mcmc.globalLogLik() << endl;

    output << "0 " << mcmc.globalLogLik() << " "; // Log likelihood
    for (size_t i = 0; i < mcmc.num_params(); i++)
        output << mcmc.parameter(i) << " 0 0 ";
    output << "0 0" << endl;

    // Chain
	for (iteration = 0; iteration < numberOfIteration; iteration++)	{
        cout << "iteration " << iteration << endl;
	    mcmc.reset_moves();

		for (iter = 0; iter < pas; iter++) {
			for (size_t selec_param = 0; selec_param < id_selec_param.size(); selec_param++) {
				parID = id_selec_param[selec_param];
				mcmc.update_parameter(parID, iteration*pas+iter+1, display);
			}
            //cout << "numberOfIteration " << iter << endl;

			// Augmented data
			for (int i = 0; i < n_iter_time_infec; i++) {
                mcmc.update_augm_infec_times(display); // Update loglik at the household level
                mcmc.update_augm_symptom_onset(display);
			}
            //cout << iteration * pas + iter << " " << mcmc.globalLogLik() << endl;
		}

		// Write log likelihood, parameter values, number of proposed/accepted move per parameter in the output file
		output << (iteration + 1) * pas << " " << mcmc.globalLogLik() << " "; // Log likelihood
		for (size_t i = 0; i < mcmc.num_params(); i++){
            output << mcmc.parameter(i) << " " << mcmc.proposed_move(i) << " " << mcmc.acceptedMove(i) << " ";
        }
        output << mcmc.proposed_move_data() << " " << mcmc.accepted_move_data() << endl;
	}

	output.close();
    //augmentedData.close();
}


// ///////////////////////////////////////////// //
// /////////////// MAIN FUNCTION /////////////// //
// ///////////////////////////////////////////// //

int main (int argc, char **argv) {

    if (argc < 8){
        // print usage here and return, since that's what you need.
    }
    
    // Arguments passed to main
    std::string tmp = argv[1];  // If we estimate contact pattern parameters
    bool inference_contact = false;
    if ((tmp == "yes") || (tmp == "true")) {
        inference_contact = true;
    }
    int numberOfIteration = std::stoi(argv[2]);  // Number of iterations in the MCMC
    std::string name_output = argv[3];  // Name of the output file
    std::string name_data = argv[4];  // Name of the date file
    std::string path_data = argv[5];  // Path data file
    std::string path_to_save = argv[6];  // Path otput file

    tmp = argv[7];
    bool display = false;
    if ((tmp == "yes") || (tmp == "true")) {
        display = true;
    }

    double max_PCR_detectability = 10.0;     
    std::string chainID = "2";

    //==========Model parameters==========
    // Initial values
    std::vector<std::string> parameterNames = {"alpha", "beta", "delta", "rU",  "r1V", "r1V90", "r2V", "r3V", "rC", "rQ", "rCC","rMC","rFC","rR"};
    //                                            0        1      2        3      4      5       6     7      8     9    10    11    12    13
    int num_params = parameterNames.size();
    //Use a lognorm prior or not for parameters (0 or 1) (useless for alpha, beta, delta)
    std::vector<int> lnormPrior(num_params);
    lnormPrior = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};

    //Standard deviation of the lognormal prior
    std::vector<double> sdLNormPrior(num_params);
    sdLNormPrior = {1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0};
    
    std::vector<double> parameter(num_params);

    // Parameters to infer according to model
    std::vector<int> selec_param(num_params, 0);
    
    //alpha
    parameter[0] = 0.001;
    selec_param[0] = 1;
    
    //beta
    parameter[1] = 0.6;
    selec_param[1] = 1;
    
    //delta 
    parameter[2] = 1.0;
    selec_param[2] = 1;

    //rU
    parameter[3] = 1;
    selec_param[3] = 1;
    
    //r1V
    parameter[4] = 1;
    selec_param[4] = 1;
    
    //r1V90
    parameter[5] = 1;
    selec_param[5] = 1;
    
    //r2V
    parameter[6] = 1;
    selec_param[6] = 1;

    //r3V
    parameter[7] = 1;
    selec_param[7] = 1;

    //rC
    parameter[8] = 1;
    selec_param[8] = 1;

    //rQ
    parameter[9] = 1;
    selec_param[9] = 1;

    //rR
    parameter[13] = 1;
    selec_param[13] = 1;

    if (inference_contact == true) {
        //CC
        parameter[10] = 1.0;
        selec_param[10] = 0;

        //rMC
        parameter[11] = 1;
        selec_param[11] = 1;

        //rFC
        parameter[12] = 1.0;
        selec_param[12] = 1;

    }else{
        //CC
        parameter[10] = 1.0;
        selec_param[10] = 0;

        //rMC
        parameter[11] = 1;
        selec_param[11] = 0;

        //rFC
        parameter[12] = 1.0;
        selec_param[12] = 0;

    }

	std::vector<int> id_selec_param(0);
    for(size_t param_index = 0; param_index < num_params; param_index++) {
        if(selec_param[param_index] == 1){
            id_selec_param.push_back(param_index);
        }
    }
    
    //Print id of selected parameters
    cout << "ID of selected parameters: "; 
    for (auto i = id_selec_param.begin(); i != id_selec_param.end(); ++i)
    	std::cout << parameterNames[*i] << ' ';
    cout << endl;

    //==========MCMC parameters==========
    /*    typedef std::chrono::high_resolution_clock myclock;
    myclock::time_point beginning = myclock::now();
    myclock::duration d = myclock::now() - beginning;
    unsigned seed = d.count();*/
    unsigned seed;
    if ( chainID == "1" ) seed = 20201208;
    if ( chainID == "2" ) seed = 20210415;
    if ( chainID == "3" ) seed = 20210601;
    cout << "Seed: " << seed << endl;

    int pas = 20; 
    //int nAdaptiveSteps = pas * numberOfIteration / 4;   // Doesn't work much...
    int nAdaptiveSteps = 0;
    int numberOfIterationTimeInfection = 1;

    // Variance of random walk
    std::vector<double> rateForRandomWalk(num_params);
    rateForRandomWalk = { 0.7,    0.1,     0.2,   0.3,  1.0,      1.0,    1.0,   0.9,  0.1,  0.2,   0.5,  0.5,  0.5,  1.0};
    //                  "alpha", "beta", "delta", "rU",  "r1V", "r1V90", "r2V", "r3V", "rC", "rQ", "rCC","rMC","rFC","rR"
    

    //==========Output files==========
    //Paths
    std::string condition;
    std::string data_file, output_file, augm_data_file,communityfile;

    //Data file
    // File structure :         0: indid, 1: hhid, 2: hhsize, 3: dds, 4: case, 5: study_period, 6: adult
    data_file = path_data+name_data+".txt";         //Input file 
    output_file = path_to_save + "mcmc_" + name_output + ".txt"; //Output file
    augm_data_file = path_to_save + name_output + ".txt"; //Output file
    communityfile = path_data+"community_infection_israel.txt";// communityfile

    cout << "Input file: " << data_file << endl;
    cout << "Output file: " << output_file << "\n\n";

    //==========Build data==========
    // Load data
    std::vector<Household> hhData = build_data(data_file, communityfile, inference_contact);

    // Initialize MCMC object
    McmcObject mcmc(
        seed, 
        numberOfIteration, 
        nAdaptiveSteps,
        hhData, 
        parameter, 
        selec_param, 
        rateForRandomWalk,
        lnormPrior,
        sdLNormPrior,
        numberOfIterationTimeInfection, 
        max_PCR_detectability          
        );

    //==========MCMC==========
    runMCMC(
        mcmc,
        output_file,
        augm_data_file,
        pas,
        id_selec_param,
        parameterNames,
        display
        );

    return 0;
}
