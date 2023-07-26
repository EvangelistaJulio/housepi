#ifndef PARAMETER__H
#define PARAMETER__H

#include <string>
#include <vector>
#include <unordered_map>
#include <cmath>
#include <random>
#include <iostream>

#include "user_model.hpp"

typedef double (*prior_pdf_ptr)(double x, std::vector<double> args);
typedef double (*prior_rand_ptr)(std::mt19937_64& rand_gen, std::vector<double> args);

//----------------------------------------------------------------------
// Parameter : Class declaration
//----------------------------------------------------------------------

class Parameter {
    // Attributes
    double p_value = 0.;                          	// value of the parameter
    std::string p_prior = "_const_";    			// choice of the prior to the parameter 
    std::vector<double> p_prior_args = {};        	// arguments to the prior function
	std::string p_metropolis_strategy = "exp";    	// strategy for the Metropolis-Hastings
	double p_metropolis_rate = 0.5;          		// rate for the Metropolis-Hastings
    std::vector<double> p_range = {};           	// range of values
	prior_pdf_ptr p_prior_log_pdf = nullptr;		// pointer to the log prior pdf
	prior_rand_ptr p_prior_rand = nullptr;			// the actual generator to execute
		
	int _metropolis_num_accept = 0;				// mcmc internal // FIXME : maybe inheritance?
	int _metropolis_num_propos = 0;
	double _metropolis_accept = 0.0;

	public:
        std::string p_name = "_no_name_";			// parameter name          
        bool p_infer = false;                     	// wheter to infer the parameter or not

			
		// FIXME prior args to be initialized
		// FIXME range to be initialized
		// Constructors
		Parameter(
			std::string name,
			double value,
			bool infer,
			std::string prior = "_const_",
			std::vector<double> args = {},
			std::string metropolis_strategy = "exp",
			double metropolis_rate = 0.,
			std::vector<double> range = {}
			);
		Parameter() = default;
		~Parameter() {};
		void init(
			std::string name,
			double value,
			bool infer,
			std::string prior = "_const_",
			std::vector<double> args = {},
			std::string metropolis_strategy = "exp",
			double metropolis_rate = 0.,
			std::vector<double> range = {}
			);

		// Getter
		double get_value () const; 
		
		// Getter (mcmc internal)
		int get_num_accept () { return _metropolis_num_accept; } 
		int get_num_propos () { return _metropolis_num_propos; } 
		double get_accept () { return _metropolis_accept; } 

		// Setter
		void set_value (double value); 
		void set_rate_metropolis (double rate); 

		// Setter (mcmc internal)
		void reset_accept_propos ();
		void set_num_accept_pp () { _metropolis_num_accept++; } 
		void set_num_propos_pp () { _metropolis_num_propos++; } 	
		void set_accept (double accept) { _metropolis_accept = accept; } 

		// Computations
		bool within_range (double value);
        double generate_from_prior(std::mt19937_64& rand_gen);
        double metropolis (std::mt19937_64& rand_gen);
        double metropolis_log_ratio_proposal();
        double metropolis_log_ratio_prior(double new_value);
};

////////////////////////////////////////////////////////////
/////////////// SHARED PARAMETERS DICTIONARY ///////////////
////////////////////////////////////////////////////////////

// Shared Parameter Dictionary
extern std::unordered_map<std::string, Parameter> Parameter_Dictionary;
// aliases
Parameter _param (std::string parID);
double _param_value (std::string parID);
bool _param_infer (std::string parID);

#endif // PARAMETER__H