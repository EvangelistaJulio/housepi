// model.cpp //
// Copyright (C) 2024 Francesco Camaglia, Institut Pasteur //

///////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////  * *** * REED-FROST SINGLE INTRODUCTION * *** *  /////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

#include "model.h"

//------------------------------------------------------------------------------------------------------------//
//---------------------------------------- - - - CHANGE-ME - - - ---------------------------------------------//
//------------------------------------------------- || -------------------------------------------------------//
//------------------------------------------------- || -------------------------------------------------------//
//------------------------------------------------- || -------------------------------------------------------//
//------------------------------------------------- \/ -------------------------------------------------------//

/* Hazard of infection (-logarithm of the person-to-person escape probability). */
double Model::hazard_of_infection()
{
	return get_value("beta");
}

//------------------------------------------------- /\ -------------------------------------------------------//
//------------------------------------------------- || -------------------------------------------------------//
//------------------------------------------------- || -------------------------------------------------------//
//------------------------------------------------- || -------------------------------------------------------//
//---------------------------------------- - - - CHANGE-ME - - - ---------------------------------------------//
//------------------------------------------------------------------------------------------------------------//

//------------------------------------------------------------------------------------------------------------//
//------------------------------------------------------------------------------------------------------------//
//------------------------------------------------------------------------------------------------------------//
//------------------------------------------------------------------------------------------------------------//
//------------------------------------------------------------------------------------------------------------//
//------------------------------------------------------------------------------------------------------------//

/////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////  WARNINGS  /////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

bool Single_Introduction_Warning = false;
bool Single_House_Member_Warning = false;

//------------------------------------------------------------------------------------------------------------//
//------------------------------------------------------------------------------------------------------------//
//------------------------------------------------------------------------------------------------------------//
//------------------------------------------------------------------------------------------------------------//
//------------------------------------------------------------------------------------------------------------//
//------------------------------------------------------------------------------------------------------------//

/* Initialize the Lookup Table. */
void Model::initialize_model()
{
	// Unknown status are not modeled
	if (report_data["n.unk.stat"] > 0)
	{
		screen::raise_error("A missing status was found : incompatible with the model.");
	}

	// households without case don't have a contribution
	if (report_data["no.HH.case"] == 0)
	{
		screen::raise_warning(
			"Hoseholds without cases are not modeled and will be ignored. Discard them to improve computation time.",
			Single_Introduction_Warning);
	}

	// households of size 1 don't have a contribution
	if (report_data["min.HH.size"] == 1)
	{
		screen::raise_warning(
			"Hoseholds with one single member are not modeled and will be ignored. Discard them to improve computation time.", 
			Single_House_Member_Warning);
	}

	//=====Initialize look-up-table to largest possible value=====//
	bax::init_generations_lookuptable(report_data["max.n.HH.cases"]);
}

// -----------------------------------------------------------------------------------------------------------//
// -----------------------------------------------------------------------------------------------------------//
// -----------------------------------------------------------------------------------------------------------//
// -----------------------------------------------------------------------------------------------------------//
// -----------------------------------------------------------------------------------------------------------//
// -----------------------------------------------------------------------------------------------------------//

/////////////////////////////////////////////////////////////////////////////////////////////////
// ///////////////////////////// Model : Methods (likelihood)  /////////////////////////////// //
/////////////////////////////////////////////////////////////////////////////////////////////////

/* Compute the log-likelihood within household. */
double Model::compute_single_house_loglike(int h)
{
	double log_likelihood = 0;

	int hh_size = househ[h].house_size;
	int n_positive = househ[h].infectious_list.size();

	double log_p2p_esc_proba = -hazard_of_infection();
	double log_p2p_infec_proba = log(1. - exp(log_p2p_esc_proba));

	//==============================================//
	//========Contribution of Sero-Negatives========//
	//==============================================//

	log_likelihood += (hh_size - n_positive) * n_positive * log_p2p_esc_proba;

	//==============================================//
	//========Contribution of Sero-Positives========//
	//==============================================//

	if (n_positive > 1) // single introduction
	{
		// sum over all possible generations
		double  sero_positive_term = 0., tmp = 0.;
		for (unsigned c = 0; c < bax::Num_Gener_LookupTable[n_positive-1]; c++)
		{
			// reference to the LUT
			bax::gentype &gen_counts = bax::Generations_LookupTable[n_positive-1][c];

			// GENERATION (0) : infected by case zero
			tmp = gen_counts[0] * log_p2p_infec_proba;
			// GENERATION (0) : escaped by case zero
			tmp += (n_positive - 1 - gen_counts[0]) * log_p2p_esc_proba;

			// GENERATION (i) : within house, i > 0
			for (size_t i = 1; i < gen_counts.size(); i++)
			{
				// GENERATION (i) : infected by previous generation i - 1
				tmp += gen_counts[i] * log(1. - exp(gen_counts[i - 1] * log_p2p_esc_proba));

				// GENERATION (i) : escaped from other generations g, i. e. 0 <= g < i - 1
				for (size_t j = 0; j < i - 1; j++)
				{
					tmp += gen_counts[i] * gen_counts[j] * log_p2p_esc_proba;
				}
			}

			// times n. of equivalent configurations returning the same kind of generation
			sero_positive_term += bax::Multipl_Gener_LookupTable[n_positive-1][c] * exp(tmp);
		}

		//log_likelihood += log(sero_positive_term * bax::binomial_coefficient(hh_size - 1, n_positive-1));
		log_likelihood += log(sero_positive_term);
	}

	return log_likelihood;
}

// -----------------------------------------------------------------------------------------------------------//
// -----------------------------------------------------------------------------------------------------------//
// -----------------------------------------------------------------------------------------------------------//
// -----------------------------------------------------------------------------------------------------------//
// -----------------------------------------------------------------------------------------------------------//
// -----------------------------------------------------------------------------------------------------------//

/////////////////////////////////////////////////////////////////////////////////////////////////
// ////////////////////////////  Model : Methods (simulation)  /////////////////////////////// //
/////////////////////////////////////////////////////////////////////////////////////////////////

/* Simulate the infection within a household. */
void Model::simulate_single_house(int h)
{
	int num_new_infectious = 0;
	int num_old_infectious = 0;
	std::vector<int> store_old_list = {};
	int hh_size = househ[h].house_size;

	double log_p2p_esc_proba = -hazard_of_infection();

	//====================//
	//=====INITIALIZE=====//
	//====================//

	// set all cases to healthy (change infectious and unknown to susceptibles)
	store_old_list = househ[h].infectious_list;
	for (int i : store_old_list)
	{
		househ[h].change_status(i, SUSCEPTIBLE);
	}

	//==============================================//
	//====GENERATION 0 : by single introduction=====//
	//==============================================//

	// pick a random introduction
	househ[h].change_status(static_cast<int>(uniform() * hh_size), INFECTIOUS);
	num_new_infectious = 1;

	//==========================================================//
	//====GENERATION >=1 : infection by previous generation=====//
	//==========================================================//

	while ((househ[h].susceptibles_list.size() > 0) && (num_new_infectious > 0))
	{
		num_old_infectious = num_new_infectious;
		num_new_infectious = 0;

		store_old_list = househ[h].susceptibles_list;
		for (int i : store_old_list)
		{
			if (uniform() > exp(num_old_infectious * log_p2p_esc_proba))
			{
				househ[h].change_status(i, INFECTIOUS);
				num_new_infectious++;
			}
		}
	}
}
