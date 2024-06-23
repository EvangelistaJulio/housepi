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

/* The logarithm of the person-to-person escape probability. */
double Model::log_person_to_person_escape_probability()
{
	return -get_value("beta");
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
bool Single_Member_House_Warning = false;

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
		raise_error_message("A missing status was found : incompatible with the model.");
	}

	// households without case don't have a contribution
	if (report_data["no.HH.case"] == 0)
	{
		raise_warning_message(
			"Hoseholds without cases are not modeled and will be ignored. Discard them to increase computation time.",
			Single_Introduction_Warning);
	}

	// households of size 1 don't have a contribution
	if (report_data["min.HH.size"] == 1)
	{
		raise_warning_message(
			"Hoseholds with one single member are not modeled and will be ignored. Discard them to increase computation time.", 
			Single_Member_House_Warning);
	}

	//=====Find maximum possible number of sero-positives=====//
	int max_num = 0, aux = 0;
	for (int h = 0; h < num_househ; h++)
	{
		aux = househ[h].infectious_list.size();
		max_num = (aux > max_num) ? aux : max_num;
	}

	//=====Initialize look-up-table to largest possible value=====//
	bax::init_generations_lookuptable(max_num);
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

	double log_p2p_esc_proba = log_person_to_person_escape_probability();
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
		int multiplicity;
		double tmp = 0., sero_positive_term = 0.;

		// sum over all possible generations
		for (bax::gentype gen_counts : bax::Generations_LookupTable.at(n_positive - 1))
		{
			// number of equivalent configurations returning the same kind of generation
			multiplicity = bax::factorial(n_positive) / bax::factorial(gen_counts[0]); // euclidean division

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

				multiplicity /= bax::factorial(gen_counts[i]);
			}

			sero_positive_term += multiplicity * exp(tmp);
		}

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

	double log_p2p_esc_proba = log_person_to_person_escape_probability();

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
