/*
 * model.cpp
 *
 * Copyright (C) 2024 Francesco Camaglia, Institut Pasteur
 *
 */

///////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////  * *** * AUGMENTED ONSET - SINGLE INTRODUCTION * *** *  /////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

#include "model.h"

//-----------------------------------------------------------------------------------------------//
//---------------------------------- - - - CHANGE-ME - - - --------------------------------------//
//------------------------------------------- || ------------------------------------------------//
//------------------------------------------- || ------------------------------------------------//
//------------------------------------------- || ------------------------------------------------//
//------------------------------------------- \/ ------------------------------------------------//

///////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////  Model : assign rate per group to individual  ////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

void Model::assign_rate_per_group_to_individual(Individual &this_indiv)
{
	//===== rate n.0 : Infectivity =====//
	if (this_indiv.get_data("vacc_status") == 0)
	{ // naive : Infectivity [REFERENCE]
		this_indiv.assign_rate_per_group("Infectivity", REFERENCE);
	}
	else
	{ // fully vaccinated : Infectivity
		this_indiv.assign_rate_per_group("Infectivity", "i_vacc");
	}

	//===== rate n.1 : Susceptibility =====//
	if (this_indiv.get_data("vacc_status") == 0)
	{ // naive : Susceptibility [REFERENCE]
		this_indiv.assign_rate_per_group("Susceptibility", REFERENCE);
	}
	else
	{ // fully vaccinated : Infectivity
		this_indiv.assign_rate_per_group("Susceptibility", "s_vacc");
	}
}

/////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////  Model : assign rate_per group to individual  //////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

/* Assign the contact patterns between a couple of individuals. */
void Model::assign_rate_per_contact_to_household(int h, int i, int k)
{
	// from `i` -> to `k`
	// children if younger than 12 y.o., adults otherwise
	bool i_is_child = (get_data(h, i, "age") < get_param("child_age"));
	bool i_is_adult = (i_is_child == false);

	bool k_is_child = (get_data(h, k, "age") < get_param("child_age"));
	bool k_is_adult = (k_is_child == false);

	if (i_is_child && k_is_child)
	{ // Child to Child
		house[h].assign_rate_per_contact(i, k, "c_kid_to_kid");
	}
	else if (i_is_child && k_is_adult)
	{ // Child to Adult
		house[h].assign_rate_per_contact(i, k, "c_kid_to_adu");
	}
	else if (i_is_adult && k_is_child)
	{ // Adult to Child
		house[h].assign_rate_per_contact(i, k, "c_adu_to_kid");
	}
	else
	{ // Adult to Adult : Contact Pattern [REFERENCE]
		house[h].assign_rate_per_contact(i, k, REFERENCE);
	}
}

///////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////  Model : transmission rate  //////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/* Time of infection */
double Model::transmission_rate(int h) const
{
	return get_param("beta") * pow(4.0 / house[h].house_size, get_param("gamma"));
}

//------------------------------------------- /\ ------------------------------------------------//
//------------------------------------------- || ------------------------------------------------//
//------------------------------------------- || ------------------------------------------------//
//------------------------------------------- || ------------------------------------------------//
//---------------------------------- - - - CHANGE-ME - - - --------------------------------------//
//-----------------------------------------------------------------------------------------------//

//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//

///////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////  Model : common initialization  ////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/* Initialize parameters and augmented variables */
void Model::init_model_for_inference()
{
	//===== Initialize Lookup Tables (infectivity profiles) =====//

	dictargs args = {{"alpha", get_param("tost_alpha")},
					 {"mu", get_param("tost_mu")},
					 {"sigma", get_param("tost_sigma")},
					 {"tau", get_param("tost_tau")}};
	infec_prof_ptr = new Infec_Prof_Skew_Logistic(args);

	// clear the Lookup Tables (in case of multiple initializations)
	LUT_inst_profile.clear();
	LUT_cumul_profile.clear();

	//===== Initialize household features =====//

	// vector with index of first infector
	LUT_house_first_case = vector<int>(num_house, -1);

	// produce the empty matrix for the infectivity profiles look-up-tables
	for (int h = 0; h < num_house; h++)
	{
		int h_size = house[h].house_size;
		LUT_inst_profile.emplace_back(h_size, vector<double>(h_size, 0.));
		LUT_cumul_profile.emplace_back(h_size, vector<double>(h_size, 0.));
	}
}

//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//

///////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////// Model : store dependencies of augmented variables  ///////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/* Store cumulative and instantaneous infectivity profiles according to an infector.
   WARNING: unpredictable behavior if i is not an infector. */
void Model::augm_store_dependecies(int h, int i)
{
	//==============================//
	//===== Infectious Members =====//
	//==============================//

	double onset_i = get_augm(h, i, "onset_date");
	double incub_i = get_augm(h, i, "incub_period");
	double t_infec_i = onset_i - incub_i;

	// init new first case to i
	double new_t0 = t_infec_i;
	LUT_house_first_case[h] = i;

	double onset_k, incub_k, t_infec_k;

	//----- Update Lookup Tables for infectious -----//
	for (int k : house[h].infectious_set)
	{
		onset_k = get_augm(h, k, "onset_date");
		incub_k = get_augm(h, k, "incub_period");
		t_infec_k = onset_k - incub_k;

		if (t_infec_i < t_infec_k)
		{ //----- `i` : infector, `k` : infectee -----//
			LUT_inst_profile[h][i][k] = infec_prof_ptr->pdf(t_infec_k - onset_i, incub_i);
			LUT_cumul_profile[h][i][k] = infec_prof_ptr->cdf(t_infec_k - onset_i, incub_i);
			LUT_inst_profile[h][k][i] = 0.;
			LUT_cumul_profile[h][k][i] = 0.;
		}
		else if (t_infec_k < t_infec_i)
		{ //----- `k` : infector, `i` : infectee -----//
			LUT_inst_profile[h][i][k] = 0.;
			LUT_cumul_profile[h][i][k] = 0.;
			LUT_inst_profile[h][k][i] = infec_prof_ptr->pdf(t_infec_i - onset_k, incub_k);
			LUT_cumul_profile[h][k][i] = infec_prof_ptr->cdf(t_infec_i - onset_k, incub_k);

			// update first case
			if (t_infec_k < new_t0)
			{
				LUT_house_first_case[h] = k;
				new_t0 = t_infec_k;
			}
		}
		// NOTE : else is the case t_infec_i == t_infec_k, i. e. nothing to update
	}

	//===============================//
	//===== Susceptible Members =====//
	//===============================//

	t_infec_k = get_augm(h, LUT_house_first_case[h], "onset_date") + get_param("study_period");
	for (int k : house[h].susceptible_set)
	{
		LUT_cumul_profile[h][i][k] = infec_prof_ptr->cdf(t_infec_k - onset_i, incub_i);
	}
}

//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//

///////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////  Model : force of infection  ///////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/* Sum over all hazards of infection. */
double Model::force_of_infec(int h, int i, const vector<squarematrix> &profiles) const
{
	double output = 0.0;

	// Risk of infection from other household members (k)
	for (int k : house[h].infectious_set)
	{
		if (k != i)
		{
			output += get_contact_rate(h, k, i) * get_group_rate(h, k, "Infectivity") * profiles[h][k][i];
		}
	}

	output *= get_group_rate(h, i, "Susceptibility") * transmission_rate(h);

	return output;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////  Model : compute single house log-likelihood  /////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/* Compute the log-likelihood within household. */
double Model::compute_house_loglike(int h) const
{
	double loglike_house = 0.0;

	// Contribution of healthy members
	for (int i : house[h].susceptible_set)
	{
		loglike_house -= force_of_infec(h, i, LUT_cumul_profile);
	}

	// Contribution of symptomatic members
	for (int i : house[h].infectious_set)
	{
		// Contribution of secondary cases
		if (i != LUT_house_first_case[h])
		{
			loglike_house += log(force_of_infec(h, i, LUT_inst_profile));
			loglike_house -= force_of_infec(h, i, LUT_cumul_profile);
		}
	}

	// model of augmented variables
	// FIXME : should be taken care by Gibbs in case...
	loglike_house += compute_augm_loglike(h);

	return loglike_house;
}

//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////// Model : Methods (simulation)  ///////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/* Poisson point process to simulate outbreak in a single household. */
void Model::simulate_single_house(int h)
{
	// note: incub periods are initialized according to option in augm file
	// note: also simulation are conditioned on the the first cases

	double Dt = 0.01;

	//=================================================//
	//===== GENERATION 0 : by single introduction =====//
	//=================================================//

	// set all cases to susceptible
	for (int i = 0; i < house[h].house_size; i++)
	{
		change_status(h, i, SUSCEPTIBLE);
	}

	// pick a random introduction and leave user
	LUT_house_first_case[h] = uniform() * house[h].house_size;
	change_status(h, LUT_house_first_case[h], INFECTIOUS);

	//=============================================================//
	//===== GENERATION >=1 : infection by previous generation =====//
	//=============================================================//

	double t = get_augm(h, LUT_house_first_case[h], "onset_date");
	double tf = t + get_param("study_period");

	// onset of secondary cases following the first case
	while ((t < tf) && (house[h].num_susceptible > 0))
	{
		t += Dt;
		for (int i = 0; i < house[h].house_size; i++)
		{
			if (house[h].is_susceptible(i))
			{
				// set new status for computations
				change_status(h, i, INFECTIOUS);
				set_augm(h, i, "onset_date", t);
				augm_store_dependecies(h, i);

				// restore status id escaped infection
				if (uniform() < exp(-force_of_infec(h, i, LUT_inst_profile) * Dt))
				{
					change_status(h, i, SUSCEPTIBLE);
				}
			}
		}
	}

	// format onset times of susceptible
	for (int i : house[h].susceptible_set)
	{
		set_augm(h, i, "incub_period", 0.);
		set_augm(h, i, "onset_date", 1000.);
	}
}