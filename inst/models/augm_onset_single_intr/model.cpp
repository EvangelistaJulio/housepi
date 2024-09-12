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
		househ[h].assign_rate_per_contact(i, k, "c_kid_to_kid");
	}
	else if (i_is_child && k_is_adult)
	{ // Child to Adult
		househ[h].assign_rate_per_contact(i, k, "c_kid_to_adu");
	}
	else if (i_is_adult && k_is_child)
	{ // Adult to Child
		househ[h].assign_rate_per_contact(i, k, "c_adu_to_kid");
	}
	else
	{ // Adult to Adult : Contact Pattern [REFERENCE]
		househ[h].assign_rate_per_contact(i, k, REFERENCE);
	}
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
/////////////////////////////////////  GLOBAL VARIABLES  //////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

//----- Attributes for Look-Up-Tables -----//
std::vector<squarematrix> LUT_inst_profile = {};
std::vector<squarematrix> LUT_cumul_profile = {};

//----- Infectivity Profile Functions -----//
Infec_Prof_Skew_Logistic LUT_infec_prof_func;

///////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////  Model : common initialization  ////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/* Initialize parameters and augmented variables */
void Model::common_init_model()
{
	//===== Initialize Lookup Tables (infectivity profiles) =====//

	sdargs_ptr args = {{"tost_alpha", get_param_address("tost_alpha")},
					   {"tost_mu", get_param_address("tost_mu")},
					   {"tost_sigma", get_param_address("tost_sigma")},
					   {"tost_tau", get_param_address("tost_tau")}};
	LUT_infec_prof_func = Infec_Prof_Skew_Logistic(args);

	// clear the Lookup Tables (in case of multiple initializations)
	LUT_inst_profile.clear();
	LUT_cumul_profile.clear();

	//===== Initialize household features =====//

	// produce the empty matrix for the infectivity profiles look-up-tables
	for (int h = 0; h < num_househ; h++)
	{
		int h_size = househ[h].house_size;

		// produce the empty matrix for the infectivity profiles look-up-tables
		LUT_inst_profile.emplace_back(h_size, std::vector<double>(h_size, 0.));
		LUT_cumul_profile.emplace_back(h_size, std::vector<double>(h_size, 0.));
	}

	// assign study period (in days) to the households
	for (int h = 0; h < num_househ; h++)
	{
		househ[h].assign_data("study_period", get_param("study_period"));
	}
}

//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////  Model : get infection time  //////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/* Time of infection */
double Model::get_infec_time(int h, int i) const
{
	return get_augm(h, i, "onset_date") - get_augm(h, i, "incub_period");
}

///////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////// Model : instantaneous infectivity profile  ///////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/* Instanteneous infectivity profile. */
double Model::compute_infec_prof_inst(int h, int i, double t_transm)
{
	// tost : time from onset of symptoms to transmission
	double tost = t_transm - get_augm(h, i, "onset_date");
	return LUT_infec_prof_func.pdf(tost, get_augm(h, i, "incub_period"));
}

///////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////// Model : cumulative infectivity profile  ///////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/* Cumulative infectivity profile. */
double Model::compute_infec_prof_cumul(int h, int i, double t_transm)
{
	// tost : time from onset of symptoms to transmission
	double tost = t_transm - get_augm(h, i, "onset_date");
	return LUT_infec_prof_func.cdf(tost, get_augm(h, i, "incub_period"));
}

///////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////// Model : store dependencies of augmented variables  ///////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/* Store cumulative and instantaneous infectivity profiles according to an infector. */
void Model::augm_store_dependecies(int h, int i)
{
	//==============================//
	//===== Infectious Members =====//
	//==============================//

	double t_infec_i = get_infec_time(h, i);

	for (int k : househ[h].infectious_set)
	{
		double t_infec_k = get_infec_time(h, k);

		if (t_infec_i < t_infec_k)
		{ // `i` : infector, `k` : infectee
			LUT_inst_profile[h][i][k] = compute_infec_prof_inst(h, i, t_infec_k);
			LUT_cumul_profile[h][i][k] = compute_infec_prof_cumul(h, i, t_infec_k);
		}
		else if (t_infec_k < t_infec_i)
		{ // `k` : infector, `i` : infectee
			LUT_inst_profile[h][k][i] = compute_infec_prof_inst(h, k, t_infec_i);
			LUT_cumul_profile[h][k][i] = compute_infec_prof_cumul(h, k, t_infec_i);
		}
		// NOTE : else is the case t_infec_i == t_infec_k, i. e. nothing to update
	}

	//===============================//
	//===== Susceptible Members =====//
	//===============================//

	double tf = get_augm(h, househ[h].first_case, "onset_date") + get_data(h, "study_period");
	for (int k : househ[h].susceptible_set)
	{
		LUT_cumul_profile[h][i][k] = compute_infec_prof_cumul(h, i, tf);
	}
}

//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//

///////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////  Model : hazards of infection  ///////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/* Sum over all hazards of infection. */
double Model::hazards_of_infec(int h, int i, std::vector<squarematrix> &profiles)
{
	double output = 0.0;

	// Risk of infection from other household members (k)
	for (int k : househ[h].infectious_set)
	{
		if (k != i)
		{
			output += househ[h].get_contact_rate(k, i) * get_group_rate(h, i, "Infectivity") * profiles[h][k][i];
		}
	}

	return output * get_param("beta") * pow(4. / househ[h].house_size, get_param("gamma"));
}

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////  Model : compute single house log-likelihood  /////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/* Compute the log-likelihood within household. */
double Model::compute_house_loglike(int h)
{
	double loglike_house = 0.0;

	// Contribution of healthy members
	for (int i : househ[h].susceptible_set)
	{
		loglike_house -= get_group_rate(h, i, "Susceptibility") * hazards_of_infec(h, i, LUT_cumul_profile);
	}

	// Contribution of symptomatic members
	for (int i : househ[h].infectious_set)
	{
		// Contribution of model associated to augmented variables
		for (std::string augmID : augm_with_model_set)
		{
			loglike_house += househ[h].indiv[i].augm_model_log_pdf(augmID);
		}

		// Contribution of secondary cases
		if (i != househ[h].first_case)
		{

			loglike_house += log(get_group_rate(h, i, "Susceptibility") * hazards_of_infec(h, i, LUT_inst_profile));
			loglike_house -= get_group_rate(h, i, "Susceptibility") * hazards_of_infec(h, i, LUT_cumul_profile);
		}
	}

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

	double t0, tf, dt = 0.01;
	double esc_rate;

	//=================================================//
	//===== GENERATION 0 : by single introduction =====//
	//=================================================//

	// set all cases to susceptible
	for (int i = 0; i < househ[h].house_size; i++)
	{
		househ[h].change_status(i, SUSCEPTIBLE);
	}

	// pick a random introduction
	househ[h].first_case = uniform() * househ[h].house_size;
	househ[h].change_status(househ[h].first_case, INFECTIOUS);

	t0 = get_infec_time(h, househ[h].first_case);
	tf = get_augm(h, househ[h].first_case, "onset_date") + get_data(h, "study_period");

	//=============================================================//
	//===== GENERATION >=1 : infection by previous generation =====//
	//=============================================================//

	for (double t = dt + t0; (t < tf) && (househ[h].num_susceptible > 0); t += dt)
	{
		for (int i = 0; i < househ[h].house_size; i++)
		{
			if (househ[h].is_susceptible[i])
			{
				// set new status for computations
				househ[h].change_status(i, INFECTIOUS);
				set_augm(h, i, "onset_date", t + uniform() * dt + get_augm(h, i, "incub_period"));
				augm_store_dependecies(h, i);

				// refuse according to Bernoulli process
				esc_rate = get_group_rate(h, i, "Susceptibility") * hazards_of_infec(h, i, LUT_inst_profile) * dt;

				if (uniform() < exp(-esc_rate))
				{
					househ[h].change_status(i, SUSCEPTIBLE);
					augm_store_dependecies(h, i);
				}
			}
		}
	}

	// format onset times of susceptible
	for (int i : househ[h].susceptible_set)
	{
		set_augm(h, i, "incub_period", 0.);
		set_augm(h, i, "onset_date", 1000.);
	}
}