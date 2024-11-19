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

/////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////  GLOBAL VARIABLES  ///////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

//----- Attributes for Look-Up-Tables -----//
std::vector<squarematrix> LUT_inst_prof;
std::vector<squarematrix> LUT_cumul_prof;
std::vector<int> LUT_hh_first_case;

//----- Infectivity Profile Functions -----//
std::shared_ptr<Infec_Prof_Ferretti20> infec_prof_ptr;

//-----------------------------------------------------------------------------------------------//
//---------------------------------- - - - CHANGE-ME - - - --------------------------------------//
//------------------------------------------- || ------------------------------------------------//
//------------------------------------------- || ------------------------------------------------//
//------------------------------------------- || ------------------------------------------------//
//------------------------------------------- \/ ------------------------------------------------//

///////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////  Model : assign rate per group to individual  ////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

void Model::assign_group_rate_to_individual(Household &house, int i)
{
	//===== rate n.0 : Infectivity =====//
	if (house.data(i, "vacc_status") == 0)
	{ // naive : Infectivity [REFERENCE]
		house.assign_group_rate(i, "Infectivity", REFERENCE);
	}
	else
	{ // fully vaccinated : Infectivity
		house.assign_group_rate(i, "Infectivity", "inf_vacc");
	}

	//===== rate n.1 : Susceptibility =====//
	if (house.data(i, "vacc_status") == 0)
	{ // naive : Susceptibility [REFERENCE]
		house.assign_group_rate(i, "Susceptibility", REFERENCE);
	}
	else
	{ // fully vaccinated : Infectivity
		house.assign_group_rate(i, "Susceptibility", "sus_vacc");
	}
}

/////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////  Model : assign rate_per group to individual  //////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

/* Assign the contact patterns between a couple of individuals. */
void Model::assign_contact_rate_to_household(Household &house, int from, int to)
{
	// child if younger than 'child_age'; adult otherwise
	bool is_from_child = (house.data(from, "age") < param("child_age"));
	bool is_to_child = (house.data(to, "age") < param("child_age"));

	if (is_from_child && is_to_child)
	{ // Child to Child
		house.assign_contact_rate(from, to, "cont_child_child");
	}
	else if (is_from_child && (!is_to_child))
	{ // Child to Adult
		house.assign_contact_rate(from, to, "cont_child_adult");
	}
	else if ((!is_from_child) && is_to_child)
	{ // Adult to Child
		house.assign_contact_rate(from, to, "cont_adult_child");
	}
	else
	{ // Adult to Adult : Contact Pattern [REFERENCE]
		house.assign_contact_rate(from, to, REFERENCE);
	}
}

///////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////  Model : transmission rate  //////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/* Time of infection */
double Model::transmission_rate(const Household &house) const
{
	return param("beta") * pow(param("size_ref") / house.size(), param("gamma"));
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
	// use infectivity profile from Ferretti et al (2020)
	infec_prof_ptr = std::make_shared<Infec_Prof_Ferretti20>();

	// clear the Lookup Tables (in case of multiple initializations)
	LUT_inst_prof.clear();
	LUT_cumul_prof.clear();

	//===== Initialize household features =====//

	// indexes of first infector
	LUT_hh_first_case = std::vector<int>(num_houses());

	// produce the empty matrix for the infectivity profiles look-up-tables
	for (size_t H = 0; H < num_houses(); ++H)
	{
		int h_size = house_size(H);
		LUT_inst_prof.emplace_back(h_size, std::vector<double>(h_size));
		LUT_cumul_prof.emplace_back(h_size, std::vector<double>(h_size));
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
void Model::update_augm_dependecies(Household &house, int i)
{
	//==============================//
	//===== Infectious Members =====//
	//==============================//

	int H = house.index();
	double t_onset_i = house.augm(i, "onset_date");
	double t_incub_i = house.augm(i, "incub_period");
	double t_infec_i = t_onset_i - t_incub_i;

	double t_onset_k, t_incub_k, t_infec_k, tost;

	// init new first case to i
	double new_t0 = t_infec_i;
	int hh_first_case = i;

	//----- Update Lookup Tables for infectious -----//
	for (int k : house.infectious_set())
	{
		if (k != i)
		{
			t_onset_k = house.augm(k, "onset_date");
			t_incub_k = house.augm(k, "incub_period");
			t_infec_k = t_onset_k - t_incub_k;

			if (t_infec_i < t_infec_k)
			{ //----- `i` : infector, `k` : infectee -----//
				tost = t_infec_k - t_onset_i;
				LUT_inst_prof[H][i][k] = infec_prof_ptr->pdf(tost, t_incub_i);
				LUT_cumul_prof[H][i][k] = infec_prof_ptr->cdf(tost, t_incub_i);
				LUT_inst_prof[H][k][i] = 0.;
				LUT_cumul_prof[H][k][i] = 0.;
			}
			else
			{ //----- `k` : infector, `i` : infectee -----//
				tost = t_infec_i - t_onset_k;
				LUT_inst_prof[H][i][k] = 0.;
				LUT_cumul_prof[H][i][k] = 0.;
				LUT_inst_prof[H][k][i] = infec_prof_ptr->pdf(tost, t_incub_k);
				LUT_cumul_prof[H][k][i] = infec_prof_ptr->cdf(tost, t_incub_k);

				//----- update first case -----//
				if (t_infec_k < new_t0)
				{
					hh_first_case = k;
					new_t0 = t_infec_k;
				}
			}
		}
	}

	//===============================//
	//===== Susceptible Members =====//
	//===============================//

	tost = house.augm(hh_first_case, "onset_date") + param("study_period") - t_onset_i;
	for (int k : house.susceptible_set())
	{
		LUT_cumul_prof[H][i][k] = infec_prof_ptr->cdf(tost, t_incub_i);
	}

	//==================================//
	//===== Store house first case =====//
	//==================================//
	LUT_hh_first_case[H] = hh_first_case;
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
double Model::foi(const Household &house, int i, const std::vector<squarematrix> &profiles) const
{
	int h = house.index();
	double output = 0.0;

	// Risk of infection from other household members (k)
	for (int k : house.infectious_set())
	{
		if (k != i)
		{
			output += house.contact_rate(k, i) * house.group_rate(k, "Infectivity") * profiles[h][k][i];
		}
	}

	output *= house.group_rate(i, "Susceptibility") * transmission_rate(house);

	return output;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////  Model : compute single house log-likelihood  /////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/* Compute the log-likelihood within household. */
double Model::compute_house_loglike(const Household &house) const
{
	int h = house.index();
	double loglike_house = 0.0;

	// Contribution of healthy members
	for (int i : house.susceptible_set())
	{
		loglike_house -= foi(house, i, LUT_cumul_prof);
	}

	// Contribution of symptomatic members
	for (int i : house.infectious_set())
	{
		// Contribution of secondary cases
		if (i != LUT_hh_first_case[h])
		{
			loglike_house += log(foi(house, i, LUT_inst_prof)) - foi(house, i, LUT_cumul_prof);
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
void Model::simulate_single_house(Household &house)
{
	double Dt = 0.01;

	//=================================================//
	//===== GENERATION 0 : by single introduction =====//
	//=================================================//

	//----- Preload incubation periods -----//
	// note: incub_period are initialized according to definitions
	// - set infectious to susceptible w/o dropping augm variables;
	// - add augmented variables to all others (also missing);
	for (int i = 0; i < house.size(); ++i)
	{
		if (house.status(i) != INFECTIOUS)
		{
			assign_all_augm(house, i);
		}
		house.change_status(i, SUSCEPTIBLE);
	}

	// pick a random introduction and change it's status
	int hh_first_case = uniform() * house.size();
	house.change_status(hh_first_case, INFECTIOUS);

	//=============================================================//
	//===== GENERATION >=1 : infection by previous generation =====//
	//=============================================================//

	double t = house.augm(hh_first_case, "onset_date");
	double tf = t + param("study_period");

	// onset of secondary cases following the first case
	while ((t < tf) && (house.num_susceptible() > 0))
	{
		t += Dt;
		for (int i = 0; i < house.size(); ++i)
		{
			if (house.status(i) == SUSCEPTIBLE)
			{
				// set new status for computations
				house.change_status(i, INFECTIOUS);
				house.set_augm_value(i, "onset_date", t);
				update_augm_dependecies(house, i);

				// restore status id escaped infection
				if (uniform() < exp(-foi(house, i, LUT_inst_prof) * Dt))
				{
					house.change_status(i, SUSCEPTIBLE);
				}
			}
		}
	}

	// format all augmented variables of susceptibles
	for (int i : house.susceptible_set())
	{
		clear_all_augm(house, i);
	}
}