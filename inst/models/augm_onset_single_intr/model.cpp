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

//----- Infectivity Profile Functions -----//
Infec_Prof_Ferretti20 *infec_prof_ptr;

//----- Attributes for Look-Up-Tables -----//
std::vector<squarematrix> LUT_inst_prof;
std::vector<squarematrix> LUT_cumul_prof;
std::vector<int> LUT_hh_first_case;

//-----------------------------------------------------------------------------------------------//
//---------------------------------- - - - CHANGE-ME - - - --------------------------------------//
//------------------------------------------- || ------------------------------------------------//
//------------------------------------------- || ------------------------------------------------//
//------------------------------------------- || ------------------------------------------------//
//------------------------------------------- \/ ------------------------------------------------//

///////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////  Model : assign rate per group to individual  ////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

void Model::assign_group_rate_to_individual(Household &this_house, int i)
{
	//===== rate n.0 : Infectivity =====//
	if (this_house.get_data(i, "vacc_status") == 0)
	{ // naive : Infectivity [REFERENCE]
		this_house.assign_group_rate(i, "Infectivity", REFERENCE);
	}
	else
	{ // fully vaccinated : Infectivity
		this_house.assign_group_rate(i, "Infectivity", "inf_vacc");
	}

	//===== rate n.1 : Susceptibility =====//
	if (this_house.get_data(i, "vacc_status") == 0)
	{ // naive : Susceptibility [REFERENCE]
		this_house.assign_group_rate(i, "Susceptibility", REFERENCE);
	}
	else
	{ // fully vaccinated : Infectivity
		this_house.assign_group_rate(i, "Susceptibility", "sus_vacc");
	}
}

/////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////  Model : assign rate_per group to individual  //////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

/* Assign the contact patterns between a couple of individuals. */
void Model::assign_contact_rate_to_household(Household &this_house, int from, int to)
{
	// from `i` -> to `k`
	// children if younger than 12 y.o., adults otherwise
	bool i_is_child = (this_house.get_data(from, "age") < get_param("child_age"));
	bool i_is_adult = (i_is_child == false);

	bool k_is_child = (this_house.get_data(to, "age") < get_param("child_age"));
	bool k_is_adult = (k_is_child == false);

	if (i_is_child && k_is_child)
	{ // Child to Child
		this_house.assign_contact_rate(from, to, "cont_child_child");
	}
	else if (i_is_child && k_is_adult)
	{ // Child to Adult
		this_house.assign_contact_rate(from, to, "cont_child_adult");
	}
	else if (i_is_adult && k_is_child)
	{ // Adult to Child
		this_house.assign_contact_rate(from, to, "cont_adult_child");
	}
	else
	{ // Adult to Adult : Contact Pattern [REFERENCE]
		this_house.assign_contact_rate(from, to, REFERENCE);
	}
}

///////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////  Model : transmission rate  //////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/* Time of infection */
double Model::transmission_rate(const Household &this_house) const
{
	return get_param("beta") * pow(4.0 / this_house.size(), get_param("gamma"));
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
	infec_prof_ptr = new Infec_Prof_Ferretti20();

	// clear the Lookup Tables (in case of multiple initializations)
	LUT_inst_prof.clear();
	LUT_cumul_prof.clear();

	//===== Initialize household features =====//

	// indexes of first infector
	LUT_hh_first_case = std::vector<int>(get_num_house(), -1);

	// produce the empty matrix for the infectivity profiles look-up-tables
	for (size_t H = 0; H < get_num_house(); H++)
	{
		int h_size = house(H).size();
		LUT_inst_prof.emplace_back(h_size, std::vector<double>(h_size, 0.));
		LUT_cumul_prof.emplace_back(h_size, std::vector<double>(h_size, 0.));
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
void Model::update_augm_dependecies(Household &this_house, int i)
{
	//==============================//
	//===== Infectious Members =====//
	//==============================//

	int H = this_house.get_index();
	double t_onset_i = this_house.get_augm(i, "onset_date");
	double t_incub_i = this_house.get_augm(i, "incub_period");
	double t_infec_i = t_onset_i - t_incub_i;

	double t_onset_k, t_incub_k, t_infec_k, tost;

	// init new first case to i
	double new_t0 = t_infec_i;
	int hh_first_case = i;

	//----- Update Lookup Tables for infectious -----//
	for (int k : this_house.get_infectious_set())
	{
		if (k != i)
		{
			t_onset_k = this_house.get_augm(k, "onset_date");
			t_incub_k = this_house.get_augm(k, "incub_period");
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

	tost = this_house.get_augm(hh_first_case, "onset_date") + get_param("study_period") - t_onset_i;
	for (int k : this_house.get_susceptible_set())
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
double Model::foi(const Household &this_house, int i, const std::vector<squarematrix> &profiles) const
{
	int h = this_house.get_index();
	double output = 0.0;

	// Risk of infection from other household members (k)
	for (int k : this_house.get_infectious_set())
	{
		if (k != i)
		{
			output += this_house.get_contact_rate(k, i) * this_house.get_group_rate(k, "Infectivity") * profiles[h][k][i];
		}
	}

	output *= this_house.get_group_rate(i, "Susceptibility") * transmission_rate(this_house);

	return output;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////  Model : compute single house log-likelihood  /////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/* Compute the log-likelihood within household. */
double Model::compute_house_loglike(const Household &this_house) const
{
	int h = this_house.get_index();
	double loglike_house = 0.0;

	// Contribution of healthy members
	for (int i : this_house.get_susceptible_set())
	{
		loglike_house -= foi(this_house, i, LUT_cumul_prof);
	}

	// Contribution of symptomatic members
	for (int i : this_house.get_infectious_set())
	{
		// Contribution of secondary cases
		if (i != LUT_hh_first_case[h])
		{
			loglike_house += log(foi(this_house, i, LUT_inst_prof)) - foi(this_house, i, LUT_cumul_prof);
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
void Model::simulate_single_house(Household &this_house)
{
	double Dt = 0.01;

	//=================================================//
	//===== GENERATION 0 : by single introduction =====//
	//=================================================//

	//----- Preload incubation periods -----//
	// note: incub_period are initialized according to definitions
	// - set infectious to susceptible w/o dropping augm variables;
	// - add augmented variables to all others (also missing);
	for (int i = 0; i < this_house.size(); ++i)
	{
		if (this_house.get_status(i) != INFECTIOUS)
		{
			assign_all_augm(this_house, i);
		}
		this_house.change_status(i, SUSCEPTIBLE);
	}

	// pick a random introduction and change it's status
	int hh_first_case = uniform() * this_house.size();
	this_house.change_status(hh_first_case, INFECTIOUS);

	//=============================================================//
	//===== GENERATION >=1 : infection by previous generation =====//
	//=============================================================//

	double t = this_house.get_augm(hh_first_case, "onset_date");
	double tf = t + get_param("study_period");

	// onset of secondary cases following the first case
	while ((t < tf) && (this_house.get_num_susceptible() > 0))
	{
		t += Dt;
		for (int i = 0; i < this_house.size(); ++i)
		{
			if (this_house.get_status(i) == SUSCEPTIBLE)
			{
				// set new status for computations
				this_house.change_status(i, INFECTIOUS);
				this_house.set_augm_value(i, "onset_date", t);
				update_augm_dependecies(this_house, i);

				// restore status id escaped infection
				if (uniform() < exp(-foi(this_house, i, LUT_inst_prof) * Dt))
				{
					this_house.change_status(i, SUSCEPTIBLE);
				}
			}
		}
	}

	// format all augmented variables of susceptibles
	for (int i : this_house.get_susceptible_set())
	{
		clear_all_augm(this_house, i);
	}
}