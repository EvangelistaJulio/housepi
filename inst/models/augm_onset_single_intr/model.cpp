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

void Model::assign_group_rate_to_individual(Individual &this_indiv)
{
	//===== rate n.0 : Infectivity =====//
	if (this_indiv.get_data("vacc_status") == 0)
	{ // naive : Infectivity [REFERENCE]
		this_indiv.assign_group_rate("Infectivity", REFERENCE);
	}
	else
	{ // fully vaccinated : Infectivity
		this_indiv.assign_group_rate("Infectivity", "inf_vacc");
	}

	//===== rate n.1 : Susceptibility =====//
	if (this_indiv.get_data("vacc_status") == 0)
	{ // naive : Susceptibility [REFERENCE]
		this_indiv.assign_group_rate("Susceptibility", REFERENCE);
	}
	else
	{ // fully vaccinated : Infectivity
		this_indiv.assign_group_rate("Susceptibility", "sus_vacc");
	}
}

/////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////  Model : assign rate_per group to individual  //////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

/* Assign the contact patterns between a couple of individuals. */
void Model::assign_contact_rate_to_household(Household &this_house, int from, int to)
{
	Individual &from_indiv = this_house.indiv[from];
	Individual &to_indiv = this_house.indiv[to];

	// from `i` -> to `k`
	// children if younger than 12 y.o., adults otherwise
	bool i_is_child = (from_indiv.get_data("age") < get_param("child_age"));
	bool i_is_adult = (i_is_child == false);

	bool k_is_child = (to_indiv.get_data("age") < get_param("child_age"));
	bool k_is_adult = (k_is_child == false);

	if (i_is_child && k_is_child)
	{ // Child to Child
		this_house.assign_rate_per_contact(from, to, "cont_child_child");
	}
	else if (i_is_child && k_is_adult)
	{ // Child to Adult
		this_house.assign_rate_per_contact(from, to, "cont_child_adult");
	}
	else if (i_is_adult && k_is_child)
	{ // Adult to Child
		this_house.assign_rate_per_contact(from, to, "cont_adult_child");
	}
	else
	{ // Adult to Adult : Contact Pattern [REFERENCE]
		this_house.assign_rate_per_contact(from, to, REFERENCE);
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
	//===== Initialize Lookup Tables (infectivity profiles) =====//
	dictargs args = {{"alpha", get_param("tost_alpha")},
					 {"mu", get_param("tost_mu")},
					 {"sigma", get_param("tost_sigma")},
					 {"tau", get_param("tost_tau")}};
	infec_prof_ptr = new Infec_Prof_Ferretti20(args);

	// clear the Lookup Tables (in case of multiple initializations)
	LUT_inst_prof.clear();
	LUT_cumul_prof.clear();

	//===== Initialize household features =====//

	// std::vector with index of first infector
	LUT_hh_first_case = std::vector<int>(num_house, -1);

	// produce the empty matrix for the infectivity profiles look-up-tables
	for (int h = 0; h < num_house; h++)
	{
		int h_size = house[h].size();
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

	int h = this_house.get_index();
	Individual &this_indiv = this_house.indiv[i];

	double onset_i = this_indiv.get_augm("onset_date");
	double incub_i = this_indiv.get_augm("incub_period");
	double t_infec_i = onset_i - incub_i;

	// init new first case to i
	double new_t0 = t_infec_i;
	int hh_first_case = i;

	double onset_k, incub_k, t_infec_k;

	//----- Update Lookup Tables for infectious -----//
	for (int k : this_house.infectious_set)
	{
		if (k != i)
		{
			Individual &other_indiv = this_house.indiv[k];
			onset_k = other_indiv.get_augm("onset_date");
			incub_k = other_indiv.get_augm("incub_period");
			t_infec_k = onset_k - incub_k;

			if (t_infec_i < t_infec_k)
			{ //----- `i` : infector, `k` : infectee -----//
				LUT_inst_prof[h][i][k] = infec_prof_ptr->pdf(t_infec_k - onset_i, incub_i);
				LUT_cumul_prof[h][i][k] = infec_prof_ptr->cdf(t_infec_k - onset_i, incub_i);
				LUT_inst_prof[h][k][i] = 0.;
				LUT_cumul_prof[h][k][i] = 0.;
			}
			else if (t_infec_k < t_infec_i)
			{ //----- `k` : infector, `i` : infectee -----//
				LUT_inst_prof[h][i][k] = 0.;
				LUT_cumul_prof[h][i][k] = 0.;
				LUT_inst_prof[h][k][i] = infec_prof_ptr->pdf(t_infec_i - onset_k, incub_k);
				LUT_cumul_prof[h][k][i] = infec_prof_ptr->cdf(t_infec_i - onset_k, incub_k);

				// update first case
				if (t_infec_k < new_t0)
				{
					hh_first_case = k;
					new_t0 = t_infec_k;
				}
			}
			// FIXME : else is the case t_infec_i == t_infec_k
		}
	}

	//===============================//
	//===== Susceptible Members =====//
	//===============================//

	t_infec_k = this_house.indiv[hh_first_case].get_augm("onset_date") + get_param("study_period");
	for (int k : this_house.susceptible_set)
	{
		LUT_cumul_prof[h][i][k] = infec_prof_ptr->cdf(t_infec_k - onset_i, incub_i);
	}

	//==================================//
	//===== Store house first case =====//
	//==================================//
	LUT_hh_first_case[h] = hh_first_case;
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
	for (int k : this_house.infectious_set)
	{
		if (k != i)
		{
			output += get_contact_rate(this_house, k, i) * get_group_rate(this_house.indiv[k], "Infectivity") * profiles[h][k][i];
		}
	}

	output *= get_group_rate(this_house.indiv[i], "Susceptibility") * transmission_rate(this_house);

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
	for (int i : this_house.susceptible_set)
	{
		loglike_house -= foi(this_house, i, LUT_cumul_prof);
	}

	// Contribution of symptomatic members
	for (int i : this_house.infectious_set)
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
	for (int i = 0; i < this_house.size(); i++)
	{
		if (!this_house.is_infectious(i))
		{
			assign_augm(this_house.indiv[i]);
		}
		this_house.change_status(i, SUSCEPTIBLE);
	}

	// pick a random introduction and change it's status
	int hh_first_case = uniform() * this_house.size();
	this_house.change_status(hh_first_case, INFECTIOUS);

	//=============================================================//
	//===== GENERATION >=1 : infection by previous generation =====//
	//=============================================================//

	double t = this_house.indiv[hh_first_case].get_augm("onset_date");
	double tf = t + get_param("study_period");

	// onset of secondary cases following the first case
	while ((t < tf) && (this_house.num_susceptible > 0))
	{
		t += Dt;
		for (int i = 0; i < this_house.size(); i++)
		{
			if (this_house.is_susceptible(i))
			{
				// set new status for computations
				this_house.change_status(i, INFECTIOUS);
				this_house.indiv[i].set_augm("onset_date", t);
				update_augm_dependecies(this_house, i);

				// restore status id escaped infection
				if (uniform() < exp(-foi(this_house, i, LUT_inst_prof) * Dt))
				{
					this_house.change_status(i, SUSCEPTIBLE);
				}
			}
		}
	}

	// remove augmented of susceptibles
	for (int i : this_house.susceptible_set)
	{
		clear_augm(this_house.indiv[i]);
	}
}