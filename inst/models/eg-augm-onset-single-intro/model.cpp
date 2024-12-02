/*
 * model.cpp
 *
 * Model template of the 'housepi' package.
 * Copyright (C) 2024 Francesco Camaglia, Institut Pasteur
 *
 */

///////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////  * *** * AUGMENTED ONSET - SINGLE INTRODUCTION * *** *  /////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

#include "model.h"

//-----------------------------------------------------------------------------------------------//
//---------------------------------------- CHANGE-ME --------------------------------------------//
//------------------------------------------- || ------------------------------------------------//
//------------------------------------------- || ------------------------------------------------//
//------------------------------------------- || ------------------------------------------------//
//------------------------------------------- \/ ------------------------------------------------//

/* [CHANGE-ME] Assign the parameters controlling infectivity for each individual. */
void Model::assign_group_rate_to_individual(Household &house, indiv i)
{
	std::string rateID;
	bool is_not_vacc = (house.data(i, "vacc_status") == 0);

	//===== Infectivity =====//
	rateID = "Infectivity";
	if (is_not_vacc)
		house.assign_group_rate(i, rateID, REFERENCE);
	else
		house.assign_group_rate(i, rateID, "inf_vacc");

	//===== Susceptibility =====//
	rateID = "Susceptibility";
	if (is_not_vacc)
		house.assign_group_rate(i, rateID, REFERENCE);
	else
		house.assign_group_rate(i, "Susceptibility", "sus_vacc");
}

/* [CHANGE-ME] Assign the contact patterns between a couple of individuals. */
void Model::assign_contact_rate_to_household(Household &house, indiv from, indiv to)
{
	// child if younger than 'child_age'; adult otherwise
	bool is_from_adult = (house.data(from, "age") < param("child_age"));
	bool is_to_adult = (house.data(to, "age") < param("child_age"));

	if (is_from_adult && is_to_adult)
		house.assign_contact_rate(from, to, REFERENCE);

	else if (is_from_adult && (!is_to_adult))
		house.assign_contact_rate(from, to, "cont_adult_child");

	else if ((!is_from_adult) && is_to_adult)
		house.assign_contact_rate(from, to, "cont_child_adult");

	else
		house.assign_contact_rate(from, to, "cont_child_child");
}

/* [CHANGE-ME] Time of infection */
double Model::transmission_rate(const Household &house) const
{
	return param("beta") * pow(param("size_ref") / house.size(), param("gamma"));
}

//------------------------------------------- /\ ------------------------------------------------//
//------------------------------------------- || ------------------------------------------------//
//------------------------------------------- || ------------------------------------------------//
//------------------------------------------- || ------------------------------------------------//
//---------------------------------------- CHANGE-ME --------------------------------------------//
//-----------------------------------------------------------------------------------------------//

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////  Model : initialization for simulation  ///////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/* Initialize parameters and augmented variables. */
void Model::init_model_for_simulation()
{
	init_model_for_inference();
}

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////  Model : initialization for inference  ////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/* Initialize parameters and augmented variables. */
void Model::init_model_for_inference()
{
	// use infectivity profile from Ferretti et al (2020)
	assign_infectivity_profile("ferretti20");

	// integer for index of introduction case
	assign_extra_integer_to_all_houses({"intro-case"}, -1);

	// matrices for the instantaneous and cumulative infectivity profiles
	assign_extra_matrix_to_all_houses({"inst", "cumul"}, 0.0);
}

///////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////  Model : followup period  ////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/* Time when followup ended. */
double Model::time_end_followup(const Household &house) const
{
	int i0 = house.extra_int("intro-case");
	return house.augm(i0, "onset_date") + param("study_period");
}

///////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////  Model : force of infection  ///////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/* Sum over all hazards of infection. */
double Model::_houseROI_(const Household &house, indiv i, const squarematrix &profiles) const
{
	double output = 0.0;

	// Risk of infection from other household members (k)
	for (indiv k : house.infectious_set())
	{
		if (k != i)
		{
			output += house.contact_rate(k, i) * house.group_rate(k, "Infectivity") * profiles[k][i];
		}
	}

	output *= house.group_rate(i, "Susceptibility") * transmission_rate(house);

	return output;
}

/* Instantaneous Risk of Infection for individual i. */
double Model::ROI_inst(const Household &house, indiv i) const
{
	return param("alpha") + _houseROI_(house, i, house.extra_matrix("inst"));
}

/* Cumulative Risk of Infection for infectious individual i. */
double Model::ROI_cumul_infec(const Household &house, indiv i) const
{
	indiv i0 = house.extra_int("intro-case");
	double t_intro = house.augm(i0, "onset_date") - house.augm(i0, "incub_period");
	double exposure = house.augm(i, "onset_date") - house.augm(i, "incub_period") - t_intro;

	return param("alpha") * exposure + _houseROI_(house, i, house.extra_matrix("cumul"));
}

/* Cumulative Risk of Infection for susceptible individual i. */
double Model::ROI_cumul_suscept(const Household &house, indiv i) const
{
	indiv i0 = house.extra_int("intro-case");
	double t_intro = house.augm(i0, "onset_date") - house.augm(i0, "incub_period");
	double exposure = time_end_followup(house) - t_intro;

	return param("alpha") * exposure + _houseROI_(house, i, house.extra_matrix("cumul"));
}

///////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////// Model : store dependencies of augmented variables  ///////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/* Store cumulative and instantaneous infectivity profiles according to an infector.
   WARNING: unpredictable behavior if i is not an infector.
   WARNING: this function should not call any mutable parameter. */
void Model::update_augm_dependecies(Household &house, indiv i)
{
	double serial;

	//==============================//
	//===== Infectious Members =====//
	//==============================//

	double onset_i = house.augm(i, "onset_date");
	double incub_i = house.augm(i, "incub_period");
	double t_infec_i = onset_i - incub_i;

	// init new first case to i
	double new_t0 = t_infec_i;
	house.set_extra_int("intro-case", i);

	//----- Update Lookup Tables for infectious -----//
	double onset_k, incub_k, t_infec_k;

	for (indiv k : house.infectious_set())
	{
		if (k != i)
		{
			onset_k = house.augm(k, "onset_date");
			incub_k = house.augm(k, "incub_period");
			t_infec_k = onset_k - incub_k;

			if (t_infec_i < t_infec_k)
			{ //----- `i` : infector, `k` : infectee -----//
				serial = t_infec_k - t_infec_i;

				house.set_extra_matrix("inst", i, k, infec_profile_inst(serial, {incub_i}));
				house.set_extra_matrix("inst", k, i, 0.0);

				house.set_extra_matrix("cumul", i, k, infec_profile_cumul(serial, {incub_i}));
				house.set_extra_matrix("cumul", k, i, 0.0);
			}
			else
			{ //----- `k` : infector, `i` : infectee -----//
				serial = t_infec_i - t_infec_k;

				house.set_extra_matrix("inst", i, k, 0.0);
				house.set_extra_matrix("inst", k, i, infec_profile_inst(serial, {incub_k}));

				house.set_extra_matrix("cumul", i, k, 0.0);
				house.set_extra_matrix("cumul", k, i, infec_profile_cumul(serial, {incub_k}));

				//----- update first case -----//
				if (t_infec_k < new_t0)
				{
					new_t0 = t_infec_k;
					house.set_extra_int("intro-case", k);
				}
			}
		}
	}

	//===============================//
	//===== Susceptible Members =====//
	//===============================//

	serial = time_end_followup(house) - t_infec_i;
	for (indiv k : house.susceptible_set())
	{
		house.set_extra_matrix("cumul", i, k, infec_profile_cumul(serial, {incub_i}));
		house.set_extra_matrix("cumul", k, i, 0.0);
	}
}

//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////  Model : compute single house log-likelihood  /////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/* Compute the log-likelihood within household. */
double Model::compute_house_loglike(Household &house) const
{
	double loglike_house = 0.0;

	//===============================//
	//===== Susceptible Members =====//
	//===============================//

	for (indiv i : house.susceptible_set())
	{
		loglike_house -= ROI_cumul_suscept(house, i);
	}

	//==============================//
	//===== Infectious Members =====//
	//==============================//
	
	for (indiv i : house.infectious_set())
	{
		// Contribution of secondary cases
		if (i != house.extra_int("intro-case"))
		{
			loglike_house += log(ROI_inst(house, i)) - ROI_cumul_infec(house, i);
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
///////////////////////////////// Model : simulate single house  //////////////////////////////////
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
	for (indiv i = 0; i < house.size(); ++i)
	{
		if (house.status(i) != INFECTIOUS)
		{
			assign_all_augm(house, i);
		}
		house.change_status(i, SUSCEPTIBLE);
	}

	// pick a random introduction and change it's status
	indiv i0 = uniform() * house.size();
	house.change_status(i0, INFECTIOUS);
	house.set_extra_int("intro-case", i0);

	//=============================================================//
	//===== GENERATION >=1 : infection by previous generation =====//
	//=============================================================//

	double t = house.augm(i0, "onset_date");
	double tf = time_end_followup(house);

	// onset of secondary cases following the first case
	while ((t < tf) && (house.n_susceptible() > 0))
	{
		t += Dt;
		for (indiv i = 0; i < house.size(); ++i)
		{
			if (house.status(i) == SUSCEPTIBLE)
			{
				// set new status for computations
				house.change_status(i, INFECTIOUS);
				house.set_augm_value(i, "onset_date", t);
				update_augm_dependecies(house, i);

				// restore status id escaped infection
				if (uniform() < exp(-ROI_inst(house, i) * Dt))
				{
					house.change_status(i, SUSCEPTIBLE);
				}
			}
		}
	}

	// format all augmented variables of susceptibles
	for (indiv i : house.susceptible_set())
	{
		clear_all_augm(house, i);
	}
}