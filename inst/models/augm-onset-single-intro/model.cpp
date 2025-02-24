/*
 * model.cpp
 *
 * Model template of the 'housepi' package.
 * Copyright (C) 2025 Francesco Camaglia, Institut Pasteur
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

/* [CHANGE-ME] Time of infection */
double Model::transmission_rate(Household &house) const
{
	return param("beta") * pow(param("size_ref") / house.size(), param("gamma"));
}

//------------------------------------------- /\ ------------------------------------------------//
//------------------------------------------- || ------------------------------------------------//
//------------------------------------------- || ------------------------------------------------//
//------------------------------------------- || ------------------------------------------------//
//---------------------------------------- CHANGE-ME --------------------------------------------//
//-----------------------------------------------------------------------------------------------//

/* Assign all group rates for each individual. */
void Model::assign_group_rate_to_individual(Household &house, indiv i)
{
	rateID r;
	UNUSED(house);
	UNUSED(i);

	//----- infectivity -----//
	r = "infectivity";
	emplace_group_rate(r, REFERENCE);

	//----- susceptibility -----//
	r = "susceptibility";
	emplace_group_rate(r, REFERENCE);
}

/* Assign the contact patterns between a couple of individuals. */
void Model::assign_contact_rate_to_household(Household &house, indiv from, indiv to)
{
	// child if younger than 'child_age'; adult otherwise
	bool is_from_adult = (house.data(from, "age") < param("child_age"));
	bool is_to_adult = (house.data(to, "age") < param("child_age"));

	if (is_from_adult && is_to_adult)
		emplace_contact_rate(REFERENCE);
	else if (is_from_adult && (!is_to_adult))
		emplace_contact_rate("cont_adult_child");
	else if ((!is_from_adult) && is_to_adult)
		emplace_contact_rate("cont_child_adult");
	else
		emplace_contact_rate("cont_child_child");
}

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////  Model : initialization  //////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/* Initialize infectivity profile. */
void Model::init_infec_profile()
{
	// default infectivity profile: Ferretti et al (2020)
	// otherwise user defined
	assign_default_infec_profile("ferretti20");
}

/* Initialize parameters and augmented variables. */
void Model::init_model_for_inference()
{
	// integer for index of introduction case
	allocate_extra_house_integer_vector("intro-case");

	// matrices for the instantaneous and cumulative infectivity profiles
	allocate_extra_house_double_matrix("inst");
	allocate_extra_house_double_matrix("cumul");
}

/* Initialize parameters and augmented variables. */
void Model::init_model_for_simulation()
{
	init_model_for_inference();
}

///////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////  Model : followup period  ////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/* Time when followup ended. */
double Model::time_end_followup(Household &house) const
{
	indiv i0 = house.extra_integer_vector("intro-case")[0];
	return house.augm(i0, "onset_date") + param("study_period");
}

///////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////  Model : force of infection  ///////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/* Sum over all hazards of infection. */
double Model::_houseROI_(Household &house, indiv i, const matrix<double> &profiles) const
{
	double output = 0.0;

	// Risk of infection from other household members (k)
	for (indiv k : house.infectious_list())
		if (k != i)
			output += house.contact_rate(k, i) * house.group_rate(k, "infectivity") * profiles[k][i];

	output *= house.group_rate(i, "susceptibility") * transmission_rate(house);

	return output;
}

/* Instantaneous Risk of Infection for individual i. */
double Model::ROI_inst(Household &house, indiv i) const
{
	return param("alpha") + _houseROI_(house, i, house.extra_double_matrix("inst"));
}

/* Cumulative Risk of Infection for infectious individual i. */
double Model::ROI_cumul_infec(Household &house, indiv i) const
{
	indiv i0 = house.extra_integer_vector("intro-case")[0];
	double t_intro = house.augm(i0, "onset_date") - house.augm(i0, "incub_period");
	double exposure = house.augm(i, "onset_date") - house.augm(i, "incub_period") - t_intro;

	return param("alpha") * exposure + _houseROI_(house, i, house.extra_double_matrix("cumul"));
}

/* Cumulative Risk of Infection for susceptible individual i. */
double Model::ROI_cumul_suscept(Household &house, indiv i) const
{
	indiv i0 = house.extra_integer_vector("intro-case")[0];
	double t_intro = house.augm(i0, "onset_date") - house.augm(i0, "incub_period");
	double exposure = time_end_followup(house) - t_intro;

	return param("alpha") * exposure + _houseROI_(house, i, house.extra_double_matrix("cumul"));
}

///////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////// Model : store dependencies of augmented variables  ///////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/* Store cumulative and instantaneous infectivity profiles according to an infector.
   WARNING: unpredictable behavior if i is not an infector.
   WARNING: this function should not call any mutable parameter. */
void Model::update_extra_dependecies(Household &house, indiv i)
{

	auto &inst_mtx = house.extra_double_matrix("inst");
	auto &cumul_mtx = house.extra_double_matrix("cumul");

	//==============================//
	//===== Infectious Members =====//
	//==============================//

	const double onset_i = house.augm(i, "onset_date");
	const double incub_i = house.augm(i, "incub_period");
	const double t_infec_i = onset_i - incub_i;

	// init new first case to i
	double new_t0 = t_infec_i;
	indiv i0 = i;

	//----- Update Lookup Tables for infectious -----//
	double onset_k, incub_k, t_infec_k;
	double gentime; // generation time

	// note: setting lut to 0 this is necessary because _houseROI_ is brute

	for (indiv k : house.infectious_list())
		if (k != i)
		{
			onset_k = house.augm(k, "onset_date");
			incub_k = house.augm(k, "incub_period");
			t_infec_k = onset_k - incub_k;

			if (t_infec_i < t_infec_k)
			{ //----- `i` : infector, `k` : infectee -----//
				gentime = t_infec_k - t_infec_i;

				inst_mtx[i][k] = infec_profile_inst(gentime, {incub_i, incub_k});
				inst_mtx[k][i] = 0.0;

				cumul_mtx[i][k] = infec_profile_cumul(gentime, {incub_i, incub_k});
				cumul_mtx[k][i] = 0.0;
			}
			else
			{ //----- `k` : infector, `i` : infectee -----//
				gentime = t_infec_i - t_infec_k;

				inst_mtx[i][k] = 0.0;
				inst_mtx[k][i] = infec_profile_inst(gentime, {incub_k, incub_i});

				cumul_mtx[i][k] = 0.0;
				cumul_mtx[k][i] = infec_profile_cumul(gentime, {incub_k, incub_i});

				//----- update first case -----//
				if (t_infec_k < new_t0)
				{
					new_t0 = t_infec_k;
					i0 = k;
				}
			}
		}

	// store intro case
	house.extra_integer_vector("intro-case")[0] = i0;

	//===============================//
	//===== Susceptible Members =====//
	//===============================//

	gentime = time_end_followup(house) - t_infec_i;
	for (indiv k : house.susceptible_list())
	{
		incub_k = house.augm(k, "incub_period");
		cumul_mtx[i][k] = infec_profile_cumul(gentime, {incub_i, incub_k});
		cumul_mtx[k][i] = 0.0;
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

	for (indiv i : house.susceptible_list())
		loglike_house -= ROI_cumul_suscept(house, i);

	//==============================//
	//===== Infectious Members =====//
	//==============================//

	for (indiv i : house.infectious_list())
		// Contribution of secondary cases
		if (i != house.extra_integer_vector("intro-case")[0])
			loglike_house += log(ROI_inst(house, i)) - ROI_cumul_infec(house, i);

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
		if (house.status(i) != INFECTED_STATUS)
			assign_all_augm(house, i);

		house.change_status(i, HEALTHY_STATUS);
	}

	// pick a random introduction and change it's status
	indiv i0 = uniform() * house.size();
	house.change_status(i0, INFECTED_STATUS);
	house.extra_integer_vector("intro-case")[0] = i0;

	double t = house.augm(i0, "onset_date");
	double tf = time_end_followup(house);

	//=============================================================//
	//===== GENERATION >=1 : infection by previous generation =====//
	//=============================================================//

	// onset of secondary cases following the first case
	for (; (t < tf) && (house.n_susceptible() > 0); t += Dt)
		for (indiv i = 0; i < house.size(); ++i)
		{
			if (house.status(i) == HEALTHY_STATUS)
			{
				// set new status for computations
				house.change_status(i, INFECTED_STATUS);
				house.set_augm_value(i, "onset_date", t);
				update_extra_dependecies(house, i);

				// restore status id escaped infection
				if (uniform() < exp(-ROI_inst(house, i) * Dt))
					house.change_status(i, HEALTHY_STATUS);
			}
		}

	// format all augmented variables of susceptibles
	for (indiv i : house.susceptible_list())
		clear_all_augm(house, i);
}
