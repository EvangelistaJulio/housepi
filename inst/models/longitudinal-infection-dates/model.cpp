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

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////  Model : initialization for simulation  ///////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/* Initialize infectivity profile. */
void Model::init_infec_profile()
{
	// default infectivity profile: negative binomial
	// otherwise user defined
	assign_default_infec_profile("nbinom");
}

/* Initialize parameters and augmented variables. */
void Model::init_model_for_simulation()
{
	// matrices for the cumulative infectivity profiles
	allocate_extra_house_double_matrix("inst");
	allocate_extra_house_double_matrix("cumul");

	// integer for index of introduction case
	allocate_extra_house_integer_vector("intro-case");
}

/* Initialize parameters and augmented variables. */
void Model::init_model_for_inference()
{
	init_model_for_simulation();

	// initialize to earliest infection date
	for (size_t H = 0; H < n_houses(); ++H)
	{
		Household &house = get_household(H);
		indiv i0 = house.size();
		for (indiv k : house.infectious_list())
			if (i0 == house.size())
				i0 = k;
			else if (house.data(k, "infec_date") < house.data(i0, "infec_date"))
				i0 = k;

		// store intro case
		house.extra_integer_vector("intro-case")[0] = i0;
	}
}

///////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////  Model : followup period  ////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/* Time when followup ended. */
double Model::time_end_followup(Household &house) const
{
	indiv i0 = house.extra_integer_vector("intro-case")[0];
	return house.data(i0, "infec_date") + param("study_period");
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
			output += house.group_rate(k, "infectivity") * profiles[k][i];

	output *= house.group_rate(i, "susceptibility") * transmission_rate(house);

	return output;
}

/* Instantaneous Risk of Infection for individual i. */
double Model::ROI_inst(Household &house, indiv i) const
{
	//----- Update Lookup Tables -----//
	auto &inst_mtx = house.extra_double_matrix("inst");
	const int t_infec_i = house.data(i, "infec_date");

	for (indiv k : house.infectious_list())
		if (k != i)
		{
			double gentime = t_infec_i - house.data(k, "infec_date");

			if (gentime > 0)
				inst_mtx[k][i] = infec_profile_inst(gentime);
		}

	return param("alpha") + _houseROI_(house, i, house.extra_double_matrix("inst"));
}

/* Cumulative Risk of Infection for infectious individual i. */
double Model::ROI_cumul_infec(Household &house, indiv i) const
{
	indiv i0 = house.extra_integer_vector("intro-case")[0];
	int t_intro = house.data(i0, "infec_date");
	int exposure = house.data(i, "infec_date") - t_intro;

	return param("alpha") * exposure + _houseROI_(house, i, house.extra_double_matrix("cumul"));
}

/* Cumulative Risk of Infection for susceptible individual i. */
double Model::ROI_cumul_suscept(Household &house, indiv i) const
{
	indiv i0 = house.extra_integer_vector("intro-case")[0];
	int t_intro = house.data(i0, "infec_date");
	int exposure = time_end_followup(house) - t_intro;

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
	auto &cumul_mtx = house.extra_double_matrix("cumul");

	//==============================//
	//===== Infectious Members =====//
	//==============================//

	const int t_infec_i = house.data(i, "infec_date");

	//----- Update Lookup Tables for infectious -----//
	for (indiv k : house.infectious_list())
		if (k != i)
		{
			double gentime = house.data(k, "infec_date") - t_infec_i;

			if (gentime > 0)
			{ // `i` : infector, `k` : infectee
				cumul_mtx[i][k] = infec_profile_inst(gentime);
				cumul_mtx[i][k] = infec_profile_cumul(gentime);
			}
			else
			{
				// `k` : infector, `i` : infectee
				cumul_mtx[k][i] = infec_profile_inst(-gentime);
				cumul_mtx[k][i] = infec_profile_cumul(-gentime);
			}
		}

	//===============================//
	//===== Susceptible Members =====//
	//===============================//

	const int exposure = time_end_followup(house) - t_infec_i;
	for (indiv k : house.susceptible_list())
		cumul_mtx[i][k] = infec_profile_cumul(exposure);
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

	// FIXME : double check loglikelihood
	for (indiv i : house.infectious_list())
		// Contribution of secondary cases
		if (i != house.extra_integer_vector("intro-case")[0])
			loglike_house += std::log(1. - std::exp(-ROI_inst(house, i))) - ROI_cumul_infec(house, i);

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
	//=================================================//
	//===== GENERATION 0 : by single introduction =====//
	//=================================================//

	for (indiv i = 0; i < house.size(); ++i)
		house.change_status(i, HEALTHY_STATUS);

	// pick a random introduction and change it's status
	indiv i0 = uniform() * house.size();
	house.change_status(i0, INFECTED_STATUS);
	house.extra_integer_vector("intro-case")[0] = i0;

	house.set_data_value(i0, "infec_date", 0);
	int tf = time_end_followup(house);

	//=============================================================//
	//===== GENERATION >=1 : infection by previous generation =====//
	//=============================================================//

	// onset of secondary cases following the first case
	for (int t = 0; (t < tf) && (house.n_susceptible() > 0); ++t)
		for (indiv i = 0; i < house.size(); ++i)
			if (house.status(i) == HEALTHY_STATUS)
			{
				// set new status for computations
				house.change_status(i, INFECTED_STATUS);
				update_extra_dependecies(house, i);
				house.set_data_value(i, "infec_date", t);

				// restore status id escaped infection
				if (uniform() < exp(-ROI_inst(house, i)))
					house.change_status(i, HEALTHY_STATUS);
			}
}
