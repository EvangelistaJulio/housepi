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
	UNUSED(house);
	UNUSED(from);
	UNUSED(to);
	emplace_contact_rate(REFERENCE);
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

/* Cumulative Risk of Infection for susceptible individual i. */
double Model::loglike_healthy(Household &house) const
{
	const indiv i0 = house.extra_integer_vector("intro-case")[0];
	const double a_cumul = param("alpha") * (param("study_period") + house.augm(i0, "incub_period"));
	const double b = transmission_rate(house);

	// load pre-compute inefctivity profiles
	const auto &cumul_mtx = house.extra_double_matrix("cumul");

	double output = 0.;

	for (indiv i : house.healthy_status_list())
	{
		double roi_cumul = 0;
		for (indiv k : house.infected_status_list())
			roi_cumul += cumul_mtx[k][i] * house.contact_rate(k, i) * house.group_rate(k, "infectivity");
		roi_cumul *= house.group_rate(i, "susceptibility") * b;
		roi_cumul += a_cumul;

		output -= roi_cumul;
	}

	return output;
}

/* Cumulative Risk of Infection for infectious individual i. */
double Model::loglike_infectious(Household &house) const
{
	double output = 0.;
	const indiv i0 = house.extra_integer_vector("intro-case")[0];
	double t_intro = house.augm(i0, "onset_date") - house.augm(i0, "incub_period");

	const auto &inst_mtx = house.extra_double_matrix("inst");
	const auto &cumul_mtx = house.extra_double_matrix("cumul");

	const double a = param("alpha");
	const double b = transmission_rate(house);

	for (indiv i : house.infected_status_list())
		if (i != i0)
		{
			// Risk of infection from other household members `k`
			double b_i = b * house.group_rate(i, "susceptibility");
			double roi_inst = 0.;
			double roi_cumul = 0.;
			for (indiv k : house.infected_status_list())
				if (cumul_mtx[k][i] != 0)
				{
					double mu_ki = house.contact_rate(k, i) * house.group_rate(k, "infectivity");
					roi_inst += inst_mtx[k][i] * mu_ki;
					roi_cumul += cumul_mtx[k][i] * mu_ki;
				}
			roi_inst *= b_i;
			roi_cumul *= b_i;

			// Risk of infection from community
			roi_inst += a;
			roi_cumul += a * (house.augm(i, "onset_date") - house.augm(i, "incub_period") - t_intro);

			output += log(roi_inst) - roi_cumul;
		}

	return output;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////// Model : store dependencies of augmented variables  ///////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/* Store cumulative and instantaneous infectivity profiles according to an infector.
   WARNING: unpredictable behavior if i is not an infector.
   WARNING: this function should not call any mutable parameter. */
void Model::update_extra_dependecies(Household &house, indiv i)
{
	const double onset_i = house.augm(i, "onset_date");
	const double incub_i = house.augm(i, "incub_period");
	const double t_infec_i = onset_i - incub_i;

	auto &inst_mtx = house.extra_double_matrix("inst");
	auto &cumul_mtx = house.extra_double_matrix("cumul");

	// init new first case to i
	double new_t0 = t_infec_i;
	indiv i0 = i;

	//==============================//
	//===== Infectious Members =====//
	//==============================//

	// note: setting lut to 0 this is necessary because _houseROI_ is brute
	for (indiv k : house.infected_status_list())
		if (k != i)
		{
			inst_mtx[i][k] = 0.0;
			cumul_mtx[i][k] = 0.0;
			inst_mtx[k][i] = 0.0;
			cumul_mtx[k][i] = 0.0;

			double incub_k = house.augm(k, "incub_period");
			double t_infec_k = house.augm(k, "onset_date") - incub_k;
			double gentime = t_infec_k - t_infec_i;

			if (gentime == 0)
				;
			else if (gentime > 0)
			{ //----- `i` : infector, `k` : infectee -----//
				inst_mtx[i][k] = infec_profile_inst(gentime, {incub_i});
				cumul_mtx[i][k] = infec_profile_cumul(gentime, {incub_i});
			}
			else
			{ //----- `k` : infector, `i` : infectee -----//
				inst_mtx[k][i] = infec_profile_inst(-gentime, {incub_k});
				cumul_mtx[k][i] = infec_profile_cumul(-gentime, {incub_k});

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

	double exposure = param("study_period") + house.augm(i0, "onset_date") - t_infec_i;
	for (indiv k : house.healthy_status_list())
	{
		cumul_mtx[i][k] = infec_profile_cumul(exposure, {incub_i});
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
	return loglike_healthy(house) + loglike_infectious(house);
}

//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//

///////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////// Model : risk of infection for simulation  ////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/* Instantaneous Risk of Infection for individual i (for simulation). */
double Model::ROI_inst_for_simulation(Household &house, indiv i, double t_infec_i) const
{
	double hh_roi_inst = 0.;

	for (indiv k : house.infected_status_list())
	{
		// `k` : infector, `i` : infectee
		double incub_k = house.augm(k, "incub_period");
		double gentime = t_infec_i - house.augm(k, "onset_date") + incub_k;
		double mu_ki = house.contact_rate(k, i) * house.group_rate(k, "infectivity");
		hh_roi_inst += infec_profile_inst(gentime, {incub_k}) * mu_ki;
	}
	hh_roi_inst *= house.group_rate(i, "susceptibility") * transmission_rate(house);

	return param("alpha") + hh_roi_inst;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////// Model : simulate single house  //////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/* Poisson point process to simulate outbreak in a single household. */
void Model::simulate_single_house(Household &house)
{
	double Dt = 0.05;

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

	// first infection occuring at time t = 0
	double tf = param("study_period") + house.augm(i0, "incub_period");
	house.set_augm_value(i0, "onset_date", house.augm(i0, "incub_period"));

	//=============================================================//
	//===== GENERATION >=1 : infection by previous generation =====//
	//=============================================================//

	// onset of secondary cases following the first case
	std::vector<indiv> new_gen_list = {};
	for (double t = 0; (t < tf) && (house.n_healthy_status() > 0); t += Dt)
	{
		for (indiv i : house.healthy_status_list())
		{
			double t_infec_i = t + Dt * uniform();
			if (uniform() >= exp(-ROI_inst_for_simulation(house, i, t_infec_i) * Dt))
			{
				house.set_augm_value(i, "onset_date",
									 t_infec_i + house.augm(i, "incub_period"));
				new_gen_list.push_back(i);
			}
		}

		// change household internal lists
		for (indiv i : new_gen_list)
			house.change_status(i, INFECTED_STATUS);
		new_gen_list.clear();
	}

	// format all augmented variables of susceptibles
	for (indiv i : house.healthy_status_list())
		clear_all_augm(house, i);
}
