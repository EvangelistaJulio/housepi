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

/* Assign all group rates to each individual. */
void Model::assign_group_rate_to_individual(Individual &member)
{
	rateID r;
	UNUSED(member);

	//----- infectivity -----//
	r = "infectivity";
	emplace_group_rate(r, REFERENCE);

	//----- susceptibility -----//
	r = "susceptibility";
	emplace_group_rate(r, REFERENCE);
}

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////  Model : initialization  //////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/* Initialize infectivity profile. */
void Model::init_infec_profile()
{
	// default infectivity profile: gamma
	// otherwise user defined
	assign_template_infec_profile({{"name", "reparam_gamma"},
								   {"arg_name", "mean/var"},
								   {"arg_value", "2/1"},
								   {"arg_at", "S/S"}});
}

/* Initialize parameters and augmented variables. */
void Model::init_model_for_inference()
{
	// integer for index of introduction case
	allocate_cache_house_integer_vector("intro_case");

	// matrices for the instantaneous and cumulative infectivity profiles
	allocate_cache_house_double_matrix("inst");
	allocate_cache_house_double_matrix("cumul");
}

/* Initialize parameters and augmented variables. */
void Model::init_model_for_simulation()
{
	init_model_for_inference();
}

/* [CHANGE-ME] Transmission rate. */
double Model::transmission_rate(Household &house) const
{
	return param("beta") * std::pow(param("size_ref") / house.size(), param("gamma"));
}

/* Cumulative Risk of Infection for susceptible individual i. */
double Model::loglike_healthy(Household &house) const
{
	const indiv i0 = house.cache_integer_vector("intro_case")[0];
	const double a_cumul = param("alpha") * (param("study_period") + house.augm(i0, "incub_period"));
	const double b = transmission_rate(house);

	// load pre-compute inefctivity profiles
	const auto &cumul_mtx = house.cache_double_matrix("cumul");

	double output = 0.;

	for (indiv i : house.get_healthy_status_hhindexes())
	{
		double roi_cumul_i = 0;
		for (indiv k : house.get_infected_status_hhindexes())
			roi_cumul_i += cumul_mtx[k][i] * house.group_rate(k, "infectivity");
		roi_cumul_i *= house.group_rate(i, "susceptibility") * b;
		output -= roi_cumul_i;
	}
	output -= a_cumul * house.n_healthy_status();

	return output;
}

/* Cumulative Risk of Infection for infectious individual i. */
double Model::loglike_infectious(Household &house) const
{
	double output = 0.;
	const indiv i0 = house.cache_integer_vector("intro_case")[0];
	const double t_intro = house.augm(i0, "infec_time");

	const auto &inst_mtx = house.cache_double_matrix("inst");
	const auto &cumul_mtx = house.cache_double_matrix("cumul");

	const double a = param("alpha");
	const double b = transmission_rate(house);

	for (indiv i : house.get_infected_status_hhindexes())
		if (i != i0)
		{
			// Risk of infection from other household members `k`
			double b_i = b * house.group_rate(i, "susceptibility");
			double roi_inst_i = 0.;
			double roi_cumul_i = 0.;
			for (indiv k : house.get_infected_status_hhindexes())
				if (cumul_mtx[k][i] != 0)
				{
					double mu_k = house.group_rate(k, "infectivity");
					roi_inst_i += inst_mtx[k][i] * mu_k;
					roi_cumul_i += cumul_mtx[k][i] * mu_k;
				}
			roi_inst_i *= b_i;
			roi_cumul_i *= b_i;

			// Risk of infection from community
			roi_inst_i += a;
			roi_cumul_i += a * (house.augm(i, "infec_time") - t_intro);

			output += log(roi_inst_i) - roi_cumul_i;
		}

	return output;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////// Model : store dependencies of augmented variables  ///////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/* Store cumulative and instantaneous infectivity profiles according to an infector.
   WARNING: unpredictable behavior if i is not an infector.
   WARNING: this function should not call any mutable parameter. */
void Model::update_cache_dependecies(Household &house, indiv i)
{
	Individual &member_i = house.get_member(i);
	const double t_infec_i = member_i.augm("infec_time");

	auto &inst_mtx = house.cache_double_matrix("inst");
	auto &cumul_mtx = house.cache_double_matrix("cumul");

	// init new first case to i
	double new_t0 = t_infec_i;
	indiv i0 = i;

	//==============================//
	//===== Infectious Members =====//
	//==============================//

	// note: setting lut to 0 this is necessary because _houseROI_ is brute
	for (indiv k : house.get_infected_status_hhindexes())
		if (k != i)
		{
			Individual &member_k = house.get_member(k);

			inst_mtx[i][k] = 0.0;
			cumul_mtx[i][k] = 0.0;
			inst_mtx[k][i] = 0.0;
			cumul_mtx[k][i] = 0.0;

			double t_infec_k = member_k.augm("infec_time");
			double gentime = t_infec_k - t_infec_i;

			if (gentime > 0)
			{ //----- `i` : infector, `k` : infectee -----//
				inst_mtx[i][k] = infec_profile_inst(member_i, gentime);
				cumul_mtx[i][k] = infec_profile_cumul(member_i, gentime);
			}
			else if (gentime < 0)
			{ //----- `k` : infector, `i` : infectee -----//
				inst_mtx[k][i] = infec_profile_inst(member_k, -gentime);
				cumul_mtx[k][i] = infec_profile_cumul(member_k, -gentime);

				//----- update first case -----//
				if (t_infec_k < new_t0)
				{
					new_t0 = t_infec_k;
					i0 = k;
				}
			}
			// case gentime == 0 is skipped
		}

	// store intro case
	house.cache_integer_vector("intro_case")[0] = i0;

	//===============================//
	//===== Susceptible Members =====//
	//===============================//

	double exposure = param("study_period") + house.augm(i0, "onset_time") - t_infec_i;
	for (indiv k : house.get_healthy_status_hhindexes())
	{
		cumul_mtx[i][k] = infec_profile_cumul(member_i, exposure);
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
double Model::compute_house_obs_loglike(Household &house) const
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
double Model::roi_inst_simulation(const std::vector<indiv> &infectors,
								  Household &house, indiv i, double t_infec_i) const
{
	double hh_roi_inst = 0.;

	for (indiv k : infectors)
	{
		// `k` : infector, `i` : infectee
		Individual &member_k = house.get_member(k);
		double gentime = t_infec_i - member_k.augm("infec_time");
		double mu_k = member_k.group_rate("infectivity");
		hh_roi_inst += infec_profile_inst(member_k, gentime) * mu_k;
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
	const double Dt = 0.05;

	//=================================================//
	//===== GENERATION 0 : by single introduction =====//
	//=================================================//

	// pick single introduction
	indiv i0;
	if (house.n_infec_status() > 0)
		// amongst pre-existing infected in data
		i0 = house.get_infected_status_hhindexes().at(uniform() * house.n_infec_status());
	else
		// amongst all members
		i0 = uniform() * house.size();

	// set all to healthy except from introduction
	// assign augmented variables to all
	for (indiv i = 0; i < house.size(); ++i)
	{
		house.change_status(i, (i == i0) ? INFECTED_STATUS : HEALTHY_STATUS);
		assign_all_augm(house, i);
	}

	//=============================================================//
	//===== GENERATION >=1 : infection by previous generation =====//
	//=============================================================//

	// first infection occuring at time t = 0
	const double t0 = house.augm(i0, "infec_time");
	const double tf = house.augm(i0, "onset_time") + param("study_period");

	// infection of secondary cases following the first case
	for (double t = t0 + Dt; (t < tf) && (house.n_healthy_status() > 0); t += Dt)
	{
		// store copies of infected/healthy individuals
		const std::vector<indiv> infected = house.get_infected_status_hhindexes();
		const std::vector<indiv> healthy = house.get_healthy_status_hhindexes();

		for (indiv i : healthy)
			if (uniform() >= exp(-roi_inst_simulation(infected, house, i, t) * Dt))
			{
				house.set_augm_value(i, "onset_time", t + house.augm(i, "incub_period"));
				house.change_status(i, INFECTED_STATUS);
			}
	}

	// format infected with unseen onset as healthy
	for (indiv i : house.get_infected_status_hhindexes())
		if (house.augm(i, "onset_time") >= tf)
			house.change_status(i, HEALTHY_STATUS);

	// format all augmented variables of susceptibles
	for (indiv i : house.get_healthy_status_hhindexes())
		clear_all_augm(house, i);
}
