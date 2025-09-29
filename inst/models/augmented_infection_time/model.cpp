/*
 *        __
 *       / /_  ____  __  __ ____ __  ____  __
 *      / __ \/ __ \/ / / / ___/ _ \/ __ \/ /
 *     / / / / /_/ / /_/ (__  )  __/ /_/ / /
 *    /_/ /_/\____/\____/____/\___/ ____/_/
 *                               /_/
 *
 * Copyright (C) 2025 Francesco Camaglia, Institut Pasteur
 *
 * model.cpp
 * Model template of the 'housepi' package.
 *
 */

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////  * *** * AUGMENT INFECTION TIMES - SINGLE INTRODUCTION * *** *  /////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

#include "model.h"

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////  Model : initialization  //////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/// Initialize infectivity profile.
void Model::init_advanced_setup()
{
	// template infectivity profile (otherwise user defined)
	init_infec_profile_template(
		{{"func", "gamma_meanvar|mean:var|{2|S}:{1|S}"}});
}

/// Initialize parameters and augmented variables.
void Model::init_model_for_inference()
{
	// integer for index of introduction case
	cache_house_integer_vector("intro_case");
	// matrices for instantaneous infectivity profiles
	cache_house_double_matrix("inst");
	// matrices for cumulative infectivity profiles
	cache_house_double_matrix("cumul");
}

/// Initialize parameters and augmented variables.
void Model::init_model_for_simulation()
{
	init_model_for_inference();
}

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////  Model : compute single house log-likelihood  /////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/// Compute the log-likelihood within household.
double Model::compute_household_obs_loglike(Household &house)
{
	return loglike_healthy(house) + loglike_infectious(house);
}

/// Cumulative Risk of Infection for susceptible individual i.
double Model::loglike_healthy(const Household &house) const
{
	const Individual &member_0 = house.get_member(house.cache_integer_vector("intro_case")[0]);
	const double alpha_cumul = param("alpha") * (param("study_period") + member_0.augm("incub_period"));
	const double beta = transmission_rate(house);

	// load pre-computed cumulative infectivity profiles
	const auto &cumul_mtx = house.cache_double_matrix("cumul");

	double output = 0.;
	for (const Individual &member_i : house.healthy())
	{
		double roi_cumul = 0;
		for (const Individual &member_k : house.infected())
			roi_cumul += cumul_mtx[member_k][member_i] * member_k.group_rate("infectivity");
		output -= alpha_cumul + member_i.group_rate("susceptibility") * beta * roi_cumul;
	}

	return output;
}

/// Cumulative Risk of Infection for infectious individual i.
double Model::loglike_infectious(const Household &house) const
{
	double output = 0.;
	const Individual &member_0 = house.get_member(house.cache_integer_vector("intro_case")[0]);
	const double t_intro = member_0.augm("infec_time");

	const auto &inst_mtx = house.cache_double_matrix("inst");
	const auto &cumul_mtx = house.cache_double_matrix("cumul");

	const double alpha = param("alpha");
	const double beta = transmission_rate(house);

	for (const Individual &member_i : house.infected())
		if (member_i != member_0)
		{
			// Risk of infection from other household members `k`
			double b_i = beta * member_i.group_rate("susceptibility");
			double roi_inst_i = 0.;
			double roi_cumul_i = 0.;
			for (const Individual &member_k : house.infected())
				if (cumul_mtx[member_k][member_i] != 0)
				{
					double mu_k = member_k.group_rate("infectivity");
					roi_inst_i += inst_mtx[member_k][member_i] * mu_k;
					roi_cumul_i += cumul_mtx[member_k][member_i] * mu_k;
				}

			roi_inst_i *= b_i;
			roi_cumul_i *= b_i;

			// Risk of infection from community
			roi_inst_i += alpha;

#ifdef _HSP_DISCRETE_TIME
			roi_cumul_i += alpha * (member_i.augm("infec_time") - 1. - t_intro);
			output += std::log(1. - std::exp(-roi_inst_i)) - roi_cumul_i;
#else
			roi_cumul_i += alpha * (member_i.augm("infec_time") - t_intro);
			output += std::log(roi_inst_i) - roi_cumul_i;
#endif
		}

	return output;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////// Model : transmission rate  ///////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/// [CHANGE-ME] Transmission rate.
double Model::transmission_rate(const Household &house) const
{
	return param("beta") * std::pow(param("size_ref") / house.size(), param("gamma"));
}

///////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////// Model : store dependencies of augmented variables  ///////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/// Store cumulative and instantaneous infectivity profiles according to an infector.
/// WARNING: unpredictable behavior if i is not an infector.
/// WARNING: this function should not call any mutable parameter.
void Model::update_cache_dependecies(Household &house, Individual &member_i)
{
	const double t_infec_i = member_i.augm("infec_time");

	matrix<double> &inst_mtx = house.cache_double_matrix("inst");
	matrix<double> &cumul_mtx = house.cache_double_matrix("cumul");

	double new_t0 = t_infec_i; // init time introduction to current case
	size_t i0 = member_i;	   // set index to current member

	///////////////////////////////////
	//////// Infectious Members ///////
	///////////////////////////////////

	// note: setting lut to 0 this is necessary because _houseROI_ is brute
	for (const Individual &member_k : house.infected())
		if (member_k != member_i)
		{
			inst_mtx[member_i][member_k] = 0.0;
			cumul_mtx[member_i][member_k] = 0.0;
			inst_mtx[member_k][member_i] = 0.0;
			cumul_mtx[member_k][member_i] = 0.0;

			const double t_infec_k = member_k.augm("infec_time");
			const double gentime = t_infec_k - t_infec_i;

			if (gentime > 0)
			{ // `i` : infector, `k` : infectee
				inst_mtx[member_i][member_k] = infec_profile_inst(member_i, gentime);
#ifdef _HSP_DISCRETE_TIME
				cumul_mtx[member_i][member_k] = infec_profile_cumul(member_i, gentime - 1);
#else
				cumul_mtx[member_i][member_k] = infec_profile_cumul(member_i, gentime);
#endif
			}
			else if (gentime < 0)
			{ // `k` : infector, `i` : infectee
				inst_mtx[member_k][member_i] = infec_profile_inst(member_k, -gentime);
#ifdef _HSP_DISCRETE_TIME
				cumul_mtx[member_k][member_i] = infec_profile_cumul(member_k, -gentime - 1);
#else
				cumul_mtx[member_k][member_i] = infec_profile_cumul(member_k, -gentime);
#endif

				// update first case
				if (t_infec_k < new_t0)
				{
					new_t0 = t_infec_k;
					i0 = member_k;
				}
			}
			// case gentime == 0 is skipped
		}

	// store intro case
	house.cache_integer_vector("intro_case")[0] = i0;

	///////////////////////////////////
	//////  Susceptible Members  //////
	///////////////////////////////////

	double exposure = param("study_period") + house.get_member(i0).augm("onset_time") - t_infec_i;
	for (const Individual &member_k : house.healthy())
	{
		cumul_mtx[member_i][member_k] = infec_profile_cumul(member_i, exposure);
		cumul_mtx[member_k][member_i] = 0.0;
	}
}

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////  Model : simulate single house  //////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/// Poisson point process to simulate outbreak in a single household.
void Model::simulate_single_household(Household &house)
{
#ifdef _HSP_DISCRETE_TIME
	const double Dt = 1;
#else
	const double Dt = HSP_TIME_STEP;
#endif

	/////////////////////////////////////////////////////////////////
	/////////////  Generation 0 : by single introduction  ///////////
	/////////////////////////////////////////////////////////////////

	// pick single introduction
	size_t i0;
	if (house.n_infec_status() > 0)
		// amongst pre-existing infectious as in data
		i0 = house.infected_hhindexes().at(uniform() * house.n_infec_status());
	else
		// amongst all members
		i0 = uniform() * house.size();
	const Individual &member_0 = house.get_member(i0);

	// set all susceptible except from introduction, then assign augmented variables
	for (Individual &member : house)
	{
		if (member == member_0)
			house.change_status(member, HSP_INFECTED_STATUS);
		else
			house.change_status(member, HSP_HEALTHY_STATUS);
		init_all_augm(member);
	}

	/////////////////////////////////////////////////////////////////
	//////  Generation >=1 : infection by previous generation  //////
	/////////////////////////////////////////////////////////////////

	const double t0 = member_0.augm("infec_time");
	const double tf = member_0.augm("onset_time") + param("study_period");

	// infection of secondary cases following the first case
	for (double t = t0 + Dt; (t < tf) && (house.n_healthy_status() > 0); t += Dt)
	{
		// store copies of infected/healthy individuals
		const std::vector<size_t> infected_hhindexes = house.infected_hhindexes();

		for (Individual &member_i : house.healthy())
			if (uniform() >= exp(-roi_inst_simulation(infected_hhindexes, house, member_i, t) * Dt))
			{
				member_i.set_augm_value("onset_time", t + member_i.augm("incub_period"));
				house.change_status(member_i, HSP_INFECTED_STATUS);
			}
	}

	// format infected with unseen onset as healthy
	for (Individual &member : house.infected())
		if (member.augm("onset_time") >= tf)
			house.change_status(member, HSP_HEALTHY_STATUS);

	// format all augmented variables of susceptibles
	for (Individual &member : house.healthy())
		clear_all_augm(member);
}

///////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////// Model : risk of infection for simulation  ////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/// Instantaneous Risk of Infection for individual i (for simulation).
double Model::roi_inst_simulation(const std::vector<size_t> &infected_hhindexes,
								  Household &house,
								  Individual &member_i,
								  double t_infec_i) const
{
	double hh_roi_inst = 0.;

	for (size_t k : infected_hhindexes)
	{
		// `k` : infector, `i` : infectee
		const Individual &member_k = house.get_member(k);
		double gentime = t_infec_i - member_k.augm("infec_time");
		double mu_k = member_k.group_rate("infectivity");
		hh_roi_inst += infec_profile_inst(member_k, gentime) * mu_k;
	}
	hh_roi_inst *= member_i.group_rate("susceptibility") * transmission_rate(house);

	return param("alpha") + hh_roi_inst;
}
