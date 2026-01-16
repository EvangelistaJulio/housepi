/* housepi - Epidemiological Analysis of Household Outbreaks.
 *
 * Copyright © 2026 Institut Pasteur (Paris), CNRS and INSERM.
 * Distributed under the terms of the GNU General Public License (GPLv3).
 * See the LICENSE file for details.
 */

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////  * *** * AUGMENT INFECTION TIMES - SINGLE INTRODUCTION * *** *  /////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

#include "tools/keys.h"
#include "tools/keys_helper.h"

#include "individual.h"
#include "household.h"
#include "model.h"

REGISTER_HOUSEPI_KEYS(AugmKey, infec_time, incub_period, onset_time)
REGISTER_HOUSEPI_KEYS(DataKey, )
REGISTER_HOUSEPI_KEYS(GroupRateKey, community, susceptibility, infectivity)
REGISTER_HOUSEPI_KEYS(HHsize_tKey, intro_case)
REGISTER_HOUSEPI_KEYS(HHdoubleKey, )
REGISTER_HOUSEPI_KEYS(HHmatrixKey, inst, cumul)
REGISTER_HOUSEPI_KEYS(ParamKey, alpha, beta, gamma, size_ref, study_period)

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////  Model : initialization  //////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/// Initialize infectivity profile.
void Model::init_advanced_setup()
{
#ifdef _HOUSEPI_DISCRETE_TIME
	request_augm_discrete("infec_time");
	request_infec_profile_discrete();
#endif

	// template infectivity profile (otherwise user defined)
	// UPDATE: this should become more user friendly
	init_infec_profile_template({{"func", "gamma_meanvar|mean:var|{2|N}:{1|N}"}});
}

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////  Model : compute single house log-likelihood  /////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/// Compute the log-likelihood within household.
double Model::compute_house_obs_log_lik(Household &house)
{
	return loglike_healthy(house) + loglike_infectious(house);
}

/// Cumulative Risk of Infection for susceptible individual i.
double Model::loglike_healthy(const Household &house) const
{
	const auto &member_0 = house.get_member(house.cache_size_t_vector(HHsize_tKey::intro_case)[0]);
	const double alpha_cumul = par(ParamKey::alpha) * (par(ParamKey::study_period) + member_0.aug("incub_period"));
	const double beta = transmission_rate(house);

	// load pre-computed cumulative infectivity profiles
	const auto &cumul_mtx = house.cache_double_matrix(HHmatrixKey::cumul);

	double output = 0.;
	for (const auto &member_i : house.healthy())
	{
		double roi_cumul = 0;
		for (const auto &member_k : house.infected())
			roi_cumul += cumul_mtx[member_k][member_i] * member_k.group_rate(GroupRateKey::infectivity);
		output -= member_i.group_rate(GroupRateKey::community) * alpha_cumul + member_i.group_rate(GroupRateKey::susceptibility) * beta * roi_cumul;
	}

	return output;
}

/// Cumulative Risk of Infection for infectious individual i.
double Model::loglike_infectious(const Household &house) const
{
	double output = 0.;
	const auto &member_0 = house.get_member(house.cache_size_t_vector(HHsize_tKey::intro_case)[0]);
	const double t_intro = member_0.aug("infec_time");

	const auto &inst_mtx = house.cache_double_matrix(HHmatrixKey::inst);
	const auto &cumul_mtx = house.cache_double_matrix(HHmatrixKey::cumul);

	const double alpha = par(ParamKey::alpha);
	const double beta = transmission_rate(house);

	for (const auto &member_i : house.infected())
		if (member_i != member_0)
		{
			// Risk of infection from other household members `k`
			double b_i = beta * member_i.group_rate(GroupRateKey::susceptibility);
			double roi_inst_i = 0.;
			double roi_cumul_i = 0.;
			for (const auto &member_k : house.infected())
			{
				const double mu_k = member_k.group_rate(GroupRateKey::infectivity);
				roi_inst_i += inst_mtx[member_k][member_i] * mu_k;
				roi_cumul_i += cumul_mtx[member_k][member_i] * mu_k;
			}

			roi_inst_i *= b_i;
			roi_cumul_i *= b_i;

			// Risk of infection from community
			roi_inst_i += member_i.group_rate(GroupRateKey::community) * alpha;

#ifdef _HOUSEPI_DISCRETE_TIME
			roi_cumul_i += member_i.group_rate(GroupRateKey::community) * alpha * (member_i.aug("infec_time") - 1. - t_intro);
			output += std::log(1. - std::exp(-roi_inst_i)) - roi_cumul_i;
#else
			roi_cumul_i += member_i.group_rate(GroupRateKey::community) * alpha * (member_i.aug("infec_time") - t_intro);
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
	return par(ParamKey::beta) * std::pow(par(ParamKey::size_ref) / house.size(), par(ParamKey::gamma));
}

///////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////// Model : store dependencies of augmented variables  ///////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/// Store cumulative and instantaneous infectivity profiles according to an infector.
/// WARNING: unpredictable behavior if i is not an infector.
/// WARNING: this function should not call any mutable parameter.
void Model::update_cache_dependecies(Household &house, Individual &member_i)
{
	const double t_infec_i = member_i.aug("infec_time");

	auto &inst_mtx = house.cache_double_matrix(HHmatrixKey::inst);
	auto &cumul_mtx = house.cache_double_matrix(HHmatrixKey::cumul);

	double new_t0 = t_infec_i; // init time introduction to current case
	size_t i0 = member_i;	   // set index to current member

	///////////////////////////////////
	//////// Infectious Members ///////
	///////////////////////////////////

	// note: setting lut to 0 this is necessary because _houseROI_ is brute
	for (const auto &member_k : house.infected())
		if (member_k != member_i)
		{
			inst_mtx[member_i][member_k] = 0.0;
			cumul_mtx[member_i][member_k] = 0.0;
			inst_mtx[member_k][member_i] = 0.0;
			cumul_mtx[member_k][member_i] = 0.0;

			const double t_infec_k = member_k.aug("infec_time");
			const double gentime = t_infec_k - t_infec_i;

			if (gentime > 0)
			{ // `i` : infector, `k` : infectee
				inst_mtx[member_i][member_k] = infec_profile_inst(member_i, gentime);
#ifdef _HOUSEPI_DISCRETE_TIME
				cumul_mtx[member_i][member_k] = infec_profile_cumul(member_i, gentime - 1);
#else
				cumul_mtx[member_i][member_k] = infec_profile_cumul(member_i, gentime);
#endif
			}
			else if (gentime < 0)
			{ // `k` : infector, `i` : infectee
				inst_mtx[member_k][member_i] = infec_profile_inst(member_k, -gentime);
#ifdef _HOUSEPI_DISCRETE_TIME
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
	house.cache_size_t_vector(HHsize_tKey::intro_case)[0] = i0;

	///////////////////////////////////
	//////  Susceptible Members  //////
	///////////////////////////////////

	const double exposure = par(ParamKey::study_period) + house.get_member(i0).aug("onset_time") - t_infec_i;
	for (const auto &member_k : house.healthy())
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
void Model::simulate_single_house(Household &house)
{
#ifdef _HOUSEPI_DISCRETE_TIME
	const double Dt = 1;
#else
	const double Dt = 0.05; // UDPDATE: it should be possible to change from outside
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
	auto &member_0 = house.get_member(i0);
	house.set_infected_status(member_0);
	init_missing_augm(member_0);

	// set all to susceptible except from introduction
	for (auto &member : house)
		if (member != member_0)
			house.set_susceptible_status(member);

	/////////////////////////////////////////////////////////////////
	//////  Generation >=1 : infection by previous generation  //////
	/////////////////////////////////////////////////////////////////

	const double t0 = member_0.aug("infec_time");
	const double tf = member_0.aug("onset_time") + par(ParamKey::study_period);

	// infection of secondary cases following the first case
	for (double t = t0 + Dt; t <= tf && house.n_healthy_status() > 0; t += Dt)
		for (auto &member_i : house.healthy())
			if (uniform() >= exp(-roi_inst_simulation(house, member_i, t) * Dt))
			{
				house.set_infected_status(member_i);
				member_i.set_augm_value("infec_time", t);
				init_missing_augm(member_i);
			}
}

///////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////// Model : risk of infection for simulation  ////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/// Instantaneous Risk of Infection for individual i (for simulation).
double Model::roi_inst_simulation(Household &house, Individual &member_i, double t_infec_i) const
{
	double roi = 0.; // instantaneous risk of infection

	// `k` : infector, `i` : infectee
	for (const Individual &member_k : house.infected())
	{
		const double gentime = t_infec_i - member_k.aug("infec_time");
		const double mu_k = member_k.group_rate(GroupRateKey::infectivity);
		roi += infec_profile_inst(member_k, gentime) * mu_k;
	}
	roi *= member_i.group_rate(GroupRateKey::susceptibility) * transmission_rate(house);

	return par(ParamKey::alpha) + roi;
}
