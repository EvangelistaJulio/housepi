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

#include "housepi_model_setup.h"
#include "model.h"

REGISTER_HOUSEPI_KEYS(AugmKey, infec_time)
REGISTER_HOUSEPI_KEYS(DataKey, end_followup)
EMPTY_HOUSEPI_KEYS(LongKey)
REGISTER_HOUSEPI_KEYS(GroupRateKey, community, susceptibility, infectivity)
REGISTER_HOUSEPI_KEYS(HHsizeKey, intro_case)
EMPTY_HOUSEPI_KEYS(HHsizevecKey)
EMPTY_HOUSEPI_KEYS(HHdoublevecKey)
REGISTER_HOUSEPI_KEYS(HHdoublemtxKey, inst, cumul)
REGISTER_HOUSEPI_KEYS(ParamKey, alpha, beta, gamma, size_ref)
HOUSEPI_MODEL_CONSTRUCTOR()

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////  Model : initialization  //////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

void Model::init_advanced_setup()
{
#ifdef _HOUSEPI_DISCRETE_TIME
	request_discrete_augm(AugmKey::infec_time);
	request_discrete_infec_profile();
#endif
	request_discrete_param(ParamKey::size_ref);

	// template infectivity profile (otherwise user defined)
	// UPDATE: this should become more user friendly
	init_infec_profile_template({{"func", "gamma_meanvar|mean:var|{2|N}:{1|N}"}});
}

void Model::init_house_for_evaluation(Household &house)
{
	size_t i0 = house.infected_hhindexes().front();
	size_t n_infec = house.n_infected_status();

	if (n_infec == 0)
		throw std::logic_error(
			"No infected member were found in " + house.err_name() +
			".\nThis model is only compatible with households with at least one infected member.");

	if (n_infec > 1)
	{
		double t_intro = house.member(i0).aug(AugmKey::infec_time);
		for (size_t a = 1; a < house.n_infected_status(); ++a)
		{
			const size_t i = house.infected_hhindexes()[a];
			double t_infec_i = house.member(i).aug(AugmKey::infec_time);
			if (t_infec_i < t_intro)
			{
				t_intro = t_infec_i;
				i0 = i;
			}
		}
	}

	cache_size(house, HHsizeKey::intro_case) = i0;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////  Model : compute single house log-likelihood  /////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

double Model::transmission_rate(size_t house_size) const
{
	return par(ParamKey::beta) * std::exp(par(ParamKey::gamma) * std::log(par(ParamKey::size_ref) / house_size));
}

double Model::compute_house_obs_log_lik(Household &house)
{
	const size_t i0 = cache_size(house, HHsizeKey::intro_case);
	const Individual &member_0 = house.member(i0);
	const double t_intro = member_0.aug(AugmKey::infec_time);

	const auto &inst_mtx = cache_double_matrix(house, HHdoublemtxKey::inst);
	const auto &cumul_mtx = cache_double_matrix(house, HHdoublemtxKey::cumul);

	const double alpha = par(ParamKey::alpha);
	const double beta = transmission_rate(house.size());

	double output = 0.;
	for (size_t k : house.susceptible_hhindexes())
	{
		const Individual &member_k = house.member(k);
		const double alpha_cumul = alpha * (member_k.data(DataKey::end_followup) - t_intro);
		output -= alpha_cumul * member_k.group_rate(GroupRateKey::community);
	}

	for (size_t i : house.infected_hhindexes())
	{
		const Individual &member_i = house.member(i);
		const double mu_i = house.member(i).group_rate(GroupRateKey::infectivity);

		for (size_t k : house.susceptible_hhindexes())
			output -= mu_i * beta * cumul_mtx[i][k] *
					  house.member(k).group_rate(GroupRateKey::susceptibility);

		if (i == i0)
			continue;

		// Risk of infection from other household members `k`
		double roi_inst_i = 0.;
		double roi_cumul_i = 0.;
		for (size_t k : house.infected_hhindexes())
		{
			if (k == i)
				continue;

			const double mu_k = house.member(k).group_rate(GroupRateKey::infectivity);
			roi_inst_i += mu_k * inst_mtx[k][i];
			roi_cumul_i += mu_k * cumul_mtx[k][i];
		}

		const double a_i = alpha * member_i.group_rate(GroupRateKey::community);
		const double b_i = beta * member_i.group_rate(GroupRateKey::susceptibility);

		roi_inst_i *= b_i;
		roi_cumul_i *= b_i;
		roi_inst_i += a_i;

#ifdef _HOUSEPI_DISCRETE_TIME
		roi_cumul_i += a_i * (member_i.aug(AugmKey::infec_time) - 1. - t_intro);
		output += std::log(1. - std::exp(-roi_inst_i)) - roi_cumul_i;
#else
		roi_cumul_i += a_i * (member_i.aug(AugmKey::infec_time) - t_intro);
		output += std::log(roi_inst_i) - roi_cumul_i;
#endif
	}

	return output;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////// Model : store dependencies of augmented variables  ///////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

inline void Model::helper_update_pairwise_cache(
	std::vector<std::vector<double>> &inst_mtx,
	std::vector<std::vector<double>> &cumul_mtx,
	const Individual &member_i, const Individual &member_k,
	double t_infec_i, double t_infec_k)
{
	const size_t i = member_i.hhindex();
	const size_t k = member_k.hhindex();
	const double gentime = t_infec_k - t_infec_i;

	// note: setting lut to 0 this is necessary because ROI is brute force
	if (gentime > 0)
	{ // `i` : infector, `k` : infectee
		inst_mtx[i][k] = infec_profile_inst(member_i, gentime);
#ifdef _HOUSEPI_DISCRETE_TIME
		cumul_mtx[i][k] = infec_profile_cumul(member_i, gentime - 1);
#else
		cumul_mtx[i][k] = infec_profile_cumul(member_i, gentime);
#endif
		inst_mtx[k][i] = 0.0;
		cumul_mtx[k][i] = 0.0;
	}
	else if (gentime < 0)
	{ // `k` : infector, `i` : infectee
		inst_mtx[i][k] = 0.0;
		cumul_mtx[i][k] = 0.0;
		inst_mtx[k][i] = infec_profile_inst(member_k, -gentime);
#ifdef _HOUSEPI_DISCRETE_TIME
		cumul_mtx[k][i] = infec_profile_cumul(member_k, -gentime - 1);
#else
		cumul_mtx[k][i] = infec_profile_cumul(member_k, -gentime);
#endif
	}
	else
	{
		inst_mtx[i][k] = 0.0;
		cumul_mtx[i][k] = 0.0;
		inst_mtx[k][i] = 0.0;
		cumul_mtx[k][i] = 0.0;
	}
}

void Model::update_house_cache_dependecies(Household &house)
{
	auto &inst_mtx = cache_double_matrix(house, HHdoublemtxKey::inst);
	auto &cumul_mtx = cache_double_matrix(house, HHdoublemtxKey::cumul);
	size_t &i0 = cache_size(house, HHsizeKey::intro_case);
	double new_t0 = 1000000;

	size_t n_infec = house.n_infected_status();
	const auto &infected_hhindexes = house.infected_hhindexes();

	for (size_t a = 0; a < n_infec; ++a)
	{
		const size_t i = infected_hhindexes[a];
		const Individual &member_i = house.member(i);
		const double t_infec_i = member_i.aug(AugmKey::infec_time);

		///////////////////////////////////
		//////  Susceptible Members  //////
		///////////////////////////////////

		for (size_t k : house.susceptible_hhindexes())
		{
			double exposure = house.member(k).data(DataKey::end_followup) - t_infec_i;
			cumul_mtx[i][k] = infec_profile_cumul(member_i, exposure);
			cumul_mtx[k][i] = 0.0;
		}

		///////////////////////////////////
		//////// Infectious Members ///////
		///////////////////////////////////

		for (size_t b = a + 1; b < n_infec; ++b)
		{
			const size_t k = infected_hhindexes[b];
			const Individual &member_k = house.member(k);
			const double t_infec_k = member_k.aug(AugmKey::infec_time);

			helper_update_pairwise_cache(inst_mtx, cumul_mtx,
										 member_i, member_k,
										 t_infec_i, t_infec_k);
		}

		// update household introduction case
		if (t_infec_i < new_t0)
		{
			new_t0 = t_infec_i;
			i0 = i;
		}
	}
}

void Model::update_member_cache_dependecies(Household &house, Individual &member_i)
{
	auto &inst_mtx = cache_double_matrix(house, HHdoublemtxKey::inst);
	auto &cumul_mtx = cache_double_matrix(house, HHdoublemtxKey::cumul);
	size_t &i0 = cache_size(house, HHsizeKey::intro_case); // set index to current member

	const size_t i = member_i.hhindex();
	const double t_infec_i = member_i.aug(AugmKey::infec_time);

	///////////////////////////////////
	//////  Susceptible Members  //////
	///////////////////////////////////

	for (size_t k : house.susceptible_hhindexes())
	{
		double exposure = house.member(k).data(DataKey::end_followup) - t_infec_i;
		cumul_mtx[i][k] = infec_profile_cumul(member_i, exposure);
		cumul_mtx[k][i] = 0.0;
	}

	///////////////////////////////////
	//////// Infectious Members ///////
	///////////////////////////////////

	if (house.n_infected_status() < 2)
		return;

	double new_t0 = t_infec_i;

	for (size_t k : house.infected_hhindexes())
	{
		if (k == i)
			continue;

		const Individual &member_k = house.member(k);
		const double t_infec_k = member_k.aug(AugmKey::infec_time);

		helper_update_pairwise_cache(inst_mtx, cumul_mtx,
									 member_i, member_k,
									 t_infec_i, t_infec_k);

		// update household introduction case
		if (t_infec_k < new_t0)
		{
			new_t0 = t_infec_k;
			i0 = k;
		}
	}
}

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////  Model : simulate single house  //////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

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
	if (house.n_infected_status() > 0)
		// amongst pre-existing infectious as in data
		i0 = *std::next(house.infected_hhindexes().begin(), uniform() * house.n_infected_status());
	else
		// amongst all members
		i0 = uniform() * house.size();
	Individual &member_0 = house.member(i0);
	house.set_infected_status(member_0);
	init_missing_augm(member_0);

	// set all to susceptible except from introduction
	for (Individual &member : house)
		if (member != member_0)
			house.set_susceptible_status(member);

	/////////////////////////////////////////////////////////////////
	//////  Generation >=1 : infection by previous generation  //////
	/////////////////////////////////////////////////////////////////

	const double alpha = par(ParamKey::alpha);
	const double beta = transmission_rate(house.size());

	const double t0 = member_0.aug(AugmKey::infec_time);
	/// FIXME: assuming internally all end_followup are the same,
	///        but this is not requested to the user!
	const double tf = member_0.data(DataKey::end_followup);

	// infection of secondary cases following the first case
	for (double t = t0 + Dt; t <= tf && house.n_susceptible_status() > 0; t += Dt)
		for (Individual &member_i : house)
		{
			if (member_i.is_infected_status())
				continue;
			if (exponential() < roi_inst_simulation(house, member_i, t, alpha, beta) * Dt)
			{
				house.set_infected_status(member_i);
				member_i.set_augm_value(AugmKey::infec_time, t);
				init_missing_augm(member_i);
			}
		}
}

double Model::roi_inst_simulation(const Household &house, const Individual &member_i,
								  double t_infec_i, double alpha, double beta) const
{
	double roi = 0.; // instantaneous risk of infection

	// `k` : infector, `i` : infectee
	for (size_t k : house.infected_hhindexes())
	{
		const Individual &member_k = house.member(k);
		const double gentime = t_infec_i - member_k.aug(AugmKey::infec_time);
		const double mu_k = member_k.group_rate(GroupRateKey::infectivity);
		roi += infec_profile_inst(member_k, gentime) * mu_k;
	}
	roi *= member_i.group_rate(GroupRateKey::susceptibility) * beta;

	return alpha + roi;
}
