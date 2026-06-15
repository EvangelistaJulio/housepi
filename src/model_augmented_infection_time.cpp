/* housepi - Epidemiological Analysis of Household Outbreaks.
 *
 * Copyright © 2024-2026 Institut Pasteur (Paris), CNRS and INSERM.
 * Distributed under the terms of the GNU General Public License (GPLv3).
 * See the LICENSE file for details.
 */

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////  * *** * AUGMENT INFECTION TIMES - SINGLE INTRODUCTION * *** *  /////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

#include "housepi.h"

#ifndef _HOUSEPI_MODULE
#define _HOUSEPI_MODULE augmented_infection_time
#endif

HOUSEPI_OPEN_MODULE_NAMESPACE()

class Model final : public ModelBase
{
public:
	HOUSEPI_DECLARE_REGISTER_KEYS(AugmKey, infec_time)
	HOUSEPI_DECLARE_REGISTER_KEYS(DataKey, start_date, study_period)
	HOUSEPI_DECLARE_EMPTY_KEYS(LongKey)
	HOUSEPI_DECLARE_REGISTER_KEYS(GroupRateKey, community, susceptibility, infectivity)
	HOUSEPI_DECLARE_REGISTER_KEYS(HHuintKey, intro_case)
	HOUSEPI_DECLARE_EMPTY_KEYS(HHuintvecKey)
	HOUSEPI_DECLARE_EMPTY_KEYS(HHdoublevecKey)
	HOUSEPI_DECLARE_REGISTER_KEYS(HHdoublemtxKey, inst, cumul)
	HOUSEPI_DECLARE_REGISTER_KEYS(ParamKey, alpha, beta, gamma, size_ref)

	/// Initialize infectivity profile.
	void init_advanced_setup() override final
	{
#ifdef _HOUSEPI_DISCRETE_TIME
		request_augm_discrete(AugmKey::infec_time);
		request_infec_profile_discrete();
#endif
		request_param_discrete(ParamKey::size_ref);
		request_data_discrete(DataKey::start_date);
		request_data_discrete(DataKey::study_period);
		request_data_positive(DataKey::study_period);

		// template infectivity profile (otherwise user defined)
		// UPDATE: this should become more user friendly
		init_infec_profile_template({{"func", "gamma|mean:var|{2|N}:{1|N}"}});
	}

	/// Allow for contact rate definitions.
	bool has_contact_rates() const override final { return true; }

	/// NOTE: missing status are not modeled
	void init_model_for_evaluation() override final { assert_no_missing_status(); }

	/// Initialize household for evaluation
	void init_house_for_evaluation(Household &house) override final
	{
		if (house.n_infected_status() == 0)
			throw std::logic_error(
				"No infected member were found in " + house.err_name() +
				".\nThis model is compatible only with households with at least one infected member.");
	}

	///////////////////////////////////////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////////////////////
	///////////////////////  Model : compute single house log-likelihood  /////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////////////////////

	/// Transmission rate.
	double transmission_rate(size_t house_size) const
	{
		return par(ParamKey::beta) *
			   std::exp(par(ParamKey::gamma) * std::log(par(ParamKey::size_ref) / house_size));
	}

	/// Compute transmission log-likelihood within household.
	double compute_house_transm_log_lik(Household &house) override final
	{
		const size_t i0 = cache_uint(house, HHuintKey::intro_case);
		const Individual &member_0 = house.member(i0);
		const double t_intro = aug(AugmKey::infec_time, member_0);

		const auto &inst_mtx = cache_double_matrix(house, HHdoublemtxKey::inst);
		const auto &cumul_mtx = cache_double_matrix(house, HHdoublemtxKey::cumul);

		const double alpha = par(ParamKey::alpha);
		const double beta = transmission_rate(house.size());

		double output = 0.;
		for (size_t s : house.susceptible_hhindexes())
		{
			const Individual &member_s = house.member(s);
			const double exposure = member_s.data(DataKey::start_date) + member_s.data(DataKey::study_period) - t_intro;
			output -= alpha * exposure * member_s.group_rate(GroupRateKey::community);
		}

		for (size_t i : house.infected_hhindexes())
		{
			const Individual &member_i = house.member(i);
			const double mu_i = house.member(i).group_rate(GroupRateKey::infectivity);

			for (size_t s : house.susceptible_hhindexes())
				output -= mu_i * beta * cumul_mtx[i][s] *
						  house.contact_rate(i, s) *
						  house.member(s).group_rate(GroupRateKey::susceptibility);

			if (i == i0)
				continue;

			// Risk of infection from other household members k
			double roi_inst_i = 0.;
			double roi_cumul_i = 0.;
			for (size_t j : house.infected_hhindexes())
			{
				if (j == i)
					continue;

				const double mu_ji = house.member(j).group_rate(GroupRateKey::infectivity) *
									 house.contact_rate(j, i);
				roi_inst_i += mu_ji * inst_mtx[j][i];
				roi_cumul_i += mu_ji * cumul_mtx[j][i];
			}

			const double a_i = alpha * member_i.group_rate(GroupRateKey::community);
			const double b_i = beta * member_i.group_rate(GroupRateKey::susceptibility);

			roi_inst_i *= b_i;
			roi_cumul_i *= b_i;
			roi_inst_i += a_i;

#ifdef _HOUSEPI_DISCRETE_TIME
			roi_cumul_i += a_i * (aug(AugmKey::infec_time, member_i) - 1. - t_intro);
			output += std::log(1. - std::exp(-roi_inst_i)) - roi_cumul_i;
#else
			roi_cumul_i += a_i * (aug(AugmKey::infec_time, member_i) - t_intro);
			output += std::log(roi_inst_i) - roi_cumul_i;
#endif
		}

		return output;
	}

	///////////////////////////////////////////////////////////////////////////////////////////////////
	//////////////////////// Model : cache dependencies of augmented variables  ///////////////////////
	///////////////////////////////////////////////////////////////////////////////////////////////////

	void helper_update_pairwise_cache(
		std::vector<std::vector<double>> &inst_mtx,
		std::vector<std::vector<double>> &cumul_mtx,
		const Individual &member_i,
		const Individual &member_j,
		double t_infec_i,
		double t_infec_j)
	{
		const size_t i = member_i.hhindex();
		const size_t j = member_j.hhindex();
		const double gentime = t_infec_j - t_infec_i;

		// note: setting lut to 0 this is necessary because ROI is brute force
		if (gentime > 0)
		{ // i : infector, j : infectee
			inst_mtx[i][j] = infec_profile_inst(member_i, gentime);
#ifdef _HOUSEPI_DISCRETE_TIME
			cumul_mtx[i][j] = infec_profile_cumul(member_i, gentime - 1);
#else
			cumul_mtx[i][j] = infec_profile_cumul(member_i, gentime);
#endif
			inst_mtx[j][i] = 0.0;
			cumul_mtx[j][i] = 0.0;
		}
		else if (gentime < 0)
		{ // k : infector, i : infectee
			inst_mtx[i][j] = 0.0;
			cumul_mtx[i][j] = 0.0;
			inst_mtx[j][i] = infec_profile_inst(member_j, -gentime);
#ifdef _HOUSEPI_DISCRETE_TIME
			cumul_mtx[j][i] = infec_profile_cumul(member_j, -gentime - 1);
#else
			cumul_mtx[j][i] = infec_profile_cumul(member_j, -gentime);
#endif
		}
#ifdef _HOUSEPI_DISCRETE_TIME
		else
		{
			inst_mtx[i][j] = 0.0;
			cumul_mtx[i][j] = 0.0;
			inst_mtx[j][i] = 0.0;
			cumul_mtx[j][i] = 0.0;
		}
#endif
	}

	/// WARNING: unpredictable behavior if i is not an infector.
	/// WARNING: these function must not call any mutable parameter.
	/// Store cumulative and instantaneous infectivity profiles for household.
	void update_house_cache_dependecies(Household &house) override final
	{
		auto &inst_mtx = cache_double_matrix(house, HHdoublemtxKey::inst);
		auto &cumul_mtx = cache_double_matrix(house, HHdoublemtxKey::cumul);
		unsigned &i0 = cache_uint(house, HHuintKey::intro_case);
		double new_t0 = 100000;

		size_t n_infec = house.n_infected_status();
		const auto &infected_hhindexes = house.infected_hhindexes();

		for (size_t a = 0; a < n_infec; ++a)
		{
			const size_t i = infected_hhindexes[a];
			const Individual &member_i = house.member(i);
			const double t_infec_i = aug(AugmKey::infec_time, member_i);

			///////////////////////////////////
			//////  Susceptible Members  //////
			///////////////////////////////////

			for (size_t s : house.susceptible_hhindexes())
			{
				const Individual &member_s = house.member(s);
				double exposure = member_s.data(DataKey::start_date) + member_s.data(DataKey::study_period) - t_infec_i;
				cumul_mtx[i][s] = infec_profile_cumul(member_i, exposure);
				cumul_mtx[s][i] = 0.0;
			}

			///////////////////////////////////
			//////// Infectious Members ///////
			///////////////////////////////////

			for (size_t b = a + 1; b < n_infec; ++b)
			{
				const size_t j = infected_hhindexes[b];
				const Individual &member_j = house.member(j);
				const double t_infec_j = aug(AugmKey::infec_time, member_j);

				helper_update_pairwise_cache(inst_mtx, cumul_mtx,
											 member_i, member_j,
											 t_infec_i, t_infec_j);
			}

			// update household introduction case
			if (t_infec_i < new_t0)
			{
				new_t0 = t_infec_i;
				i0 = i;
			}
		}
	}

	/// Store cumulative and instantaneous infectivity profiles according to an infector.
	void update_member_cache_dependecies(Household &house, Individual &member_i) override final
	{
		auto &inst_mtx = cache_double_matrix(house, HHdoublemtxKey::inst);
		auto &cumul_mtx = cache_double_matrix(house, HHdoublemtxKey::cumul);
		unsigned &i0 = cache_uint(house, HHuintKey::intro_case);

		const size_t i = member_i.hhindex();
		const double t_infec_i = aug(AugmKey::infec_time, member_i);

		///////////////////////////////////
		//////  Susceptible Members  //////
		///////////////////////////////////

		for (size_t s : house.susceptible_hhindexes())
		{
			const Individual &member_s = house.member(s);
			double exposure = member_s.data(DataKey::start_date) + member_s.data(DataKey::study_period) - t_infec_i;
			cumul_mtx[i][s] = infec_profile_cumul(member_i, exposure);
			cumul_mtx[s][i] = 0.0;
		}

		///////////////////////////////////
		//////// Infectious Members ///////
		///////////////////////////////////

		if (house.n_infected_status() < 2)
			return;

		// commit introduction to current member
		i0 = i;
		double new_t0 = t_infec_i;

		for (size_t j : house.infected_hhindexes())
		{
			if (j == i)
				continue;

			const Individual &member_j = house.member(j);
			const double t_infec_j = aug(AugmKey::infec_time, member_j);

			helper_update_pairwise_cache(inst_mtx, cumul_mtx,
										 member_i, member_j,
										 t_infec_i, t_infec_j);

			// update household introduction case
			if (t_infec_j < new_t0)
			{
				new_t0 = t_infec_j;
				i0 = j;
			}
		}
	}

	///////////////////////////////////////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////////////////////
	////////////////////////////////  Model : simulate single house  //////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////////////////////

	/// Poisson point process to simulate outbreak in a single household.
	void simulate_single_house(Household &house) override final
	{
#ifdef _HOUSEPI_DISCRETE_TIME
		const double Dt = 1;
#else
		const double Dt = 0.05; // UDPDATE: it should be possible to change from outside
#endif
		/////////////////////////////////////////////////////////////////
		/////////////  Generation 0 : by single introduction  ///////////
		/////////////////////////////////////////////////////////////////

		size_t i0 = helper_rand_intro(house);
		Individual &member_0 = house.member(i0);

		/////////////////////////////////////////////////////////////////
		//////  Generation >=1 : infection by previous generation  //////
		/////////////////////////////////////////////////////////////////

		const double alpha = par(ParamKey::alpha);
		const double beta = transmission_rate(house.size());

		// largest study period
		double study_period = 0.;
		for (const Individual &member_i : house)
			study_period = std::max(study_period, member_i.data(DataKey::study_period));

		// FIXME: first approximation
		// this is assuming the anchor at time 0; wrong because:
		// 1 - anaware of first case infection time definition;
		// 2 - anchor cannot be specified by user;
		// 3 - a secondary case could actually have been the index;
		
		const double t0 = aug(AugmKey::infec_time, member_0);
		const double end_date = study_period; // FIXME: temporary patch

		// infection of secondary cases following the first case
		for (double t = t0 + Dt; t <= end_date && house.n_susceptible_status() > 0; t += Dt)
			for (Individual &member_i : house)
			{
				if (member_i.is_infected_status() || t > member_i.data(DataKey::study_period))
					continue;
				if (exponential() < roi_inst_simulation(house, member_i, t, alpha, beta) * Dt)
				{
					house.set_infected_status(member_i);
					set_augm_value(AugmKey::infec_time, member_i, t);
					init_missing_augm(member_i);
				}
			}
	}

	/// Draw introduction uniformly among all members.
	size_t helper_rand_intro(Household &house)
	{
		size_t i0 = uniform() * house.size();
		Individual &member_0 = house.member(i0);
		house.set_infected_status(member_0);
		init_missing_augm(member_0);

		// set all others members to susceptible
		for (Individual &member : house)
			if (member != member_0)
				house.set_susceptible_status(member);

		return i0;
	}

	/// Instantaneous risk of Infection for individual i. Used in simulation.
	double roi_inst_simulation(
		const Household &house,
		const Individual &member_i,
		double t_infec_i,
		double alpha,
		double beta) const
	{
		size_t i = member_i.hhindex();
		double roi = 0.; // instantaneous risk of infection

		// j : infector, i : infectee
		for (size_t j : house.infected_hhindexes())
		{
			const Individual &member_j = house.member(j);
			const double gentime = t_infec_i - aug(AugmKey::infec_time, member_j);
			roi += member_j.group_rate(GroupRateKey::infectivity) *
				   house.contact_rate(j, i) *
				   infec_profile_inst(member_j, gentime);
		}

		return alpha * member_i.group_rate(GroupRateKey::community) +
			   beta * member_i.group_rate(GroupRateKey::susceptibility) * roi;
	}
};
HOUSEPI_CLOSE_MODULE_NAMESPACE()
