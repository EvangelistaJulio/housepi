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
#define _HOUSEPI_MODULE case_ascertained
#endif

HOUSEPI_OPEN_MODULE_NAMESPACE()

class Model final : public ModelBase
{
public:
	HOUSEPI_DECLARE_REGISTER_KEYS(Augm, infec_time)
	HOUSEPI_DECLARE_REGISTER_KEYS(Data, start_date, end_date)
	HOUSEPI_DECLARE_EMPTY_KEYS(Long)
	HOUSEPI_DECLARE_REGISTER_KEYS(GroupRisk, community, susceptibility, infectivity)
	HOUSEPI_DECLARE_REGISTER_KEYS(HHmember, intro_case)
	HOUSEPI_DECLARE_EMPTY_KEYS(HHmembervec)
	HOUSEPI_DECLARE_EMPTY_KEYS(HHdoublevec)
	HOUSEPI_DECLARE_REGISTER_KEYS(HHdoublemtx, inst, cumul)
	HOUSEPI_DECLARE_REGISTER_KEYS(Param, alpha, beta, gamma, hhsize_ref)

	/// Initialize infectivity profile.
	void init_model_requirements() override final
	{
		request_augm_transmission(
			{Augm::infec_time});
		request_param_discrete(
			{Param::hhsize_ref});
		request_data_discrete(
			{Data::start_date, Data::end_date});
		request_risk_identifiability(
			Param::alpha, {GroupRisk::community});
		request_risk_identifiability(
			Param::beta, {GroupRisk::susceptibility, GroupRisk::infectivity}, true);
#ifdef _HOUSEPI_DISCRETE_TIME
		request_augm_discrete(
			{Augm::infec_time});
		request_infec_profile_discrete();
#endif
	}

	/// Allow for contact risk definitions.
	bool has_contact_risks() const override final { return true; }

	/// Transmission is driven by an infectivity profile.
	bool has_infec_profile() const override final { return true; }

	/// NOTE: missing status are not modeled
	void init_model_for_evaluation() override final { assert_no_missing_status(); }

	/// Initialize household for evaluation
	void init_house_for_evaluation(const Household &house) override final
	{
		if (house.n_infected_status() == 0)
			throw std::logic_error(
				"No infected members were found in " + err_name(house) +
				".\nCase-ascertained model assumes at least one infected member per household.");

		// study period must be non-negative
		for (const Individual &indiv : house)
			if (data(Data::start_date, indiv) > data(Data::end_date, indiv))
				throw std::logic_error(
					"Found negative follow-up period for " + err_name(indiv) +
					".\nCheck " + msg::arg("start_date") + " and " + msg::arg("end_date") + ".");
	}

	///////////////////////////////////////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////////////////////
	///////////////////////  Model : compute single house log-likelihood  /////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////////////////////

	/// Transmission rate.
	double transmission_rate(MemberIndex house_size) const
	{
		return par(Param::beta) *
			   std::exp(par(Param::gamma) * std::log(par(Param::hhsize_ref) / house_size));
	}

	/// End of follow-up `end_date` is the last fully-observed day (inclusive)
	double follow_up_end(const Individual &indiv) const
	{
#ifdef _HOUSEPI_DISCRETE_TIME
		return data(Data::end_date, indiv);
#else
		return data(Data::end_date, indiv) + 1.;
#endif
	}

	/// Compute transmission log-likelihood within household.
	double compute_house_transm_log_lik(const Household &house) override final
	{
		const MemberIndex i0 = cache_member(house, HHmember::intro_case);
		const Individual &indiv_0 = house.at(i0);
		const double t_intro = aug(Augm::infec_time, indiv_0);

		const auto &inst_mtx = cache_double_matrix(house, HHdoublemtx::inst);
		const auto &cumul_mtx = cache_double_matrix(house, HHdoublemtx::cumul);

		const double alpha = par(Param::alpha);
		const double beta = transmission_rate(house.size());

		double output = 0.;
		for (const Individual &indiv_s : house.escaped_indivs())
		{
			const double exposure = follow_up_end(indiv_s) - t_intro;
			output -= alpha * exposure * group_risk(GroupRisk::community, indiv_s);
		}

		for (MemberIndex i : house.infected_members())
		{
			const Individual &indiv_i = house.at(i);
			const double mu_i = group_risk(GroupRisk::infectivity, indiv_i);

			for (MemberIndex s : house.escaped_members())
				output -= mu_i * beta * cumul_mtx[i][s] *
						  contact_risk(house, i, s) *
						  group_risk(GroupRisk::susceptibility, house.at(s));

			if (i == i0)
				continue;

			// Risk of infection from other household members k
			double roi_inst_i = 0.;
			double roi_cumul_i = 0.;
			for (MemberIndex j : house.infected_members())
			{
				if (j == i)
					continue;

				const double mu_ji = group_risk(GroupRisk::infectivity, house.at(j)) *
									 contact_risk(house, j, i);
				roi_inst_i += mu_ji * inst_mtx[j][i];
				roi_cumul_i += mu_ji * cumul_mtx[j][i];
			}

			const double a_i = alpha * group_risk(GroupRisk::community, indiv_i);
			const double b_i = beta * group_risk(GroupRisk::susceptibility, indiv_i);

			roi_inst_i *= b_i;
			roi_cumul_i *= b_i;
			roi_inst_i += a_i;

#ifdef _HOUSEPI_DISCRETE_TIME
			// days escaped infection before transmission:
			// a member infected the day after the introduction survived none
			// NOTE: co-primary cannot survive fewer
			const double survived = aug(Augm::infec_time, indiv_i) - 1. - t_intro;
			if (survived > 0.)
				roi_cumul_i += a_i * survived;
			output += std::log(1. - std::exp(-roi_inst_i)) - roi_cumul_i;
#else
			roi_cumul_i += a_i * (aug(Augm::infec_time, indiv_i) - t_intro);
			output += std::log(roi_inst_i) - roi_cumul_i;
#endif
		}

		return output;
	}

	///////////////////////////////////////////////////////////////////////////////////////////////////
	//////////////////// Model : cached quantities consuming augmented variables  ////////////////////
	///////////////////////////////////////////////////////////////////////////////////////////////////

	void helper_update_pairwise_cache(
		std::vector<std::vector<double>> &inst_mtx,
		std::vector<std::vector<double>> &cumul_mtx,
		const Individual &indiv_i,
		const Individual &indiv_j,
		double t_infec_i,
		double t_infec_j)
	{
		const MemberIndex i = Index::member(indiv_i);
		const MemberIndex j = Index::member(indiv_j);
		const double gentime = t_infec_j - t_infec_i;

		// note: setting lut to 0 this is necessary because ROI is brute force
		if (gentime > 0)
		{ // i : infector, j : infectee
			inst_mtx[i][j] = infec_profile_inst(indiv_i, gentime);
#ifdef _HOUSEPI_DISCRETE_TIME
			cumul_mtx[i][j] = infec_profile_cumul(indiv_i, gentime - 1);
#else
			cumul_mtx[i][j] = infec_profile_cumul(indiv_i, gentime);
#endif
			inst_mtx[j][i] = 0.0;
			cumul_mtx[j][i] = 0.0;
		}
		else if (gentime < 0)
		{ // k : infector, i : infectee
			inst_mtx[i][j] = 0.0;
			cumul_mtx[i][j] = 0.0;
			inst_mtx[j][i] = infec_profile_inst(indiv_j, -gentime);
#ifdef _HOUSEPI_DISCRETE_TIME
			cumul_mtx[j][i] = infec_profile_cumul(indiv_j, -gentime - 1);
#else
			cumul_mtx[j][i] = infec_profile_cumul(indiv_j, -gentime);
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

	/// Introduction rule: earliest infection, ties to smallest hosuehold index.
	static void helper_commit_intro(MemberIndex &i0, double &t0, MemberIndex i, double t)
	{
		if (t < t0 || (t == t0 && i < i0))
		{
			t0 = t;
			i0 = i;
		}
	}

	/// Store cumulative and instantaneous infectivity profiles for household.
	/// WARNING: unpredictable behavior if i is not an infector.
	/// WARNING: these function must not call any mutable parameter.
	void update_house_transmission_cache(const Household &house) override final
	{
		auto &inst_mtx = cache_double_matrix(house, HHdoublemtx::inst);
		auto &cumul_mtx = cache_double_matrix(house, HHdoublemtx::cumul);
		MemberIndex &i0 = cache_member(house, HHmember::intro_case);
		double new_t0 = 100000;

		const MemberIndex n_infec = house.n_infected_status();
		const auto &infected_members = house.infected_members();

		for (MemberIndex a = 0; a < n_infec; ++a)
		{
			const MemberIndex i = infected_members[a];
			const Individual &indiv_i = house.at(i);
			const double t_infec_i = aug(Augm::infec_time, indiv_i);

			///////////////////////////////
			//////  Escaped Members  //////
			///////////////////////////////

			for (MemberIndex s : house.escaped_members())
			{
				const Individual &indiv_s = house.at(s);
				const double exposure = follow_up_end(indiv_s) - t_infec_i;
				cumul_mtx[i][s] = infec_profile_cumul(indiv_i, exposure);
				cumul_mtx[s][i] = 0.0;
			}

			///////////////////////////////////
			//////// Infectious Members ///////
			///////////////////////////////////

			for (MemberIndex b = a + 1; b < n_infec; ++b)
			{
				const MemberIndex j = infected_members[b];
				const Individual &indiv_j = house.at(j);
				const double t_infec_j = aug(Augm::infec_time, indiv_j);

				helper_update_pairwise_cache(inst_mtx, cumul_mtx,
											 indiv_i, indiv_j,
											 t_infec_i, t_infec_j);
			}

			// update household introduction case
			helper_commit_intro(i0, new_t0, i, t_infec_i);
		}
	}

	/// Store cumulative and instantaneous infectivity profiles according to an infector.
	void update_member_transmission_cache(const Household &house,
										  const Individual &indiv_i) override final
	{
		auto &inst_mtx = cache_double_matrix(house, HHdoublemtx::inst);
		auto &cumul_mtx = cache_double_matrix(house, HHdoublemtx::cumul);
		MemberIndex &i0 = cache_member(house, HHmember::intro_case);

		const MemberIndex i = Index::member(indiv_i);
		const double t_infec_i = aug(Augm::infec_time, indiv_i);

		///////////////////////////////
		//////  Escaped Members  //////
		///////////////////////////////

		for (MemberIndex s : house.escaped_members())
		{
			const Individual &indiv_s = house.at(s);
			double exposure = follow_up_end(indiv_s) - t_infec_i;
			cumul_mtx[i][s] = infec_profile_cumul(indiv_i, exposure);
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

		for (MemberIndex j : house.infected_members())
		{
			if (j == i)
				continue;

			const Individual &indiv_j = house.at(j);
			const double t_infec_j = aug(Augm::infec_time, indiv_j);

			helper_update_pairwise_cache(inst_mtx, cumul_mtx,
										 indiv_i, indiv_j,
										 t_infec_i, t_infec_j);

			// update household introduction case
			helper_commit_intro(i0, new_t0, j, t_infec_j);
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
		const double Dt = real_time_step();
#endif
		/////////////////////////////////////////////////////////////////
		/////////////  Generation 0 : by single introduction  ///////////
		/////////////////////////////////////////////////////////////////

		const MemberIndex i0 = helper_rand_intro(house);
		Individual &indiv_0 = house.at(i0);

		/////////////////////////////////////////////////////////////////
		//////  Generation >=1 : infection by previous generation  //////
		/////////////////////////////////////////////////////////////////

		const double alpha = par(Param::alpha);
		const double beta = transmission_rate(house.size());

		const double t0 = aug(Augm::infec_time, indiv_0);
		for (double t = t0 + Dt; house.n_escaped_status() > 0; t += Dt)
		{
			bool any_at_risk = false;
			for (Individual &indiv_i : house)
			{
				if (indiv_i.is_infected_status())
					continue;
				if (t > follow_up_end(indiv_i))
					continue;
				any_at_risk = true;

				if (exponential() < roi_inst_simulation(house, indiv_i, t, alpha, beta) * Dt)
				{
					set_infected(house, indiv_i);
					set_augm_value(Augm::infec_time, indiv_i, t);
					init_missing_augm(indiv_i);
				}
			}
			if (!any_at_risk)
				break;
		}
	}

	/// Draw introduction uniformly among all members.
	MemberIndex helper_rand_intro(Household &house)
	{
		const MemberIndex i0 = uniform() * house.size();
		Individual &indiv_0 = house.at(i0);
		set_infected(house, indiv_0);
		init_missing_augm(indiv_0);

		// set all others members to escaped status
		for (Individual &indiv : house)
			if (indiv != indiv_0)
				set_escaped(house, indiv);

		return i0;
	}

	/// Instantaneous risk of Infection for individual i. Used in simulation.
	double roi_inst_simulation(
		const Household &house,
		const Individual &indiv_i,
		double t_infec_i,
		double alpha,
		double beta) const
	{
		MemberIndex i = Index::member(indiv_i);
		double roi = 0.; // instantaneous risk of infection

		// j : infector, i : infectee
		for (MemberIndex j : house.infected_members())
		{
			const Individual &indiv_j = house.at(j);
			const double gentime = t_infec_i - aug(Augm::infec_time, indiv_j);
			roi += group_risk(GroupRisk::infectivity, indiv_j) *
				   contact_risk(house, j, i) *
				   infec_profile_inst(indiv_j, gentime);
		}

		return alpha * group_risk(GroupRisk::community, indiv_i) +
			   beta * group_risk(GroupRisk::susceptibility, indiv_i) * roi;
	}
};
HOUSEPI_CLOSE_MODULE_NAMESPACE()
