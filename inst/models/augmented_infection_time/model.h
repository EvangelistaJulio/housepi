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

#ifndef HOUSEPI__MODEL__H
#define HOUSEPI__MODEL__H

#include "mcmc.h"

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////// Model : class //////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

class Model final : public MCMC
{
public:
	Model();

	/// Initialize infectivity profile.
	void init_advanced_setup() override final;
	/// Initialize household introduction case evaluation (fixed for single infector).
	void init_house_for_evaluation(Household &house) override final;

	void helper_update_pairwise_cache(
		std::vector<std::vector<double>> &inst_mtx,
		std::vector<std::vector<double>> &cumul_mtx,
		const Individual &member_i,
		const Individual &member_k,
		double t_infec_i,
		double t_infec_k);

	/// WARNING: unpredictable behavior if i is not an infector.
	/// WARNING: these function must not call any mutable parameter.
	/// Store cumulative and instantaneous infectivity profiles for household.
	void update_house_cache_dependecies(Household &house) override final;
	/// Store cumulative and instantaneous infectivity profiles according to an infector.
	void update_member_cache_dependecies(Household &house, Individual &member_i) override final;

	/// Compute the log-likelihood within household.
	double compute_house_obs_log_lik(Household &house) override final
	{
		return log_lik_susceptible(house) + log_lik_infectious(house);
	}
	/// Transmission rate.
	double transmission_rate(size_t house_size) const;
	/// Susceptible members contribution to log-likelihood.
	double log_lik_susceptible(const Household &house) const;
	/// Infectious members contribution to log-likelihood.
	double log_lik_infectious(const Household &house) const;

	/// Poisson point process to simulate outbreak in a single household.
	void simulate_single_house(Household &house) override final;
	/// Instantaneous risk of Infection for individual i. Used in simulation.
	double roi_inst_simulation(Household &house, Individual &member_i, double t_infec_i) const;
};

#endif // HOUSEPI__MODEL__H