/* housepi - Epidemiological Analysis of Household Outbreaks.
 *
 * Copyright © 2026 Institut Pasteur (Paris), CNRS and INSERM.
 * Distributed under the terms of the GNU General Public License (GPLv3).
 * See the LICENSE file for details.
 */

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////  * *** * CHAIN BINOMIAL HETEROGENEOUS * *** *  /////////////////////////
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

	/// Initialize the LookUp Table.
	void init_model_for_evaluation() override final;
	/// Initialize internal order for house escape matrix
	/// ordered members : 1) infectious, 2) susceptibles
	/// NOTE: this works only if size = n_infected_status + n_susceptible_status
	void init_house_for_evaluation(Household &house) override final;

	/// Compute the log-likelihood for a single household.
	double compute_house_obs_log_lik(Household &house) override final;
	/// Risk of infection from the community.
	double community_risk() const;
	/// Risk of infection from another household member.
	double person_to_person_risk() const;
	/// The person-to-person escape probabilities stored inside a matrix.
	/// Important: it must be re-ordered to have seropositives first.
	void fill_house_escape_matrix(Household &house);

	/// Simulate the infection within a household.
	void simulate_single_house(Household &house) override final;
};

#endif // HOUSEPI__MODEL__H