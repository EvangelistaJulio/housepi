/* housepi - Epidemiological Analysis of Household Outbreaks.
 *
 * Copyright © 2024-2026 Institut Pasteur (Paris), CNRS and INSERM.
 * Distributed under the terms of the GNU General Public License (GPLv3).
 * See the LICENSE file for details.
 */

//////////////////////////////  test  //////////////////////////////

#include "housepi.h"

#ifndef _HOUSEPI_MODULE
#define _HOUSEPI_MODULE test
#endif

HOUSEPI_OPEN_MODULE_NAMESPACE()
/// Blank model: unit likelihood, latents initialized in simulation.
class Model final : public ModelBase
{
public:
	HOUSEPI_DECLARE_EMPTY_KEYS(Latent)
	HOUSEPI_DECLARE_EMPTY_KEYS(Data)
	HOUSEPI_DECLARE_EMPTY_KEYS(Long)
	HOUSEPI_DECLARE_EMPTY_KEYS(GroupRisk)
	HOUSEPI_DECLARE_EMPTY_KEYS(HHmember)
	HOUSEPI_DECLARE_EMPTY_KEYS(HHmembervec)
	HOUSEPI_DECLARE_EMPTY_KEYS(HHdoublevec)
	HOUSEPI_DECLARE_EMPTY_KEYS(HHdoublemtx)
	HOUSEPI_DECLARE_EMPTY_KEYS(Param)

	/// Return a null log-likelihood: unit likelihood.
	double compute_house_transm_log_lik(const Household &) override final { return 0.; }

	/// Initialize the latent variables of infected members.
	void simulate_single_house(Household &house) override final
	{
		for (MemberIdx i = 0; i < house.size(); ++i)
		{
			Individual &indiv = house.at(i);
			if (indiv.is_infected_status())
				init_missing_latent(indiv);
		}
	}
};
HOUSEPI_CLOSE_MODULE_NAMESPACE()
