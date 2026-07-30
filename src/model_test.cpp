/* housepi - Epidemiological Analysis of Household Outbreaks.
 *
 * Copyright © 2024-2026 Institut Pasteur (Paris), CNRS and INSERM.
 * Distributed under the terms of the GNU General Public License (GPLv3).
 * See the LICENSE file for details.
 */

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////  * *** * BLANK MODEL * *** *  /////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

#include "housepi.h"

#ifndef _HOUSEPI_MODULE
#define _HOUSEPI_MODULE test
#endif

HOUSEPI_OPEN_MODULE_NAMESPACE()
class Model final : public ModelBase
{
public:
	HOUSEPI_DECLARE_EMPTY_KEYS(AugmKey)
	HOUSEPI_DECLARE_EMPTY_KEYS(DataKey)
	HOUSEPI_DECLARE_EMPTY_KEYS(LongKey)
	HOUSEPI_DECLARE_EMPTY_KEYS(GroupRiskKey)
	HOUSEPI_DECLARE_EMPTY_KEYS(HHuintKey)
	HOUSEPI_DECLARE_EMPTY_KEYS(HHuintvecKey)
	HOUSEPI_DECLARE_EMPTY_KEYS(HHdoublevecKey)
	HOUSEPI_DECLARE_EMPTY_KEYS(HHdoublemtxKey)
	HOUSEPI_DECLARE_EMPTY_KEYS(ParamKey)

	/// Compute transmission log-likelihood within household. Unitary likelihood.
	double compute_house_transm_log_lik(Household &) override final { return 0.; }

	/// Augmented variable initialization.
	void simulate_single_house(Household &house) override final
	{
		for (size_t i : house.infected_members())
			init_missing_augm(house.member(i));
	}
};
HOUSEPI_CLOSE_MODULE_NAMESPACE()
