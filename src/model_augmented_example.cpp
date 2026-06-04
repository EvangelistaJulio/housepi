/* housepi - Epidemiological Analysis of Household Outbreaks.
 *
 * Copyright © 2024-2026 Institut Pasteur (Paris), CNRS and INSERM.
 * Distributed under the terms of the GNU General Public License (GPLv3).
 * See the LICENSE file for details.
 */

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////  * *** * DATA AUGMENTATION EXAMPLE * *** *  ///////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

#include "housepi.h"

#ifndef _HOUSEPI_MODULE
#define _HOUSEPI_MODULE augmented_example
#endif

HOUSEPI_OPEN_MODULE_NAMESPACE()
class Model final : public ModelBase
{
public:
	HOUSEPI_DECLARE_REGISTER_KEYS(AugmKey, proba)
	HOUSEPI_DECLARE_EMPTY_KEYS(DataKey)
	HOUSEPI_DECLARE_EMPTY_KEYS(LongKey)
	HOUSEPI_DECLARE_EMPTY_KEYS(GroupRateKey)
	HOUSEPI_DECLARE_EMPTY_KEYS(HHsizeKey)
	HOUSEPI_DECLARE_EMPTY_KEYS(HHsizevecKey)
	HOUSEPI_DECLARE_EMPTY_KEYS(HHdoublevecKey)
	HOUSEPI_DECLARE_EMPTY_KEYS(HHdoublemtxKey)
	HOUSEPI_DECLARE_EMPTY_KEYS(ParamKey)

	/// Compute transmission log-likelihood within household.
	double compute_house_transm_log_lik(Household &house) override final
	{
		double output = 0.0;
		for (size_t i : house.infected_hhindexes())
			output += std::log(house.member(i).aug(AugmKey::proba));
		return output;
	}

	/// Simulate augmented infectious members augmented variables.
	void simulate_single_house(Household &house) override final
	{
		for (size_t i : house.infected_hhindexes())
			init_missing_augm(house.member(i));
	}
};
HOUSEPI_CLOSE_MODULE_NAMESPACE()
