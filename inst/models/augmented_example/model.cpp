/* housepi - Epidemiological Analysis of Household Outbreaks.
 *
 * Copyright © 2026 Institut Pasteur (Paris), CNRS and INSERM.
 * Distributed under the terms of the GNU General Public License (GPLv3).
 * See the LICENSE file for details.
 */

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////  * *** * DATA AUGMENTATION EXAMPLE * *** *  ///////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

#include "tools/keys.h"
#include "tools/keys_helper.h"

#include "individual.h"
#include "household.h"
#include "model.h"

REGISTER_HOUSEPI_KEYS(AugmKey, proba)
REGISTER_HOUSEPI_KEYS(DataKey, )
REGISTER_HOUSEPI_KEYS(GroupRateKey, )
REGISTER_HOUSEPI_KEYS(HHsize_tKey, )
REGISTER_HOUSEPI_KEYS(HHdoubleKey, )
REGISTER_HOUSEPI_KEYS(HHmatrixKey, )
REGISTER_HOUSEPI_KEYS(ParamKey, )

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////  Model : compute single house log-likelihood  /////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/// Compute the log-likelihood within household.
double Model::compute_house_obs_log_lik(Household &house)
{
	double output = 0.0;
	for (const auto &member : house.infected())
		output += std::log(member.aug("proba"));
	return output;
}

/// Augmented variable initialization.
void Model::simulate_single_house(Household &house)
{
	for (auto &member : house.infected())
		init_missing_augm(member);
}
