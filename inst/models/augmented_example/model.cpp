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

#include "housepi_model_setup.h"
#include "model.h"

REGISTER_HOUSEPI_KEYS(AugmKey, proba)
EMPTY_HOUSEPI_KEYS(DataKey)
EMPTY_HOUSEPI_KEYS(LongKey)
EMPTY_HOUSEPI_KEYS(GroupRateKey)
EMPTY_HOUSEPI_KEYS(HHsizeKey)
EMPTY_HOUSEPI_KEYS(HHsizevecKey)
EMPTY_HOUSEPI_KEYS(HHdoublevecKey)
EMPTY_HOUSEPI_KEYS(HHdoublemtxKey)
EMPTY_HOUSEPI_KEYS(ParamKey)
HOUSEPI_MODEL_CONSTRUCTOR()

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////  Model : compute single house log-likelihood  /////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

double Model::compute_house_obs_log_lik(Household &house)
{
	double output = 0.0;
	for (size_t i : house.infected_hhindexes())
		output += std::log(house.member(i).aug(AugmKey::proba));
	return output;
}

void Model::simulate_single_house(Household &house)
{
	for (size_t i : house.infected_hhindexes())
		init_missing_augm(house.member(i));
}
