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

HOUSEPI_DECLARE_REGISTER_KEYS(AugmKey, proba)
HOUSEPI_DECLARE_EMPTY_KEYS(DataKey)
HOUSEPI_DECLARE_EMPTY_KEYS(LongKey)
HOUSEPI_DECLARE_EMPTY_KEYS(GroupRateKey)
HOUSEPI_DECLARE_EMPTY_KEYS(HHsizeKey)
HOUSEPI_DECLARE_EMPTY_KEYS(HHsizevecKey)
HOUSEPI_DECLARE_EMPTY_KEYS(HHdoublevecKey)
HOUSEPI_DECLARE_EMPTY_KEYS(HHdoublemtxKey)
HOUSEPI_DECLARE_EMPTY_KEYS(ParamKey)
HOUSEPI_DEFINE_MODEL_CONSTRUCTOR()
HOUSEPI_DECLARE_RCPP_MODULE(HOUSEPI_MODULE_NAME)

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
