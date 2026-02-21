/* housepi - Epidemiological Analysis of Household Outbreaks.
 *
 * Copyright © 2026 Institut Pasteur (Paris), CNRS and INSERM.
 * Distributed under the terms of the GNU General Public License (GPLv3).
 * See the LICENSE file for details.
 */

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////  * *** * BLANK MODEL * *** *  /////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

#include "housepi_model_setup.h"
#include "tools/messages.h"
#include "model.h"

#ifndef HOUSEPI_MODULE_NAME
#define HOUSEPI_MODULE_NAME housepiblank
#endif

HOUSEPI_DECLARE_EMPTY_KEYS(AugmKey)
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

/// Raise error when "blank" model.
void Model::init_model_for_evaluation() { hsp::msg::warning(
	"Computed log-likelihood with 'blank' model."); }

/// Raise error when "blank" model.
void Model::init_model_for_simulation() { hsp::msg::warning(
	"Performed simulation with 'blank' model."); }

/// Unitary likelihood.
double Model::compute_house_obs_log_lik(Household &) { return 0.; }

/// Augmented variable initialization.
void Model::simulate_single_house(Household &house)
{
	for (size_t i : house.infected_hhindexes())
		init_missing_augm(house.member(i));
}
