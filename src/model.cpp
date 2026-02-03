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

EMPTY_HOUSEPI_KEYS(AugmKey)
EMPTY_HOUSEPI_KEYS(DataKey)
EMPTY_HOUSEPI_KEYS(LongKey)
EMPTY_HOUSEPI_KEYS(GroupRateKey)
EMPTY_HOUSEPI_KEYS(HHsizeKey)
EMPTY_HOUSEPI_KEYS(HHsizevecKey)
EMPTY_HOUSEPI_KEYS(HHdoublevecKey)
EMPTY_HOUSEPI_KEYS(HHdoublemtxKey)
EMPTY_HOUSEPI_KEYS(ParamKey)
HOUSEPI_MODEL_CONSTRUCTOR()

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
