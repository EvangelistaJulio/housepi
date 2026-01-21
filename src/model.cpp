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

REGISTER_HOUSEPI_KEYS(AugmKey, )
REGISTER_HOUSEPI_KEYS(DataKey, )
REGISTER_HOUSEPI_KEYS(GroupRateKey, )
REGISTER_HOUSEPI_KEYS(HHsize_tKey, )
REGISTER_HOUSEPI_KEYS(HHdoubleKey, )
REGISTER_HOUSEPI_KEYS(HHmatrixKey, )
REGISTER_HOUSEPI_KEYS(ParamKey, )

/// Raise error when "blank" model.
void Model::init_model_for_inference() { hsp::msg::warning(
	"Computed log-likelihood with 'blank' model."); }

/// Raise error when "blank" model.
void Model::init_model_for_simulation() { hsp::msg::warning(
	"Performed simulation with 'blank' model."); }

/// Unitary likelihood.
double Model::compute_house_obs_log_lik(Household &) { return 0.; }

/// Augmented variable initialization.
void Model::simulate_single_house(Household &house)
{
	for (Individual &member : house.infected())
		init_missing_augm(member);
}
