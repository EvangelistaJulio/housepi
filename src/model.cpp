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

#include "housepi.h"
#include "tools/messages.h"
#include "mcmc.h"

HOUSEPI_OPEN_NAMED_MODULE_NAMESPACE()
HOUSEPI_DECLARE_EMPTY_KEYS(AugmKey)
HOUSEPI_DECLARE_EMPTY_KEYS(DataKey)
HOUSEPI_DECLARE_EMPTY_KEYS(LongKey)
HOUSEPI_DECLARE_EMPTY_KEYS(GroupRateKey)
HOUSEPI_DECLARE_EMPTY_KEYS(HHsizeKey)
HOUSEPI_DECLARE_EMPTY_KEYS(HHsizevecKey)
HOUSEPI_DECLARE_EMPTY_KEYS(HHdoublevecKey)
HOUSEPI_DECLARE_EMPTY_KEYS(HHdoublemtxKey)
HOUSEPI_DECLARE_EMPTY_KEYS(ParamKey)

class Model final : public MCMC
{
public:
	HOUSEPI_MODEL_CONSTRUCTOR()

	/// Raise error when "blank" model.
	void init_model_for_evaluation() override final { hsp::msg::warning(
		"Computed log-likelihood with 'blank' model."); }

	/// Raise error when "blank" model.
	void init_model_for_simulation() override final { hsp::msg::warning(
		"Performed simulation with 'blank' model."); }

	/// Unitary likelihood.
	double compute_house_obs_log_lik(Household &) override final { return 0.; }

	/// Augmented variable initialization.
	void simulate_single_house(Household &house) override final
	{
		for (size_t i : house.infected_hhindexes())
			init_missing_augm(house.member(i));
	}
};
HOUSEPI_CLOSE_NAMED_MODULE_NAMESPACE()
HOUSEPI_REGISTER_RCPP_MODULE_INTERFACE()
