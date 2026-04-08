/* housepi - Epidemiological Analysis of Household Outbreaks.
 *
 * Copyright © 2026 Institut Pasteur (Paris), CNRS and INSERM.
 * Distributed under the terms of the GNU General Public License (GPLv3).
 * See the LICENSE file for details.
 */

#pragma once

#include <memory>
#include <string>
#include <vector>

#include "tools/atomic.h"
#include "tools/loglike.h"
#include "tools/nan.h"

#include "setup.h"
#include "data_store.h"
#include "parameter.h"
#include "group_rate.h"
#include "contact_rate.h"
#include "infec_profile.h"
#include "augmented.h"
#include "posterior.h"

/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////  *** MCMCBase : class ***  //////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

namespace hsp
{
	class MCMCBase
	{
	public:
		atomic::Error atomic_error;
		LogLike tot_log_lik;
		Posterior posterior;

		std::shared_ptr<DataStore> pDataStore;

		Setup setup;
		Parameter param;
		GroupRate group_rate;
		ContactRate contact_rate;
		InfecProfile infec_profile;
		Augmented augm;

		///////////////////////////
		//////  CONSTRUCTOR  //////
		///////////////////////////

		MCMCBase();
		virtual ~MCMCBase();

		////////////////////////////
		//////  DATA ACCESS  //////
		////////////////////////////

		DataWide &data_wide() { return pDataStore->data_wide; }
		const DataWide &data_wide() const { return pDataStore->data_wide; }
		DataLong &data_long() { return pDataStore->data_long; }
		const DataLong &data_long() const { return pDataStore->data_long; }

		void set_data_store(std::shared_ptr<DataStore> ds);

		void init_registry(const std::vector<std::string> &data_wide_keys,
						   const std::vector<std::string> &hh_size_keys,
						   const std::vector<std::string> &hh_size_vec_keys,
						   const std::vector<std::string> &hh_double_vec_keys,
						   const std::vector<std::string> &hh_double_mtx_keys,
						   const std::vector<std::string> &data_long_keys,
						   const std::vector<std::string> &param_keys,
						   const std::vector<std::string> &group_rate_keys,
						   const std::vector<std::string> &augm_keys);

	protected:
		///////////////////////////////////
		//////  METROPOLIS-HASTINGS  //////
		///////////////////////////////////

		/// Assert setup definitions are ok to proceed with mcmc sampling.
		void init_monte_carlo_();
		/// Intialize a single Markov chain for Monte Carlo.
		void init_monte_carlo_chain_();

		/// Recover the old parameter value.
		void recover_parameter_(size_t p, bool is_cache_depen);
	};
}
