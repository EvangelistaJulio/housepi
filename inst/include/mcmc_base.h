/* housepi - Epidemiological Analysis of Household Outbreaks.
 *
 * Copyright © 2024-2026 Institut Pasteur (Paris), CNRS and INSERM.
 * Distributed under the terms of the GNU General Public License (GPLv3).
 * See the LICENSE file for details.
 */

#pragma once

#include <memory>
#include <string>
#include <vector>

#include "tools/logproba.h"
#include "tools/nan.h"

#include "setup.h"
#include "data_store.h"
#include "parameter.h"
#include "group_rate.h"
#include "contact_rate.h"
#include "infec_profile.h"
#include "augmented.h"
#include "posterior.h"
#include "translator.h"

/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////  *** MCMCBase : class ***  //////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

namespace housepi
{
	class MCMCBase
	{
	public:
		Setup setup;
		LogProba tot_obs_log_lik;	 /// observation log-likelihood
		LogProba tot_transm_log_lik; /// transmission log-likelihood
		LogProba tot_augm_log_prior; /// augmented-variable log-priors
		Posterior posterior;

		/// Pointer to data storage.
		std::shared_ptr<DataStore> pDataStore;

		/// Key holders.
		Parameter param;
		GroupRate group_rate;
		ContactRate contact_rate;
		InfecProfile infec_profile;
		Augmented augm;
		Translator translator;

		///////////////////////////
		//////  CONSTRUCTOR  //////
		///////////////////////////

		MCMCBase();
		virtual ~MCMCBase();

		/// Initialize all key-holders.
		void init_key_holder(
			const std::vector<std::string> &data_wide_keys,
			const std::vector<std::string> &hh_size_keys,
			const std::vector<std::string> &hh_size_vec_keys,
			const std::vector<std::string> &hh_double_vec_keys,
			const std::vector<std::string> &hh_double_mtx_keys,
			const std::vector<std::string> &data_long_keys,
			const std::vector<std::string> &param_keys,
			const std::vector<std::string> &group_rate_keys,
			const std::vector<std::string> &augm_keys);

		///////////////////////////
		//////  DATA ACCESS  //////
		///////////////////////////

		/// Associate data storage.
		void set_data_store(std::shared_ptr<DataStore> ds);

		DataWide &data_wide() noexcept { return pDataStore->data_wide; }
		const DataWide &data_wide() const noexcept { return pDataStore->data_wide; }
		DataLong &data_long() noexcept { return pDataStore->data_long; }
		const DataLong &data_long() const noexcept { return pDataStore->data_long; }

		/// Initialize equivalence class compression.
		/// NOTE: assume unique if augmented definitions required, but not imported yet.
		void init_equivalence_class_compression();

		////////////////////////////
		//////  AUGMENTATION  //////
		////////////////////////////

		double aug(size_t a, const Individual &m) const { return augm.at(a, m).value(); }
		double data(size_t d, const Individual &m) const noexcept { return data_wide().data(d, m); }
		Variable &get_augm_varia(size_t a, const Individual &m) noexcept { return augm.at(a, m); }
		const Variable &get_augm_varia(size_t a, const Individual &m) const noexcept { return augm.at(a, m); }

		//////////////////////
		//////  GETTER  //////
		//////////////////////

		/// Whether two household members belong to same partition according to rates.
		bool is_same_rate_partition(const Household &house, size_t i, size_t j) const;

		/// Column names of data shared across simulations.
		std::vector<std::string> get_static_data_names() const;
		/// Column name vector for the per-iteration shallow list.
		std::vector<std::string> get_col_names() const;
		// Pre-compute per-household row offsets. Used to fix OpenMP race in simulations.
		std::vector<size_t> get_row_offsets() const;
		/// Sum per-infected-member observation log-likelihhood.
		double house_log_obs(const Household &house) const;

	protected:
		///////////////////////////////////
		//////  METROPOLIS-HASTINGS  //////
		///////////////////////////////////

		/// Assert and initialize configuration to proceed with Monte Carlo.
		void init_monte_carlo_();
		/// Intialize a single Markov chain for Monte Carlo.
		void init_monte_carlo_chain_();

		/// Recover the old total log-likelihood, augmented log-priors and observation.
		/// NOTE: called after a refused parameter move.
		void recover_total_log_density_(bool is_cache_depen,
										bool with_augm_prior = true,
										bool with_obs = true);
	};
}
