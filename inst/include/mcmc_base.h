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
#include "risks.h"
#include "infec_profile.h"
#include "augmented.h"
#include "posterior.h"
#include "translator.h"

//////////////////////////////  MCMCBase  //////////////////////////////

namespace housepi
{
	class MCMCBase
	{
	public:
		Setup setup;
		LogProba tot_obs_log_lik;	 /// Observation log-likelihood.
		LogProba tot_transm_log_lik; /// Transmission log-likelihood.
		LogProba tot_latent_log_lik; /// Latent-model log-likelihood.
		Posterior posterior;

		std::shared_ptr<DataStore> pDataStore; /// Data storage.

		//////  Key holders  //////

		Parameter param;
		Risks risks;
		InfecProfile infec_profile;
		Augmented augm;
		Translator translator;

		//////  Constructor  //////

		MCMCBase();
		virtual ~MCMCBase();

		/// Initialize every key holder.
		void init_key_holder(const std::vector<std::string> &data_wide_keys,
							 const std::vector<std::string> &hh_member_keys,
							 const std::vector<std::string> &hh_member_vec_keys,
							 const std::vector<std::string> &hh_double_vec_keys,
							 const std::vector<std::string> &hh_double_mtx_keys,
							 const std::vector<std::string> &data_long_keys,
							 const std::vector<std::string> &param_keys,
							 const std::vector<std::string> &group_risk_keys,
							 const std::vector<std::string> &latent_keys);

		//////  Data access  //////

		/// Attach data store `ds`.
		void set_data_store(std::shared_ptr<DataStore> ds);

		DataWide &data_wide() noexcept { return pDataStore->data_wide; }
		const DataWide &data_wide() const noexcept { return pDataStore->data_wide; }
		DataLong &data_long() noexcept { return pDataStore->data_long; }
		const DataLong &data_long() const noexcept { return pDataStore->data_long; }

		/// Initialize equivalence-class compression.
		/// - latent definitions required, not imported yet: one class per household.
		void init_equivalence_class_compression();

		/// Throw if the generic profile transmits at generation time 0.
		/// - status-specific profile: checked per member, with its variable.
		void assert_valid_infectivity_profile() const;

		//////  Getter  //////

		/// Latent `a` of member `m`.
		double aug(size_t a, const Individual &m) const { return augm.at(a, m).value(); }
		/// Datum `d` of member `m`.
		double data(size_t d, const Individual &m) const noexcept { return data_wide().data(d, m); }

		/// Data column names shared across simulations.
		std::vector<std::string> get_static_data_names() const;
		/// Column names of the per-iteration shallow list.
		std::vector<std::string> get_col_names() const;
		/// Per-household row offsets, so parallel simulations write disjoint rows.
		std::vector<size_t> get_row_offsets() const;
		/// Sum the observation log-likelihood over infected members of `house`.
		double house_obs_log_lik(const Household &house) const;
		/// Sum the latent log-likelihood over infected members of `house`.
		double house_latent_log_lik(const Household &house) const;

	protected:
		//////  Log-likelihood  //////

		/// Observation-model log-likelihood of `m`.
		LogProba &obs_log_lik_(const Individual &m) noexcept
		{
			return augm.state.of(m).obs_log_lik;
		}
		const LogProba &obs_log_lik_(const Individual &m) const noexcept
		{
			return augm.state.of(m).obs_log_lik;
		}

		/// Latent-model log-likelihood of `m`.
		LogProba &latent_log_lik_(const Individual &m) noexcept
		{
			return augm.state.of(m).latent_log_lik;
		}
		const LogProba &latent_log_lik_(const Individual &m) const noexcept
		{
			return augm.state.of(m).latent_log_lik;
		}

		//////  Metropolis-Hastings  //////

		/// Assert configuration and initialize Monte Carlo.
		void init_monte_carlo_();
		/// Initialize one Markov chain.
		void init_monte_carlo_chain_();

		/// Recover previous latent log-likelihood of `m`, after a rejected move.
		void recover_latent_log_lik_(const Individual &m) noexcept
		{
			latent_log_lik_(m).recover();
		}
		/// Recover previous observation log-likelihood of `m`, after a rejected move.
		void recover_obs_log_lik_(const Individual &m) noexcept
		{
			obs_log_lik_(m).recover();
		}
		/// Recover previous transmission log-likelihood of `house`, after a rejected move.
		void recover_transm_log_lik_(const Household &house) noexcept
		{
			data_wide().transm_log_lik(house).recover();
		}
		/// Recover previous total log-likelihood.
		/// @param with_cache Whether to recover cached values too.
		void recover_total_log_lik_(bool with_cache,
									bool with_latent_ll = true,
									bool with_obs_ll = true);
	};
}
