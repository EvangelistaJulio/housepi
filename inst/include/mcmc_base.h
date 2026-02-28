/* housepi - Epidemiological Analysis of Household Outbreaks.
 *
 * Copyright © 2026 Institut Pasteur (Paris), CNRS and INSERM.
 * Distributed under the terms of the GNU General Public License (GPLv3).
 * See the LICENSE file for details.
 */

#pragma once

#include <string>
#include <vector>

#include "tools/atomic.h"
#include "tools/loglike.h"
#include "tools/nan.h"

#include "setup.h"
#include "data_wide.h"
#include "data_long.h"
#include "parameter.h"
#include "group_rate.h"
#include "infec_profile.h"
#include "augmented.h"

/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////  *** MCMCBase : class ***  //////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

namespace hsp
{
	struct TotalLogLike : public LogLike
	{
	public:
		std::vector<double> smpl_values = {}; /// total log-likelihood samples

		void init_draws(unsigned long sample_size)
		{
			count_ = 0;
			smpl_values.clear();
			smpl_values.assign(sample_size, na::DOUBLE);
		}
		void store_draws() { smpl_values[count_++] = value_; }
		unsigned long count() const noexcept { return count_; }

	private:
		unsigned long count_ = 0; /// count of stored samples
	};

	class MCMCBase
	{
	public:
		atomic::Error atomic_error;

		Setup setup;
		DataWide data_wide;
		DataLong data_long;
		Parameter param;
		GroupRate group_rate;
		InfecProfile infec_profile;
		Augmented augm;

		TotalLogLike tot_log_lik;

		///////////////////////////
		//////  CONSTRUCTOR  //////
		///////////////////////////

		MCMCBase();
		virtual ~MCMCBase();

		void init_registry(const std::vector<std::string> &data_wide_keys,
						   const std::vector<std::string> &hh_size_keys,
						   const std::vector<std::string> &hh_size_vec_keys,
						   const std::vector<std::string> &hh_double_vec_keys,
						   const std::vector<std::string> &hh_double_mtx_keys,
						   const std::vector<std::string> &data_long_keys,
						   const std::vector<std::string> &param_keys,
						   const std::vector<std::string> &group_rate_keys,
						   const std::vector<std::string> &augm_keys);

		/// Auxiliary method to initialize simulations.
		/// UPDATE: once I allow for multi-treading and a flag other than n_chains = 0, I can drop this.
		void helper_init_simulation()
		{
			data_wide.assert_non_empty(); // check if data is empty
			setup.n_chains = 0;			  // use 0 chains as a state variable
			setup.n_threads = 1;		  // impose no parallelization
			setup.set_n_threads();		  // change accordingly
		}

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

		/// Init posterior sampling.
		//// UPDATE: I would like to store also power posterior sample.
		/// NOTE: it will depend on chosen rungs in that case.
		void init_posterior_sample_();
		/// Store current Bayesian posterior sample.
		void store_posterior_sample_();
	};
}
