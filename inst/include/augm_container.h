/* housepi - Epidemiological Analysis of Household Outbreaks.
 *
 * Copyright © 2026 Institut Pasteur (Paris), CNRS and INSERM.
 * Distributed under the terms of the GNU General Public License (GPLv3).
 * See the LICENSE file for details.
 */

#pragma once

#include <string>
#include <vector>

#include "tools/loglike.h"

#include "variable.h"

namespace hsp::status
{
	constexpr int INFECTED = 1;	   /// Infected status (>=)
	constexpr int SUSCEPTIBLE = 0; /// Susceptible status (==)
	constexpr int MISSING = -1;	   /// Missing status (<=)

	/// Wheter individual status is infected.
	inline bool is_infected(int s) { return s >= INFECTED; }
	/// Wheter individual status is susceptible.
	inline bool is_susceptible(int s) { return s == SUSCEPTIBLE; }
	/// Wheter individual status is missing.
	inline bool is_missing(int s) { return s <= MISSING; }
}

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////// *** AugmContainer : class *** ///////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

class AugmContainer
{
	friend struct Augmented;
	friend struct DataWide;
	friend class MCMCBase;
	friend class MCMC;

public:
	///////////////////////////
	//////  CONSTRUCTOR  //////
	///////////////////////////

	virtual ~AugmContainer() = default;

	//////////////////////
	//////  SETTER  //////
	//////////////////////

	/// Unique identifier as provided in data.
	const std::string &id() const noexcept { return id_; }
	/// Associated status.
	int status() const noexcept { return status_; }

	/// Set numeric value of augmented variable. If not allowed, raise exception.
	void set_augm_value(size_t a, double x) { augm_vec_[a].set(x); }

	//////////////////////
	//////  GETTER  //////
	//////////////////////

	double aug(size_t a) const { return augm_vec_[a].value(); }
	Variable &get_augm_variable(size_t a) { return augm_vec_[a]; }
	const Variable &get_augm_variable(size_t a) const { return augm_vec_[a]; }

	/////////////////////
	//////  FLAGS  //////
	/////////////////////

	bool is_infected_status() const noexcept { return hsp::status::is_infected(status_); }
	bool is_healthy_status() const noexcept { return hsp::status::is_susceptible(status_); }
	bool is_missing_status() const noexcept { return hsp::status::is_missing(status_); }

protected:
	std::string id_ = "";				  /// user id
	int status_ = hsp::status::MISSING;	  /// status of the individual
	size_t augm_idx_;					  /// index of status specific augmented definition
	std::vector<Variable> augm_vec_ = {}; /// dictionary of augmented definitions
	hsp::LogLike augm_log_lik_;			  /// augmented model log-likelihood

	AugmContainer() = default;
};
