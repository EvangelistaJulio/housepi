/* housepi - Epidemiological Analysis of Household Outbreaks.
 *
 * Copyright © 2026 Institut Pasteur (Paris), CNRS and INSERM.
 * Distributed under the terms of the GNU General Public License (GPLv3).
 * See the LICENSE file for details.
 */

#pragma once

#include <string>
#include <vector>

#include "tools/keys.h"
#include "variable.h"

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////// *** AugmentedContainer : class *** /////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/// Class to deal with multiple augmented defintions.
class AugmentedContainer
{
	friend struct Augmented;

public:
	virtual ~AugmentedContainer() = default;

	/// Set numeric value of augmented variable. If not allowed, raise exception.
	void set_augm_value(AugmKey a, double x) { augm_vec_[static_cast<size_t>(a)].set(x); }

	virtual std::string err_name() const = 0;
	double aug(AugmKey a) const { return augm_vec_[static_cast<size_t>(a)].get(); }
	Variable &get_augm_variable(size_t a) { return augm_vec_[a]; }
	const Variable &get_augm_variable(size_t a) const { return augm_vec_[a]; }
	std::string str_augm_summary() const;

	double get_augm_model_log_lik() const { return augm_log_lik_; }
	void recover_augm_model_log_lik() { augm_log_lik_ = buffer_; };

protected:
	size_t augm_idx_status_;			  /// index for augmented model
	std::vector<Variable> augm_vec_ = {}; /// dictionary of augmented definitions
	double augm_log_lik_ = 0;			  /// current augmented model log-likelihood
	double buffer_ = 0;					  /// buffer augmented model log-likelihood

	AugmentedContainer() = default;
};
