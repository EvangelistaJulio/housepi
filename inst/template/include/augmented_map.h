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

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////// *** AugmentedContainer : class *** /////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/// Class to deal with multiple augmented defintions.
class AugmentedContainer
{
	friend struct Augmented;
	friend struct DataWide;
	friend class MCMCImpl;

public:
	virtual ~AugmentedContainer() = default;

	/// Set numeric value of augmented variable. If not allowed, raise exception.
	void set_augm_value(size_t a, double x) { augm_vec_[a].set(x); }

	double aug(size_t a) const { return augm_vec_[a].get(); }
	Variable &get_augm_variable(size_t a) { return augm_vec_[a]; }
	const Variable &get_augm_variable(size_t a) const { return augm_vec_[a]; }
	std::string str_augm_summary() const;

protected:
	size_t augm_idx_status_;			  /// index for augmented model
	std::vector<Variable> augm_vec_ = {}; /// dictionary of augmented definitions
	hsp::LogLike augm_log_lik_; /// augmented model log-likelihood

	AugmentedContainer() = default;
};
