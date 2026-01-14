/* housepi - Epidemiological Analysis of Household Outbreaks.
 *
 * Copyright © 2026 Institut Pasteur (Paris), CNRS and INSERM.
 * Distributed under the terms of the GNU General Public License (GPLv3).
 * See the LICENSE file for details.
 */

#ifndef HOUSEPI__AUGMENTED__MAP__H
#define HOUSEPI__AUGMENTED__MAP__H

#include <string>
#include <vector>

#include "tools/operators.h"
#include "tools/types.h"

#include "variable.h"

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////// *** AugmentedMap : class *** ////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

namespace hsp::parse
{
	struct Definition;
}

/// Class to deal with multiple augmented defintions.
/// UPDATE: to be extended to Household.
class AugmentedMap
{
public:
	///////////////////////////////////////////////////////////////////////////////////////////////
	/////////////////////////////////  AugmentedMap : constructor  ////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////////////////

	AugmentedMap() {}
	virtual ~AugmentedMap() = default;

	///////////////////////////////////////////////////////////////////////////////////////////////
	///////////////////////////////////  AugmentedMap : setter  ///////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////////////////

	/// Clear augmented definition (also used for allocation).
	void clear_augm(const std::string &a)
	{
		augm_map_[a] = Variable();
		augm_map_[a].set_name(a, "augmented");
	}

	/// Assign definitions for the augmented variable model.
	/// NOTE: a Variable must already be allocated with clear_augm.
	void init_augm(const std::string &a,
				   double value,
				   bool is_const,
				   const hsp::parse::Definition &func_name,
				   bool is_prim,
				   const hsp::parse::Definition &model_name)
	{
		auto &this_augm = hsp::safe_at(augm_map_, a, "augmented variable");
		this_augm.init(value, is_const, func_name, is_prim);
		this_augm.init_distr(model_name);
	}

	/// Set numeric value of augmented variable. If not allowed, throw exception.
	void set_augm_value(const std::string &a, double x) { hsp::safe_at(augm_map_, a, "augmented").set(x); }

	///////////////////////////////////////////////////////////////////////////////////////////////
	///////////////////////////////////  AugmentedMap : getter  ///////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////////////////

	virtual std::string err_name() const = 0;
	double aug(const std::string &a) const { return hsp::safe_at(augm_map_, a, "augmented").get(); }
	Variable &get_augm_variable(const std::string &a) { return hsp::safe_at(augm_map_, a, "augmented"); }
	const Variable &get_augm_variable(const std::string &a) const { return hsp::safe_at(augm_map_, a, "augmented"); }
	std::string str_augm_summary() const;

	///////////////////////////////////////////////////////////////////////////////////////////////
	/////////////////////////////////  Individual : log-likelihood  ///////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////////////////

	void init_augm_model_logratios(const hsp::str_table &models_ordered);
	double get_augm_model_log_lik() const { return augm_log_lik_; }
	double materialize_augm_model_log_lik();
	double materialize_augm_model_logratio(const std::string &a);
	void recover_augm_model_log_lik() { augm_log_lik_ = buffer_; };
	void assert_log_lik_not_NaN(const std::string &a = "") const;

private:
	hsp::dict<Variable> augm_map_ = {}; /// dictionary of augmented definitions

	double augm_log_lik_ = 0; // current augmented model log-likelihood
	double buffer_ = 0;		  // buffer augmented model log-likelihood

	// container for the logratio lists according to dependencies
	hsp::dict<std::vector<const Variable *>> model_logratios_ = {};
};

#endif // HOUSEPI__AUGMENTED__MAP__H
