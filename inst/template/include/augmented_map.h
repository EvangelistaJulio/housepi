/* housepi - Epidemiological Analysis of Household Outbreaks.
 *
 * Copyright © 2026 Institut Pasteur (Paris), CNRS and INSERM.
 * Distributed under the terms of the GNU General Public License (GPLv3).
 * See the LICENSE file for details.
 */

#ifndef HOUSEPI__AUGMENTED__MAP__H
#define HOUSEPI__AUGMENTED__MAP__H

#include "tools/operators.h"
#include "tools/parser.h"
#include "tools/types.h"
#include "variable.h"

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////// *** AugmentedMap : class *** ////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

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
		Variable &this_augm = safe_at(augm_map_, a, "augmented variable");
		this_augm.init(value, is_const, func_name, is_prim);
		this_augm.init_distr(model_name);
	}

	/// Set numeric value of augmented variable. If not allowed, throw exception.
	void set_augm_value(const std::string &a, double x) { safe_at(augm_map_, a, "augmented").set(x); }

	///////////////////////////////////////////////////////////////////////////////////////////////
	///////////////////////////////////  AugmentedMap : getter  ///////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////////////////

	virtual std::string err_name() const { return ""; }
	double aug(const std::string &a) const { return safe_at(augm_map_, a, "augmented").get(); }
	Variable &get_augm_variable(const std::string &a) { return safe_at(augm_map_, a, "augmented"); }
	const Variable &get_augm_variable(const std::string &a) const { return safe_at(augm_map_, a, "augmented"); }
	std::string str_augm_summary() const;

	///////////////////////////////////////////////////////////////////////////////////////////////
	/////////////////////////////////  Individual : log-likelihood  ///////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////////////////

	void init_augm_model_logratios(const u_dict<str_vec> &models_ordered);
	double get_augm_model_log_lik() const { return augm_log_lik_; }
	double materialize_augm_model_log_lik();
	double materialize_augm_model_logratio(const std::string &a);
	void recover_augm_model_log_lik() { augm_log_lik_ = buffer_; };
	void assert_log_lik_not_NaN(const std::string &a = hsp::na::STRING) const
	{
		if (hsp::na::is(augm_log_lik_))
			throw std::runtime_error(
				"Augmented log-likelihood returned NaN for " + err_name() +
				(hsp::na::is(a) ? "" : " after change of " + augm_map_.find(a)->second.err_name()) + ".");
	}

private:
	u_dict<Variable> augm_map_ = {}; /// dictionary of augmented definitions

	double augm_log_lik_ = 0; // current augmented model log-likelihood
	double buffer_ = 0;		  // buffer augmented model log-likelihood

	// container for the logratio lists according to dependencies
	u_dict<std::vector<const Variable *>> model_logratios_ = {};
};

#endif // HOUSEPI__AUGMENTED__MAP__H
