/* housepi - Epidemiological Analysis of Household Outbreaks.
 *
 * Copyright © 2026 Institut Pasteur (Paris), CNRS and INSERM.
 * Distributed under the terms of the GNU General Public License (GPLv3).
 * See the LICENSE file for details.
 */

#pragma once

#include <string>
#include <string_view>
#include <vector>

#include "attribute.h"
#include "individual.h"

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////  *** Household : class ***  ////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/// Household class.
class Household
{
	friend struct DataWide;

public:
	///////////////////////////
	//////  CONSTRUCTOR  //////
	///////////////////////////

	Household(const std::string &house_id = "", size_t index = SIZE_MAX)
	{
		house_id_ = house_id;
		index_ = index;
		multiplicity_class_ = index;
	}
	~Household() = default;

	//////////////////////
	//////  SETTER  //////
	//////////////////////

	/// Assign a multiplicity to the household.
	void set_multiplicity(unsigned long m, size_t class_idx)
	{
		house_multiplicity_ = m;
		multiplicity_class_ = class_idx;
	}
	/// Emplace a new individual in the household.
	void init_member(size_t index, const std::string &id,
					 int status, const std::vector<double> &data);

	void set_susceptible_status(Individual &member) { helper_set_status_(member, hsp::status::SUSCEPTIBLE); }
	void set_infected_status(Individual &member) { helper_set_status_(member, hsp::status::INFECTED); }
	void set_infected_status(Individual &member, int new_status);
	void init_contact_rates(size_t from, size_t to, const Attribute &x);

	/////////////////////
	//////  CACHE  //////
	/////////////////////

	// Cache size_t vector
	std::vector<size_t> &cache_size_t_vector(size_t e) { return cache_size_t_vec_[e]; }
	const std::vector<size_t> &cache_size_t_vector(size_t e) const { return cache_size_t_vec_[e]; }

	// Cached double vector
	std::vector<double> &cache_double_vector(size_t e) { return cache_dbl_vec_[e]; }
	const std::vector<double> &cache_double_vector(size_t e) const { return cache_dbl_vec_[e]; }

	// Cache double matrix
	std::vector<std::vector<double>> &cache_double_matrix(size_t e) { return cache_dbl_mtx_[e]; }
	const std::vector<std::vector<double>> &cache_double_matrix(size_t e) const { return cache_dbl_mtx_[e]; }

	//////////////////////////////
	//////  LOG-LIKELIHOOD  //////
	//////////////////////////////

	/// Return value of the household observation log-likelihood.
	double get_obs_log_lik() const noexcept { return obs_log_lik_; }
	/// Store the log-likelihood of household observation.
	void set_obs_log_lik(double x)
	{
		obs_buffer_ = obs_log_lik_;
		obs_log_lik_ = x;
		assert_log_lik_not_NaN();
	}
	/// Assert that current log-likelihood id not a NaN.
	void assert_log_lik_not_NaN() const;

	//////////////////////
	//////  GETTER  //////
	//////////////////////

	Individual &member(size_t i) { return members_.at(i); }
	const Individual &member(size_t i) const { return members_.at(i); }

	unsigned long multiplicity() const noexcept { return house_multiplicity_; }
	size_t multiplicity_class_index() const noexcept { return multiplicity_class_; }
	bool is_multiplicity_head_class() const noexcept { return index_ == multiplicity_class_; }

	const std::vector<size_t> &infected_hhindexes() const noexcept { return infected_list_; }
	const std::vector<size_t> &susceptible_hhindexes() const noexcept { return healthy_list_; }
	const std::vector<size_t> &missing_hhindexes() const noexcept { return missing_list_; }

	size_t n_infected_status() const noexcept { return infected_list_.size(); }
	size_t n_susceptible_status() const noexcept { return healthy_list_.size(); }
	size_t n_missing_status() const noexcept { return missing_list_.size(); }

	double contact_rate(size_t from, size_t to) const { return contact_rates_[from][to].get(); }

	/////////////////////////
	//////  ITERATORS  //////
	/////////////////////////

	auto begin() noexcept { return members_.begin(); }
	auto end() noexcept { return members_.end(); }
	auto begin() const noexcept { return members_.begin(); }
	auto end() const noexcept { return members_.end(); }

	////////////////////
	//////  SHOW  //////
	////////////////////

	const std::string &house_id() const noexcept { return house_id_; }
	std::string err_name() const { return "household '" + house_id_ + "'"; }
	size_t house_index() const noexcept { return index_; }
	size_t size() const noexcept { return members_.size(); }

private:
	std::string house_id_ = "";			   /// User house id
	size_t index_ = -1;					   /// Internal unique index
	unsigned long house_multiplicity_ = 1; /// N. of equivalent houses
	size_t multiplicity_class_ = -1;	   /// Index of multiplicity class

	std::vector<Individual> members_ = {};	 /// Household members
	std::vector<size_t> infected_list_ = {}; /// hh indexes of infected members
	std::vector<size_t> healthy_list_ = {};	 /// hh indexes of susceptible members
	std::vector<size_t> missing_list_ = {};	 /// hh indexes of members with missing status

	std::vector<std::vector<Attribute>> contact_rates_; /// contact rates
	std::vector<std::vector<size_t>> cache_size_t_vec_ = {};
	std::vector<std::vector<double>> cache_dbl_vec_ = {};
	std::vector<std::vector<std::vector<double>>> cache_dbl_mtx_ = {};

	std::vector<std::vector<size_t>> buffer_size_t_vec_ = {};
	std::vector<std::vector<double>> buffer_dbl_vec_ = {};
	std::vector<std::vector<std::vector<double>>> buffer_dbl_mtx_ = {};

	double obs_log_lik_ = 0;				/// current observation log-likelihood
	double obs_buffer_ = 0;					/// buffer observation log-likelihood
	std::vector<double> smpl_log_lik_ = {}; /// samples of household pointwise log-likelihood

	/// Change_status of individual i to new status (without controls).
	void helper_set_status_(Individual &member, int new_status);
};
