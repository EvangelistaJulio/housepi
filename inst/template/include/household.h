/* housepi - Epidemiological Analysis of Household Outbreaks.
 *
 * Copyright © 2026 Institut Pasteur (Paris), CNRS and INSERM.
 * Distributed under the terms of the GNU General Public License (GPLv3).
 * See the LICENSE file for details.
 */

#ifndef HOUSEPI__HOUSEHOLD__H
#define HOUSEPI__HOUSEHOLD__H

#include "tools/types.h"
#include "tools/na_utils.h"
#include "tools/operators.h"
#include "tools/global.h"

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
public:
	///////////////////////////////////////////////////////////////////////////////////////////////
	//////////////////////////////////  Household : constructor  //////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////////////////

	Household(const std::string &house_id = hsp::na::STRING, size_t index = SIZE_MAX)
	{
		house_id_ = house_id;
		index_ = index;
		house_multiplicity_ = 1;
		multiplicity_class_ = index;
		house_size_ = 0;
	}
	~Household() = default;

	///////////////////////////////////////////////////////////////////////////////////////////////
	///////////////////////////////////  Household : iterators  ///////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////////////////

	auto begin() { return members_.begin(); }
	auto end() { return members_.end(); }
	auto begin() const { return members_.begin(); }
	auto end() const { return members_.end(); }

	// T should be Individual for mutable access
	// or const Individual for const access.
	/// NOTE: if Individual is constant, then indexes are immutable as well;
	///		  otherwise, indexes are copied and access is safe.
	template <typename T>
	class View
	{
		// conditional vector reference, if constant or not
		using VectorRef = std::conditional_t<std::is_const<T>::value,
											 const std::vector<Individual> &,
											 std::vector<Individual> &>;							 
		// conditional index, reference or copy
		using IndexRef = std::conditional_t<std::is_const<T>::value,
											   const std::vector<size_t> &,
											   std::vector<size_t>>;

	public:
		VectorRef vec;
		IndexRef indexes; 

		View(VectorRef m, const std::vector<size_t> &i)
			: vec(m), indexes(i) {}

		class iterator
		{
		public:
			iterator(View *v, size_t p) : view_(v), pos_(p) {}

			T &operator*() const { return view_->vec[view_->indexes[pos_]]; }
			T *operator->() const { return &view_->vec[view_->indexes[pos_]]; }

			iterator &operator++()
			{
				++pos_;
				return *this;
			}
			iterator operator++(int)
			{
				iterator tmp = *this;
				++(*this);
				return tmp;
			}

			bool operator==(const iterator &other) const { return pos_ == other.pos_; }
			bool operator!=(const iterator &other) const { return !(*this == other); }

		private:
			View *view_;
			size_t pos_;
		};

		iterator begin() { return iterator(this, 0); }
		iterator end() { return iterator(this, indexes.size()); }
	};

	View<Individual> infected() { return View<Individual>(members_, infected_list_); }
	View<Individual> healthy() { return View<Individual>(members_, healthy_list_); }
	View<Individual> missing() { return View<Individual>(members_, missing_list_); }
	View<const Individual> infected() const { return View<const Individual>(members_, infected_list_); }
	View<const Individual> healthy() const { return View<const Individual>(members_, healthy_list_); }
	View<const Individual> missing() const { return View<const Individual>(members_, missing_list_); }

	const std::vector<size_t> &infected_hhindexes() const { return infected_list_; }
	const std::vector<size_t> &healthy_hhindexes() const { return healthy_list_; }
	const std::vector<size_t> &missing_hhindexes() const { return missing_list_; }

	size_t n_infec_status() const { return infected_list_.size(); }
	size_t n_healthy_status() const { return healthy_list_.size(); }
	size_t n_missing_status() const { return missing_list_.size(); }

	///////////////////////////////////////////////////////////////////////////////////////////////
	////////////////////////////////////  Household : getter  /////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////////////////

	// Individual Getter
	Individual &get_member(size_t i) { return members_.at(i); }
	const Individual &get_member(size_t i) const { return members_.at(i); }

	double contact_rate(size_t from, size_t to) const { return contact_rates_[from][to].get(); }

	// Household Getter
	const std::string &house_id() const { return house_id_; }
	std::string err_name() const { return "household '" + house_id_ + "'"; }
	size_t house_index() const { return index_; }
	size_t size() const { return house_size_; }

	unsigned long multiplicity() const { return house_multiplicity_; }
	size_t multiplicity_class_index() const { return multiplicity_class_; }
	bool is_multiplicity_head_class() const { return (index_ == multiplicity_class_); }

	///////////////////////////////////////////////////////////////////////////////////////////////
	////////////////////////////////////  Household : setter  /////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////////////////

	/// Assign a multiplicity to the household.
	void set_multiplicity(unsigned long m, size_t class_idx)
	{
		house_multiplicity_ = m;
		multiplicity_class_ = class_idx;
	}
	void init_member(size_t index, const std::string &id, int status, const u_dict<double> &data);
	void set_healthy_status(Individual &member) { set_status_(member, hsp::status::SUSCEPTIBLE); }
	void set_infected_status(Individual &member, int new_status = hsp::status::INFECTED)
	{
		if (!hsp::status::is_infected(new_status))
			throw std::invalid_argument(
				"Tried to set as infected status " + std::to_string(new_status) +
				" to " + err_name() + ".\nInfected status must be an integer greater than 0.");
		set_status_(member, new_status);
	}
	void init_contact_rates(size_t from, size_t to, const Attribute &x);

	///////////////////////////////////////////////////////////////////////////////////////////////
	///////////////////////////////  Household : cache attributes  ////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////////////////

	/// Check if any cache quantity is associated to the household.
	bool is_any_cache() { return !(cache_int_vec_.empty() && cache_dbl_vec_.empty() && cache_dbl_mtx_.empty()); }
	/// Deep copies of cache attributes in the buffer
	void buffer_all_cache()
	{
		buffer_cache_int_vec_ = cache_int_vec_;
		buffer_cache_size_t_vec_ = cache_size_t_vec_;
		buffer_cache_dbl_vec_ = cache_dbl_vec_;
		buffer_cache_dbl_mtx_ = cache_dbl_mtx_;
	}

	// Cache integer vector
	void init_cache_integer_vector(const std::string &e) { cache_int_vec_[e] = int_vec(size(), 0); }
	int_vec &cache_integer_vector(const std::string &e) { return safe_at(cache_int_vec_, e, "cached object"); }
	const int_vec &cache_integer_vector(const std::string &e) const { return safe_at(cache_int_vec_, e, "cached object"); }

	// Cache size_t vector
	void init_cache_size_t_vector(const std::string &e) { cache_size_t_vec_[e] = std::vector<size_t>(size(), 0); }
	std::vector<size_t> &cache_size_t_vector(const std::string &e) { return safe_at(cache_size_t_vec_, e, "cached object"); }
	const std::vector<size_t> &cache_size_t_vector(const std::string &e) const { return safe_at(cache_size_t_vec_, e, "cached object"); }

	// Cached double vector
	void init_cache_double_vector(const std::string &e) { cache_dbl_vec_[e] = dbl_vec(size(), 0.); }
	dbl_vec &cache_double_vector(const std::string &e) { return safe_at(cache_dbl_vec_, e, "cached object"); }
	const dbl_vec &cache_double_vector(const std::string &e) const { return safe_at(cache_dbl_vec_, e, "cached object"); }

	// Cache double matrix
	void init_cache_double_matrix(const std::string &e) { cache_dbl_mtx_[e] = matrix<double>(size(), dbl_vec(size(), 0.)); }
	matrix<double> &cache_double_matrix(const std::string &e) { return safe_at(cache_dbl_mtx_, e, "cached object"); }
	const matrix<double> &cache_double_matrix(const std::string &e) const { return safe_at(cache_dbl_mtx_, e, "cached object"); }

	///////////////////////////////////////////////////////////////////////////////////////////////
	/////////////////////////////////  Household : log-likelihood  ////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////////////////

	/// Return value of the household observation log-likelihood.
	double get_obs_log_lik() const { return obs_log_lik_; }
	/// Store the log-likelihood of household observation.
	void set_obs_log_lik(double x)
	{
		obs_buffer_ = obs_log_lik_;
		obs_log_lik_ = x;
		assert_log_lik_not_NaN();
	}
	/// Assert that current log-likelihood id not a NaN.
	void assert_log_lik_not_NaN() const
	{
		if (hsp::na::is(obs_log_lik_))
			throw std::runtime_error(
				"Observation log-likelihood returned NaN.");
	}
	void recover_buffer_values(bool is_cache_depen);

private:
	std::string house_id_ = hsp::na::STRING; /// User house id
	size_t index_ = -1;						 /// Internal unique index
	unsigned short house_size_ = 0;			 /// N. of house members
	unsigned long house_multiplicity_ = 1;	 /// N. of equivalent houses
	size_t multiplicity_class_ = -1;		 /// Index of multiplicity class

	std::vector<Individual> members_ = {};
	std::vector<size_t> infected_list_ = {};
	std::vector<size_t> healthy_list_ = {};
	std::vector<size_t> missing_list_ = {};

	// contact rate pointers
	matrix<Attribute> contact_rates_;

	/// Change_status of individual i to new status (without controls).
	void set_status_(Individual &member, int new_status);

	///////////////////////////////////////////////////////////////////////////////////////////////
	////////////////////////////  Household : log-likelihood [private]  ///////////////////////////
	///////////////////////////////////////////////////////////////////////////////////////////////

	double obs_log_lik_ = 0; // current observation log-likelihood
	double obs_buffer_ = 0;	 // buffer observation log-likelihood

	dbl_vec smpl_log_lik_ = {}; // diagnostics

	///////////////////////////////////////////////////////////////////////////////////////////////
	///////////////////////////  Household : cache attributes [private]  //////////////////////////
	///////////////////////////////////////////////////////////////////////////////////////////////

	u_dict<int_vec> cache_int_vec_ = {};
	u_dict<std::vector<size_t>> cache_size_t_vec_ = {};
	u_dict<dbl_vec> cache_dbl_vec_ = {};
	u_dict<matrix<double>> cache_dbl_mtx_ = {};

	u_dict<int_vec> buffer_cache_int_vec_ = {};
	u_dict<std::vector<size_t>> buffer_cache_size_t_vec_ = {};
	u_dict<dbl_vec> buffer_cache_dbl_vec_ = {};
	u_dict<matrix<double>> buffer_cache_dbl_mtx_ = {};
};

#endif // HOUSEPI__HOUSEHOLD__H
