/*
 *        __
 *       / /_  ____  __  __ ____ __  ____  __
 *      / __ \/ __ \/ / / / ___/ _ \/ __ \/ /
 *     / / / / /_/ / /_/ (__  )  __/ /_/ / /
 *    /_/ /_/\____/\____/____/\___/ ____/_/
 *                               /_/
 *
 * housepi - Epidemiological Analysis of Household Outbreaks.
 *
 * Author: Francesco Camaglia
 * Copyright © 2025 Institut Pasteur (Paris), CNRS and INSERM.
 *
 * Distributed under the terms of the GNU General Public License (GPLv3).
 * See the LICENSE file for details.
 *
 * household.h
 *
 */

#ifndef HOUSEPI__HOUSEHOLD__H
#define HOUSEPI__HOUSEHOLD__H

#include <functional>

#include "tools/types.h"
#include "tools/na_utils.h"
#include "tools/operators.h"
#include "tools/status.h"
#include "tools/rates.h"
#include "function.h"
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
	u_dict<std::vector<std::function<double()>>> model_logratios_ = {};
};

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////// *** Individual : class *** /////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/// Forward declaration.
class Household;

/// Class for individual household study member.
class Individual : public AugmentedMap
{
public:
	///////////////////////////////////////////////////////////////////////////////////////////////
	//////////////////////////////////  Individual : constructor  /////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////////////////

	Individual() {}
	Individual(size_t index,
			   size_t hhindex,
			   const std::string &id,
			   int status,
			   const u_dict<double> &data)
	{
		index_ = index;
		hhindex_ = hhindex;
		id_ = id;
		status_ = status;
		data_map_ = data;
	}

	/// Sugar to tell if the corresponding members are the same.
	bool operator==(const Individual &other) const { return this->index_ == other.index_; }
	/// Sugar to tell if the corresponding members are different
	bool operator!=(const Individual &other) const { return this->index_ != other.index_; }
	/// Sugar to call an Individual object as house index.
	operator size_t() const { return hhindex_; }

	///////////////////////////////////////////////////////////////////////////////////////////////
	////////////////////////////////////  Individual : getter  ////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////////////////

	const std::string &id() const { return id_; }
	std::string err_name() const override final { return "individual '" + id_ + "'"; }
	int status() const { return status_; }
	bool is_infected_status() const { return hsp::status::is_infected(status_); }
	bool is_healthy_status() const { return hsp::status::is_healthy(status_); }
	bool is_missing_status() const { return hsp::status::is_missing(status_); }
	double data(const std::string &d) const { return safe_at(data_map_, d, "data"); }
	const dbl_vec &data_long(const std::string &d) const
	{
		// return a reference to a stable empty vector when no datalong is attached
		static const dbl_vec EMPTY_DBL_VEC = dbl_vec();
		return (pdatalong_ == nullptr) ? EMPTY_DBL_VEC : pdatalong_->at(d);
	}
	double group_rate(const std::string &g) const { return safe_at(group_rate_map_, g, "group rate").get(); }
	const Attribute &get_group_rate_attribute(const std::string &g) const { return safe_at(group_rate_map_, g, "group rate"); }

	///////////////////////////////////////////////////////////////////////////////////////////////
	////////////////////////////////////  Individual : setter  ////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////////////////

	/// Initialize pointer to data long format.
	void init_pointer_to_data_long(u_dict<dbl_vec> *y) { pdatalong_ = y; }
	/// Initialize group rate according to given attribute.
	void init_group_rate(const std::string &g, const Attribute &attr)
	{
		group_rate_map_[g].transfer(attr);
		group_rate_map_[g].set_name(g, "group rate");
	}
	/// Initialize group rate to reference.
	void init_group_rate_to_reference(const std::string &g)
	{
		group_rate_map_[g].transfer(Attribute(hsp::rates::VALUE_REFERENCE));
		group_rate_map_[g].set_name(g, "group rate");
	}

private:
	std::string id_ = hsp::na::STRING;	// user id
	size_t index_ = 0;					// housepi internal index
	size_t hhindex_ = 0;				// household relative index
	int status_ = hsp::status::MISSING; // status of the individual

	// NOTE: allowd data variables are numerical only
	u_dict<double> data_map_ = {};			// data
	u_dict<Attribute> group_rate_map_ = {}; // group rates

	u_dict<dbl_vec> *pdatalong_ = nullptr; // data long format associated

	friend class Household;
};

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////  Household : class  ////////////////////////////////////////
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
	void set_healthy_status(Individual &member) { set_status_(member, hsp::status::HEALTHY); }
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
