/*
 *        __
 *       / /_  ____  __  __ ____ __  ____  __
 *      / __ \/ __ \/ / / / ___/ _ \/ __ \/ /
 *     / / / / /_/ / /_/ (__  )  __/ /_/ / /
 *    /_/ /_/\____/\____/____/\___/ ____/_/
 *                               /_/
 *
 * Copyright (C) 2025 Francesco Camaglia, Institut Pasteur
 *
 * household.h
 *
 */

#ifndef HOUSEPI__HOUSEHOLD__H
#define HOUSEPI__HOUSEHOLD__H

#include <memory>
#include <string>
#include <vector>

#include "tools.h"
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

	void init_augm(const std::string &a,
				   double value,
				   bool is_const,
				   bool is_discr,
				   const DefinitionParser &func_name,
				   bool is_prim,
				   const DefinitionParser &model_name);

	void set_augm_value(const std::string &a, double x) { safe_at(augm_map_, a, "augmented").set(x); }
	void clear_augm(const std::string &a) { augm_map_[a] = Variable(); }
	void user_reset(rnd_gen &prng);

	///////////////////////////////////////////////////////////////////////////////////////////////
	///////////////////////////////////  AugmentedMap : getter  ///////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////////////////

	virtual std::string err_name() const { return ""; }
	double augm(const std::string &a) const { return safe_at(augm_map_, a, "augmented").get(); }
	Variable &get_augm_variable(const std::string &a) { return safe_at(augm_map_, a, "augmented"); }
	const Variable &get_augm_variable(const std::string &a) const { return safe_at(augm_map_, a, "augmented"); }

	///////////////////////////////////////////////////////////////////////////////////////////////
	/////////////////////////////////  Individual : loglikelihood  ////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////////////////

	void init_augm_model_logratios(const u_dict<str_vec> &models_ordered);
	double get_augm_model_loglike() const { return loglike_; }
	double compute_and_set_augm_model_loglike();
	double compute_and_set_augm_model_logratio(const std::string &a);
	void recover_augm_model_loglike() { loglike_ = buffer_; };

private:
	u_dict<Variable> augm_map_ = {}; /// dictionary of augmented definitions
	
	double loglike_ = 0; // current augmented model loglikelihood
	double buffer_ = 0;  // buffer augmented model loglikelihood
	
	// container for the logratio lists according to dependencies 
	u_dict<std::vector<std::function<double()>>> model_logratios_ = {};
};

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////// *** Individual : class *** /////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

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
			   const u_dict<double> &data);

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
	std::string err_name() const override final { return "individual \"" + id_ + "\""; }
	int status() const { return status_; }
	bool is_infected_status() const { return HSP_IS_INFECTED(status_); }
	bool is_healthy_status() const { return HSP_IS_HEALTHY(status_); }
	bool is_missing_status() const { return HSP_IS_MISSING(status_); }
	double data(const std::string &d) const { return safe_at(data_map_, d, "data"); }
	const dbl_vec &datalong(const std::string &d) const
	{
		if (pdatalong_ == nullptr)
			throw std::invalid_argument(
				"No datalong entries for " + err_name() + ".");
		return pdatalong_->at(d);
	}
	double group_rate(const std::string &r) const { return safe_at(rate_map_, r, "group rate").get(); }

	///////////////////////////////////////////////////////////////////////////////////////////////
	////////////////////////////////////  Individual : setter  ////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////////////////

	void set_status(int s) { status_ = s; }
	void init_pointer_to_datalong(u_dict<dbl_vec> *y) { pdatalong_ = y; }
	void init_all_group_rates(const u_dict<std::vector<observ_dbl_ptr>> &ptrs);

private:
	///////////////////////////////////////////////////////////////////////////////////////////////
	//////////////////////////////  Individual : GroupRate_ : struct  /////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////////////////

	struct GroupRate_
	{
	public:
		GroupRate_(const std::string &name = HSP_NA_STRING,
				   const std::vector<observ_dbl_ptr> &ptrs = {});

		double get() const
		{
			double output = 1.0;
			for (observ_dbl_ptr ptr : pointers_)
				output *= *ptr;
			return output;
		}

	private:
		std::string name_ = "";
		std::vector<observ_dbl_ptr> pointers_ = {};
	};

	std::string id_ = HSP_NA_STRING;  // user id
	size_t index_ = 0;				  // housepi internal index
	size_t hhindex_ = 0;			  // household relative index
	int status_ = HSP_MISSING_STATUS; // status of the individual

	// NOTE: allowd data variables are numerical only
	u_dict<double> data_map_ = {};	   // data
	u_dict<GroupRate_> rate_map_ = {}; // group rates

	u_dict<dbl_vec> *pdatalong_ = nullptr; // datalong to be associated
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

	Household(const std::string &house_id = HSP_NA_STRING,
			  size_t index = SIZE_MAX);
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
	template <typename T>
	class ViewIterator
	{
		// conditional to select the vector type if constant or not
		using VectorType = std::conditional_t<std::is_const<T>::value, const std::vector<Individual>, std::vector<Individual>>;

		VectorType &pVector;
		const std::vector<size_t> &indexes;
		size_t pos;

	public:
		ViewIterator(VectorType &m, const std::vector<size_t> &i, size_t p)
			: pVector(m), indexes(i), pos(p) {}

		T &operator*() const { return pVector[indexes[pos]]; }
		T *operator->() const { return &pVector[indexes[pos]]; }

		ViewIterator &operator++()
		{
			++pos;
			return *this;
		}
		ViewIterator operator++(int)
		{
			ViewIterator tmp = *this;
			++(*this);
			return tmp;
		}

		bool operator==(const ViewIterator &other) const { return pos == other.pos; }
		bool operator!=(const ViewIterator &other) const { return !(*this == other); }
	};

	template <typename T>
	class View
	{
		// conditional to select the vector type if constant or not
		using VectorType = std::conditional_t<std::is_const<T>::value, const std::vector<Individual>, std::vector<Individual>>;
		using iterator = ViewIterator<T>;

	public:
		VectorType &vec;
		const std::vector<size_t> &indexes;

		View(VectorType &m, const std::vector<size_t> &i)
			: vec(m), indexes(i) {}

		iterator begin() { return iterator(vec, indexes, 0); }
		iterator end() { return iterator(vec, indexes, indexes.size()); }
	};

	View<Individual> infected() { return View<Individual>(members_, infected_list_); }
	View<Individual> healthy() { return View<Individual>(members_, healthy_list_); }
	View<Individual> missing() { return View<Individual>(members_, missing_list_); }
	View<Individual> select(const std::vector<size_t> lst) { return View<Individual>(members_, lst); }
	View<const Individual> infected() const { return View<const Individual>(members_, infected_list_); }
	View<const Individual> healthy() const { return View<const Individual>(members_, healthy_list_); }
	View<const Individual> missing() const { return View<const Individual>(members_, missing_list_); }
	View<const Individual> select(const std::vector<size_t> lst) const { return View<const Individual>(members_, lst); }

	const std::vector<size_t> &infected_hhindexes() const { return infected_list_; }
	const std::vector<size_t> &healthy_hhindexes() const { return healthy_list_; }
	const std::vector<size_t> &missing_hhindexes() const { return missing_list_; }

	///////////////////////////////////////////////////////////////////////////////////////////////
	////////////////////////////////////  Household : getter  /////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////////////////

	// Individual Getter
	Individual &get_member(size_t i) { return members_.at(i); }
	const Individual &get_member(size_t i) const { return members_.at(i); }

	double contact_rate(size_t from, size_t to) const { return *contact_rates_[from][to]; }

	// Household Getter
	const std::string &house_id() const { return house_id_; }
	std::string err_name() const { return "household \"" + house_id_ + "\""; }
	size_t house_index() const { return index_; }
	size_t size() const { return house_size_; }

	size_t n_infec_status() const { return infected_list_.size(); }
	size_t n_healthy_status() const { return healthy_list_.size(); }
	size_t n_missing_status() const { return missing_list_.size(); }

	unsigned long multiplicity() const { return house_multiplicity_; }
	size_t multiplicity_class_index() const { return multiplicity_class_; }
	bool is_multiplicity_head_class() const { return (index_ == multiplicity_class_); }

	///////////////////////////////////////////////////////////////////////////////////////////////
	////////////////////////////////////  Household : setter  /////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////////////////

	void set_multiplicity(unsigned long m, size_t class_idx);
	void init_member(size_t index, const std::string &id, int status, const u_dict<double> &data);
	bool change_status(Individual &member, int new_status);
	void init_contact_rates(size_t from, size_t to, observ_dbl_ptr x);

	///////////////////////////////////////////////////////////////////////////////////////////////
	///////////////////////////////  Household : cache attributes  ////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////////////////

	/// Check if any cache quantity is associated to the household.
	bool is_any_cache() { return !(cache_int_vec_.empty() && cache_dbl_vec_.empty() && cache_dbl_mtx_.empty()); }
	void buffer_all_cache();

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
	/////////////////////////////////  Household : loglikelihood  /////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////////////////

	double get_obs_loglike() const { return obs_loglike_; }
	void set_obs_loglike(double x);
	void recover_buffer_values(bool is_cache_depen);

private:
	std::string house_id_ = HSP_NA_STRING; // user house id
	size_t index_ = -1;					   // housepi internal index
	unsigned short house_size_ = 0;		   // n. of house members
	unsigned long house_multiplicity_ = 1; // n. of equivalent houses
	size_t multiplicity_class_ = -1;	   // index of multiplicity class

	std::vector<Individual> members_ = {}; // members

	// list of members for each status
	std::vector<size_t> infected_list_ = {};
	std::vector<size_t> healthy_list_ = {};
	std::vector<size_t> missing_list_ = {};

	// contact rate pointers
	matrix<observ_dbl_ptr> contact_rates_;

	///////////////////////////////////////////////////////////////////////////////////////////////
	////////////////////////////  Household : loglikelihood [private]  ////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////////////////

	double obs_loglike_ = 0; // current observation loglikelihood
	double obs_buffer_ = 0;	 // buffer observation loglikelihood

	dbl_vec smpl_loglike_ = {}; // diagnostics

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
