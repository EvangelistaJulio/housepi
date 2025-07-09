/*      __
 *     / /_  ____  __  __ ____ __  ____  __
 *    / __ \/ __ \/ / / / ___/ _ \/ __ \/ /
 *   / / / / /_/ / /_/ (__  )  __/ /_/ / /
 *  /_/ /_/\____/\____/____/\___/ ____/_/
 *                             /_/
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

//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//

// status aliases
#define INFECTED_STATUS 1 // infected status (>=)
#define HEALTHY_STATUS 0  // healthy status (==)
#define MISSING_STATUS -1 // missing status (<=)

///////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////  Individual : class  ////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/* Class for individual household study member. */
class Individual
{
public:
	//=======================//
	//===== CONSTRUCTOR =====//
	//=======================//

	Individual() {}
	Individual(size_t index,
			   size_t hhindex,
			   const std::string &id,
			   int status,
			   const data_key_map<double> &data);
	~Individual();

	bool operator==(const Individual &other) const { return this->index_ == other.index_; }
	bool operator!=(const Individual &other) const { return this->index_ != other.index_; }

	//==================//
	//===== GETTER =====//
	//==================//

	const std::string &id() const { return id_; }
	operator size_t() const { return hhindex_; }

	int status() const { return status_; }
	bool is_infected_status() const { return status_ >= INFECTED_STATUS; }
	bool is_healthy_status() const { return status_ == HEALTHY_STATUS; }
	bool is_missing_status() const { return status_ <= MISSING_STATUS; }
	double data(const data_key &d) const { return data_map_.at(d); }
	const dbl_vec &datalist(const data_key &d) const
	{
		if (pdatalist_ == nullptr)
			throw std::invalid_argument(
				"Individual \"" + id_ + "\" has no datalist entries.");
		else
			return pdatalist_->at(d);
	}
	double group_rate(const rate_key &r) const { return rate_map_.at(r).get(); }
	double augm(const augm_key &a) const { return augm_map_.at(a).get(); }
	Variable &get_augm_variable(const augm_key &a) { return augm_map_.at(a); }
	const Variable &get_augm_variable(const augm_key &a) const { return augm_map_.at(a); }

	//==================//
	//===== SETTER =====//
	//==================//

	void set_status(int s) { status_ = s; }

	void set_augm_value(const augm_key &a, double x) { augm_map_.at(a).set_value(x); }
	void assign_pointer_to_datalist(data_key_map<dbl_vec> *y) { pdatalist_ = y; }
	void clear_augm(const augm_key &a) { augm_map_[a] = Variable(); }
	void assign_all_group_rates(const rate_key_mapvec<const double *> &ptrs);
	void assign_augm(const augm_key &a,
					 const std::string &name,
					 double value,
					 bool is_discrete,
					 bool is_constant,
					 const std::string &func_name,
					 const std::string &model_name);
	void init_augm_model(const std::vector<augm_key> &augm_model_list);

	//=========================//
	//===== LOGLIKELIHOOD =====//
	//=========================//

	double get_augm_model_loglike() const { return augm_loglike_; }
	double compute_and_set_augm_model_loglike();
	void recover_augm_model_loglike() { augm_loglike_ = augm_buffer_; };

private:
	// individual group rate nested class
	struct GroupRate
	{
		std::string group_rate_name = "";
		ushort n_rates = 0;
		const double **pointers = {};

		GroupRate(const std::string &name = HSP_NA_STRING,
				  const std::vector<const double *> &pointers = {});
		~GroupRate();
		double get() const;
	};

	std::string id_ = HSP_NA_STRING; // user id
	size_t index_ = 0;				 // housepi internal index
	size_t hhindex_ = 0;			 // household relative index
	int status_ = MISSING_STATUS;	 // status of the individual

	// NOTE: I am forced to store numerical variables
	data_key_map<double> data_map_ = {};	// dataset
	augm_key_map<Variable> augm_map_ = {};	// augmented
	rate_key_map<GroupRate> rate_map_ = {}; // group rates

	data_key_map<dbl_vec> *pdatalist_ = nullptr; // datalist to be associated

	std::vector<std::function<double()>> augm_model_logpdfs_ = {};
	double augm_loglike_ = 0; // current augmented model loglikelihood
	double augm_buffer_ = 0;  // buffer augmented model loglikelihood
};

//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//

///////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////  Household : class  ////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/* Household class */
class Household
{
public:
	//=======================//
	//===== CONSTRUCTOR =====//
	//=======================//

	Household(const std::string &house_id = HSP_NA_STRING,
			  size_t index = SIZE_MAX);
	~Household();

	//=====================//
	//===== ITERATORS =====//
	//=====================//

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
	View<const Individual> infected() const { return View<const Individual>(members_, infected_list_); }
	View<const Individual> healthy() const { return View<const Individual>(members_, healthy_list_); }
	View<const Individual> missing() const { return View<const Individual>(members_, missing_list_); }

	const std::vector<size_t> &infected_hhindexes() const { return infected_list_; }
	const std::vector<size_t> &healthy_hhindexes() const { return healthy_list_; }
	const std::vector<size_t> &missing_hhindexes() const { return missing_list_; }

	//==================//
	//===== GETTER =====//
	//==================//

	// Individual Getter
	const std::string &id(size_t i) const { return members_.at(i).id(); }
	int status(size_t i) const { return members_.at(i).status(); }
	Individual &get_member(size_t i) { return members_.at(i); }
	const Individual &get_member(size_t i) const { return members_.at(i); }

	double contact_rate(size_t from, size_t to) const { return *contact_rates_[from][to]; }

	// Household Getter
	const std::string &house_id() const { return house_id_; }
	size_t house_index() const { return index_; }
	size_t size() const { return house_size_; }

	size_t n_infec_status() const { return infected_list_.size(); }
	size_t n_healthy_status() const { return healthy_list_.size(); }
	size_t n_missing_status() const { return missing_list_.size(); }

	ulong multiplicity() const { return house_multiplicity_; }
	size_t multiplicity_class_index() const { return multiplicity_class_; }
	bool is_multiplicity_head_class() const { return (index_ == multiplicity_class_); }

	//==================//
	//===== SETTER =====//
	//==================//

	void set_multiplicity(ulong m, size_t class_idx);
	void assign_member(size_t index,
					   const std::string &id,
					   int status,
					   const data_key_map<double> &data);
	bool change_status(Individual &member, int new_status);
	bool change_status(size_t i, int new_status) { return change_status(members_.at(i), new_status); }
	void assign_contact_rate(size_t from, size_t to, const double *x);

	//============================//
	//===== CACHE ATTRIBUTES =====//
	//============================//

	/* Check if any cache quantity is associated to the household. */
	bool is_any_cache() { return !(cache_int_vec_.empty() && cache_dbl_vec_.empty() && cache_dbl_mtx_.empty()); }
	void buffer_all_cache();

	//----- Cached integer vector -----//
	void assign_cache_integer_vector(const cache_key &e) { cache_int_vec_[e] = int_vec(size(), 0); }
	int_vec &cache_integer_vector(const cache_key &e) { return cache_int_vec_.at(e); }
	const int_vec &cache_integer_vector(const cache_key &e) const { return cache_int_vec_.at(e); }

	//----- Cached double vector -----//
	void assign_cache_double_vector(const cache_key &e) { cache_dbl_vec_[e] = dbl_vec(size(), 0.); }
	dbl_vec &cache_double_vector(const cache_key &e) { return cache_dbl_vec_.at(e); }
	const dbl_vec &cache_double_vector(const cache_key &e) const { return cache_dbl_vec_.at(e); }

	//----- Cached double matrix -----//
	void assign_cache_double_matrix(const cache_key &e) { cache_dbl_mtx_[e] = matrix<double>(size(), dbl_vec(size(), 0.)); }
	matrix<double> &cache_double_matrix(const cache_key &e) { return cache_dbl_mtx_.at(e); }
	const matrix<double> &cache_double_matrix(const cache_key &e) const { return cache_dbl_mtx_.at(e); }

	//=========================//
	//===== LOGLIKELIHOOD =====//
	//=========================//

	double get_obs_loglike() const { return obs_loglike_; }
	void set_obs_loglike(double x);

	void recover_buffer_values(bool is_cache_depen);

private:
	std::string house_id_ = HSP_NA_STRING; // user house id
	size_t index_ = -1;					   // housepi internal index
	ushort house_size_ = 0;				   // n. of house members
	ulong house_multiplicity_ = 1;		   // n. of equivalent houses
	size_t multiplicity_class_ = -1;	   // index of multiplicity class

	std::vector<Individual> members_ = {}; // members

	// list of members for each status
	std::vector<size_t> infected_list_ = {};
	std::vector<size_t> healthy_list_ = {};
	std::vector<size_t> missing_list_ = {};

	// contact rate pointers
	matrix<const double *> contact_rates_;

	//===== LOGLIKELIHOOD =====//

	double obs_loglike_ = 0; // current observation loglikelihood
	double obs_buffer_ = 0;	 // buffer observation loglikelihood

	dbl_vec smpl_loglike_ = {}; // diagnostics

	//===== EXTRA ATTRIBUTES =====//

	//----- Extra Attributes -----//
	cache_key_map<int_vec> cache_int_vec_ = {};
	cache_key_map<dbl_vec> cache_dbl_vec_ = {};
	cache_key_map<matrix<double>> cache_dbl_mtx_ = {};

	//----- Buffer EXTRA Attributes -----//
	cache_key_map<int_vec> buffer_cache_int_vec_ = {};
	cache_key_map<dbl_vec> buffer_cache_dbl_vec_ = {};
	cache_key_map<matrix<double>> buffer_cache_dbl_mtx_ = {};
};

#endif // HOUSEPI__HOUSEHOLD__H
