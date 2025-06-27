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

	bool operator==(const Individual &other) const { return this->_index == other._index; }
	bool operator!=(const Individual &other) const { return this->_index != other._index; }

	//==================//
	//===== GETTER =====//
	//==================//

	const std::string &id() const { return _id; }
	operator size_t() const { return _hhindex; }

	int status() const { return _status; }
	double data(const data_key &d) const { return _data_map.at(d); }
	const dbl_vec &datalist(const data_key &d) const
	{
		if (_datalist_ptr == nullptr)
			throw std::invalid_argument(
				"Individual \"" + _id + "\" has no datalist entries.");
		else
			return _datalist_ptr->at(d);
	}
	double group_rate(const rate_key &r) const { return _rate_map.at(r).get(); }
	double augm(const augm_key &a) const { return _augm_map.at(a).get(); }
	Variable &get_augm_variable(const augm_key &a) { return _augm_map.at(a); }
	const Variable &get_augm_variable(const augm_key &a) const { return _augm_map.at(a); }

	const std::vector<data_key> keys_dataset() const;
	const std::vector<data_key> keys_datalist() const;
	const std::vector<rate_key> keys_group_rate() const;
	const std::vector<augm_key> keys_augm() const;

	//==================//
	//===== SETTER =====//
	//==================//

	void set_status(int s) { _status = s; }

	void set_augm_value(const augm_key &a, double x) { _augm_map.at(a).set_value(x); }
	void assign_pointer_to_datalist(data_key_map<dbl_vec> *y) { _datalist_ptr = y; }
	void clear_augm(const augm_key &a) { _augm_map[a] = Variable(); }
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

	double get_augm_model_loglike() const { return _augm_loglike; }
	double compute_and_set_augm_model_loglike();
	void recover_augm_model_loglike() { _augm_loglike = _augm_buffer; };

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

	std::string _id = HSP_NA_STRING; // user id
	size_t _index = 0;				 // housepi internal index
	size_t _hhindex = 0;			 // household relative index
	int _status = MISSING_STATUS;	 // status of the individual

	// NOTE: I am forced to store numerical variables
	data_key_map<double> _data_map = {};	// dataset
	augm_key_map<Variable> _augm_map = {};	// augmented
	rate_key_map<GroupRate> _rate_map = {}; // group rates

	data_key_map<dbl_vec> *_datalist_ptr = nullptr; // datalist to be associated

	std::vector<std::function<double()>> _augm_model_logpdfs = {};
	double _augm_loglike = 0; // current augmented model loglikelihood
	double _augm_buffer = 0;  // buffer augmented model loglikelihood
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

	//====================//
	//===== ITERATOR =====//
	//====================//

	auto begin() { return _members.begin(); }
	auto end() { return _members.end(); }
	auto begin() const { return _members.begin(); }
	auto end() const { return _members.end(); }

	// T should be Individual for mutable access
	// or const Individual for const access.
	template <typename T>
	class ViewIterator
	{
		// conditional to select the vector type if constant or not
		using VectorType = std::conditional_t<std::is_const<T>::value, const std::vector<Individual>, std::vector<Individual>>;

		VectorType &members_ptr;			// Either vector<Individual>& or const vector<Individual>&
		const std::vector<size_t> &indices; // indices always read-only
		size_t pos;

	public:
		ViewIterator(VectorType &m,
					 const std::vector<size_t> &i,
					 size_t p)
			: members_ptr(m), indices(i), pos(p) {}

		T &operator*() const { return members_ptr[indices[pos]]; }
		T *operator->() const { return &members_ptr[indices[pos]]; }

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
		VectorType &members;
		const std::vector<size_t> &indices;

		View(VectorType &m, const std::vector<size_t> &i)
			: members(m), indices(i) {}

		iterator begin() { return iterator(members, indices, 0); }
		iterator end() { return iterator(members, indices, indices.size()); }
	};

	View<Individual> infected() { return View<Individual>(_members, _infected_list); }
	View<Individual> healthy() { return View<Individual>(_members, _healthy_list); }
	View<Individual> missing() { return View<Individual>(_members, _missing_list); }
	View<const Individual> infected() const { return View<const Individual>(_members, _infected_list); }
	View<const Individual> healthy() const { return View<const Individual>(_members, _healthy_list); }
	View<const Individual> missing() const { return View<const Individual>(_members, _missing_list); }

	const std::vector<size_t> &infected_hhindexes() const { return _infected_list; }
	const std::vector<size_t> &healthy_hhindexes() const { return _healthy_list; }
	const std::vector<size_t> &missing_hhindexes() const { return _missing_list; }

	//==================//
	//===== GETTER =====//
	//==================//

	// Individual Getter
	const std::string &id(size_t i) const { return _members.at(i).id(); }
	int status(size_t i) const { return _members.at(i).status(); }
	Individual &get_member(size_t i) { return _members.at(i); }
	const Individual &get_member(size_t i) const { return _members.at(i); }

	double contact_rate(size_t from, size_t to) const { return *_contact_rates[from][to]; }

	// Household Getter
	const std::string &house_id() const { return _house_id; }
	size_t house_index() const { return _index; }
	size_t size() const { return _house_size; }

	size_t n_infec_status() const { return _infected_list.size(); }
	size_t n_healthy_status() const { return _healthy_list.size(); }
	size_t n_missing_status() const { return _missing_list.size(); }

	ulong multiplicity() const { return _house_multiplicity; }
	size_t multiplicity_class_index() const { return _multiplicity_class; }
	bool is_multiplicity_head_class() const { return (_index == _multiplicity_class); }

	//==================//
	//===== SETTER =====//
	//==================//

	void set_multiplicity(ulong m, size_t class_idx);
	void assign_member(size_t index,
					   const std::string &id,
					   int status,
					   const data_key_map<double> &data);
	bool change_status(Individual &member, int new_status);
	bool change_status(size_t i, int new_status) { return change_status(_members.at(i), new_status); }
	void assign_contact_rate(size_t from, size_t to, const double *x);

	//============================//
	//===== CACHE ATTRIBUTES =====//
	//============================//

	/* Check if any cache quantity is associated to the household. */
	bool is_any_cache() { return !(_cache_int_vec.empty() && _cache_dbl_vec.empty() && _cache_dbl_mtx.empty()); }
	void buffer_all_cache();

	//----- Cached integer vector -----//
	void assign_cache_integer_vector(const cache_key &e) { _cache_int_vec[e] = int_vec(size(), 0); }
	int_vec &cache_integer_vector(const cache_key &e) { return _cache_int_vec.at(e); }
	const int_vec &cache_integer_vector(const cache_key &e) const { return _cache_int_vec.at(e); }

	//----- Cached double vector -----//
	void assign_cache_double_vector(const cache_key &e) { _cache_dbl_vec[e] = dbl_vec(size(), 0.); }
	dbl_vec &cache_double_vector(const cache_key &e) { return _cache_dbl_vec.at(e); }
	const dbl_vec &cache_double_vector(const cache_key &e) const { return _cache_dbl_vec.at(e); }

	//----- Cached double matrix -----//
	void assign_cache_double_matrix(const cache_key &e) { _cache_dbl_mtx[e] = matrix<double>(size(), dbl_vec(size(), 0.)); }
	matrix<double> &cache_double_matrix(const cache_key &e) { return _cache_dbl_mtx.at(e); }
	const matrix<double> &cache_double_matrix(const cache_key &e) const { return _cache_dbl_mtx.at(e); }

	//=========================//
	//===== LOGLIKELIHOOD =====//
	//=========================//

	double get_obs_loglike() const { return _obs_loglike; }
	void set_obs_loglike(double x);

	void recover_buffer_values(bool is_cache_depen);

private:
	std::string _house_id = HSP_NA_STRING; // user house id
	size_t _index = -1;					   // housepi internal index
	ushort _house_size = 0;				   // n. of house members
	ulong _house_multiplicity = 1;		   // n. of equivalent houses
	size_t _multiplicity_class = -1;	   // index of multiplicity class

	std::vector<Individual> _members = {}; // members

	// list of members for each status
	std::vector<size_t> _infected_list = {};
	std::vector<size_t> _healthy_list = {};
	std::vector<size_t> _missing_list = {};

	// contact rate pointers
	matrix<const double *> _contact_rates;

	//===== LOGLIKELIHOOD =====//

	double _obs_loglike = 0; // current observation loglikelihood
	double _obs_buffer = 0;	 // buffer observation loglikelihood

	dbl_vec _smpl_loglike = {}; // diagnostics

	//===== EXTRA ATTRIBUTES =====//

	//----- Extra Attributes -----//
	cache_key_map<int_vec> _cache_int_vec = {};
	cache_key_map<dbl_vec> _cache_dbl_vec = {};
	cache_key_map<matrix<double>> _cache_dbl_mtx = {};

	//----- Buffer EXTRA Attributes -----//
	cache_key_map<int_vec> _buffer_cache_int_vec = {};
	cache_key_map<dbl_vec> _buffer_cache_dbl_vec = {};
	cache_key_map<matrix<double>> _buffer_cache_dbl_mtx = {};
};

#endif // HOUSEPI__HOUSEHOLD__H
