/*
 * household.cpp
 *
 * Copyright (C) 2025 Francesco Camaglia, Institut Pasteur
 */

#include <algorithm>

#include "variable.h"
#include "household.h"

//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//

///////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////  Individual : class  ////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/* Struct for individual household study member. */
class Individual
{
public:
	//===== CONSTRUCTOR =====//
	Individual(const std::string &id = NULL_STRING,
			   int status = MISSING_STATUS,
			   const dataIDmap<double> &data = {});
	Individual(Individual &&other) noexcept; // Move constructor
	~Individual();

	//===== PUBLIC ATTRIBUTES =====//
	std::string id = NULL_STRING; // user id
	int status = MISSING_STATUS;  // status of the individual

	//===== DATA =====//
	// UPDATE: if I had the chance to use boost I would do:
	// NOTE: BH R library doesn't contain variant
	// > #include <boost/variant.hpp>
	// > ...
	// > using ValueType = boost::variant<double, std::string>;
	// This way I could have string data (nice group_rates names)
	// The the getter would become something like:
	// template<typename T>
	// T get(const dataID& d) const { return boost::get<T>(data.find(d)->second); }
	//
	dataIDmap<double> data = {}; // database associated

	//===== AUGMENTED VARIABLES =====//
	augmIDmap<Variable> augm_map = {}; // augmented variables

	//===== GROUP RATES =====//
	void assign_all_group_rates(const rateIDmapvec<const double *> &ptrs);

	struct GroupRate
	{
		std::string name = "";
		ushort n_rates = 0;
		const double **pointers = {};

		GroupRate(const std::string &name = NULL_STRING,
				  const std::vector<const double *> &pointers = {});
		~GroupRate();
		double get() const;
	};
	rateIDmap<GroupRate> rate_map; // group rates associated
};

//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////  Individual : constructors  ///////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/* Individual constructor */
Individual::Individual(const std::string &id,
					   int status,
					   const dataIDmap<double> &data)
{
	// associate the new id
	this->id = id;
	// associate the new status
	this->status = (status < MISSING_STATUS) ? MISSING_STATUS : status;
	// associate the new data
	this->data = data;
}

/* Individual move constructor. */
// UPDATE : is this really necessary?
Individual::Individual(Individual &&other) noexcept
	: id(std::move(other.id)),
	  status(std::move(other.status)),
	  data(std::move(other.data)),
	  augm_map(std::move(other.augm_map)),
	  rate_map(std::move(other.rate_map)) {}

/* Individual deconstructor */
Individual::~Individual() = default;

///////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////  Individual : assign all group rates  //////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

void Individual::assign_all_group_rates(const rateIDmapvec<const double *> &ptrs)
{
	for (const auto &pair : ptrs)
	{
		rateID r = pair.first;
		std::vector<const double *> vec_par_ptr = pair.second;

		if (vec_par_ptr.empty())
			std::invalid_argument(
				"Missing rate \"" + r + "\" for individual " + id + ".");
		else
			this->rate_map[r] = GroupRate(r, vec_par_ptr);
	}
}

//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//

///////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////  Individual::Group Rate : class  /////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

//===== CONSTRUCTOR	=====//
Individual::GroupRate::GroupRate(const std::string &name,
								 const std::vector<const double *> &pointers)
{
	if (!tools::is_na(name))
	{
		this->name = name;

		if (pointers.empty())
			throw std::invalid_argument(
				"No parameters associated to group rate \"" + name + "\".");

		this->n_rates = pointers.size();
		this->pointers = new const double *[pointers.size()];
		for (ushort c = 0; c < this->n_rates; ++c)
			this->pointers[c] = pointers.at(c);
	}
}
Individual::GroupRate::~GroupRate() = default;

//===== GETTER =====//
double Individual::GroupRate::get() const
{
	double output = 1.0;
	for (ushort r = 0; r < n_rates; ++r)
		output *= *pointers[r];

	return output;
}

//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//

class Household::Impl
{
public:
	std::string house_id = NULL_STRING; // user id
	size_t index = -1;					// index of household
	ushort house_size = 0;				// number of house members

	std::vector<Individual> members = {}; //

	double house_loglike = NULL_DOUBLE; // current loglikelihood value
	double buffer = NULL_DOUBLE;		// buffer loglikelihood value

	dbl_vec dgn_loglike = {}; // diagnostics

	void assign_individual(const std::string &id,
						   int status,
						   const dataIDmap<double> &data);
	bool change_status(indiv i, int new_status);

	// contact rate pointers
	matrix<const double *> contact_rates;

	ushort n_infectious = 0;
	ushort n_susceptible = 0;
	ushort n_missing_status = 0;

	// wheter members correspond to a given status or not
	std::vector<indiv> infectious_list = {};
	std::vector<indiv> susceptible_list = {};
	std::vector<indiv> missing_status_list = {};
};


//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//

/* Link an individual to the household. */
void Household::Impl::assign_individual(const std::string &id,
										int status,
										const dataIDmap<double> &data)
{
	//----- update status lists within the household -----//
	if (status >= INFECTED_STATUS)
	{
		++n_infectious;
		infectious_list.push_back(house_size);
	}
	else if (status == HEALTHY_STATUS)
	{
		++n_susceptible;
		susceptible_list.push_back(house_size);
	}
	else
	{
		++n_missing_status;
		missing_status_list.push_back(house_size);
	}

	//----- emplace new individual -----//
	members.emplace_back(id, status, data);
	++house_size;
}

/* Change_status of individual i to new status (without controls).
	NOTE : this introduces augmented variables -> it requires RJMCMC in case of inference. */
bool Household::Impl::change_status(indiv i, int new_status)
{
	Individual &this_indiv = members[i];
	int old_status = this_indiv.status;

	if (new_status == old_status)
		return false;

	else
	{
		if (old_status >= INFECTED_STATUS)
		{
			--n_infectious;
			infectious_list.erase(std::find(infectious_list.begin(), infectious_list.end(), i));
		}
		else if (old_status == HEALTHY_STATUS)
		{
			--n_susceptible;
			susceptible_list.erase(std::find(susceptible_list.begin(), susceptible_list.end(), i));
		}
		else if (old_status == MISSING_STATUS)
		{
			--n_missing_status;
			missing_status_list.erase(std::find(missing_status_list.begin(), missing_status_list.end(), i));
		}
		else
			throw std::runtime_error(
				"Unknown status \"" + std::to_string(old_status) + "\" for individual " + this_indiv.id + ".");

		if (new_status >= INFECTED_STATUS)
		{
			++n_infectious;
			infectious_list.push_back(i);
		}
		else if (new_status == HEALTHY_STATUS)
		{
			++n_susceptible;
			susceptible_list.push_back(i);
		}
		else if (new_status == MISSING_STATUS)
		{
			++n_missing_status;
			missing_status_list.push_back(i);
		}
		else
			throw std::runtime_error(
				"Unknown new status \"" + std::to_string(new_status) + "\" for individual " + this_indiv.id + ".");

		this_indiv.status = new_status;
		return true;
	}
}

//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//

//===== CONSTRUCTOR =====//
Household::Household(const std::string &house_id,
					 size_t H) : pImpl(std::make_unique<Impl>())
{
	pImpl->house_id = house_id;
	pImpl->index = H;
}
Household::Household(Household &&other) noexcept = default;
Household &Household::operator=(Household &&other) noexcept = default;
Household::~Household() = default;

void Household::assign_individual(const std::string &id,
								  int status,
								  const dataIDmap<double> &data)
{
	pImpl->assign_individual(id, status, data);
}

/* Change_status of individual i to new status (without controls).
	NOTE : this introduces augmented variables -> it requires RJMCMC in case of inference. */
bool Household::change_status(indiv i, int new_status) { return pImpl->change_status(i, new_status); }

/* Associate a pointer to parameter value for the contact from->to.
Allocate memory if it  */
void Household::assign_contact_rate(indiv from, indiv to, const double *x)
{
	// initialize contact rate matrix if necessary
	if (pImpl->contact_rates.empty())
		pImpl->contact_rates.resize(size(), std::vector<const double *>(size(), &REFERENCE_VALUE));

	pImpl->contact_rates[from][to] = x;
}

void Household::set_loglike(double x)
{
	pImpl->buffer = pImpl->house_loglike;
	pImpl->house_loglike = x;
}
void Household::restore_loglike() { pImpl->house_loglike = pImpl->buffer; }

//===== GETTER =====//
const std::string &Household::house_id() const { return pImpl->house_id; }
const std::string &Household::id(indiv i) const { return pImpl->members[i].id; }
size_t Household::index() const { return pImpl->index; }
ushort Household::size() const { return pImpl->house_size; }

int Household::status(indiv i) const { return pImpl->members[i].status; }
double Household::contact_rate(indiv from, indiv to) const { return *pImpl->contact_rates[from][to]; }

double Household::data(indiv i, const dataID &d) const { return pImpl->members[i].data.at(d); }
double Household::group_rate(indiv i, const rateID &r) const { return pImpl->members[i].rate_map.at(r).get(); }
double Household::augm(indiv i, const augmID &a) const { return pImpl->members[i].augm_map.at(a).get_value(); }

double Household::augm_model_logpdf(indiv i, const augmID &a) const { return pImpl->members[i].augm_map.at(a).distr_logpdf(); }
double Household::augm_model_logratio(indiv i, const augmID &a) const { return pImpl->members[i].augm_map.at(a).distr_logratio(); }
double Household::augm_propose_logratio(indiv i, const augmID &a) const { return pImpl->members[i].augm_map.at(a).propose_logratio(); }
double Household::augm_move_logratio(indiv i, const augmID &a) const { return pImpl->members[i].augm_map.at(a).get_move_logratio(); }

ushort Household::n_infectious() const { return pImpl->n_infectious; }
ushort Household::n_susceptible() const { return pImpl->n_susceptible; }
ushort Household::n_missing_status() const { return pImpl->n_missing_status; }

double Household::get_loglike() const { return pImpl->house_loglike; }
double Household::get_delta_loglike() const { return pImpl->house_loglike - pImpl->buffer; }

const std::vector<indiv> & Household::infectious_list() const { return pImpl->infectious_list; }
const std::vector<indiv> & Household::susceptible_list() const { return pImpl->susceptible_list; }
const std::vector<indiv> & Household::missing_status_list() const { return pImpl->missing_status_list; }

///////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////  Household : setter  ///////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

void Household::set_data_value(indiv i, const dataID &d, double x)
{
	pImpl->members[i].data[d] = x;
}

void Household::assign_all_group_rates(indiv i, const rateIDmapvec<const double *> &x)
{
	pImpl->members[i].assign_all_group_rates(x);
}

/* Assign definitions for the augmented variables. */
void Household::assign_augm(indiv i, const augmID &a, const std::string &name, double value,
							bool is_discrete, const std::string &model_name,
							const std::string &rule_name, double rule_scale, double rule_from)
{
	Variable &this_augm = pImpl->members[i].augm_map.at(a);
	this_augm.init(name, value, is_discrete);
	this_augm.assign_distribution(model_name);
	this_augm.assign_proposer(rule_name, rule_scale, rule_from);
}

/* Assign scalar argument for augmented variable model. */
void Household::assign_augm_model_scalar_arg(indiv i, const augmID &a, const std::string &arg_name, double x)
{
	pImpl->members[i].augm_map.at(a).assign_distr_scalar_arg(arg_name, x);
}
/* Assign pointer to argument for augmented variable model. */
void Household::assign_augm_model_pointer_arg(indiv i, const augmID &a, const std::string &arg_name, const double *y)
{
	pImpl->members[i].augm_map.at(a).assign_distr_pointer_arg(arg_name, y);
}

/* Assign scalar range minimum for augmented variable model. */
void Household::assign_augm_model_scalar_range_min(indiv i, const augmID &a, double x)
{
	pImpl->members[i].augm_map.at(a).assign_distr_scalar_range_min(x);
}
/* Assign pointer range minimum for augmented variable model. */
void Household::assign_augm_model_pointer_range_min(indiv i, const augmID &a, const double *y)
{
	pImpl->members[i].augm_map.at(a).assign_distr_pointer_range_min(y);
}

/* Assign scalar range maximum for augmented variable model. */
void Household::assign_augm_model_scalar_range_max(indiv i, const augmID &a, double x)
{
	pImpl->members[i].augm_map.at(a).assign_distr_scalar_range_max(x);
}
/* Assign pointer range maximum for augmented variable model. */
void Household::assign_augm_model_pointer_range_max(indiv i, const augmID &a, const double *y)
{
	pImpl->members[i].augm_map.at(a).assign_distr_pointer_range_max(y);
}

void Household::init_augm_model_reshaping(indiv i, const augmID &a)
{
	pImpl->members[i].augm_map.at(a).init_distr_reshaping();
}

/* Set a value for the augmented variable. */
void Household::set_augm_value(indiv i, const augmID &a, double x)
{
	pImpl->members[i].augm_map.at(a).set_value(x);
}

/* Draw a value from the augmented variable distribution. */
void Household::augm_model_rand(indiv i, const augmID &a, rnd_gen &prng)
{
	pImpl->members[i].augm_map.at(a).rand(prng);
}

/* True if the propose augm move is valid. */
bool Household::propose_valid_augm_move(indiv i, const augmID &a, rnd_gen &prng)
{
	Variable &this_augm = pImpl->members[i].augm_map.at(a);
	this_augm.propose(prng); // propose new value and buffer old

	if (this_augm.valid_move()) // proceed only if valid
		return true;
	else
		this_augm.restore_value(); // reassign buffered value

	return false;
}

/* Restore buffer valies of augmented variable, extra and loglikelihood. */
void Household::restore_augm_value(indiv i, const augmID &a)
{
	pImpl->members[i].augm_map.at(a).restore_value();
	restore_all_extra();
	restore_loglike();
}

/* Assign empty augmented variable. */
void Household::clear_augm(indiv i, const augmID &a) { pImpl->members[i].augm_map[a] = Variable(); }

//===== DIAGNOSTICS =====//

/* Update the values for the Bayesian diagnostics. */
void Household::init_diagnostics(ulong sample_size)
{
	pImpl->dgn_loglike.clear();
	pImpl->dgn_loglike = dbl_vec(sample_size);
}

/* Update the values for the Bayesian diagnostics. */
void Household::update_diagnostics(ulong count) { pImpl->dgn_loglike[count] = pImpl->house_loglike; }

/* Return the loglikelihood sample vector. */
const dbl_vec &Household::get_dgn_loglike() const { return pImpl->dgn_loglike; }

///////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////  Household : extra attributes  //////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

bool Household::is_any_extra()
{
	return !(extra_int_vec.empty() && extra_dbl_vec.empty() && extra_dbl_mtx.empty());
}

/* Deep copies of extra attributes in the buffer */
void Household::buffer_all_extra()
{
	_buffer_extra_int_vec_ = extra_int_vec;
	_buffer_extra_dbl_vec_ = extra_dbl_vec;
	_buffer_extra_dbl_mtx_ = extra_dbl_mtx;
}

/* Swap buffer into  extra attributes */
void Household::restore_all_extra()
{
	extra_int_vec.swap(_buffer_extra_int_vec_);
	extra_dbl_vec.swap(_buffer_extra_dbl_vec_);
	extra_dbl_mtx.swap(_buffer_extra_dbl_mtx_);
}

//----- Extra integer vector -----//
void Household::allocate_extra_integer_vector(const extraID &e) { extra_int_vec[e] = int_vec(size(), 0); }
int_vec &Household::extra_integer_vector(const extraID &e) { return extra_int_vec[e]; }

//----- Extra double vector -----//
void Household::allocate_extra_double_vector(const extraID &e) { extra_dbl_vec[e] = dbl_vec(size(), 0.); }
dbl_vec &Household::extra_double_vector(const extraID &e) { return extra_dbl_vec[e]; }

//----- Extra double matrix -----//
void Household::allocate_extra_double_matrix(const extraID &e) { extra_dbl_mtx[e] = matrix<double>(size(), dbl_vec(size(), 0.)); }
matrix<double> &Household::extra_double_matrix(const extraID &e) { return extra_dbl_mtx[e]; }
