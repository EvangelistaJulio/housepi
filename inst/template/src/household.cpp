/*
 * household.cpp
 *
 * Copyright (C) 2025 Francesco Camaglia, Institut Pasteur
 */

#include "individual.h"
#include "household.h"

//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//

class Household_Impl
{
public:
	//===== CONSTRUCTOR =====//
	Household_Impl() {};
	~Household_Impl() {};

	//===== ATTRIBUTES =====//
	std::vector<Individual> members = {};

	std::string ID = NULL_STRING; // user ID
	size_t index = -1;			// index of household
	ushort house_size = 0;		// number of house members

	double house_loglike = NULL_DOUBLE; // current loglikelihood value
	double buffer = NULL_DOUBLE;		// buffer loglikelihood value

	ushort n_infectious = 0;
	ushort n_susceptible = 0;
	ushort n_missing_status = 0;

	// wheter members correspond to a given status or not
	std::vector<bool> _is_infectious_ = {};
	std::vector<bool> _is_susceptible_ = {};
	std::vector<bool> _is_missing_status_ = {};

	void assign_individual(const std::string &ID,
						   int status,
						   const dataIDmap<double> &data);
	bool change_status(indiv i, int new_status);

	// diagnostics
	dbl_vec dgn_loglike = {};

	// contact rate pointers
	matrix<const double *> contact_rates;

	// list with ONLY indexes of members with a given status
	std::set<indiv> infectious_set = {};
	std::set<indiv> susceptible_set = {};
	std::set<indiv> missing_status_set = {};

	//===== EXTRA ATTRIBUTES =====//

	//----- Extra Attributes -----//
	extraIDmap<int> extra_integer = {};
	extraIDmap<double> extra_double = {};
	extraIDmap<dbl_vec> extra_vec = {};
	extraIDmap<matrix<double>> extra_mtx = {};

	//----- Buffer EXTRA Attributes -----//
	extraIDmap<int> _buffer_extra_integer_ = {};
	extraIDmap<double> _buffer_extra_double_ = {};
	extraIDmap<dbl_vec> _buffer_extra_vec_ = {};
	extraIDmap<matrix<double>> _buffer_extra_mtx_ = {};

	void buffer_all_extra();
	void restore_all_extra();

	// UPDATE: assertion vs deafult
	// pointer-to-method is a drawback for all getters
	// NOTE: assertion should be done elsewhere and knocked off elsewhere
	// NOTE : in any case temporary if dropping maps

	void disable_assert();

	//===== Setter (auxiliary methods) =====//
	void (Household_Impl::*_func_set_augm_)(indiv, const augmID &, double) = &Household_Impl::_set_augm_assert_;
	void _set_augm_assert_(indiv i, const augmID &a, double xnew);
	void _set_augm_default_(indiv i, const augmID &a, double xnew);

	double get_data(indiv i, const dataID &d) const;
	double get_augm(indiv i, const augmID &a) const;
	double get_rate(indiv i, const rateID &r) const;

	//===== Getter (auxiliary methods) =====//
	double (Household_Impl::*_func_get_augm_)(indiv, const augmID &) const = &Household_Impl::_get_augm_assert_;
	double _get_augm_assert_(indiv i, const augmID &a) const;
	double _get_augm_default_(indiv i, const augmID &a) const;

	double (Household_Impl::*_func_get_data_)(indiv, const dataID &) const = &Household_Impl::_get_data_assert_;
	double _get_data_assert_(indiv i, const dataID &d) const;
	double _get_data_default_(indiv i, const dataID &d) const;

	double (Household_Impl::*_func_get_rate_)(indiv, const rateID &) const = &Household_Impl::_get_rate_assert_;
	double _get_rate_assert_(indiv i, const rateID &d) const;
	double _get_rate_default_(indiv i, const rateID &d) const;
};

//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//

/* Link an individual to the household. */
void Household_Impl::assign_individual(const std::string &ID,
									   int status,
									   const dataIDmap<double> &data)
{
	//----- update status lists within the household -----//
	if (status >= INFECTED_STATUS)
	{
		++n_infectious;
		infectious_set.insert(house_size);

		_is_infectious_.push_back(true);
		_is_susceptible_.push_back(false);
		_is_missing_status_.push_back(false);
	}
	else if (status == HEALTHY_STATUS)
	{
		++n_susceptible;
		susceptible_set.insert(house_size);

		_is_infectious_.push_back(false);
		_is_susceptible_.push_back(true);
		_is_missing_status_.push_back(false);
	}
	else
	{
		if (status < MISSING_STATUS)
		{
			std::ostringstream oss;
			oss << "Status '" << status << "' of individual " << ID << " interpreted as missing."
				<< " Use (=)" << HEALTHY_STATUS << " for susceptible, (>=)" << INFECTED_STATUS << " for infectious,"
				<< " and (=) " << MISSING_STATUS << " for missing status.";
			tools::warning(oss.str());
		}
		status = MISSING_STATUS;

		++n_missing_status;
		missing_status_set.insert(house_size);

		_is_infectious_.push_back(false);
		_is_susceptible_.push_back(false);
		_is_missing_status_.push_back(true);
	}

	//----- update household size -----//
	++house_size;

	//----- emplace new individual -----//

	members.emplace_back(ID, status, data);
}

/* Change_status of individual i to new status (without controls).
	NOTE : this introduces augmented variables -> it requires RJMCMC in case of inference. */
bool Household_Impl::change_status(indiv i, int new_status)
{
	Individual &this_indiv = members[i];
	int old_status = this_indiv.status;

	if (new_status == old_status)
	{
		return false;
	}
	else
	{
		if (old_status >= INFECTED_STATUS)
		{
			--n_infectious;
			infectious_set.erase(i);
			_is_infectious_[i] = false;
		}
		else if (old_status == HEALTHY_STATUS)
		{
			--n_susceptible;
			susceptible_set.erase(i);
			_is_susceptible_[i] = false;
		}
		else if (old_status == MISSING_STATUS)
		{
			--n_missing_status;
			missing_status_set.erase(i);
			_is_missing_status_[i] = false;
		}
		else
		{
			throw std::runtime_error(
				"Unknown status '" + std::to_string(old_status) + "' for individual " + this_indiv.ID + ".");
		}

		if (new_status >= INFECTED_STATUS)
		{
			++n_infectious;
			infectious_set.insert(i);
			_is_infectious_[i] = true;
		}
		else if (new_status == HEALTHY_STATUS)
		{
			++n_susceptible;
			susceptible_set.insert(i);
			_is_susceptible_[i] = true;
		}
		else if (new_status == MISSING_STATUS)
		{
			++n_missing_status;
			missing_status_set.insert(i);
			_is_missing_status_[i] = true;
		}
		else
		{
			throw std::runtime_error(
				"Unknown new status '" + std::to_string(new_status) + "' for individual " + this_indiv.ID + ".");
		}

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

double Household_Impl::get_data(indiv i, const dataID &d) const { return (this->*_func_get_data_)(i, d); }
double Household_Impl::get_augm(indiv i, const augmID &a) const { return (this->*_func_get_augm_)(i, a); }
double Household_Impl::get_rate(indiv i, const rateID &r) const { return (this->*_func_get_rate_)(i, r); }

/* Disable data and augmented variable assertions. */
void Household_Impl::disable_assert()
{
	this->_func_set_augm_ = &Household_Impl::_set_augm_default_;
	this->_func_get_augm_ = &Household_Impl::_get_augm_default_;
	this->_func_get_data_ = &Household_Impl::_get_data_default_;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////  Household_Impl : functional methods  //////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/* Default method to get the value of data. */
double Household_Impl::_get_data_default_(indiv i, const dataID &d) const
{
	return members[i].data.at(d);
}

/* Safe method to get the value of data with previous assertion. */
double Household_Impl::_get_data_assert_(indiv i, const dataID &d) const
{
	if (members[i].assert_data(d))
		return _get_data_default_(i, d);
	else
		return NULL_DOUBLE;
}

/* Default method to get the value of an augmented variable. */
double Household_Impl::_get_augm_default_(indiv i, const augmID &a) const
{
	return members[i].augm_map.at(a).get_value();
}

/* Safe method to get the value of an augmented variable with previous assertion. */
double Household_Impl::_get_augm_assert_(indiv i, const augmID &a) const
{
	if (members[i].assert_augm(a))
		return _get_augm_default_(i, a);
	else
		return NULL_DOUBLE;
}

/* Default method to get the value of rate. */
double Household_Impl::_get_rate_default_(indiv i, const rateID &r) const
{
	return members[i].rate_map.at(r).get();
}

/* Safe method to get the rate with previous assertion. */
double Household_Impl::_get_rate_assert_(indiv i, const rateID &r) const
{
	if (members[i].assert_rate(r))
		return _get_rate_default_(i, r);
	else
		return NULL_DOUBLE;
}

/* Default method to set the value of an augmented variable. */
void Household_Impl::_set_augm_default_(indiv i, const augmID &a, double xnew)
{
	members[i].augm_map.at(a).set_value(xnew);
}

/* Safe method to set the value of an augmented variable with previous assertion. */
void Household_Impl::_set_augm_assert_(indiv i, const augmID &a, double xnew)
{
	if (members[i].assert_augm(a))
		_set_augm_default_(i, a, xnew);
}

//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//

/* Deep copies of extra attributes in the buffer */
void Household_Impl::buffer_all_extra()
{
	_buffer_extra_integer_ = extra_integer;
	_buffer_extra_double_ = extra_double;
	_buffer_extra_vec_ = extra_vec;
	_buffer_extra_mtx_ = extra_mtx;
}

/* Swap buffer into  extra attributes */
void Household_Impl::restore_all_extra()
{
	extra_integer.swap(_buffer_extra_integer_);
	extra_double.swap(_buffer_extra_double_);
	extra_vec.swap(_buffer_extra_vec_);
	extra_mtx.swap(_buffer_extra_mtx_);
}

//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//

//===== CONSTRUCTOR =====//
Household::Household() : pImpl(std::make_unique<Household_Impl>()) {};
Household::Household(const std::string &ID, size_t H) : pImpl(std::make_unique<Household_Impl>())
{
	pImpl->ID = ID;
	pImpl->index = H;
}
Household::Household(Household &&other) noexcept = default;
Household &Household::operator=(Household &&other) noexcept = default;
Household::~Household() = default;

void Household::assign_individual(const std::string &ID,
								  int status,
								  const dataIDmap<double> &data)
{
	pImpl->assign_individual(ID, status, data);
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
const std::string &Household::house_ID() const { return pImpl->ID; }
const std::string &Household::ID(indiv i) const { return pImpl->members[i].ID; }
size_t Household::index() const { return pImpl->index; }
ushort Household::size() const { return pImpl->house_size; }

int Household::status(indiv i) const { return pImpl->members[i].status; }
double Household::contact_rate(indiv from, indiv to) const { return *pImpl->contact_rates[from][to]; }

// UPDATE : assertion is a drawback
// NOTE : temporary before change into vectors
double Household::augm(indiv i, const augmID &a) const { return pImpl->get_augm(i, a); }
double Household::data(indiv i, const dataID &d) const { return pImpl->get_data(i, d); }
double Household::group_rate(indiv i, const rateID &r) const { return pImpl->get_rate(i, r); }

double Household::augm_model_logpdf(indiv i, const augmID &a) const { return pImpl->members[i].augm_map.at(a).distr_logpdf(); }
double Household::augm_model_logratio(indiv i, const augmID &a) const { return pImpl->members[i].augm_map.at(a).distr_logratio(); }
double Household::augm_propose_logratio(indiv i, const augmID &a) const { return pImpl->members[i].augm_map.at(a).propose_logratio(); }
double Household::augm_move_logratio(indiv i, const augmID &a) const { return pImpl->members[i].augm_map.at(a).get_move_logratio(); }

ushort Household::n_infectious() const { return pImpl->n_infectious; }
ushort Household::n_susceptible() const { return pImpl->n_susceptible; }
ushort Household::n_missing_status() const { return pImpl->n_missing_status; }

const std::set<indiv> &Household::infectious_set() const { return pImpl->infectious_set; }
const std::set<indiv> &Household::susceptible_set() const { return pImpl->susceptible_set; }
const std::set<indiv> &Household::missing_status_set() const { return pImpl->missing_status_set; }

double Household::get_loglike() const { return pImpl->house_loglike; }
double Household::get_delta_loglike() const { return pImpl->house_loglike - pImpl->buffer; }

//===== SETTER =====//

/* Disable data and augmented variable assertions. */
void Household::disable_assert() { pImpl->disable_assert(); }

void Household::assign_all_group_rates(indiv i, const rateIDmapvec<const double *> &x)
{
	pImpl->members[i].assign_all_group_rates(x);
}

/* Assign definitions for the augmented variables. */
void Household::assign_augm(indiv i, const augmID &a, const std::string &name, double value,
							bool is_discrete, const std::string &model_name,
							const std::string &rule_name, double rule_scale, double rule_from)
{
	Base_Variable &this_augm = pImpl->members[i].augm_map.at(a);
	this_augm.init(name, value, is_discrete);
	this_augm.assign_distribution(model_name);
	this_augm.assign_proposer(rule_name, rule_scale, rule_from);
}

/* Assign scalar argument for augmented variable model. */
void Household::assign_augm_model_scalar_arg(indiv i, const augmID &a,
											 const std::string &arg_name,
											 double x)
{
	pImpl->members[i].augm_map.at(a).assign_distr_scalar_arg(arg_name, x);
}

/* Assign pointer to argument for augmented variable model. */
void Household::assign_augm_model_pointer_arg(indiv i, const augmID &a,
											  const std::string &arg_name,
											  const double *y)
{
	pImpl->members[i].augm_map.at(a).assign_distr_pointer_arg(arg_name, y);
}

/* Assign range for augmented variable model. */
void Household::assign_augm_model_range(indiv i, const augmID &a,
										double lower, double upper)
{
	pImpl->members[i].augm_map.at(a).assign_distr_range(lower, upper);
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

/**/
bool Household::propose_valid_augm_move(indiv i, const augmID &a, rnd_gen &prng)
{
	Base_Variable &this_augm = pImpl->members[i].augm_map.at(a);

	this_augm.propose(prng); // propose new value and buffer old

	if (this_augm.valid_move()) // proceed only if valid
		return true;
	else
		this_augm.restore_value(); // reassign buffered value

	return false;
}

void Household::restore_augm_value(indiv i, const augmID &a)
{
	pImpl->members[i].augm_map.at(a).restore_value();
}

void Household::clear_augm(indiv i, const augmID &a)
{
	pImpl->members[i].augm_map[a] = Base_Variable();
}

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

//===== EXTRA ATTRIBUTES =====//
void Household::_buffer_all_extra() { pImpl->buffer_all_extra(); }
void Household::_restore_all_extra() { pImpl->restore_all_extra(); }

//----- Single integers -----//
void Household::allocate_extra_integer(const std::vector<extraID> &names, int xinit)
{
	for (const extraID &e : names)
		pImpl->extra_integer[e] = xinit;
}
void Household::set_extra_integer(const extraID &e, int xnew)
{
	pImpl->extra_integer.at(e) = xnew;
}

int Household::extra_integer(const extraID &e) const
{
	return pImpl->extra_integer.at(e);
}

void Household::allocate_extra_double(const std::vector<extraID> &names, double xinit)
{
	for (const extraID &e : names)
		pImpl->extra_double[e] = xinit;
}
void Household::set_extra_double(const extraID &e, double xnew)
{
	pImpl->extra_double.at(e) = xnew;
}

double Household::extra_double(const extraID &e) const
{
	return pImpl->extra_double.at(e);
}

void Household::allocate_extra_vector(const std::vector<extraID> &names, double xinit)
{
	for (const extraID &e : names)
		pImpl->extra_vec[e] = dbl_vec(size(), xinit);
}
void Household::set_extra_vector(const extraID &e, size_t i, double xnew)
{
	pImpl->extra_vec.at(e)[i] = xnew;
}
double Household::extra_vector(const extraID &e, size_t i) const
{
	return pImpl->extra_vec.at(e)[i];
}
const dbl_vec &Household::extra_vector(const extraID &e) const
{
	return pImpl->extra_vec.at(e);
}

void Household::allocate_extra_matrix(const std::vector<extraID> &names, double xinit)
{
	for (const extraID &e : names)
		pImpl->extra_mtx[e] = matrix<double>(size(), dbl_vec(size(), xinit));
}

void Household::set_extra_matrix(const extraID &e, size_t from, size_t to, double xnew)
{
	pImpl->extra_mtx.at(e)[from][to] = xnew;
}

double Household::extra_matrix(const extraID &e, size_t from, size_t to) const
{
	return pImpl->extra_mtx.at(e)[from][to];
}

const matrix<double> &Household::extra_matrix(const extraID &e) const
{
	return pImpl->extra_mtx.at(e);
}
