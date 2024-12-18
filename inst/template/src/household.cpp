/*
 * household.cpp
 *
 * Copyright (C) 2024 Francesco Camaglia, Institut Pasteur
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
	std::vector<Individual> indiv_vector = {};

	std::string ID = UNDEFINED; // user ID
	size_t index = -1;			// index of household
	ushort house_size = 0;		// number of house members

	double house_loglike = -INF_DOUBLE; // current loglikelihood value
	double buffer = -INF_DOUBLE;		// buffer loglikelihood value

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
	ulong dgn_count = 0;
	std::vector<double> dgn_loglike = {};

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
	extraIDmap<std::vector<double>> extra_vec = {};
	extraIDmap<matrix<double>> extra_mtx = {};

	//----- Buffer EXTRA Attributes -----//
	extraIDmap<int> _buffer_extra_integer_ = {};
	extraIDmap<double> _buffer_extra_double_ = {};
	extraIDmap<std::vector<double>> _buffer_extra_vec_ = {};
	extraIDmap<matrix<double>> _buffer_extra_mtx_ = {};

	// deep copies of extra attributes in the buffer
	void buffer_all_extra()
	{
		_buffer_extra_integer_ = extra_integer;
		_buffer_extra_double_ = extra_double;
		_buffer_extra_vec_ = extra_vec;
		_buffer_extra_mtx_ = extra_mtx;
	}

	// swap buffer into  extra attributes
	void restore_all_extra()
	{
		extra_integer.swap(_buffer_extra_integer_);
		extra_double.swap(_buffer_extra_double_);
		extra_vec.swap(_buffer_extra_vec_);
		extra_mtx.swap(_buffer_extra_mtx_);
	}
};

/* Link an individual to the household */
void Household_Impl::assign_individual(const std::string &ID,
									   int status,
									   const dataIDmap<double> &data)
{
	//----- update status lists within the household -----//
	if (status >= INFECTIOUS)
	{
		++n_infectious;
		infectious_set.insert(house_size);

		_is_infectious_.push_back(true);
		_is_susceptible_.push_back(false);
		_is_missing_status_.push_back(false);
	}
	else if (status == SUSCEPTIBLE)
	{
		++n_susceptible;
		susceptible_set.insert(house_size);

		_is_infectious_.push_back(false);
		_is_susceptible_.push_back(true);
		_is_missing_status_.push_back(false);
	}
	else if (status == MISSING_STATUS)
	{
		++n_missing_status;
		missing_status_set.insert(house_size);

		_is_infectious_.push_back(false);
		_is_susceptible_.push_back(false);
		_is_missing_status_.push_back(true);
	}
	else
	{
		throw std::invalid_argument(
			"Unknown status '" + std::to_string(status) + "' for individual " + ID + ".");
	}

	//----- update household size -----//
	++house_size;

	//----- emplace new individual -----//
	indiv_vector.emplace_back(ID, status, data);
}

/* Change_status of individual i to new status (without controls).
	NOTE : this introduces augmented variables -> it requires RJMCMC in case of inference. */
bool Household_Impl::change_status(indiv i, int new_status)
{
	Individual &this_indiv = indiv_vector[i];
	int old_status = this_indiv.get_status();

	if (new_status == old_status)
	{
		return false;
	}
	else
	{
		if (old_status >= INFECTIOUS)
		{
			--n_infectious;
			infectious_set.erase(i);
			_is_infectious_[i] = false;
		}
		else if (old_status == SUSCEPTIBLE)
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
			std::ostringstream oss;
			oss << "Unknown status '"
				<< old_status << "' for individual " << this_indiv.get_ID() << ".";
			throw std::runtime_error(oss.str());
		}

		if (new_status >= INFECTIOUS)
		{
			++n_infectious;
			infectious_set.insert(i);
			_is_infectious_[i] = true;
		}
		else if (new_status == SUSCEPTIBLE)
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
			std::ostringstream oss;
			oss << "While changing status: unknown new status '"
				<< new_status << "' for individual " << this_indiv.get_ID() << ".";
			throw std::runtime_error(oss.str());
		}

		this_indiv.set_status(new_status);
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
const std::string &Household::ID(indiv i) const { return pImpl->indiv_vector[i].get_ID(); }
size_t Household::index() const { return pImpl->index; }
ushort Household::size() const { return pImpl->house_size; }

double Household::group_rate(indiv i, const rateID &r) const { return pImpl->indiv_vector[i].get_group_rate(r); }
double Household::contact_rate(indiv from, indiv to) const { return *pImpl->contact_rates[from][to]; }

int Household::status(indiv i) const { return pImpl->indiv_vector[i].get_status(); }
double Household::augm(indiv i, const augmID &a) const { return pImpl->indiv_vector[i].get_augm(a); }
double Household::data(indiv i, const dataID &d) const { return pImpl->indiv_vector[i].get_data(d); }
const double *Household::data_pointer(indiv i, const dataID &d) const { return pImpl->indiv_vector[i].get_data_pointer(d); }

ushort Household::n_infectious() const { return pImpl->n_infectious; }
ushort Household::n_susceptible() const { return pImpl->n_susceptible; }
ushort Household::n_missing_status() const { return pImpl->n_missing_status; }

const std::set<indiv> &Household::infectious_set() const { return pImpl->infectious_set; }
const std::set<indiv> &Household::susceptible_set() const { return pImpl->susceptible_set; }
const std::set<indiv> &Household::missing_status_set() const { return pImpl->missing_status_set; }

double Household::get_loglike() const { return pImpl->house_loglike; }
double Household::get_delta_loglike() const { return pImpl->house_loglike - pImpl->buffer; }

double Household::augm_model_logpdf(indiv i, const augmID &a) const { return pImpl->indiv_vector[i].augm_model_logpdf(a); }
double Household::augm_model_logratio(indiv i, const augmID &a) const { return pImpl->indiv_vector[i].augm_model_logratio(a); }
double Household::augm_propose_logratio(indiv i, const augmID &a) const { return pImpl->indiv_vector[i].augm_propose_logratio(a); }
double Household::augm_move_logratio(indiv i, const augmID &a) const { return pImpl->indiv_vector[i].augm_move_logratio(a); }

//===== SETTER =====//

void Household::disable_assert()
{
	for (indiv i = 0; i < size(); ++i)
		pImpl->indiv_vector[i].disable_assert();
}

void Household::assign_all_group_rates(indiv i, const rateIDmapvec<const double *> &x) { pImpl->indiv_vector[i].assign_all_group_rates(x); }

void Household::assign_augm(indiv i, const augmID &a, const std::string &name, double value, bool is_discrete,
							const std::string &model_name, const dict_args_ptr &model_args_ptr, const dict_args &range,
							const std::string &rule_name, double rule_scale, double rule_from)
{
	pImpl->indiv_vector[i].assign_augm(a, name, value, is_discrete,
									   model_name, model_args_ptr, range,
									   rule_name, rule_scale, rule_from);
}
void Household::set_augm_value(indiv i, const augmID &a, double x) { pImpl->indiv_vector[i].set_augm_value(a, x); }
void Household::augm_model_rand(indiv i, const augmID &a, rnd_gen &prng) { pImpl->indiv_vector[i].augm_model_rand(a, prng); }
bool Household::propose_valid_augm_move(indiv i, const augmID &a, rnd_gen &prng) { return pImpl->indiv_vector[i].propose_valid_augm_move(a, prng); }
void Household::restore_augm_value(indiv i, const augmID &a) { pImpl->indiv_vector[i].restore_augm_value(a); }
void Household::clear_augm(indiv i, const augmID &a) { pImpl->indiv_vector[i].clear_augm(a); }

//===== DIAGNOSTICS =====//

/* Update the values for the Bayesian diagnostics. */
void Household::init_diagnostics(ulong sample_size)
{
	pImpl->dgn_count = 0;
	pImpl->dgn_loglike.clear();
	pImpl->dgn_loglike = std::vector<double>(sample_size);
}

/* Update the values for the Bayesian diagnostics. */
void Household::update_diagnostics() { pImpl->dgn_loglike[pImpl->dgn_count++] = pImpl->house_loglike; }

/* Return the loglikelihood sample vector. */
const std::vector<double> &Household::get_dgn_loglike() const { return pImpl->dgn_loglike; }

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
		pImpl->extra_vec[e] = std::vector<double>(size(), xinit);
}
void Household::set_extra_vector(const extraID &e, size_t i, double xnew)
{
	pImpl->extra_vec.at(e)[i] = xnew;
}
double Household::extra_vector(const extraID &e, size_t i) const
{
	return pImpl->extra_vec.at(e)[i];
}
const std::vector<double> &Household::extra_vector(const extraID &e) const
{
	return pImpl->extra_vec.at(e);
}

void Household::allocate_extra_matrix(const std::vector<extraID> &names, double xinit)
{
	for (const extraID &e : names)
		pImpl->extra_mtx[e] = matrix<double>(size(), std::vector<double>(size(), xinit));
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
