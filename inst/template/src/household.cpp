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
	std::vector<Individual> indiv = {};

	std::string ID = UNDEFINED; // user ID
	int index = -1;				// index of household
	int house_size = 0;			// number of house members

	double house_loglike = -INF_DOUBLE; // current loglikelihood value
	double buffer = -INF_DOUBLE;		// buffer loglikelihood value

	int n_infectious = 0;
	int n_susceptible = 0;
	int n_missing_status = 0;

	// wheter members correspond to a given status or not
	std::vector<bool> _is_infectious_ = {};
	std::vector<bool> _is_susceptible_ = {};
	std::vector<bool> _is_missing_status_ = {};

	void assign_individual(const std::string &ID, int status, const dictargs &data);
	bool change_status(int i, int new_status);

	// diagnostics
	ulong dgn_count = 0;
	std::vector<double> dgn_loglike = {};

	// contact rate pointers
	std::vector<std::vector<const double *>> contact_rates;
	std::vector<std::vector<std::string>> _which_rate_;

	// list with ONLY indexes of members with a given status
	std::set<int> infectious_set = {};
	std::set<int> susceptible_set = {};
	std::set<int> missing_status_set = {};

	//===== EXTRA ATTRIBUTES =====//

	//----- Extra Attributes -----//
	std::unordered_map<std::string, int> extra_int = {};
	std::unordered_map<std::string, double> extra_double = {};
	std::unordered_map<std::string, std::vector<double>> extra_vec = {};
	std::unordered_map<std::string, squarematrix> extra_mtx = {};

	//----- Buffer EXTRA Attributes -----//
	std::unordered_map<std::string, int> _buffer_extra_int_ = {};
	std::unordered_map<std::string, std::vector<double>> _buffer_extra_vec_ = {};
	std::unordered_map<std::string, squarematrix> _buffer_extra_mtx_ = {};

	// deep copies of extra attributes in the buffer
	void buffer_all_extra()
	{
		_buffer_extra_int_ = extra_int;
		_buffer_extra_vec_ = extra_vec;
		_buffer_extra_mtx_ = extra_mtx;
	}

	// swap buffer into  extra attributes
	void restore_all_extra()
	{
		extra_int.swap(_buffer_extra_int_);
		extra_vec.swap(_buffer_extra_vec_);
		extra_mtx.swap(_buffer_extra_mtx_);
	}
};

/* Link an individual to the household */
void Household_Impl::assign_individual(const std::string &ID, int status, const dictargs &data)
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
	indiv.emplace_back(ID, status, data);
}

/* Change_status of individual i to new status (without controls).
	NOTE : this introduces augmented variables -> it requires RJMCMC in case of inference. */
bool Household_Impl::change_status(int i, int new_status)
{
	Individual &this_indiv = indiv[i];
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
Household::Household(const std::string &ID, int H) : pImpl(std::make_unique<Household_Impl>())
{
	pImpl->ID = ID;
	pImpl->index = H;
}
Household::Household(Household &&other) noexcept = default;
Household &Household::operator=(Household &&other) noexcept = default;
Household::~Household() = default;

void Household::assign_individual(const std::string &ID, int status, const dictargs &data)
{
	pImpl->assign_individual(ID, status, data);
}

/* Change_status of individual i to new status (without controls).
	NOTE : this introduces augmented variables -> it requires RJMCMC in case of inference. */
bool Household::change_status(int i, int new_status) { return pImpl->change_status(i, new_status); }

/* Initializes matrix of contact rates. */
void Household::init_contact_matrix()
{
	pImpl->_which_rate_.resize(size(), std::vector<std::string>(size(), UNDEFINED));
	pImpl->contact_rates.resize(size(), std::vector<const double *>(size(), &REFERENCE_VALUE));
}

void Household::assign_contact_rate(int from, int to, const std::string &parID) { pImpl->_which_rate_[from][to] = parID; }
void Household::set_contact_rate(int from, int to, const double *x) { pImpl->contact_rates[from][to] = x; }
void Household::set_loglike(double x)
{
	pImpl->buffer = pImpl->house_loglike;
	pImpl->house_loglike = x;
}
void Household::restore_loglike() { pImpl->house_loglike = pImpl->buffer; }

//===== GETTER =====//
const std::string &Household::house_ID() const { return pImpl->ID; }
const std::string &Household::ID(int i) const { return pImpl->indiv[i].get_ID(); }
int Household::index() const { return pImpl->index; }
int Household::size() const { return pImpl->house_size; }

const std::unordered_map<std::string, std::string> &Household::get_group_rate_map(int i) const { return pImpl->indiv[i].get_group_rate_map(); }
double Household::group_rate(int i, const std::string &rateID) const { return pImpl->indiv[i].get_group_rate(rateID); }

const std::string &Household::get_contact_rate_map(int from, int to) const { return pImpl->_which_rate_[from][to]; }
double Household::contact_rate(int from, int to) const { return *pImpl->contact_rates[from][to]; }

int Household::status(int i) const { return pImpl->indiv[i].get_status(); }
double Household::augm(int i, const std::string &augmID) const { return pImpl->indiv[i].get_augm(augmID); }
double Household::data(int i, const std::string &dataID) const { return pImpl->indiv[i].get_data(dataID); }
const double *Household::data_pointer(int i, const std::string &dataID) const { return pImpl->indiv[i].get_data_pointer(dataID); }

int Household::n_infectious() const { return pImpl->n_infectious; }
int Household::n_susceptible() const { return pImpl->n_susceptible; }
int Household::n_missing_status() const { return pImpl->n_missing_status; }

const std::set<int> &Household::infectious_set() const { return pImpl->infectious_set; }
const std::set<int> &Household::susceptible_set() const { return pImpl->susceptible_set; }
const std::set<int> &Household::missing_status_set() const { return pImpl->missing_status_set; }

double Household::get_loglike() const { return pImpl->house_loglike; }
double Household::get_delta_loglike() const { return pImpl->house_loglike - pImpl->buffer; }

double Household::augm_model_logpdf(int i, const std::string &augmID) const { return pImpl->indiv[i].augm_model_logpdf(augmID); }
double Household::augm_model_logratio(int i, const std::string &augmID) const { return pImpl->indiv[i].augm_model_logratio(augmID); }
double Household::augm_propose_logratio(int i, const std::string &augmID) const { return pImpl->indiv[i].augm_propose_logratio(augmID); }
double Household::augm_move_logratio(int i, const std::string &augmID) const { return pImpl->indiv[i].augm_move_logratio(augmID); }

//===== SETTER =====//
void Household::disable_assert()
{
	for (int i = 0; i < size(); ++i)
	{
		pImpl->indiv[i].disable_assert();
	}
}

void Household::assign_augm(int i, const std::string &augmID, double value,
							const std::string &model_name, dictargs_ptr model_args_ptr, dictargs range,
							const std::string &rule_name, dictargs rule_args, double rule_from)
{
	pImpl->indiv[i].assign_augm(augmID, value,
								model_name, model_args_ptr, range,
								rule_name, rule_args, rule_from);
}
void Household::set_augm_value(int i, const std::string &augmID, double x) { pImpl->indiv[i].set_augm_value(augmID, x); }
void Household::augm_model_rand(int i, const std::string &augmID, sys_rnd_gen &prng) { pImpl->indiv[i].augm_model_rand(augmID, prng); }
bool Household::propose_valid_augm_move(int i, const std::string &augmID, sys_rnd_gen &prng) { return pImpl->indiv[i].propose_valid_augm_move(augmID, prng); }
void Household::restore_augm_value(int i, const std::string &augmID) { pImpl->indiv[i].restore_augm_value(augmID); }

void Household::assign_group_rate(int i, const std::string &rateID, const std::string &parID) { pImpl->indiv[i].assign_group_rate(rateID, parID); }
void Household::set_rate(int i, const std::string &rateID, const double *ptr) { pImpl->indiv[i].set_rate(rateID, ptr); }
void Household::clear_augm(int i, const std::string &augmID) { pImpl->indiv[i].clear_augm(augmID); }

//===== DIAGNOSTICS =====//

/* Update the values for the Bayesian diagnostics. */
void Household::init_diagnostics(ulong sample_size)
{
	pImpl->dgn_count = 0;
	pImpl->dgn_loglike.clear();
	pImpl->dgn_loglike = std::vector<double>(sample_size);
}

/* Update the values for the Bayesian diagnostics. */
void Household::update_diagnostics()
{
	pImpl->dgn_loglike[pImpl->dgn_count++] = pImpl->house_loglike;
}

/* Return the loglikelihood sample vector. */
const std::vector<double> &Household::get_dgn_loglike() const
{
	return pImpl->dgn_loglike;
}

//===== EXTRA ATTRIBUTES =====//
void Household::_buffer_all_extra() { pImpl->buffer_all_extra(); }
void Household::_restore_all_extra() { pImpl->restore_all_extra(); }

void Household::allocate_extra_int(const std::vector<std::string> names, int init_value)
{
	for (const std::string &intID : names)
	{
		pImpl->extra_int.emplace(intID, init_value);
	}
}
void Household::set_extra_int(const std::string &intID, int new_value)
{
	pImpl->extra_int.at(intID) = new_value;
}

int Household::extra_int(const std::string &intID) const
{
	return pImpl->extra_int.at(intID);
}


void Household::allocate_extra_double(const std::vector<std::string> names, double init_value)
{
	for (const std::string &dblID : names)
	{
		pImpl->extra_double.emplace(dblID, init_value);
	}
}
void Household::set_extra_double(const std::string &dblID, double new_value)
{
	pImpl->extra_double.at(dblID) = new_value;
}

double Household::extra_double(const std::string &dblID) const
{
	return pImpl->extra_double.at(dblID);
}

void Household::allocate_extra_vector(const std::vector<std::string> names, double init_value)
{
	for (const std::string &vecID : names)
	{
		pImpl->extra_vec.emplace(vecID, std::vector<double>(size(), init_value));
	}
}
void Household::set_extra_vector(const std::string &vecID, size_t i, double new_value)
{
	pImpl->extra_vec.at(vecID)[i] = new_value;
}
double Household::extra_vector(const std::string &vecID, size_t i) const
{
	return pImpl->extra_vec.at(vecID)[i];
}
const std::vector<double> &Household::extra_vector(const std::string &vecID) const
{
	return pImpl->extra_vec.at(vecID);
}

void Household::allocate_extra_matrix(const std::vector<std::string> names, double init_value)
{
	for (const std::string &mtxID : names)
	{
		pImpl->extra_mtx.emplace(mtxID, std::vector<std::vector<double>>(size(), std::vector<double>(size(), init_value)));
	}
}

void Household::set_extra_matrix(const std::string &mtxID, size_t from, size_t to, double new_value)
{
	pImpl->extra_mtx.at(mtxID)[from][to] = new_value;
}

double Household::extra_matrix(const std::string &mtxID, size_t from, size_t to) const
{
	return pImpl->extra_mtx.at(mtxID)[from][to];
}

const squarematrix &Household::extra_matrix(const std::string &mtxID) const
{
	return pImpl->extra_mtx.at(mtxID);
}
