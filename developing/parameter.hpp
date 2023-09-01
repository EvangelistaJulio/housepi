// parameter.hpp //

#ifndef PARAMETER__H
#define PARAMETER__H

#include <string>
#include <vector>
#include <unordered_map>
#include <cmath>
#include <random>
#include <iostream>

#include "utils.hpp"

///////////////////////////////////////////////////////////////////////
// ////////////////// Metropolis_Hastings : class ////////////////// //
///////////////////////////////////////////////////////////////////////

// FIXME : certain argument must be equal to 0.

class Metropolis_Hastings
{
protected:
	// acceptation
	int num_accept = 0;
	int num_propos = 0;
	double accept_rate = 0.0;

	double _adapter_(int step, int num_adapt_iter);

public:
	// Constructors
	Metropolis_Hastings() = default;
	~Metropolis_Hastings() {}

	double get_accept_rate() const { return accept_rate; };
	void update_accept_rate(bool is_accepted, int step);
	void set_num_accept_pp() { num_accept++; };
	void set_num_propos_pp() { num_propos++; };
	void reset_moves();

	void set_range(kwargs new_range);
	bool valid(double new_value) { return (this->*_valid_func_)(new_value); }

	// --- functions (virtual) --- //
	virtual double propose(std::mt19937_64 &rng, double old_value)
	{
		std::uniform_real_distribution<double> distr_gen(-0.5, 0.5);
		return old_value + distr_gen(rng);
	}
	virtual double log_ratio() { return 0.0; } // NOTE : the options at the moment do not depend on old/new values
	virtual void adapt(int step, int num_adapt_iter)
	{
		(void)step;
		(void)num_adapt_iter;
	};

private:
	// --- Validation : functional methods --- //
	double _range_low_ = -INFINITY;
	double _range_up_ = INFINITY;
	bool (Metropolis_Hastings::*_valid_func_)(double new_value);

	bool _valid_dummy_(double new_value)
	{
		(void)new_value;
		return true;
	}
	bool _valid_range_low_(double new_value) { return (new_value > _range_low_); }
	bool _valid_range_up_(double new_value) { return (new_value < _range_up_); }
	bool _valid_range_both_(double new_value) { return _valid_range_low_(new_value) && _valid_range_up_(new_value); }
};

///////////////////////////////////////////////////////////////////////
// ///////////////////// Normal_Proposal : class /////////////////// //
///////////////////////////////////////////////////////////////////////

class Normal_Proposal : public Metropolis_Hastings
{
private:
	Normal_Distr distrib;
	std::vector<std::string> args_to_adapt = {"sd"};

public:
	// Constructors
	Normal_Proposal(kwargs new_args = {{"mean", (double)0.0}, {"sd", (double)1.0}}) : distrib(new_args){};
	~Normal_Proposal() {}

	// --- functions (override) --- //
	double propose(std::mt19937_64 &rng, double old_value) override { return old_value + distrib.rand(rng); }
	double log_ratio() override { return 0.0; }
	void adapt(int step, int num_adapt_iter) override
	{
		for (std::string argID : args_to_adapt)
		{
			distrib.set_arg(argID, distrib.get_arg(argID) * _adapter_(step, num_adapt_iter));
		}
	}
};

///////////////////////////////////////////////////////////////////////
// ///////////////////// Lognorm_Proposal : class ////////////////// //
///////////////////////////////////////////////////////////////////////

class Lognorm_Proposal : public Metropolis_Hastings
{
private:
	Lognorm_Distr distrib;
	std::vector<std::string> args_to_adapt = {"sd"};

public:
	// Constructors
	Lognorm_Proposal(kwargs new_args = {{"mean", (double)0.0}, {"sd", (double)1.0}}) : distrib(new_args){};
	~Lognorm_Proposal() {}

	// --- functions (override) --- //
	double propose(std::mt19937_64 &rng, double old_value) override { return old_value * distrib.rand(rng); }
	double log_ratio() override { return distrib.get_arg("sd"); }
	void adapt(int step, int num_adapt_iter) override
	{
		for (std::string argID : args_to_adapt)
		{
			distrib.set_arg(argID, distrib.get_arg(argID) * _adapter_(step, num_adapt_iter));
		}
	}
};

///////////////////////////////////////////////////////////////////////
// ////////////////////// Prior_Proposal : class /////////////////// //
///////////////////////////////////////////////////////////////////////

// FIXME : missing...

/////////////////////////////////////////////////////////////////
// /////////////////// Base_Variable : class ///////////////// //
/////////////////////////////////////////////////////////////////

template <typename double_or_ptr_1, typename double_or_ptr_2>
class Base_Variable
{
private:
	// objects for inference
	Distribution_1D _distribution_; // FIXME: this 
	Metropolis_Hastings _proposer_;

protected:
	// general attributes
	double_or_ptr_1 value; // double or pointer to the external value
	double value_0 = 0.0;  // initial value
	bool infer = false;	   // inference mode

public:
	// Constructor // UPDATE: templates are messy in header + constructor to be cleaned
	Base_Variable(
		double_or_ptr_1 value_double_or_ptr,
		std::string distr_name,
		std::unordered_map<std::string, double_or_ptr_2> distr_args,
		std::string strat_name,
		kwargs strat_args,
		kwargs value_range) : _distribution_(), _proposer_()
	{
		value = value_double_or_ptr;

		// --- metropolis strategy --- //

		if (strat_name == UNDEF)
		{
			infer = false;
		}
		else if (strat_name == "normal")
		{
			infer = true;
			_proposer_ = Normal_Proposal(strat_args);
		}
		else if (strat_name == "lognorm")
		{
			infer = true;
			_proposer_ = Lognorm_Proposal(strat_args);
		}
		else
		{
			std::cerr << "ERROR: Unknown metropolis strategy " << strat_name << std::endl;
			std::cerr << "Exit...\n"
					  << std::endl;
			exit(EXIT_FAILURE);
		}
		_proposer_.set_range(value_range);

		// --- prior function --- //
		if ((infer == true) && (distr_name != UNDEF))
		{
			if (distr_name == "unif")
			{
				_distribution_ = Uniform_Distr(distr_args);
			}
			else if (distr_name == "normal")
			{
				_distribution_ = Normal_Distr(distr_args);
			}
			else if (distr_name == "lognorm")
			{
				_distribution_ = Lognorm_Distr(distr_args);
			}
			else if (distr_name == "gamma")
			{
				//FLAG
				std::mt19937_64 rng;
				rng.seed(235543);
				std::cerr << "before Gamma_Distr : " << std::endl;
				std::cerr << _distribution_.rand(rng) << std::endl;
				//FLAG
				_distribution_ = Gamma_Distr(distr_args);
				//FLAG
				rng.seed(235543);
				std::cerr << "after Gamma_Distr : " << std::endl;
				std::cerr << _distribution_.rand(rng) << std::endl;
				//FLAG
			}
			else if (distr_name == "exp")
			{
				_distribution_ = Exponential_Distr(distr_args);
			}
			else
			{
				std::cerr << "\nERROR: Unknown distribution " << distr_name << std::endl;
				std::cerr << "Exit...\n"
						  << std::endl;
				exit(EXIT_FAILURE);
			}
		}
		else
		{
			distr_name = UNDEF;
			distr_args = {};
		}
	}
	~Base_Variable() {}

	// Getter
	bool get_infer() const { return infer; }

	// Computations
	bool valid(double new_value) { return _proposer_.valid(new_value); }

	// Methods for _distribution_
	double distr_rand(std::mt19937_64 &rng) { return _distribution_.rand(rng); }
	double distr_log_pdf(double new_value) { return _distribution_.log_pdf(new_value); }
	double distr_log_ratio(double old_value, double new_value) { return _distribution_.log_ratio(old_value, new_value); }

	// Methods for _proposer_
	double propose(std::mt19937_64 &rng, double old_value) { return (_proposer_.propose)(rng, old_value); }
	double propose_log_ratio() { return (_proposer_.log_ratio)(); }

	double get_accept_rate() const { return _proposer_.get_accept_rate(); }
	void set_num_propos_pp() { _proposer_.set_num_propos_pp(); }
	void set_num_accept_pp() { _proposer_.set_num_accept_pp(); }
	void update_accept_rate(bool is_accepted, int step) { _proposer_.update_accept_rate(is_accepted, step); }
	void reset_moves() { _proposer_.reset_moves(); }
	void adapt(int step, int num_adapt_steps) { (_proposer_.adapt)(step, num_adapt_steps); }
};

/////////////////////////////////////////////////////////////////
// ///////////////////// Parameter : class /////////////////// //
/////////////////////////////////////////////////////////////////

class Parameter : public Base_Variable<double *, double>
{
public:
	Parameter(
		double *new_value_ptr,
		std::string prior_name = UNDEF,
		kwargs prior_args = {},
		std::string strat_name = UNDEF,
		kwargs strat_args = {},
		kwargs value_range = {}) : Base_Variable(new_value_ptr,
												 prior_name,
												 prior_args,
												 strat_name,
												 strat_args,
												 value_range) { value_0 = get_value(); }
	// Getter
	double get_value() const { return *(value); }
	double *get_address() const { return value; }
	double get_value_0() const { return value_0; }
	
	// Setter
	void set_value(double new_value) { *(value) = new_value; }

	// Aliases
	double prior_rand(std::mt19937_64 &rng) { return distr_rand(rng); }
	double prior_log_ratio(double old_value, double new_value) { return distr_log_ratio(old_value, new_value); }
};

/////////////////////////////////////////////////////////////////
// ///////////////////// Augmented : class /////////////////// //
/////////////////////////////////////////////////////////////////

class Augmented : public Base_Variable<double, double *>
{
public:
	Augmented(
		double new_value,
		std::string model_name = UNDEF,
		kwargs_ext model_args_ext = {},
		std::string strat_name = UNDEF,
		kwargs strat_args = {},
		kwargs value_range = {}) : Base_Variable(new_value,
												 model_name,
												 model_args_ext,
												 strat_name,
												 strat_args,
												 value_range) { value_0 = get_value(); }
	// Getter
	double get_value() const { return value; }
	double get_value_0() const { return value_0; }

	// Setter
	void set_value(double new_value) { value = new_value; }

	// Aliases
	double model_rand(std::mt19937_64 &rng) { return distr_rand(rng); }
	double model_log_pdf(double new_value) { return distr_log_pdf(new_value); }
	double model_log_ratio(double old_value, double new_value) { return distr_log_ratio(old_value, new_value); }
};

///////////////////////////////////////////////////////////
/////// Parameter_Dictionary : class declaration //////////
///////////////////////////////////////////////////////////

// This is a wrapper class containing aliases for the global params.
class Parameter_Dictionary
{
private:
	// Attributes
	std::vector<double> param_values = {};
	std::unordered_map<std::string, Parameter> params = {};
	size_t _counter_ = 0;

public:
	// Attributes
	std::vector<std::string> to_be_inferred = {};
	std::vector<std::string> not_to_be_inferred = {};

	void resize(size_t num) { param_values.resize(num); }
	void add_param(
		const std::string k,
		double new_value,
		std::string prior_name = UNDEF,
		kwargs prior_args = {},
		std::string strat_name = UNDEF,
		kwargs strat_args = {},
		kwargs value_range = {});

	// Getter
	double *get_address(const std::string k) const { return params.at(k).get_address(); }
	double get_value(const std::string k) const { return params.at(k).get_value(); }
	bool get_infer(const std::string k) const { return params.at(k).get_infer(); }

	// Setter
	void set_value(const std::string k, double new_value) { params.at(k).set_value(new_value); }
	void init_value(const std::string k, std::mt19937_64 &rng) { set_value(k, params.at(k).prior_rand(rng)); }

	double get_accept_rate(const std::string k) const { return params.at(k).get_accept_rate(); }
	void set_num_propos_pp(const std::string k) { params.at(k).set_num_propos_pp(); }
	void set_num_accept_pp(const std::string k) { params.at(k).set_num_accept_pp(); }
	void update_accept_rate(const std::string k, bool is_accepted, int step) { params.at(k).update_accept_rate(is_accepted, step); }
	void reset_moves(const std::string k) { params.at(k).reset_moves(); }

	void adapt(const std::string k, int step, int m_num_adapt_iter) { params.at(k).adapt(step, m_num_adapt_iter); }

	// Computations
	bool valid_param(const std::string k, double new_value) { return params.at(k).valid(new_value); }
	double propose(const std::string k, std::mt19937_64 &rng, double old_value) { return params.at(k).propose(rng, old_value); }
	double log_ratio_proposal(const std::string k) { return params.at(k).propose_log_ratio(); }
	double log_ratio_prior(const std::string k, double old_value, double new_value) { return params.at(k).prior_log_ratio(old_value, new_value); }

	// Display
	void print();
};

/////////////////////////////////////////////////////////
//////////// global parameters declaration //////////////
/////////////////////////////////////////////////////////

extern Parameter_Dictionary glob_params;

#endif // PARAMETER__H