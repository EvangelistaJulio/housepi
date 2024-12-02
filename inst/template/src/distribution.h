/*
 * distribution.h
 *
 * Copyright (C) 2024 Francesco Camaglia, Institut Pasteur
 */

#ifndef HOUSEPI__DISTRIBUTION__H
#define HOUSEPI__DISTRIBUTION__H

#include <memory>
#include <random>
#include <string>
#include <vector>
#include <map>
#include <unordered_map>
#include <cmath>

#include "tools.h"
#include "function.h"

//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//

/////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////  Base Distribution : class  ////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

/* A parent class to all distributions. */
class Base_Distribution : public Base_Function
{
public:
	using Base_Function::Base_Function;

	double pdf(double x) const;
	double cdf(double x) const;
	double logpdf(double x) const;
	double logratio(double xold, double xnew) const;
	double rand(sys_rnd_gen &prng);

protected:
	//----- Distribution specific virtual methods -----//
	virtual double _pdf_(double x) const;
	virtual double _cdf_(double x) const;
	virtual double _logpdf_(double x) const;
	virtual double _logratio_(double xold, double xnew) const;
	virtual double _rand_(sys_rnd_gen &prng);

	//----- Adjust methods for truncation and/or discretizat -----//
	void _assign_reshaping_(const dictargs &range, const std::string &discr);
	void _check_rand_attempts_(int &n_attempts) const;

	//----- Functional methods for public methods -----//
	double (Base_Distribution::*_func_pdf_)(double) const = &Base_Distribution::_default_pdf_;
	double (Base_Distribution::*_func_cdf_)(double) const = &Base_Distribution::_default_cdf_;
	double (Base_Distribution::*_func_logpdf_)(double) const = &Base_Distribution::_default_logpdf_;
	double (Base_Distribution::*_func_logratio_)(double, double) const = &Base_Distribution::_default_logratio_;
	double (Base_Distribution::*_func_rand_)(sys_rnd_gen &) = &Base_Distribution::_default_rand_;

	//----- Functional methods for truncation -----//
	double (Base_Distribution::*_func_trunc_)(double) const = &Base_Distribution::_trunc_dummy_;
	double (Base_Distribution::*_func_logtrunc_)(double) const = &Base_Distribution::_logtrunc_dummy_;
	double (Base_Distribution::*_func_shift_)(double) const = &Base_Distribution::_shift_dummy_;

private:
	double _undefined_call_error_(const std::string &txt = UNDEFINED, double x = 0, double y = 0) const;

	//----- Deafult functional methods -----//
	double _default_pdf_(double x) const;
	double _default_cdf_(double x) const;
	double _default_logpdf_(double x) const;
	double _default_logratio_(double xold, double xnew) const;
	double _default_rand_(sys_rnd_gen &prng);

	//----- Functional methods for the center adjusted discretization -----//
	double _discr_center_pdf_(double x) const;
	double _discr_center_cdf_(double x) const;
	double _discr_center_logpdf_(double x) const;
	double _discr_center_logratio_(double xold, double xnew) const;
	double _discr_center_rand_(sys_rnd_gen &prng);

	//----- Functional methods for the left adjusted discretization -----//
	double _discr_left_pdf_(double x) const;
	double _discr_left_cdf_(double x) const;
	double _discr_left_logpdf_(double x) const;
	double _discr_left_logratio_(double xold, double xnew) const;
	double _discr_left_rand_(sys_rnd_gen &prng);

	//----- Functional methods for the truncation and/or discretization -----//
	double _trunc_dummy_(double x) const;
	double _trunc_min_(double x) const;
	double _trunc_max_(double x) const;
	double _trunc_both_sides_(double x) const;

	//----- Functional methods for the (log) truncation -----//
	double _logtrunc_dummy_(double x) const;
	double _logtrunc_min_(double x) const;
	double _logtrunc_max_(double x) const;
	double _logtrunc_both_sides_(double x) const;

	//----- Functional methods to shift cumulative -----//
	double _shift_dummy_(double x) const;
	double _shift_min_(double x) const;
};

/////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////  distribution switchboard  ///////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

/* overloaded distribution with static values. */
std::shared_ptr<Base_Distribution> distribution_switchboard(const std::string &name = UNDEFINED,
															const dictargs &args = {},
															const dictargs &range = {});

/* overloaded distribution with pointer to values. */
std::shared_ptr<Base_Distribution> distribution_switchboard(const std::string &name = UNDEFINED,
															const dictargs_ptr &args = {},
															const dictargs &range = {});

#endif // HOUSEPI__DISTRIBUTION__H
