/*
 * distribution.h
 *
 * Copyright (C) 2025 Francesco Camaglia, Institut Pasteur
 */

#ifndef HOUSEPI__DISTRIBUTION__H
#define HOUSEPI__DISTRIBUTION__H

#include <memory>
#include <string>
#include <vector>

#include "function.h"

#define D_DISCR_TYPE_NONE 0 // none
#define D_DISCR_TYPE_LEFT 1 // left
#define D_DISCR_TYPE_CNTR 2 // center
#define D_DISCR_TYPE_DFND 3 // defined

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
class Distribution : public Function
{
public:
	using Function::Function;

	//===== SETTER =====//
	virtual void init_reshaping(bool discrete);

	//===== GETTER =====//
	bool is_discrete() const;
	double pdf(double x, const dbl_vec &opts = {}) const;
	double cdf(double x, const dbl_vec &opts = {}) const;
	double logpdf(double x, const dbl_vec &opts = {}) const;
	double logratio(double xold, double xnew, const dbl_vec &opts = {}) const;
	double rand(rnd_gen &prng);

	//===== SHOW =====//
    const std::map<double, dbl_vec> showcase(dbl_vec opts = {},
											 dbl_vec xlim = {-20., 20.}) const;

protected:
	//===== SETTER =====//
	void init_discretization(ushort discr_type);
	void init_truncation(ushort discr_type);

	void _init_min_zero_reshaping_(bool discrete, bool is_open_zero);

	//===== GETTER =====//
	// note: distribution specific virtual methods
	virtual double _pdf_(double x, const dbl_vec &opts) const;
	virtual double _cdf_(double x, const dbl_vec &opts) const;
	virtual double _logpdf_(double x, const dbl_vec &opts) const;
	virtual double _logratio_(double xold, double xnew, const dbl_vec &opts) const;
	virtual void _update_args_();
	virtual double _rand_(rnd_gen &prng);

private:
	//===== Functional methods for public getter methods =====//
	double (Distribution::*_func_pdf_)(double, const dbl_vec &) const = &Distribution::_default_pdf_;
	double (Distribution::*_func_cdf_)(double, const dbl_vec &) const = &Distribution::_default_cdf_;
	double (Distribution::*_func_logpdf_)(double, const dbl_vec &) const = &Distribution::_default_logpdf_;
	double (Distribution::*_func_logratio_)(double, double, const dbl_vec &) const = &Distribution::_default_logratio_;
	double (Distribution::*_func_rand_)(rnd_gen &) = &Distribution::_default_rand_;

	//----- Default functional methods -----//
	double _default_pdf_(double x, const dbl_vec &opts) const;
	double _default_cdf_(double x, const dbl_vec &opts) const;
	double _default_logpdf_(double x, const dbl_vec &opts) const;
	double _default_logratio_(double xold, double xnew, const dbl_vec &opts) const;
	double _default_rand_(rnd_gen &prng);

	//----- Functional methods for the discrete defined distributions -----//
	double _discr_defined_pdf_(double x, const dbl_vec &opts) const;
	double _discr_defined_cdf_(double x, const dbl_vec &opts) const;
	double _discr_defined_logpdf_(double x, const dbl_vec &opts) const;
	double _discr_defined_logratio_(double xold, double xnew, const dbl_vec &opts) const;
	// note : _discr_defined_rand_ == _default_rand_

	//----- Functional methods for the center adjusted discretization -----//
	double _discr_center_pdf_(double x, const dbl_vec &opts) const;
	double _discr_center_cdf_(double x, const dbl_vec &opts) const;
	double _discr_center_logpdf_(double x, const dbl_vec &opts) const;
	double _discr_center_logratio_(double xold, double xnew, const dbl_vec &opts) const;
	double _discr_center_rand_(rnd_gen &prng);

	//----- Functional methods for the left adjusted discretization -----//
	double _discr_left_pdf_(double x, const dbl_vec &opts) const;
	double _discr_left_cdf_(double x, const dbl_vec &opts) const;
	double _discr_left_logpdf_(double x, const dbl_vec &opts) const;
	double _discr_left_logratio_(double xold, double xnew, const dbl_vec &opts) const;
	double _discr_left_rand_(rnd_gen &prng);

	//===== Functional methods for truncation =====//
	double (Distribution::*_func_trunc_)(double, const dbl_vec &opts) const = &Distribution::_trunc_dummy_;
	double (Distribution::*_func_logtrunc_)(double, const dbl_vec &opts) const = &Distribution::_logtrunc_dummy_;
	double (Distribution::*_func_shift_)(double, const dbl_vec &opts) const = &Distribution::_shift_dummy_;

	//----- Functional methods for the truncation and/or discretization -----//
	double _trunc_dummy_(double x, const dbl_vec &opts) const;
	double _trunc_min_(double x, const dbl_vec &opts) const;
	double _trunc_max_(double x, const dbl_vec &opts) const;
	double _trunc_both_sides_(double x, const dbl_vec &opts) const;

	// case of center discretization
	double _trunc_min_discr_center_(double x, const dbl_vec &opts) const;
	double _trunc_max_discr_center_(double x, const dbl_vec &opts) const;
	double _trunc_both_sides_discr_center_(double x, const dbl_vec &opts) const;

	// case of left discretization
	double _trunc_min_discr_left_(double x, const dbl_vec &opts) const;
	double _trunc_max_discr_left_(double x, const dbl_vec &opts) const;
	double _trunc_both_sides_discr_left_(double x, const dbl_vec &opts) const;

	//----- Functional methods for the (log) truncation -----//
	double _logtrunc_dummy_(double x, const dbl_vec &opts) const;
	double _logtrunc_min_(double x, const dbl_vec &opts) const;
	double _logtrunc_max_(double x, const dbl_vec &opts) const;
	double _logtrunc_both_sides_(double x, const dbl_vec &opts) const;

	// case of center discretization
	double _logtrunc_min_discr_center_(double x, const dbl_vec &opts) const;
	double _logtrunc_max_discr_center_(double x, const dbl_vec &opts) const;
	double _logtrunc_both_sides_discr_center_(double x, const dbl_vec &opts) const;

	// case of left discretization
	double _logtrunc_min_discr_left_(double x, const dbl_vec &opts) const;
	double _logtrunc_max_discr_left_(double x, const dbl_vec &opts) const;
	double _logtrunc_both_sides_discr_left_(double x, const dbl_vec &opts) const;

	//----- Functional methods to shift cumulative -----//
	double _shift_dummy_(double x, const dbl_vec &opts) const;
	double _shift_min_(double x, const dbl_vec &opts) const;
	double _shift_min_discr_center_(double x, const dbl_vec &opts) const;
	double _shift_min_discr_left_(double x, const dbl_vec &opts) const;
};

/////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////  distribution switchboard  ///////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

/* Distribution with default scalar arguments. */
std::shared_ptr<Distribution> distribution_switchboard(
	const std::string &name = NULL_STRING);

#endif // HOUSEPI__DISTRIBUTION__H
