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

#include "tools.h"
#include "function.h"

//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//

#define D_DISCR_TYPE_NONE 0 // none
#define D_DISCR_TYPE_LEFT 1 // left
#define D_DISCR_TYPE_CNTR 2 // center
#define D_DISCR_TYPE_DFND 3 // defined

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
	const std::vector<dbl_vec> showcase(const dbl_vec &opts = {},
										int n_bins = 2000,
										double ytol = 1.e-5) const;

protected:
	//===== SETTER =====//
	void init_functional_discretization(ushort discr_type);
	void init_functional_truncation(ushort discr_type);

	void init_min_zero_reshaping(bool discrete, bool is_open_zero);
	void init_min_zero_discr_defined_reshaping();

	//===== GETTER =====//
	// note: distribution specific virtual methods
	virtual double raw_pdf(double x, const dbl_vec &opts) const;
	virtual double raw_cdf(double x, const dbl_vec &opts) const;
	virtual double raw_logpdf(double x, const dbl_vec &opts) const;
	virtual double raw_logratio(double xold, double xnew, const dbl_vec &opts) const;
	virtual void update_args();
	virtual double raw_rand(rnd_gen &prng);

private:
	double _rejection_sampling_(rnd_gen &prng, double (Distribution::*func)(rnd_gen &));

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

	// default case
	double _trunc_dummy_(double x, const dbl_vec &opts) const;

	// continuos case
	double _trunc_min_(double x, const dbl_vec &opts) const;
	double _trunc_max_(double x, const dbl_vec &opts) const;
	double _trunc_both_sides_(double x, const dbl_vec &opts) const;
	double _trunc_rand_(rnd_gen &prng);

	// case of discrete defined
	double _trunc_min_discr_defined_(double x, const dbl_vec &opts) const;
	double _trunc_max_discr_defined_(double x, const dbl_vec &opts) const;
	double _trunc_both_sides_discr_defined_(double x, const dbl_vec &opts) const;
	// note : _trunx_rand_discr_defined == _trunc_rand_

	// case of center discretization
	double _trunc_min_discr_center_(double x, const dbl_vec &opts) const;
	double _trunc_max_discr_center_(double x, const dbl_vec &opts) const;
	double _trunc_both_sides_discr_center_(double x, const dbl_vec &opts) const;
	double _trunc_rand_discr_center_(rnd_gen &prng);

	// case of left discretization
	double _trunc_min_discr_left_(double x, const dbl_vec &opts) const;
	double _trunc_max_discr_left_(double x, const dbl_vec &opts) const;
	double _trunc_both_sides_discr_left_(double x, const dbl_vec &opts) const;
	double _trunc_rand_discr_left_(rnd_gen &prng);

	//----- Functional methods for the log-truncation -----//

	// default case
	double _logtrunc_dummy_(double x, const dbl_vec &opts) const;

	// continuos case
	double _logtrunc_min_(double x, const dbl_vec &opts) const;
	double _logtrunc_max_(double x, const dbl_vec &opts) const;
	double _logtrunc_both_sides_(double x, const dbl_vec &opts) const;

	// case of discrete defined
	double _logtrunc_min_discr_defined_(double x, const dbl_vec &opts) const;
	double _logtrunc_max_discr_defined_(double x, const dbl_vec &opts) const;
	double _logtrunc_both_sides_discr_defined_(double x, const dbl_vec &opts) const;

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
	double _shift_min_discr_defined_(double x, const dbl_vec &opts) const;
	double _shift_min_discr_center_(double x, const dbl_vec &opts) const;
	double _shift_min_discr_left_(double x, const dbl_vec &opts) const;

	//===== DEVELOPER CHECK =====//
	double __developer_check_pdf_norm__(const dbl_vec &opts = {},
										int n_bins = 5000,
										double ytol = 1.e-8) const;
};

/////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////  distribution switchboard  ///////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

/* Distribution with default scalar arguments. */
std::shared_ptr<Distribution> distribution_switchboard(
	const std::string &name = NULL_STRING);

#endif // HOUSEPI__DISTRIBUTION__H
