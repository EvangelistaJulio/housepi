/*
 *        __
 *       / /_  ____  __  __ ____ __  ____  __
 *      / __ \/ __ \/ / / / ___/ _ \/ __ \/ /
 *     / / / / /_/ / /_/ (__  )  __/ /_/ / /
 *    /_/ /_/\____/\____/____/\___/ ____/_/
 *                               /_/
 *
 * housepi - Epidemiological Analysis of Household Outbreaks.
 *
 * Author: Francesco Camaglia
 * Copyright © 2025 Institut Pasteur (Paris), CNRS and INSERM.
 *
 * Distributed under the terms of the GNU General Public License (GPLv3).
 * See the LICENSE file for details.
 *
 * distribution.h
 *
 */

#ifndef HOUSEPI__DISTRIBUTION__H
#define HOUSEPI__DISTRIBUTION__H

#include "tools/types.h"
#include "tools/na_utils.h"
#include "tools/random.h"
#include "function.h"

/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////// *** Distribution : class *** ///////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

/// Parent class of all distributions.
class Distribution : public Function
{
public:
	/////////////////////////////////////////////////////////////////////////////////////////////
	///////////////////////////////  Distribution : constructor  ////////////////////////////////
	/////////////////////////////////////////////////////////////////////////////////////////////

	Distribution(const std::string &name = hsp::na::STRING,
				 const str_vec &arg_name = {}, const dbl_vec &arg_value = {},
				 double domain_min = hsp::na::DOUBLE, bool is_min_open = true,
				 double domain_max = hsp::na::DOUBLE, bool is_max_open = true)
		: Function(name, arg_name, arg_value,
				   domain_min, is_min_open,
				   domain_max, is_max_open)
	{
		err_type_ = "distribution";
		set_eval_as_distr_call();
	}
	virtual ~Distribution() = default;

	/////////////////////////////////////////////////////////////////////////////////////////////
	//////////////////////////////////  Distribution : getter  //////////////////////////////////
	/////////////////////////////////////////////////////////////////////////////////////////////

	double rand(rnd_gen &prng) { return (this->*rand_)(prng); }
	/// Return value of probability density/mass function.
	double pdf(double x) const { return eval_function_as_distr_(x); }
	/// Return value of log-probability density/mass function.
	double logpdf(double x) const { return (this->*log_func_)(x); }
	/// Return value of cumulative distribution.
	double cdf(double x) const { return eval_primitive_as_distr_(x); }
	/// Return value of log-ratio of probability density/mass functions.
	double logratio(double xo, double xn) const { return (this->*log_ratio_)(xo, xn); }

protected:
	/////////////////////////////////////////////////////////////////////////////////////////////
	////////////////////////////  Distribution : setter [protected] /////////////////////////////
	/////////////////////////////////////////////////////////////////////////////////////////////

	/// Initialize functional discretization.
	void init_discr_helper_(unsigned short discr_type) override final;
	void init_trunc_helper_();

	/////////////////////////////////////////////////////////////////////////////////////////////
	////////////////////////////  Distribution : virtual [protected] ////////////////////////////
	/////////////////////////////////////////////////////////////////////////////////////////////

	/// Raw random sampling method of the distribution.
	/// @param prng : pseudo-random number generator.
	virtual double raw_rand_(rnd_gen &) { throw std::invalid_argument(
		"Called undefined method 'rand' of " + err_name() + "."); }

	/// Raw probability density/mass function of the distribution.
	virtual double raw_pdf_(double) const { throw std::invalid_argument(
		"Called to undefined method 'pdf' of " + err_name() + "."); }

	/// Raw log-probability density/mass function of the distribution.
	virtual double raw_logpdf_(double) const { throw std::invalid_argument(
		"Called to undefined method 'logpdf' of " + err_name() + "."); }

	/// Raw cumulative distribution function.
	virtual double raw_cdf_(double) const { throw std::invalid_argument(
		"Called to undefined method 'cdf' of " + err_name() + "."); }

	/// Raw quantile.
	virtual double raw_quantile_(double) const { throw std::invalid_argument(
		"Called to undefined method 'quantile' of " + err_name() + "."); }

	/////////////////////////////////////////////////////////////////////////////////////////////
	////////////////////////////  Distribution : getter [protected] /////////////////////////////
	/////////////////////////////////////////////////////////////////////////////////////////////

	static double cast_center_(double x) { return static_cast<double>(std::lround(x)); }
	static double cast_left_(double x) { return std::floor(x); }

	double (Distribution::*rand_)(rnd_gen &) = &Distribution::raw_rand_;
	double do_discr_center_rand_(rnd_gen &prng) { return cast_center_(raw_rand_(prng)); }
	double do_discr_left_rand_(rnd_gen &prng) { return cast_left_(raw_rand_(prng)); }

	double raw_function_(double x) const override final { return (this->*trunc_)(raw_pdf_(x)); }
	double do_function_discr_center_(double x) const override final
	{
		const double i = std::round(x);
		return (this->*trunc_)(raw_cdf_(i + .5) - raw_cdf_(i - .5));
	}
	double do_function_discr_left_(double x) const override final
	{
		const double i = std::floor(x);
		return (this->*trunc_)(raw_cdf_(i + 1.) - raw_cdf_(i));
	}

	double raw_logfunc_(double x) const override final { return (this->*logtrunc_)(raw_logpdf_(x)); }

	double raw_primitive_(double x) const override final { return (this->*trunc_)((this->*shift_)(raw_cdf_(x))); }
	double do_primitive_discr_center_(double x) const override final { return (this->*trunc_)((this->*shift_)(raw_cdf_(std::round(x) + .5))); }
	double do_primitive_discr_left_(double x) const override final { return (this->*trunc_)((this->*shift_)(raw_cdf_(std::floor(x) + 1.))); }

	double do_logratio_discr_center_(double xo, double xn) const override final
	{
		const double io = std::round(xo), in = std::round(xn);
		return std::log((raw_cdf_(in + .5) - raw_cdf_(in - .5)) / (raw_cdf_(io + .5) - raw_cdf_(io - .5)));
	}
	double do_logratio_discr_left_(double xo, double xn) const override final
	{
		const double io = std::floor(xo), in = std::floor(xn);
		return std::log((raw_cdf_(in + 1.) - raw_cdf_(in)) / (raw_cdf_(io + 1.) - raw_cdf_(io)));
	}

	/////////////////////////////////////////////////////////////////////////////////////////////
	//////////////////////  Distribution : functional truncation [private] //////////////////////
	/////////////////////////////////////////////////////////////////////////////////////////////

	// UPDATE: this static isn't elegant
	double inverse_sampling_(rnd_gen &prng)
	{
		static std::uniform_real_distribution<double> unif_gen(0., 1.); 
		const double r = unif_gen(prng);
		const double sx = range_lower_.is_empty() ? 0. : (1. - r) * raw_cdf_(lower());
		const double dx = range_upper_.is_empty() ? r : r * raw_cdf_(upper());
		return raw_quantile_(sx + dx);
	}

	double (Distribution::*trunc_)(double) const = &Distribution::do_trunc_dummy_;

	// default case
	double do_trunc_dummy_(double x) const { return x; }

	// continuos case
	double do_trunc_min_(double x) const { return x / (1. - raw_cdf_(lower())); }
	double do_trunc_max_(double x) const { return x / (raw_cdf_(upper())); }
	double do_trunc_both_sides_(double x) const { return x / (raw_cdf_(upper()) - raw_cdf_(lower())); }
	double do_trunc_rand_(rnd_gen &prng) { return inverse_sampling_(prng); }

	// case of discrete defined
	double do_trunc_min_discr_defined_(double x) const { return x / (1. - raw_cdf_(lower() - 1)); }
	double do_trunc_max_discr_defined_(double x) const { return x / (raw_cdf_(upper() + 1)); }
	double do_trunc_both_sides_discr_defined_(double x) const { return x / (raw_cdf_(upper() + 1) - raw_cdf_(lower() - 1)); }
	// note : do_trunx_rand_discr_defined == do_trunc_rand_

	// case of center discretization
	double do_trunc_min_discr_center_(double x) const { return x / (1. - raw_cdf_(std::ceil(lower()) - .5)); }
	double do_trunc_max_discr_center_(double x) const { return x / (raw_cdf_(std::floor(upper()) + .5)); }
	double do_trunc_both_sides_discr_center_(double x) const { return x / (raw_cdf_(std::floor(upper()) + .5) - raw_cdf_(std::ceil(lower()) - .5)); }
	double do_trunc_rand_discr_center_(rnd_gen &prng) { return cast_center_(do_trunc_rand_(prng)); }

	// case of left discretization
	double do_trunc_min_discr_left_(double x) const { return x / (1. - raw_cdf_(std::floor(lower()))); }
	double do_trunc_max_discr_left_(double x) const { return x / (raw_cdf_(std::ceil(upper() + 1))); }
	double do_trunc_both_sides_discr_left_(double x) const { return x / (raw_cdf_(std::ceil(upper() + 1)) - raw_cdf_(std::floor(lower()))); }
	double do_trunc_rand_discr_left_(rnd_gen &prng) { return cast_left_(do_trunc_rand_(prng)); }

	/////////////////////////////////////////////////////////////////////////////////////////////
	////////////////////  Distribution : functional log-truncation [private] ////////////////////
	/////////////////////////////////////////////////////////////////////////////////////////////

	/// note : this is not improving the case of discretization with repspect to simply log of pdf
	double (Distribution::*logtrunc_)(double) const = &Distribution::do_logtrunc_dummy_;

	// default case
	double do_logtrunc_dummy_(double x) const { return x; }

	// continuos case
	double do_logtrunc_min_(double x) const { return x - std::log(1. - raw_cdf_(lower())); }
	double do_logtrunc_max_(double x) const { return x - std::log(raw_cdf_(upper())); }
	double do_logtrunc_both_sides_(double x) const { return x - std::log(raw_cdf_(upper()) - raw_cdf_(lower())); }

	// case of discrete defined
	double do_logtrunc_min_discr_defined_(double x) const { return x - std::log(1. - raw_cdf_(lower() - 1)); }
	double do_logtrunc_max_discr_defined_(double x) const { return x - std::log(raw_cdf_(upper() + 1)); }
	double do_logtrunc_both_sides_discr_defined_(double x) const { return x - std::log(raw_cdf_(upper() + 1) - raw_cdf_(lower() - 1)); }

	// case of center discretization
	double do_logtrunc_min_discr_center_(double x) const { return x - std::log(1. - raw_cdf_(std::ceil(lower()) - .5)); }
	double do_logtrunc_max_discr_center_(double x) const { return x - std::log(raw_cdf_(std::floor(upper()) + .5)); }
	double do_logtrunc_both_sides_discr_center_(double x) const { return x - std::log(raw_cdf_(std::floor(upper()) + .5) - raw_cdf_(std::ceil(lower()) - .5)); }

	// case of left discretization
	double do_logtrunc_min_discr_left_(double x) const { return x - std::log(1. - raw_cdf_(std::floor(lower()))); }
	double do_logtrunc_max_discr_left_(double x) const { return x - std::log(raw_cdf_(std::ceil(upper()))); }
	double do_logtrunc_both_sides_discr_left_(double x) const { return x - std::log(raw_cdf_(std::ceil(upper())) - raw_cdf_(std::floor(lower()))); }

	/////////////////////////////////////////////////////////////////////////////////////////////
	/////////////////////////  Distribution : functional shift [private] ////////////////////////
	/////////////////////////////////////////////////////////////////////////////////////////////

	double (Distribution::*shift_)(double) const = &Distribution::do_shift_dummy_;
	double do_shift_dummy_(double x) const { return x; }
	double do_shift_min_(double x) const { return x - raw_cdf_(lower()); }
	double do_shift_min_discr_defined_(double x) const { return x - raw_cdf_(lower() - 1); }
	double do_shift_min_discr_center_(double x) const { return x - raw_cdf_(std::floor(lower()) - .5); }
	double do_shift_min_discr_left_(double x) const { return x - raw_cdf_(std::floor(lower())); }
};

#endif // HOUSEPI__DISTRIBUTION__H
