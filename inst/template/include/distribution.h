/* housepi - Epidemiological Analysis of Household Outbreaks.
 *
 * Copyright © 2026 Institut Pasteur (Paris), CNRS and INSERM.
 * Distributed under the terms of the GNU General Public License (GPLv3).
 * See the LICENSE file for details.
 */

#pragma once

#include <memory>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

#include "tools/nan.h"
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
protected:
	enum class Trunc : uint8_t
	{
		None,  /// No truncation
		Active /// Active truncation
	} trunc_mode_ = Trunc::None;

public:
	///////////////////////////
	//////  CONSTRUCTOR  //////
	///////////////////////////

	Distribution(const std::vector<std::pair<std::string, double>> &args = {},
				 double domain_min = hsp::na::DOUBLE, bool is_min_open = true,
				 double domain_max = hsp::na::DOUBLE, bool is_max_open = true)
		: Function(args,
				   domain_min, is_min_open,
				   domain_max, is_max_open) { set_eval_as_distr_call(); }
	virtual ~Distribution() = default;

	//////////////////////
	//////  GETTER  //////
	//////////////////////

	double rand(hsp::rnd_gen &prng) const
	{
		switch (trunc_mode_)
		{
		case Trunc::None:
			return relocate_(cast_discr_(raw_rand_(prng)));
		case Trunc::Active:
			return relocate_(cast_discr_(inverse_sampling_(prng)));
		}
		__builtin_unreachable();
	}
	/// Return value of log-probability density/mass function.
	double logpdf(double x) const { return log_func_(delocate_(x)); }
	/// Return value of log-ratio of probability density/mass functions.
	double logratio(double xo, double xn) const { return log_ratio_(delocate_(xo), delocate_(xn)); }

	std::string err_name() const override { return "distribution '" + std::string(name()) + "'"; }

protected:
	/////////////////////////////////////////////////////////////////////////////////////////////
	////////////////////////////  Distribution : virtual [protected] ////////////////////////////
	/////////////////////////////////////////////////////////////////////////////////////////////

	/// Raw random sampling method of the distribution.
	/// @param prng : pseudo-random number generator.
	virtual double raw_rand_(hsp::rnd_gen &) const { throw std::invalid_argument(
		"Called undefined method 'rand' of " + err_name() + "."); }

	/// Raw probability density/mass function of the distribution.
	virtual double raw_pdf_(double) const = 0;

	/// Raw log-probability density/mass function of the distribution.
	virtual double raw_logpdf_(double) const { throw std::invalid_argument(
		"Called to undefined method 'logpdf' of " + err_name() + "."); }

	/// Raw cumulative distribution function.
	virtual double raw_cdf_(double) const = 0;

	/// Raw quantile.
	virtual double raw_quantile_(double) const { throw std::invalid_argument(
		"Called to undefined method 'quantile' of " + err_name() + "."); }

	/// Raw method to call.
	double raw_function_(double) const override final { throw std::runtime_error(
		"Forbbiden call to raw_function_ for " + err_name() + "\n"); }

	/// Polymorphic function method.
	double function_(double x) const override final
	{
		double y;
		switch (discr_mode_)
		{
		case Mode::Discr::None:
			y = raw_pdf_(x);
			break;
		case Mode::Discr::Center:
		{
			const double i = std::round(x);
			y = raw_cdf_(i + .5) - raw_cdf_(i - .5);
			break;
		}
		case Mode::Discr::Left:
		{
			const double i = std::floor(x);
			y = raw_cdf_(i + 1.) - raw_cdf_(i);
			break;
		}
		case Mode::Discr::Defined:
			y = raw_pdf_(std::floor(x));
			break;
		default:
			__builtin_unreachable();
		}

		return y / (this->*trunc_)();
	}

	/// Polymorphic log-function method.
	double log_func_(double x) const override final
	{
		double y;
		switch (discr_mode_)
		{
		case Mode::Discr::None:
			y = raw_logpdf_(x);
			break;
		case Mode::Discr::Center:
		{
			const double i = std::round(x);
			y = raw_cdf_(i + 0.5) - raw_cdf_(i - 0.5);
			break;
		}
		case Mode::Discr::Left:
		{
			const double i = std::floor(x);
			y = raw_cdf_(i + 1.) - raw_cdf_(i);
			break;
		}
		case Mode::Discr::Defined:
			y = raw_logpdf_(std::floor(x));
			break;
		default:
			__builtin_unreachable();
		}

		return y - logtrunc_();
	}

	/// Polymorphic primitive method.
	double primitive_(double x) const override final
	{
		double adjusted_x;
		switch (discr_mode_)
		{
		case Mode::Discr::None:
			adjusted_x = x;
			break;
		case Mode::Discr::Center:
			adjusted_x = std::round(x) + 0.5;
			break;
		case Mode::Discr::Left:
			adjusted_x = std::floor(x) + 1.0;
			break;
		case Mode::Discr::Defined:
			adjusted_x = std::floor(x);
			break;
		default:
			__builtin_unreachable();
		}

		return (raw_cdf_(adjusted_x) - (this->*shift_)()) / (this->*trunc_)();
	}

	/// Polymorphic log_ratio method.
	double log_ratio_(double xo, double xn) const override final
	{
		switch (discr_mode_)
		{
		case Mode::Discr::None:
			return raw_logratio_(xo, xn);
		case Mode::Discr::Center:
		{
			const double io = std::round(xo), in = std::round(xn);
			if (io == in)
				return 0.;
			return std::log(raw_cdf_(in + .5) - raw_cdf_(in - .5)) -
				   std::log(raw_cdf_(io + .5) - raw_cdf_(io - .5));
		}
		case Mode::Discr::Left:
		{
			const double io = std::floor(xo), in = std::floor(xn);
			if (io == in)
				return 0.;
			return std::log(raw_cdf_(in + 1.) - raw_cdf_(in)) -
				   std::log(raw_cdf_(io + 1.) - raw_cdf_(io));
		}
		case Mode::Discr::Defined:
			return raw_logratio_(std::floor(xo), std::floor(xn));
		}
		__builtin_unreachable();
	}

	/////////////////////////////////////////////////////////////////////////////////////////////
	////////////////////////////  Distribution : getter [protected] /////////////////////////////
	/////////////////////////////////////////////////////////////////////////////////////////////

	double cast_discr_(double x) const
	{
		switch (discr_mode_)
		{
		case Mode::Discr::Center:
			return static_cast<double>(std::lround(x));
		case Mode::Discr::Left:
			return std::floor(x);
		case Mode::Discr::Defined:
		case Mode::Discr::None:
			return x;
		}
		__builtin_unreachable();
	}

	double inverse_sampling_(hsp::rnd_gen &prng) const
	{
		static std::uniform_real_distribution<double> unif_gen(0., 1.);
		const double r = unif_gen(prng);
		const double sx = lower_arg_.is_empty() ? 0. : (1. - r) * raw_cdf_(lower());
		const double dx = upper_arg_.is_empty() ? r : r * raw_cdf_(upper());
		return raw_quantile_(sx + dx);
	}

	/////////////////////////////////////////////////////////////////////////////////////////////
	//////////////////////  Distribution : functional truncation [private] //////////////////////
	/////////////////////////////////////////////////////////////////////////////////////////////

	void init_trunc_helper_();
	/// UPDATE: if lower/upper values don't change, I could just cached them.

	double (Distribution::*trunc_)() const = &Distribution::do_trunc_dummy_;

	// default case
	double do_trunc_dummy_() const { return 1.; }

	// continuos case
	double do_trunc_min_() const { return 1. - raw_cdf_(lower()); }
	double do_trunc_max_() const { return raw_cdf_(upper()); }
	double do_trunc_both_sides_() const { return raw_cdf_(upper()) - raw_cdf_(lower()); }

	// case of discrete defined
	double do_trunc_min_discr_defined_() const { return 1. - raw_cdf_(lower() - 1); }
	double do_trunc_max_discr_defined_() const { return raw_cdf_(upper() + 1); }
	double do_trunc_both_sides_discr_defined_() const { return raw_cdf_(upper() + 1) - raw_cdf_(lower() - 1); }

	// case of center discretization
	double do_trunc_min_discr_center_() const { return 1. - raw_cdf_(std::ceil(lower()) - 0.5); }
	double do_trunc_max_discr_center_() const { return raw_cdf_(std::floor(upper()) + 0.5); }
	double do_trunc_both_sides_discr_center_() const { return raw_cdf_(std::floor(upper()) + 0.5) - raw_cdf_(std::ceil(lower()) - 0.5); }

	// case of left discretization
	double do_trunc_min_discr_left_() const { return 1. - raw_cdf_(std::floor(lower())); }
	double do_trunc_max_discr_left_() const { return raw_cdf_(std::ceil(upper() + 1)); }
	double do_trunc_both_sides_discr_left_() const { return raw_cdf_(std::ceil(upper() + 1)) - raw_cdf_(std::floor(lower())); }

	double logtrunc_() const {
		switch (trunc_mode_)
		{
		case Trunc::None:
			return 0;
		case Trunc::Active:
			return std::log((this->*trunc_)());
		}
		__builtin_unreachable();
	} 

	/////////////////////////////////////////////////////////////////////////////////////////////
	/////////////////////////  Distribution : functional shift [private] ////////////////////////
	/////////////////////////////////////////////////////////////////////////////////////////////

	double (Distribution::*shift_)() const = &Distribution::do_shift_dummy_;
	double do_shift_dummy_() const { return 0.; }
	double do_shift_min_() const { return raw_cdf_(lower()); }
	double do_shift_min_discr_defined_() const { return raw_cdf_(lower() - 1); }
	double do_shift_min_discr_center_() const { return raw_cdf_(std::floor(lower()) - 0.5); }
	double do_shift_min_discr_left_() const { return raw_cdf_(std::floor(lower())); }
};

/// Compose registries to get a distribution from them.
std::unique_ptr<Distribution> compose_distr_registries(std::string_view opt,
													   const std::string &err_type,
													   const std::vector<std::string> &registries,
													   bool is_cumul = false);
