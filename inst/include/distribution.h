/* housepi - Epidemiological Analysis of Household Outbreaks.
 *
 * Copyright © 2026 Institut Pasteur (Paris), CNRS and INSERM.
 * Distributed under the terms of the GNU General Public License (GPLv3).
 * See the LICENSE file for details.
 */

#pragma once

#include <cstdint>
#include <memory>
#include <random>
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

namespace hsp
{
	/// Parent class of all distributions.
	class Distribution : public Function
	{
	public:
		///////////////////////////
		//////  CONSTRUCTOR  //////
		///////////////////////////

		Distribution(const std::vector<std::pair<std::string, double>> &args = {},
					 double domain_min = na::DOUBLE,
					 double domain_max = na::DOUBLE)
			: Function(args, domain_min, domain_max) { set_eval_as_distr_call(); }
		virtual ~Distribution() = default;

		virtual std::unique_ptr<Distribution> clone_distrib() const = 0;
		std::unique_ptr<Function> clone() const override final { return clone_distrib(); };

		std::string err_name() const override { return "distribution '" + std::string(name()) + "'"; }

		//////////////////////
		//////  GETTER  //////
		//////////////////////

		double rand(rnd_gen &prng) const
		{
			switch (trunc_mode_)
			{
			case Trunc::None:
				return relocate_(cast_discr_(raw_rand_(prng)));
			case Trunc::Active:
				return relocate_(cast_discr_(inverse_sampling_(prng)));
			default:
				__builtin_unreachable();
			}
		}
		/// Return value of log-probability density/mass function.
		double log_pdf(double x) const { return log_func_(delocate_(x)); }
		/// Return value of log-ratio of probability density/mass functions.
		double log_ratio(double xo, double xn) const { return log_ratio_(delocate_(xo), delocate_(xn)); }

	protected:
		/////////////////////////////////////////////////////////////////////////////////////////////
		////////////////////////////  Distribution : virtual [protected] ////////////////////////////
		/////////////////////////////////////////////////////////////////////////////////////////////

		/// Raw random sampling method of the distribution.
		/// @param prng : pseudo-random number generator.
		virtual double raw_rand_(rnd_gen &) const { raise_undefined_call_err_("rand"); }

		/// Raw probability density/mass function of the distribution.
		virtual double raw_pdf_(double) const = 0;

		/// Raw log-probability density/mass function of the distribution.
		virtual double raw_log_pdf_(double) const { raise_undefined_call_err_("log_pdf"); }

		/// Raw cumulative distribution function.
		virtual double raw_cdf_(double) const = 0;

		/// Raw quantile.
		virtual double raw_quantile_(double) const { raise_undefined_call_err_("quantile"); }

		/// Raw method to call.
		double raw_function_(double) const override final { raise_undefined_call_err_("function"); }

		/// Polymorphic function method.
		double function_(double x) const override final
		{
			double y;
			switch (discr_mode_)
			{
			case Mode::Discrete::None:
				y = raw_pdf_(x);
				break;
			case Mode::Discrete::Center:
			{
				const double i = std::round(x);
				y = raw_cdf_(i + .5) - raw_cdf_(i - .5);
				break;
			}
			case Mode::Discrete::Left:
			{
				const double i = std::floor(x);
				y = raw_cdf_(i + 1.) - raw_cdf_(i);
				break;
			}
			case Mode::Discrete::Defined:
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
			case Mode::Discrete::None:
				y = raw_log_pdf_(x);
				break;
			case Mode::Discrete::Center:
			{
				const double i = std::round(x);
				y = raw_cdf_(i + 0.5) - raw_cdf_(i - 0.5);
				break;
			}
			case Mode::Discrete::Left:
			{
				const double i = std::floor(x);
				y = raw_cdf_(i + 1.) - raw_cdf_(i);
				break;
			}
			case Mode::Discrete::Defined:
				y = raw_log_pdf_(std::floor(x));
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
			case Mode::Discrete::None:
				adjusted_x = x;
				break;
			case Mode::Discrete::Center:
				adjusted_x = std::round(x) + 0.5;
				break;
			case Mode::Discrete::Left:
				adjusted_x = std::floor(x) + 1.0;
				break;
			case Mode::Discrete::Defined:
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
			case Mode::Discrete::None:
				return raw_log_ratio_(xo, xn);
			case Mode::Discrete::Center:
			{
				const double io = std::round(xo), in = std::round(xn);
				if (io == in)
					return 0.;
				return std::log(raw_cdf_(in + .5) - raw_cdf_(in - .5)) -
					   std::log(raw_cdf_(io + .5) - raw_cdf_(io - .5));
			}
			case Mode::Discrete::Left:
			{
				const double io = std::floor(xo), in = std::floor(xn);
				if (io == in)
					return 0.;
				return std::log(raw_cdf_(in + 1.) - raw_cdf_(in)) -
					   std::log(raw_cdf_(io + 1.) - raw_cdf_(io));
			}
			case Mode::Discrete::Defined:
				return raw_log_ratio_(std::floor(xo), std::floor(xn));
			default:
				__builtin_unreachable();
			}
		}

		/////////////////////////////////////////////////////////////////////////////////////////////
		////////////////////////////  Distribution : getter [protected] /////////////////////////////
		/////////////////////////////////////////////////////////////////////////////////////////////

		double cast_discr_(double x) const
		{
			switch (discr_mode_)
			{
			case Mode::Discrete::Center:
				return static_cast<double>(std::lround(x));
			case Mode::Discrete::Left:
				return std::floor(x);
			case Mode::Discrete::Defined:
			case Mode::Discrete::None:
				return x;
			default:
				__builtin_unreachable();
			}
		}

		double inverse_sampling_(rnd_gen &prng) const
		{
			thread_local static std::uniform_real_distribution<double> unif_gen(0., 1.);
			const double r = unif_gen(prng);
			const double sx = lower_.is_empty() ? 0. : (1. - r) * raw_cdf_(lower_.value());
			const double dx = upper_.is_empty() ? r : r * raw_cdf_(upper_.value());
			return raw_quantile_(sx + dx);
		}

		/////////////////////////////////////////////////////////////////////////////////////////////
		//////////////////////  Distribution : functional truncation [private] //////////////////////
		/////////////////////////////////////////////////////////////////////////////////////////////

		enum class Trunc : std::uint8_t
		{
			None,  /// No truncation
			Active /// Active truncation
		} trunc_mode_ = Trunc::None;

		/// Select functional method in case of truncation.
		void helper_init_distr_truncation_();

		/// UPDATE: if values don't change, I should just cache!
		/// UPDATE: truncation should become a switch too

		double (Distribution::*trunc_)() const = &Distribution::do_trunc_dummy_;

		// default case
		double do_trunc_dummy_() const { return 1.; }

		// continuos case
		double do_trunc_min_() const { return 1. - raw_cdf_(lower_.value()); }
		double do_trunc_max_() const { return raw_cdf_(upper_.value()); }
		double do_trunc_both_sides_() const { return raw_cdf_(upper_.value()) - raw_cdf_(lower_.value()); }

		// case of discrete defined
		double do_trunc_min_discr_defined_() const { return 1. - raw_cdf_(lower_.value() - 1); }
		double do_trunc_max_discr_defined_() const { return raw_cdf_(upper_.value() + 1); }
		double do_trunc_both_sides_discr_defined_() const { return raw_cdf_(upper_.value() + 1) - raw_cdf_(lower_.value() - 1); }

		// case of center discretization
		double do_trunc_min_discr_center_() const { return 1. - raw_cdf_(std::ceil(lower_.value()) - 0.5); }
		double do_trunc_max_discr_center_() const { return raw_cdf_(std::floor(upper_.value()) + 0.5); }
		double do_trunc_both_sides_discr_center_() const { return raw_cdf_(std::floor(upper_.value()) + 0.5) - raw_cdf_(std::ceil(lower_.value()) - 0.5); }

		// case of left discretization
		double do_trunc_min_discr_left_() const { return 1. - raw_cdf_(std::floor(lower_.value())); }
		double do_trunc_max_discr_left_() const { return raw_cdf_(std::ceil(upper_.value() + 1)); }
		double do_trunc_both_sides_discr_left_() const { return raw_cdf_(std::ceil(upper_.value() + 1)) - raw_cdf_(std::floor(lower_.value())); }

		double logtrunc_() const
		{
			switch (trunc_mode_)
			{
			case Trunc::None:
				return 0;
			case Trunc::Active:
				return std::log((this->*trunc_)());
			default:
				__builtin_unreachable();
			}
		}

		/////////////////////////////////////////////////////////////////////////////////////////////
		/////////////////////////  Distribution : functional shift [private] ////////////////////////
		/////////////////////////////////////////////////////////////////////////////////////////////

		double (Distribution::*shift_)() const = &Distribution::do_shift_dummy_;
		double do_shift_dummy_() const { return 0.; }
		double do_shift_min_() const { return raw_cdf_(lower_.value()); }
		double do_shift_min_discr_defined_() const { return raw_cdf_(lower_.value() - 1); }
		double do_shift_min_discr_center_() const { return raw_cdf_(std::floor(lower_.value()) - 0.5); }
		double do_shift_min_discr_left_() const { return raw_cdf_(std::floor(lower_.value())); }
	};

	/// Compose registries to get a distribution from them.
	std::unique_ptr<Distribution> compose_distr_registries(std::string_view opt,
														   const std::string &err_type,
														   const std::vector<std::string> &registries,
														   bool is_cumul = false);
}
