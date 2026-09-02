/* housepi - Epidemiological Analysis of Household Outbreaks.
 *
 * Copyright © 2024-2026 Institut Pasteur (Paris), CNRS and INSERM.
 * Distributed under the terms of the GNU General Public License (GPLv3).
 * See the LICENSE file for details.
 */

#pragma once

#include <cmath>
#include <cstdint>
#include <memory>
#include <random>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

#include "tools/messages.h"
#include "tools/nan.h"
#include "tools/random.h"

#include "function.h"

/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////// *** Distribution : class *** ///////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

namespace housepi
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
			: Function(args, domain_min, domain_max)
		{
			set_eval_as_pdf();
			err_family_ = ErrType::Distribution;
		}
		virtual ~Distribution() = default;

		virtual std::unique_ptr<Distribution> clone_distrib() const = 0;
		std::unique_ptr<Function> clone() const override final { return clone_distrib(); };

		//////////////////////
		//////  GETTER  //////
		//////////////////////

		double rand(rnd_gen &prng) const
		{
			switch (trunc_mode_)
			{
			case Trunc::None:
				return relocate_(cast_discr_(raw_rand_(prng)));
			default:
				return relocate_(cast_discr_(truncated_rand_(prng)));
			}
		}
		/// Return value of log-probability density/mass function.
		double log_pdf(double x) const { return log_func_(delocate_(x)); }
		/// Check for validity, then return value of log-probability density/mass function.
		double guarded_log_pdf(double x) const { return is_valid(x) ? log_pdf(x) : -inf::DOUBLE; }

		/// Check for validity, then  return the quantile (inverse CDF) at probability p.
		double quantile(double p) const
		{
			if (p < 0 || p > 1)
				return na::DOUBLE;
			if (p <= 0.)
				return domain_.min.is_empty() ? -inf::DOUBLE : relocate_(domain_.min.value());
			if (p >= 1.)
				return domain_.max.is_empty() ? inf::DOUBLE : relocate_(domain_.max.value());
			return relocate_(cast_discr_(raw_quantile_(shift_() + p * trunc_())));
		}

		/// Values the distribution draws can take.
		Codomain support() const noexcept { return domain_.relocated(lattice_()); }

		/// Values returned by evaluate().
		/// NOTE: pdf() and cdf() return zero off the domain.
		/// NOTE: density is unbounded above, a mass and a cumulative are not.
		Codomain codomain() const noexcept override
		{
			Codomain c;
			c.inf = 0.;
			c.below = Interval::Closed;
			if (eval_mode_ == Mode::Evaluate::Probability && !is_discretized())
				return c;
			c.sup = 1.;
			c.above = Interval::Closed;
			return c;
		}

	protected:
		/// Integer-value structure of the draws.
		/// NOTE: location off the integers moves the support off as well.
		Lattice lattice_() const noexcept
		{
			if (discr_mode_ == Mode::Discrete::None)
				return Lattice::Doubles;
			if (domain_.loc.is_empty() || domain_.loc.codomain().on_integers())
				return Lattice::Integers;
			return Lattice::Doubles;
		}

		/////////////////////////////////////////////////////////////////////////////////////////////
		////////////////////////////  Distribution : virtual [protected] ////////////////////////////
		/////////////////////////////////////////////////////////////////////////////////////////////

		/// Raw random sampling method of the distribution.
		/// @param prng : pseudo-random number generator.
		virtual double raw_rand_(rnd_gen & /*prng*/) const { raise_undefined_call_err_("rand"); }

		/// Raw probability density/mass function of the distribution.
		virtual double raw_pdf_(double /*x*/) const = 0;

		/// Raw log-probability density/mass function of the distribution.
		virtual double raw_log_pdf_(double /*x*/) const { raise_undefined_call_err_("log_pdf"); }

		/// Raw cumulative distribution function.
		virtual double raw_cdf_(double /*x*/) const = 0;

		/// Raw complementary CDF (1 - cdf).
		/// NOTE: used in discretization when a tail-accurate closed form exists.
		virtual double raw_ccdf_(double x) const { return 1. - raw_cdf_(x); }

		/// Raw method to call.
		double raw_function_(double /*x*/) const override final { raise_undefined_call_err_("function"); }

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
				y = discr_dist_center_(x);
				break;
			case Mode::Discrete::Left:
				y = discr_dist_left_(x);
				break;
			default: // Mode::Discrete::Defined:
				y = raw_pdf_(std::floor(x));
				break;
			}
			return y / trunc_();
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
				y = std::log(discr_dist_center_(x));
				break;
			case Mode::Discrete::Left:
				y = std::log(discr_dist_left_(x));
				break;
			default: // Mode::Discrete::Defined:
				y = raw_log_pdf_(std::floor(x));
				break;
			}
			return y - logtrunc_();
		}

		/// Polymorphic primitive method.
		double primitive_(double x) const override final { return (raw_cdf_(primitive_x_helper_(x)) - shift_()) / trunc_(); }

		/// Polymorphic log_ratio method.
		double log_ratio_(double xo, double xn) const override final
		{
			switch (discr_mode_)
			{
			case Mode::Discrete::None:
				return raw_log_ratio_(xo, xn);
			case Mode::Discrete::Center:
				return std::log(discr_dist_center_(xn)) - std::log(discr_dist_center_(xo));
			case Mode::Discrete::Left:
				return std::log(discr_dist_left_(xn)) - std::log(discr_dist_left_(xo));
			default: // Mode::Discrete::Defined:
				return raw_log_ratio_(std::floor(xo), std::floor(xn));
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
			default: // cases Mode::Discrete::Defined, Mode::Discrete::None
				return x;
			}
		}

		/// Minimum CDF window for rejection sampling.
		static constexpr double MIN_REJECT_WINDOW = 0.5;
		/// Rejection attempts before falling back to inverse sampling.
		/// NOTE: succeeds with > 99% proba
		static constexpr int MAX_REJECT_DRAWS = 8;

		/// Sample from the truncated distribution.
		/// When truncation window is large, uses rejection-sampling.
		/// Otherwise falls back to inverse sampling.
		double truncated_rand_(rnd_gen &prng) const
		{
			const double cdf_at_min = trunc_mode_ == Trunc::Max ? 0. : cdf_min_();
			const double cdf_at_max = trunc_mode_ == Trunc::Min ? 1. : cdf_max_();
			const double window = cdf_at_max - cdf_at_min;

			if (window > MIN_REJECT_WINDOW)
			{
				const double min_val = trunc_mode_ == Trunc::Max ? -inf::DOUBLE : domain_.min.value();
				const double max_val = trunc_mode_ == Trunc::Min ? inf::DOUBLE : domain_.max.value();
				for (int k = 0; k < MAX_REJECT_DRAWS; ++k)
				{
					const double x = cast_discr_(raw_rand_(prng));
					if (x >= min_val && x <= max_val)
						return x;
				}
			}
			return raw_quantile_(cdf_at_min + standard_uniform(prng) * window);
		}

		/////////////////////////////////////////////////////////////////////////////////////////////
		//////////////////////  Distribution : functional truncation [private] //////////////////////
		/////////////////////////////////////////////////////////////////////////////////////////////

		enum class Trunc : std::uint8_t
		{
			None, /// No truncation
			Min,  /// Lower-bound truncation
			Max,  /// Upper-bound truncation
			Both  /// Two-sided truncation
		} trunc_mode_ = Trunc::None;

		/// Select truncation mode.
		void helper_init_distr_truncation_()
		{
			const bool has_min = domain_.below != Mode::Interval::Unbounded;
			const bool has_max = domain_.above != Mode::Interval::Unbounded;

			if (!has_min && !has_max)
				trunc_mode_ = Trunc::None;
			else if (has_min && !has_max)
				trunc_mode_ = Trunc::Min;
			else if (!has_min && has_max)
				trunc_mode_ = Trunc::Max;
			else
				trunc_mode_ = Trunc::Both;
		}

		/// CDF at the minimum, with discretization correction.
		double cdf_min_() const
		{
			switch (discr_mode_)
			{
			case Mode::Discrete::Defined:
				return raw_cdf_(std::ceil(domain_.min.value()) - 1.);
			case Mode::Discrete::Center:
				return raw_cdf_(std::ceil(domain_.min.value()) - 0.5);
			case Mode::Discrete::Left:
				return raw_cdf_(std::ceil(domain_.min.value()));
			default: // case Mode::Discrete::None:
				return raw_cdf_(domain_.min.value());
			}
		}

		/// CDF at the maximum, with discretization correction.
		double cdf_max_() const
		{
			switch (discr_mode_)
			{
			case Mode::Discrete::Defined:
				return raw_cdf_(std::floor(domain_.max.value()));
			case Mode::Discrete::Center:
				return raw_cdf_(std::floor(domain_.max.value()) + 0.5);
			case Mode::Discrete::Left:
				return raw_cdf_(std::floor(domain_.max.value()) + 1.);
			default: // case Mode::Discrete::None:
				return raw_cdf_(domain_.max.value());
			}
		}

		/// Normalizing constant for the truncated distribution.
		double trunc_() const
		{
			switch (trunc_mode_)
			{
			case Trunc::Min:
				return 1. - cdf_min_();
			case Trunc::Max:
				return cdf_max_();
			case Trunc::Both:
				return cdf_max_() - cdf_min_();
			default: // case Trunc::None:
				return 1.;
			}
		}

		/// CDF offset subtracted by the primitive for truncated distributions.
		double shift_() const
		{
			switch (trunc_mode_)
			{
			case Trunc::Min:
			case Trunc::Both:
				return cdf_min_();
			default: // cases Trunc::None, Trunc::Max
				return 0.;
			}
		}

		double logtrunc_() const
		{
			switch (trunc_mode_)
			{
			case Trunc::None:
				return 0.;
			default:
				return std::log(trunc_());
			}
		}

	private:
		/// Helper for Center discretization.
		double discr_dist_center_(double x) const
		{
			const double i = std::round(x);
			const double p = raw_cdf_(i + 0.5) - raw_cdf_(i - 0.5);
			// fall back to the complementary CDF only when the plain difference
			// cancels in the right tail (keeps the pmf finite); the bulk pays nothing.
			return p > 0. ? p : raw_ccdf_(i - 0.5) - raw_ccdf_(i + 0.5);
		}
		/// Helper for Left discretization.
		double discr_dist_left_(double x) const
		{
			const double i = std::floor(x);
			const double p = raw_cdf_(i + 1.) - raw_cdf_(i);
			return p > 0. ? p : raw_ccdf_(i) - raw_ccdf_(i + 1.);
		}
	};

}
