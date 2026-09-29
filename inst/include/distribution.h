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

//////////////////////////////  Distribution  //////////////////////////////

namespace housepi
{
	/// Base of every distribution.
	class Distribution : public Function
	{
	public:
		//////  Constructor  //////

		Distribution(const default_args &args = {},
					 double domain_min = na::DOUBLE,
					 double domain_max = na::DOUBLE)
			: Function(args, domain_min, domain_max)
		{
			set_eval_as_pdf();
			err_family_ = ErrType::Distribution;
		}
		virtual ~Distribution() = default;

		/// Return a copy, typed as a distribution.
		virtual std::unique_ptr<Distribution> clone_distrib() const = 0;
		/// Return a copy.
		std::unique_ptr<Function> clone() const override final { return clone_distrib(); };

		/// Finalize, then store the truncation normalizer when nothing it reads can move.
		void finalize() override
		{
			Function::finalize();
			if (trunc_mode_ == Trunc::None)
				logtrunc_mode_ = LogTrunc::None;
			else if (has_constant_args_())
				logtrunc_mode_ = LogTrunc::LookUp, lookup_log_trunc_ = std::log(trunc_());
			else
				logtrunc_mode_ = LogTrunc::Compute;
		}

		//////  Getter  //////

		/// Draw a value, truncation included.
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
		/// Return log-density or log-mass at `x`.
		double log_pdf(double x) const { return log_func_(delocate_(x)); }
		/// Return log-density or log-mass at `x`, minus infinity outside the domain.
		double guarded_log_pdf(double x) const { return is_valid(x) ? log_pdf(x) : -inf::DOUBLE; }

		/// Return the quantile at probability `p`, NA outside [0, 1].
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
		///
		/// - pdf() and cdf() return zero off the domain;
		/// - a density is unbounded above, a mass and a cumulative are not.
		Codomain codomain() const noexcept override
		{
			Codomain c;
			c.inf = 0.;
			c.below = Interval::Closed;
			if (eval_mode_ == Evaluate::Probability && !is_discretized())
				return c;
			c.sup = 1.;
			c.above = Interval::Closed;
			return c;
		}

	protected:
		//////  Virtual  //////

		/// Raw random draw.
		virtual double raw_rand_(rnd_gen & /*prng*/) const = 0;
		/// Raw density or mass.
		virtual double raw_pdf_(double /*x*/) const = 0;
		/// Raw log-density or log-mass.
		virtual double raw_log_pdf_(double /*x*/) const = 0;
		/// Raw cumulative distribution function.
		virtual double raw_cdf_(double /*x*/) const = 0;

		/// Raw complementary cumulative distribution, 1 - cdf.
		/// Discretization reads it where a tail-accurate closed form exists.
		virtual double raw_ccdf_(double x) const { return 1. - raw_cdf_(x); }

		/// Raw function: undefined for a distribution.
		double raw_function_(double /*x*/) const override final { raise_undefined_call_err_("function"); }

		/// Evaluate density or mass under the discretization scheme.
		double function_(double x) const override final
		{
			double y;
			switch (discr_mode_)
			{
			case Discrete::None:
				y = raw_pdf_(x);
				break;
			case Discrete::Center:
				y = discr_dist_center_(x);
				break;
			case Discrete::Left:
				y = discr_dist_left_(x);
				break;
			default: // Discrete::Defined:
				y = raw_pdf_(std::floor(x));
				break;
			}
			return y / trunc_();
		}

		/// Evaluate log-density or log-mass under the discretization scheme.
		double log_func_(double x) const override final
		{
			double y;
			switch (discr_mode_)
			{
			case Discrete::None:
				y = raw_log_pdf_(x);
				break;
			case Discrete::Center:
				y = std::log(discr_dist_center_(x));
				break;
			case Discrete::Left:
				y = std::log(discr_dist_left_(x));
				break;
			default: // Discrete::Defined:
				y = raw_log_pdf_(std::floor(x));
				break;
			}
			return y - logtrunc_();
		}

		/// Evaluate the truncated cumulative under the discretization scheme.
		double primitive_(double x) const override final { return (raw_cdf_(primitive_x_helper_(x)) - shift_()) / trunc_(); }

		/// Evaluate the log-ratio of densities under the discretization scheme.
		double log_ratio_(double xo, double xn) const override final
		{
			switch (discr_mode_)
			{
			case Discrete::None:
				return raw_log_ratio_(xo, xn);
			case Discrete::Center:
				return std::log(discr_dist_center_(xn)) - std::log(discr_dist_center_(xo));
			case Discrete::Left:
				return std::log(discr_dist_left_(xn)) - std::log(discr_dist_left_(xo));
			default: // Discrete::Defined:
				return raw_log_ratio_(std::floor(xo), std::floor(xn));
			}
		}

		//////  Getter  //////

		/// Lattice of the draws.
		/// - location off the integers: the support moves off with it.
		Lattice lattice_() const noexcept
		{
			if (discr_mode_ == Discrete::None)
				return Lattice::Doubles;
			if (domain_.loc.is_empty() || domain_.loc.codomain().on_integers())
				return Lattice::Integers;
			return Lattice::Doubles;
		}

		/// Round draw `x` onto the lattice of the discretization scheme.
		double cast_discr_(double x) const
		{
			switch (discr_mode_)
			{
			case Discrete::Center:
				return static_cast<double>(std::lround(x));
			case Discrete::Left:
				return std::floor(x);
			default: // cases Discrete::Defined, Discrete::None
				return x;
			}
		}

		/// Minimum cumulative window for rejection sampling.
		static constexpr double MIN_REJECT_WINDOW = 0.5;
		/// Rejection attempts before falling back to inverse sampling, succeeding with probability above 99%.
		static constexpr int MAX_REJECT_DRAWS = 8;

		/// Draw from the truncated distribution.
		///
		/// - window above MIN_REJECT_WINDOW: rejection sampling;
		/// - otherwise, or after MAX_REJECT_DRAWS rejections: inverse sampling.
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

		//////  Truncation  //////

		/// Truncated sides.
		enum class Trunc : std::uint8_t
		{
			None, /// No truncation.
			Min,  /// Lower-bound truncation.
			Max,  /// Upper-bound truncation.
			Both  /// Two-sided truncation.
		} trunc_mode_ = Trunc::None;
		/// Source of the log-normalizer.
		enum class LogTrunc : std::uint8_t
		{
			None,	/// No truncation, null log-normalizer.
			LookUp, /// Constant log-normalizer, stored by finalize().
			Compute /// Log-normalizer computed per evaluation.
		} logtrunc_mode_ = LogTrunc::Compute;
		/// Log-normalizer read in LogTrunc::LookUp mode.
		double lookup_log_trunc_ = na::DOUBLE;

		/// Select truncation mode.
		void helper_init_distr_truncation_()
		{
			const bool has_min = domain_.below != Interval::Unbounded;
			const bool has_max = domain_.above != Interval::Unbounded;

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
			case Discrete::Defined:
				return raw_cdf_(std::ceil(domain_.min.value()) - 1.);
			case Discrete::Center:
				return raw_cdf_(std::ceil(domain_.min.value()) - 0.5);
			case Discrete::Left:
				return raw_cdf_(std::ceil(domain_.min.value()));
			default: // case Discrete::None:
				return raw_cdf_(domain_.min.value());
			}
		}

		/// CDF at the maximum, with discretization correction.
		double cdf_max_() const
		{
			switch (discr_mode_)
			{
			case Discrete::Defined:
				return raw_cdf_(std::floor(domain_.max.value()));
			case Discrete::Center:
				return raw_cdf_(std::floor(domain_.max.value()) + 0.5);
			case Discrete::Left:
				return raw_cdf_(std::floor(domain_.max.value()) + 1.);
			default: // case Discrete::None:
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

		/// Log-normalizer of the truncated distribution.
		double logtrunc_() const
		{
			switch (logtrunc_mode_)
			{
			case LogTrunc::None:
				return 0.;
			case LogTrunc::LookUp:
				return lookup_log_trunc_;
			default: // case LogTrunc::Compute:
				return std::log(trunc_());
			}
		}

	private:
		//////  Discretization  //////

		/// Mass of the unit bin centered on `x`.
		double discr_dist_center_(double x) const
		{
			const double i = std::round(x);
			const double p = raw_cdf_(i + 0.5) - raw_cdf_(i - 0.5);
			// complementary cumulative only where the plain difference cancels in the right tail
			// the pmf stays finite, the bulk pays nothing
			return p > 0. ? p : raw_ccdf_(i - 0.5) - raw_ccdf_(i + 0.5);
		}
		/// Mass of the unit bin starting at floor(`x`).
		double discr_dist_left_(double x) const
		{
			const double i = std::floor(x);
			const double p = raw_cdf_(i + 1.) - raw_cdf_(i);
			return p > 0. ? p : raw_ccdf_(i) - raw_ccdf_(i + 1.);
		}
	};

	//////////////////////////////  NullDistribution  //////////////////////////////

	/// Point `x0` carrying no density: the model of a variable pinned to `x0`.
	///
	/// - log-density and log-ratio: 0, the variable never moves;
	/// - draws: `x0`.
	/// @note Unregistered: built in place of a collapsed model, never named by the user.
	class NullDistribution final : public Distribution
	{
	public:
		NullDistribution(double x0) : x0_(x0) {}
		std::unique_ptr<Distribution> clone_distrib() const override final { return std::make_unique<NullDistribution>(*this); }
		std::string name() const noexcept override final { return "null"; }

	protected:
		double raw_rand_(rnd_gen & /*prng*/) const override final { return x0_; }
		double raw_pdf_(double /*x*/) const override final { return 1.; }
		double raw_log_pdf_(double /*x*/) const override final { return 0.; }
		double raw_cdf_(double x) const override final { return x < x0_ ? 0. : 1.; }

	private:
		double x0_; /// Pinned point.
	};
}
