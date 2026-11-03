/* housepi - Epidemiological Analysis of Household Outbreaks.
 * Copyright 2024-2026 Institut Pasteur (Paris), France.
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#pragma once

#include <cmath>
#include <cstdint>
#include <memory>
#include <string>
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
        ~Distribution() override = default;

        /// Return a copy, typed as a distribution.
        std::unique_ptr<Distribution> clone_distrib() const
        {
            return std::unique_ptr<Distribution>(static_cast<Distribution *>(clone().release()));
        }

        //////  Getter  //////

        /// Draw a value, truncation included.
        double rand(RndGen &prng) const
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
        /// Return log-ratio of densities or masses from `xo` to `xn`.
        /// * null move -> 0, whatever the density at that point.
        double log_ratio(double xo, double xn) const
        {
            return xo == xn ? 0. : log_ratio_(delocate_(xo), delocate_(xn));
        }
        /// Return log-ratio of densities or masses from `xo` to `xn`, validity checked first.
        ///
        /// * only `xn` valid -> infinity;
        /// * only `xo` valid -> minus infinity;
        /// * neither -> NA.
        double guarded_log_ratio(double xo, double xn) const
        {
            switch (is_valid(xo) + 2 * is_valid(xn))
            {
            case 3:
                return log_ratio(xo, xn);
            case 2:
                return inf::DOUBLE;
            case 1:
                return -inf::DOUBLE;
            default:
                return na::DOUBLE;
            }
        }

        /// Return the quantile at probability `p`, NA outside [0, 1].
        double quantile(double p) const
        {
            if (p < 0 || p > 1)
                return na::DOUBLE;
            if (p <= 0.)
                return domain_.min.is_empty() ? -inf::DOUBLE : relocate_(domain_.min.value());
            if (p >= 1.)
                return domain_.max.is_empty() ? inf::DOUBLE : relocate_(domain_.max.value());
            if (trunc_mode_ != Trunc::None)
                p = trunc_mode_ == Trunc::Max ? p * cdf_max_() : cdf_min_() + p * trunc_();
            return relocate_(cast_discr_(raw_quantile_(p)));
        }

        /// Values the distribution draws can take.
        virtual Codomain support() const noexcept { return domain_.relocated(lattice_()); }

        /// Values returned by evaluate().
        ///
        /// * pdf() and cdf() return zero off the domain;
        /// * a density is unbounded above, a mass and a cumulative are not.
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
        //////  Arguments  //////

        /// Require arguments `idx` of the active parametrization to be positive.
        ArgRequirements positive_args_(std::initializer_list<size_t> idx) const
        {
            ArgRequirements reqs;
            for (size_t i : idx)
                reqs.emplace_back(&arguments_(i), Codomain::positive());
            return reqs;
        }

        //////  Virtual  //////

        /// Raw quantile at cumulative probability `p`.
        virtual double raw_quantile_(double /*p*/) const = 0;
        /// Raw log-ratio of densities or masses from `xo` to `xn`.
        virtual double raw_log_ratio_(double /*xo*/, double /*xn*/) const = 0;
        /// Raw random draw.
        virtual double raw_rand_(RndGen & /*prng*/) const = 0;

        /// Evaluate the log-ratio of densities or masses under the discretization scheme.
        double log_ratio_(double xo, double xn) const
        {
            switch (discr_mode_)
            {
            case Discrete::None:
                return raw_log_ratio_(xo, xn);
            case Discrete::Center:
                return std::log(discr_func_center_(xn)) - std::log(discr_func_center_(xo));
            case Discrete::Left:
                return std::log(discr_func_left_(xn)) - std::log(discr_func_left_(xo));
            default: // Discrete::Defined:
                return raw_log_ratio_(std::floor(xo), std::floor(xn));
            }
        }

        //////  Getter  //////

        /// Lattice of the draws.
        /// * location off the integers -> the support moves off with it.
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
        /// * window above MIN_REJECT_WINDOW -> rejection sampling;
        /// * otherwise, or after MAX_REJECT_DRAWS rejections -> inverse sampling.
        double truncated_rand_(RndGen &prng) const
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
    };

    //////////////////////////////  NullDistribution  //////////////////////////////

    /// Closed range carrying no density: the model of a variable pinned within it.
    ///
    /// * log-density and log-ratio -> 0, the variable never moves;
    /// * support -> the range, lattice included;
    /// * draws -> the lower end, the pinned point of a collapsed model.
    /// @note Unregistered: built in place of a collapsed model or of a pinned prior, never named by the user.
    class NullDistribution final : public Cloneable<NullDistribution, Distribution>
    {
    public:
        static constexpr const char *NAME = "null";
        NullDistribution(Codomain range) : range_(range) {}

        Codomain support() const noexcept override final { return range_; }

    protected:
        double raw_rand_(RndGen & /*prng*/) const override final { return range_.inf; }
        double raw_function_(double /*x*/) const override final { return 1.; }
        double raw_log_func_(double /*x*/) const override final { return 0.; }
        double raw_primitive_(double x) const override final { return x < range_.inf ? 0. : 1.; }
        double raw_quantile_(double /*p*/) const override final { return range_.inf; }
        double raw_log_ratio_(double /*xo*/, double /*xn*/) const override final { return 0.; }

    private:
        Codomain range_; /// Values the pinned variable takes.
    };
}
