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
        /// Require arguments `idx` of the active parametrization to be positive.
        ArgRequirements positive_args_(std::initializer_list<size_t> idx) const;
        /// Whether arguments `idx` of the active parametrization are constant.
        bool are_args_constant_(std::initializer_list<size_t> idx) const;

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
            case Discrete::Left:
                return std::log(discr_func_left_(xn)) - std::log(discr_func_left_(xo));
            // Discrete::Defined:
            default:
                return raw_log_ratio_(std::floor(xo), std::floor(xn));
            }
        }

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
            case Discrete::Left:
                return std::floor(x);
            // cases Discrete::Defined, Discrete::None
            default:
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
        double truncated_rand_(RndGen &prng) const;

        /// Select truncation mode from the sides that cut mass.
        /// * side unbounded -> no cut;
        /// * side on a support edge -> bounds values, no cut.
        /// @param min_edge Whether the minimum sits on the support floor.
        /// @param max_edge Whether the maximum sits on the support ceiling.
        void helper_init_distr_truncation_(bool min_edge = false, bool max_edge = false);
    };

    //////////////////////////////  PinnedDistribution  //////////////////////////////

    /// Model of a variable pinned within a closed range, which it reports as its values.
    ///
    /// * single point -> point mass, as for a collapsed model;
    /// * interval -> support only, as for a parameter read from a column of values per repetition;
    /// * log-density and log-ratio -> 0, the variable never moves;
    /// * domain and support -> the range, lattice included;
    /// * value outside the range -> guarded log-density is minus infinity;
    /// * draws and quantiles -> the lower end.
    /// @note Unregistered: built in place of a collapsed model or of a pinned prior, never named by the user.
    class PinnedDistribution final : public Cloneable<PinnedDistribution, Distribution>
    {
    public:
        static constexpr const char *NAME = "pinned";
        PinnedDistribution(Codomain range)
            : Cloneable({}, range.inf, range.sup), range_(range) {}

        Codomain support() const noexcept override final { return range_; }

    protected:
        double raw_rand_(RndGen & /*prng*/) const override final { return range_.inf; }
        double raw_function_(double /*x*/) const override final { return 1.; }
        double raw_log_func_(double /*x*/) const override final { return 0.; }
        double raw_primitive_(double x) const override final { return x < range_.inf ? 0. : 1.; }
        double raw_quantile_(double /*p*/) const override final { return range_.inf; }
        double raw_log_ratio_(double /*xo*/, double /*xn*/) const override final { return 0.; }

    private:
        /// Values the pinned variable takes.
        Codomain range_;
    };
}
