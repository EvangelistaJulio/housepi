/* housepi - Epidemiological Analysis of Household Outbreaks.
 * Copyright 2024-2026 Institut Pasteur (Paris), France.
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#pragma once

#include <cmath>
#include <cstdint>
#include <memory>
#include <string>
#include <string_view>
#include <type_traits>
#include <utility>
#include <vector>

#include "tools/messages.h"
#include "tools/nan.h"

namespace housepi
{
    class Function;

    /// Index for the global name pool of the process.
    using NameIdx = std::uint32_t;

    /// Openness of an interval side.
    enum class Interval : std::uint8_t
    {
        Unbounded,
        Open,
        Closed
    };

    /// Integer-value structure of the values returned.
    enum class Lattice : std::uint8_t
    {
        Doubles,
        Integers
    };

    /// Values a definition returns, in the absolute frame.
    struct Codomain
    {
        /// Infimum.
        double inf = -inf::DOUBLE;
        /// Supremum.
        double sup = inf::DOUBLE;
        /// Openness of the lower side.
        Interval below = Interval::Unbounded;
        /// Openness of the upper side.
        Interval above = Interval::Unbounded;
        /// Lattice of the values.
        Lattice lattice = Lattice::Doubles;

        /// Values greater than zero.
        static Codomain positive(Lattice l = Lattice::Doubles) noexcept
        {
            return {0., inf::DOUBLE, Interval::Open, Interval::Unbounded, l};
        }
        /// Values not below zero.
        static Codomain non_negative(Lattice l = Lattice::Doubles) noexcept
        {
            return {0., inf::DOUBLE, Interval::Closed, Interval::Unbounded, l};
        }
        /// Values between zero and one, each side open or closed.
        static Codomain unit(Interval below, Interval above) noexcept
        {
            return {0., 1., below, above};
        }
        /// Finite values, `-Inf` and `Inf` excluded.
        static Codomain finite(Lattice l = Lattice::Doubles) noexcept
        {
            return {-inf::DOUBLE, inf::DOUBLE, Interval::Open, Interval::Open, l};
        }

        /// Whether every value returned is an integer.
        bool on_integers() const noexcept { return lattice == Lattice::Integers; }
        /// Whether every value lies within `outer`, its lattice included.
        /// An unbounded side holds finite values: open at -Inf or Inf.
        bool is_within(const Codomain &outer) const noexcept;

        /// Interval notation, in the absolute frame.
        std::string str() const;

        /// Whether every value returned is provably greater than zero.
        /// @note One-sided: false does not mean negative.
        bool is_positive() const noexcept
        {
            return below != Interval::Unbounded &&
                   (inf > 0. || (inf == 0. && below == Interval::Open));
        }
    };

    //////////////////////////////  Attribute  //////////////////////////////

    /// Hold one value of a definition: empty, a number, a pointer to another attribute, or a function.
    class Attribute
    {
    public:
        //////  Constructor  //////

        Attribute();
        Attribute(double x) { init_numeric(x); }
        virtual ~Attribute();

        /// Copy constructor.
        Attribute(const Attribute &other);
        /// Move constructor.
        Attribute(Attribute &&) noexcept;
        /// Copy assignment.
        Attribute &operator=(Attribute other);

        /// Swap contents with `other`.
        void swap(Attribute &other) noexcept;

        //////  Operators  //////

        /// Whether both attributes are equal.
        bool operator==(const Attribute &other) const noexcept;
        /// Whether the attributes differ.
        bool operator!=(const Attribute &other) const noexcept { return !(*this == other); }

        //////  Setter  //////

        /// Set attribute name.
        void name(const std::string &txt);
        void name(const std::string &txt, ErrType err_type)
        {
            this->name(txt);
            err_type_ = err_type;
        }

        /// Initialize attribute as numeric.
        void init_numeric(double x);
        /// Initialize attribute as pointer.
        void init_pointer(const Attribute *y);
        /// Initialize attribute as function.
        void init_function(std::unique_ptr<Function> f);
        /// Initialize attribute as the product of the attributes in `other_ptrs`.
        void init_product(const std::vector<const Attribute *> &other_ptrs);
        /// Transfer definition from `other`.
        void transfer(const Attribute &other);
        /// Finalize: prune indirections, collapse constants to numeric.
        /// @note Call only when every read flagged constant keeps its value for the whole run.
        virtual void finalize();

        /// Set attribute value.
        void set(double x)
        {
            switch (mode_)
            {
            case Mode::Empty:
                mode_ = Mode::Numeric;
                value_ = x;
                return;
            case Mode::Numeric:
                value_ = x;
                return;
            // cases Mode::Pointer, Mode::Function
            default:
                msg::stop("Tried to change non numeric " + str() + ".");
            }
        }

        /// Initialize argument `key` of the function from `attr`.
        /// @note Call only in function mode.
        void init_func_attr(std::string_view key, const Attribute &attr);

        //////  Getter  //////

        /// Values the attribute takes.
        virtual Codomain codomain() const noexcept;

        /// Get numeric value.
        /// @note Defined after Function, which function mode evaluates.
        double value() const;

        /// Get the attribute holding the value, pointers followed.
        const Attribute *pointer() const noexcept
        {
            switch (mode_)
            {
            case Mode::Pointer:
                return attr_->pointer();
            // cases Mode::Empty, Mode::Numeric, Mode::Function
            default:
                return this;
            }
        }
        /// Get the function, nullptr unless in function mode.
        /// @note Pointers are not followed: call pointer() first.
        const Function *function() const noexcept { return func_.get(); }

        //////  Flags  //////

        /// Whether it holds no value.
        /// * function -> never missing, by convention;
        /// * reduction function -> may still return NA.
        bool is_na() const noexcept;
        /// Whether it is constant.
        bool is_constant() const noexcept;
        /// Whether every value it returns is an integer.
        bool is_discrete() const noexcept { return codomain().on_integers(); }
        /// Whether every value it returns is positive.
        bool is_positive() const noexcept { return codomain().is_positive(); }
        /// Whether it can return a value.
        bool is_callable() const noexcept;

        /// Whether it has no definition.
        bool is_empty() const noexcept { return mode_ == Mode::Empty; }
        /// Whether it is a number.
        bool is_numeric() const noexcept { return mode_ == Mode::Numeric; }
        /// Whether it is a pointer.
        bool is_pointer() const noexcept { return mode_ == Mode::Pointer; }
        /// Whether it is a function.
        bool is_function() const noexcept { return mode_ == Mode::Function; }

        //////  Show  //////

        /// Get attribute name.
        std::string name() const;
        /// Name formatted for error messages.
        std::string str() const;
        /// Source formatted for error messages.
        /// * number -> its value;
        /// * pointer -> the attribute it reads, e.g. `data "lo"`;
        /// * function, or pointer to one -> the definition of the attribute holding it, e.g. "`sd` definition".
        std::string str_source() const;

    protected:
        /// Storage mode.
        enum class Mode : std::uint8_t
        {
            Empty,
            Numeric,
            Pointer,
            Function
        } mode_ = Mode::Empty;

        /// Type name displayed in errors.
        ErrType err_type_ = ErrType::Empty;
        /// Whether the user declared it constant.
        bool is_user_const_ = false;
        /// Interned name id, 0 when empty.
        NameIdx name_id_ = 0;
        double value_ = na::DOUBLE;
        /// Attribute read through, owned by a module of the same chain.
        /// * chain rebuilt -> pointed again by its build, never copied.
        const Attribute *attr_ = nullptr;
        std::unique_ptr<Function> func_ = nullptr;
    };

    //////////////////////////////  Domain  //////////////////////////////

    /// Values a definition may be called on.
    /// Bounds are live attributes.
    struct Domain
    {
        /// User-defined minimum.
        Attribute min = Attribute();
        /// User-defined maximum.
        Attribute max = Attribute();
        /// Shift between the frames.
        Attribute loc = Attribute();
        /// Openness of the lower side.
        Interval below = Interval::Unbounded;
        /// Openness of the upper side.
        Interval above = Interval::Unbounded;

        /// Cast `x` into the frame the bounds and all raw math live in.
        /// @note Read per evaluation, so the guard stays a mode compare.
        double delocate(double x) const { return loc.is_empty() ? x : x - loc.value(); }
        /// Cast `x` back into the frame the caller sees.
        double relocate(double x) const { return loc.is_empty() ? x : x + loc.value(); }

        /// Whether `x` falls below the minimum.
        bool is_below(double x) const
        {
            switch (below)
            {
            case Interval::Open:
                return delocate(x) <= min.value();
            case Interval::Closed:
                return delocate(x) < min.value();
            // case Interval::Unbounded:
            default:
                return false;
            }
        }

        /// Whether `x` falls above the maximum.
        bool is_above(double x) const
        {
            switch (above)
            {
            case Interval::Open:
                return delocate(x) >= max.value();
            case Interval::Closed:
                return delocate(x) > max.value();
            // case Interval::Unbounded:
            default:
                return false;
            }
        }

        /// Whether `x` falls inside.
        bool is_valid(double x) const { return !is_below(x) && !is_above(x); }
        /// The interval the bounds delimit, in the absolute frame.
        /// Bounds read through their codomains: an undrawn bound still places the interval it reaches.
        Codomain relocated() const;
        /// Snapshot into the absolute frame, for the support of a distribution.
        /// A dense support puts no mass on its floor: a closed bound reads open there.
        /// Boundary attainability is decided here.
        Codomain relocated(Lattice lattice) const;

        /// Set the openness of both sides from the bounds provided.
        /// @param owner Name of the definition, for the error it may raise.
        void init_standard(bool is_min_open, bool is_max_open, const std::string &owner);
        /// Lower the minimum onto `min_val`, returning whether truncation is needed.
        bool init_min(double min_val, bool is_min_open, const std::string &owner);

        /// Whether both bounds read and delimit a non-empty interval.
        /// @todo Account for the codomains and constantness of the bounds.
        bool is_consistent() const;
        /// Whether both bounds read and coincide on a single point.
        bool is_degenerate() const { return min.is_callable() && max.is_callable() && min.value() == max.value(); }

        /// Raise a uniform error if the domain is not consistent.
        /// @param owner Name of the definition, for the error it may raise.
        void assert_consistent(const std::string &owner) const;

        /// Interval notation of the values the bounds reach, in the absolute frame.
        std::string str() const { return relocated().str(); }
    };

    //////////////////////////////  Function  //////////////////////////////

    /// Base of every function and distribution, with arguments stored in arrays.
    class Function
    {
    protected:
        /// Method an implicit evaluate() dispatches to.
        enum class Evaluate : std::uint8_t
        {
            Function,
            Primitive,
            Probability,
            Cumulative,
            Fold
        };
        /// Discretization scheme.
        enum class Discrete : std::uint8_t
        {
            None,
            Left,
            Defined
        };

    public:
        //////  Constructor  //////

        using default_args = std::vector<std::pair<std::string, double>>;
        Function(const default_args &args = {},
                 double domain_min = na::DOUBLE,
                 double domain_max = na::DOUBLE);
        virtual ~Function() = default;

        /// Return a copy.
        virtual std::unique_ptr<Function> clone() const = 0;

        //////  Setter  //////

        /// Initialize discretization.
        /// @todo Make discretize() and truncate() protected.
        virtual void discretize() { raise_undefined_call_err_("discretize"); }
        /// Initialize truncation.
        /// @note Call only after discretize().
        virtual void truncate() { helper_init_standard_domain_(); };

        /// Discretize, truncate, then assert consistency.
        /// @note Call once every argument and the domain are assigned.
        void reshape(bool discr)
        {
            if (discr)
                discretize();
            truncate();

            assert_is_callable_explicitly();
            assert_location_matches_discrete();
            assert_arg_requirements();
        }

        /// Select primitive as implicit evaluate.
        void set_eval_as_prim() noexcept { eval_mode_ = Evaluate::Primitive; }
        /// Select pdf as implicit evaluate, the default for distributions.
        void set_eval_as_pdf() noexcept { eval_mode_ = Evaluate::Probability; }
        /// Select cdf as implicit evaluate.
        void set_eval_as_cdf() noexcept { eval_mode_ = Evaluate::Cumulative; }
        /// Select an NA-skipping fold, for reductions over household members.
        void set_eval_as_fold() noexcept { eval_mode_ = Evaluate::Fold; }

        /// Set argument `i` to `x`.
        /// @todo Narrow access to set_arg(), set_min(), set_max(): AdaptableProposal is their only caller.
        void set_arg(size_t i, double x) { alts_[alt_idx_].args[i].set(x); }
        /// Set the minimum.
        void set_min(double x) { domain_.min.set(x); }
        /// Set the maximum.
        void set_max(double x) { domain_.max.set(x); }

        /// Transfer live bound `a` into the minimum or maximum truncation.
        /// Pointers and functions stay live: the domain tracks the source per evaluation.
        void transfer_min(const Attribute &a) { domain_.min.transfer(a); }
        void transfer_max(const Attribute &a) { domain_.max.transfer(a); }

        /// Initialize argument `key` from `attr`.
        void init_attr(std::string_view key, const Attribute &attr);
        /// Leave argument `key` unbound, so no requirement is judged on it.
        void init_attr_unbound(std::string_view key) { find_arg_(key).init_numeric(na::DOUBLE); }

        /// Finalize: prune indirections, collapse constants to numeric,
        /// store the truncation normalizer when nothing it reads can move.
        virtual void finalize();

        //////  Getter  //////

        /// Values returned, composed from the operands.
        ///
        /// * computed live, never stored: finalize() may fold the operands after reshape();
        /// * defaulted, so an out-of-tree custom function still compiles.
        virtual Codomain codomain() const noexcept { return Codomain(); }

        /// Number of arguments in the *active* parametrization.
        size_t n_args() const noexcept { return arguments_().size(); }
        /// Read value of argument `i` of the *active* parametrization.
        double arg(size_t i) const { return alts_[alt_idx_].args[i].value(); }

        /// Location, minimum and maximum of the domain.
        const Attribute &loc() const noexcept { return domain_.loc; }
        const Attribute &min() const noexcept { return domain_.min; }
        const Attribute &max() const noexcept { return domain_.max; }
        /// Values `x` may take, in the de-located frame.
        const Domain &domain() const noexcept { return domain_; }

        /// Evaluate at the bound argument `x`.
        double evaluate() const { return evaluate(x_arg_.value()); }

        /// Evaluate at `x`, according to the selected mode.
        double evaluate(double x) const
        {
            switch (eval_mode_)
            {
            case Evaluate::Primitive:
                return primitive_(delocate_(x));
            case Evaluate::Probability:
                return pdf(x);
            case Evaluate::Cumulative:
                return cdf(x);
            // cases Evaluate::Function, Evaluate::Fold
            default:
                return function_(delocate_(x));
            }
        }

        /// Return function value at `x`, delocated.
        double func(double x) const { return function_(delocate_(x)); }
        /// Return log-function value at `x`, delocated.
        double log_func(double x) const { return log_func_(delocate_(x)); }
        /// Return primitive value at `x`, delocated.
        double prim(double x) const { return primitive_(delocate_(x)); }

        /// Return density or mass at `x`, 0 outside the domain.
        double pdf(double x) const
        {
            return is_valid(x) ? function_(delocate_(x)) : 0.;
        }
        /// Return cumulative distribution at `x`.
        /// * above the domain -> primitive at the maximum, which lives in the de-located frame.
        double cdf(double x) const
        {
            return is_below(x) ? 0 : primitive_(is_above(x) ? domain_.max.value() : delocate_(x));
        }

        //////  Flags  //////

        /// Whether the return value is constant.
        /// * mutable domain -> not constant.
        bool is_constant() const noexcept
        {
            return x_arg_.is_constant() && has_constant_args_();
        }

        /// Whether discretize() put the value on a unit lattice, of unspecified offset.
        ///
        /// * location off the integers -> the lattice shifts with it;
        /// * a discretized definition may return values off the integers.
        /// @note For "the values are integers", ask codomain().on_integers().
        bool is_discretized() const noexcept { return discr_mode_ != Discrete::None; }
        /// Whether the return value is always strictly positive.
        bool is_positive() const noexcept { return codomain().is_positive(); }
        /// Whether no argument is missing, see get_missing_argument().
        bool is_callable_explicitly() const noexcept;
        /// Whether it can be called implicitly, with `x` bound too.
        bool is_callable_implicitly() const noexcept { return !x_arg_.is_empty() && is_callable_explicitly(); }

        /// Collect names of missing arguments.
        /// * argument of the active parametrization -> missing when empty;
        /// * domain bound -> missing when empty and listed by arg_requirements_().
        std::vector<std::string> get_missing_argument() const;
        /// Throw if not callable explicitly.
        void assert_is_callable_explicitly() const;
        /// Throw if discretization and location do not match.
        void assert_location_matches_discrete() const;
        /// Throw if an argument can take a value outside its requirement.
        /// * constant -> its value is judged;
        /// * otherwise -> every value its source can take, from its codomain.
        /// @note Unbound and empty arguments are not judged.
        void assert_arg_requirements() const;

        /// Whether `x` is below the domain minimum.
        bool is_below(double x) const { return domain_.is_below(x); }
        /// Whether `x` is above the domain maximum.
        bool is_above(double x) const { return domain_.is_above(x); }
        /// Whether `x` is inside the domain.
        bool is_valid(double x) const { return domain_.is_valid(x); }

        //////  Show  //////

        /// Get function name.
        virtual std::string name() const noexcept = 0;
        /// Get name of argument `i`, from the active parametrization.
        std::string arg_name(unsigned short i) const { return alts_[alt_idx_].args[i].name(); }
        /// Name formatted for error messages.
        virtual std::string str() const;
        /// Set the family word shown in errors, see str().
        void set_err_family(ErrType family) noexcept { err_family_ = family; }

    protected:
        //////  Domain  //////

        /// Values `x` may take, and their frame.
        Domain domain_ = Domain();

        /// Cast into the de-located frame.
        double delocate_(double x) const { return domain_.delocate(x); }
        /// Cast back into the absolute frame.
        double relocate_(double x) const { return domain_.relocate(x); }

        /// Whether the active arguments and the bounded sides of the domain are constant.
        bool has_constant_args_() const noexcept;

        /// Set the standard domain, open or closed on each side.
        /// @todo Avoid computing str() on every call.
        void helper_init_standard_domain_(bool is_min_open = false, bool is_max_open = false)
        {
            domain_.init_standard(is_min_open, is_max_open, str());
        }
        /// Set a domain with minimum `min_val` and no maximum.
        bool helper_init_domain_min_(double min_val = 0., bool is_min_open = false)
        {
            return domain_.init_min(min_val, is_min_open, str());
        }

        //////  Arguments  //////

        Attribute x_arg_ = Attribute();

        /// Arguments of one parametrization.
        struct AltArgSet
        {
            /// Parametrization name.
            std::string name = {};
            std::vector<Attribute> args = {};
        };
        /// Alternative parametrizations, 0 canonical.
        std::vector<AltArgSet> alts_;
        /// Active parametrization.
        std::uint8_t alt_idx_ = 0;
        /// Whether the user picked a re-parametrization.
        bool is_user_alt_ = false;
        /// Family word shown in errors.
        ErrType err_family_ = ErrType::Function;
        /// Method an implicit evaluate() dispatches to.
        Evaluate eval_mode_ = Evaluate::Function;
        /// Discretization scheme.
        Discrete discr_mode_ = Discrete::None;

        /// Register an alternative parametrization.
        void register_alt_(const std::string &name, const default_args &args);

        /// Assert argument `key` accepts `attr`.
        void helper_assert_valid_attribute_(std::string_view key, const Attribute &attr);

        /// Arguments of the active parametrization.
        std::vector<Attribute> &arguments_() noexcept { return alts_[alt_idx_].args; }
        const std::vector<Attribute> &arguments_() const noexcept { return alts_[alt_idx_].args; }
        Attribute &arguments_(size_t i) noexcept { return alts_[alt_idx_].args[i]; }
        const Attribute &arguments_(size_t i) const noexcept { return alts_[alt_idx_].args[i]; }

        /// Return argument `arg_name`, an empty attribute if absent.
        virtual Attribute &find_arg_(std::string_view arg_name);

        /// Append argument `key`, for variadic parametrizations.
        /// * standard arguments `x`, `loc`, `min`, `max` -> returned, not appended.
        /// @note Call once per key: repeated names are refused R-side.
        Attribute &append_arg_(std::string_view key);

        //////  Virtual  //////

        /// Throw on a call to virtual method `call_name`, left undefined.
        [[noreturn]] void raise_undefined_call_err_(const std::string &call_name) const;

        /// Arguments, each with the values it may take.
        using ArgRequirements = std::vector<std::pair<const Attribute *, Codomain>>;
        /// Arguments of the active parametrization, each with the values it may take.
        /// Read once, when the definition is reshaped.
        virtual ArgRequirements arg_requirements_() const { return {}; }

        /// Raw function, density or mass for a distribution.
        virtual double raw_function_(double /*x*/) const = 0;
        /// Raw log-function, log-density or log-mass for a distribution.
        virtual double raw_log_func_(double /*x*/) const { raise_undefined_call_err_("log_func"); }
        /// Raw primitive, cumulative distribution for a distribution.
        virtual double raw_primitive_(double /*x*/) const { raise_undefined_call_err_("primitive"); }
        /// Raw complementary cumulative distribution, 1 - cdf.
        ///
        /// * continuous distribution -> a tail-accurate closed form, read by Left discretization;
        /// * any other -> undefined.
        virtual double raw_ccdf_(double /*x*/) const { raise_undefined_call_err_("ccdf"); }

        /// Evaluate the function under the discretization scheme, truncation normalized.
        double function_(double x) const
        {
            double y;
            switch (discr_mode_)
            {
            case Discrete::None:
                y = raw_function_(x);
                break;
            case Discrete::Left:
                y = discr_func_left_(x);
                break;
            // Discrete::Defined:
            default:
                y = raw_function_(std::floor(x));
                break;
            }
            return trunc_mode_ == Trunc::None ? y : y / trunc_();
        }

        /// Evaluate the log-function under the discretization scheme, truncation normalized.
        double log_func_(double x) const
        {
            double y;
            switch (discr_mode_)
            {
            case Discrete::None:
                y = raw_log_func_(x);
                break;
            case Discrete::Left:
                y = std::log(discr_func_left_(x));
                break;
            // Discrete::Defined:
            default:
                y = raw_log_func_(std::floor(x));
                break;
            }
            return trunc_mode_ == Trunc::None ? y : y - logtrunc_();
        }

        /// Evaluate the primitive under the discretization scheme, truncation normalized.
        double primitive_(double x) const
        {
            const double y = raw_primitive_(primitive_x_helper_(x));
            switch (trunc_mode_)
            {
            case Trunc::None:
                return y;
            case Trunc::Max:
                return y / trunc_();
            // cases Trunc::Min, Trunc::Both
            default:
            {
                const double lo = cdf_min_();
                if (trunc_norm_mode_ == TruncNorm::LookUp)
                    return (y - lo) / lookup_trunc_;
                return (y - lo) / ((trunc_mode_ == Trunc::Min ? 1. : cdf_max_()) - lo);
            }
            }
        }

        //////  Truncation  //////

        /// Truncated sides, set by a distribution only.
        enum class Trunc : std::uint8_t
        {
            None,
            Min,
            Max,
            Both
        } trunc_mode_ = Trunc::None;
        /// Source of the normalizer.
        enum class TruncNorm : std::uint8_t
        {
            /// Constant normalizer, stored by finalize().
            LookUp,
            /// Normalizer computed per evaluation.
            Compute
        } trunc_norm_mode_ = TruncNorm::Compute;
        /// Log-normalizer read in TruncNorm::LookUp mode.
        double lookup_log_trunc_ = na::DOUBLE;
        /// Normalizer read in TruncNorm::LookUp mode.
        double lookup_trunc_ = na::DOUBLE;

        /// Primitive at the minimum, with discretization correction.
        /// @note Call only with a minimum.
        double cdf_min_() const
        {
            switch (discr_mode_)
            {
            case Discrete::Defined:
                return raw_primitive_(std::ceil(domain_.min.value()) - 1.);
            case Discrete::Left:
                return raw_primitive_(std::ceil(domain_.min.value()));
            // case Discrete::None:
            default:
                return raw_primitive_(domain_.min.value());
            }
        }

        /// Primitive at the maximum, with discretization correction.
        /// @note Call only with a maximum.
        double cdf_max_() const
        {
            switch (discr_mode_)
            {
            case Discrete::Defined:
                return raw_primitive_(std::floor(domain_.max.value()));
            case Discrete::Left:
                return raw_primitive_(std::floor(domain_.max.value()) + 1.);
            // case Discrete::None:
            default:
                return raw_primitive_(domain_.max.value());
            }
        }

        /// Normalizing constant of the truncation.
        /// @note Call only when truncated.
        double trunc_() const
        {
            if (trunc_norm_mode_ == TruncNorm::LookUp)
                return lookup_trunc_;
            switch (trunc_mode_)
            {
            case Trunc::Min:
                return 1. - cdf_min_();
            case Trunc::Max:
                return cdf_max_();
            // case Trunc::Both:
            default:
                return cdf_max_() - cdf_min_();
            }
        }

        /// Log-normalizer of the truncation.
        /// @note Call only when truncated.
        double logtrunc_() const { return trunc_norm_mode_ == TruncNorm::LookUp ? lookup_log_trunc_ : std::log(trunc_()); }

        //////  Discretization  //////

        /// Map `x` to the point the primitive reads, per discretization scheme.
        double primitive_x_helper_(double x) const
        {
            switch (discr_mode_)
            {
            case Discrete::Defined:
                return std::floor(x);
            case Discrete::Left:
                return std::floor(x) + 1.;
            default:
                return x;
            }
        }

        /// Mass of the unit bin starting at floor(`x`).
        double discr_func_left_(double x) const { return discr_bin_(std::floor(x)); }

    private:
        /// Upper-tail mass below which cumulatives lose precision.
        static constexpr double TAIL_MASS = 0x1p-20;

        /// Mass of the unit bin starting at `a`.
        ///
        /// * default -> difference of primitives;
        /// * bin above the median with mass below `TAIL_MASS` -> difference of complementary cumulatives.
        double discr_bin_(double a) const
        {
            const double hi = raw_primitive_(a + 1.), lo = raw_primitive_(a), p = hi - lo;
            return p > TAIL_MASS || lo < 0.5 ? p : raw_ccdf_(a) - raw_ccdf_(a + 1.);
        }
    };

    inline double Attribute::value() const
    {
        switch (mode_)
        {
        case Mode::Numeric:
            return value_;
        case Mode::Pointer:
            return attr_->value();
        case Mode::Function:
            return func_->evaluate();
        // case Mode::Empty
        default:
            return na::DOUBLE;
        }
    }

    //////////////////////////////  Cloneable  //////////////////////////////

    /// Implement clone() and name() of a concrete function class `D` deriving from `B`.
    ///
    /// * `D` derives from `Cloneable<D, B>` and inherits the constructors of `B`;
    /// * `D` is `final`, so no subclass clones as `D`;
    /// * `D` declares `static constexpr const char *NAME`, the name it registers under.
    template <class D, class B>
    struct Cloneable : B
    {
        using B::B;

        std::unique_ptr<Function> clone() const override final
        {
            static_assert(std::is_final_v<D>, "Class deriving from Cloneable must be final.");
            return std::make_unique<D>(static_cast<const D &>(*this));
        }
        std::string name() const noexcept override final { return D::NAME; }
    };
}
