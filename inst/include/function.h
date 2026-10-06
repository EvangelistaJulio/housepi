/* housepi - Epidemiological Analysis of Household Outbreaks.
 *
 * Copyright © 2024-2026 Institut Pasteur (Paris), CNRS and INSERM.
 * Distributed under the terms of the GNU General Public License (GPLv3).
 * See the LICENSE file for details.
 */

#pragma once

#include <cstdint>
#include <cmath>
#include <memory>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

#include "tools/messages.h"
#include "tools/nan.h"

namespace housepi
{
    class Function;

    /// Index for the global name pool of the process.
    using NameIndex = std::uint32_t;

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
        double inf = -inf::DOUBLE;            /// Infimum.
        double sup = inf::DOUBLE;             /// Supremum.
        Interval below = Interval::Unbounded; /// Openness of the lower side.
        Interval above = Interval::Unbounded; /// Openness of the upper side.
        Lattice lattice = Lattice::Doubles;   /// Lattice of the values.

        /// Whether every value returned is an integer.
        bool on_integers() const noexcept { return lattice == Lattice::Integers; }

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
        Attribute(Attribute &&);
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
        /// @note Evaluation only, never run for simulations.
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
            default: // cases Mode::Pointer, Mode::Function
                raise_set_error_();
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
                return pAttr_->pointer();
            default: // cases Mode::Empty, Mode::Numeric, Mode::Function
                return this;
            }
        }
        /// Get the function, nullptr unless in function mode.
        /// @note Pointers are not followed: call pointer() first.
        const Function *function() const noexcept { return pFunc_.get(); }

        //////  Flags  //////

        /// Whether it holds no value.
        /// - function: never missing, by convention;
        /// - reduction function: may still return NA.
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
        /// Value formatted for error messages.
        std::string str_val() const;

    protected:
        /// Storage mode.
        enum class Mode : std::uint8_t
        {
            Empty,
            Numeric,
            Pointer,
            Function
        } mode_ = Mode::Empty;

        ErrType err_type_ = ErrType::Empty;         /// Type name displayed in errors.
        bool is_user_const_ = false;                /// Whether the user declared it constant.
        NameIndex name_id_ = 0;                     /// Interned name id, 0 when empty.
        double value_ = na::DOUBLE;                 /// Numeric value.
        const Attribute *pAttr_ = nullptr;          /// Pointed attribute.
        std::unique_ptr<Function> pFunc_ = nullptr; /// Owned function.

        /// Throw on setting a value into an attribute that is not numeric.
        [[noreturn]] void raise_set_error_() const;
    };

    //////////////////////////////  Domain  //////////////////////////////

    /// Values a definition may be called on.
    /// Bounds are live attributes.
    struct Domain
    {
        Attribute min = Attribute();          /// User-defined minimum.
        Attribute max = Attribute();          /// User-defined maximum.
        Attribute loc = Attribute();          /// Shift between the frames.
        Interval below = Interval::Unbounded; /// Openness of the lower side.
        Interval above = Interval::Unbounded; /// Openness of the upper side.

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
            default: // case Interval::Unbounded:
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
            default: // case Interval::Unbounded:
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

        /// Interval notation, in the absolute frame.
        std::string str() const;
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
        } eval_mode_ = Evaluate::Function;
        /// Discretization scheme.
        enum class Discrete : std::uint8_t
        {
            None,
            Center,
            Left,
            Defined
        } discr_mode_ = Discrete::None;

    public:
        //////  Constructor  //////

        using default_args = std::vector<std::pair<std::string, double>>;
        Function(const default_args &args = {},
                 double domain_min = na::DOUBLE,
                 double domain_max = na::DOUBLE);
        virtual ~Function() = default;

        //////  Operators  //////

        /// Call evaluate() at `x`.
        double operator()(double x) const { return evaluate(x); }

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
            assert_discrete_args();
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

        /// Finalize: prune indirections, collapse constants to numeric.
        virtual void finalize();

        //////  Getter  //////

        /// Values returned, composed from the operands.
        ///
        /// - computed live, never stored: finalize() may fold the operands after reshape();
        /// - defaulted, so an out-of-tree custom function still compiles.
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
            default: // cases Evaluate::Function, Evaluate::Fold
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
        /// - above the domain: primitive at the maximum, which lives in the de-located frame.
        double cdf(double x) const
        {
            return is_below(x) ? 0 : primitive_(is_above(x) ? domain_.max.value() : delocate_(x));
        }

        //////  Flags  //////

        /// Whether the return value is constant.
        /// - mutable domain: not constant.
        bool is_constant() const noexcept
        {
            return x_arg_.is_constant() && has_constant_args_();
        }

        /// Whether discretize() put the value on a unit lattice, of unspecified offset.
        ///
        /// - location off the integers: the lattice shifts with it;
        /// - a discretized definition may return values off the integers.
        /// @note For "the values are integers", ask codomain().on_integers().
        bool is_discretized() const noexcept { return discr_mode_ != Discrete::None; }
        /// Whether the return value is always strictly positive.
        bool is_positive() const noexcept { return codomain().is_positive(); }
        /// Whether no argument is empty.
        bool is_callable_explicitly() const noexcept;
        /// Whether it can be called implicitly, with `x` bound too.
        bool is_callable_implicitly() const noexcept { return !x_arg_.is_empty() && is_callable_explicitly(); }

        /// Collect names of missing arguments.
        std::vector<std::string> get_missing_argument() const;
        /// Throw if not callable explicitly.
        void assert_is_callable_explicitly() const;
        /// Throw if discretization and location do not match.
        void assert_location_matches_discrete() const;
        /// Throw if a declared argument does not hold an integer value.
        void assert_discrete_args() const;

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

        Domain domain_ = Domain(); /// Values `x` may take, and their frame.

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

        Attribute x_arg_ = Attribute(); /// Argument `x`.

        /// Arguments of one parametrization.
        struct AltArgSet
        {
            std::string name = {};            /// Parametrization name.
            std::vector<Attribute> args = {}; /// Arguments.
        };
        std::vector<AltArgSet> alts_;            /// Alternative parametrizations, 0 canonical.
        std::uint8_t alt_idx_ = 0;               /// Active parametrization.
        bool is_user_alt_ = false;               /// Whether the user picked a re-parametrization.
        ErrType err_family_ = ErrType::Function; /// Family word shown in errors.

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
        /// - standard arguments `x`, `loc`, `min`, `max`: returned, not appended.
        /// @note Call once per key: repeated names are refused R-side.
        Attribute &append_arg_(std::string_view key);

        //////  Virtual  //////

        /// Throw on a call to virtual method `call_name`, left undefined.
        [[noreturn]] void raise_undefined_call_err_(const std::string &call_name) const;

        /// Arguments that must hold integer values.
        /// Read once, when the definition is reshaped.
        virtual std::vector<const Attribute *> discrete_args_() const { return {}; }

        /// Raw function.
        virtual double raw_function_(double /*x*/) const = 0;
        /// Raw log-function.
        virtual double raw_log_func_(double /*x*/) const { raise_undefined_call_err_("log_func"); }
        /// Raw primitive.
        virtual double raw_primitive_(double /*x*/) const { raise_undefined_call_err_("primitive"); }

        /// Evaluate the function under the discretization scheme.
        virtual double function_(double x) const
        {
            switch (discr_mode_)
            {
            case Discrete::None:
                return raw_function_(x);
            case Discrete::Center:
                return discr_func_center_(x);
            case Discrete::Left:
                return discr_func_left_(x);
            default: // Discrete::Defined:
                return raw_function_(std::floor(x));
            }
        }

        /// Evaluate the log-function under the discretization scheme.
        virtual double log_func_(double x) const
        {
            switch (discr_mode_)
            {
            case Discrete::None:
                return raw_log_func_(x);
            case Discrete::Center:
                return std::log(discr_func_center_(x));
            case Discrete::Left:
                return std::log(discr_func_left_(x));
            default: // Discrete::Defined:
                return raw_log_func_(std::floor(x));
            }
        }

        /// Evaluate the primitive under the discretization scheme.
        virtual double primitive_(double x) const { return raw_primitive_(primitive_x_helper_(x)); }

        //////  Discretization  //////

        /// Map `x` to the point the primitive reads, per discretization scheme.
        double primitive_x_helper_(double x) const
        {
            switch (discr_mode_)
            {
            case Discrete::Defined:
                return std::floor(x);
            case Discrete::Center:
                return std::round(x) + 0.5;
            case Discrete::Left:
                return std::floor(x) + 1.;
            default:
                return x;
            }
        }

        /// Mass of the unit bin centered on `x`.
        double discr_func_center_(double x) const
        {
            const double i = std::round(x);
            return raw_primitive_(i + 0.5) - raw_primitive_(i - 0.5);
        }
        /// Mass of the unit bin starting at floor(`x`).
        double discr_func_left_(double x) const
        {
            const double i = std::floor(x);
            return raw_primitive_(i + 1.) - raw_primitive_(i);
        }
    };

    inline double Attribute::value() const
    {
        switch (mode_)
        {
        case Mode::Numeric:
            return value_;
        case Mode::Pointer:
            return pAttr_->value();
        case Mode::Function:
            return pFunc_->evaluate();
        default: // case Mode::Empty
            return na::DOUBLE;
        }
    }

}
