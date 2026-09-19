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
        double inf = -inf::DOUBLE;
        double sup = inf::DOUBLE;
        Interval below = Interval::Unbounded;
        Interval above = Interval::Unbounded;
        Lattice lattice = Lattice::Doubles;

        /// Wheter every value returned is an integer.
        bool on_integers() const noexcept { return lattice == Lattice::Integers; }

        /// Interval notation, in the absolute frame.
        std::string str() const;

        /// Wheter every value returned is provably greater than zero.
        /// NOTE: one-sided, "false" does not mean "negative".
        bool is_positive() const noexcept
        {
            return below != Interval::Unbounded &&
                   (inf > 0. || (inf == 0. && below == Interval::Open));
        }
    };

    /////////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////// *** Attribute : class *** ///////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////

    class Attribute
    {
    public:
        /////////////////////
        ///  CONSTRUCTOR  ///
        /////////////////////

        Attribute();
        Attribute(double x) { init_numeric(x); }
        virtual ~Attribute();

        /// Copy constructor
        Attribute(const Attribute &other);
        /// Move constructor
        Attribute(Attribute &&);
        /// Copy assignment
        Attribute &operator=(Attribute other);

        void swap(Attribute &other) noexcept;

        ///////////////////
        ///  OPERATORS  ///
        ///////////////////

        /// Equal to operator.
        /// UPDATE: are these still necessary?
        bool operator==(const Attribute &other) const noexcept;
        /// Not equal to operator.
        bool operator!=(const Attribute &other) const noexcept { return !(*this == other); }

        ////////////////
        ///  SETTER  ///
        ////////////////

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
        /// Productory of pointed attributes values.
        void init_product(const std::vector<const Attribute *> &other_ptrs);
        /// Transfer from other Attribute.
        void transfer(const Attribute &other);
        /// Finalizing by pruning indirections and setting constants to numeric.
        /// NOTE: evaluation-only, never run for simulations.
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

        /// Initialize function attribute by assigning from another one.
        void init_func_attr(std::string_view key, const Attribute &attr);

        ////////////////
        ///  GETTER  ///
        ////////////////

        /// Values the attribute takes.
        virtual Codomain codomain() const noexcept;

        /// Get attribute numeric value.
        // UPDATE: Mode::Pointer forces value out-of-line.
        double value() const;

        /// Get pointer to attribute value.
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
        /// Get pointer to attribute function.
        /// FIXME: this should not be a constant method.
        const Function *function() const noexcept
        {
            switch (mode_)
            {
            case Mode::Pointer:
                return pAttr_->function();
            case Mode::Function:
                return pFunc_.get();
            default: // cases Mode::Empty, Mode::Numeric
                return nullptr;
            }
        }

        ///////////////
        ///  FLAGS  ///
        ///////////////

        /// Wheter it holds no value.
        /// NOTE: by convention function is never missing.
        /// NOTE: reduction function can still return NA.
        bool is_na() const noexcept;
        /// Wheter it's constant.
        bool is_constant() const noexcept;
        /// Wheter every value it returns is a discrete number.
        bool is_discrete() const noexcept { return codomain().on_integers(); }
        /// Wheter every value it returns is a positive number.
        bool is_positive() const noexcept { return codomain().is_positive(); }
        /// Wheter it's able to return.
        bool is_callable() const noexcept;

        /// Wheter it has no definitions.
        bool is_empty() const noexcept { return mode_ == Mode::Empty; }
        /// Wheter it's a number.
        bool is_numeric() const noexcept { return mode_ == Mode::Numeric; }
        /// Wheter it's a pointer.
        bool is_pointer() const noexcept { return mode_ == Mode::Pointer; }
        /// Wheter it's a function.
        bool is_function() const noexcept { return mode_ == Mode::Function; }

        //////////////
        ///  SHOW  ///
        //////////////

        /// Get attribute name.
        std::string name() const;
        /// Name formatted for error messages.
        std::string str() const;
        /// Value formatted for error messages.
        std::string str_val() const;

    protected:
        enum class Mode : std::uint8_t
        {
            Empty,
            Numeric,
            Pointer,
            Function
        } mode_ = Mode::Empty;

        ErrType err_type_ = ErrType::Empty;         /// type name displayed in errors
        bool is_user_const_ = false;                /// user flag for constant
        NameIndex name_id_ = 0;                     /// interned name id (0 == empty)
        double value_ = na::DOUBLE;                 /// internal value
        const Attribute *pAttr_ = nullptr;          /// pointer to attribute
        std::unique_ptr<Function> pFunc_ = nullptr; /// pointer to function

        [[noreturn]] void raise_set_error_() const;
    };

    /////////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////// *** Domain : struct *** /////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////

    /// Values a definition may be called on.
    /// NOTE: bounds are live attributes.
    struct Domain
    {
        Attribute min = Attribute(); /// user defined minimum
        Attribute max = Attribute(); /// user defined maximum
        Attribute loc = Attribute(); /// shift between the frames
        Interval below = Interval::Unbounded;
        Interval above = Interval::Unbounded;

        /// Into the frame the bounds and all raw math live in.
        /// NOTE: read per evaluation, so the guard stays a mode compare.
        double delocate(double x) const { return loc.is_empty() ? x : x - loc.value(); }
        /// Back into the frame the caller sees.
        double relocate(double x) const { return loc.is_empty() ? x : x + loc.value(); }

        /// Wheter x falls below the minimum.
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

        /// Wheter x falls above the maximum.
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

        /// Wheter x falls inside.
        bool is_valid(double x) const { return !is_below(x) && !is_above(x); }
        /// The interval the bounds delimit, in the absolute frame.
        /// NOTE: composed from the codomains of the bounds, so an undrawn
        /// bound still places the interval it can reach.
        Codomain relocated() const;
        /// Snapshot into the absolute frame, for the support of a distribution.
        /// NOTE: a dense support puts no mass on its floor, so a closed bound
        /// reads open there - that is where boundary attainability is decided.
        Codomain relocated(Lattice lattice) const;

        /// Set the openness of both sides from the bounds provided.
        /// @param owner : name of the definition, for the error it may raise.
        void init_standard(bool is_min_open, bool is_max_open, const std::string &owner);
        /// Lower the minimum onto 'min_val'. Returns whether truncation is needed.
        bool init_min(double min_val, bool is_min_open, const std::string &owner);

        /// UPDATE: consistency check now are anaware of codomains/constantness of attributes.
        /// Wheter both bounds read and delimit a non-empty interval.
        bool is_consistent() const;
        /// Wheter both bounds read and coincide on a single point.
        bool is_degenerate() const { return min.is_callable() && max.is_callable() && min.value() == max.value(); }

        /// Raise a uniform error if the domain is not consistent.
        /// @param owner : name of the definition, for the error it may raise.
        void assert_consistent(const std::string &owner) const;

        /// Interval notation, in the absolute frame.
        std::string str() const;
    };

    /////////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////// *** Function : class *** /////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////

    /// Function parent class. Arguments stored in arrays.
    class Function
    {
    protected:
        enum class Evaluate : std::uint8_t
        {
            Function,
            Primitive,
            Probability,
            Cumulative
        } eval_mode_ = Evaluate::Function;
        enum class Discrete : std::uint8_t
        {
            None,
            Center,
            Left,
            Defined
        } discr_mode_ = Discrete::None;

    public:
        /////////////////////
        ///  CONSTRUCTOR  ///
        /////////////////////

        Function(const std::vector<std::pair<std::string, double>> &args = {},
                 double domain_min = na::DOUBLE,
                 double domain_max = na::DOUBLE);
        virtual ~Function() = default;

        ///////////////////
        ///  OPERATORS  ///
        ///////////////////

        /// Sugar operator to call evaluate.
        double operator()(double x) const { return evaluate(x); }

        virtual std::unique_ptr<Function> clone() const = 0;

        ////////////////
        ///  SETTER  ///
        ////////////////

        /// UPDATE: these two methods can turn into protected probably.
        /// Virtual function to initialize discretization.
        virtual void discretize() { raise_undefined_call_err_("discretize"); }
        /// Virtual function to initialize truncation.
        /// NOTE: always to be called ONLY after discretize.
        virtual void truncate() { helper_init_standard_domain_(); };

        /// Finalize discretization, truncation and consistency asserts.
        /// NOTE: to be called once all arguments and domain are assigned.
        void reshape(bool discr)
        {
            if (discr)
                discretize();
            truncate();

            assert_is_callable_explicitly();
            assert_location_matches_discrete();
            assert_discrete_args();
        }

        /// Select primitive method for implicit evaluate.
        void set_eval_as_prim() noexcept { eval_mode_ = Evaluate::Primitive; }
        /// Select pdf as implicit evaluate (default for distributions).
        void set_eval_as_pdf() noexcept { eval_mode_ = Evaluate::Probability; }
        /// Select cdf as implicit evaluate.
        void set_eval_as_cdf() noexcept { eval_mode_ = Evaluate::Cumulative; }

        /// Set argument value.
        /// UPDATE: just required by AdaptableProposal.
        void set_arg(size_t i, double x) { alts_[alt_idx_].args[i].set(x); }
        /// Set the minimum.
        /// UPDATE: just required by AdaptableProposal.
        void set_min(double x) { domain_.min.set(x); }
        /// Set the maximum.
        /// UPDATE: just required by AdaptableProposal.
        void set_max(double x) { domain_.max.set(x); }

        /// Transfer a live bound attribute into the min/max truncation, preserving
        /// pointer/function liveness so the domain tracks the source per evaluation.
        void transfer_min(const Attribute &a) { domain_.min.transfer(a); }
        void transfer_max(const Attribute &a) { domain_.max.transfer(a); }

        /// Initialize function attribute by assigning from another one.
        void init_attr(std::string_view key, const Attribute &attr);

        /// Finalize by pruning indirections and setting constants to numeric.
        virtual void finalize();

        ////////////////
        ///  GETTER  ///
        ////////////////

        /// Values returned, composed from the operands.
        /// NOTE: computed live, never stored - the algebraic nodes read their
        /// operands, which finalize() may fold after reshape() has run.
        /// NOTE: defaulted, so an out-of-tree custom function still compiles.
        virtual Codomain codomain() const noexcept { return Codomain(); }

        /// Number of arguments in the *active* parametrization.
        size_t n_args() const noexcept { return arguments_().size(); }
        /// Read the i-th argument value of the *active* parametrization.
        double arg(size_t i) const { return alts_[alt_idx_].args[i].value(); }

        const Attribute &loc() const noexcept { return domain_.loc; }
        const Attribute &min() const noexcept { return domain_.min; }
        const Attribute &max() const noexcept { return domain_.max; }
        /// Values 'x' may take, in the de-located frame.
        const Domain &domain() const noexcept { return domain_; }

        /// Evaluate implicitly at 'x' value.
        double evaluate() const { return evaluate(x_arg_.value()); }

        /// Central method to return function's value at given x.
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
            default: // case  Evaluate::Function:
                return function_(delocate_(x));
            }
        }

        /// Return function value at x (with delocation applied).
        double func(double x) const { return function_(delocate_(x)); }
        /// Return log-function value at x (with delocation applied).
        double log_func(double x) const { return log_func_(delocate_(x)); }
        /// Return primitive value at x (with delocation applied).
        double prim(double x) const { return primitive_(delocate_(x)); }

        /// Return value as for a probability density/mass function.
        double pdf(double x) const
        {
            return is_valid(x) ? function_(delocate_(x)) : 0.;
        }
        /// Return value as for a cumulative distribution function.
        /// NOTE: the maximum already lives in the de-located frame.
        double cdf(double x) const
        {
            return is_below(x) ? 0 : primitive_(is_above(x) ? domain_.max.value() : delocate_(x));
        }

        /// Return value of log-ratio of probability density/mass functions.
        /// NOTE: a null move carries no density change, whatever the density at that point.
        double log_ratio(double xo, double xn) const
        {
            return xo == xn ? 0. : log_ratio_(delocate_(xo), delocate_(xn));
        }
        /// Check for validity, then return value of log-ratio of probability density/mass functions.
        double guarded_log_ratio(double xo, double xn) const
        {
            bool ov = is_valid(xo), nv = is_valid(xn);
            if (ov && nv)
                return log_ratio(xo, xn);
            else if (!ov && nv)
                return inf::DOUBLE;
            else if (ov && !nv)
                return -inf::DOUBLE;
            return na::DOUBLE;
        }

        ///////////////
        ///  FLAGS  ///
        ///////////////

        /// Wheter return value is always constant.
        /// NOTE: at the moment also a mutable domain makes the function non constant.
        bool is_constant() const noexcept
        {
            return x_arg_.is_constant() && has_constant_args_();
        }

        /// Wheter discretize() has put the value on a unit lattice, of unspecified
        /// offset. This is the machinery: a location off the integers carries the
        /// lattice with it, so a discretized definition need not return integer values.
        /// NOTE: for "the values are integers" ask codomain().on_integers() instead.
        bool is_discretized() const noexcept { return discr_mode_ != Discrete::None; }
        /// Wheter return value is always strictly positive.
        bool is_positive() const noexcept { return codomain().is_positive(); }
        /// Wheter all arguments are not empty.
        bool is_callable_explicitly() const noexcept;
        /// Wheter it can be called implicitely, i. e. if also 'x' has been linked.
        bool is_callable_implicitly() const noexcept { return !x_arg_.is_empty() && is_callable_explicitly(); }

        /// Collect names of missing arguments.
        std::vector<std::string> get_missing_argument() const;
        /// Throw exception if not callable explicitly.
        void assert_is_callable_explicitly() const;
        /// Throw exception if discrete and location do not match.
        void assert_location_matches_discrete() const;
        /// Throw exception if a declared argument does not hold an integer value.
        void assert_discrete_args() const;

        /// Returns TRUE if value x is below the domain minimum.
        bool is_below(double x) const { return domain_.is_below(x); }
        /// Returns TRUE if value x is above the domain maximum.
        bool is_above(double x) const { return domain_.is_above(x); }
        /// Returns TRUE if value x is inside the domain.
        bool is_valid(double x) const { return domain_.is_valid(x); }

        //////////////
        ///  SHOW  ///
        //////////////

        /// Get function name.
        virtual std::string name() const noexcept = 0;
        /// Get argument name (from the active parametrization).
        std::string arg_name(unsigned short i) const { return alts_[alt_idx_].args[i].name(); }
        /// Name formatted for error messages.
        virtual std::string str() const;
        /// Set the family word shown in errors (see str()).
        void set_err_family(ErrType family) noexcept { err_family_ = family; }

    protected:
        ////////////////
        ///  DOMAIN  ///
        ////////////////

        Domain domain_ = Domain(); /// values 'x' may take, and the frame they live in

        /// Cast into the de-located frame.
        double delocate_(double x) const { return domain_.delocate(x); }
        /// Cast back into the absolute frame.
        double relocate_(double x) const { return domain_.relocate(x); }

        /// Whether the active arguments and the bounded sides of the domain are constant.
        bool has_constant_args_() const noexcept;

        /// Set the domain of validity of the value, and the corresponding validation function.
        /// UPDATE: passing str() forces to compute it at runtime...
        void helper_init_standard_domain_(bool is_min_open = false, bool is_max_open = false)
        {
            domain_.init_standard(is_min_open, is_max_open, str());
        }
        /// Set a domain for given minimum and no maximum.
        bool helper_init_domain_min_(double min_val = 0., bool is_min_open = false)
        {
            return domain_.init_min(min_val, is_min_open, str());
        }

        ///////////////////
        ///  ARGUMENTS  ///
        ///////////////////

        Attribute x_arg_ = Attribute(); /// 'x' argument

        /// Alternative arguments for re-parametrization
        struct AltArgSet
        {
            std::string name = {};            /// Alternative set name
            std::vector<Attribute> args = {}; /// Alternative arguments
        };
        std::vector<AltArgSet> alts_;            /// Alternative parametrizations (canonical is 0)
        std::uint8_t alt_idx_ = 0;               /// Current select arguments
        bool is_user_alt_ = false;               /// Whether user has picked a re-parametrization
        ErrType err_family_ = ErrType::Function; /// family word shown in errors

        /// Register an alternative parametrization.
        void register_alt_(const std::string &name,
                           const std::vector<std::pair<std::string, double>> &args);

        /// Assert that proposed function attribute is valid.
        void helper_assert_valid_attribute_(std::string_view key, const Attribute &attr);

        std::vector<Attribute> &arguments_() noexcept { return alts_[alt_idx_].args; }
        const std::vector<Attribute> &arguments_() const noexcept { return alts_[alt_idx_].args; }
        Attribute &arguments_(size_t i) noexcept { return alts_[alt_idx_].args[i]; }
        const Attribute &arguments_(size_t i) const noexcept { return alts_[alt_idx_].args[i]; }

        /// Find attribute and return its reference.
        /// If doesn't exists, return an empty attribute.
        virtual Attribute &find_arg_(std::string_view arg_name);

        /// Find-or-append a named argument, for variadic parametrizations.
        /// If attribute is absent at initialization, create it.
        Attribute &find_or_append_arg_(std::string_view key);

        /////////////////
        ///  VIRTUAL  ///
        /////////////////

        /// Helper to raise call if undefined virtual methods are called.
        [[noreturn]] void raise_undefined_call_err_(const std::string &call_name) const;

        /// Attributes that must return integer values.
        /// NOTE: read once, when the definition is reshaped.
        virtual std::vector<const Attribute *> discrete_args_() const { return {}; }

        /// Raw method to call.
        virtual double raw_function_(double /*x*/) const = 0;
        /// Raw method to call log-function.
        virtual double raw_log_func_(double /*x*/) const { raise_undefined_call_err_("log_func"); }
        /// Raw primitive method.
        virtual double raw_primitive_(double /*x*/) const { raise_undefined_call_err_("primitive"); }
        /// Raw quantile method.
        virtual double raw_quantile_(double /*p*/) const { raise_undefined_call_err_("quantile"); }
        /// Raw log-ratio method. NOTE: empty function return 0 by convention.
        virtual double raw_log_ratio_(double /*xo*/, double /*xn*/) const { return 0.; }

        /// Polymorphic function method.
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

        /// Polymorphic log-function method.
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

        /// Polymorphic primitive method.
        virtual double primitive_(double x) const { return raw_primitive_(primitive_x_helper_(x)); }

        /// Polymorphic log_ratio method.
        virtual double log_ratio_(double xo, double xn) const
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

        ////////////////////////
        ///  DISCRETIZATION  ///
        ////////////////////////

        /// Helper for the primitive method.
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

        /// Helper for Center discretization.
        double discr_func_center_(double x) const
        {
            const double i = std::round(x);
            return raw_primitive_(i + 0.5) - raw_primitive_(i - 0.5);
        }
        /// Helper for Left discretization.
        double discr_func_left_(double x) const
        {
            const double i = std::floor(x);
            return raw_primitive_(i + 1.) - raw_primitive_(i);
        }
    };

}
