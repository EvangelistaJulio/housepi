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

#include "tools/nan.h"

namespace housepi
{
    class Function;

    /////////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////// *** Attribute : class *** ///////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////

    class Attribute
    {
    public:
        ///////////////////////////
        //////  CONSTRUCTOR  //////
        ///////////////////////////

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

        /////////////////////////
        //////  OPERATORS  //////
        /////////////////////////

        /// Equal to operator.
        bool operator==(const Attribute &other) const noexcept;
        /// Not equal to operator.
        bool operator!=(const Attribute &other) const noexcept { return !(*this == other); }

        //////////////////////
        //////  SETTER  //////
        //////////////////////

        /// Set attribute name.
        void name(const std::string &name) { name_ = name; }
        void name(const std::string &name, const std::string &err_type) {
            name_ = name;
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
        /// Optimizing by pruning indirections and setting constants to numeric.
        /// NOTE: evaluation-only, never run for simulations.
        virtual void optimize();

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

        //////////////////////
        //////  GETTER  //////
        //////////////////////

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

        /////////////////////
        //////  FLAGS  //////
        /////////////////////

        /// Wheter it's constant.
        bool is_constant() const noexcept;
        /// Wheter it's discrete.
        virtual bool is_discrete() const noexcept;

        /// Wheter it has no definitions.
        bool is_empty() const noexcept { return mode_ == Mode::Empty; }
        /// Wheter it's a number.
        bool is_numeric() const noexcept { return mode_ == Mode::Numeric; }
        /// Wheter it's a pointer.
        bool is_pointer() const noexcept { return mode_ == Mode::Pointer; }
        /// Wheter it's a function.
        bool is_function() const noexcept { return mode_ == Mode::Function; }

        /// Wheter it's able to return.
        bool is_callable() const noexcept;

        ////////////////////
        //////  SHOW  //////
        ////////////////////

        /// Get attribute name.
        std::string_view name() const noexcept { return name_; }
        /// Get attribute name to display in errors.
        std::string err_name() const { return err_type_ + " '" + std::string(name()) + "'"; }
        /// Get attribute value to display in errors.
        std::string str_val() const;
#ifndef NDEBUG
        /// Get a summary of the attribute.
        virtual std::string str_summary() const;
#endif

    protected:
        enum class Mode : std::uint8_t
        {
            Empty,
            Numeric,
            Pointer,
            Function
        } mode_ = Mode::Empty;

        std::string err_type_ = "empty";            /// type name displayed in errors
        std::string name_ = "";                     /// name of argument
        bool is_user_const_ = false;                /// user flag for constant
        double value_ = na::DOUBLE;                 /// internal value
        const Attribute *pAttr_ = nullptr;          /// pointer to attribute
        std::unique_ptr<Function> pFunc_ = nullptr; /// pointer to function

        [[noreturn]] void raise_set_error_() const;
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
        struct Mode
        {
            enum class Discrete : std::uint8_t
            {
                None,
                Center,
                Left,
                Defined
            };
            enum class Interval : std::uint8_t
            {
                None,
                Open,
                Closed
            };
            enum class Location : std::uint8_t
            {
                None,
                Active
            };
            enum class Evaluate : std::uint8_t
            {
                Function,
                Primitive,
                Probability,
                Cumulative
            };
        };
        Mode::Evaluate eval_mode_ = Mode::Evaluate::Function;
        Mode::Discrete discr_mode_ = Mode::Discrete::None;
        Mode::Location deloc_mode_ = Mode::Location::None;
        Mode::Interval below_mode_ = Mode::Interval::None;
        Mode::Interval above_mode_ = Mode::Interval::None;

    public:
        ///////////////////////////
        //////  CONSTRUCTOR  //////
        ///////////////////////////

        Function(const std::vector<std::pair<std::string, double>> &args = {},
                 double domain_min = na::DOUBLE,
                 double domain_max = na::DOUBLE);
        virtual ~Function() = default;

        /////////////////////////
        //////  OPERATORS  //////
        /////////////////////////

        /// Sugar operator to call evaluate.
        double operator()(double x) const { return evaluate(x); }

        virtual std::unique_ptr<Function> clone() const = 0;

        //////////////////////
        //////  SETTER  //////
        //////////////////////

        /// Virtual function to initialize discretization.
        virtual void init_discretization() { raise_undefined_call_err_("discretization"); }

        /// Virtual function to initialize truncation.
        /// NOTE: always to be called ONLY after init_discretization.
        virtual void init_truncation() { helper_init_standard_range_(); };

        /// Select primitive method for implicit evaluate.
        void set_eval_as_prim() noexcept { eval_mode_ = Mode::Evaluate::Primitive; }
        /// Select pdf as implicit evaluate (default for distributions).
        void set_eval_as_pdf() noexcept { eval_mode_ = Mode::Evaluate::Probability; }
        /// Select cdf as implicit evaluate.
        void set_eval_as_cdf() noexcept { eval_mode_ = Mode::Evaluate::Cumulative; }

        /// Set argument value.
        /// UPDATE: just required by AdaptableProposal.
        void set_arg(size_t i, double x) { alts_[alt_idx_].args[i].set(x); }
        /// Set lower value.
        /// UPDATE: just required by AdaptableProposal.
        void set_lower(double x) { lower_.set(x); }
        /// Set upper value.
        /// UPDATE: just required by AdaptableProposal.
        void set_upper(double x) { upper_.set(x); }

        /// Initialize function attribute by assigning from another one.
        void init_attr(std::string_view key, const Attribute &attr);

        /// Optimize by pruning indirections and setting constants to numeric.
        virtual void optimize();

        //////////////////////
        //////  GETTER  //////
        //////////////////////

        /// Number of arguments in the *active* parametrization.
        size_t n_args() const noexcept { return arguments_().size(); }
        /// Read the i-th argument value of the *active* parametrization.
        double arg(size_t i) const { return alts_[alt_idx_].args[i].value(); }

        const Attribute &loc() const noexcept { return loc_; }
        const Attribute &lower() const noexcept { return lower_; }
        const Attribute &upper() const noexcept { return upper_; }

        /// Evaluate implicitly at 'x' value.
        double evaluate() const { return evaluate(x_arg_.value()); }

        /// Central method to return function's value at given x.
        double evaluate(double x) const
        {
            switch (eval_mode_)
            {
            case Mode::Evaluate::Primitive:
                return primitive_(delocate_(x));
            case Mode::Evaluate::Probability:
                return pdf(x);
            case Mode::Evaluate::Cumulative:
                return cdf(x);
            default: // case  Mode::Evaluate::Function:
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
        double pdf(double x) const { return is_valid(x) ? function_(delocate_(x)) : 0.; }
        /// Return value as for a cumulative distribution function.
        double cdf(double x) const { return is_below(x) ? 0 : primitive_(delocate_(is_above(x) ? upper_.value() : x)); }

        /// Return value of log-ratio of probability density/mass functions.
        double log_ratio(double xo, double xn) const { return log_ratio_(delocate_(xo), delocate_(xn)); }
        /// Check for validity, then return value of log-ratio of probability density/mass functions.
        double guarded_log_ratio(double xo, double xn) const
        {
            bool ov = is_valid(xo), nv = is_valid(xn);
            if (ov && nv)
                return log_ratio_(delocate_(xo), delocate_(xn));
            else if (!ov && nv)
                return std::numeric_limits<double>::infinity();
            else if (ov && !nv)
                return -std::numeric_limits<double>::infinity();
            return na::DOUBLE;
        }

        /////////////////////
        //////  FLAGS  //////
        /////////////////////

        /// Wheter return value is always constant.
        bool is_constant() const noexcept;
        /// Wheter return value is always discrete.
        virtual bool is_discrete() const noexcept { return discr_mode_ != Mode::Discrete::None; }
        /// Wheter all arguments are not empty.
        bool is_callable_explicitly() const noexcept;
        /// Wheter it can be called implicitely, i. e. if also 'x' has been linked.
        bool is_callable_implicitly() const noexcept { return !x_arg_.is_empty() && is_callable_explicitly(); }

        /// Collect names of missing arguments.
        std::vector<std::string> get_missing_argument() const;
        /// Throw exception if not callable explicitly.
        void assert_is_callable_explicitly() const;

        /// Returns TRUE if value x is below range minimum.
        bool is_below(double x) const
        {
            switch (below_mode_)
            {
            case Mode::Interval::Open:
                return delocate_(x) <= lower_.value();
            case Mode::Interval::Closed:
                return delocate_(x) < lower_.value();
            default: // case Mode::Interval::None:
                return false;
            }
        }

        /// Returns TRUE if value x is above range maximum.
        bool is_above(double x) const
        {
            switch (above_mode_)
            {
            case Mode::Interval::Open:
                return delocate_(x) >= upper_.value();
            case Mode::Interval::Closed:
                return delocate_(x) > upper_.value();
            default: // case Mode::Interval::None:
                return false;
            }
        }

        /// Returns TRUE if value x is within range.
        bool is_valid(double x) const { return !is_below(x) && !is_above(x); }

        ////////////////////
        //////  SHOW  //////
        ////////////////////

        /// Get function name.
        virtual std::string_view name() const noexcept = 0;
        /// Get argument name (from the active parametrization).
        std::string_view arg_name(unsigned short i) const noexcept { return alts_[alt_idx_].args[i].name(); }
        /// Get function name complete with type to display in errors.
        virtual std::string err_name() const { return "function '" + std::string(name()) + "'"; }

        /// Get a summary of the range.
        std::string str_range() const;
#ifndef NDEBUG
        /// Get a summary of the function.
        std::string str_summary(bool show_x_arg = true) const;
#endif

    protected:
        /////////////////////////
        //////  ARGUMENTS  //////
        /////////////////////////

        Attribute x_arg_ = Attribute(); /// 'x' argument
        Attribute loc_ = Attribute();   /// 'loc' argument
        Attribute lower_ = Attribute(); /// lower bound
        Attribute upper_ = Attribute(); /// upper bound

        /// Alternative arguments for re-parametrization
        struct AltArgSet
        {
            std::string name = {};            /// Alternative set name
            std::vector<Attribute> args = {}; /// Alternative arguments
        };
        std::vector<AltArgSet> alts_; /// Alternative parametrizations (canonical is 0)
        unsigned short alt_idx_ = 0;  /// Current select arguments
        bool is_user_alt_ = false;    /// Whether user has picked a re-parametrization

        /// Register an alternative parametrization.
        void register_alt_(const std::string &name,
                           const std::vector<std::pair<std::string, double>> &args);

        std::vector<Attribute> &arguments_() noexcept { return alts_[alt_idx_].args; }
        const std::vector<Attribute> &arguments_() const noexcept { return alts_[alt_idx_].args; }
        Attribute &arguments_(size_t i) noexcept { return alts_[alt_idx_].args[i]; }
        const Attribute &arguments_(size_t i) const noexcept { return alts_[alt_idx_].args[i]; }

        /// Find attribute and return its reference.
        /// If doesn't exists, return an empty attribute.
        virtual Attribute &find_arg_(std::string_view arg_name);

        ///////////////////////
        //////  VIRTUAL  //////
        ///////////////////////

        /// Helper to raise call if undefined virtual methods are called.
        [[noreturn]] void raise_undefined_call_err_(const std::string &call_name) const;

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
            case Mode::Discrete::None:
                return raw_function_(x);
            case Mode::Discrete::Center:
                return discr_func_center_(x);
            case Mode::Discrete::Left:
                return discr_func_left_(x);
            default: // Mode::Discrete::Defined:
                return raw_function_(std::floor(x));
            }
        }

        /// Polymorphic log-function method.
        virtual double log_func_(double x) const
        {
            switch (discr_mode_)
            {
            case Mode::Discrete::None:
                return raw_log_func_(x);
            case Mode::Discrete::Center:
                return std::log(discr_func_center_(x));
            case Mode::Discrete::Left:
                return std::log(discr_func_left_(x));
            default: // Mode::Discrete::Defined:
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
            case Mode::Discrete::None:
                return raw_log_ratio_(xo, xn);
            case Mode::Discrete::Center:
                return std::log(discr_func_center_(xn)) - std::log(discr_func_center_(xo));
            case Mode::Discrete::Left:
                return std::log(discr_func_left_(xn)) - std::log(discr_func_left_(xo));
            default: // Mode::Discrete::Defined:
                return raw_log_ratio_(std::floor(xo), std::floor(xn));
            }
        }

        //////////////////////////////
        //////  DISCRETIZATION  //////
        //////////////////////////////

        /// Helper for the primitive method.
        double primitive_x_helper_(double x) const
        {
            switch (discr_mode_)
            {
            case Mode::Discrete::Defined:
                return std::floor(x);
            case Mode::Discrete::Center:
                return std::round(x) + 0.5;
            case Mode::Discrete::Left:
                return std::floor(x) + 1.;
            default:
                return x;
            }
        }

        ////////////////////////
        //////  LOCATION  //////
        ////////////////////////

        /// Delocate according to mode.
        /// FIXME: this hasn't been tested.
        double delocate_(double x) const
        {
            switch (deloc_mode_)
            {
            case Mode::Location::Active:
                return x - loc_.value();
            default: // case Mode::Location::None:
                return x;
            }
        }

        /// Relocate according to mode. Used for random number generation.
        double relocate_(double x) const
        {
            switch (deloc_mode_)
            {
            case Mode::Location::Active:
                return x + loc_.value();
            default: // case Mode::Location::None:
                return x;
            }
        }

        ///////////////////////
        //////  HELPERS  //////
        ///////////////////////

        /// Set the range of validity of the value, and the corresponding validation function.
        void helper_init_standard_range_(bool is_min_open = false, bool is_max_open = false);
        /// Functional range for lower bound and no upper bound.
        bool helper_init_range_min_(double min_val = 0., bool is_min_open = false);

    private:
        /// Assert that proposed function attribute is valid.
        void helper_assert_valid_attribute_(std::string_view key, const Attribute &attr);

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

    /// Compose registries to get a function from them.
    std::unique_ptr<Function> compose_func_registries(
        std::string_view opt,
        const std::string &err_type,
        const std::vector<std::string> &registries,
        bool is_prim = false);
}
