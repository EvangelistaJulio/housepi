/* housepi - Epidemiological Analysis of Household Outbreaks.
 *
 * Copyright © 2026 Institut Pasteur (Paris), CNRS and INSERM.
 * Distributed under the terms of the GNU General Public License (GPLv3).
 * See the LICENSE file for details.
 */

#pragma once

#include <cmath>
#include <memory>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

#include "tools/nan.h"

/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////// *** Attribute : class *** ///////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

class Function;

/// Attribute parent class.
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

    //////////////////////
    //////  SETTER  //////
    //////////////////////

    /// Set attribute name.
    void name(const std::string &name) { name_ = name; }
    /// Set error type.
    void set_err_type(const std::string &err_type) { err_type_ = err_type; }

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
    /// Optimizing pruning indirections and setting constants to numeric.
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
        case Mode::Pointer:
        case Mode::Function:
            throw std::runtime_error(
                "Tried to change non numeric " + err_name() + ".");
        }
        __builtin_unreachable();
    }

    /// Initialize function attribute by assigning from another one.
    void init_func_attr(std::string_view key, const Attribute &attr);

    //////////////////////
    //////  GETTER  //////
    //////////////////////

    /// Get attribute numeric value.
    double value() const;

    /// Get pointer to attribute value.
    const Attribute *pointer() const noexcept
    {
        switch (mode_)
        {
        case Mode::Empty:
        case Mode::Numeric:
            return this;
        case Mode::Pointer:
            return pAttr_->pointer();
        case Mode::Function:
            return this;
        }
        __builtin_unreachable();
    }
    /// Get pointer to attribute function.
    /// FIXME: this should not be a constant method.
    const Function *function() const noexcept
    {
        switch (mode_)
        {
        case Mode::Empty:
        case Mode::Numeric:
            return nullptr;
        case Mode::Pointer:
            return pAttr_->function();
        case Mode::Function:
            return pFunc_.get();
        }
        __builtin_unreachable();
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
    /// Get a summary of the attribute.
    virtual std::string str_summary() const;

protected:
    enum class Mode : uint8_t
    {
        Empty,
        Numeric,
        Pointer,
        Function
    } mode_ = Mode::Empty;

    std::string err_type_ = "empty";            /// type name displayed in errors
    std::string name_ = "";                     /// name of argument
    bool is_user_const_ = false;                /// user flag for constant
    double value_ = hsp::na::DOUBLE;            /// internal value
    const Attribute *pAttr_ = nullptr;          /// pointer to attribute
    std::unique_ptr<Function> pFunc_ = nullptr; /// pointer to function
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
        enum class Discrete : uint8_t
        {
            None,
            Center,
            Left,
            Defined
        };
        enum class Interval : uint8_t
        {
            None,
            Open,
            Closed
        };
        enum class Location : uint8_t
        {
            None,
            Active
        };
        enum class Evaluate : uint8_t
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
             double domain_min = hsp::na::DOUBLE,
             double domain_max = hsp::na::DOUBLE);
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
    virtual void init_discretization() { throw std::invalid_argument(
        "Called undefined method discretization of " + err_name() + "."); }

    /// Virtual function to initialize truncation.
    /// NOTE: always to be called ONLY after init_discretization.
    virtual void init_truncation() { helper_init_standard_range_(); };

    /// Select function method for return value.
    void set_eval_call() noexcept { eval_mode_ = Mode::Evaluate::Function; }
    /// Select primitive method for return value.
    void set_eval_primitive() noexcept { eval_mode_ = Mode::Evaluate::Primitive; }
    /// Select call as a distributin method for return value.
    void set_eval_as_distr_call() noexcept { eval_mode_ = Mode::Evaluate::Probability; }
    /// Select primitive as a distributin method for return value.
    void set_eval_as_distr_primitive() noexcept { eval_mode_ = Mode::Evaluate::Cumulative; }

    /// Set argument value.
    /// UPDATE: just required by AdaptableProposal.
    void set_arg(size_t i, double x) { arguments_[i].set(x); }
    /// Set lower value.
    /// UPDATE: just required by AdaptableProposal.
    void set_lower(double x) { lower_.set(x); }
    /// Set upper value.
    /// UPDATE: just required by AdaptableProposal.
    void set_upper(double x) { upper_.set(x); }

    /// Initialize function attribute by assigning from another one.
    void init_attr(std::string_view key, const Attribute &attr);
    
    /// Optimizing pruning indirections and setting constants to numeric.
    void optimize();

    //////////////////////
    //////  GETTER  //////
    //////////////////////

    size_t n_args() const noexcept { return arguments_.size(); }
    double arg(size_t i) const { return arguments_[i].value(); }

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
        case Mode::Evaluate::Function:
            return function_(delocate_(x));
        case Mode::Evaluate::Primitive:
            return primitive_(delocate_(x));
        case Mode::Evaluate::Probability:
            return pdf(x);
        case Mode::Evaluate::Cumulative:
            return cdf(x);
        }
        __builtin_unreachable();
    }

    /// Return value as for a probability density/mass function.
    double pdf(double x) const { return is_valid(x) ? function_(delocate_(x)) : 0.; }

    /// Return value as for a cumulative distribution function.
    double cdf(double x) const { return is_below(x) ? 0 : primitive_(delocate_(is_above(x) ? upper_.value() : x)); }

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
        case Mode::Interval::None:
            return false;
        case Mode::Interval::Open:
            return delocate_(x) <= lower_.value();
        case Mode::Interval::Closed:
            return delocate_(x) < lower_.value();
        }
        __builtin_unreachable();
    }

    /// Returns TRUE if value x is above range maximum.
    bool is_above(double x) const
    {
        switch (above_mode_)
        {
        case Mode::Interval::None:
            return false;
        case Mode::Interval::Open:
            return delocate_(x) >= upper_.value();
        case Mode::Interval::Closed:
            return delocate_(x) > upper_.value();
        }
        __builtin_unreachable();
    }

    /// Returns TRUE if value x is within range.
    bool is_valid(double x) const { return !is_below(x) && !is_above(x); }

    ////////////////////
    //////  SHOW  //////
    ////////////////////

    /// Get function name.
    virtual std::string_view name() const noexcept = 0;
    /// Get argument name.
    std::string_view arg_name(unsigned short i) const noexcept { return arguments_[i].name(); }
    /// Get function name complete with type to display in errors.
    virtual std::string err_name() const { return "function '" + std::string(name()) + "'"; }

    std::string str_range() const;
    std::string str_summary(bool show_x_arg = true) const;

protected:
    std::vector<Attribute> arguments_ = {}; /// array with arguments
    Attribute x_arg_ = Attribute();         /// 'x' argument
    Attribute loc_ = Attribute();           /// 'loc' argument
    Attribute lower_ = Attribute();         /// lower bound
    Attribute upper_ = Attribute();         /// upper bound

    /// Find attribute and return. If doesn't exists, return an empty attribute.
    virtual Attribute &find_arg_(std::string_view arg_name);

    /////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////  Function : virtual methods [protected]  ///////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////

    /// Raw method to call.
    virtual double raw_function_(double) const = 0;
    /// Raw method to call log-function.
    virtual double raw_logfunc_(double) const { throw std::invalid_argument(
        "Called undefined method 'logfunc' of " + err_name() + "."); }
    /// Raw primitive method.
    virtual double raw_primitive_(double) const { throw std::invalid_argument(
        "Called undefined method 'primitive' of " + err_name() + "."); }
    /// Raw log-ratio method. NOTE: empty function return 0 by convention.
    virtual double raw_logratio_(double, double) const { return 0.; }

    /////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////  Function : discretization [protected]  ///////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////

    /// Polymorphic function method.
    virtual double function_(double x) const
    {
        switch (discr_mode_)
        {
        case Mode::Discrete::None:
            return raw_function_(x);
        case Mode::Discrete::Center:
        {
            const double i = std::round(x);
            return raw_primitive_(i + 0.5) - raw_primitive_(i - 0.5);
        }
        case Mode::Discrete::Left:
        {
            const double i = std::floor(x);
            return raw_primitive_(i + 1.) - raw_primitive_(i);
        }
        case Mode::Discrete::Defined:
            return raw_function_(std::floor(x));
        }
        __builtin_unreachable();
    }

    /// Polymorphic log-function method.
    virtual double log_func_(double x) const
    {
        switch (discr_mode_)
        {
        case Mode::Discrete::None:
            return raw_logfunc_(x);
        case Mode::Discrete::Center:
        {
            const double i = std::round(x);
            return std::log(raw_primitive_(i + 0.5) - raw_primitive_(i - 0.5));
        }
        case Mode::Discrete::Left:
        {
            const double i = std::floor(x);
            return std::log(raw_primitive_(i + 1.) - raw_primitive_(i));
        }
        case Mode::Discrete::Defined:
            return raw_logfunc_(std::floor(x));
        }
        __builtin_unreachable();
    }

    /// Polymorphic primitive method.
    virtual double primitive_(double x) const
    {
        switch (discr_mode_)
        {
        case Mode::Discrete::None:
            return raw_primitive_(x);
        case Mode::Discrete::Center:
            return raw_primitive_(std::round(x) + 0.5);
        case Mode::Discrete::Left:
            return raw_primitive_(std::floor(x) + 1.);
        case Mode::Discrete::Defined:
            return raw_primitive_(std::floor(x));
        }
        __builtin_unreachable();
    }

    /// Polymorphic log_ratio method.
    virtual double log_ratio_(double xo, double xn) const
    {
        switch (discr_mode_)
        {
        case Mode::Discrete::None:
            return raw_logratio_(xo, xn);
        case Mode::Discrete::Center:
        {
            const double io = std::round(xo), in = std::round(xn);
            if (io == in)
                return 0.;
            return std::log(raw_primitive_(in + .5) - raw_primitive_(in - .5)) -
                   std::log(raw_primitive_(io + .5) - raw_primitive_(io - .5));
        }
        case Mode::Discrete::Left:
        {
            const double io = std::floor(xo), in = std::floor(xn);
            if (io == in)
                return 0.;
            return std::log(raw_primitive_(in + 1.) - raw_primitive_(in)) -
                   std::log(raw_primitive_(io + 1.) - raw_primitive_(io));
        }
        case Mode::Discrete::Defined:
            return raw_logratio_(std::floor(xo), std::floor(xn));
        }
        __builtin_unreachable();
    }

    /////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////  Function : lcoation [protected]  //////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////

    /// Delocate according to mode.
    /// FIXME: this hasn't been tested.
    double delocate_(double x) const
    {
        switch (deloc_mode_)
        {
        case Mode::Location::None:
            return x;
        case Mode::Location::Active:
            return x - loc_.value();
        }
        __builtin_unreachable();
    }

    /// Relocate according to mode. Used for random number generation.
    double relocate_(double x) const
    {
        switch (deloc_mode_)
        {
        case Mode::Location::None:
            return x;
        case Mode::Location::Active:
            return x + loc_.value();
        }
        __builtin_unreachable();
    }

    /////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////  Function : helpers [protected]  ///////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////

    /// Assert that proposed function attribute is valid.
    void helper_assert_valid_attribute_(std::string_view key, const Attribute &attr);

    /// Set the range of validity of the value, and the corresponding validation function.
    void helper_init_standard_range_(bool is_min_open = false, bool is_max_open = false);

    /// Functional range for lower bound and no upper bound.
    bool helper_init_range_min_(double min_val = 0., bool is_min_open = false);
};

inline double Attribute::value() const
{
    switch (mode_)
    {
    case Mode::Empty:
        throw std::runtime_error(
            "Impossible to return value for " + err_name() + ".");
    case Mode::Numeric:
        return value_;
    case Mode::Pointer:
        return pAttr_->value();
    case Mode::Function:
        return pFunc_->evaluate();
    default:
        __builtin_unreachable();
    }
}

/// Compose registries to get a function from them.
std::unique_ptr<Function> compose_func_registries(std::string_view opt,
                                                  const std::string &err_type,
                                                  const std::vector<std::string> &registries,
                                                  bool is_prim = false);
