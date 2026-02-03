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

#include "attribute.h"

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
        enum class Discr : uint8_t
        {
            None,   /// No discretization
            Center, /// Centered discretization
            Left,   /// Left discretization
            Defined /// Discrete defined
        };
        enum class Interv : uint8_t
        {
            None,  /// No interval
            Open,  /// Open interval
            Closed /// Closed interval
        };
        enum class Deloc : uint8_t
        {
            None,  /// No delocalization
            Active /// Active delocalization
        };
        enum class Eval : uint8_t
        {
            Func,  /// as function
            Prim,  /// as primitive
            Proba, /// as probability mass/density
            Cumul  /// as cumulative
        };
    };
    Mode::Eval eval_mode_ = Mode::Eval::Func;      /// evaluation mode
    Mode::Discr discr_mode_ = Mode::Discr::None;   /// discretization mode
    Mode::Deloc deloc_mode_ = Mode::Deloc::None;   /// delocation mode
    Mode::Interv below_mode_ = Mode::Interv::None; /// lower interval below mode
    Mode::Interv above_mode_ = Mode::Interv::None; /// upper interval above mode

public:
    ///////////////////////////
    //////  CONSTRUCTOR  //////
    ///////////////////////////

    Function(const std::vector<std::pair<std::string, double>> &args = {},
             double domain_min = hsp::na::DOUBLE, bool is_min_open = false,
             double domain_max = hsp::na::DOUBLE, bool is_max_open = false);
    virtual ~Function() = default;

    /////////////////////////
    //////  OPERATORS  //////
    /////////////////////////

    /// Sugar operator to call evaluate.
    double operator()(double x) const { return evaluate(x); }

    //////////////////////
    //////  SETTER  //////
    //////////////////////

    /// Initialize given function attribute.
    void init_attr(std::string_view key, const Attribute &attr);

    /// Virtual function to initialize discretization.
    virtual void init_discretization() { throw std::invalid_argument(
        "Called undefined method discretization of " + err_name() + "."); }

    /// Virtual function to initialize truncation.
    /// NOTE: always to be call after init_discretization.
    virtual void init_truncation() { throw std::invalid_argument(
        "Called to undefined method truncation of " + err_name() + "."); }

    /// Select function method for return value.
    void set_eval_call() noexcept { eval_mode_ = Mode::Eval::Func; }
    /// Select primitive method for return value.
    void set_eval_primitive() noexcept { eval_mode_ = Mode::Eval::Prim; }
    /// Select call as a distributin method for return value.
    void set_eval_as_distr_call() noexcept { eval_mode_ = Mode::Eval::Proba; }
    /// Select primitive as a distributin method for return value.
    void set_eval_as_distr_primitive() noexcept { eval_mode_ = Mode::Eval::Cumul; }

    /// Set argument value.
    /// UPDATE: just required by AdaptableProposal.
    void set_arg(size_t i, double x) { arguments_[i].set(x); }
    /// Set lower value.
    /// UPDATE: just required by AdaptableProposal.
    void set_lower(double x) { lower_arg_.set(x); }
    /// Set upper value.
    /// UPDATE: just required by AdaptableProposal.
    void set_upper(double x) { upper_arg_.set(x); }

    //////////////////////
    //////  GETTER  //////
    //////////////////////

    size_t n_args() const noexcept { return arguments_.size(); }
    double arg(size_t i) const { return arguments_[i].get(); }
    double lower() const { return lower_arg_.get(); }
    double upper() const { return upper_arg_.get(); }

    const Attribute &get_loc() const noexcept { return loc_arg_; }
    const Attribute &get_lower() const noexcept { return lower_arg_; }
    const Attribute &get_upper() const noexcept { return upper_arg_; }

    /// Evaluate implicitely at 'x' value.
    double evaluate() const { return evaluate(x_arg_.get()); }

    /// Central method to return function's value at given x.
    double evaluate(double x) const
    {
        switch (eval_mode_)
        {
        case Mode::Eval::Func:
            return function_(delocate_(x));
        case Mode::Eval::Prim:
            return primitive_(delocate_(x));
        case Mode::Eval::Proba:
            return pdf(x);
        case Mode::Eval::Cumul:
            return cdf(x);
        }
        __builtin_unreachable();
    }

    /// Return value as for a probability density/mass function.
    double pdf(double x) const { return is_valid(x) ? function_(delocate_(x)) : 0.; }

    /// Return value as for a cumulative distribution function.
    double cdf(double x) const { return is_below(x) ? 0 : primitive_(delocate_(is_above(x) ? upper() : x)); }

    /////////////////////
    //////  FLAGS  //////
    /////////////////////

    /// Wheter return value is always constant.
    bool is_constant() const noexcept;
    /// Wheter return value is always discrete.
    virtual bool is_discrete() const noexcept { return discr_mode_ != Mode::Discr::None; }
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
        case Mode::Interv::None:
            return false;
        case Mode::Interv::Open:
            return delocate_(x) <= lower_arg_.get();
        case Mode::Interv::Closed:
            return delocate_(x) < lower_arg_.get();
        }
        __builtin_unreachable();
    }

    /// Returns TRUE if value x is above range maximum.
    bool is_above(double x) const
    {
        switch (above_mode_)
        {
        case Mode::Interv::None:
            return false;
        case Mode::Interv::Open:
            return delocate_(x) >= upper_arg_.get();
        case Mode::Interv::Closed:
            return delocate_(x) > upper_arg_.get();
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
    Attribute loc_arg_ = Attribute();       /// 'loc' argument
    Attribute lower_arg_ = Attribute();     /// lower bound
    Attribute upper_arg_ = Attribute();     /// upper bound

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
        case Mode::Discr::None:
            return raw_function_(x);
        case Mode::Discr::Center:
        {
            const double i = std::round(x);
            return raw_primitive_(i + 0.5) - raw_primitive_(i - 0.5);
        }
        case Mode::Discr::Left:
        {
            const double i = std::floor(x);
            return raw_primitive_(i + 1.) - raw_primitive_(i);
        }
        case Mode::Discr::Defined:
            return raw_function_(std::floor(x));
        }
        __builtin_unreachable();
    }

    /// Polymorphic log-function method.
    virtual double log_func_(double x) const
    {
        switch (discr_mode_)
        {
        case Mode::Discr::None:
            return raw_logfunc_(x);
        case Mode::Discr::Center:
        {
            const double i = std::round(x);
            return std::log(raw_primitive_(i + 0.5) - raw_primitive_(i - 0.5));
        }
        case Mode::Discr::Left:
        {
            const double i = std::floor(x);
            return std::log(raw_primitive_(i + 1.) - raw_primitive_(i));
        }
        case Mode::Discr::Defined:
            return raw_logfunc_(std::floor(x));
        }
        __builtin_unreachable();
    }

    /// Polymorphic primitive method.
    virtual double primitive_(double x) const
    {
        switch (discr_mode_)
        {
        case Mode::Discr::None:
            return raw_primitive_(x);
        case Mode::Discr::Center:
            return raw_primitive_(std::round(x) + 0.5);
        case Mode::Discr::Left:
            return raw_primitive_(std::floor(x) + 1.);
        case Mode::Discr::Defined:
            return raw_primitive_(std::floor(x));
        }
        __builtin_unreachable();
    }

    /// Polymorphic log_ratio method.
    virtual double log_ratio_(double xo, double xn) const
    {
        switch (discr_mode_)
        {
        case Mode::Discr::None:
            return raw_logratio_(xo, xn);
        case Mode::Discr::Center:
        {
            const double io = std::round(xo), in = std::round(xn);
            if (io == in)
                return 0.;
            return std::log(raw_primitive_(in + .5) - raw_primitive_(in - .5)) -
                   std::log(raw_primitive_(io + .5) - raw_primitive_(io - .5));
        }
        case Mode::Discr::Left:
        {
            const double io = std::floor(xo), in = std::floor(xn);
            if (io == in)
                return 0.;
            return std::log(raw_primitive_(in + 1.) - raw_primitive_(in)) -
                   std::log(raw_primitive_(io + 1.) - raw_primitive_(io));
        }
        case Mode::Discr::Defined:
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
        case Mode::Deloc::None:
            return x;
        case Mode::Deloc::Active:
            return x - loc_arg_.get();
        }
        __builtin_unreachable();
    }

    /// Relocate according to mode. Used for random number generation.
    double relocate_(double x) const
    {
        switch (deloc_mode_)
        {
        case Mode::Deloc::None:
            return x;
        case Mode::Deloc::Active:
            return x + loc_arg_.get();
        }
        __builtin_unreachable();
    }

    /////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////  Function : range [protected]  ////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////

    /// Internal code to define the type of range.
    bool init_range_min_(double min_val = 0., bool is_min_open = false);
    void init_range_helper_(bool is_min_open = false, bool is_max_open = false);
};

/// Compose registries to get a function from them.
std::unique_ptr<Function> compose_func_registries(std::string_view opt,
                                                  const std::string &err_type,
                                                  const std::vector<std::string> &registries,
                                                  bool is_prim = false);
