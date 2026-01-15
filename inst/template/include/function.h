/* housepi - Epidemiological Analysis of Household Outbreaks.
 *
 * Copyright © 2026 Institut Pasteur (Paris), CNRS and INSERM.
 * Distributed under the terms of the GNU General Public License (GPLv3).
 * See the LICENSE file for details.
 */

#ifndef HOUSEPI__FUNCTION__H
#define HOUSEPI__FUNCTION__H

#include <string>
#include <vector>

#include "tools/nan.h"
#include "tools/types.h"
#include "attribute.h"

/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////// *** Function : class *** /////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

/// Function parent class. Arguments stored in arrays.
class Function
{
public:
    /////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////  Function : constructor  ///////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////

    Function(const std::string &name,
             const hsp::named_args &args = {},
             double domain_min = hsp::na::DOUBLE, bool is_min_open = false,
             double domain_max = hsp::na::DOUBLE, bool is_max_open = false);
    virtual ~Function() = default;

    /////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////  Function : return value  //////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////

    /// Sugar operator to call evaluate.
    double operator()(double x) const { return evaluate(x); }

    /// Central method to return function's value at given x.
    double evaluate(double x) const { return (this->*evaluate_)(x); }
    /// By convention, the argument 'x' is automatically the last one.
    double evaluate() const { return (this->*evaluate_)(x_arg_.get()); }

    /// Select call method for return value.
    void set_eval_call() { evaluate_ = &Function::eval_function_; }
    /// Select primitive method for return value.
    void set_eval_primitive() { evaluate_ = &Function::eval_primitive_; }
    /// Select call as a distributin method for return value.
    void set_eval_as_distr_call() { evaluate_ = &Function::eval_function_as_distr_; }
    /// Select primitive as a distributin method for return value.
    void set_eval_as_distr_primitive() { evaluate_ = &Function::eval_primitive_as_distr_; }

    /////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////  Function : getter  /////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////

    size_t n_args() const { return arguments_.size(); }
    double arg(unsigned short i) const { return arguments_[i].get(); }
    double lower() const { return lower_arg_.get(); }
    double upper() const { return upper_arg_.get(); }

    const Attribute &get_arg(unsigned short i) const { return arguments_[i]; }
    const Attribute &get_loc() const { return loc_arg_; }
    const Attribute &get_lower() const { return lower_arg_; }
    const Attribute &get_upper() const { return upper_arg_; }

    /// Wheter it is empty.
    bool is_empty() const;
    /// Wheter return value is always constant.
    bool is_constant() const;
    /// Wheter return value is always discrete.
    virtual bool is_discrete() const { return function_ != &Function::raw_function_; }
    /// Wheter all arguments are not empty.
    bool is_callable_explicitly() const;
    /// Wheter it can be called implicitely, i. e. if also 'x' has been linked.
    bool is_callable_implicitly() const { return !x_arg_.is_empty() && is_callable_explicitly(); }

    /// Collect names of missing arguments.
    std::vector<std::string> get_missing_argument() const;
    /// Throw exception if not callable explicitly.
    void assert_is_callable_explicitly() const;

    /// Returns TRUE if value x is below range minimum.
    bool is_below(double x) const { return (this->*is_below_)((this->*delocate_)(x)); }
    /// Returns TRUE if value x is above range maximum.
    bool is_above(double x) const { return (this->*is_above_)((this->*delocate_)(x)); }
    /// Returns TRUE if value x is within range.
    bool is_valid(double x) const { return (this->*is_valid_)((this->*delocate_)(x)); }

    /////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////  Function : setter  /////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////

    /// Initialize given function attribute.
    void init_attr(const std::string &key, const Attribute &attr);

    /// Set argument value.
    void set_arg(unsigned short i, double x) { arguments_[i].set(x); }

    /// Virtual function to initialize discretization.
    virtual void init_discretization() { throw std::invalid_argument(
        "Called undefined method discretization of " + err_name() + "."); }

    /// Virtual function to initialize truncation.
    /// NOTE: always to be call after init_discretization.
    virtual void init_truncation() { throw std::invalid_argument(
        "Called to undefined method truncation of " + err_name() + "."); }

    /////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////  Function : show  //////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////

    /// Get function name.
    std::string name() const { return name_; }
    /// Get argument name.
    std::string arg_name(unsigned short i) const { return arguments_[i].name(); }
    /// Get function name complete with type to display in errors.
    std::string err_name() const { return err_type_ + " '" + name_ + "'"; }

    std::string str_range() const;
    std::string str_summary(bool show_x_arg = true) const;

protected:
    std::string err_type_ = "function";                                      // type name displayed in errors
    std::string name_;                                                       /// function name
    std::vector<Attribute> arguments_ = {};                                  /// array with arguments
    Attribute x_arg_ = Attribute();                                          /// 'x' argument
    Attribute loc_arg_ = Attribute();                                        /// 'loc' argument
    Attribute lower_arg_ = Attribute();                                      /// lower bound
    Attribute upper_arg_ = Attribute();                                      /// upper bound
    double (Function::*evaluate_)(double) const = &Function::eval_function_; /// return value method

    /// Find attribute and return. If doesn't exists, return an empty attribute.
    virtual Attribute &find_arg_(const std::string &arg_name);

    /// Return value of call.
    double eval_function_(double x) const { return (this->*function_)((this->*delocate_)(x)); }

    /// Return value as for a probability density/mass function.
    double eval_function_as_distr_(double x) const
    {
        if (is_valid(x))
            return eval_function_(x);
        return 0.;
    }

    /// Return value of primitive.
    double eval_primitive_(double x) const { return (this->*primitive_)((this->*delocate_)(x)); }

    /// Return value as for a cumulative distribution function.
    double eval_primitive_as_distr_(double x) const
    {
        if (is_valid(x))
            return eval_primitive_(x);
        if (is_above(x))
            return eval_primitive_((this->*delocate_)(upper()));
        return 0;
    }

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

    /// Internal code to define the type of discretization.
    unsigned short which_discrete_case_() const;
    /// Initialize functional discretization.
    virtual void init_discr_helper_(unsigned short discr_type);

    /// Polymorphic function method.
    double (Function::*function_)(double) const = &Function::raw_function_;
    /// Call for discrete defined function.
    double do_function_discrete_(double x) const { return raw_function_(std::floor(x)); }
    /// Call for center discretization.
    virtual double do_function_discr_center_(double x) const
    {
        const double i = std::round(x);
        return raw_primitive_(i + .5) - raw_primitive_(i - .5);
    }
    /// Call for left discretization.
    virtual double do_function_discr_left_(double x) const
    {
        const double i = std::floor(x);
        return raw_primitive_(i + 1.) - raw_primitive_(i);
    }

    /// Polymorphic function method.
    double (Function::*log_func_)(double) const = &Function::raw_logfunc_;
    double do_logfunc_discrete_(double x) const { return raw_logfunc_(std::floor(x)); }
    double do_logfunc_discr_center_(double x) const { return std::log(do_function_discr_center_(x)); }
    double do_logfunc_discr_left_(double x) const { return std::log(do_function_discr_left_(x)); }

    double (Function::*primitive_)(double) const = &Function::raw_primitive_;
    double do_primitive_discrete_(double x) const { return raw_primitive_(std::floor(x)); }
    virtual double do_primitive_discr_center_(double x) const { return raw_primitive_(std::round(x) + .5); }
    virtual double do_primitive_discr_left_(double x) const { return raw_primitive_(std::floor(x) + 1.); }

    double (Function::*log_ratio_)(double, double) const = &Function::raw_logratio_;
    double do_logratio_discrete_(double xo, double xn) const { return raw_logratio_(std::floor(xo), std::floor(xn)); }
    virtual double do_logratio_discr_center_(double xo, double xn) const
    {
        const double io = std::round(xo), in = std::round(xn);
        if (io == in)
            return 0.;
        return std::log((raw_primitive_(in + .5) - raw_primitive_(in - .5)) / (raw_primitive_(io + .5) - raw_primitive_(io - .5)));
    }
    virtual double do_logratio_discr_left_(double xo, double xn) const
    {
        const double io = std::floor(xo), in = std::floor(xn);
        if (io == in)
            return 0.;
        return std::log((raw_primitive_(in + 1.) - raw_primitive_(in)) / (raw_primitive_(io + 1.) - raw_primitive_(io)));
    }

    /////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////  Function : lcoation [protected]  //////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////

    double (Function::*delocate_)(double) const = &Function::do_delocate_none_;
    double do_delocate_none_(double x) const { return x; }
    double do_delocate_(double x) const { return x - loc_arg_.get(); }

    /////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////  Function : range [protected]  ////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////

    /// Internal code to define the type of range.
    bool init_range_min_(double min_val = 0., bool is_min_open = false);
    void init_range_helper_(bool is_min_open = false, bool is_max_open = false);

    //---- alternative valid range methods -----//
    bool (Function::*is_below_)(double) const = &Function::do_is_below_never_;
    bool (Function::*is_above_)(double) const = &Function::do_is_above_never_;
    bool (Function::*is_valid_)(double) const = &Function::do_is_valid_always_;

    //----- range : none ------//
    bool do_is_below_never_(double) const { return false; }
    bool do_is_above_never_(double) const { return false; }
    bool do_is_valid_always_(double) const { return true; }

    //----- range : minimum -----//
    bool do_is_below_open_min_(double x) const { return (x <= lower_arg_.get()); }
    bool do_is_below_close_min_(double x) const { return (x < lower_arg_.get()); }
    bool do_is_valid_open_min_(double x) const { return (x > lower_arg_.get()); }
    bool do_is_valid_close_min_(double x) const { return (x >= lower_arg_.get()); }

    //----- range : maximum -----//
    bool do_is_above_open_max_(double x) const { return (x >= upper_arg_.get()); }
    bool do_is_above_close_max_(double x) const { return (x > upper_arg_.get()); }
    bool do_is_valid_open_max_(double x) const { return (x < upper_arg_.get()); }
    bool do_is_valid_close_max_(double x) const { return (x <= upper_arg_.get()); }

    //----- range : both -----//
    bool do_is_valid_both_open_(double x) const { return (x > lower_arg_.get()) && (x < upper_arg_.get()); }
    bool do_is_valid_open_min_close_max_(double x) const { return (x > lower_arg_.get()) && (x <= upper_arg_.get()); }
    bool do_is_valid_close_min_open_max_(double x) const { return (x >= lower_arg_.get()) && (x < upper_arg_.get()); }
    bool do_is_valid_both_close_(double x) const { return (x >= lower_arg_.get()) && (x <= upper_arg_.get()); }
};

#endif // HOUSEPI__FUNCTION__H
