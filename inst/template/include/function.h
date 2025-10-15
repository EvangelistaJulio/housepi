/*
 *        __
 *       / /_  ____  __  __ ____ __  ____  __
 *      / __ \/ __ \/ / / / ___/ _ \/ __ \/ /
 *     / / / / /_/ / /_/ (__  )  __/ /_/ / /
 *    /_/ /_/\____/\____/____/\___/ ____/_/
 *                               /_/
 *
 * Copyright (C) 2025 Francesco Camaglia, Institut Pasteur
 *
 * function.h
 *
 */

#ifndef HOUSEPI__FUNCTION__H
#define HOUSEPI__FUNCTION__H

#include <memory>
#include <functional>

#include "tools/types.h"
#include "tools/na_utils.h"

#ifdef _HOUSEPI_RCPP
#include <Rcpp.h>
#endif

/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////// *** Attribute : class *** ///////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

// forward declaration
class Function;

/// Attribute parent class.
class Attribute
{
public:
    /////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////  Attribute : constructor  ///////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////

    /// Empty attribute.
    Attribute() {};
    /// Numeric attribute.
    Attribute(double x);
    /// Pointer to double attribute.
    Attribute(observ_dbl_ptr y);
    /// Function attribute.
    Attribute(const std::shared_ptr<Function> f);
    /// Function attribute by registry choice.
    Attribute(const std::string &s, bool is_primitive = false);
    /// Deconstructor.
    virtual ~Attribute() = default;

    /////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////  Attribute : setter  //////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////

    /// Set attribute value.
    void set(double x) { (this->*do_set_)(x); }
    /// Set attribute name. Eventually also err type.
    void name(const std::string &s, const std::string &e = hsp::na::STRING)
    {
        name_ = s;
        if (!hsp::na::is(e))
            err_type_ = e;
    }

    /// Initialize attribute as numeric.
    void init_numeric(double x);
    /// Initialize attribute as pointer.
    void init_pointer(observ_dbl_ptr y);
    /// Initialize attribute as function.
    void init_function(const std::shared_ptr<Function> f);
    /// Initialize attribute as function given the function name.
    void init_function(const std::string &s, bool is_primitive = false);
    /// Transfer from other Attribute.
    void transfer(const Attribute &other);

    /////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////  Attribute : getter  //////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////

    /// Get attribute name.
    std::string name() const { return name_; }
    /// Get attribute name to display in errors.
    std::string err_name() const { return err_type_ + " \"" + name_ + "\""; }
    /// Get a description about the functional status.
    virtual std::string str_summary() const;

    /// Get attribute numeric value.
    double get() const { return (this->*get_)(); }
    /// Get pointer to attribute value.
    observ_dbl_ptr get_ptr() const { return (this->*get_ptr_)(); }
    /// Get pointer to attribute function.
    const std::shared_ptr<Function> get_function() const { return pFunc_; }

    /// Wheter it's constant.
    virtual bool is_constant() const;
    /// Wheter it's discrete.
    virtual bool is_discrete() const;

    /// Wheter it's a number.
    bool is_numeric() const { return !hsp::na::is(value_); }
    /// Wheter it's a pointer.
    bool is_pointer() const { return ptr_to_value_ != nullptr; }
    /// Wheter it's a function.
    bool is_function() const { return pFunc_ != nullptr; }
    /// Wheter it has no definitions.
    bool is_empty() const { return !is_numeric() && !is_pointer() && !is_function(); }

protected:
    std::string err_type_ = "attribute";        // type name displayed in errors
    std::string name_ = hsp::na::STRING;        /// name of argument
    double value_ = hsp::na::DOUBLE;            /// internal value
    observ_dbl_ptr ptr_to_value_ = nullptr;     /// pointer to value // UPDATE: nullptr might be dangerous
    std::shared_ptr<Function> pFunc_ = nullptr; /// pointer to function // UPDATE: nullptr might be dangerous

    // Functional getter
    double (Attribute::*get_)() const = &Attribute::do_get_default_;
    /// Get attribute value from internal value.
    double do_get_default_() const { return value_; }
    /// Get attribute value from external value.
    double do_get_from_pointer_() const { return *(ptr_to_value_); }
    /// Get attribute value from external function.
    double do_get_from_function_() const;
    /// Get attribute error for not callable functions.
    double do_get_error_() const;

    // Functional pointer getter
    observ_dbl_ptr (Attribute::*get_ptr_)() const = &Attribute::do_get_ptr_default_;
    /// Get pointer to attribute error.
    observ_dbl_ptr do_get_pointer_error_() const { throw std::runtime_error(
        "Tried to return a pointer to " + err_name() + "."); }
    /// Get pointer to attribute value internal value.
    observ_dbl_ptr do_get_ptr_default_() const { return &value_; }
    /// Get pointer to attribute value from external value.
    observ_dbl_ptr do_get_pointer_from_pointer_() const { return ptr_to_value_; }

    // Functional setter
    void (Attribute::*do_set_)(double) = &Attribute::do_set_default_;
    /// Set attribute internal value.
    virtual void do_set_default_(double x) { value_ = x; }
    /// Set error. Non numeric attributes can't be set, but just assigned.
    void do_set_error_(double) { throw std::runtime_error(
        "Tried to change non numeric " + err_name() + "."); }
};

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

    Function(const std::string &name = hsp::na::STRING,
             const str_vec &arg_name = {}, const dbl_vec &arg_value = {},
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
    void set_eval_function() { evaluate_ = &Function::eval_function_; }
    /// Select primitive method for return value.
    void set_eval_primitive() { evaluate_ = &Function::eval_primitive_; }
    /// Select call as a distributin method for return value.
    void set_eval_as_distr_call() { evaluate_ = &Function::eval_function_as_distr_; }
    /// Select primitive as a distributin method for return value.
    void set_eval_as_distr_primitive() { evaluate_ = &Function::eval_primitive_as_distr_; }

    /////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////  Function : getter  /////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////

    double arg(unsigned short i) const { return arguments_[i].get(); }
    double get_min() const { return range_min_.get(); }
    double get_max() const { return range_max_.get(); }

    observ_dbl_ptr get_ptr_to_arg(unsigned short i) const { return arguments_[i].get_ptr(); }
    observ_dbl_ptr get_ptr_to_min() const { return range_min_.get_ptr(); }
    observ_dbl_ptr get_ptr_to_max() const { return range_max_.get_ptr(); }

    /// Wheter it is empty.
    bool is_empty() const { return hsp::na::is(name()); }
    /// Wheter return value is always constant.
    bool is_constant() const;
    /// Wheter return value is always discrete.
    virtual bool is_discrete() const { return function_ != &Function::raw_function_; }
    /// Wheter it can be called implicitely, i. e. if a 'x' has been linked.
    bool is_implicit_callable() const { return !x_arg_.is_empty(); }

    /// Returns TRUE if value x is below range minimum.
    bool is_below(double x) const { return (this->*is_below_)(x); }
    /// Returns TRUE if value x is above range maximum.
    bool is_above(double x) const { return (this->*is_above_)(x); }
    /// Returns TRUE if value x is within range.
    bool is_valid(double x) const { return (this->*is_valid_)(x); }

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

    std::string name() const { return name_; }
    /// Get attribute name to display in errors.
    std::string err_name() const { return err_type_ + " \"" + name_ + "\""; }

    std::string str_range() const;
    std::string str_summary(bool show_x_arg = true) const;

#ifdef _HOUSEPI_RCPP
    Rcpp::NumericVector rcpp_export_args() const;
    Rcpp::NumericVector rcpp_export_range() const;
#endif

protected:
    std::string err_type_ = "function";                                      // type name displayed in errors
    std::string name_ = hsp::na::STRING;                                     /// function name
    std::vector<Attribute> arguments_ = {};                                  /// array with arguments
    Attribute x_arg_ = Attribute();                                          /// 'x' argument
    Attribute range_min_ = Attribute();                                      /// lower bound
    Attribute range_max_ = Attribute();                                      /// upper bound
    double (Function::*evaluate_)(double) const = &Function::eval_function_; /// return value method

    /// Find attribute and return. If doesn't exists, return an empty attribute.
    virtual Attribute &find_arg_(const std::string &arg_name);

    /// Return value of call.
    double eval_function_(double x) const { return (this->*function_)(x); }
    /// Return value as for a probability density/mass function.
    double eval_function_as_distr_(double x) const { return is_valid(x) ? eval_function_(x) : 0.; }
    /// Return value of primitive.
    double eval_primitive_(double x) const { return (this->*primitive_)(x); }
    /// Return value as for a cumulative distribution function.
    double eval_primitive_as_distr_(double x) const { return is_valid(x) ? eval_primitive_(x) : (is_above(x) ? eval_primitive_(get_max()) : 0.); }

    /////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////  Function : virtual methods [protected]  ///////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////

    /// Raw call method.
    virtual double raw_function_(double) const { throw std::invalid_argument(
        "Called undefined method \"function\" of " + err_name() + "."); }
    /// Raw log-call method.
    virtual double raw_logfunc_(double) const { throw std::invalid_argument(
        "Called undefined method \"logfunc\" of " + err_name() + "."); }
    /// Raw primitive method.
    virtual double raw_primitive_(double) const { throw std::invalid_argument(
        "Called undefined method \"primitive\" of " + err_name() + "."); }
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
    /////////////////////////////  Function : range [protected]  ////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////

    /// Internal code to define the type of range.
    unsigned short which_valid_case_() const;
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
    bool do_is_below_open_min_(double x) const { return (x <= range_min_.get()); }
    bool do_is_below_close_min_(double x) const { return (x < range_min_.get()); }
    bool do_is_valid_open_min_(double x) const { return (x > range_min_.get()); }
    bool do_is_valid_close_min_(double x) const { return (x >= range_min_.get()); }

    //----- range : maximum -----//
    bool do_is_above_open_max_(double x) const { return (x >= range_max_.get()); }
    bool do_is_above_close_max_(double x) const { return (x > range_max_.get()); }
    bool do_is_valid_open_max_(double x) const { return (x < range_max_.get()); }
    bool do_is_valid_close_max_(double x) const { return (x <= range_max_.get()); }

    //----- range : both -----//
    bool do_is_valid_both_open_(double x) const { return (x > range_min_.get()) && (x < range_max_.get()); }
    bool do_is_valid_open_min_close_max_(double x) const { return (x > range_min_.get()) && (x <= range_max_.get()); }
    bool do_is_valid_close_min_open_max_(double x) const { return (x >= range_min_.get()) && (x < range_max_.get()); }
    bool do_is_valid_both_close_(double x) const { return (x >= range_min_.get()) && (x <= range_max_.get()); }
};

/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////// *** FunctionRegistry : class *** //////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

class FunctionRegistry
{
public:
    /////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////  FunctionRegistry : constructor  ////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////

    FunctionRegistry(const std::string &name)
        : name_(name) {}
    virtual ~FunctionRegistry() = default;

    /////////////////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////  FunctionRegistry : getter  //////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////

    const std::string &name() const { return name_; }
    virtual str_vec get_options() const;
    virtual std::shared_ptr<Function> get(const std::string &opt, bool is_primitive) const;

    /////////////////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////  FunctionRegistry : setter  //////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////

    template <typename T>
    void register_func()
    {
        static_assert(std::is_base_of<Function, T>::value,
                      "T must be derived from Function");
        factory_map_[T().name()] = []()
        { return std::make_shared<T>(); };
    }

private:
    std::string name_;
    u_dict<std::function<std::shared_ptr<Function>()>> factory_map_ = {};
};

/// Separate functions registries.
#define DECLARE_HOUSEPI_FUNCTION_REGISTRY(ClassName, Description) \
    class ClassName : public FunctionRegistry                     \
    {                                                             \
    public:                                                       \
        static ClassName &inst();                                 \
                                                                  \
    private:                                                      \
        ClassName() : FunctionRegistry(Description) {}            \
    };
DECLARE_HOUSEPI_FUNCTION_REGISTRY(ArithmeticRegistry, "arithmetic function")
DECLARE_HOUSEPI_FUNCTION_REGISTRY(TranscendentalRegistry, "transcendental function")
DECLARE_HOUSEPI_FUNCTION_REGISTRY(UserFuncRegistry, "user function")

/// Register a Function derived class T to FungRegistry derived class R.
template <typename R, typename T>
struct StaticFunctionRegistatror
{
    StaticFunctionRegistatror()
    {
        static_assert(std::is_base_of<FunctionRegistry, R>::value,
                      "R must be derived from FunctionRegistry.");
        R::inst().template register_func<T>();
    }
};
#define REGISTER_HOUSEPI_FUNCTION(R, T)        \
    namespace                                  \
    {                                          \
        static StaticFunctionRegistatror<R, T> \
            housepi_register_##T##_to_##R;     \
    }

#define REGISTER_USER_FUNCTION(T) REGISTER_HOUSEPI_FUNCTION(UserFuncRegistry, T)

/// Compose registries to get a function from them.
std::shared_ptr<Function> compose_func_registries(const std::string &opt,
                                                  const std::string &name,
                                                  const std::vector<const FunctionRegistry *> &registries,
                                                  bool is_primitive = false);

#endif // HOUSEPI__FUNCTION__H
