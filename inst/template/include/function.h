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

#include <string>
#include <vector>

#include "tools.h"

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
    /////////////////////////////////////  Attribute : getter  //////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////

    /// Get attribute name.
    std::string name() const { return attribute_name_; }
    /// Get attribute numeric value.
    double get() const { return (this->*get_)(); }
    /// Get pointer to attribute value.
    observer_dbl_ptr get_ptr() const { return (this->*get_ptr_)(); }
    /// Get pointer to attribute function.
    const std::shared_ptr<Function> get_function() const { return pFunc_; }
    /// Get formatted value.
    std::string str() const { return tools::str::format(get()); }

    /// Wheter it's constant.
    bool is_constant() const;
    /// Wheter it's a scalar.
    bool is_scalar() const { return (!is_pointer()) && (!is_function()); }
    /// Wheter it's a pointer.
    bool is_pointer() const { return ptr_to_value_ != nullptr; }
    /// Wheter it's a function.
    bool is_function() const { return pFunc_ != nullptr; }

    /////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////  Attribute : setter  //////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////

    /// Set attribute name.
    void name(const std::string &name) { attribute_name_ = name; }
    /// Set attribute value.
    void set(double x) { (this->*do_set_)(x); }

    /// Initialize attribute as scalar.
    void init_scalar(double x);
    /// Initialize attribute as pointer.
    void init_pointer(observer_dbl_ptr y);
    /// Initialize attribute as function.
    void init_function(const std::shared_ptr<Function> f);

private:
    std::string attribute_name_ = HSP_NA_STRING; /// name of argument
    double value_ = HSP_NA_DOUBLE;               /// internal value
    observer_dbl_ptr ptr_to_value_ = nullptr;    /// pointer to value
    std::shared_ptr<Function> pFunc_ = nullptr;  /// pointer to function

    // Functional getter
    double (Attribute::*get_)() const = &Attribute::do_get_default_;
    /// Get attribute value from internal value.
    double do_get_default_() const { return value_; }
    /// Get attribute value from external value.
    double do_get_from_pointer_() const { return *(ptr_to_value_); }
    /// Get attribute value from external function.
    double do_get_from_function_() const;

    // Functional pointer getter
    observer_dbl_ptr (Attribute::*get_ptr_)() const = &Attribute::do_get_ptr_default_;
    /// Get pointer to attribute value internal value.
    observer_dbl_ptr do_get_ptr_default_() const { return &value_; }
    /// Get pointer to attribute value from external value.
    observer_dbl_ptr do_get_pointer_from_pointer_() const { return ptr_to_value_; }
    /// Get pointer to attribute error.
    observer_dbl_ptr do_get_pointer_error_() const { throw std::runtime_error(
        "Cannot return a pointer to attribute \"" + attribute_name_ + "\"."); }

    // Functional setter
    void (Attribute::*do_set_)(double) = &Attribute::do_set_default_;
    /// Set attribute internal value.
    void do_set_default_(double x) { value_ = x; }
    /// Set error.
    /// note: so far set_arg is only used for adaptation.
    /// Pointer and fucntion attributes can't be set, but just assigned.
    void do_set_error_(double) { throw std::runtime_error(
        "Non scalar attribute \"" + attribute_name_ + "\" cannot be changed."); }
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

    Function(const std::string &name = HSP_NA_STRING,
             const str_vec &arg_name = {}, const dbl_vec &arg_value = {},
             double domain_min = HSP_NA_DOUBLE, bool is_min_open = false,
             double domain_max = HSP_NA_DOUBLE, bool is_max_open = false);
    virtual ~Function() = default;

    /////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////  Function : return value  //////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////

    /// Sugar operator to call return_value.
    double operator()(double x) const { return return_value(x); }

    /// Central method to return function's value at given x.
    double return_value(double x) const { return (this->*return_value_)(x); }
    /// By convention, the argument 'x' is automatically the last one.
    double return_value() const { return (this->*return_value_)(arg(n_args_)); }

    /// Select call method for return value.
    void set_return_call() { return_value_ = call_; }
    /// Select primitive method for return value.
    void set_return_primitive() { return_value_ = primitive_; }
    /// Select call as a distributin method for return value.
    void set_return_as_distr_call() { return_value_ = &Function::do_as_distr_call_; }
    /// Select primitive as a distributin method for return value.
    void set_return_as_distr_primitive() { return_value_ = &Function::do_as_distr_primitive_; }

    /////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////  Function : getter  /////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////

    unsigned short assert_arg(const std::string &arg_name) const;
    double arg(unsigned short i) const { return arguments_[i].get(); }
    double get_min() const { return range_min_.get(); }
    double get_max() const { return range_max_.get(); }

    observer_dbl_ptr get_ptr_to_arg(unsigned short i) const { return arguments_[i].get_ptr(); }
    observer_dbl_ptr get_ptr_to_min() const { return range_min_.get_ptr(); }
    observer_dbl_ptr get_ptr_to_max() const { return range_max_.get_ptr(); }

    /// Wheter it is discretized.
    bool is_discrete() const { return call_ != &Function::do_call_; }
    /// Wheter all arguments are actually constants.
    bool is_constant() const;
    /// Returns TRUE if value x is below range minimum.
    bool is_below(double x) const { return (this->*is_below_)(x); }
    /// Returns TRUE if value x is above range maximum.
    bool is_above(double x) const { return (this->*is_above_)(x); }
    /// Returns TRUE if value x is within range.
    bool is_valid(double x) const { return (this->*is_valid_)(x); }

    /////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////  Function : setter  /////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////

    void init_attr(const std::string &arg_name, bool is_arg, double value);
    void init_attr(const std::string &arg_name, bool is_arg, observer_dbl_ptr ptr_to_value);
    void init_attr(const std::string &arg_name, bool is_arg, std::shared_ptr<Function> f);

    /// Set argument value.
    void set_arg(unsigned short i, double x) { arguments_[i].set(x); }

    /// Virtual function to initialize reshaping.
    virtual void init_reshaping(bool) {}

    /////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////  Function : show  //////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////

    std::string name() const { return function_name_; }
    std::string str_arg_names() const;
    std::string str_arg_values() const;

    std::string str_min() const { return range_min_.str(); }
    std::string str_max() const { return range_max_.str(); }

    std::string str_range() const;
    std::string str_summary() const;

#ifdef _HOUSEPI_RCPP
    Rcpp::NumericVector rcpp_export_args() const;
    Rcpp::NumericVector rcpp_export_range() const;
#endif

protected:
    std::string function_name_ = HSP_NA_STRING;              /// function name
    unsigned short n_args_ = 0;                              /// n. of arguments
    Attribute *arguments_ = nullptr;                         /// array with arguments
    Attribute range_min_ = Attribute();                      /// lower bound
    Attribute range_max_ = Attribute();                      /// upper bound
    double (Function::*return_value_)(double) const = call_; /// return value method

    /////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////  Function : virtual methods [protected]  ///////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////

    /// Raw call method.
    virtual double do_call_(double) const { throw std::invalid_argument(
        "Called undefined method \"call\" of function \"" + name() + "\"."); }
    /// Raw log-call method.
    virtual double do_log_call_(double) const { throw std::invalid_argument(
        "Called undefined method \"log_call\" of function \"" + name() + "\"."); }
    /// Raw primitive method.
    virtual double do_primitive_(double) const { throw std::invalid_argument(
        "Called undefined method \"primitive\" of function \"" + name() + "\"."); }
    /// Raw log-ratio method.
    virtual double do_log_ratio_(double, double) const { throw std::invalid_argument(
        "Called undefined method \"log_ratio\" of function \"" + name() + "\"."); }

    /// Return value as for a probability density/mass function.
    double do_as_distr_call_(double x) const { return is_valid(x) ? (this->*call_)(x) : 0.; }
    /// Return value as for a cumulative distribution function.
    double do_as_distr_primitive_(double x) const { return is_valid(x) ? (this->*primitive_)(x) : (is_above(x) ? (this->*primitive_)(get_max()) : 0.); }

    /////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////  Function : discretization [protected]  ///////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////

    /// Internal code to define the type of discretization.
    unsigned short which_discrete_case_() const;
    /// Initialize functional discretization.
    void init_functional_discretization_(unsigned short discr_type);

    /// Polymorphic call method.
    double (Function::*call_)(double) const = &Function::do_call_;
    /// Call for discrete defined function.
    double do_call_discrete_(double x) const { return do_call_(std::floor(x)); }
    /// Call for center discretization.
    virtual double do_call_discr_center_(double x) const
    {
        const double i = std::round(x);
        return do_primitive_(i + .5) - do_primitive_(i - .5);
    }
    /// Call for left discretization.
    virtual double do_call_discr_left_(double x) const
    {
        const double i = std::floor(x);
        return do_primitive_(i + 1.) - do_primitive_(i);
    }

    double (Function::*log_call_)(double) const = &Function::do_log_call_;
    double do_log_call_discrete_(double x) const { return do_log_call_(std::floor(x)); }
    double do_log_call_discr_center_(double x) const { return std::log(do_call_discr_center_(x)); }
    double do_log_call_discr_left_(double x) const { return std::log(do_call_discr_left_(x)); }

    double (Function::*primitive_)(double) const = &Function::do_primitive_;
    double do_primitive_discrete_(double x) const { return do_primitive_(std::floor(x)); }
    virtual double do_primitive_discr_center_(double x) const { return do_primitive_(std::round(x) + .5); }
    virtual double do_primitive_discr_left_(double x) const { return do_primitive_(std::floor(x) + 1.); }

    double (Function::*log_ratio_)(double, double) const = &Function::do_log_ratio_;
    double do_log_ratio_discrete_(double xo, double xn) const { return do_log_ratio_(std::floor(xo), std::floor(xn)); }
    virtual double do_log_ratio_discr_center_(double xo, double xn) const
    {
        const double io = std::round(xo), in = std::round(xn);
        return std::log((do_primitive_(in + .5) - do_primitive_(in - .5)) / (do_primitive_(io + .5) - do_primitive_(io - .5)));
    }
    virtual double do_log_ratio_discr_left_(double xo, double xn) const
    {
        const double io = std::floor(xo), in = std::floor(xn);
        return std::log((do_primitive_(in + 1.) - do_primitive_(in)) / (do_primitive_(io + 1.) - do_primitive_(io)));
    }

    /////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////  Function : arguments [protected]  //////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////

    void init_arg_scalar_(const std::string &arg_name, double value);
    void init_arg_pointer_(const std::string &arg_name, observer_dbl_ptr ptr_to_value);
    void init_arg_function_(const std::string &arg_name, const std::shared_ptr<Function> f);

    /////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////  Function : range [protected]  ////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////

    void init_range_scalar_(const std::string &range_name, double lower);
    void init_range_pointer_(const std::string &range_name, observer_dbl_ptr ptr_to_lower);
    void init_range_function_(const std::string &range_name, const std::shared_ptr<Function> f);

    /// Internal code to define the type of range.
    unsigned short which_valid_case_() const;
    void init_functional_range_(bool is_min_open, bool is_max_open);
    void init_uniform_range_();

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
    using FactoryFunc = std::function<std::shared_ptr<Function>()>;

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
    bool is_valid_option(const std::string &opt) const { return get_options() % opt; }
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
    u_dict<FactoryFunc> factory_map_ = {};
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
                                                  bool is_primitive,
                                                  const std::vector<const FunctionRegistry *> &registries);

#endif // HOUSEPI__FUNCTION__H
