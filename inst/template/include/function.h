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
    //////////////////////////////////  Attribute : constructor  ////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////

    Attribute() {};
    virtual ~Attribute() = default;

    /////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////  Attribute : getter  //////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////

    /// Get attribute value
    double get() const { return (this->*get_)(); }
    /// Get pointer to attribute value.
    const double *get_ptr() const { return (this->*get_ptr_)(); }
    /// Get pointer to attribute function.
    const std::shared_ptr<Function> get_function() const { return pFunc_; }

    bool is_constant() const;
    bool is_scalar() const { return (!is_pointer()) && (!is_function()); }
    bool is_pointer() const { return _ptr_to_value != nullptr; }
    bool is_function() const { return pFunc_ != nullptr; }

    std::string str() const { return tools::str::format(get()); }
    std::string name() const { return _attribute_name; }

    /////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////  Attribute : setter  //////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////

    /// Set attribute value.
    void set(double x) { (this->*_do_set)(x); }
    void name(const std::string &name) { _attribute_name = name; }

    void init_scalar(double x);
    void init_pointer(const double *y);
    void init_function(const std::shared_ptr<Function> f);

private:
    std::string _attribute_name = HSP_NA_STRING; /// name of argument
    double value_ = HSP_NA_DOUBLE;               /// internal value
    const double *_ptr_to_value = nullptr;       /// pointer to value
    std::shared_ptr<Function> pFunc_ = nullptr;  /// pointer to function

    // Functional getter
    double (Attribute::*get_)() const = &Attribute::do_get_default_;
    /// Get attribute value from internal value.
    double do_get_default_() const { return value_; }
    /// Get attribute value from external value.
    double __get_from_pointer() const { return *(_ptr_to_value); }
    /// Get attribute value from external function.
    double __get_from_function() const;

    // Functional pointer getter
    const double *(Attribute::*get_ptr_)() const = &Attribute::do_get_ptr_default_;
    /// Get pointer to attribute value internal value.
    const double *do_get_ptr_default_() const { return &value_; }
    /// Get pointer to attribute value from external value.
    const double *__get_pointer_from_pointer() const { return _ptr_to_value; }
    /// Get pointer to attribute error.
    const double *__get_pointer_error() const { throw std::runtime_error(
        "Cannot return a pointer to attribute \"" + _attribute_name + "\"."); }

    // Functional setter
    void (Attribute::*_do_set)(double) = &Attribute::__set_default;
    /// Set attribute internal value.
    void __set_default(double x) { value_ = x; }
    /// Set error.
    /// note: so far set_arg is only used for adaptation.
    /// Pointer and fucntion attributes can't be set, but just assigned.
    void __set_error(double) { throw std::runtime_error(
        "Non scalar attribute \"" + _attribute_name + "\" cannot be changed."); }
};

/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////// *** Function : class *** /////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

/// Function parent class. Arguements stored in arrays.
class Function
{
public:
    /////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////  Function : constructor  ///////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////

    Function(const std::string &name = HSP_NA_STRING,
             const str_vec &arg_name = {},
             const dbl_vec &arg_value = {},
             double domain_min = HSP_NA_DOUBLE,
             double domain_max = HSP_NA_DOUBLE,
             bool is_min_open = true,
             bool is_max_open = true);
    virtual ~Function() = default;

    /// Sugar operator to call return_value.
    double operator()(double x) const { return return_value(x); }

    /////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////  Function : getter  /////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////

    /// Virtual function to return a value.
    virtual double return_value(double) const { return 0.; }
    /// By convention, the argument 'x' is automatically the last one.
    double return_value() const { return return_value(arg(n_args_)); }

    unsigned short assert_arg(const std::string &arg_name) const;
    double arg(unsigned short i) const { return arguments_[i].get(); }
    double get_min() const { return range_min_.get(); }
    double get_max() const { return range_max_.get(); }

    const double *get_ptr_to_arg(unsigned short i) const { return arguments_[i].get_ptr(); }
    const double *get_ptr_to_min() const { return range_min_.get_ptr(); }
    const double *get_ptr_to_max() const { return range_max_.get_ptr(); }

    /// Wheter all arguments are actually constants.
    bool is_constant() const;
    /// Returns TRUE if value x is below range minimum.
    bool is_below(double x) const { return (this->*is_below_)(x); }
    /// Returns TRUE if value x is above range maximum.
    bool is_above(double x) const { return (this->*is_above_)(x); }
    /// Returns TRUE if value x is within range.
    bool is_valid(double x) const { return (this->*is_valid_)(x); }

    /// Internal code to define the type of range.
    unsigned short which_valid_case() const;

    /////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////  Function : arguments & range  ////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////

    void init_attr(const std::string &arg_name, bool is_arg, double value);
    void init_attr(const std::string &arg_name, bool is_arg, const double *ptr_to_value);
    void init_attr(const std::string &arg_name, bool is_arg, std::shared_ptr<Function> f);

    /// Set argument value.
    void set_arg(unsigned short i, double x) { arguments_[i].set(x); }

    void init_uniform_range();
    void init_functional_range(bool is_min_open, bool is_max_open);

    /// Virtual function to initialize reshaping in Distribution, etc.
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

private:
    std::string function_name_ = HSP_NA_STRING; /// function name

    Attribute *arguments_ = nullptr; /// array with arguments
    unsigned short n_args_ = 0;      /// n. of arguments

    Attribute range_min_ = Attribute(); /// lower bound
    Attribute range_max_ = Attribute(); /// upper bound

    /////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////  Function : arguments & range [private]  ///////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////

    void init_arg_scalar_(const std::string &arg_name, double value);
    void init_arg_pointer_(const std::string &arg_name, const double *ptr_to_value);
    void init_arg_function_(const std::string &arg_name, const std::shared_ptr<Function> f);

    void init_range_scalar_(const std::string &range_name, double lower);
    void init_range_pointer_(const std::string &range_name, const double *ptr_to_lower);
    void init_range_function_(const std::string &range_name, const std::shared_ptr<Function> f);

    //===== VALID =====//

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
        : _name(name) {}
    virtual ~FunctionRegistry() = default;

    /////////////////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////  FunctionRegistry : getter  //////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////

    const std::string &name() const { return _name; }
    virtual str_vec get_options() const;
    virtual bool is_valid_option(const std::string &opt) const { return get_options() % opt; }
    virtual std::shared_ptr<Function> get(const std::string &opt) const;

    /////////////////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////  FunctionRegistry : setter  //////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////

    template <typename T>
    void register_func()
    {
        static_assert(std::is_base_of<Function, T>::value,
                      "T must be derived from Function");
        _factory_map[T().name()] = []()
        { return std::make_shared<T>(); };
    }

private:
    std::string _name;
    std::unordered_map<std::string, FactoryFunc> _factory_map = {};
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
                                                  const std::vector<const FunctionRegistry *> &registries);

#endif // HOUSEPI__FUNCTION__H
