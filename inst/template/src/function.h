/*      __
 *     / /_  ____  __  __ ____ __  ____  __
 *    / __ \/ __ \/ / / / ___/ _ \/ __ \/ /
 *   / / / / /_/ / /_/ (__  )  __/ /_/ / /
 *  /_/ /_/\____/\____/____/\___/ ____/_/
 *                             /_/
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

//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//

// forward declaration
class Function;

/////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////  Attribute : class ////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

/* Attribute parent class. */
class Attribute
{
public:
    //=======================//
    //===== CONSTRUCTOR =====//
    //=======================//

    Attribute() {};
    virtual ~Attribute() = default;

    //==================//
    //===== GETTER =====//
    //==================//

    /* Get attribute value */
    double get() const { return (this->*_do_get)(); }
    /* Get pointer to attribute value. */
    const double *get_pointer() const { return (this->*_do_get_pointer)(); }
    /* Get pointer to attribute function. */
    const std::shared_ptr<Function> get_function() const { return _func_ptr; }

    bool is_scalar() const { return (!is_pointer()) && (!is_function()); }
    bool is_pointer() const { return _ptr_to_value != nullptr; }
    bool is_function() const { return _func_ptr != nullptr; }

    const std::string show() const { return tools::na_to_str(get()); }
    const std::string show_name() const { return _attribute_name; }

    //==================//
    //===== SETTER =====//
    //==================//

    /* Set attribute value. */
    void set(double x) { (this->*_do_set)(x); }

    void set_name(const std::string &name) { _attribute_name = name; }

    void assign_scalar(double x);
    void assign_pointer(const double *y);
    void assign_function(const std::shared_ptr<Function> f);

private:
    std::string _attribute_name = HSP_NA_STRING;   // name of argument
    double _value = HSP_NA_DOUBLE;                 // internal value
    const double *_ptr_to_value = nullptr;         // pointer to value
    std::shared_ptr<Function> _func_ptr = nullptr; // pointer to function

    // Functional getter //
    double (Attribute::*_do_get)() const = &Attribute::__get_default;
    /* Get attribute value from internal value. */
    double __get_default() const { return _value; }
    /* Get attribute value from external value. */
    double __get_from_pointer() const { return *(_ptr_to_value); }
    /* Get attribute value from external function. */
    double __get_from_function() const;

    // Functional pointer getter //
    const double *(Attribute::*_do_get_pointer)() const = &Attribute::__get_pointer_default;
    /* Get pointer to attribute value internal value. */
    const double *__get_pointer_default() const { return &_value; }
    /* Get pointer to attribute value from external value. */
    const double *__get_pointer_from_pointer() const { return _ptr_to_value; }
    /* Get pointer to attribute error. */
    const double *__get_pointer_error() const { throw std::runtime_error(
        "Cannot return a pointer to attribute \"" + show_name() + "\"."); }

    // Functional setter //
    void (Attribute::*_do_set)(double) = &Attribute::__set_default;
    /* Set attribute internal value. */
    void __set_default(double x) { _value = x; }
    /* Set error.
     * note: so far set_arg is only used for adaptation.
     * Pointer and fucntion attributes can't be set, but just assigned.
     */
    void __set_error(double) { throw std::runtime_error(
        "Non scalar attribute \"" + show_name() + "\" cannot be changed."); }
};

/////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////  Base Function : class //////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

/* Function parent class. Arguements stored in arrays. */
class Function
{
public:
    //=======================//
    //===== CONSTRUCTOR =====//
    //=======================//

    Function(const std::string &name = HSP_NA_STRING,
             const str_vec &arg_name = {},
             const dbl_vec &arg_value = {},
             double domain_min = HSP_NA_DOUBLE,
             double domain_max = HSP_NA_DOUBLE,
             bool is_min_open = true,
             bool is_max_open = true);
    virtual ~Function();

    double operator()(double x) const { return return_value(x); }

    //==================//
    //===== GETTER =====//
    //==================//

    /* Virtual function to return a value */
    virtual double return_value(double) const { return 0.; }
    /* By convention, the argument 'x' is automatically the last one. */
    double return_value() const { return return_value(arg(_n_args)); }

    bool is_able_to_showcase() const;
    ushort assert_arg(const std::string &arg_name) const;
    double arg(ushort i) const { return _arguments[i].get(); }
    double get_min() const { return _range_min.get(); }
    double get_max() const { return _range_max.get(); }

    const double *get_ptr_to_arg(ushort i) const { return _arguments[i].get_pointer(); }
    const double *get_ptr_to_min() const { return _range_min.get_pointer(); }
    const double *get_ptr_to_max() const { return _range_max.get_pointer(); }

    /* Returns TRUE if value x is below range minimum. */
    bool below(double x) const { return (this->*_do_below)(x); }
    /* Returns TRUE if value x is above range maximum. */
    bool above(double x) const { return (this->*_do_above)(x); }
    /* Returns TRUE if value x is within range. */
    bool valid(double x) const { return (this->*_do_valid)(x); }

    /* Internal code to define the type of range. */
    ushort which_valid_case() const;

    //=============================//
    //===== ARGUMENTS & RANGE =====//
    //=============================//

    void assign_attr(const std::string &arg_name, bool is_arg, double value);
    void assign_attr(const std::string &arg_name, bool is_arg, const double *ptr_to_value);
    void assign_attr(const std::string &arg_name, bool is_arg, std::shared_ptr<Function> f);

    /* Set argument value. */
    void set_arg(ushort i, double x) { _arguments[i].set(x); }

    void init_uniform_range();
    void init_functional_range(bool is_min_open, bool is_max_open);

    /* Virtual function to initialize reshaping in Distribution, etc. */
    virtual void init_reshaping(bool) {}

    //===== SHOW =====//
    /* Return function name. */
    const std::string show_name() const { return _function_name; }
    const std::string show_arg_names() const;
    const std::string show_arg_values() const;
    /* Return formatted range minimum. */
    const std::string show_min() const { return _range_min.show(); }
    /* Return formatted range maximum. */
    const std::string show_max() const { return _range_max.show(); }
    const std::string show_range() const;
    const std::string show_summary() const;

private:
    std::string _function_name = HSP_NA_STRING;

    //===== Arguments =====//
    Attribute *_arguments = nullptr;
    ushort _n_args = 0;

    //===== Domain =====//
    Attribute _range_min = Attribute();
    Attribute _range_max = Attribute();

    //=============================//
    //===== ARGUMENTS & RANGE =====//
    //=============================//
    void _assign_arg_scalar(const std::string &arg_name, double value);
    void _assign_arg_pointer(const std::string &arg_name, const double *ptr_to_value);
    void _assign_arg_function(const std::string &arg_name, const std::shared_ptr<Function> f);

    void _assign_range_scalar(const std::string &range_name, double lower);
    void _assign_range_pointer(const std::string &range_name, const double *ptr_to_lower);
    void _assign_range_function(const std::string &range_name, const std::shared_ptr<Function> f);

    //===== VALID =====//
    // note: for the moment range extremes can't be one open and other close

    //---- alternative valid range methods -----//
    bool (Function::*_do_below)(double) const = &Function::__below_never;
    bool (Function::*_do_above)(double) const = &Function::__above_never;
    bool (Function::*_do_valid)(double) const = &Function::__valid_always;

    //----- range : none ------//
    bool __below_never(double) const { return false; }
    bool __above_never(double) const { return false; }
    bool __valid_always(double) const { return true; }

    //----- range : minimum -----//
    bool __below_open_min(double x) const { return (x <= _range_min.get()); }
    bool __below_close_min(double x) const { return (x < _range_min.get()); }
    bool __valid_open_min(double x) const { return (x > _range_min.get()); }
    bool __valid_close_min(double x) const { return (x >= _range_min.get()); }

    //----- range : maximum -----//
    bool __above_open_max(double x) const { return (x >= _range_max.get()); }
    bool __above_close_max(double x) const { return (x > _range_max.get()); }
    bool __valid_open_max(double x) const { return (x < _range_max.get()); }
    bool __valid_close_max(double x) const { return (x <= _range_max.get()); }

    //----- range : both -----//
    bool __valid_both_open(double x) const { return (x > _range_min.get()) && (x < _range_max.get()); }
    bool __valid_open_min_close_max(double x) const { return (x > _range_min.get()) && (x <= _range_max.get()); }
    bool __valid_close_min_open_max(double x) const { return (x >= _range_min.get()) && (x < _range_max.get()); }
    bool __valid_both_close(double x) const { return (x >= _range_min.get()) && (x <= _range_max.get()); }
};

/////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////  Registry : class  /////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

class FunctionRegistry
{
    using FactoryFunc = std::function<std::shared_ptr<Function>()>;

public:
    //=======================//
    //===== CONSTRUCTOR =====//
    //=======================//

    FunctionRegistry(const std::string &name)
        : _name(name) {}
    virtual ~FunctionRegistry() = default;

    //==================//
    //===== GETTER =====//
    //==================//

    const std::string &name() const { return _name; }
    virtual const str_vec get_options() const;
    virtual bool is_valid_option(const std::string &opt) const { return contains(get_options(), opt); }
    virtual std::shared_ptr<Function> get(const std::string &opt) const;

    //==================//
    //===== SETTER =====//
    //==================//

    template <typename T>
    void register_func()
    {
        static_assert(std::is_base_of<Function, T>::value,
                      "T must be derived from Function");
        _factory_map[T().show_name()] = []()
        { return std::make_shared<T>(); };
    }

private:
    std::string _name;
    std::unordered_map<std::string, FactoryFunc> _factory_map = {};
};

/* Separate functions registries. */
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

/* Register a Function derived class T to FungRegistry derived class R. */
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

/* Compose registries to get a function from them. */
std::shared_ptr<Function> compose_func_registries(const std::string &opt,
                                                  const std::string &name,
                                                  const std::vector<const FunctionRegistry *> &registries);

#endif // HOUSEPI__FUNCTION__H
