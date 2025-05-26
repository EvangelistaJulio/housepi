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

    double get() const;
    const double *get_pointer() const;
    const std::shared_ptr<Function> get_function() const;
    bool is_scalar() const;
    bool is_pointer() const;
    bool is_function() const;
    const std::string show() const;
    const std::string show_name() const;

    //==================//
    //===== SETTER =====//
    //==================//

    void set(double x);
    void set_name(const std::string &name);
    void assign_scalar(double x);
    void assign_pointer(const double *y);
    void assign_function(const std::shared_ptr<Function> f);

private:
    std::string _attribute_name = HOUSEPI_NA_STRING; // name of argument
    double _value = HOUSEPI_NA_DOUBLE;               // internal value
    const double *_ptr_to_value = nullptr;           // pointer to value
    std::shared_ptr<Function> _func_ptr = nullptr;   // pointer to function

    // Functional getter //
    double (Attribute::*_do_get)() const = &Attribute::__get_default;
    double __get_default() const;
    double __get_from_pointer() const;
    double __get_from_function() const;

    // Functional pointer getter //
    const double *(Attribute::*_do_get_pointer)() const = &Attribute::__get_pointer_default;
    const double *__get_pointer_default() const;
    const double *__get_pointer_from_pointer() const;
    const double *__get_pointer_error() const;

    // Functional setter //
    void (Attribute::*_do_set)(double) = &Attribute::__set_default;
    void __set_default(double x);
    void __set_error(double x);
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

    Function(const std::string &name = HOUSEPI_NA_STRING,
             const str_vec &arg_name = {},
             const dbl_vec &arg_value = {},
             double domain_min = HOUSEPI_NA_DOUBLE,
             double domain_max = HOUSEPI_NA_DOUBLE,
             bool is_min_open = true,
             bool is_max_open = true);
    virtual ~Function();

    //==================//
    //===== GETTER =====//
    //==================//

    virtual double return_value() const;

    ushort n_args() const;
    ushort find(const std::string &arg_name) const;
    ushort assert_arg(const std::string &arg_name) const;
    double get_arg(ushort i) const;
    double get_min() const;
    double get_max() const;
    const double *get_ptr_to_arg(ushort i) const;
    const double *get_ptr_to_min() const;
    const double *get_ptr_to_max() const;
    bool below(double x) const;
    bool above(double x) const;
    bool valid(double x) const;
    ushort which_valid_case() const;

    //=============================//
    //===== ARGUMENTS & RANGE =====//
    //=============================//

    void assign_attr(const std::string &arg_name, bool is_arg, double value);
    void assign_attr(const std::string &arg_name, bool is_arg, const double *ptr_to_value);
    void assign_attr(const std::string &arg_name, bool is_arg, std::shared_ptr<Function> f);

    void set_arg(ushort i, double new_value);
    void init_uniform_range();
    void init_functional_range(bool is_min_open, bool is_max_open);

    //===== SHOW =====//
    const std::string show_name() const;
    const std::string show_arg_names() const;
    const std::string show_arg_values() const;
    const std::string show_min() const;
    const std::string show_max() const;
    const std::string show_range() const;
    const std::string show_summary() const;

private:
    std::string _function_name = HOUSEPI_NA_STRING;

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
    bool __below_never(double x) const
    {
        UNUSED(x);
        return false;
    }
    bool __above_never(double x) const
    {
        UNUSED(x);
        return false;
    }
    bool __valid_always(double x) const
    {
        UNUSED(x);
        return true;
    }

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

class FuncRegistry
{
public:
    //=======================//
    //===== CONSTRUCTOR =====//
    //=======================//

    FuncRegistry(const std::string &name, const str_vec &options)
        : _name(name), _options(options) {}
    virtual ~FuncRegistry() = default;

    //==================//
    //===== GETTER =====//
    //==================//

    bool valid_option(const std::string &opt) const { return contains(_options, opt); }
    const str_vec &get_options() const { return _options; }
    virtual std::shared_ptr<Function> get(const std::string &opt) const = 0;

protected:
    std::string _name;
    str_vec _options;
};

/////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////  CompositeRegistry : class  /////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

class CompositeFuncRegistry : public FuncRegistry
{
private:
    std::vector<std::shared_ptr<FuncRegistry>> _registries;

public:
    CompositeFuncRegistry(const std::string &name,
                          std::vector<std::shared_ptr<FuncRegistry>> registries)
        : FuncRegistry(name, {}), _registries(std::move(registries))
    {
        for (const auto &r : _registries)
            _options += r->get_options();
    }

    std::shared_ptr<Function> get(const std::string &opt) const override final
    {
        if (tools::is_na(opt))
            return std::make_shared<Function>();
        else
            for (const auto &r : _registries)
                if (r->valid_option(opt))
                    return r->get(opt);

        throw std::invalid_argument(
            "Unknown " + _name + " \"" + opt + "\"." + tools::find_closest(opt, _options));

        return nullptr;
    }
};

/////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////  ArithmeticRegistry : class  ////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

class ArithmeticRegistry : public FuncRegistry
{
public:
    ArithmeticRegistry()
        : FuncRegistry("algebraic",
                       {"eq", "add", "sub", "mult", "div"}) {}

    std::shared_ptr<Function> get(const std::string &name) const override final;
};

/////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////  TranscendentalRegistry : class  ////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

class TranscendentalRegistry : public FuncRegistry
{
public:
    TranscendentalRegistry()
        : FuncRegistry("algebraic",
                       {"pow", "exp", "log"}) {}

    std::shared_ptr<Function> get(const std::string &name) const override final;
};

#endif // HOUSEPI__FUNCTION__H
