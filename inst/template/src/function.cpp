/*      __
 *     / /_  ____  __  __ ____ __  ____  __
 *    / __ \/ __ \/ / / / ___/ _ \/ __ \/ /
 *   / / / / /_/ / /_/ (__  )  __/ /_/ / /
 *  /_/ /_/\____/\____/____/\___/ ____/_/
 *                             /_/
 *
 * Copyright (C) 2025 Francesco Camaglia, Institut Pasteur
 *
 * function.cpp
 *
 */

#include "function.h"

//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//

/////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////  Attribute : getter //////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

bool Attribute::is_scalar() const { return (!is_pointer()) && (!is_function()); }
bool Attribute::is_pointer() const { return _ptr_to_value != nullptr; }
bool Attribute::is_function() const { return _func_ptr != nullptr; }
const std::string Attribute::show() const { return tools::na_to_str(get()); }
const std::string Attribute::show_name() const { return _attribute_name; }

/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////  Attribute : get (functional) /////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

/* Get attribute value */
inline double Attribute::get() const { return (this->*_do_get)(); }

/* Get attribute value internal value. */
double Attribute::__get_default() const { return _value; }

/* Get attribute value from external value. */
double Attribute::__get_from_pointer() const { return *(_ptr_to_value); }

/* Set attribute value internal value. */
// FIXME: return_value should never be dependent on anything
double Attribute::__get_from_function() const { return _func_ptr->return_value({}); }

/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////  Attribute : get_pointer (functional) /////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

/* Get pointer to attribute value. */
inline const double *Attribute::get_pointer() const { return (this->*_do_get_pointer)(); }

/* Get pointer to attribute value internal value. */
const double *Attribute::__get_pointer_default() const { return &_value; }

/* Get pointer to attribute value from external value. */
const double *Attribute::__get_pointer_from_pointer() const { return _ptr_to_value; }

/* Get pointer to attribute error. */
const double *Attribute::__get_pointer_error() const
{
    throw std::runtime_error(
        "Cannot return a pointer to attribute \"" + _attribute_name + "\".");
    return nullptr;
}

/////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////  Attribute : get_function ////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

inline const std::shared_ptr<Function> Attribute::get_function() const { return _func_ptr; }

/////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////  Attribute : setter //////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

inline void Attribute::set_name(const std::string &name) { _attribute_name = name; }
void Attribute::assign_scalar(double x)
{
    _do_get = &Attribute::__get_default;
    _do_get_pointer = &Attribute::__get_pointer_default;
    _do_set = &Attribute::__set_default;
    _value = x;
    _ptr_to_value = nullptr;
    _func_ptr = nullptr;
}
void Attribute::assign_pointer(const double *y)
{
    _do_get = &Attribute::__get_from_pointer;
    _do_get_pointer = &Attribute::__get_pointer_from_pointer;
    _do_set = &Attribute::__set_error;
    _value = HOUSEPI_NA_DOUBLE;
    _ptr_to_value = y;
    _func_ptr = nullptr;
}
void Attribute::assign_function(const std::shared_ptr<Function> f)
{
    _do_get = &Attribute::__get_from_function;
    _do_get_pointer = &Attribute::__get_pointer_error;
    _do_set = &Attribute::__set_error;
    _value = HOUSEPI_NA_DOUBLE;
    _ptr_to_value = nullptr;
    _func_ptr = f;
}

/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////  Attribute : set (functional) /////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

/* Set attribute value. */
inline void Attribute::set(double x) { (this->*_do_set)(x); }

/* Set attribute internal value. */
void Attribute::__set_default(double x) { _value = x; }

/* Set error. */
void Attribute::__set_error(double x)
{
    // note: so far set_arg is only used for adaptation.
    // Pointer and fucntion attributes can't be set, but just assigned.
    UNUSED(x);
    throw std::runtime_error(
        "Non scalar attribute \"" + _attribute_name + "\" cannot be changed.");
}

//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//

/////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////  Base Function : Getter  ///////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

/* Virtual function to return a value */
double Function::return_value(const dbl_vec &input) const
{
    UNUSED(input);
    return 0.;
}

/* Auxiliary function to retrieve the index of a given attribute by name. */
inline ushort Function::find(const std::string &arg_name) const
{
    for (ushort c = 0; c < _n_args; ++c)
        if (_arguments[c].show_name() == arg_name)
            return c;
    return _n_args;
}

/* Assert argument exists at assignment. Returns it's index. */
ushort Function::assert_arg(const std::string &arg_name) const
{
    ushort i = this->find(arg_name);
    if (i == _n_args)
    {
        str_vec options = {};
        for (ushort c = 0; c < _n_args; ++c)
            options.push_back(_arguments[c].show_name());
        throw std::invalid_argument(
            "Unknown argument \"" + arg_name + "\" for \"" + _function_name + "\"." + tools::find_closest(arg_name, options));
    }
    return i;
}

double Function::get_arg(ushort i) const { return _arguments[i].get(); }
double Function::get_min() const { return _range_min.get(); }
double Function::get_max() const { return _range_max.get(); }

const double *Function::get_ptr_to_arg(ushort i) const { return _arguments[i].get_pointer(); }
const double *Function::get_ptr_to_min() const { return _range_min.get_pointer(); }
const double *Function::get_ptr_to_max() const { return _range_max.get_pointer(); }

/* Returns TRUE if value x is below range minimum. */
bool Function::below(double x) const { return (this->*_do_below)(x); }
/* Returns TRUE if value x is above range maximum. */
bool Function::above(double x) const { return (this->*_do_above)(x); }
/* Returns TRUE if value x is within range. */
bool Function::valid(double x) const { return (this->*_do_valid)(x); }

/////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////  Base Function : set range  ////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

/* Returns the case of the range of validity. */
ushort Function::which_valid_case() const
{
    // note: bounded if not pointing internally to NA
    bool is_lower_bound = !(_range_min.is_scalar() && tools::is_na(_range_min.get()));
    bool is_upper_bound = !(_range_max.is_scalar() && tools::is_na(_range_max.get()));
    return is_lower_bound + 2 * is_upper_bound;
}

/* Ensure uniform range (min, max) match argument (min, max). */
void Function::init_uniform_range()
{
    // impose range minimum equal to argument minimum
    Attribute &min_arg = _arguments[this->assert_arg("min")];
    if (min_arg.is_function())
        _range_min.assign_function(min_arg.get_function());
    else
        _range_min.assign_pointer(min_arg.get_pointer());

    // impose range maximum equal to argument maximum
    Attribute &max_arg = _arguments[this->assert_arg("max")];
    if (max_arg.is_function())
        _range_max.assign_function(max_arg.get_function());
    else
        _range_max.assign_pointer(max_arg.get_pointer());
}

/* Set the range of validity of the value, and the corresponding validation function. */
void Function::init_functional_range(bool is_min_open, bool is_max_open)
{
    switch (which_valid_case())
    {
    case 0: //----- unbounded -----//
        _do_below = &Function::__below_never;
        _do_above = &Function::__above_never;
        _do_valid = &Function::__valid_always;
        break;

    case 1: //----- lower bound -----//
        _do_above = &Function::__above_never;
        if (is_min_open)
        {
            _do_below = &Function::__below_open_min;
            _do_valid = &Function::__valid_open_min;
        }
        else
        {
            _do_below = &Function::__below_close_min;
            _do_valid = &Function::__valid_close_min;
        }
        break;

    case 2: //----- upper bound -----//
        _do_below = &Function::__below_never;
        if (is_max_open)
        {
            _do_above = &Function::__above_open_max;
            _do_valid = &Function::__valid_open_max;
        }
        else
        {
            _do_above = &Function::__above_close_max;
            _do_valid = &Function::__valid_close_max;
        }
        break;

    case 3: //----- within an interval -----//

        if (_range_min.get() >= _range_max.get())
            throw std::invalid_argument(
                "Inconsistent range " + show_range() + ".");

        if ((is_min_open) && (is_max_open))
        {
            _do_below = &Function::__below_open_min;
            _do_above = &Function::__above_open_max;
            _do_valid = &Function::__valid_both_open;
        }
        else if ((is_min_open) && (!is_max_open))
        {
            _do_below = &Function::__below_open_min;
            _do_above = &Function::__above_close_max;
            _do_valid = &Function::__valid_open_min_close_max;
        }
        else if ((!is_min_open) && (is_max_open))
        {
            _do_below = &Function::__below_close_min;
            _do_above = &Function::__above_open_max;
            _do_valid = &Function::__valid_close_min_open_max;
        }
        else
        {
            _do_below = &Function::__below_close_min;
            _do_above = &Function::__above_close_max;
            _do_valid = &Function::__valid_both_close;
        }
        break;
    }
}

const std::string Function::show_min() const { return _range_min.show(); }
const std::string Function::show_max() const { return _range_max.show(); }

/* Return formatted range. */
const std::string Function::show_range() const
{
    if (_do_valid == &Function::__valid_always)
        return "";

    else if (_do_valid == &Function::__valid_open_min)
        return "(" + show_min() + ", Inf)";

    else if (_do_valid == &Function::__valid_close_min)
        return "[" + show_min() + ", Inf)";

    else if (_do_valid == &Function::__valid_open_max)
        return "(-Inf, " + show_max() + ")";

    else if (_do_valid == &Function::__valid_close_max)
        return "(-Inf, " + show_max() + "]";

    else if (_do_valid == &Function::__valid_both_open)
        return "(" + show_min() + ", " + show_max() + ")";

    else if (_do_valid == &Function::__valid_open_min_close_max)
        return "(" + show_min() + ", " + show_max() + "]";

    else if (_do_valid == &Function::__valid_close_min_open_max)
        return "[" + show_min() + ", " + show_max() + ")";

    else // (_do_valid == &Function::__valid_both_close)
        return "[" + show_min() + ", " + show_max() + "]";
}

//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//

/////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////  Base Function : Arguments  /////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

/* Set argument value. */
void Function::set_arg(ushort i, double x) { _arguments[i].set(x); }

//--- Assign scalar argument or range value. ---//
void Function::assign_attr_scalar(const std::string &arg_name, bool is_arg, double x)
{
    if (is_arg == ATTR_ARGUMENT)
        _assign_arg_scalar(arg_name, x);
    else
        _assign_range_scalar(arg_name, x);
}
//* Assign pointer argument or range value. *//
void Function::assign_attr_pointer(const std::string &arg_name, bool is_arg, const double *y)
{
    // check the function is different from the current one
    if (y == nullptr)
        throw std::invalid_argument(
            "Invalid pointer argument \"" + arg_name + "\" for \"" + _function_name + "\".");

    if (is_arg == ATTR_ARGUMENT)
        _assign_arg_pointer(arg_name, y);
    else
        _assign_range_pointer(arg_name, y);
}
//* Assign fucntion argument or range value. *//
void Function::assign_attr_function(const std::string &arg_name, bool is_arg, std::shared_ptr<Function> f)
{
    // check the function is different from the current one
    if (f == nullptr)
        throw std::invalid_argument(
            "Invalid function argument \"" + arg_name + "\" for \"" + _function_name + "\".");
    else if (f.get() == this)
        throw std::invalid_argument(
            "Cannot assign function \"" + _function_name + "\" to its own argument \"" + arg_name + "\".");

    if (is_arg == ATTR_ARGUMENT)
        _assign_arg_function(arg_name, f);
    else
        _assign_range_function(arg_name, f);
}

/* Assignment scalar argument. */
void Function::_assign_arg_scalar(const std::string &arg_name, double x)
{
    ushort i = this->assert_arg(arg_name);
    _arguments[i].assign_scalar(x);
}
/* Assignment pointer argument. */
void Function::_assign_arg_pointer(const std::string &arg_name, const double *y)
{
    ushort i = this->assert_arg(arg_name);
    _arguments[i].assign_pointer(y);
}
/* Assignment function argument. */
void Function::_assign_arg_function(const std::string &arg_name, std::shared_ptr<Function> f)
{
    ushort i = this->assert_arg(arg_name);
    _arguments[i].assign_function(f);
}
/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////  Base Function : Range  ///////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

/* Assign scalar range. */
void Function::_assign_range_scalar(const std::string &range_name, double x)
{
    if (_function_name == "unif")
    {
        if ((range_name != "min") && (range_name != "max"))
            throw std::invalid_argument(
                "Unknown range name \"" + range_name + "\" for \"" + _function_name + "\".");
        else if (!tools::is_na(x))
            _assign_arg_scalar(range_name, x);
    }
    else
    {
        if (range_name == "min")
            _range_min.assign_scalar(x);
        else if (range_name == "max")
            _range_max.assign_scalar(x);
        else
            throw std::invalid_argument(
                "Unknown range name \"" + range_name + "\" for \"" + _function_name + "\".");
    }
}

/* Assign pointer range. */
void Function::_assign_range_pointer(const std::string &range_name, const double *y)
{
    if (_function_name == "unif")
    {
        if ((range_name != "min") && (range_name != "max"))
            throw std::invalid_argument(
                "Unknown range name \"" + range_name + "\" for \"" + _function_name + "\".");
        else if (y != nullptr)
            _assign_arg_pointer(range_name, y);
    }
    else
    {
        if (range_name == "min")
            _range_min.assign_pointer(y);
        else if (range_name == "max")
            _range_max.assign_pointer(y);
        else
            throw std::invalid_argument(
                "Unknown range name \"" + range_name + "\" for \"" + _function_name + "\".");
    }
}

/* Assign function range. */
void Function::_assign_range_function(const std::string &range_name, const std::shared_ptr<Function> f)
{
    if (_function_name == "unif")
    {
        if ((range_name != "min") && (range_name != "max"))
            throw std::invalid_argument(
                "Unknown range name \"" + range_name + "\" for \"" + _function_name + "\".");
        else if (f != nullptr)
            _assign_arg_function(range_name, f);
    }
    else
    {
        if (range_name == "min")
            _range_min.assign_function(f);
        else if (range_name == "max")
            _range_max.assign_function(f);
        else
            throw std::invalid_argument(
                "Unknown range name \"" + range_name + "\" for \"" + _function_name + "\".");
    }
}

//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//

/////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////  Base Function : Show  /////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

/* Return function name. */
const std::string Function::show_name() const { return _function_name; }

/* Return formatted argument names of the distribution. */
const std::string Function::show_arg_names() const
{
    if (_n_args > 0)
    {
        std::string output = _arguments[0].show_name();
        for (ushort c = 1; c < _n_args; ++c)
            output += C_LIST_DELIM + _arguments[c].show_name();
        return output;
    }

    return HOUSEPI_NA_STRING;
}

/* Return formatted argument names of the distribution. */
const std::string Function::show_arg_values() const
{
    if (_n_args > 0)
    {
        std::string output = _arguments[0].show();
        for (ushort c = 1; c < _n_args; ++c)
            output += C_LIST_DELIM + _arguments[c].show();
        return output;
    }

    return HOUSEPI_NA_STRING;
}

/* Return formatted function summary. */
const std::string Function::show_summary() const
{
    std::string output = show_name();
    if (_n_args > 0)
    {
        output += "(" + _arguments[0].show_name() + "=" + _arguments[0].show();
        for (ushort c = 1; c < _n_args; ++c)
            output += ", " + _arguments[c].show_name() + "=" + _arguments[c].show();
        output += ")";
    }
    if (_do_valid != &Function::__valid_always)
        output += " : "+ show_range();

    return output;
}

//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//

/* Initial defaults for function. */
void Function::init_defaults(const std::string &name,
                             const str_vec &arg_name,
                             const dbl_vec &arg_value,
                             double domain_min,
                             double domain_max,
                             bool is_min_open,
                             bool is_max_open)
{
    // name of the function
    _function_name = name;

    // number of arguments
    if (arg_name.size() != arg_value.size())
        throw std::length_error(
            "Internal inconsistency for \"" + _function_name + "\" arguments.");
    _n_args = arg_name.size();

    // assign default arguments
    _arguments = new Attribute[arg_name.size()];
    for (ushort c = 0; c < _n_args; ++c)
    {
        _arguments[c].set_name(arg_name[c]);
        _arguments[c].set(arg_value[c]);
    }

    //----- assign domain -----//
    // NOTE: at the moment only place where to set open ranges
    // NOTE: not sure if setting pointer to smtg exiting the domain I still have control
    _assign_range_scalar("min", domain_min);
    _assign_range_scalar("max", domain_max);
    init_functional_range(is_min_open, is_max_open);
}

//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//

/* Identity Function. */
class Equality_Function : public Function
{
public:
    //===== CONSTRUCTOR =====//
    Equality_Function()
    {
        init_defaults("identity",
                      {"x"}, {0.},
                      HOUSEPI_NA_DOUBLE, HOUSEPI_NA_DOUBLE,
                      false, false);
    }
    ~Equality_Function() = default;

    double return_value(const dbl_vec &input) const override final
    {
        UNUSED(input);
        return get_arg(0);
    }
};

//--- Addition Function. ---//
class Addition_Function : public Function
{
public:
    //===== CONSTRUCTOR =====//
    Addition_Function()
    {
        init_defaults("add",
                      {"x", "y"}, {0., 0.},
                      HOUSEPI_NA_DOUBLE, HOUSEPI_NA_DOUBLE,
                      false, false);
    }
    ~Addition_Function() = default;

    //==================//
    //===== GETTER =====//
    //==================//
    double return_value(const dbl_vec &input) const override final
    {
        UNUSED(input);
        return get_arg(0) + get_arg(1);
    }
};

//--- Subtraction Function. ---//
class Subtraction_Function : public Function
{
public:
    //===== CONSTRUCTOR =====//
    Subtraction_Function()
    {
        init_defaults("sum",
                      {"x", "y"}, {0., 0.},
                      HOUSEPI_NA_DOUBLE, HOUSEPI_NA_DOUBLE,
                      false, false);
    }
    ~Subtraction_Function() = default;

    //==================//
    //===== GETTER =====//
    //==================//
    double return_value(const dbl_vec &input) const override final
    {
        UNUSED(input);
        return get_arg(0) - get_arg(1);
    }
};

//--- Multiplication Function. ---//
class Multiplication_Function : public Function
{
public:
    //===== CONSTRUCTOR =====//
    Multiplication_Function()
    {
        init_defaults("mult",
                      {"x", "y"}, {0., 0.},
                      HOUSEPI_NA_DOUBLE, HOUSEPI_NA_DOUBLE,
                      false, false);
    }
    ~Multiplication_Function() = default;

    //==================//
    //===== GETTER =====//
    //==================//
    double return_value(const dbl_vec &input) const override final
    {
        UNUSED(input);
        return get_arg(0) * get_arg(1);
    }
};
//--- Division Function. ---//
class Division_Function : public Function
{
public:
    //===== CONSTRUCTOR =====//
    Division_Function()
    {
        init_defaults("div",
                      {"x", "y"}, {0., 0.},
                      HOUSEPI_NA_DOUBLE, HOUSEPI_NA_DOUBLE,
                      false, false);
    }
    ~Division_Function() = default;

    //==================//
    //===== GETTER =====//
    //==================//
    double return_value(const dbl_vec &input) const override final
    {
        UNUSED(input);
        return get_arg(0) / get_arg(1);
    }
};
//--- Power Function. ---//
class Power_Function : public Function
{
public:
    //===== CONSTRUCTOR =====//
    Power_Function()
    {
        init_defaults("pow",
                      {"x", "y"}, {0., 0.},
                      HOUSEPI_NA_DOUBLE, HOUSEPI_NA_DOUBLE,
                      false, false);
    }
    ~Power_Function() = default;

    //==================//
    //===== GETTER =====//
    //==================//
    double return_value(const dbl_vec &input) const override final
    {
        UNUSED(input);
        return std::pow(get_arg(0), get_arg(1));
    }
};

//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//

/////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////  Function Switchboard  /////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

/* Auxiliary hidden template to return the chosen distribution. */
std::shared_ptr<Function> function_switch(const std::string &name)
{
    const str_vec &options = {"eq", "add", "sub", "mult", "div", "pow"};

    if (tools::is_na(name))
        return std::make_shared<Function>();

    else if (name == "eq")
        return std::make_shared<Equality_Function>();

    else if (name == "add")
        return std::make_shared<Addition_Function>();

    else if (name == "sub")
        return std::make_shared<Subtraction_Function>();

    else if (name == "mult")
        return std::make_shared<Multiplication_Function>();

    else if (name == "div")
        return std::make_shared<Division_Function>();

    else if (name == "pow")
        return std::make_shared<Power_Function>();

    else
        throw std::invalid_argument(
            "Unknown function \"" + name + "\"." + tools::find_closest(name, options));
}
