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
//////////////////////////////////////  Argument : getter ///////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

inline double Argument::get() const { return *(_ptr_to_value); }
inline const double *Argument::get_ptr() const { return _ptr_to_value; }
bool Argument::is_internal_ptr() const { return _ptr_to_value == &_value; }
const std::string Argument::show() const { return tools::natostring(get()); }
const std::string Argument::show_name() const { return _argument_name; }

/////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////  Argument : setter ///////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

inline void Argument::set_name(const std::string &name) { _argument_name = name; }
inline void Argument::set(double x) { (this->*_do_set_arg)(x); }
void Argument::assign_ptr(const double *y)
{
    _do_set_arg = &Argument::__set_arg_pointer;
    _ptr_to_value = y;
}
void Argument::restore_internal_ptr()
{
    _do_set_arg = &Argument::__set_arg_default;
    _ptr_to_value = &_value;
}

/////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////  Argument : functional /////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

void Argument::__set_arg_default(double x) { _value = x; }

/* note : the set_arg is only used for the adaptive
 * for this reason a pointer argument can't be set, but just assigned
 */
void Argument::__set_arg_pointer(double x)
{
    UNUSED(x);
    throw std::runtime_error(
        "Pointer argument \"" + _argument_name + "\" can not be changed.");
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

/* Auxiliary function to retrieve the index of a given argument by name. */
inline ushort Function::find(const std::string &arg_name) const
{
    for (ushort c = 0; c < _n_args; ++c)
        if (_arguments[c].show_name() == arg_name)
            return c;
    return _n_args;
}

/* Assert argument exists at assignment. */
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

const double *Function::get_ptr_to_arg(ushort i) const { return _arguments[i].get_ptr(); }
const double *Function::get_ptr_to_min() const { return _range_min.get_ptr(); }
const double *Function::get_ptr_to_max() const { return _range_max.get_ptr(); }

/* Returns TRUE if value x is below range minimum. */
bool Function::below(double x) const { return (this->*_do_below)(x); }
/* Returns TRUE if value x is above range maximum. */
bool Function::above(double x) const { return (this->*_do_above)(x); }
/* Returns TRUE if value x is within range. */
bool Function::valid(double x) const { return (this->*_do_valid)(x); }

/////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////  Base Function : set range  ////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

ushort Function::which_valid_case() const
{
    // bounded if not pointing internally to NA
    bool is_lower_bound = !(_range_min.is_internal_ptr() && tools::is_na(_range_min.get()));
    bool is_upper_bound = !(_range_max.is_internal_ptr() && tools::is_na(_range_max.get()));
    return is_lower_bound + 2 * is_upper_bound;
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
        return "(" + show_min() + ", infty)";

    else if (_do_valid == &Function::__valid_close_min)
        return "[" + show_min() + ", infty)";

    else if (_do_valid == &Function::__valid_open_max)
        return "(-infty, " + show_max() + ")";

    else if (_do_valid == &Function::__valid_close_max)
        return "(-infty, " + show_max() + "]";

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
void Function::assign_scalar(const std::string &arg_name, bool is_arg, double value)
{
    if (is_arg)
        _assign_scalar_arg(arg_name, value);
    else
        _assign_scalar_range(arg_name, value);
}
//--- Assign pointer argument or range value. ---//
void Function::assign_pointer(const std::string &arg_name, bool is_arg, const double *ptr_to_value)
{
    if (is_arg)
        _assign_pointer_arg(arg_name, ptr_to_value);
    else
        _assign_pointer_range(arg_name, ptr_to_value);
}

/* Assignment of non default POINTER arguments. */
void Function::_assign_scalar_arg(const std::string &arg_name, double x)
{
    ushort i = this->assert_arg(arg_name);
    _arguments[i].set(x);
    _arguments[i].restore_internal_ptr();
}

/* Assignment of non default POINTER arguments. */
void Function::_assign_pointer_arg(const std::string &arg_name, const double *y)
{
    ushort i = this->assert_arg(arg_name);
    _arguments[i].assign_ptr(y);
}

/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////  Base Function : Range  ///////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

/* Assign scalar range. */
void Function::_assign_scalar_range(const std::string &range_name, double x)
{
    if (_function_name == "unif")
    {
        // Ensure uniform range (min, max) match  == argument (min,max)
        // UPDATE: this is a necessary step and should be moved elsewhere
        if (range_name == "min")
            _range_min.assign_ptr(_arguments[this->find("min")].get_ptr());
        else if (range_name == "max")
            _range_max.assign_ptr(_arguments[this->find("max")].get_ptr());
        else
            throw std::invalid_argument(
                "Unknown range name \"" + range_name + "\" for \"" + _function_name + "\".");

        if (!tools::is_na(x))
        {
            tools::warning(
                "Depracated: overwriting \"unif\" argument \"" + range_name + "\" with range value.");
            _assign_scalar_arg(range_name, x);
        }
    }
    else
    {
        if (range_name == "min")
        {
            _range_min.set(x);
            _range_min.restore_internal_ptr();
        }
        else if (range_name == "max")
        {
            _range_max.set(x);
            _range_max.restore_internal_ptr();
        }
        else
            throw std::invalid_argument(
                "Unknown range name \"" + range_name + "\" for \"" + _function_name + "\".");
    }
}

/* Assign pointer range. */
void Function::_assign_pointer_range(const std::string &range_name, const double *y)
{
    if (_function_name == "unif")
    {
        // Ensure uniform range (min, max) match  == argument (min,max)
        // UPDATE: this is a necessary step and should be moved elsewhere
        if (range_name == "min")
            _range_min.assign_ptr(_arguments[this->find("min")].get_ptr());
        else if (range_name == "max")
            _range_max.assign_ptr(_arguments[this->find("max")].get_ptr());
        else
            throw std::invalid_argument(
                "Unknown range name \"" + range_name + "\" for \"" + _function_name + "\".");

        if (y != nullptr)
        {
            tools::warning(
                "Depracated: overwriting \"unif\" argument \"" + range_name + "\" with range value.");
            _assign_pointer_arg(range_name, y);
        }
    }
    else
    {
        if (range_name == "min")
            _range_min.assign_ptr(y);
        else if (range_name == "max")
            _range_max.assign_ptr(y);
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
    _arguments = new Argument[arg_name.size()];
    for (ushort c = 0; c < _n_args; ++c)
    {
        _arguments[c].set_name(arg_name[c]);
        _arguments[c].set(arg_value[c]);
    }

    //----- assign domain -----//
    // NOTE: at the moment only place where to set open ranges
    // NOTE: not sure if setting pointer to smtg exiting the domain I still have control
    _assign_scalar_range("min", domain_min);
    _assign_scalar_range("max", domain_max);
    init_functional_range(is_min_open, is_max_open);
}

//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//

/* Identity Function. */
class Identity_Function : public Function
{
public:
    //===== CONSTRUCTOR =====//
    Identity_Function()
    {
        init_defaults("identity",
                      {"x"}, {0.},
                      HOUSEPI_NA_DOUBLE, HOUSEPI_NA_DOUBLE,
                      false, false);
    }
    ~Identity_Function() = default;

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
    const str_vec &options = {"identity", "add", "sub"};

    if (tools::is_na(name))
        return std::make_shared<Function>();

    else if (name == "add")
        return std::make_shared<Addition_Function>();

    else if (name == "sub")
        return std::make_shared<Subtraction_Function>();

    else if (name == "identity")
        return std::make_shared<Identity_Function>();

    else
        throw std::invalid_argument(
            "Unknown function \"" + name + "\"." + tools::find_closest(name, options));
}
