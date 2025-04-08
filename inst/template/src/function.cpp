/*
 * function.cpp
 *
 * Copyright (C) 2025 Francesco Camaglia, Institut Pasteur
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

inline double Argument::get() const { return *(ptr_to_value); }
inline const double *Argument::get_ptr() const { return ptr_to_value; }
bool Argument::is_internal_ptr() const { return ptr_to_value == &value; }
const std::string Argument::show() const { return (tools::natostring(get())); }

/////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////  Argument : setter ///////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

inline void Argument::set(double x) { (this->*_func_set_)(x); }
void Argument::assign_ptr(const double *y)
{
    _func_set_ = &Argument::_set_arg_pointer_;
    ptr_to_value = y;
}
void Argument::restore_internal_ptr()
{
    _func_set_ = &Argument::_set_arg_default_;
    ptr_to_value = &value;
}

/////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////  Argument : functional /////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

void Argument::_set_arg_default_(double x) { value = x; }

/* note : the set_arg is only used for the adaptive
 * for this reason a pointer argument can't be set, but just assigned
 */
void Argument::_set_arg_pointer_(double x)
{
    UNUSED(x);
    throw std::runtime_error(
        "Pointer argument \"" + name + "\" can not be changed.");
}

//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////  Base Function : Constructors  ////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

Function::Function() {};
Function::~Function() = default;

/////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////  Base Function : Getter  ///////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

/* Auxiliary function to retrieve the index of a given argument by name. */
inline ushort Function::find(const std::string &arg_name) const
{
    for (ushort c = 0; c < n_args; ++c)
        if (arguments[c].name == arg_name)
            return c;
    return n_args;
}

/* Assert argument exists at assignment. */
ushort Function::assert_arg(const std::string &arg_name) const
{
    ushort i = this->find(arg_name);
    if (i == n_args)
    {
        str_vec options = {};
        for (ushort c = 0; c < n_args; ++c)
            options.push_back(arguments[c].name);
        throw std::invalid_argument(
            "Unknown argument \"" + arg_name + "\" for \"" + name + "\"." + tools::find_closest(arg_name, options));
    }
    return i;
}

double Function::get_arg(ushort i) const { return arguments[i].get(); }
double Function::get_min() const { return range_min.get(); }
double Function::get_max() const { return range_max.get(); }

const double *Function::get_ptr_to_arg(ushort i) const { return arguments[i].get_ptr(); }
const double *Function::get_ptr_to_min() const { return range_min.get_ptr(); }
const double *Function::get_ptr_to_max() const { return range_max.get_ptr(); }

/* Returns TRUE if value x is below range minimum. */
bool Function::below(double x) const { return (this->*_func_below_)(x); }
/* Returns TRUE if value x is above range maximum. */
bool Function::above(double x) const { return (this->*_func_above_)(x); }
/* Returns TRUE if value x is within range. */
bool Function::valid(double x) const { return (this->*_func_valid_)(x); }

/////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////  Base Function : set range  ////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

ushort Function::which_valid_case() const
{
    // bounded if not pointing internally to NA
    bool is_lower_bound = !(range_min.is_internal_ptr() && tools::is_na(range_min.get()));
    bool is_upper_bound = !(range_max.is_internal_ptr() && tools::is_na(range_max.get()));
    return is_lower_bound + 2 * is_upper_bound;
}

/* Set the range of validity of the value, and the corresponding validation function. */
void Function::init_functional_range(bool is_min_open, bool is_max_open)
{
    switch (which_valid_case())
    {
    case 0: //----- unbounded -----//
        _func_below_ = &Function::_below_never_;
        _func_above_ = &Function::_above_never_;
        _func_valid_ = &Function::_valid_always_;
        break;

    case 1: //----- lower bound -----//
        _func_above_ = &Function::_above_never_;
        if (is_min_open)
        {
            _func_below_ = &Function::_below_open_min_;
            _func_valid_ = &Function::_valid_open_min_;
        }
        else
        {
            _func_below_ = &Function::_below_close_min_;
            _func_valid_ = &Function::_valid_close_min_;
        }
        break;

    case 2: //----- upper bound -----//
        _func_below_ = &Function::_below_never_;
        if (is_max_open)
        {
            _func_above_ = &Function::_above_open_max_;
            _func_valid_ = &Function::_valid_open_max_;
        }
        else
        {
            _func_above_ = &Function::_above_close_max_;
            _func_valid_ = &Function::_valid_close_max_;
        }
        break;

    case 3: //----- within an interval -----//

        if (range_min.get() >= range_max.get())
            throw std::invalid_argument(
                "Inconsistent range " + show_range() + ".");

        if ((is_min_open) && (is_max_open))
        {
            _func_below_ = &Function::_below_open_min_;
            _func_above_ = &Function::_above_open_max_;
            _func_valid_ = &Function::_valid_both_open_;
        }
        else if ((is_min_open) && (!is_max_open))
        {
            _func_below_ = &Function::_below_open_min_;
            _func_above_ = &Function::_above_close_max_;
            _func_valid_ = &Function::_valid_open_min_close_max_;
        }
        else if ((!is_min_open) && (is_max_open))
        {
            _func_below_ = &Function::_below_close_min_;
            _func_above_ = &Function::_above_open_max_;
            _func_valid_ = &Function::_valid_close_min_open_max_;
        }
        else
        {
            _func_below_ = &Function::_below_close_min_;
            _func_above_ = &Function::_above_close_max_;
            _func_valid_ = &Function::_valid_both_close_;
        }
        break;
    }
}

const std::string Function::show_min() const { return range_min.show(); }
const std::string Function::show_max() const { return range_max.show(); }

/* Return formatted range. */
const std::string Function::show_range() const
{
    if (_func_valid_ == &Function::_valid_always_)
        return "";

    else if (_func_valid_ == &Function::_valid_open_min_)
        return "(" + show_min() + ", infty)";

    else if (_func_valid_ == &Function::_valid_close_min_)
        return "[" + show_min() + ", infty)";

    else if (_func_valid_ == &Function::_valid_open_max_)
        return "(-infty, " + show_max() + ")";

    else if (_func_valid_ == &Function::_valid_close_max_)
        return "(-infty, " + show_max() + "]";

    else if (_func_valid_ == &Function::_valid_both_open_)
        return "(" + show_min() + ", " + show_max() + ")";

    else if (_func_valid_ == &Function::_valid_open_min_close_max_)
        return "(" + show_min() + ", " + show_max() + "]";

    else if (_func_valid_ == &Function::_valid_close_min_open_max_)
        return "[" + show_min() + ", " + show_max() + ")";

    else // (_func_valid_ == &Function::_valid_both_close_)
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
void Function::set_arg(ushort i, double x) { arguments[i].set(x); }

/* Assignment of non default POINTER arguments. */
void Function::assign_scalar_arg(const std::string &arg_name, double x)
{
    ushort i = this->assert_arg(arg_name);
    arguments[i].set(x);
    arguments[i].restore_internal_ptr();
}

/* Assignment of non default POINTER arguments. */
void Function::assign_pointer_arg(const std::string &arg_name, const double *y)
{
    ushort i = this->assert_arg(arg_name);
    arguments[i].assign_ptr(y);
}

/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////  Base Function : Range  ///////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

/* Assign scalar range minimum. */
void Function::assign_scalar_range_min(double x)
{
    if (name == "unif")
    {
        // Ensure uniform range min == argument "min"
        range_min.assign_ptr(arguments[this->find("min")].get_ptr());
        if (!tools::is_na(x))
        {
            tools::warning(
                "Depracated: overwriting \"unif\" argument min with range min.");
            assign_scalar_arg("min", x);
        }
    }
    else
    {
        range_min.set(x);
        range_min.restore_internal_ptr();
    }
}

/* Assign pointer range minimum. */
void Function::assign_pointer_range_min(const double *y)
{
    // Ensure uniform range min == argument "min"
    if (name == "unif")
    {
        // Ensure uniform range min == argument "min"
        range_min.assign_ptr(arguments[this->find("min")].get_ptr());
        if (y != nullptr)
        {
            tools::warning(
                "Depracated: overwriting \"unif\" argument min with range min.");
            assign_pointer_arg("min", y);
        }
    }
    else
    {
        range_min.assign_ptr(y);
    }
}

/* Assign scalar range maximum. */
void Function::assign_scalar_range_max(double x)
{
    if (name == "unif")
    {
        // Ensure uniform range max == argument "max"
        range_max.assign_ptr(arguments[this->find("max")].get_ptr());
        if (!tools::is_na(x))
        {
            tools::warning(
                "Depracated: overwriting \"unif\" argument max with range max.");
            assign_scalar_arg("max", x);
        }
    }
    else
    {
        range_max.set(x);
        range_max.restore_internal_ptr();
    }
}

/* Assign pointer range maximum. */
void Function::assign_pointer_range_max(const double *y)
{
    // Ensure uniform range max == argument "max"
    if (name == "unif")
    {
        // Ensure uniform range max == argument "max"
        range_max.assign_ptr(arguments[this->find("max")].get_ptr());
        if (y != nullptr)
        {
            tools::warning(
                "Depracated: overwriting \"unif\" argument max with range max.");
            assign_pointer_arg("max", y);
        }
    }
    else
    {
        range_max.assign_ptr(y);
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
const std::string Function::show_name() const { return name; }

/* Return formatted argument names of the distribution. */
const std::string Function::show_arg_names() const
{
    if (n_args > 0)
    {
        std::string output = arguments[0].name;
        for (ushort c = 1; c < n_args; ++c)
            output += C_LIST_DELIM + arguments[c].name;
        return output;
    }

    return NULL_STRING;
}

/* Return formatted argument names of the distribution. */
const std::string Function::show_arg_values() const
{
    if (n_args > 0)
    {
        std::string output = arguments[0].show();
        for (ushort c = 1; c < n_args; ++c)
            output += C_LIST_DELIM + arguments[c].show();
        return output;
    }

    return NULL_STRING;
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
    this->name = name;

    // number of arguments
    if (arg_name.size() != arg_value.size())
        throw std::length_error(
            "Internal inconsistency for \"" + name + "\" arguments.");
    n_args = arg_name.size();

    // assign default arguments
    arguments = new Argument[arg_name.size()];
    for (ushort c = 0; c < n_args; ++c)
    {
        arguments[c].name = arg_name[c];
        arguments[c].set(arg_value[c]);
    }

    //----- assign domain -----//
    // NOTE: at the moment only place where to set open ranges
    // NOTE: not sure if setting pointer to smtg exiting the domain I still have control
    assign_scalar_range_min(domain_min);
    assign_scalar_range_max(domain_max);
    init_functional_range(is_min_open, is_max_open);
}
