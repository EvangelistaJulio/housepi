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

class Argument
{
public:
    std::string name = NULL_STRING;

    //===== GETTER =====//
    inline double get() const { return *(ptr_to_value); }
    inline const double *get_ptr() const { return ptr_to_value; }
    bool is_internal_ptr() const { return ptr_to_value == &value; }
    const std::string show() const { return (tools::natostring(get())); }

    //===== SETTER =====//
    inline void set(double x) { (this->*_func_set_)(x); }
    void assign_ptr(const double *y)
    {
        _func_set_ = &Argument::_set_arg_pointer_;
        ptr_to_value = y;
    }
    void restore_internal_ptr()
    {
        _func_set_ = &Argument::_set_arg_default_;
        ptr_to_value = &value;
    }

private:
    double value = NULL_DOUBLE;
    const double *ptr_to_value = &value;

    // note : the set_arg is only used for the adaptive
    // for this reason a pointer argument can't be set, but just assigned
    void (Argument::*_func_set_)(double) = &Argument::_set_arg_default_;
    void _set_arg_default_(double x) { value = x; }
    void _set_arg_pointer_(double x)
    {
        UNUSED(x);
        throw std::runtime_error(
            "Pointer argument \"" + name + "\" can not be changed.");
    }
};

//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//

class Function::Impl
{
public:
    std::string name = NULL_STRING;

    //===== Arguments =====//

    Argument *arguments = nullptr;
    ushort n_args = 0;

    //----- Domain -----//
    Argument range_min = Argument();
    Argument range_max = Argument();

    //===== Non-default arguments =====//
    ushort find(const std::string &arg_name) const;

    //===== Range =====//
    void init_functional_range(bool is_min_open, bool is_max_open);

    //===== Getter =====//
    bool valid(double x) const { return (this->*_func_valid_range_)(x); }
    ushort which_valid_case() const;

    //===== SHOW =====//
    const std::string show_min() const;
    const std::string show_max() const;
    const std::string show_range() const;

private:
    //===== VALID =====//
    // note: for the moment range extremes can't be one open and other close

    //---- alternative valid range methods -----//
    bool (Impl::*_func_valid_range_)(double) const = &Impl::_valid_range_always_;
    bool _valid_range_always_(double x) const
    {
        UNUSED(x);
        return true;
    }

    //----- range : minimum -----//
    bool _valid_range_open_min_(double x) const { return (x > range_min.get()); }
    bool _valid_range_close_min_(double x) const { return (x >= range_min.get()); }

    //----- range : maximum -----//
    bool _valid_range_open_max_(double x) const { return (x < range_max.get()); }
    bool _valid_range_close_max_(double x) const { return (x <= range_max.get()); }

    //----- range : both -----//
    bool _valid_range_both_open_(double x) const { return (x > range_min.get()) && (x < range_max.get()); }
    bool _valid_range_open_min_close_max_(double x) const { return (x > range_min.get()) && (x <= range_max.get()); }
    bool _valid_range_close_min_open_max_(double x) const { return (x >= range_min.get()) && (x < range_max.get()); }
    bool _valid_range_both_close_(double x) const { return (x >= range_min.get()) && (x <= range_max.get()); }
};

/////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////  Base Function : methods  //////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

inline ushort Function::Impl::find(const std::string &arg_name) const
{
    for (ushort c = 0; c < n_args; ++c)
        if (arguments[c].name == arg_name)
            return c;
    return n_args;
}

/////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////  Base Function : set range  ////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

ushort Function::Impl::which_valid_case() const
{
    // bounded if not pointing internally to NA
    bool is_lower_bound = !(range_min.is_internal_ptr() && tools::is_na(range_min.get()));
    bool is_upper_bound = !(range_max.is_internal_ptr() && tools::is_na(range_max.get()));
    return is_lower_bound + 2 * is_upper_bound;
}

/* Set the range of validity of the value, and the corresponding validation function. */
void Function::Impl::init_functional_range(bool is_min_open, bool is_max_open)
{
    switch (which_valid_case())
    {
    case 0: //----- unbounded -----//
        _func_valid_range_ = &Impl::_valid_range_always_;
        break;

    case 1: //----- lower bound -----//
        if (is_min_open)
            _func_valid_range_ = &Impl::_valid_range_open_min_;
        else
            _func_valid_range_ = &Impl::_valid_range_close_min_;
        break;

    case 2: //----- upper bound -----//
        if (is_max_open)
            _func_valid_range_ = &Impl::_valid_range_open_max_;
        else
            _func_valid_range_ = &Impl::_valid_range_close_max_;
        break;

    case 3: //----- within an interval -----//
        if ((is_min_open) && (is_max_open))
            _func_valid_range_ = &Impl::_valid_range_both_open_;

        else if ((is_min_open) && (!is_max_open))
            _func_valid_range_ = &Impl::_valid_range_open_min_close_max_;

        else if ((!is_min_open) && (is_max_open))
            _func_valid_range_ = &Impl::_valid_range_close_min_open_max_;

        else
            _func_valid_range_ = &Impl::_valid_range_both_close_;

        if (range_min.get() >= range_max.get())
            throw std::invalid_argument(
                "Inconsistent range " + show_range() + ".");

        break;
    }
}

const std::string Function::Impl::show_min() const { return range_min.show(); }
const std::string Function::Impl::show_max() const { return range_max.show(); }

/* Return formatted maximum range. */
const std::string Function::Impl::show_range() const
{
    if (_func_valid_range_ == &Impl::_valid_range_always_)
        return "";

    else if (_func_valid_range_ == &Impl::_valid_range_open_min_)
        return "(" + show_min() + ", infty)";

    else if (_func_valid_range_ == &Impl::_valid_range_close_min_)
        return "[" + show_min() + ", infty)";

    else if (_func_valid_range_ == &Impl::_valid_range_open_max_)
        return "(-infty, " + show_max() + ")";

    else if (_func_valid_range_ == &Impl::_valid_range_close_max_)
        return "(-infty, " + show_max() + "]";

    else if (_func_valid_range_ == &Impl::_valid_range_both_open_)
        return "(" + show_min() + ", " + show_max() + ")";

    else if (_func_valid_range_ == &Impl::_valid_range_open_min_close_max_)
        return "(" + show_min() + ", " + show_max() + "]";

    else if (_func_valid_range_ == &Impl::_valid_range_close_min_open_max_)
        return "[" + show_min() + ", " + show_max() + ")";

    else if (_func_valid_range_ == &Impl::_valid_range_both_close_)
        return "[" + show_min() + ", " + show_max() + "]";
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

Function::Function() : pImpl(std::make_unique<Impl>()) {};
Function::~Function() = default;

/////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////  Base Function : Getter  ///////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

ushort Function::find(const std::string &arg_name) const { return pImpl->find(arg_name); }
double Function::get_arg(ushort i) const { return pImpl->arguments[i].get(); }
double Function::get_min() const { return pImpl->range_min.get(); }
double Function::get_max() const { return pImpl->range_max.get(); }
const double *Function::get_ptr_to_arg(ushort i) const { return pImpl->arguments[i].get_ptr(); }
const double *Function::get_ptr_to_min() const { return pImpl->range_min.get_ptr(); }
const double *Function::get_ptr_to_max() const { return pImpl->range_max.get_ptr(); }
bool Function::valid(double x) const { return pImpl->valid(x); }

/////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////  Base Function : Arguments  /////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

/* Set argument value. */
void Function::set_arg(ushort i, double x) { pImpl->arguments[i].set(x); }

/* Assignment of non default POINTER arguments. */
void Function::assign_scalar_arg(const std::string &arg_name, double x)
{
    ushort i = this->find(arg_name);
    if (i == pImpl->n_args)
        throw std::invalid_argument(
            "Unknown argument \"" + arg_name + "\" for \"" + pImpl->name + "\".");
    else
    {
        pImpl->arguments[i].set(x);
        pImpl->arguments[i].restore_internal_ptr();
    }
}

/* Assignment of non default POINTER arguments. */
void Function::assign_pointer_arg(const std::string &arg_name, const double *y)
{
    ushort i = this->find(arg_name);
    if (i == pImpl->n_args)
        throw std::invalid_argument(
            "Unknown argument \"" + arg_name + "\" for \"" + pImpl->name + "\".");
    else
        pImpl->arguments[i].assign_ptr(y);
}

/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////  Base Function : Range  ///////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

/* Assign scalar range minimum. */
void Function::assign_scalar_range_min(double x)
{
    if (pImpl->name == "unif")
    {
        // Ensure uniform range min == argument "min"
        pImpl->range_min.assign_ptr(pImpl->arguments[this->find("min")].get_ptr());
        if (!tools::is_na(x))
        {
            tools::warning(
                "Depracated: overwriting \"unif\" argument min with range min.");
            assign_scalar_arg("min", x);
        }
    }
    else
    {
        pImpl->range_min.set(x);
        pImpl->range_min.restore_internal_ptr();
    }
}

/* Assign pointer range minimum. */
void Function::assign_pointer_range_min(const double *y)
{
    // Ensure uniform range min == argument "min"
    if (pImpl->name == "unif")
    {
        // Ensure uniform range min == argument "min"
        pImpl->range_min.assign_ptr(pImpl->arguments[this->find("min")].get_ptr());
        if (y != nullptr)
        {
            tools::warning(
                "Depracated: overwriting \"unif\" argument min with range min.");
            assign_pointer_arg("min", y);
        }
    }
    else
    {
        pImpl->range_min.assign_ptr(y);
    }
}

/* Assign scalar range maximum. */
void Function::assign_scalar_range_max(double x)
{
    if (pImpl->name == "unif")
    {
        // Ensure uniform range max == argument "max"
        pImpl->range_max.assign_ptr(pImpl->arguments[this->find("max")].get_ptr());
        if (!tools::is_na(x))
        {
            tools::warning(
                "Depracated: overwriting \"unif\" argument max with range max.");
            assign_scalar_arg("max", x);
        }
    }
    else
    {
        pImpl->range_max.set(x);
        pImpl->range_max.restore_internal_ptr();
    }
}

/* Assign pointer range maximum. */
void Function::assign_pointer_range_max(const double *y)
{
    // Ensure uniform range max == argument "max"
    if (pImpl->name == "unif")
    {
        // Ensure uniform range max == argument "max"
        pImpl->range_max.assign_ptr(pImpl->arguments[this->find("max")].get_ptr());
        if (y != nullptr)
        {
            tools::warning(
                "Depracated: overwriting \"unif\" argument max with range max.");
            assign_pointer_arg("max", y);
        }
    }
    else
    {
        pImpl->range_max.assign_ptr(y);
    }
}

/* Initialize open range functional methods. */
void Function::init_functional_range(bool is_min_open, bool is_max_open) { pImpl->init_functional_range(is_min_open, is_max_open); }

/* Return which valid case. */
ushort Function::which_valid_case() const { return pImpl->which_valid_case(); }

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
const std::string Function::show_name() const { return pImpl->name; }

/* Return formatted argument names of the distribution. */
const std::string Function::show_arg_names() const
{
    if (pImpl->n_args > 0)
    {
        std::string output = pImpl->arguments[0].name;
        for (ushort c = 1; c < pImpl->n_args; ++c)
            output += C_LIST_DELIM + pImpl->arguments[c].name;
        return output;
    }

    return NULL_STRING;
}

/* Return formatted argument names of the distribution. */
const std::string Function::show_arg_values() const
{
    if (pImpl->n_args > 0)
    {
        std::string output = pImpl->arguments[0].show();
        for (ushort c = 1; c < pImpl->n_args; ++c)
            output += C_LIST_DELIM + pImpl->arguments[c].show();
        return output;
    }

    return NULL_STRING;
}

/* Return formatted minimum value. */
const std::string Function::show_min() const { return pImpl->show_min(); }

/* Return formatted maximum value. */
const std::string Function::show_max() const { return pImpl->show_max(); }

/* Return formatted maximum range. */
const std::string Function::show_range() const { return pImpl->show_range(); }

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
    pImpl->name = name;

    // number of arguments
    if (arg_name.size() != arg_value.size())
        throw std::length_error(
            "Internal inconsistency for \"" + name + "\" arguments.");
    pImpl->n_args = arg_name.size();

    // assign default arguments
    pImpl->arguments = new Argument[arg_name.size()];
    for (ushort c = 0; c < pImpl->n_args; ++c)
    {
        pImpl->arguments[c].name = arg_name[c];
        pImpl->arguments[c].set(arg_value[c]);
    }

    //----- assign domain -----//
    // NOTE: at the moment only place where to set open ranges
    // NOTE: not sure if setting pointer to smtg exiting the domain I still have control
    assign_scalar_range_min(domain_min);
    assign_scalar_range_max(domain_max);
    pImpl->init_functional_range(is_min_open, is_max_open);
}
