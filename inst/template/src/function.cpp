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

class Function_Argument
{
public:
    std::string name = NULL_STRING;

    //===== GETTER =====//
    inline double get() const { return *(ptr_to_value); }
    inline const double *get_ptr() const { return &value; }
    bool is_internal_ptr() const { return ptr_to_value == get_ptr(); }

    //===== SETTER =====//
    inline void set(double x) { (this->*_func_set_)(x); }
    void assign_ptr(const double *y)
    {
        _func_set_ = &Function_Argument::_set_arg_pointer_;
        ptr_to_value = y;
    }
    void restore_internal_ptr()
    {
        _func_set_ = &Function_Argument::_set_arg_default_;
        ptr_to_value = &value;
    }

private:
    double value = NULL_DOUBLE;
    const double *ptr_to_value = &value;

    // note : the set_arg is only used for the adaptive
    // for this reason a pointer argument can't be set, but just assigned
    void (Function_Argument::*_func_set_)(double) = &Function_Argument::_set_arg_default_;
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

class Base_Function::Impl
{
public:
    std::string name = NULL_STRING;

    //===== Arguments =====//

    Function_Argument *arguments = nullptr;
    ushort n_args = 0;

    //----- Domain -----//
    Function_Argument range_min = Function_Argument();
    Function_Argument range_max = Function_Argument();

    //===== Non-default arguments =====//
    ushort find(const std::string &arg_name);

    //===== Range =====//
    void init_range(bool is_open);

    //===== Getter =====//
    bool valid(double x) const { return (this->*_func_valid_range_)(x); }
    ushort which_valid_case() const;

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

    //----- close range -----//
    bool _valid_range_close_min_(double x) const { return (x >= range_min.get()); }
    bool _valid_range_close_max_(double x) const { return (x <= range_max.get()); }
    bool _valid_range_close_both_(double x) const { return (x >= range_min.get()) && (x <= range_max.get()); }

    //----- open range -----//
    bool _valid_range_open_min_(double x) const { return (x > range_min.get()); }
    bool _valid_range_open_max_(double x) const { return (x < range_max.get()); }
    bool _valid_range_open_both_(double x) const { return (x > range_min.get()) && (x < range_max.get()); }
};

/////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////  Base Function : methods  //////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

inline ushort Base_Function::Impl::find(const std::string &arg_name)
{
    for (ushort c = 0; c < n_args; ++c)
        if (arguments[c].name == arg_name)
            return c;
    return n_args;
}

/////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////  Base Function : set range  ////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

ushort Base_Function::Impl::which_valid_case() const
{
    // bounded if not pointing internally to NA
    bool is_lower_bound = !(range_min.is_internal_ptr() && tools::is_na(range_min.get()));
    bool is_upper_bound = !(range_max.is_internal_ptr() && tools::is_na(range_max.get()));
    return is_lower_bound + 2 * is_upper_bound;
}

/* Set the range of validity of the value, and the corresponding (close) validation function. */
void Base_Function::Impl::init_range(bool is_open)
{
    switch (which_valid_case())
    {
    case 0:
        _func_valid_range_ = &Impl::_valid_range_always_;
        break;
    case 1:
        if (is_open)
            _func_valid_range_ = &Impl::_valid_range_open_min_;
        else
            _func_valid_range_ = &Impl::_valid_range_close_min_;
        break;
    case 2:
        if (is_open)
            _func_valid_range_ = &Impl::_valid_range_open_max_;
        else
            _func_valid_range_ = &Impl::_valid_range_close_max_;
        break;
    case 3:
        if (range_min.get() >= range_max.get())
        {
            std::ostringstream oss;
            oss << "Inconsistent range {min:" << range_min.get()
                << ", max:" << range_max.get() << "}.";
            throw std::invalid_argument(oss.str());
        }
        if (is_open)
            _func_valid_range_ = &Impl::_valid_range_open_both_;
        else
            _func_valid_range_ = &Impl::_valid_range_close_both_;
        break;
    }
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

Base_Function::Base_Function() : pImpl(std::make_unique<Impl>()) {};
Base_Function::~Base_Function() = default;

/////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////  Base Function : Getter  ///////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

const std::string &Base_Function::get_name() const { return pImpl->name; }
double Base_Function::get_arg(ushort i) const { return pImpl->arguments[i].get(); }
double Base_Function::get_range_min() const { return pImpl->range_min.get(); }
double Base_Function::get_range_max() const { return pImpl->range_max.get(); }
bool Base_Function::valid(double x) const { return pImpl->valid(x); }

/////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////  Base Function : Arguments  /////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

/* Set argument value. */
void Base_Function::set_arg(ushort i, double x) { pImpl->arguments[i].set(x); }

/* Assignment of non default POINTER arguments. */
void Base_Function::assign_scalar_arg(const std::string &arg_name, double x)
{
    ushort it = pImpl->find(arg_name);
    if (it == pImpl->n_args)
        throw std::invalid_argument(
            "Unknown argument \"" + arg_name + "\" for \"" + pImpl->name + "\".");
    else
    {
        pImpl->arguments[it].set(x);
        pImpl->arguments[it].restore_internal_ptr();
    }
}

/* Assignment of non default POINTER arguments. */
void Base_Function::assign_pointer_arg(const std::string &arg_name, const double *y)
{
    ushort it = pImpl->find(arg_name);
    if (it == pImpl->n_args)
        throw std::invalid_argument(
            "Unknown argument \"" + arg_name + "\" for \"" + pImpl->name + "\".");
    else
        pImpl->arguments[it].assign_ptr(y);
}

/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////  Base Function : Range  ///////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

/* Assign scalar range minimum. */
void Base_Function::assign_scalar_range_min(double x)
{
    if (pImpl->name == "unif")
    {
        // Ensure uniform range min == argument "min"
        pImpl->range_min.assign_ptr(pImpl->arguments[pImpl->find("min")].get_ptr());
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
void Base_Function::assign_pointer_range_min(const double *y)
{
    // Ensure uniform range min == argument "min"
    if (pImpl->name == "unif")
    {
        // Ensure uniform range min == argument "min"
        pImpl->range_min.assign_ptr(pImpl->arguments[pImpl->find("min")].get_ptr());
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
void Base_Function::assign_scalar_range_max(double x)
{
    if (pImpl->name == "unif")
    {
        // Ensure uniform range max == argument "max"
        pImpl->range_max.assign_ptr(pImpl->arguments[pImpl->find("max")].get_ptr());
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
void Base_Function::assign_pointer_range_max(const double *y)
{
    // Ensure uniform range max == argument "max"
    if (pImpl->name == "unif")
    {
        // Ensure uniform range max == argument "max"
        pImpl->range_max.assign_ptr(pImpl->arguments[pImpl->find("max")].get_ptr());
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
void Base_Function::init_open_range() { pImpl->init_range(true); }

/* Initialize open range functional methods. */
void Base_Function::init_close_range() { pImpl->init_range(false); }

/* Return which valid case. */
ushort Base_Function::which_valid_case() const { return pImpl->which_valid_case(); }

//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//

/////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////  Base Function : Show  /////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

/* Return formatted argument names of the distribution. */
const std::string Base_Function::show_arg_names() const
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
const std::string Base_Function::show_arg_values() const
{
    if (pImpl->n_args > 0)
    {
        std::string output = std::to_string(pImpl->arguments[0].get());
        for (ushort c = 1; c < pImpl->n_args; ++c)
            output += C_LIST_DELIM + std::to_string(pImpl->arguments[c].get());
        return output;
    }

    return NULL_STRING;
}

/* Return formatted minimum value. */
const std::string Base_Function::show_range_min() const
{
    if (tools::is_na(pImpl->range_min.get()))
        return NULL_STRING;
    else
        return std::to_string(pImpl->range_min.get());
}

/* Return formatted maximum value. */
const std::string Base_Function::show_range_max() const
{
    if (tools::is_na(pImpl->range_max.get()))
        return NULL_STRING;
    else
        return std::to_string(pImpl->range_max.get());
}

//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//

/* Initial defaults for function. */
void Base_Function::_init_defaults_(const std::string &name,
                                    const str_vec &arg_name,
                                    const dbl_vec &arg_value,
                                    double domain_min,
                                    double domain_max,
                                    bool is_open)
{
    // name of the function
    pImpl->name = name;

    // number of arguments
    if (arg_name.size() != arg_value.size())
        throw std::length_error(
            "Internal inconsistency for \"" + name + "\" arguments.");
    pImpl->n_args = arg_name.size();

    // assign default arguments
    pImpl->arguments = new Function_Argument[arg_name.size()];
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
    pImpl->init_range(is_open);
}
