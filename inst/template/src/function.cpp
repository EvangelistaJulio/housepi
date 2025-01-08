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
    std::string name = UNDEFINED_STRING;

    Function_Argument() = default;
    ~Function_Argument() = default;

    inline double get() const { return *(ptr_to_value); }
    inline void set(double x) { (this->*_func_set_)(x); }
    void assign_ptr(const double *y)
    {
        _func_set_ = &Function_Argument::_set_arg_pointer_;
        ptr_to_value = y;
    }

private:
    double value = UNDEFINED_DOUBLE;
    const double *ptr_to_value = &value;

    // note : the set_arg is only used for the adaptive
    // for this reason a pointer argument can't be set, but just assigned
    void (Function_Argument::*_func_set_)(double) = &Function_Argument::_set_arg_value_;
    void _set_arg_value_(double x) { value = x; }
    void _set_arg_pointer_(double x)
    {
        UNUSED(x);
        throw std::runtime_error(
            "Pointer argument '" + name + "' can not be changed.");
    }
};

//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//

class Base_Function_Impl
{
public:
    //===== Arguments =====//
    std::string name = UNDEFINED_STRING;
    ushort _n_args_ = 0;
    Function_Argument *_arguments_ = nullptr;

    //----- Domain -----//
    double _range_min_ = -UNDEFINED_DOUBLE; // UPDATE: not elegant minus...
    double _range_max_ = UNDEFINED_DOUBLE;

    //===== CONSTRUCTOR =====//
    Base_Function_Impl() = default;
    ~Base_Function_Impl() = default;

    //===== Initial default arguments =====//
    void init_defaults(const std::string &name,
                       const str_vec &arg_name,
                       const dbl_vec &arg_value,
                       double domain_min,
                       double domain_max);

    //===== Range =====//
    ushort set_range(double minimum = UNDEFINED_DOUBLE,
                     double maximum = UNDEFINED_DOUBLE);

    //===== Getter =====//
    bool valid(double x) const { return (this->*_func_valid_range_)(x); }

    //===== Show =====//
    const std::string show_arg_names() const;
    const std::string show_arg_values() const;
    const std::string show_range_min() const;
    const std::string show_range_max() const;

    //===== Non-default arguments =====//
    ushort find(const std::string &arg_name);
    void assign_scalar_arg(const std::string &arg_name, double value);
    void assign_pointer_arg(const std::string &arg_name, const double *ptr_to_value);

private:
    //===== VALID =====//
    // alternative valid range methods
    bool (Base_Function_Impl::*_func_valid_range_)(double) const = &Base_Function_Impl::_valid_range_always_;
    bool _valid_range_always_(double x) const
    {
        UNUSED(x);
        return true;
    }
    bool _valid_range_min_(double x) const { return (x > _range_min_); }
    bool _valid_range_max_(double x) const { return (x < _range_max_); }
    bool _valid_range_both_(double x) const { return _valid_range_min_(x) && _valid_range_max_(x); }
};

/////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////  Base Function : methods  //////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

ushort Base_Function_Impl::find(const std::string &arg_name)
{
    for (ushort c = 0; c < _n_args_; ++c)
        if (_arguments_[c].name == arg_name)
            return c;
    return _n_args_;
}

/* Overload assignment of non default POINTER arguments. */
void Base_Function_Impl::assign_scalar_arg(const std::string &arg_name, double x)
{
    ushort it = find(arg_name);
    if (it == _n_args_)
        throw std::invalid_argument(
            "Unknown argument '" + arg_name + "' for '" + name + "'.");
    else
        _arguments_[it].set(x);
}

/* Overload assignment of non default POINTER arguments. */
void Base_Function_Impl::assign_pointer_arg(const std::string &arg_name, const double *y)
{
    ushort it = find(arg_name);
    if (it == _n_args_)
        throw std::invalid_argument(
            "Unknown argument '" + arg_name + "' for '" + name + "'.");
    else
        _arguments_[it].assign_ptr(y);
}

/* Overload constructor for argument. */
void Base_Function_Impl::init_defaults(const std::string &name,
                                       const str_vec &arg_name,
                                       const dbl_vec &arg_value,
                                       double domain_min,
                                       double domain_max)
{
    this->name = name;
    this->_n_args_ = arg_name.size();

    if (this->_n_args_ != arg_value.size())
        throw std::length_error(
            "Internal inconsistency for '" + name + "' arguments.");

    // assign default arguments
    this->_arguments_ = new Function_Argument[arg_name.size()];
    for (ushort c = 0; c < this->_n_args_; ++c)
    {
        _arguments_[c].name = arg_name[c];
        _arguments_[c].set(arg_value[c]);
    }

    // assign domain
    set_range(domain_min, domain_max);
}

//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//

/////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////  Base Function : set range  ////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

/* Set the range of validity of the value, and the corresponding validation function. */
ushort Base_Function_Impl::set_range(double lower, double upper)
{
    // Ignore lower if 'min' exists as argument (e .g. uniform distribution)
    if (find("min") < _n_args_)
    {
        if (!tools::is_na(lower))
            tools::warning(
                "Overwriting range min with argument min in '" + name + "'.");
        lower = _arguments_[find("min")].get();
    }
    // Ignore upper if 'max' exists as argument
    if (find("max") < _n_args_)
    {
        if (!tools::is_na(upper))
            tools::warning(
                "Overwriting range max with argument max in '" + name + "'.");
        upper = _arguments_[find("max")].get();
    }

    // assign range and respective validation method
    // ... if at least one from previous range (domain) should be finite
    bool is_lower_bound = !(tools::is_na(lower) && tools::is_na(_range_min_));
    bool is_upper_bound = !(tools::is_na(upper) && tools::is_na(_range_max_));

    if ((!is_lower_bound) && (!is_upper_bound))
    {
        _func_valid_range_ = &Base_Function_Impl::_valid_range_always_;
        return 0;
    }
    else if ((is_lower_bound) && (!is_upper_bound))
    {
        if (!tools::is_na(lower))
            _range_min_ = lower;
        _func_valid_range_ = &Base_Function_Impl::_valid_range_min_;
        return 1;
    }
    else if ((!is_lower_bound) && (is_upper_bound))
    {
        if (!tools::is_na(upper))
            _range_max_ = upper;
        _func_valid_range_ = &Base_Function_Impl::_valid_range_max_;
        return 2;
    }
    else // ((is_lower_bound) && (is_upper_bound))
    {
        if (!tools::is_na(lower))
            _range_min_ = lower;
        if (!tools::is_na(upper))
            _range_max_ = upper;
        _func_valid_range_ = &Base_Function_Impl::_valid_range_both_;

        if (_range_min_ >= _range_max_)
            throw std::invalid_argument(
                "Inconsistent range {min:" + std::to_string(_range_min_) + ", max:" + std::to_string(_range_max_) + "}.");

        return 3;
    }
}

/////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////  Base Function : show args  ////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

/* Return formatted argument names of the distribution. */
const std::string Base_Function_Impl::show_arg_names() const
{
    if (_n_args_ > 0)
    {
        std::string output = _arguments_[0].name;
        for (ushort c = 1; c < _n_args_; ++c)
            output += C_LIST_DELIM + _arguments_[c].name;
        return output;
    }

    return UNDEFINED_STRING;
}

/* Return formatted argument names of the distribution. */
const std::string Base_Function_Impl::show_arg_values() const
{
    if (_n_args_ > 0)
    {
        std::string output = std::to_string(_arguments_[0].get());
        for (ushort c = 1; c < _n_args_; ++c)
            output += C_LIST_DELIM + std::to_string(_arguments_[c].get());
        return output;
    }

    return UNDEFINED_STRING;
}

/* */
const std::string Base_Function_Impl::show_range_min() const
{
    return tools::is_na(_range_min_) ? UNDEFINED_STRING : std::to_string(_range_min_);
}

/* */
const std::string Base_Function_Impl::show_range_max() const
{
    return tools::is_na(_range_max_) ? UNDEFINED_STRING : std::to_string(_range_max_);
}

//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//

/////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////  Base Function : implementations  //////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

//----- CONSTRUCTORS -----//
Base_Function::Base_Function() : pImpl(std::make_unique<Base_Function_Impl>()) {};
Base_Function::~Base_Function() = default;

//----- FUNCTION NAME -----//
const std::string &Base_Function::get_name() const { return pImpl->name; }

//----- ARGUMENTS -----//
double Base_Function::get_arg(ushort i) const { return pImpl->_arguments_[i].get(); }
void Base_Function::set_arg(ushort i, double x) { pImpl->_arguments_[i].set(x); }
const std::string Base_Function::show_arg_names() const { return pImpl->show_arg_names(); }
const std::string Base_Function::show_arg_values() const { return pImpl->show_arg_values(); }

//----- RANGE -----//
const std::string Base_Function::show_range_min() const { return pImpl->show_range_min(); }
const std::string Base_Function::show_range_max() const { return pImpl->show_range_max(); }

//----- VALID -----//
bool Base_Function::valid(double x) const { return pImpl->valid(x); }

//----- ARGUMENTS -----//
void Base_Function::assign_scalar_arg(const std::string &arg_name, double x)
{
    pImpl->assign_scalar_arg(arg_name, x);
}

void Base_Function::assign_pointer_arg(const std::string &arg_name, const double *y)
{
    pImpl->assign_pointer_arg(arg_name, y);
}

//----- RANGE -----//
ushort Base_Function::set_range(double lower, double upper)
{
    return pImpl->set_range(lower, upper);
}
double Base_Function::get_range_min() const { return pImpl->_range_min_; }
double Base_Function::get_range_max() const { return pImpl->_range_max_; }

/* Overload constructor for argument VALUES. */
void Base_Function::_init_defaults_(const std::string &name,
                                    const str_vec &arg_name,
                                    const dbl_vec &arg_value,
                                    double domain_min,
                                    double domain_max)
{
    pImpl->init_defaults(name, arg_name, arg_value, domain_min, domain_max);
}
