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
private:
    double value = INF_DOUBLE;
    const double *ptr_to_value = &value; // UPDATE: these const type is not convincing me

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

public:
    std::string name = UNDEFINED;

    Function_Argument() = default;
    ~Function_Argument() = default;

    /* Assign argument getter and setter for pointer mode. */
    void assign_ptr(const double *y)
    {
        _func_set_ = &Function_Argument::_set_arg_pointer_;
        ptr_to_value = y;
    }
    inline double get() const { return *(ptr_to_value); }
    inline void set(double x) { (this->*_func_set_)(x); }
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
    std::string name = UNDEFINED;
    ushort _n_args_ = 0;
    Function_Argument *_arguments_ = nullptr;

    //----- Domain -----//
    double _range_min_ = -INF_DOUBLE; // init to -infty
    double _range_max_ = INF_DOUBLE;  // init to +infty

    //===== CONSTRUCTOR =====//
    Base_Function_Impl() = default;
    ~Base_Function_Impl() = default;

    //===== Getter =====//
    bool valid(double x) const { return (this->*_func_valid_range_)(x); }

    //===== Show =====//
    const std::string show_arg_names() const;
    const std::string show_arg_values() const;
    const std::string show_range_min() const;
    const std::string show_range_max() const;

    //===== Initial default arguments =====//
    void init_defaults(const std::string &name,
                       const dict_args &args_val,
                       const dict_args &domain = {});

    //===== Non-default arguments =====//
    ushort find(const std::string &arg_name);
    void assign_argument(const std::string &arg_name, double value);
    void assign_argument(const std::string &arg_name, const double *ptr_to_value);

    //===== Range =====//
    ushort set_range(const dict_args &new_range);

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
void Base_Function_Impl::assign_argument(const std::string &arg_name, double value)
{
    ushort it = find(arg_name);
    if (it == _n_args_)
        throw std::invalid_argument(
            "Unknown argument '" + arg_name + "' for '" + name + "'.");
    else
        _arguments_[it].set(value);
}

/* Overload assignment of non default POINTER arguments. */
void Base_Function_Impl::assign_argument(const std::string &arg_name, const double *ptr_to_value)
{
    ushort it = find(arg_name);
    if (it == _n_args_)
        throw std::invalid_argument(
            "Unknown argument '" + arg_name + "' for '" + name + "'.");
    else
        _arguments_[it].assign_ptr(ptr_to_value);
}

/* Overload constructor for argument. */
void Base_Function_Impl::init_defaults(const std::string &name,
                                       const dict_args &args_val,
                                       const dict_args &domain)
{
    this->name = name;
    this->_n_args_ = args_val.size();
    this->_arguments_ = new Function_Argument[args_val.size()];

    ushort c = 0;
    for (const auto &pair : args_val)
    {
        _arguments_[c].name = pair.first;
        _arguments_[c].set(pair.second);
        ++c;
    }

    set_range(domain); // assign domain
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
ushort Base_Function_Impl::set_range(const dict_args &new_range)
{
    double min_val = _range_min_;
    double max_val = _range_max_;
    bool is_new_min = false;
    bool is_new_max = false;
    ushort output = 0;

    //===== Read the new range =====//
    for (auto const &x : new_range)
    {
        if (x.first == "min")
            is_new_min = true;
        else if (x.first == "max")
            is_new_max = true;
        else
            throw std::invalid_argument(
                "Invalid range keyword '" + x.first + "'. Provide only 'min' and/or 'max'.");

        // load the new values
        if (is_new_min)
            min_val = x.second;
        else if (is_new_max)
            max_val = x.second;
    }

    //===== Adjust new min value and domain =====//
    if (find("min") < _n_args_)
    {
        // Overwrite min if it exists as argument (e .g. uniform distribution)
        min_val = _arguments_[find("min")].get();
        if (is_new_min)
            tools::warning(
                "Overwriting range min with argument min " + std::to_string(min_val) + " in '" + name + "'.");
        is_new_min = true;
    }
    else if (is_new_min)
    {
        // Check compatibility with min from domain (e .g. uniform distribution)
        if (min_val < _range_min_)
            throw std::invalid_argument(
                "Invalid min value " + std::to_string(min_val) + " for distribution '" + name + "'.");
    }

    //===== Adjust new max value and domain =====//
    if (find("max") < _n_args_)
    {
        // Overwrite max if it exists as argument
        max_val = _arguments_[find("max")].get();
        if (is_new_max)
            tools::warning(
                "Overwriting range max with argument max " + std::to_string(max_val) + " in '" + name + "'.");
        is_new_max = true;
    }
    else if (is_new_max)
    {
        // Check compatibility with max from domain
        if (max_val > _range_max_)
            throw std::invalid_argument(
                "Invalid max value " + std::to_string(max_val) + " for distribution '" + name + "'.");
    }

    //===== Assign range and respective validation method =====//

    // at least one from previous range (domain) should be finite
    // UPDATE : to be moved to NA...
    is_new_min = (is_new_min != -INF_DOUBLE) || (_range_min_ != -INF_DOUBLE);
    is_new_max = (is_new_max != INF_DOUBLE) || (_range_max_ != INF_DOUBLE);

    if ((!is_new_min) && (!is_new_max))
    {
        output = 0;
        _func_valid_range_ = &Base_Function_Impl::_valid_range_always_;
    }
    else if ((is_new_min) && (!is_new_max))
    {
        output = 1;
        _range_min_ = min_val;
        _func_valid_range_ = &Base_Function_Impl::_valid_range_min_;
    }
    else if ((!is_new_min) && (is_new_max))
    {
        output = 2;
        _range_max_ = max_val;
        _func_valid_range_ = &Base_Function_Impl::_valid_range_max_;
    }
    else // ((is_new_min) && (is_new_max))
    {
        output = 3;
        if (_range_min_ >= _range_max_)
        {
            throw std::invalid_argument(
                "Inconsistent range {min:" + std::to_string(_range_min_) + ", max:" + std::to_string(_range_max_) + "}.");
        }
        else
        {
            _range_min_ = min_val;
            _range_max_ = max_val;
            _func_valid_range_ = &Base_Function_Impl::_valid_range_both_;
        }
    }

    return output;
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

    return UNDEFINED;
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

    return UNDEFINED;
}

/* */
const std::string Base_Function_Impl::show_range_min() const
{
    return (_range_min_ > -INF_DOUBLE) ? std::to_string(_range_min_) : UNDEFINED;
}

/* */
const std::string Base_Function_Impl::show_range_max() const
{
    return (_range_max_ < INF_DOUBLE) ? std::to_string(_range_min_) : UNDEFINED;
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
// note: just for initial non-default arguments
void Base_Function::assign_args(const dict_args &new_args_val)
{
    for (const auto &pair : new_args_val)
        pImpl->assign_argument(pair.first, pair.second);
}
void Base_Function::assign_args(const dict_args_ptr &new_args_ptr)
{
    for (const auto &pair : new_args_ptr)
        pImpl->assign_argument(pair.first, pair.second);
}

//----- RANGE -----//
ushort Base_Function::set_range(const dict_args &new_range) { return pImpl->set_range(new_range); }
double Base_Function::get_range_min() const { return pImpl->_range_min_; }
double Base_Function::get_range_max() const { return pImpl->_range_max_; }

/* Overload constructor for argument VALUES. */
void Base_Function::init_defaults(const std::string &name,
                                  const dict_args &args_val,
                                  const dict_args &domain)
{
    pImpl->init_defaults(name, args_val, domain);
}
