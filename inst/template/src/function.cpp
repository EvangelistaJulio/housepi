/*
 * function.cpp
 *
 * Copyright (C) 2024 Francesco Camaglia, Institut Pasteur
 */

#include "function.h"

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
    std::string _name_ = UNDEFINED;
    std::vector<std::string> _args_names_ = {}; // ordered list of argument names
    double *_args_val_ = nullptr;               // internal arguments (fixed at initialization)
    const double **_args_ptr_ = nullptr;        // external arguments (pointer to variable)

    //===== Constructors =====//
    Base_Function_Impl() {};
    ~Base_Function_Impl() {};

    void assign_val_mode();
    void assign_ptr_mode();

    //===== Getter =====//
    bool valid(double x) const { return (this->*_func_valid_range_)(x); }
    double get_arg(int i) const { return (this->*_func_get_arg_)(i); }

    //===== Setter (for internal arguments) =====//
    void set_arg(int i, double new_value) { (this->*_func_set_arg_)(i, new_value); }

    //===== Show =====//
    std::string show_arg_names() const;
    std::string show_arg_values() const;
    std::string show_range_min() const { return (_range_min_ > -INF_DOUBLE) ? std::to_string(_range_min_) : UNDEFINED; }
    std::string show_range_max() const { return (_range_max_ < INF_DOUBLE) ? std::to_string(_range_min_) : UNDEFINED; }

    //===== Initial non-default arguments =====//
    void assign_args(const dictargs &new_args_val);
    void assign_args(const dictargs_ptr &new_args_ptr);

    //===== Range =====//
    ushort set_range(const dictargs &new_range);
    double get_range_min() const { return _range_min_; }
    double get_range_max() const { return _range_max_; }

private:
    //===== Argument getter : functional methods =====//
    double (Base_Function_Impl::*_func_get_arg_)(int) const;
    double _get_arg_val_(int i) const { return _args_val_[i]; }
    double _get_arg_ptr_(int i) const { return *(_args_ptr_[i]); }

    //===== Argument setter : functional methods =====//
    // note : the set_arg is only used for the adaptive
    // for this reason a pointer argument can't be set, but just assigned
    void (Base_Function_Impl::*_func_set_arg_)(int, double);
    void _set_arg_val_(int i, double new_value) { _args_val_[i] = new_value; }
    void _set_arg_ptr_(int i, double new_value)
    {
        UNUSED(new_value);
        std::ostringstream oss;
        oss << "Inconsistent change of pointer argument'" << _args_names_[i]
            << "' within '" << _name_ << "'.";
        tools::error(oss.str());
    }

    //===== Valid : functional methods =====//
    double _range_min_ = -INF_DOUBLE; // init to -infty
    double _range_max_ = INF_DOUBLE;  // init to +infty

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

/* Assign argument getter and setter for scalar value mode. */
void Base_Function_Impl::assign_val_mode()
{
    // assign getter
    _func_get_arg_ = &Base_Function_Impl::_get_arg_val_;
    // assign setter
    _func_set_arg_ = &Base_Function_Impl::_set_arg_val_;
}

/* Assign argument getter and setter for pointer mode. */
void Base_Function_Impl::assign_ptr_mode()
{
    // assign getter
    _func_get_arg_ = &Base_Function_Impl::_get_arg_ptr_;
    // assign setter
    _func_set_arg_ = &Base_Function_Impl::_set_arg_ptr_;
}

/* Overload assignment of non default POINTER arguments. */
void Base_Function_Impl::assign_args(const dictargs &args_val)
{
    for (const auto &pair : args_val)
    {
        auto it = std::find(_args_names_.begin(), _args_names_.end(), pair.first);
        if (it == _args_names_.end())
        {
            std::ostringstream oss;
            oss << "Unknown argument '" << pair.first << "' for '" << _name_ << "'.";
            throw std::invalid_argument(oss.str());
        }
        else
        {
            _args_val_[std::distance(_args_names_.begin(), it)] = pair.second;
        }
    }
}

/* Overload assignment of non default POINTER arguments. */
void Base_Function_Impl::assign_args(const dictargs_ptr &args_ptr)
{
    for (const auto &pair : args_ptr)
    {
        auto it = std::find(_args_names_.begin(), _args_names_.end(), pair.first);
        if (it == _args_names_.end())
        {
            std::ostringstream oss;
            oss << "Unknown argument '" << pair.first << "' for '" << _name_ << "'.";
            throw std::invalid_argument(oss.str());
        }
        else
        {
            _args_ptr_[std::distance(_args_names_.begin(), it)] = pair.second;
        }
    }
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
ushort Base_Function_Impl::set_range(const dictargs &new_range)
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
        {
            if (is_new_min)
            {
                throw std::invalid_argument("Multiple definition for 'min' in '" + _name_ + "'.");
            }
            else
            {
                min_val = x.second;
                is_new_min = true;
            }
        }
        else if (x.first == "max")
        {
            if (is_new_max)
            {
                throw std::invalid_argument("Multiple definition for 'max' in '" + _name_ + "'.");
            }
            else
            {
                max_val = x.second;
                is_new_max = true;
            }
        }
        else
        {
            std::ostringstream oss;
            oss << "Invalid range keyword '" << x.first << "'. Provide only 'min' and/or 'max'.";
            throw std::invalid_argument(oss.str());
        }
    }

    //===== Adjust new min value and domain =====//
    if (std::find(_args_names_.begin(), _args_names_.end(), "min") != _args_names_.end())
    {
        //----- Overwrite min if it exists as argument (e .g. uniform distribution) -----//
        min_val = get_arg(std::find(_args_names_.begin(), _args_names_.end(), "min") - _args_names_.begin());
        if (is_new_min)
        {
            std::ostringstream oss;
            oss << "Overwriting range min with argument min " << min_val << " in '" << _name_ << "'.";
            tools::warning(oss.str());
        }
        is_new_min = true;
    }
    else if (is_new_min)
    {
        //----- Check compatibility with min from domain (e .g. uniform distribution) -----//
        if (min_val < _range_min_)
        {
            std::ostringstream oss;
            oss << "Invalid min value " << min_val << " for distribution '" << _name_ << "'.";
            throw std::invalid_argument(oss.str());
        }
    }

    //===== Adjust new max value and domain =====//
    if (std::find(_args_names_.begin(), _args_names_.end(), "max") != _args_names_.end())
    {
        //----- Overwrite max if it exists as argument -----//
        max_val = get_arg(std::find(_args_names_.begin(), _args_names_.end(), "max") - _args_names_.begin());
        if (is_new_max)
        {
            std::ostringstream oss;
            oss << "Overwriting range max with argument max " << max_val << " in '" << _name_ << "'.";
            tools::warning(oss.str());
        }
        is_new_max = true;
    }
    else if (is_new_max)
    {
        //----- Check compatibility with max from domain -----//
        if (max_val > _range_max_)
        {
            std::ostringstream oss;
            oss << "Invalid max value " << max_val << " for distribution '" << _name_ << "'.";
            throw std::invalid_argument(oss.str());
        }
    }

    //===== Assign range and respective validation method =====//

    // at least one from previous range (domain) should be finite
    // FIXME : to be moved to NA...
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
            std::ostringstream oss;
            oss << "Provided an inconsistent range {min:" << _range_min_ << ", max:" << _range_max_ << "}.";
            throw std::invalid_argument(oss.str());
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
std::string Base_Function_Impl::show_arg_names() const
{
    if (_args_names_.empty())
    {
        return UNDEFINED;
    }
    else
    {
        std::ostringstream oss;
        oss << _args_names_[0];
        for (size_t i = 1; i < _args_names_.size(); ++i)
        {
            oss << C_LIST_DELIM << _args_names_[i];
        }
        return oss.str();
    }
}

/* Return formatted argument names of the distribution. */
std::string Base_Function_Impl::show_arg_values() const
{
    if (_args_names_.empty())
    {
        return UNDEFINED;
    }
    else
    {
        std::ostringstream oss;
        oss << get_arg(0);
        for (size_t i = 1; i < _args_names_.size(); ++i)
        {
            oss << C_LIST_DELIM << get_arg(i);
        }
        return oss.str();
    }
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
std::string Base_Function::get_name() const { return pImpl->_name_; }

//----- ARGUMENTS -----//
double Base_Function::get_arg(int i) const { return pImpl->get_arg(i); }
void Base_Function::set_arg(int i, double x) { pImpl->set_arg(i, x); }
std::string Base_Function::show_arg_names() const { return pImpl->show_arg_names(); }
std::string Base_Function::show_arg_values() const { return pImpl->show_arg_values(); }

//----- RANGE -----//
std::string Base_Function::show_range_min() const { return pImpl->show_range_min(); }
std::string Base_Function::show_range_max() const { return pImpl->show_range_max(); }

//----- VALID -----//
bool Base_Function::valid(double x) const { return pImpl->valid(x); }

//----- ARGUMENTS -----//
// note: just for initial non-default arguments
void Base_Function::assign_args(const dictargs &new_args_val) { pImpl->assign_args(new_args_val); }
void Base_Function::assign_args(const dictargs_ptr &new_args_ptr) { pImpl->assign_args(new_args_ptr); }

//----- RANGE -----//
ushort Base_Function::set_range(const dictargs &new_range) { return pImpl->set_range(new_range); }
double Base_Function::get_range_min() const { return pImpl->get_range_min(); }
double Base_Function::get_range_max() const { return pImpl->get_range_max(); }

/* Overload constructor for argument VALUES. */
void Base_Function::init(const std::string &name,
                         const dictargs &args_val,
                         const std::vector<std::string> &order,
                         const dictargs &domain)
{
    pImpl->_name_ = name;
    pImpl->_args_names_ = order;
    pImpl->_args_val_ = new double[args_val.size()];
    size_t i = 0;
    for (const std::string &argID : order)
    {
        pImpl->_args_val_[i++] = args_val.at(argID);
    }
    // assign getter and setter for scalar value mode
    pImpl->assign_val_mode();
    // assign domain
    pImpl->set_range(domain);
}

/* Overload constructor for argument POINTERS. */
void Base_Function::init(const std::string &name,
                         const dictargs_ptr &args_ptr,
                         const std::vector<std::string> &order,
                         const dictargs &domain)
{
    pImpl->_name_ = name;
    pImpl->_args_names_ = order;
    pImpl->_args_ptr_ = new const double *[args_ptr.size()];
    size_t i = 0;
    for (const std::string &argID : order)
    {
        pImpl->_args_ptr_[i++] = args_ptr.at(argID);
    }
    // assign getter and setter for pointer mode
    pImpl->assign_ptr_mode();
    // assign domain
    pImpl->set_range(domain);
}