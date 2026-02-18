/* housepi - Epidemiological Analysis of Household Outbreaks.
 *
 * Copyright © 2026 Institut Pasteur (Paris), CNRS and INSERM.
 * Distributed under the terms of the GNU General Public License (GPLv3).
 * See the LICENSE file for details.
 */

#include <cmath>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

#include "tools/messages.h"
#include "tools/nan.h"
#include "tools/string.h"
#include "tools/types.h"

#include "function.h"

/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////  *** Attribute : class ***  /////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

Attribute::Attribute() = default;
Attribute::~Attribute() = default;

Attribute::Attribute(const Attribute &other)
    : mode_(other.mode_),
      err_type_(other.err_type_),
      name_(other.name_),
      is_user_const_(other.is_user_const_),
      value_(other.value_),
      pAttr_(other.pAttr_),
      pFunc_(other.pFunc_ ? other.pFunc_->clone() : nullptr) {}

Attribute::Attribute(Attribute &&) = default;

Attribute &Attribute::operator=(Attribute other)
{
    swap(other);
    return *this;
}

void Attribute::swap(Attribute &other) noexcept
{
    std::swap(mode_, other.mode_);
    std::swap(err_type_, other.err_type_);
    std::swap(name_, other.name_);
    std::swap(is_user_const_, other.is_user_const_);
    std::swap(value_, other.value_);
    std::swap(pAttr_, other.pAttr_);
    std::swap(pFunc_, other.pFunc_);
}

std::string Attribute::str_val() const
{
    switch (mode_)
    {
    case Mode::Empty:
        return hsp::str::format(hsp::na::DOUBLE);
    case Mode::Numeric:
        return hsp::str::format(value_);
    case Mode::Pointer:
        return pAttr_->str_val();
    case Mode::Function:
        if (pFunc_->is_callable_implicitly())
            return hsp::str::format(pFunc_->evaluate());
        return "<>";
    }
    __builtin_unreachable();
}

std::string Attribute::str_summary() const
{
    std::string sms = err_name() + " ";

    sms += (is_function() && !pFunc_->is_callable_implicitly()) ? "<> | " : str_val() + " | ";
    sms += is_discrete() ? "discrete" : "real";
    sms += is_constant() ? " constant" : " variable";

    switch (mode_)
    {
    case Mode::Empty:
        return sms + " empty";
    case Mode::Numeric:
        return sms + " numeric";
    case Mode::Pointer:
        return sms + " pointer -> " + pAttr_->err_name();
    case Mode::Function:
        return sms + " function " + pFunc_->str_summary();
    }
    __builtin_unreachable();
}

bool Attribute::is_constant() const noexcept
{
    switch (mode_)
    {
    case Mode::Empty:
        return false;
    case Mode::Numeric:
        return is_user_const_;
    case Mode::Pointer:
        return pAttr_->is_constant();
    case Mode::Function:
        return pFunc_->is_constant();
    }
    __builtin_unreachable();
}

bool Attribute::is_discrete() const noexcept
{
    switch (mode_)
    {
    case Mode::Empty:
        return false;
    case Mode::Numeric:
        return value_ == std::floor(value_);
    case Mode::Pointer:
        return pAttr_->is_discrete();
    case Mode::Function:
        return pFunc_->is_discrete();
    }
    __builtin_unreachable();
}

bool Attribute::is_callable() const noexcept
{
    switch (mode_)
    {
    case Mode::Empty:
        return false;
    case Mode::Numeric:
        return !hsp::na::is(value_);
    case Mode::Pointer:
        return pAttr_->is_callable();
    case Mode::Function:
        return pFunc_->is_callable_implicitly();
    }
    __builtin_unreachable();
}

/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////  Attribute : setter //////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

void Attribute::init_numeric(double x)
{
    if (err_type_ == "empty")
        set_err_type("numeric");

    mode_ = Mode::Numeric;
    value_ = x;
    pAttr_ = nullptr;
    pFunc_ = nullptr;
    is_user_const_ = true;
}

void Attribute::init_pointer(const Attribute *y)
{
    if (!y)
        throw std::invalid_argument(
            "Tried to assign empty pointer to " + err_name() + ".");

    if (err_type_ == "empty")
        set_err_type("pointer");

    mode_ = Mode::Pointer;
    value_ = hsp::na::DOUBLE;
    pAttr_ = y;
    pFunc_ = nullptr;
    is_user_const_ = false;
}

void Attribute::init_function(std::unique_ptr<Function> f)
{
    // check the function is different from the current one
    if (!f)
        throw std::invalid_argument(
            "Tried to assign empty function to " + err_name() + ".");
    if (f == pFunc_)
        throw std::invalid_argument(
            "Tried to assign " + err_name() + " as its own argument.");

    if (err_type_ == "empty")
        err_type_ = "function";
    if (hsp::na::is(name_))
        name_ = f->name();
    mode_ = Mode::Function;
    value_ = hsp::na::DOUBLE;
    pAttr_ = nullptr;
    pFunc_ = std::move(f);
    is_user_const_ = false;
}

void Attribute::transfer(const Attribute &other)
{
    switch (other.mode_)
    {
    case Mode::Empty:
        // NOTE: nothing to transfer if other is pure empty
        if (other.err_type_ != "empty")
            init_pointer(other.pointer());
        return;
    case Mode::Numeric:
        if (other.is_user_const_)
            init_numeric(other.value());
        else
            init_pointer(other.pointer());
        return;
    case Mode::Pointer:
        init_pointer(other.pointer());
        return;
    case Mode::Function:
        init_function(other.function()->clone());
        return;
    }
    __builtin_unreachable();
}
void Attribute::optimize()
{
    if (is_constant())
        init_numeric(value());
    else if (is_user_const_)
        throw std::invalid_argument(
            "Impossible to set " + err_name() + " as constant. Check its dependencies.");

    // reduce function indirection
    if (is_pointer() && pointer()->is_function())
        init_function(pointer()->function()->clone());
    if (is_function())
        pFunc_->optimize();
}

void Attribute::init_func_attr(std::string_view key, const Attribute &attr)
{
    switch (mode_)
    {
    case Mode::Function:
        pFunc_->init_attr(key, attr);
        return;
    default:
        throw std::runtime_error(
            "Attempted to initialize function attribute of " + err_name() + ".");
    }
}
