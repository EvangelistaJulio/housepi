/* housepi - Epidemiological Analysis of Household Outbreaks.
 *
 * Copyright © 2026 Institut Pasteur (Paris), CNRS and INSERM.
 * Distributed under the terms of the GNU General Public License (GPLv3).
 * See the LICENSE file for details.
 */

#include <memory>
#include <stdexcept>
#include <string>

#include "tools/nan.h"
#include "tools/string.h"

#include "attribute.h"
#include "function.h"

/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////  *** Attribute : class ***  /////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

Attribute::Attribute() = default;
Attribute::~Attribute() = default;

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

void Attribute::init_function(std::shared_ptr<Function> f)
{
    // check the function is different from the current one
    if (f == nullptr)
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
    pFunc_ = f;
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
        init_function(other.function());
        return;
    }
    __builtin_unreachable();
}

inline double Attribute::value() const
{
    switch (mode_)
    {
    case Mode::Empty:
        throw std::runtime_error(
            "Impossible to return value for " + err_name() + ".");
    case Mode::Numeric:
        return value_;
    case Mode::Pointer:
        return pAttr_->value();
    case Mode::Function:
        return pFunc_->evaluate();
    default:
        __builtin_unreachable();
    }
}
