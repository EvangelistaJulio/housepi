/* housepi - Epidemiological Analysis of Household Outbreaks.
 *
 * Copyright © 2026 Institut Pasteur (Paris), CNRS and INSERM.
 * Distributed under the terms of the GNU General Public License (GPLv3).
 * See the LICENSE file for details.
 */

#pragma once

#include <memory>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

#include "tools/nan.h"

/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////// *** Attribute : class *** ///////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

class Function;

/// Attribute parent class.
class Attribute
{
protected:
    enum class Mode : uint8_t
    {
        Empty,
        Numeric,
        Pointer,
        Function
    } mode_ = Mode::Empty;

public:
    ///////////////////////////
    //////  CONSTRUCTOR  //////
    ///////////////////////////

    Attribute();
    Attribute(double x) { init_numeric(x); }
    virtual ~Attribute();

    //////////////////////
    //////  SETTER  //////
    //////////////////////

    /// Set attribute name.
    void name(const std::string &name) { name_ = name; }
    /// Set error type.
    void set_err_type(const std::string &err_type) { err_type_ = err_type; }

    /// Initialize attribute as numeric.
    void init_numeric(double x)
    {
        if (err_type_ == "empty")
            set_err_type("numeric");

        mode_ = Mode::Numeric;
        value_ = x;
        pAttr_ = nullptr;
        pFunc_ = nullptr;
        is_user_const_ = true;
    }
    /// Initialize attribute as pointer.
    void init_pointer(const Attribute *y)
    {
        if (y == nullptr)
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
    /// Initialize attribute as function.
    void init_function(std::shared_ptr<Function> f);
    /// Productory of pointed attributes values.
    void init_productory(const std::vector<const Attribute *> &other_ptrs);
    /// Transfer from other Attribute.
    void transfer(const Attribute &other);

    /// Set attribute value.
    void set(double x)
    {
        switch (mode_)
        {
        case Mode::Empty:
            mode_ = Mode::Numeric;
            value_ = x;
            return;
        case Mode::Numeric:
            value_ = x;
            return;
        case Mode::Pointer:
        case Mode::Function:
            throw std::runtime_error(
                "Tried to change non numeric " + err_name() + ".");
        }
        __builtin_unreachable();
    }

    //////////////////////
    //////  GETTER  //////
    //////////////////////

    /// Get attribute numeric value.
    /// NOTE: this can't become noexcept because of special function calls.
    double value() const;

    /// Get pointer to attribute value.
    const Attribute *pointer() const noexcept
    {
        switch (mode_)
        {
        case Mode::Empty:
        case Mode::Numeric:
            return this;
        case Mode::Pointer:
            return pAttr_->pointer();
        case Mode::Function:
            return this;
        }
        __builtin_unreachable();
    }
    /// Get pointer to attribute function.
    /// FIXME: this should not be a constant method.
    std::shared_ptr<Function> function() const noexcept
    {
        switch (mode_)
        {
        case Mode::Empty:
        case Mode::Numeric:
            return nullptr;
        case Mode::Pointer:
            return pAttr_->function();
        case Mode::Function:
            return pFunc_;
        }
        __builtin_unreachable();
    }

    /////////////////////
    //////  FLAGS  //////
    /////////////////////

    /// Wheter it's constant.
    bool is_constant() const noexcept;
    /// Wheter it's discrete.
    virtual bool is_discrete() const noexcept;

    /// Wheter it has no definitions.
    bool is_empty() const noexcept { return mode_ == Mode::Empty; }
    /// Wheter it's a number.
    bool is_numeric() const noexcept { return mode_ == Mode::Numeric; }
    /// Wheter it's a pointer.
    bool is_pointer() const noexcept { return mode_ == Mode::Pointer; }
    /// Wheter it's a function.
    bool is_function() const noexcept { return mode_ == Mode::Function; }

    /// Wheter it's able to return. 
    bool is_callable() const noexcept;

    ////////////////////
    //////  SHOW  //////
    ////////////////////

    /// Get attribute name.
    std::string_view name() const noexcept { return name_; }
    /// Get attribute name to display in errors.
    std::string err_name() const { return err_type_ + " '" + std::string(name()) + "'"; }
    /// Get attribute value to display in errors.
    std::string str_val() const;
    /// Get a summary of the attribute.
    virtual std::string str_summary() const;

protected:
    std::string err_type_ = "empty";            /// type name displayed in errors
    std::string name_ = "";                     /// name of argument
    bool is_user_const_ = false;                /// user flag for constant
    double value_ = hsp::na::DOUBLE;            /// internal value
    const Attribute *pAttr_ = nullptr;          /// pointer to attribute
    std::shared_ptr<Function> pFunc_ = nullptr; /// pointer to function
};
