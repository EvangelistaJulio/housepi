/* housepi - Epidemiological Analysis of Household Outbreaks.
 *
 * Copyright © 2026 Institut Pasteur (Paris), CNRS and INSERM.
 * Distributed under the terms of the GNU General Public License (GPLv3).
 * See the LICENSE file for details.
 */

#ifndef HOUSEPI__ATTRIBUTE__H
#define HOUSEPI__ATTRIBUTE__H

#include <memory>

#include "tools/types.h"

/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////// *** Attribute : class *** ///////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

// forward declaration
class Function;

/// Attribute parent class.
class Attribute
{
public:
    /////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////  Attribute : constructor  ///////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////

    /// Empty attribute constructor.
    Attribute();
    /// Numeric attribute.
    explicit Attribute(double x);
    /// Pointer to double attribute.
    explicit Attribute(const Attribute *y);
    /// Function attribute.
    explicit Attribute(const std::shared_ptr<Function> f);
    /// Function attribute by registry choice.
    explicit Attribute(const std::string &s, bool is_primitive = false);
    /// Attribute deconstructor.
    virtual ~Attribute();

    /////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////  Attribute : setter  //////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////

    /// Set attribute value.
    void set(double x) { (this->*do_set_)(x); }
    /// Set attribute name. Eventually also err type.
    void set_name(const std::string &name, const std::string &err_type = hsp::na::STRING)
    {
        this->name_ = name;
        if (!hsp::na::is(err_type))
            this->err_type_ = err_type;
    }

    /// Initialize attribute as numeric.
    void init_numeric(double x);
    /// Initialize attribute as pointer.
    void init_pointer(const Attribute *y);
    /// Initialize attribute as function.
    void init_function(const std::shared_ptr<Function> f);
    /// Initialize attribute as function given the function name.
    void init_function(const std::string &s, bool is_primitive = false);
    /// Transfer from other Attribute.
    void transfer(const Attribute &other);

    /////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////  Attribute : getter  //////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////

    /// Get attribute numeric value.
    double get() const { return (this->*get_)(); }
    /// Get pointer to attribute value.
    const Attribute *get_ptr() const { return (this->*get_ptr_)(); }
    /// Get pointer to attribute function.
    const std::shared_ptr<Function> get_function() const { return pFunc_; }

    /// Wheter it's constant.
    virtual bool is_constant() const;
    /// Wheter it's discrete.
    virtual bool is_discrete() const;

    /// Wheter it's a number.
    bool is_numeric() const { return !hsp::na::is(value_); }
    /// Wheter it's a pointer.
    bool is_pointer() const { return pAttr_ != nullptr; }
    /// Wheter it's a function.
    bool is_function() const { return pFunc_ != nullptr; }
    /// Wheter it has no definitions.
    bool is_empty() const { return !is_numeric() && !is_pointer() && !is_function(); }

    /// Get attribute name.
    std::string get_name() const { return name_; }
    /// Get attribute name to display in errors.
    std::string err_name() const { return err_type_ + " '" + name_ + "'"; }
    /// Get a description about the functional status.
    virtual std::string str_summary() const;

protected:
    std::string err_type_ = "attribute";        // type name displayed in errors
    std::string name_ = hsp::na::STRING;        /// name of argument
    double value_ = hsp::na::DOUBLE;            /// internal value
    const Attribute *pAttr_ = nullptr;          /// pointer to attribute
    std::shared_ptr<Function> pFunc_ = nullptr; /// pointer to function

    // Functional getter
    double (Attribute::*get_)() const = &Attribute::do_get_default_;
    /// Get attribute value from internal value.
    double do_get_default_() const { return value_; }
    /// Get attribute value from external value.
    double do_get_from_pointer_() const { return pAttr_->get(); }
    /// Get attribute value from external function.
    double do_get_from_function_() const;
    /// Get attribute error for not callable functions.
    double do_get_error_() const;

    // Functional pointer getter
    const Attribute *(Attribute::*get_ptr_)() const = &Attribute::do_get_ptr_default_;
    /// Get pointer to attribute error.
    const Attribute *do_get_pointer_error_() const { throw std::runtime_error(
        "Tried to return a pointer to " + err_name() + "."); }
    /// Get pointer to attribute value internal value.
    const Attribute *do_get_ptr_default_() const { return this; }
    /// Get pointer to attribute value from external value.
    const Attribute *do_get_pointer_from_pointer_() const { return pAttr_; }

    // Functional setter
    void (Attribute::*do_set_)(double) = &Attribute::do_set_default_;
    /// Set attribute internal value.
    virtual void do_set_default_(double x) { value_ = x; }
    /// Set error. Non numeric attributes can't be set, but just assigned.
    void do_set_error_(double) { throw std::runtime_error(
        "Tried to change non numeric " + err_name() + "."); }
};

#endif // HOUSEPI__ATTRIBUTE__H
