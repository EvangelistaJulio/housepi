/*
 *        __
 *       / /_  ____  __  __ ____ __  ____  __
 *      / __ \/ __ \/ / / / ___/ _ \/ __ \/ /
 *     / / / / /_/ / /_/ (__  )  __/ /_/ / /
 *    /_/ /_/\____/\____/____/\___/ ____/_/
 *                               /_/
 *
 * housepi - Epidemiological Analysis of Household Outbreaks.
 * 
 * Author: Francesco Camaglia
 * Copyright © 2025 Institut Pasteur (Paris), CNRS and INSERM.
 * See the COPYRIGHT file for details.
 * 
 * Distributed under the terms of the GNU General Public License (GPLv3).
 * See the COPYING file for details.
 *
 * function.cpp
 *
 */

#include "tools/discrete.h"
#include "tools/params.h"
#include "tools/messages.h"
#include "tools/string.h"
#include "tools/operators.h"
#include "function.h"
#include "function_registry.h"

/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////  Attribute : constructor ///////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

Attribute::Attribute(double x)
{
    set_name(hsp::na::STRING, "numeric");
    init_numeric(x);
}
Attribute::Attribute(const Attribute *y)
{
    set_name(hsp::na::STRING, "pointer");
    init_pointer(y);
}
Attribute::Attribute(const std::shared_ptr<Function> f)
{
    set_name(f->name(), "function");
    init_function(f);
}
Attribute::Attribute(const std::string &s, bool is_primitive)
{
    set_name(s, "function");
    init_function(s, is_primitive);
}

/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////  Attribute : getter //////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

std::string Attribute::str_summary() const
{
    std::string sms = "";
    if (!is_function() || pFunc_->is_implicit_callable())
        sms = std::to_string(get()) + " | ";

    if (is_discrete())
        sms += "discrete";
    else
        sms += "real";

    if (is_constant())
        sms += " constant";
    else
        sms += " variable";

    if (is_numeric())
        sms += " numeric";
    else if (is_pointer())
        sms += " pointer -> " + pAttr_->err_name();
    else if (is_function())
        sms += " function " + pFunc_->str_summary();
    else
        sms += " empty";
    return sms;
}

/// Get attribute value from external function.
///
/// NOTE: these can't be moved in the class because of forward declaration of Function
double Attribute::do_get_from_function_() const { return pFunc_->evaluate(); }
/// Get attribute error for not callable functions.
double Attribute::do_get_error_() const { throw std::runtime_error(
    "Impossible to return value for " + err_name() +
    ". Argument 'x' of " + pFunc_->err_name() + " cannot be missing."); }

/// Wheter it's constant.
/// NOTE: by construction convention, a pointer to constant must be initialized to numeric instead.
/// NOTE: if empty, then it's constant by convention.
bool Attribute::is_constant() const
{
    if (is_numeric())
        return true;
    if (is_pointer())
        return pAttr_->is_constant();
    if (is_function())
        return pFunc_->is_constant();
    return true;
}

/// Wheter it's discrete.
/// NOTE: if empty, then it's discrete by convention.
bool Attribute::is_discrete() const
{
    if (is_numeric())
        return value_ == std::floor(value_);
    if (is_pointer())
        return pAttr_->is_discrete();
    if (is_function())
        return pFunc_->is_discrete();
    return true;
}

/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////  Attribute : setter //////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

/// Initialize attribute as numeric.
void Attribute::init_numeric(double x)
{
    get_ = &Attribute::do_get_default_;
    get_ptr_ = &Attribute::do_get_ptr_default_;
    do_set_ = &Attribute::do_set_default_;
    value_ = x;
    pAttr_ = nullptr;
    pFunc_ = nullptr;
}
/// Initialize attribute as pointer.
void Attribute::init_pointer(const Attribute *y)
{
    if (y == nullptr)
        throw std::invalid_argument(
            "Tried to assign empty pointer to " + err_name() + ".");

    get_ = &Attribute::do_get_from_pointer_;
    get_ptr_ = &Attribute::do_get_pointer_from_pointer_;
    do_set_ = &Attribute::do_set_error_;
    value_ = hsp::na::DOUBLE;
    pAttr_ = y;
    pFunc_ = nullptr;
}
/// Initialize attribute as function.
void Attribute::init_function(const std::shared_ptr<Function> f)
{
    // check the function is different from the current one
    if (f == nullptr)
        throw std::invalid_argument(
            "Tried to assign empty function to " + err_name() + ".");
    if (f == pFunc_)
        throw std::invalid_argument(
            "Tried to assign " + err_name() + " as its own argument.");

    get_ = &Attribute::do_get_from_function_;
    get_ptr_ = &Attribute::do_get_pointer_error_;
    do_set_ = &Attribute::do_set_error_;
    value_ = hsp::na::DOUBLE;
    pAttr_ = nullptr;
    pFunc_ = f;
}

/// Transfer from other Attribute.
void Attribute::transfer(const Attribute &other)
{
    if (other.is_constant())
        init_numeric(other.get());
    else if (other.is_function())
        init_function(other.get_function());
    else
        init_pointer(other.get_ptr());
}

/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////  Function : constructor  /////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

/// Initial defaults for function.
Function::Function(const std::string &name,
                   const str_vec &arg_name, const dbl_vec &arg_value,
                   double domain_min, bool is_min_open,
                   double domain_max, bool is_max_open)
{
    // initialize default 'x' argument
    x_arg_.set_name(hsp::params::X);
    x_arg_.set(hsp::na::DOUBLE);

    // name of the function
    this->name_ = name;

    // number of arguments
    if (arg_name.size() != arg_value.size())
        hsp::msg::error(
            "Inconsistent arguments for " + err_name() + "." +
            "Found " + std::to_string(arg_name.size()) + " default names but " +
            std::to_string(arg_value.size()) + " values.");

    // allocate attribute for remaining arguments
    size_t n_args = arg_name.size() - (arg_name % hsp::params::X);
    arguments_ = std::vector<Attribute>(n_args);

    for (size_t k = 0; k < arg_name.size(); ++k)
        if (arg_name[k] == hsp::params::X)
        {
            x_arg_.set(arg_value[k]);
            --k;
        }
        else
        {
            arguments_[k].set_name(arg_name[k]);
            arguments_[k].set(arg_value[k]);
        }

    // Assign Domain
    // UPDATE: at the moment only place where to set open ranges
    // UPDATE: not sure if setting pointer to smtg exiting the domain I still have control
    init_attr(hsp::params::MIN, Attribute(domain_min));
    init_attr(hsp::params::MAX, Attribute(domain_max));
    init_range_helper_(is_min_open, is_max_open);
}

/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////  Function : Getter  /////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

/// Assert argument exists at assignment. Returns it's index.
Attribute &Function::find_arg_(const std::string &arg_name)
{
    if (arg_name == hsp::params::X)
        return x_arg_;

    for (Attribute &attr : arguments_)
        if (attr.get_name() == arg_name)
            return attr;

    // else format the error showing options
    str_vec options = {hsp::params::X};
    for (const Attribute &attr : arguments_)
        options.push_back(attr.get_name());
    throw std::invalid_argument(
        "Unknown argument '" + arg_name + "' for " + err_name() + "." +
        hsp::str::find_closest(arg_name, options));
}

/// Wheter return value is always constant.
bool Function::is_constant() const
{
    // check X value
    // NOTE: a function with undefined 'x' is not constant
    if (x_arg_.is_empty() || !x_arg_.is_constant())
        return false;

    // check arguments
    for (const Attribute &attr : arguments_)
        if (!attr.is_constant())
            return false;

    // check range
    if (!range_lower_.is_constant())
        return false;
    if (!range_upper_.is_constant())
        return false;

    // otherwise constant
    return true;
}

/////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////  Function : set range  //////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

/// Functional range for lower bound and no upper bound.
bool Function::init_range_min_(double min_val, bool is_min_open)
{
    if (lower() < min_val)
        throw std::invalid_argument(
            "Invoked invalid min value " + hsp::str::format(lower()) + " for " + err_name() + ".");

    if (lower() == min_val && hsp::na::is(upper()))
        return false;

    init_range_helper_(is_min_open ? lower() == min_val : false, false);
    return true;
}

/// Set the range of validity of the value, and the corresponding validation function.
void Function::init_range_helper_(bool is_min_open, bool is_max_open)
{
    if (range_lower_.is_empty() && range_upper_.is_empty()) // Unbounded
    {
        is_below_ = &Function::do_is_below_never_;
        is_above_ = &Function::do_is_above_never_;
        is_valid_ = &Function::do_is_valid_always_;
    }
    else if (!range_lower_.is_empty() && range_upper_.is_empty()) // Lower Bound
    {
        is_above_ = &Function::do_is_above_never_;
        if (is_min_open)
        {
            is_below_ = &Function::do_is_below_open_min_;
            is_valid_ = &Function::do_is_valid_open_min_;
        }
        else
        {
            is_below_ = &Function::do_is_below_close_min_;
            is_valid_ = &Function::do_is_valid_close_min_;
        }
    }
    else if (range_lower_.is_empty() && !range_upper_.is_empty()) // Upper Bound
    {
        is_below_ = &Function::do_is_below_never_;
        if (is_max_open)
        {
            is_above_ = &Function::do_is_above_open_max_;
            is_valid_ = &Function::do_is_valid_open_max_;
        }
        else
        {
            is_above_ = &Function::do_is_above_close_max_;
            is_valid_ = &Function::do_is_valid_close_max_;
        }
    }
    else // Within an Interval
    {
        if (is_min_open && is_max_open)
        {
            is_below_ = &Function::do_is_below_open_min_;
            is_above_ = &Function::do_is_above_open_max_;
            is_valid_ = &Function::do_is_valid_both_open_;
        }
        else if (is_min_open && !is_max_open)
        {
            is_below_ = &Function::do_is_below_open_min_;
            is_above_ = &Function::do_is_above_close_max_;
            is_valid_ = &Function::do_is_valid_open_min_close_max_;
        }
        else if (!is_min_open && is_max_open)
        {
            is_below_ = &Function::do_is_below_close_min_;
            is_above_ = &Function::do_is_above_open_max_;
            is_valid_ = &Function::do_is_valid_close_min_open_max_;
        }
        else
        {
            is_below_ = &Function::do_is_below_close_min_;
            is_above_ = &Function::do_is_above_close_max_;
            is_valid_ = &Function::do_is_valid_both_close_;
        }

        // check for consistency
        bool wrong_open = (is_valid_ != &Function::do_is_valid_both_close_) && (range_lower_.get() >= range_upper_.get());
        bool wrong_close = (is_valid_ == &Function::do_is_valid_both_close_) && (range_lower_.get() > range_upper_.get());
        if (wrong_open || wrong_close)
            throw std::invalid_argument(
                "Generated an inconsistent range " + str_range() + " for " + err_name() + ".");
    }
}

/// Return formatted range.
std::string Function::str_range() const
{
    const std::string str_min = hsp::str::format(range_lower_.get());
    const std::string str_max = hsp::str::format(range_upper_.get());

    if (is_valid_ == &Function::do_is_valid_always_)
        return "(-Inf, Inf)";

    if (is_valid_ == &Function::do_is_valid_open_min_)
        return "(" + str_min + ", Inf)";

    if (is_valid_ == &Function::do_is_valid_close_min_)
        return "[" + str_min + ", Inf)";

    if (is_valid_ == &Function::do_is_valid_open_max_)
        return "(-Inf, " + str_max + ")";

    if (is_valid_ == &Function::do_is_valid_close_max_)
        return "(-Inf, " + str_max + "]";

    if (is_valid_ == &Function::do_is_valid_both_open_)
        return "(" + str_min + ", " + str_max + ")";

    if (is_valid_ == &Function::do_is_valid_open_min_close_max_)
        return "(" + str_min + ", " + str_max + "]";

    if (is_valid_ == &Function::do_is_valid_close_min_open_max_)
        return "[" + str_min + ", " + str_max + ")";

    // (is_valid_ == &Function::do_is_valid_both_close_)
    return "[" + str_min + ", " + str_max + "]";
}

/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////  Function : Arguments  ///////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

/// Assign attribute to function.
void Function::init_attr(const std::string &key, const Attribute &attr)
{
    // assert that function attribute are implicitely callable
    if (attr.is_function() && !attr.get_function()->is_implicit_callable())
        throw std::invalid_argument(
            "Impossible to associate " + attr.err_name() + " as argument '" + key + "' of " + err_name() + ".\n" +
            "Argument 'x' of " + attr.err_name() + " cannot be missing.");

    // distinguish between range keys and standard keys
    if (key == hsp::params::MIN || key == hsp::params::MAX)
    {
        // skip empty constants for range
        // NOTE: Variable overwrites 'is_constant' method,
        // thus empty does not always imply constant
        if (attr.is_empty() && attr.is_constant())
            return;

        if (name_ == "uniform")
            find_arg_(key).transfer(attr);
        else if (key == hsp::params::MIN)
            range_lower_.transfer(attr);
        else // key == hsp::params::MAX
            range_upper_.transfer(attr);
    }
    else
        find_arg_(key).transfer(attr);
}

/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////  Function : discretization  //////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

/// Select functional method for discretization.
void Function::init_discr_helper_(unsigned short discr_type)
{
    if (discr_type == hsp::discr::DEFINED)
    {
        function_ = &Function::do_function_discrete_;
        log_func_ = &Function::do_logfunc_discrete_;
        primitive_ = &Function::do_primitive_discrete_;
        log_ratio_ = &Function::do_logratio_discrete_;
    }
    else if (discr_type == hsp::discr::CENTER)
    {
        function_ = &Function::do_function_discr_center_;
        log_func_ = &Function::do_logfunc_discr_center_;
        primitive_ = &Function::do_primitive_discr_center_;
        log_ratio_ = &Function::do_logratio_discr_center_;
    }
    else if (discr_type == hsp::discr::LEFT)
    {
        function_ = &Function::do_function_discr_left_;
        log_func_ = &Function::do_logfunc_discr_left_;
        primitive_ = &Function::do_primitive_discr_left_;
        log_ratio_ = &Function::do_logratio_discr_left_;
    }
    else // no discretization (default)
    {
        function_ = &Function::raw_function_;
        log_func_ = &Function::raw_logfunc_;
        primitive_ = &Function::raw_primitive_;
        log_ratio_ = &Function::raw_logratio_;
    }
}

/// Returns the case of the discretization.
unsigned short Function::which_discrete_case_() const
{
    if (function_ == &Function::do_function_discrete_)
        return hsp::discr::DEFINED;
    if (function_ == &Function::do_function_discr_center_)
        return hsp::discr::CENTER;
    if (function_ == &Function::do_function_discr_left_)
        return hsp::discr::LEFT;
    return hsp::discr::NONE;
}

/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////  Function : show  ///////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

/// Return formatted function summary.
std::string Function::str_summary(bool show_x_arg) const
{
    std::string output = name_ + "(";
    if (show_x_arg)
        output += x_arg_.get_name() + "=" + hsp::str::format(x_arg_.get());
    for (size_t c = 0; c < arguments_.size(); ++c)
    {
        if (c > 0 || show_x_arg)
            output += ", ";
        output += arguments_[c].get_name() + "=" + hsp::str::format(arguments_[c].get());
    }
    return output + ") : " + str_range();
}

/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////  FunctionRegistry : entries  //////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

// FIXME:
// verify if is_discrete should include Function::is_discrete() || (...)

/// Minus Function.
class Minus_Function : public Function
{
public:
    Minus_Function() : Function("minus") {}
    void init_truncation() override final {}
    bool is_discrete() const override final { return x_arg_.is_discrete(); }

protected:
    double raw_function_(double x) const override final { return -x; }
};

/// Subtraction function.
class Subtraction_Function : public Function
{
public:
    Subtraction_Function() : Function("sub", {"x2"}, {0.}) {}
    void init_truncation() override final {}
    bool is_discrete() const override final { return get_arg(0).is_discrete() && x_arg_.is_discrete(); }

protected:
    double raw_function_(double x) const override final { return x - arg(0); }
};

/// Division Function.
class Division_Function : public Function
{
public:
    Division_Function() : Function("div", {"x2"}, {1}) {}
    void init_truncation() override final {}

protected:
    double raw_function_(double x) const override final { return x / arg(0); }
};

/// Power Function.
class Power_Function : public Function
{
public:
    Power_Function() : Function("pow", {"x2"}, {1.}) {}
    void init_truncation() override final { init_range_helper_(); }
    bool is_discrete() const override final { return Function::is_discrete() || (get_arg(0).is_discrete() && x_arg_.is_discrete()); }

protected:
    double raw_function_(double x) const override final { return std::pow(x, arg(0)); }
    double raw_logfunc_(double x) const override final { return arg(0) * std::log(x); }
    double raw_primitive_(double x) const override final
    {
        if (std::abs(arg(0) + 1.) < 1.e-10)
            return std::log(std::abs(x));
        return std::pow(x, arg(0) + 1.) / (arg(0) + 1.);
    }
    double raw_logratio_(double xo, double xn) const override final { return arg(0) * std::log(xn / xo); }
};

/// SquareRoot Function.
class SquareRoot_Function : public Function
{
public:
    SquareRoot_Function() : Function("sqrt", {}, {},
                                     0., false, hsp::na::DOUBLE, false) {}
    void init_truncation() override final { init_range_min_(); }

protected:
    double raw_function_(double x) const override final { return std::sqrt(x); }
    double raw_logfunc_(double x) const override final { return 0.5 * std::log(x); }
    double raw_primitive_(double x) const override final { return std::pow(x, 1.5) / 1.5; }
    double raw_logratio_(double xo, double xn) const override final { return 0.5 * std::log(xn / xo); }
};

/// Exponential Function.
class Exponential_Function : public Function
{
public:
    Exponential_Function() : Function("exp") {}
    void init_truncation() override final { init_range_helper_(); }

protected:
    double raw_function_(double x) const override final { return std::exp(x); }
    double raw_logfunc_(double x) const override final { return x; }
    double raw_primitive_(double x) const override final { return std::exp(x); }
    double raw_logratio_(double xo, double xn) const override final { return xn - xo; }
};

/// Logarithm Function.
class Logarithm_Function : public Function
{
public:
    Logarithm_Function() : Function("log", {}, {},
                                    0., true, hsp::na::DOUBLE, false) {}
    void init_truncation() override final { init_range_min_(0., true); }

protected:
    double raw_function_(double x) const override final { return std::log(x); }
    double raw_logfunc_(double x) const override final { return std::log(std::log(x)); }
    double raw_primitive_(double x) const override final { return x * std::log(x) - x; }
    double raw_logratio_(double xo, double xn) const override final { return std::log(std::log(xn) / std::log(xo)); }
};

/// Sine Function.
class Sine_Function : public Function
{
public:
    Sine_Function() : Function("sin") {}
    void init_truncation() override final { init_range_helper_(); }

protected:
    double raw_function_(double x) const override final { return std::sin(x); }
    double raw_logfunc_(double x) const override final { return std::log(std::sin(x)); }
    double raw_primitive_(double x) const override final { return -std::cos(x); }
    double raw_logratio_(double xo, double xn) const override final { return std::log(std::sin(xn) / std::sin(xo)); }
};

/// Cosine Function.
class Cosine_Function : public Function
{
public:
    Cosine_Function() : Function("cos") {}
    void init_truncation() override final { init_range_helper_(); }

protected:
    double raw_function_(double x) const override final { return std::cos(x); }
    double raw_logfunc_(double x) const override final { return std::log(std::cos(x)); }
    double raw_primitive_(double x) const override final { return std::sin(x); }
    double raw_logratio_(double xo, double xn) const override final { return std::log(std::cos(xn) / std::cos(xo)); }
};

/// Tangent Function.
class Tangent_Function : public Function
{
public:
    Tangent_Function() : Function("tan") {}
    void init_truncation() override final { init_range_helper_(); }

protected:
    double raw_function_(double x) const override final { return std::tan(x); }
    double raw_logfunc_(double x) const override final { return std::log(std::tan(x)); }
    double raw_logratio_(double xo, double xn) const override final { return std::log(std::tan(xn) / std::tan(xo)); }
};

/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////  Variadic_Function : class  ///////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

/// Auxialiary class to deal with variadic functions.
class Variadic_Function : public Function
{
public:
    Variadic_Function(const std::string &name = hsp::na::STRING) : Function(name) {}
    virtual ~Variadic_Function() = default;

    // UPDATE: this might change for future variadic functions
    void init_truncation() override final {}
    bool is_discrete() const override final
    {
        if (!x_arg_.is_discrete())
            return false;
        for (const Attribute &attr : arguments_)
            if (!attr.is_discrete())
                return false;
        return true;
    }

protected:
    /// Variadic function version of the find argument method.
    Attribute &find_arg_(const std::string &arg_name) override final
    {
        if (arg_name == hsp::params::X)
            return x_arg_;

        for (Attribute &attr : arguments_)
            if (attr.get_name() == arg_name)
                return attr;

        // if attribute doesn't exist at initialization, create it
        arguments_.emplace_back();
        return arguments_.back();
    }
};

#define MAKE_VARIADIC_ALGEBRAIC(ClassName, FunctionName, Symbol) \
    class ClassName : public Variadic_Function                   \
    {                                                            \
    public:                                                      \
        ClassName() : Variadic_Function(FunctionName) {}         \
                                                                 \
    protected:                                                   \
        double raw_function_(double x) const override final      \
        {                                                        \
            for (const Attribute &attr : arguments_)             \
                x Symbol attr.get();                             \
            return x;                                            \
        }                                                        \
    };
MAKE_VARIADIC_ALGEBRAIC(Summation_Function, "sum", +=)
MAKE_VARIADIC_ALGEBRAIC(Product_Function, "prod", *=)

#define MAKE_VARIADIC_MIN_MAX(ClassName, FunctionName, Symbol) \
    class ClassName : public Variadic_Function                 \
    {                                                          \
    public:                                                    \
        ClassName() : Variadic_Function(FunctionName) {}       \
                                                               \
    protected:                                                 \
        double raw_function_(double x) const override final    \
        {                                                      \
            for (const Attribute &attr : arguments_)           \
            {                                                  \
                const double new_val = attr.get();             \
                if (new_val Symbol x)                          \
                    x = new_val;                               \
            }                                                  \
            return x;                                          \
        }                                                      \
    };
MAKE_VARIADIC_MIN_MAX(Minimum_Function, "min", <)
MAKE_VARIADIC_MIN_MAX(Maximum_Function, "max", >)

/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////  FunctionRegistry : getter  //////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

str_vec FunctionRegistry::get_options() const
{
    str_vec keys;
    keys.reserve(factory_map_.size());
    for (const auto &pair : factory_map_)
        keys.push_back(pair.first);
    return keys;
}

std::unique_ptr<Function> FunctionRegistry::get(const std::string &opt, bool is_primitive) const
{
    auto output = factory_map_.find(opt)->second();
    if (is_primitive)
        output->set_eval_primitive();
    else
        output->set_eval_call();
    return output;
}

std::unique_ptr<Function> compose_func_registries(const std::string &opt,
                                                  const std::string &name,
                                                  const std::vector<const FunctionRegistry *> &registries,
                                                  bool is_prim)
{
    str_vec options = {};

    if (registries.empty())
        throw std::invalid_argument(
            "No registries were provided for '" + name + "'.");

    if (hsp::na::is(opt))
        return std::make_unique<Function>();

    for (const auto &r : registries)
    {
        if (r->get_options() % opt)
            return r->get(opt, is_prim);
        options += r->get_options();
    }

    throw std::invalid_argument(
        "Unknown " + name + " '" + opt + "'." + hsp::str::find_closest(opt, options));
}

/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////  Register functions  //////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

#define IMPLEMENT_HOUSEPI_FUNCTION_REGISTRY(ClassName) \
    ClassName &ClassName::inst()                       \
    {                                                  \
        static ClassName _inst;                        \
        return _inst;                                  \
    }

IMPLEMENT_HOUSEPI_FUNCTION_REGISTRY(ArithmeticRegistry)
REGISTER_HOUSEPI_FUNCTION(ArithmeticRegistry, Minimum_Function)
REGISTER_HOUSEPI_FUNCTION(ArithmeticRegistry, Maximum_Function)
REGISTER_HOUSEPI_FUNCTION(ArithmeticRegistry, Minus_Function)
REGISTER_HOUSEPI_FUNCTION(ArithmeticRegistry, Subtraction_Function)
REGISTER_HOUSEPI_FUNCTION(ArithmeticRegistry, Division_Function)
REGISTER_HOUSEPI_FUNCTION(ArithmeticRegistry, Summation_Function)
REGISTER_HOUSEPI_FUNCTION(ArithmeticRegistry, Product_Function)

IMPLEMENT_HOUSEPI_FUNCTION_REGISTRY(TranscendentalRegistry)
REGISTER_HOUSEPI_FUNCTION(TranscendentalRegistry, Power_Function)
REGISTER_HOUSEPI_FUNCTION(TranscendentalRegistry, SquareRoot_Function)
REGISTER_HOUSEPI_FUNCTION(TranscendentalRegistry, Exponential_Function)
REGISTER_HOUSEPI_FUNCTION(TranscendentalRegistry, Logarithm_Function)
REGISTER_HOUSEPI_FUNCTION(TranscendentalRegistry, Sine_Function)
REGISTER_HOUSEPI_FUNCTION(TranscendentalRegistry, Cosine_Function)
REGISTER_HOUSEPI_FUNCTION(TranscendentalRegistry, Tangent_Function)

IMPLEMENT_HOUSEPI_FUNCTION_REGISTRY(UserFuncRegistry)
