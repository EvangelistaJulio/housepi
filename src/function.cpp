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
#include "tools/string.h"
#include "tools/types.h"

#include "function.h"
#include "registry.h"

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
    {
        output += std::string(x_arg_.name());
        output += "=";
        output += hsp::str::format(x_arg_.get());
    }
    for (size_t c = 0; c < arguments_.size(); ++c)
    {
        if (c > 0 || show_x_arg)
            output += ", ";
        output += std::string(arguments_[c].name());
        output += "=";
        output += hsp::str::format(arguments_[c].get());
    }
    return output + ") : " + str_range();
}

/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////// *** FunctionRegistry : class *** //////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

std::vector<std::string> FunctionRegistry::get_options() const
{
    std::vector<std::string> entries;
    entries.reserve(factory_map_.size());
    for (const auto &[key, val] : factory_map_)
        entries.push_back(key);
    return entries;
}

std::unique_ptr<Function> FunctionRegistry::create(std::string_view opt, bool is_primitive) const
{
    auto it = factory_map_.find(opt);
    if (it == factory_map_.end())
        return nullptr;

    auto output = it->second();
    if (is_primitive)
        output->set_eval_primitive();
    else
        output->set_eval_call();
    return output;
}

std::unique_ptr<Function> compose_func_registries(std::string_view opt,
                                                  const std::string &err_type,
                                                  const std::vector<std::string> &registries,
                                                  bool is_prim)
{
    if (registries.empty())
        throw std::invalid_argument(
            "No registries were provided for '" + err_type + "'.");

    if (hsp::na::is(opt))
        return nullptr;

    std::vector<const FunctionRegistry *> registry_ptrs = {};
    for (const std::string &r : registries)
        if (r == "user")
            registry_ptrs.push_back(&UserFuncRegistry::instance());
        else if (r == "transcendental")
            registry_ptrs.push_back(&TranscendentalRegistry::instance());
        else if (r == "arithmetic")
            registry_ptrs.push_back(&ArithmeticRegistry::instance());
        else if (r == "continuos")
            registry_ptrs.push_back(&ContinuosDistrRegistry::instance());
        else if (r == "discrete")
            registry_ptrs.push_back(&DiscreteDistrRegistry::instance());
        else if (r == "infecprofile")
            registry_ptrs.push_back(&InfecProfileRegistry::instance());
        else
            throw std::invalid_argument("Unknown " + err_type + " registry '" + r + "'.");

    std::vector<std::string> options = {};
    for (const auto &ptr : registry_ptrs)
    {
        if (ptr->count(opt))
            return ptr->create(opt, is_prim);
        const std::vector<std::string> &y = ptr->get_options();
        options.reserve(options.size() + y.size());
        options.insert(options.end(), y.begin(), y.end());
    }

    throw std::invalid_argument("Unknown " + err_type + " '" + std::string(opt) + "'." +
                                hsp::str::find_closest(opt, options));
}

/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////  Attribute : constructor ///////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

/// Empty attribute constructor.
Attribute::Attribute()
{
    init_numeric(hsp::na::DOUBLE);
    name("", "numeric");
}
Attribute::Attribute(double x)
{
    init_numeric(x);
    name("", "numeric");
}
Attribute::Attribute(const Attribute *y)
{
    init_pointer(y);
    name("", "pointer");
}
Attribute::Attribute(const std::shared_ptr<Function> f)
{
    init_function(f);
    name(std::string(f->name()), "function");
}
Attribute::Attribute(std::string_view s, bool is_primitive)
{
    init_function(s, is_primitive);
    name(std::string(pFunc_->name()), "function");
}
/// Attribute deconstructor.
Attribute::~Attribute() = default;

/// Wheter it's a number.
bool Attribute::is_numeric() const noexcept { return !hsp::na::is(value_); }

std::string Attribute::str_summary() const
{
    std::string sms = "";
    if (!is_function() || pFunc_->is_callable_implicitly())
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
bool Attribute::is_constant() const noexcept
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
bool Attribute::is_discrete() const noexcept
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

/// Set attribute name. Eventually also err type.
void Attribute::name(const std::string &name, const std::string &err_type)
{
    name_ = name;
    if (!hsp::na::is(err_type))
        err_type_ = err_type;
}

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
/// Initialize attribute as function given the function name.
void Attribute::init_function(std::string_view s, bool is_primitive)
{
    init_function(compose_func_registries(
        s, "function",
        {"user", "transcendental", "arithmetic",
         "continuos", "discrete", "infecprofile"},
        is_primitive));
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
/// NOTE: to access arguments by index, I must provide separate name/value.
Function::Function(const std::string &name,
                   const std::vector<std::pair<std::string, double>> &args,
                   double domain_min, bool is_min_open,
                   double domain_max, bool is_max_open)
{
    // initialize default arguments
    x_arg_.name(hsp::params::X);
    x_arg_.set(hsp::na::DOUBLE);
    loc_arg_.name(hsp::params::LOC);
    loc_arg_.set(hsp::na::DOUBLE);

    // name of the function
    name_ = name;

    // allocate attribute for remaining arguments
    size_t n_args = args.size();
    for (const auto &[key, val] : args)
        if (key == hsp::params::X || key == hsp::params::LOC)
            --n_args;
    arguments_ = std::vector<Attribute>(n_args);

    size_t counter = 0;
    for (const auto &[key, val] : args)
        if (key == hsp::params::X)
        {
            x_arg_.set(val);
        }
        else if (key == hsp::params::LOC)
        {
            loc_arg_.set(val);
        }
        else
        {
            arguments_[counter].name(key);
            arguments_[counter++].set(val);
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
Attribute &Function::find_arg_(std::string_view arg_name)
{
    if (hsp::na::is(arg_name))
        throw std::invalid_argument(
            "Provided an unnamed argument to " + err_name() + ".");

    if (arg_name == hsp::params::X)
        return x_arg_;
    if (arg_name == hsp::params::LOC)
        return loc_arg_;

    for (Attribute &attr : arguments_)
        if (attr.name() == arg_name)
            return attr;

    // else format the error showing options
    std::vector<std::string> options = {hsp::params::X, hsp::params::LOC};
    for (const Attribute &attr : arguments_)
        options.emplace_back(attr.name());
    throw std::invalid_argument(
        "Unknown argument '" + std::string(arg_name) + "' for " + err_name() + "." +
        hsp::str::find_closest(arg_name, options));
}

/// Wheter return value is always constant.
bool Function::is_constant() const noexcept
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
    if (!lower_arg_.is_constant())
        return false;
    if (!upper_arg_.is_constant())
        return false;

    // otherwise constant
    return true;
}

/// Wheter all arguments are not empty.
bool Function::is_callable_explicitly() const noexcept
{
    for (const Attribute &attr : arguments_)
        if (attr.is_empty())
            return false;
    return true;
}

/// Collect names of missing arguments.
std::vector<std::string> Function::get_missing_argument() const
{
    std::vector<std::string> undefined_args = {};
    for (const Attribute &attr : arguments_)
        if (attr.is_empty())
            undefined_args.emplace_back(attr.name());
    return undefined_args;
}

/// Throw exception if not callable explicitly.
void Function::assert_is_callable_explicitly() const
{
    const std::vector<std::string> undefined_args = get_missing_argument();
    if (!undefined_args.empty())
        throw std::invalid_argument(
            "Impossible to initialize " + err_name() + ": missing non default argument" +
            (undefined_args.size() > 1 ? "s " : " ") + hsp::str::format(undefined_args));
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
    if (lower_arg_.is_empty() && upper_arg_.is_empty()) // Unbounded
    {
        is_below_ = &Function::do_is_below_never_;
        is_above_ = &Function::do_is_above_never_;
        is_valid_ = &Function::do_is_valid_always_;
    }
    else if (!lower_arg_.is_empty() && upper_arg_.is_empty()) // Lower Bound
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
    else if (lower_arg_.is_empty() && !upper_arg_.is_empty()) // Upper Bound
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
        bool wrong_open = (is_valid_ != &Function::do_is_valid_both_close_) && (lower_arg_.get() >= upper_arg_.get());
        bool wrong_close = (is_valid_ == &Function::do_is_valid_both_close_) && (lower_arg_.get() > upper_arg_.get());
        if (wrong_open || wrong_close)
            throw std::invalid_argument(
                "Generated an inconsistent range " + str_range() + " for " + err_name() + ".");
    }
}

/// Return formatted range.
std::string Function::str_range() const
{
    const std::string str_min = hsp::str::format(lower_arg_.get());
    const std::string str_max = hsp::str::format(upper_arg_.get());

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
void Function::init_attr(std::string_view key, const Attribute &attr)
{
    // assert that function attribute are implicitely callable
    if (attr.is_function())
    {
        const std::vector<std::string> undefined_args = attr.get_function()->get_missing_argument();
        if (!undefined_args.empty())
            throw std::invalid_argument(
                "Impossible to associate " + attr.err_name() + " as argument '" + std::string(key) + "' of " + err_name() + ".\n" +
                "Non default argument" + (undefined_args.size() > 1 ? "s " : " ") + hsp::str::format(undefined_args) +
                " of " + attr.err_name() + " cannot be missing.");

        if (!attr.get_function()->is_callable_implicitly())
            throw std::invalid_argument(
                "Impossible to associate " + attr.err_name() + " as argument '" + std::string(key) + "' of " + err_name() + ".\n" +
                "Argument 'x' of " + attr.err_name() + " cannot be missing.");
    }

    // distinguish between range keys and standard keys
    if (key == hsp::params::MIN || key == hsp::params::MAX)
    {
        // skip empty constants for range
        // NOTE: Variable overwrites 'is_constant' method,
        // thus empty does not always imply constant
        if (attr.is_empty() && attr.is_constant())
            return;

        if (name() == "uniform")
            find_arg_(key).transfer(attr);
        else if (key == hsp::params::MIN)
            lower_arg_.transfer(attr);
        else // key == hsp::params::MAX
            upper_arg_.transfer(attr);
    }
    else
        find_arg_(key).transfer(attr);

    // Apply location argument shift.
    if (key == hsp::params::LOC)
        delocate_ = &Function::do_delocate_;
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

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////  *** Function : implementations ***  /////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

namespace
{
    using std::log;

    /// Minus Function.
    class Minus_ final : public Function
    {
    public:
        Minus_() : Function("minus") {}
        void init_truncation() override final {}
        bool is_discrete() const noexcept override final { return x_arg_.is_discrete(); }

    protected:
        double raw_function_(double x) const override final { return -x; }
    };

    /// Subtraction function.
    class Subtraction_ final : public Function
    {
    public:
        Subtraction_() : Function("sub", {{"x2", 0.}}) {}
        void init_truncation() override final {}
        bool is_discrete() const noexcept override final { return get_arg(0).is_discrete() && x_arg_.is_discrete(); }

    protected:
        double raw_function_(double x) const override final { return x - arg(0); }
    };

    /// Division Function.
    class Division_ final : public Function
    {
    public:
        Division_() : Function("div", {{"x2", 1}}) {}
        void init_truncation() override final {}

    protected:
        double raw_function_(double x) const override final { return x / arg(0); }
    };

    /// Power Function.
    class Power_ final : public Function
    {
    public:
        Power_() : Function("pow", {{"x2", 1.}}) {}
        void init_truncation() override final { init_range_helper_(); }
        bool is_discrete() const noexcept override final { return Function::is_discrete() || (get_arg(0).is_discrete() && x_arg_.is_discrete()); }

    protected:
        double raw_function_(double x) const override final { return std::pow(x, arg(0)); }
        double raw_logfunc_(double x) const override final { return arg(0) * log(x); }
        double raw_primitive_(double x) const override final
        {
            if (std::abs(arg(0) + 1.) < 1.e-10)
                return log(std::abs(x));
            return std::pow(x, arg(0) + 1.) / (arg(0) + 1.);
        }
        double raw_logratio_(double xo, double xn) const override final { return arg(0) * log(xn / xo); }
    };

    /// SquareRoot Function.
    class SquareRoot_ final : public Function
    {
    public:
        SquareRoot_() : Function("sqrt", {},
                                 0., false, hsp::na::DOUBLE, false) {}
        void init_truncation() override final { init_range_min_(); }

    protected:
        double raw_function_(double x) const override final { return std::sqrt(x); }
        double raw_logfunc_(double x) const override final { return 0.5 * log(x); }
        double raw_primitive_(double x) const override final { return std::pow(x, 1.5) / 1.5; }
        double raw_logratio_(double xo, double xn) const override final { return 0.5 * log(xn / xo); }
    };

    /// Exponential Function.
    class Exponential_ final : public Function
    {
    public:
        Exponential_() : Function("exp") {}
        void init_truncation() override final { init_range_helper_(); }

    protected:
        double raw_function_(double x) const override final { return std::exp(x); }
        double raw_logfunc_(double x) const override final { return x; }
        double raw_primitive_(double x) const override final { return std::exp(x); }
        double raw_logratio_(double xo, double xn) const override final { return xn - xo; }
    };

    /// Logarithm Function.
    class Logarithm_ final : public Function
    {
    public:
        Logarithm_() : Function("log", {},
                                0., true, hsp::na::DOUBLE, false) {}
        void init_truncation() override final { init_range_min_(0., true); }

    protected:
        double raw_function_(double x) const override final { return log(x); }
        double raw_logfunc_(double x) const override final { return log(log(x)); }
        double raw_primitive_(double x) const override final { return x * log(x) - x; }
        double raw_logratio_(double xo, double xn) const override final { return log(log(xn) / log(xo)); }
    };

    /// Sine Function.
    class Sine_ final : public Function
    {
    public:
        Sine_() : Function("sin") {}
        void init_truncation() override final { init_range_helper_(); }

    protected:
        double raw_function_(double x) const override final { return std::sin(x); }
        double raw_logfunc_(double x) const override final { return log(std::sin(x)); }
        double raw_primitive_(double x) const override final { return -std::cos(x); }
        double raw_logratio_(double xo, double xn) const override final { return log(std::sin(xn) / std::sin(xo)); }
    };

    /// Cosine Function.
    class Cosine_ final : public Function
    {
    public:
        Cosine_() : Function("cos") {}
        void init_truncation() override final { init_range_helper_(); }

    protected:
        double raw_function_(double x) const override final { return std::cos(x); }
        double raw_logfunc_(double x) const override final { return log(std::cos(x)); }
        double raw_primitive_(double x) const override final { return std::sin(x); }
        double raw_logratio_(double xo, double xn) const override final { return log(std::cos(xn) / std::cos(xo)); }
    };

    /// Tangent Function.
    class Tangent_ final : public Function
    {
    public:
        Tangent_() : Function("tan") {}
        void init_truncation() override final { init_range_helper_(); }

    protected:
        double raw_function_(double x) const override final { return std::tan(x); }
        double raw_logfunc_(double x) const override final { return log(std::tan(x)); }
        double raw_logratio_(double xo, double xn) const override final { return log(std::tan(xn) / std::tan(xo)); }
    };

    /////////////////////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////  VariadicFunction : class  ///////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////

    /// Auxialiary class to deal with variadic functions.
    class VariadicFunction : public Function
    {
    public:
        VariadicFunction(const std::string &name = "") : Function(name) {}
        virtual ~VariadicFunction() = default;

        // UPDATE: this might change for future variadic functions
        void init_truncation() override final {}
        bool is_discrete() const noexcept override final
        {
            if (!x_arg_.is_discrete())
                return false;
            for (const Attribute &attr : arguments_)
                if (!attr.is_discrete())
                    return false;
            return true;
        }

    protected:
        /// Override the find argument method at assigment.
        /// If attribute doesn't exist at initialization, create it.
        Attribute &find_arg_(std::string_view arg_name) override final
        {
            if (arg_name == hsp::params::X)
                return x_arg_;
            if (arg_name == hsp::params::LOC)
                return loc_arg_;
            for (Attribute &attr : arguments_)
                if (!hsp::na::is(attr.name()) && attr.name() == arg_name)
                    return attr;

            arguments_.emplace_back();
            return arguments_.back();
        }
    };

#define MAKE_VARIADIC_ALGEBRAIC(ClassName, FunctionName, Symbol) \
    class ClassName final : public VariadicFunction              \
    {                                                            \
    public:                                                      \
        ClassName() : VariadicFunction(FunctionName) {}          \
                                                                 \
    protected:                                                   \
        double raw_function_(double x) const override final      \
        {                                                        \
            for (const Attribute &attr : arguments_)             \
                x Symbol attr.get();                             \
            return x;                                            \
        }                                                        \
    };
    MAKE_VARIADIC_ALGEBRAIC(Summation_, "sum", +=)
    MAKE_VARIADIC_ALGEBRAIC(Product_, "prod", *=)

#define MAKE_VARIADIC_MIN_MAX(ClassName, FunctionName, Symbol) \
    class ClassName final : public VariadicFunction            \
    {                                                          \
    public:                                                    \
        ClassName() : VariadicFunction(FunctionName) {}        \
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
    MAKE_VARIADIC_MIN_MAX(Minimum_, "min", <)
    MAKE_VARIADIC_MIN_MAX(Maximum_, "max", >)
}

/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////  Register functions  //////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

REGISTER_HOUSEPI_FUNCTION(Minimum_, ArithmeticRegistry)
REGISTER_HOUSEPI_FUNCTION(Maximum_, ArithmeticRegistry)
REGISTER_HOUSEPI_FUNCTION(Minus_, ArithmeticRegistry)
REGISTER_HOUSEPI_FUNCTION(Subtraction_, ArithmeticRegistry)
REGISTER_HOUSEPI_FUNCTION(Division_, ArithmeticRegistry)
REGISTER_HOUSEPI_FUNCTION(Summation_, ArithmeticRegistry)
REGISTER_HOUSEPI_FUNCTION(Product_, ArithmeticRegistry)

REGISTER_HOUSEPI_FUNCTION(Power_, TranscendentalRegistry)
REGISTER_HOUSEPI_FUNCTION(SquareRoot_, TranscendentalRegistry)
REGISTER_HOUSEPI_FUNCTION(Exponential_, TranscendentalRegistry)
REGISTER_HOUSEPI_FUNCTION(Logarithm_, TranscendentalRegistry)
REGISTER_HOUSEPI_FUNCTION(Sine_, TranscendentalRegistry)
REGISTER_HOUSEPI_FUNCTION(Cosine_, TranscendentalRegistry)
REGISTER_HOUSEPI_FUNCTION(Tangent_, TranscendentalRegistry)

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////  Attribute : extras  /////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/// Productory of pointed attributes values.
Attribute::Attribute(const std::vector<const Attribute *> other_ptrs)
{
    if (other_ptrs.size() == 1)
    {
        init_pointer(other_ptrs[0]);
        return;
    }

    auto fun = std::make_shared<Product_>();
    // initialize attributes to pointers
    fun->init_attr(hsp::params::X, Attribute(other_ptrs[0]));
    for (size_t i = 1; i < other_ptrs.size(); ++i)
        fun->init_attr("", Attribute(other_ptrs[i]));
    init_function(fun);
}
