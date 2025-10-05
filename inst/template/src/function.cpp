/*
 *        __
 *       / /_  ____  __  __ ____ __  ____  __
 *      / __ \/ __ \/ / / / ___/ _ \/ __ \/ /
 *     / / / / /_/ / /_/ (__  )  __/ /_/ / /
 *    /_/ /_/\____/\____/____/\___/ ____/_/
 *                               /_/
 *
 * Copyright (C) 2025 Francesco Camaglia, Institut Pasteur
 *
 * function.cpp
 *
 */

#include "function.h"

/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////  Attribute : constructor ///////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

Attribute::Attribute(double x)
{
    name(HSP_NA_STRING, "numeric");
    init_numeric(x);
}
Attribute::Attribute(observ_dbl_ptr y)
{
    name(HSP_NA_STRING, "pointer");
    init_pointer(y);
}
Attribute::Attribute(const std::shared_ptr<Function> f)
{
    name(f->name(), "function");
    init_function(f);
}
Attribute::Attribute(const std::string &s, bool is_primitive)
{
    name(s, "function");
    init_function(s, is_primitive);
}

/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////  Attribute : getter //////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

std::string Attribute::what() const
{
    std::string sms = is_constant() ? "constant" : "variable";
    if (is_numeric())
        sms += " numeric";
    else if (is_pointer())
        sms += " pointer";
    else if (is_function())
        sms += " function";
    else
        sms += " empty";
    return sms;
}

/// Get attribute value from external function.
/// note: this can't be moved in the class because of forward declaration of Function
double Attribute::do_get_from_function_() const { return pFunc_->evaluate(); }

/// Wheter it's constant.
/// NOTE: by construction convention, a pointer to constant must be initialized to numeric instead.
/// NOTE: if empty, then it's constant by convention.
bool Attribute::is_constant() const
{
    if (is_numeric())
        return true;
    else if (is_pointer())
        return false;
    else if (is_function())
        return pFunc_->is_constant();
    else // NOTE: is_empty(), constant by convention
        return true;
}

/// Wheter it's discrete.
/// NOTE: if empty, then it's discrete by convention.
bool Attribute::is_discrete() const
{
    if (is_numeric())
        return value_ == std::floor(value_);
    else if (is_pointer())
        return false; // FIXME: here it becomes problematic..
    else if (is_function())
        return pFunc_->is_discrete();
    else
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
    ptr_to_value_ = nullptr;
    pFunc_ = nullptr;
}
/// Initialize attribute as pointer.
void Attribute::init_pointer(observ_dbl_ptr y)
{
    if (y == nullptr)
        throw std::invalid_argument(
            "Attempt to assign empty pointer to " + err_name() + ".");

    get_ = &Attribute::do_get_from_pointer_;
    get_ptr_ = &Attribute::do_get_pointer_from_pointer_;
    do_set_ = &Attribute::do_set_error_;
    value_ = HSP_NA_DOUBLE;
    ptr_to_value_ = y;
    pFunc_ = nullptr;
}
/// Initialize attribute as function.
void Attribute::init_function(const std::shared_ptr<Function> f)
{
    // check the function is different from the current one
    if (f == nullptr)
        throw std::invalid_argument(
            "Attempt to assign empty function to " + err_name() + ".");
    else if (f == pFunc_)
        throw std::invalid_argument(
            "Cannot assign " + err_name() + " as its own argument.");

    get_ = &Attribute::do_get_from_function_;
    get_ptr_ = &Attribute::do_get_pointer_error_;
    do_set_ = &Attribute::do_set_error_;
    value_ = HSP_NA_DOUBLE;
    ptr_to_value_ = nullptr;
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
    // name of the function
    this->name_ = name;

    // number of arguments
    if (arg_name.size() != arg_value.size())
        tools::error(
            "Inconsistent arguments for " + err_name() + "." +
            "Found " + std::to_string(arg_name.size()) + " default names but " +
            std::to_string(arg_value.size()) + " values.");

    // initialize default 'x' argument
    x_arg_.name(HSP_PARAM_X);
    x_arg_.set(HSP_NA_DOUBLE);

    // allocate attribute for remaining arguments
    size_t n_args = arg_name.size() - (arg_name % HSP_PARAM_X);
    arguments_ = std::vector<Attribute>(n_args);

    size_t c = 0;
    for (size_t k = 0; k < arg_name.size(); ++k)
        if (arg_name[k] == HSP_PARAM_X)
            x_arg_.set(arg_value[k]);
        else
        {
            arguments_[c].name(arg_name[k]);
            arguments_[c++].set(arg_value[k]);
        }

    // Assign Domain
    // UPDATE: at the moment only place where to set open ranges
    // UPDATE: not sure if setting pointer to smtg exiting the domain I still have control
    init_attr(HSP_PARAM_MIN, domain_min);
    init_attr(HSP_PARAM_MAX, domain_max);
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
    if (arg_name == HSP_PARAM_X)
        return x_arg_;
    else
        for (Attribute &attr : arguments_)
            if (attr.name() == arg_name)
                return attr;

    // else throw error showing options
    str_vec options = {HSP_PARAM_X};
    for (const Attribute &attr : arguments_)
        options.push_back(attr.name());

    // throw error showing options
    throw std::invalid_argument(
        "Unknown argument \"" + arg_name + "\" for " + err_name() + "." +
        tools::str::find_closest(arg_name, options, HSP_MIN_STR_DIST));
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
    if (!range_min_.is_constant())
        return false;
    if (!range_max_.is_constant())
        return false;

    // else is constant
    return true;
}

/////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////  Function : set range  //////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

/// Returns the case of the range of validity.
unsigned short Function::which_valid_case_() const
{
    return (!range_min_.is_empty()) + 2 * (!range_max_.is_empty());
}

/// Ensure uniform range (min, max) match argument (min, max).
void Function::init_uniform_range_()
{
    // impose range minimum equal to argument minimum
    const Attribute &min_arg = find_arg_(HSP_PARAM_MIN);
    if (min_arg.is_function())
        range_min_.init_function(min_arg.get_function());
    else
        range_min_.init_pointer(min_arg.get_ptr());

    // impose range maximum equal to argument maximum
    const Attribute &max_arg = find_arg_(HSP_PARAM_MAX);
    if (max_arg.is_function())
        range_max_.init_function(max_arg.get_function());
    else
        range_max_.init_pointer(max_arg.get_ptr());
}

/// Functional range for lower bound and no upper bound.
bool Function::init_range_min_(double min_val, bool is_min_open)
{
    if (get_min() < min_val)
        throw std::invalid_argument(
            "Invalid min value " + str_min() + " for " + err_name() + ".");
    else if ((get_min() > min_val) || (!tools::is_na(get_max())))
    {
        init_range_helper_(((is_min_open) ? (get_min() == min_val) : false), false);
        return true;
    }
    return false;
}

/// Set the range of validity of the value, and the corresponding validation function.
void Function::init_range_helper_(bool is_min_open, bool is_max_open)
{
    switch (which_valid_case_())
    {
    case 0: // Unbounded
        is_below_ = &Function::do_is_below_never_;
        is_above_ = &Function::do_is_above_never_;
        is_valid_ = &Function::do_is_valid_always_;
        break;

    case 1: // Lower Bound
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
        break;

    case 2: // Upper Bound
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
        break;

    case 3: // Within an Interval
        if ((is_min_open) && (is_max_open))
        {
            is_below_ = &Function::do_is_below_open_min_;
            is_above_ = &Function::do_is_above_open_max_;
            is_valid_ = &Function::do_is_valid_both_open_;
        }
        else if ((is_min_open) && (!is_max_open))
        {
            is_below_ = &Function::do_is_below_open_min_;
            is_above_ = &Function::do_is_above_close_max_;
            is_valid_ = &Function::do_is_valid_open_min_close_max_;
        }
        else if ((!is_min_open) && (is_max_open))
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
        bool wrong_open = (is_valid_ != &Function::do_is_valid_both_close_) && (range_min_.get() >= range_max_.get());
        bool wrong_close = (is_valid_ == &Function::do_is_valid_both_close_) && (range_min_.get() > range_max_.get());
        if (wrong_open || wrong_close)
            throw std::invalid_argument(
                "Inconsistent range " + str_range() + " for " + err_name() + ".");
        break;
    }
}

/// Return formatted range.
std::string Function::str_range() const
{
    if (is_valid_ == &Function::do_is_valid_always_)
        return "(-Inf, Inf)";

    else if (is_valid_ == &Function::do_is_valid_open_min_)
        return "(" + str_min() + ", Inf)";

    else if (is_valid_ == &Function::do_is_valid_close_min_)
        return "[" + str_min() + ", Inf)";

    else if (is_valid_ == &Function::do_is_valid_open_max_)
        return "(-Inf, " + str_max() + ")";

    else if (is_valid_ == &Function::do_is_valid_close_max_)
        return "(-Inf, " + str_max() + "]";

    else if (is_valid_ == &Function::do_is_valid_both_open_)
        return "(" + str_min() + ", " + str_max() + ")";

    else if (is_valid_ == &Function::do_is_valid_open_min_close_max_)
        return "(" + str_min() + ", " + str_max() + "]";

    else if (is_valid_ == &Function::do_is_valid_close_min_open_max_)
        return "[" + str_min() + ", " + str_max() + ")";

    else // (is_valid_ == &Function::do_is_valid_both_close_)
        return "[" + str_min() + ", " + str_max() + "]";
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
    if (attr.is_function())
        if (!attr.get_function()->is_implicit_callable())
            throw std::invalid_argument(
                "Impossible to associate " + attr.err_name() + " as argument \"" + key + "\" of " + err_name() + ".\n" +
                "Value of argument \"x\" of " + attr.err_name() + " cannot be missing.");

    // distinguish between range keys and standard keys
    if (HSP_RANGE_KEYS % key)
    {
        if (attr.is_empty())
            ; // skip NA for range
        else if (name_ == "uniform")
            find_arg_(key).transfer(attr);
        else if (key == HSP_PARAM_MIN)
            range_min_.transfer(attr);
        else // key == HSP_PARAM_MAX
            range_max_.transfer(attr);
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
    switch (discr_type)
    {
    case HSP_DISCR_DEFINED: // discrete defined
        function_ = &Function::do_function_discrete_;
        log_func_ = &Function::do_logfunc_discrete_;
        primitive_ = &Function::do_primitive_discrete_;
        log_ratio_ = &Function::do_logratio_discrete_;
        break;

    case HSP_DISCR_CENTER: // centered discretization
        function_ = &Function::do_function_discr_center_;
        log_func_ = &Function::do_logfunc_discr_center_;
        primitive_ = &Function::do_primitive_discr_center_;
        log_ratio_ = &Function::do_logratio_discr_center_;
        break;

    case HSP_DISCR_LEFT: // left discretization
        function_ = &Function::do_function_discr_left_;
        log_func_ = &Function::do_logfunc_discr_left_;
        primitive_ = &Function::do_primitive_discr_left_;
        log_ratio_ = &Function::do_logratio_discr_left_;
        break;

    default: // no discretization (default)
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
        return HSP_DISCR_DEFINED;
    else if (function_ == &Function::do_function_discr_center_)
        return HSP_DISCR_CENTER;
    else if (function_ == &Function::do_function_discr_left_)
        return HSP_DISCR_LEFT;
    return HSP_DISCR_NONE;
}

/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////  Function : show  ///////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

/// Return formatted argument names of the distribution.
std::string Function::str_arg_names() const
{
    if (!arguments_.empty())
    {
        std::string output = arguments_[0].name();
        for (size_t c = 1; c < arguments_.size(); ++c)
            output += HSP_LIST_DELIM + arguments_[c].name();
        return output;
    }
    return HSP_NA_STRING;
}
/// Return formatted argument names of the distribution.
std::string Function::str_arg_values() const
{
    if (!arguments_.empty())
    {
        std::string output = arguments_[0].str();
        for (size_t c = 1; c < arguments_.size(); ++c)
            output += HSP_LIST_DELIM + arguments_[c].str();
        return output;
    }
    return HSP_NA_STRING;
}
/// Return formatted function summary.
std::string Function::str_summary(bool show_x_arg) const
{
    std::string output = name_ + "(";
    if (show_x_arg)
        output += x_arg_.name() + "=" + x_arg_.str();
    for (size_t c = 0; c < arguments_.size(); ++c)
        output += ", " + arguments_[c].name() + "=" + arguments_[c].str();
    output += ") : " + str_range();
    return output;
}

#ifdef _HOUSEPI_RCPP
/// Export argument list to R.
Rcpp::NumericVector Function::rcpp_export_args() const
{
    using namespace Rcpp;

    size_t n_args = arguments_.size();
    NumericVector vals(n_args);
    CharacterVector keys(n_args);
    for (size_t k = 0; k < n_args; ++k)
    {
        vals[k] = arguments_[k].get();
        keys[k] = arguments_[k].name();
    }
    vals.names() = std::move(keys);
    return vals;
}
/// Export range list to R.
Rcpp::NumericVector Function::rcpp_export_range() const
{
    using namespace Rcpp;
    return NumericVector::create(_[HSP_PARAM_MIN] = range_min_.get(),
                                 _[HSP_PARAM_MAX] = range_max_.get());
}
#endif

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

std::shared_ptr<Function> FunctionRegistry::get(const std::string &opt, bool is_primitive) const
{
    auto output = factory_map_.find(opt)->second();
    if (is_primitive)
        output->set_eval_primitive();
    else
        output->set_eval_function();
    return output;
}

std::shared_ptr<Function> compose_func_registries(const std::string &opt,
                                                  const std::string &name,
                                                  const std::vector<const FunctionRegistry *> &registries,
                                                  bool is_primitive)
{
    str_vec options = {};

    if (registries.empty())
        throw std::invalid_argument("No registries provided for \"" + name + "\".");
    else if (tools::is_na(opt))
        return std::make_shared<Function>();
    else
        for (const auto &r : registries)
            if (r->is_valid_option(opt))
                return r->get(opt, is_primitive);
            else
                options += r->get_options();

    throw std::invalid_argument(
        "Unknown " + name + " \"" + opt + "\"." +
        tools::str::find_closest(opt, options, HSP_MIN_STR_DIST));
}

/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////  FunctionRegistry : entries  //////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

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

/// Arithmetic function factory.
#define MAKE_ARITHMETIC_FUNCTION(ClassName, FunctionName, Symbol, DefaultY)             \
    class ClassName : public Function                                                   \
    {                                                                                   \
    public:                                                                             \
        ClassName() : Function(FunctionName, {"x2"}, {DefaultY}) {}                     \
                                                                                        \
        void init_truncation() override final {}                                        \
                                                                                        \
        bool is_discrete() const override final                                         \
        {                                                                               \
            return arguments_[0].is_discrete() && x_arg_.is_discrete();                 \
        }                                                                               \
                                                                                        \
    protected:                                                                          \
        double raw_function_(double x) const override final { return x Symbol arg(0); } \
    };
MAKE_ARITHMETIC_FUNCTION(Addition_Function, "add", +, 0.)
MAKE_ARITHMETIC_FUNCTION(Subtraction_Function, "sub", -, 0.)
MAKE_ARITHMETIC_FUNCTION(Multiplication_Function, "mult", *, 1.)

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
    Power_Function()
        : Function("pow", {"x2"}, {1.}) {}

    void init_truncation() override final { init_range_helper_(); }

    bool is_discrete() const override final
    {
        return Function::is_discrete() || (arguments_[0].is_discrete() && x_arg_.is_discrete());
    }

protected:
    double raw_function_(double x) const override final { return std::pow(x, arg(0)); }
    double raw_logfunc_(double x) const override final { return arg(0) * std::log(x); }
    double raw_primitive_(double x) const override final
    {
        if (std::abs(arg(0) + 1.) < 1.e-10)
            return std::log(std::abs(x));
        else
            return std::pow(x, arg(0) + 1.) / (arg(0) + 1.);
    }
    double raw_logratio_(double xo, double xn) const override final { return arg(0) * std::log(xn / xo); }
};

/// SquareRoot Function.
class SquareRoot_Function : public Function
{
public:
    SquareRoot_Function() : Function("sqrt", {}, {},
                                     0., false, HSP_NA_DOUBLE, false) {}

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
                                    0., true, HSP_NA_DOUBLE, false) {}

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
    Variadic_Function(const std::string &name = HSP_NA_STRING) : Function(name) {}
    virtual ~Variadic_Function() = default;

protected:
    /// Variadic function version of the find argument method.
    Attribute &find_arg_(const std::string &arg_name) override final
    {
        if (arg_name == HSP_PARAM_X)
            return x_arg_;
        else
            for (Attribute &attr : arguments_)
                if (attr.name() == arg_name)
                    return attr;

        // if attribute doesn't exist at initialization, create it
        arguments_.emplace_back();
        return arguments_.back();
    }
};

/// VARIADIC function factory.
#define MAKE_VARIADIC_FUNCTION(ClassName, FunctionName, Symbol) \
    class ClassName : public Variadic_Function                  \
    {                                                           \
    public:                                                     \
        ClassName() : Variadic_Function(FunctionName) {}        \
                                                                \
        void init_truncation() override final {}                \
        bool is_discrete() const override final                 \
        {                                                       \
            if (!x_arg_.is_discrete())                          \
                return false;                                   \
            for (const Attribute &attr : arguments_)            \
                if (!attr.is_discrete())                        \
                    return false;                               \
            return true;                                        \
        }                                                       \
                                                                \
    protected:                                                  \
        double raw_function_(double x) const override final     \
        {                                                       \
            for (const Attribute &attr : arguments_)            \
                x Symbol attr.get();                            \
            return x;                                           \
        }                                                       \
    };

MAKE_VARIADIC_FUNCTION(Summation_Function, "sum", +=)
MAKE_VARIADIC_FUNCTION(Product_Function, "prod", *=)

/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////  Variadic_Function : class  ///////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

#define IMPLEMENT_HOUSEPI_FUNCTION_REGISTRY(ClassName) \
    ClassName &ClassName::inst()                       \
    {                                                  \
        static ClassName _inst;                        \
        return _inst;                                  \
    }

IMPLEMENT_HOUSEPI_FUNCTION_REGISTRY(ArithmeticRegistry)
REGISTER_HOUSEPI_FUNCTION(ArithmeticRegistry, Minus_Function)
REGISTER_HOUSEPI_FUNCTION(ArithmeticRegistry, Addition_Function)
REGISTER_HOUSEPI_FUNCTION(ArithmeticRegistry, Subtraction_Function)
REGISTER_HOUSEPI_FUNCTION(ArithmeticRegistry, Multiplication_Function)
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
