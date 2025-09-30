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
//////////////////////////////////////  Attribute : getter //////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

/// Get attribute value from external function.
/// note: this can't be moved in the class because of forward declaration of Function
double Attribute::do_get_from_function_() const { return pFunc_->evaluate(); }

/// Wheter it's constant.
/// NOTE: by construction convention, a pointer to constant must initialize to scalar instead.
/// NOTE: if empty, then it's constant by convention.
bool Attribute::is_constant() const
{
    if (is_scalar())
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
    if (is_scalar())
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

/// Initialize attribute as scalar.
void Attribute::init_scalar(double x)
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
        init_scalar(other.get());
    else if (other.is_pointer())
        init_pointer(other.get_ptr());
    else if (other.is_function())
        init_function(other.get_function());
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
    bool is_x_default = arg_name % "x";

    // name of the function
    this->name_ = name;

    // number of arguments
    if (arg_name.size() != arg_value.size())
    {
        std::ostringstream oss;
        oss << "Inconsistent arguments for " << err_name() << "."
            << "Found " << arg_name.size() << " default names but "
            << arg_value.size() << " values.";
        tools::error(oss.str());
    }
    n_args_ = arg_name.size() - (is_x_default);

    // assign default arguments
    args_ = new Attribute[n_args_ + 1];
    // by convention, last argument is 'x' and NA by default
    args_[n_args_].name("x");
    args_[n_args_].set(HSP_NA_DOUBLE);

    unsigned short c = 0;
    for (unsigned short k = 0; k < arg_name.size(); ++k)
        if (arg_name[k] == "x")
            args_[n_args_].set(arg_value[k]);
        else
        {
            args_[c].name(arg_name[k]);
            args_[c++].set(arg_value[k]);
        }

    // Assign Domain
    // NOTE: at the moment only place where to set open ranges
    // NOTE: not sure if setting pointer to smtg exiting the domain I still have control
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
unsigned short Function::assert_arg(const std::string &arg_name) const
{
    if (arg_name == "x")
        return n_args_; // by convention, last argument is 'x'
    else
        for (unsigned short i = 0; i < n_args_; ++i)
            if (args_[i].name() == arg_name)
                return i;

    // throw error showing options
    str_vec options = {};
    for (unsigned short i = 0; i < n_args_; ++i)
        options.push_back(args_[i].name());
    throw std::invalid_argument(
        "Unknown argument \"" + arg_name + "\" for " + err_name() + "." +
        tools::str::find_closest(arg_name, options, HSP_MIN_STR_DIST));
}

/// Wheter return value is always constant.
bool Function::is_constant() const
{
    // check arguments
    for (unsigned short i = 0; i < n_args_; ++i)
        if (!args_[i].is_constant())
            return false;
    // check X value
    // NOTE: a function with undefined 'x' is not constant
    if (args_[n_args_].is_empty() || !args_[n_args_].is_constant())
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
    const Attribute &min_arg = args_[this->assert_arg(HSP_PARAM_MIN)];
    if (min_arg.is_function())
        range_min_.init_function(min_arg.get_function());
    else
        range_min_.init_pointer(min_arg.get_ptr());

    // impose range maximum equal to argument maximum
    const Attribute &max_arg = args_[this->assert_arg(HSP_PARAM_MAX)];
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
void Function::init_attr(const std::string &what, const Attribute &attr)
{
    // assert that function attribute are implicitely callable
    if (attr.is_function())
        if (!attr.get_function()->is_implicit_callable())
            throw std::invalid_argument(
                "Undefined attribute \"" + what + "\" for " + err_name() + "    .");

    // distinguish between range keys and standard keys
    if (HSP_RANGE_KEYS % what)
    {
        if (attr.is_empty())
            ; // skip NA for range
        else if (name_ == "unif")
            args_[assert_arg(what)].transfer(attr);
        else if (what == HSP_PARAM_MIN)
            range_min_.transfer(attr);
        else // what == HSP_PARAM_MAX
            range_max_.transfer(attr);
    }
    else
        args_[assert_arg(what)].transfer(attr);
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
    if (n_args_ > 0)
    {
        std::string output = args_[0].name();
        for (unsigned short c = 1; c < n_args_; ++c)
            output += HSP_LIST_DELIM + args_[c].name();
        return output;
    }
    return HSP_NA_STRING;
}
/// Return formatted argument names of the distribution.
std::string Function::str_arg_values() const
{
    if (n_args_ > 0)
    {
        std::string output = args_[0].str();
        for (unsigned short c = 1; c < n_args_; ++c)
            output += HSP_LIST_DELIM + args_[c].str();
        return output;
    }
    return HSP_NA_STRING;
}

/// Return formatted function summary.
std::string Function::str_summary() const
{
    std::string output = name_;
    output += "(" + args_[n_args_].name() + "=" + args_[n_args_].str();
    for (unsigned short c = 0; c < n_args_; ++c)
        output += ", " + args_[c].name() + "=" + args_[c].str();
    output += ") : " + str_range();
    return output;
}

#ifdef _HOUSEPI_RCPP
/// Export argument list to R.
Rcpp::NumericVector Function::rcpp_export_args() const
{
    using namespace Rcpp;

    NumericVector vals(n_args_);
    CharacterVector keys(n_args_);
    for (size_t k = 0; k < n_args_; ++k)
    {
        vals[k] = args_[k].get();
        keys[k] = args_[k].name();
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
                                                  bool is_primitive,
                                                  const std::vector<const FunctionRegistry *> &registries)
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

class Minus_Function : public Function
{
public:
    Minus_Function() : Function("minus", {}, {}) {}

    void init_truncation() override final {}
    bool is_discrete() const override final { return args_[n_args_].is_discrete(); }

protected:
    double raw_function_(double x) const override final { return -x; }
};

#define DEFINE_ARITHMETIC_FUNCTION(ClassName, FunctionName, Symbol, DefaultY)           \
    class ClassName : public Function                                                   \
    {                                                                                   \
    public:                                                                             \
        ClassName() : Function(FunctionName, {"x2"}, {DefaultY}) {}                     \
                                                                                        \
        void init_truncation() override final {}                                        \
                                                                                        \
        bool is_discrete() const override final                                         \
        {                                                                               \
            return args_[0].is_discrete() && args_[n_args_].is_discrete();              \
        }                                                                               \
                                                                                        \
    protected:                                                                          \
        double raw_function_(double x) const override final { return x Symbol arg(0); } \
    };
DEFINE_ARITHMETIC_FUNCTION(Addition_Function, "add", +, 0.)
DEFINE_ARITHMETIC_FUNCTION(Subtraction_Function, "sub", -, 0.)
DEFINE_ARITHMETIC_FUNCTION(Multiplication_Function, "mult", *, 1.)

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
        return Function::is_discrete() || (args_[0].is_discrete() && args_[n_args_].is_discrete());
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

IMPLEMENT_HOUSEPI_FUNCTION_REGISTRY(TranscendentalRegistry)
REGISTER_HOUSEPI_FUNCTION(TranscendentalRegistry, Power_Function)
REGISTER_HOUSEPI_FUNCTION(TranscendentalRegistry, SquareRoot_Function)
REGISTER_HOUSEPI_FUNCTION(TranscendentalRegistry, Exponential_Function)
REGISTER_HOUSEPI_FUNCTION(TranscendentalRegistry, Logarithm_Function)
REGISTER_HOUSEPI_FUNCTION(TranscendentalRegistry, Sine_Function)
REGISTER_HOUSEPI_FUNCTION(TranscendentalRegistry, Cosine_Function)
REGISTER_HOUSEPI_FUNCTION(TranscendentalRegistry, Tangent_Function)

IMPLEMENT_HOUSEPI_FUNCTION_REGISTRY(UserFuncRegistry)
