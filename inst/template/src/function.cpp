/*      __
 *     / /_  ____  __  __ ____ __  ____  __
 *    / __ \/ __ \/ / / / ___/ _ \/ __ \/ /
 *   / / / / /_/ / /_/ (__  )  __/ /_/ / /
 *  /_/ /_/\____/\____/____/\___/ ____/_/
 *                             /_/
 *
 * Copyright (C) 2025 Francesco Camaglia, Institut Pasteur
 *
 * function.cpp
 *
 */

#include "function.h"

//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//

/////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////  Attribute : setter //////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

void Attribute::assign_scalar(double x)
{
    get_ = &Attribute::do_get_default_;
    get_ptr_ = &Attribute::do_get_ptr_default_;
    _do_set = &Attribute::__set_default;
    value_ = x;
    _ptr_to_value = nullptr;
    pFunc_ = nullptr;
}
void Attribute::assign_pointer(const double *y)
{
    get_ = &Attribute::__get_from_pointer;
    get_ptr_ = &Attribute::__get_pointer_from_pointer;
    _do_set = &Attribute::__set_error;
    value_ = HSP_NA_DOUBLE;
    _ptr_to_value = y;
    pFunc_ = nullptr;
}
void Attribute::assign_function(const std::shared_ptr<Function> f)
{
    get_ = &Attribute::__get_from_function;
    get_ptr_ = &Attribute::__get_pointer_error;
    _do_set = &Attribute::__set_error;
    value_ = HSP_NA_DOUBLE;
    _ptr_to_value = nullptr;
    pFunc_ = f;
}

/* Get attribute value from external function. */
double Attribute::__get_from_function() const { return pFunc_->return_value(); }

//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//

/////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////  Function : constructor  /////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

const std::string HSP_X_ = "x";

/* Initial defaults for function. */
Function::Function(const std::string &name,
                   const str_vec &arg_name,
                   const dbl_vec &arg_value,
                   double domain_min,
                   double domain_max,
                   bool is_min_open,
                   bool is_max_open)
{
    bool is_x_default = arg_name % HSP_X_;

    // name of the function
    function_name_ = name;

    // number of arguments
    if (arg_name.size() != arg_value.size())
    {
        std::ostringstream oss;
        oss << "Inconsistent arguments for \"" << function_name_ << "\"."
            << " Found " << arg_name.size() << " default names but "
            << arg_value.size() << " values.";
        throw std::length_error(oss.str());
    }
    n_args_ = arg_name.size() - (is_x_default);

    // assign default arguments
    arguments_ = new Attribute[n_args_ + 1];
    // by convention, last argument is "x" and NA by default
    arguments_[n_args_].set_name(HSP_X_);
    arguments_[n_args_].set(HSP_NA_DOUBLE);

    ushort c = 0;
    for (ushort k = 0; k < arg_name.size(); ++k)
        if (arg_name[k] == HSP_X_)
            arguments_[n_args_].set(arg_value[k]);
        else
        {
            arguments_[c].set_name(arg_name[k]);
            arguments_[c++].set(arg_value[k]);
        }

    //----- assign domain -----//
    // NOTE: at the moment only place where to set open ranges
    // NOTE: not sure if setting pointer to smtg exiting the domain I still have control
    assign_range_scalar_("min", domain_min);
    assign_range_scalar_("max", domain_max);
    init_functional_range(is_min_open, is_max_open);
}

Function::~Function() = default;

//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//

/////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////  Base Function : Getter  ///////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

/* Tell if all arguments are defined, except for 'x' argument. */
bool Function::is_able_to_showcase() const
{
    for (ushort i = 0; i < n_args_; ++i)
        if (tools::is_na(arg(i)))
            return false;
    return true;
}

/* Assert argument exists at assignment. Returns it's index. */
ushort Function::assert_arg(const std::string &arg_name) const
{
    if (arg_name == HSP_X_)
        return n_args_; // by convention, last argument is "x"
    else
        for (ushort c = 0; c < n_args_; ++c)
            if (arguments_[c].show_name() == arg_name)
                return c;

    str_vec options = {};
    for (ushort c = 0; c < n_args_; ++c)
        options.push_back(arguments_[c].show_name());
    throw std::invalid_argument(
        "Unknown argument \"" + arg_name + "\" for \"" + show_name() + "\"." + tools::find_closest(arg_name, options));
}

/////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////  Base Function : set range  ////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

/* Returns the case of the range of validity. */
ushort Function::which_valid_case() const
{
    // note: bounded if not pointing internally to NA
    bool is_lower_bound = !(range_min_.is_scalar() && tools::is_na(range_min_.get()));
    bool is_upper_bound = !(range_max_.is_scalar() && tools::is_na(range_max_.get()));
    return is_lower_bound + 2 * is_upper_bound;
}

/* Ensure uniform range (min, max) match argument (min, max). */
void Function::init_uniform_range()
{
    // impose range minimum equal to argument minimum
    Attribute &min_arg = arguments_[this->assert_arg("min")];
    if (min_arg.is_function())
        range_min_.assign_function(min_arg.get_function());
    else
        range_min_.assign_pointer(min_arg.get_ptr());

    // impose range maximum equal to argument maximum
    Attribute &max_arg = arguments_[this->assert_arg("max")];
    if (max_arg.is_function())
        range_max_.assign_function(max_arg.get_function());
    else
        range_max_.assign_pointer(max_arg.get_ptr());
}

/* Set the range of validity of the value, and the corresponding validation function. */
void Function::init_functional_range(bool is_min_open, bool is_max_open)
{
    switch (which_valid_case())
    {
    case 0: //----- unbounded -----//
        below_ = &Function::do_below_never_;
        above_ = &Function::do_above_never_;
        valid_ = &Function::do_valid_always_;
        break;

    case 1: //----- lower bound -----//
        above_ = &Function::do_above_never_;
        if (is_min_open)
        {
            below_ = &Function::do_below_open_min_;
            valid_ = &Function::do_valid_open_min_;
        }
        else
        {
            below_ = &Function::do_below_close_min_;
            valid_ = &Function::do_valid_close_min_;
        }
        break;

    case 2: //----- upper bound -----//
        below_ = &Function::do_below_never_;
        if (is_max_open)
        {
            above_ = &Function::do_above_open_max_;
            valid_ = &Function::do_valid_open_max_;
        }
        else
        {
            above_ = &Function::do_above_close_max_;
            valid_ = &Function::do_valid_close_max_;
        }
        break;

    case 3: //----- within an interval -----//
        if ((is_min_open) && (is_max_open))
        {
            below_ = &Function::do_below_open_min_;
            above_ = &Function::do_above_open_max_;
            valid_ = &Function::do_valid_both_open_;
        }
        else if ((is_min_open) && (!is_max_open))
        {
            below_ = &Function::do_below_open_min_;
            above_ = &Function::do_above_close_max_;
            valid_ = &Function::do_valid_open_min_close_max_;
        }
        else if ((!is_min_open) && (is_max_open))
        {
            below_ = &Function::do_below_close_min_;
            above_ = &Function::do_above_open_max_;
            valid_ = &Function::do_valid_close_min_open_max_;
        }
        else
        {
            below_ = &Function::do_below_close_min_;
            above_ = &Function::do_above_close_max_;
            valid_ = &Function::do_valid_both_close_;
        }

        // check for consistency
        bool wrong_open = (valid_ != &Function::do_valid_both_close_) && (range_min_.get() >= range_max_.get());
        bool wrong_close = (valid_ == &Function::do_valid_both_close_) && (range_min_.get() > range_max_.get());
        if (wrong_open || wrong_close)
            throw std::invalid_argument(
                "Inconsistent range " + show_range() + ".");
        break;
    }
}

/* Return formatted range. */
const std::string Function::show_range() const
{
    if (valid_ == &Function::do_valid_always_)
        return "";

    else if (valid_ == &Function::do_valid_open_min_)
        return "(" + show_min() + ", Inf)";

    else if (valid_ == &Function::do_valid_close_min_)
        return "[" + show_min() + ", Inf)";

    else if (valid_ == &Function::do_valid_open_max_)
        return "(-Inf, " + show_max() + ")";

    else if (valid_ == &Function::do_valid_close_max_)
        return "(-Inf, " + show_max() + "]";

    else if (valid_ == &Function::do_valid_both_open_)
        return "(" + show_min() + ", " + show_max() + ")";

    else if (valid_ == &Function::do_valid_open_min_close_max_)
        return "(" + show_min() + ", " + show_max() + "]";

    else if (valid_ == &Function::do_valid_close_min_open_max_)
        return "[" + show_min() + ", " + show_max() + ")";

    else // (valid_ == &Function::do_valid_both_close_)
        return "[" + show_min() + ", " + show_max() + "]";
}

//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//

/////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////  Base Function : Arguments  /////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

//--- Assign scalar argument or range value. ---//
void Function::assign_attr(const std::string &arg_name, bool is_arg, double x)
{
    if (is_arg == HSP_ATTR_ARGUMENT)
        assign_arg_scalar_(arg_name, x);
    else
        assign_range_scalar_(arg_name, x);
}
//* Assign pointer argument or range value. *//
void Function::assign_attr(const std::string &arg_name, bool is_arg, const double *y)
{
    // check the function is different from the current one
    if (y == nullptr)
        throw std::invalid_argument(
            "Invalid pointer argument \"" + arg_name + "\" for \"" + function_name_ + "\".");

    if (is_arg == HSP_ATTR_ARGUMENT)
        assign_arg_pointer_(arg_name, y);
    else
        assign_range_pointer_(arg_name, y);
}
//* Assign fucntion argument or range value. *//
void Function::assign_attr(const std::string &arg_name, bool is_arg, std::shared_ptr<Function> f)
{
    // check the function is different from the current one
    if (f == nullptr)
        throw std::invalid_argument(
            "Invalid function argument \"" + arg_name + "\" for \"" + function_name_ + "\".");
    else if (f.get() == this)
        throw std::invalid_argument(
            "Cannot assign function \"" + function_name_ + "\" to its own argument \"" + arg_name + "\".");

    if (is_arg == HSP_ATTR_ARGUMENT)
        assign_arg_function_(arg_name, f);
    else
        assign_range_function_(arg_name, f);
}

/* Assignment scalar argument. */
void Function::assign_arg_scalar_(const std::string &arg_name, double x)
{
    ushort i = this->assert_arg(arg_name);
    arguments_[i].assign_scalar(x);
}
/* Assignment pointer argument. */
void Function::assign_arg_pointer_(const std::string &arg_name, const double *y)
{
    ushort i = this->assert_arg(arg_name);
    arguments_[i].assign_pointer(y);
}
/* Assignment function argument. */
void Function::assign_arg_function_(const std::string &arg_name, std::shared_ptr<Function> f)
{
    ushort i = this->assert_arg(arg_name);
    arguments_[i].assign_function(f);
}
/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////  Base Function : Range  ///////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

/* Assign scalar range. */
void Function::assign_range_scalar_(const std::string &range_name, double x)
{
    if (function_name_ == "unif")
    {
        if ((range_name != "min") && (range_name != "max"))
            throw std::invalid_argument(
                "Unknown range name \"" + range_name + "\" for \"" + function_name_ + "\".");
        else if (!tools::is_na(x))
            assign_arg_scalar_(range_name, x);
    }
    else
    {
        if (range_name == "min")
            range_min_.assign_scalar(x);
        else if (range_name == "max")
            range_max_.assign_scalar(x);
        else
            throw std::invalid_argument(
                "Unknown range name \"" + range_name + "\" for \"" + function_name_ + "\".");
    }
}

/* Assign pointer range. */
void Function::assign_range_pointer_(const std::string &range_name, const double *y)
{
    if (function_name_ == "unif")
    {
        if ((range_name != "min") && (range_name != "max"))
            throw std::invalid_argument(
                "Unknown range name \"" + range_name + "\" for \"" + function_name_ + "\".");
        else if (y != nullptr)
            assign_arg_pointer_(range_name, y);
    }
    else
    {
        if (range_name == "min")
            range_min_.assign_pointer(y);
        else if (range_name == "max")
            range_max_.assign_pointer(y);
        else
            throw std::invalid_argument(
                "Unknown range name \"" + range_name + "\" for \"" + function_name_ + "\".");
    }
}

/* Assign function range. */
void Function::assign_range_function_(const std::string &range_name, const std::shared_ptr<Function> f)
{
    if (function_name_ == "unif")
    {
        if ((range_name != "min") && (range_name != "max"))
            throw std::invalid_argument(
                "Unknown range name \"" + range_name + "\" for \"" + function_name_ + "\".");
        else if (f != nullptr)
            assign_arg_function_(range_name, f);
    }
    else
    {
        if (range_name == "min")
            range_min_.assign_function(f);
        else if (range_name == "max")
            range_max_.assign_function(f);
        else
            throw std::invalid_argument(
                "Unknown range name \"" + range_name + "\" for \"" + function_name_ + "\".");
    }
}

//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//

/////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////  Base Function : Show  /////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

/* Return formatted argument names of the distribution. */
const std::string Function::show_arg_names() const
{
    if (n_args_ > 0)
    {
        std::string output = arguments_[0].show_name();
        for (ushort c = 1; c < n_args_; ++c)
            output += C_LIST_DELIM + arguments_[c].show_name();
        return output;
    }

    return HSP_NA_STRING;
}

/* Return formatted argument names of the distribution. */
const std::string Function::show_arg_values() const
{
    if (n_args_ > 0)
    {
        std::string output = arguments_[0].show();
        for (ushort c = 1; c < n_args_; ++c)
            output += C_LIST_DELIM + arguments_[c].show();
        return output;
    }

    return HSP_NA_STRING;
}

/* Return formatted function summary. */
const std::string Function::show_summary() const
{
    std::string output = show_name();
    if (n_args_ > 0)
    {
        output += "(" + arguments_[0].show_name() + "=" + arguments_[0].show();
        for (ushort c = 1; c < n_args_; ++c)
            output += ", " + arguments_[c].show_name() + "=" + arguments_[c].show();
        output += ")";
    }
    if (valid_ != &Function::do_valid_always_)
        output += " : " + show_range();

    return output;
}

#ifdef _HOUSEPI_RCPP
/* Export argument list to R.*/
Rcpp::List Function::export_arg_list() const
{
    using namespace Rcpp;

    List output = List::create();

    for (size_t k = 0; k < n_args_; ++k)
        output.push_back(arguments_[k].get(), arguments_[k].show_name());

    return (output);
}

/* Export range list to R.*/
Rcpp::List Function::export_range_list() const
{
    using namespace Rcpp;

    List output = List::create(_["min"] = range_min_.get(),
                               _["max"] = range_max_.get());

    return (output);
}
#endif

//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//

const str_vec FunctionRegistry::get_options() const
{
    str_vec keys;
    keys.reserve(_factory_map.size());
    for (const auto &pair : _factory_map)
        keys.push_back(pair.first);
    return keys;
}

std::shared_ptr<Function> FunctionRegistry::get(const std::string &opt) const
{
    auto it = _factory_map.find(opt);
    if (it != _factory_map.end())
        return it->second();
    return nullptr;
}

std::shared_ptr<Function> compose_func_registries(const std::string &opt, const std::string &name,
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
                return r->get(opt);
            else
                options += r->get_options();

    throw std::invalid_argument(
        "Unknown " + name + " \"" + opt + "\"." + tools::find_closest(opt, options));
}

//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//

#define IMPLEMENT_HOUSEPI_FUNCTION_REGISTRY(ClassName) \
    ClassName &ClassName::inst()                       \
    {                                                  \
        static ClassName _inst;                        \
        return _inst;                                  \
    }

IMPLEMENT_HOUSEPI_FUNCTION_REGISTRY(ArithmeticRegistry)

/* Equality Function. */
class Equality_Function : public Function
{
public:
    //===== CONSTRUCTOR =====//
    Equality_Function()
        : Function("eq") {}
    ~Equality_Function() = default;

    //===== GETTER =====//
    double return_value(double x) const override final { return x; }
};

// transcendental_functions.cpp
#define DEFINE_ARITHMETIC_FUNCTION(ClassName, FunctionName, Symbol, DefaultY)          \
    class ClassName : public Function                                                  \
    {                                                                                  \
    public:                                                                            \
        ClassName() : Function(FunctionName, {"y"}, {DefaultY}) {}                     \
        ~ClassName() = default;                                                        \
        double return_value(double x) const override final { return x Symbol arg(0); } \
    };
DEFINE_ARITHMETIC_FUNCTION(Addition_Function, "add", +, 0.0)
DEFINE_ARITHMETIC_FUNCTION(Subtraction_Function, "sub", -, 0.0)
DEFINE_ARITHMETIC_FUNCTION(Multiplication_Function, "mult", *, 1.0)
DEFINE_ARITHMETIC_FUNCTION(Division_Function, "div", /, 1.0)

REGISTER_HOUSEPI_FUNCTION(ArithmeticRegistry, Equality_Function)
REGISTER_HOUSEPI_FUNCTION(ArithmeticRegistry, Addition_Function)
REGISTER_HOUSEPI_FUNCTION(ArithmeticRegistry, Subtraction_Function)
REGISTER_HOUSEPI_FUNCTION(ArithmeticRegistry, Multiplication_Function)
REGISTER_HOUSEPI_FUNCTION(ArithmeticRegistry, Division_Function)

IMPLEMENT_HOUSEPI_FUNCTION_REGISTRY(TranscendentalRegistry)

//* Power Function. *//
class Power_Function : public Function
{
public:
    //===== CONSTRUCTOR =====//
    Power_Function()
        : Function("pow", {"y"}, {1.}) {}
    ~Power_Function() = default;

    //===== GETTER =====//
    double return_value(double x) const override final { return std::pow(x, arg(0)); }
};

#define DEFINE_TRANSCENDENTAL_FUNCTION(ClassName, FunctionName)                             \
    class ClassName : public Function                                                       \
    {                                                                                       \
    public:                                                                                 \
        ClassName() : Function(#FunctionName) {}                                            \
        ~ClassName() = default;                                                             \
        double return_value(double x) const override final { return std::FunctionName(x); } \
    };
DEFINE_TRANSCENDENTAL_FUNCTION(SquareRoot_Function, sqrt)
DEFINE_TRANSCENDENTAL_FUNCTION(Exponential_Function, exp)
DEFINE_TRANSCENDENTAL_FUNCTION(Logarithm_Function, log)
DEFINE_TRANSCENDENTAL_FUNCTION(Sine_Function, sin)
DEFINE_TRANSCENDENTAL_FUNCTION(Cosine_Function, cos)
DEFINE_TRANSCENDENTAL_FUNCTION(Tangent_Function, tan)

REGISTER_HOUSEPI_FUNCTION(TranscendentalRegistry, Power_Function)
REGISTER_HOUSEPI_FUNCTION(TranscendentalRegistry, SquareRoot_Function)
REGISTER_HOUSEPI_FUNCTION(TranscendentalRegistry, Exponential_Function)
REGISTER_HOUSEPI_FUNCTION(TranscendentalRegistry, Logarithm_Function)
REGISTER_HOUSEPI_FUNCTION(TranscendentalRegistry, Sine_Function)
REGISTER_HOUSEPI_FUNCTION(TranscendentalRegistry, Cosine_Function)
REGISTER_HOUSEPI_FUNCTION(TranscendentalRegistry, Tangent_Function)

IMPLEMENT_HOUSEPI_FUNCTION_REGISTRY(UserFuncRegistry)
