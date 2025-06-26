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
    _do_get = &Attribute::__get_default;
    _do_get_pointer = &Attribute::__get_pointer_default;
    _do_set = &Attribute::__set_default;
    _value = x;
    _ptr_to_value = nullptr;
    _func_ptr = nullptr;
}
void Attribute::assign_pointer(const double *y)
{
    _do_get = &Attribute::__get_from_pointer;
    _do_get_pointer = &Attribute::__get_pointer_from_pointer;
    _do_set = &Attribute::__set_error;
    _value = HSP_NA_DOUBLE;
    _ptr_to_value = y;
    _func_ptr = nullptr;
}
void Attribute::assign_function(const std::shared_ptr<Function> f)
{
    _do_get = &Attribute::__get_from_function;
    _do_get_pointer = &Attribute::__get_pointer_error;
    _do_set = &Attribute::__set_error;
    _value = HSP_NA_DOUBLE;
    _ptr_to_value = nullptr;
    _func_ptr = f;
}

/* Get attribute value from external function. */
double Attribute::__get_from_function() const { return _func_ptr->return_value(); }

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
    bool is_x_default = contains(arg_name, HSP_X_);

    // name of the function
    _function_name = name;

    // number of arguments
    if (arg_name.size() != arg_value.size())
    {
        std::ostringstream oss;
        oss << "Inconsistent arguments for \"" << _function_name << "\"."
            << " Found " << arg_name.size() << " default names but "
            << arg_value.size() << " values.";
        throw std::length_error(oss.str());
    }
    _n_args = arg_name.size() - (is_x_default);

    // assign default arguments
    _arguments = new Attribute[_n_args + 1];
    // by convention, last argument is "x" and NA by default
    _arguments[_n_args].set_name(HSP_X_);
    _arguments[_n_args].set(HSP_NA_DOUBLE);

    ushort c = 0;
    for (ushort k = 0; k < arg_name.size(); ++k)
        if (arg_name[k] == HSP_X_)
            _arguments[_n_args].set(arg_value[k]);
        else
        {
            _arguments[c].set_name(arg_name[k]);
            _arguments[c++].set(arg_value[k]);
        }

    //----- assign domain -----//
    // NOTE: at the moment only place where to set open ranges
    // NOTE: not sure if setting pointer to smtg exiting the domain I still have control
    _assign_range_scalar("min", domain_min);
    _assign_range_scalar("max", domain_max);
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
    for (ushort i = 0; i < _n_args; ++i)
        if (tools::is_na(arg(i)))
            return false;
    return true;
}

/* Assert argument exists at assignment. Returns it's index. */
ushort Function::assert_arg(const std::string &arg_name) const
{
    if (arg_name == HSP_X_)
        return _n_args; // by convention, last argument is "x"
    else
        for (ushort c = 0; c < _n_args; ++c)
            if (_arguments[c].show_name() == arg_name)
                return c;

    str_vec options = {};
    for (ushort c = 0; c < _n_args; ++c)
        options.push_back(_arguments[c].show_name());
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
    bool is_lower_bound = !(_range_min.is_scalar() && tools::is_na(_range_min.get()));
    bool is_upper_bound = !(_range_max.is_scalar() && tools::is_na(_range_max.get()));
    return is_lower_bound + 2 * is_upper_bound;
}

/* Ensure uniform range (min, max) match argument (min, max). */
void Function::init_uniform_range()
{
    // impose range minimum equal to argument minimum
    Attribute &min_arg = _arguments[this->assert_arg("min")];
    if (min_arg.is_function())
        _range_min.assign_function(min_arg.get_function());
    else
        _range_min.assign_pointer(min_arg.get_pointer());

    // impose range maximum equal to argument maximum
    Attribute &max_arg = _arguments[this->assert_arg("max")];
    if (max_arg.is_function())
        _range_max.assign_function(max_arg.get_function());
    else
        _range_max.assign_pointer(max_arg.get_pointer());
}

/* Set the range of validity of the value, and the corresponding validation function. */
void Function::init_functional_range(bool is_min_open, bool is_max_open)
{
    switch (which_valid_case())
    {
    case 0: //----- unbounded -----//
        _do_below = &Function::__below_never;
        _do_above = &Function::__above_never;
        _do_valid = &Function::__valid_always;
        break;

    case 1: //----- lower bound -----//
        _do_above = &Function::__above_never;
        if (is_min_open)
        {
            _do_below = &Function::__below_open_min;
            _do_valid = &Function::__valid_open_min;
        }
        else
        {
            _do_below = &Function::__below_close_min;
            _do_valid = &Function::__valid_close_min;
        }
        break;

    case 2: //----- upper bound -----//
        _do_below = &Function::__below_never;
        if (is_max_open)
        {
            _do_above = &Function::__above_open_max;
            _do_valid = &Function::__valid_open_max;
        }
        else
        {
            _do_above = &Function::__above_close_max;
            _do_valid = &Function::__valid_close_max;
        }
        break;

    case 3: //----- within an interval -----//
        if ((is_min_open) && (is_max_open))
        {
            _do_below = &Function::__below_open_min;
            _do_above = &Function::__above_open_max;
            _do_valid = &Function::__valid_both_open;
        }
        else if ((is_min_open) && (!is_max_open))
        {
            _do_below = &Function::__below_open_min;
            _do_above = &Function::__above_close_max;
            _do_valid = &Function::__valid_open_min_close_max;
        }
        else if ((!is_min_open) && (is_max_open))
        {
            _do_below = &Function::__below_close_min;
            _do_above = &Function::__above_open_max;
            _do_valid = &Function::__valid_close_min_open_max;
        }
        else
        {
            _do_below = &Function::__below_close_min;
            _do_above = &Function::__above_close_max;
            _do_valid = &Function::__valid_both_close;
        }

        // check for consistency
        bool wrong_open = (_do_valid != &Function::__valid_both_close) && (_range_min.get() >= _range_max.get());
        bool wrong_close = (_do_valid == &Function::__valid_both_close) && (_range_min.get() > _range_max.get());
        if (wrong_open || wrong_close)
            throw std::invalid_argument(
                "Inconsistent range " + show_range() + ".");
        break;
    }
}

/* Return formatted range. */
const std::string Function::show_range() const
{
    if (_do_valid == &Function::__valid_always)
        return "";

    else if (_do_valid == &Function::__valid_open_min)
        return "(" + show_min() + ", Inf)";

    else if (_do_valid == &Function::__valid_close_min)
        return "[" + show_min() + ", Inf)";

    else if (_do_valid == &Function::__valid_open_max)
        return "(-Inf, " + show_max() + ")";

    else if (_do_valid == &Function::__valid_close_max)
        return "(-Inf, " + show_max() + "]";

    else if (_do_valid == &Function::__valid_both_open)
        return "(" + show_min() + ", " + show_max() + ")";

    else if (_do_valid == &Function::__valid_open_min_close_max)
        return "(" + show_min() + ", " + show_max() + "]";

    else if (_do_valid == &Function::__valid_close_min_open_max)
        return "[" + show_min() + ", " + show_max() + ")";

    else // (_do_valid == &Function::__valid_both_close)
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
        _assign_arg_scalar(arg_name, x);
    else
        _assign_range_scalar(arg_name, x);
}
//* Assign pointer argument or range value. *//
void Function::assign_attr(const std::string &arg_name, bool is_arg, const double *y)
{
    // check the function is different from the current one
    if (y == nullptr)
        throw std::invalid_argument(
            "Invalid pointer argument \"" + arg_name + "\" for \"" + _function_name + "\".");

    if (is_arg == HSP_ATTR_ARGUMENT)
        _assign_arg_pointer(arg_name, y);
    else
        _assign_range_pointer(arg_name, y);
}
//* Assign fucntion argument or range value. *//
void Function::assign_attr(const std::string &arg_name, bool is_arg, std::shared_ptr<Function> f)
{
    // check the function is different from the current one
    if (f == nullptr)
        throw std::invalid_argument(
            "Invalid function argument \"" + arg_name + "\" for \"" + _function_name + "\".");
    else if (f.get() == this)
        throw std::invalid_argument(
            "Cannot assign function \"" + _function_name + "\" to its own argument \"" + arg_name + "\".");

    if (is_arg == HSP_ATTR_ARGUMENT)
        _assign_arg_function(arg_name, f);
    else
        _assign_range_function(arg_name, f);
}

/* Assignment scalar argument. */
void Function::_assign_arg_scalar(const std::string &arg_name, double x)
{
    ushort i = this->assert_arg(arg_name);
    _arguments[i].assign_scalar(x);
}
/* Assignment pointer argument. */
void Function::_assign_arg_pointer(const std::string &arg_name, const double *y)
{
    ushort i = this->assert_arg(arg_name);
    _arguments[i].assign_pointer(y);
}
/* Assignment function argument. */
void Function::_assign_arg_function(const std::string &arg_name, std::shared_ptr<Function> f)
{
    ushort i = this->assert_arg(arg_name);
    _arguments[i].assign_function(f);
}
/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////  Base Function : Range  ///////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

/* Assign scalar range. */
void Function::_assign_range_scalar(const std::string &range_name, double x)
{
    if (_function_name == "unif")
    {
        if ((range_name != "min") && (range_name != "max"))
            throw std::invalid_argument(
                "Unknown range name \"" + range_name + "\" for \"" + _function_name + "\".");
        else if (!tools::is_na(x))
            _assign_arg_scalar(range_name, x);
    }
    else
    {
        if (range_name == "min")
            _range_min.assign_scalar(x);
        else if (range_name == "max")
            _range_max.assign_scalar(x);
        else
            throw std::invalid_argument(
                "Unknown range name \"" + range_name + "\" for \"" + _function_name + "\".");
    }
}

/* Assign pointer range. */
void Function::_assign_range_pointer(const std::string &range_name, const double *y)
{
    if (_function_name == "unif")
    {
        if ((range_name != "min") && (range_name != "max"))
            throw std::invalid_argument(
                "Unknown range name \"" + range_name + "\" for \"" + _function_name + "\".");
        else if (y != nullptr)
            _assign_arg_pointer(range_name, y);
    }
    else
    {
        if (range_name == "min")
            _range_min.assign_pointer(y);
        else if (range_name == "max")
            _range_max.assign_pointer(y);
        else
            throw std::invalid_argument(
                "Unknown range name \"" + range_name + "\" for \"" + _function_name + "\".");
    }
}

/* Assign function range. */
void Function::_assign_range_function(const std::string &range_name, const std::shared_ptr<Function> f)
{
    if (_function_name == "unif")
    {
        if ((range_name != "min") && (range_name != "max"))
            throw std::invalid_argument(
                "Unknown range name \"" + range_name + "\" for \"" + _function_name + "\".");
        else if (f != nullptr)
            _assign_arg_function(range_name, f);
    }
    else
    {
        if (range_name == "min")
            _range_min.assign_function(f);
        else if (range_name == "max")
            _range_max.assign_function(f);
        else
            throw std::invalid_argument(
                "Unknown range name \"" + range_name + "\" for \"" + _function_name + "\".");
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
    if (_n_args > 0)
    {
        std::string output = _arguments[0].show_name();
        for (ushort c = 1; c < _n_args; ++c)
            output += C_LIST_DELIM + _arguments[c].show_name();
        return output;
    }

    return HSP_NA_STRING;
}

/* Return formatted argument names of the distribution. */
const std::string Function::show_arg_values() const
{
    if (_n_args > 0)
    {
        std::string output = _arguments[0].show();
        for (ushort c = 1; c < _n_args; ++c)
            output += C_LIST_DELIM + _arguments[c].show();
        return output;
    }

    return HSP_NA_STRING;
}

/* Return formatted function summary. */
const std::string Function::show_summary() const
{
    std::string output = show_name();
    if (_n_args > 0)
    {
        output += "(" + _arguments[0].show_name() + "=" + _arguments[0].show();
        for (ushort c = 1; c < _n_args; ++c)
            output += ", " + _arguments[c].show_name() + "=" + _arguments[c].show();
        output += ")";
    }
    if (_do_valid != &Function::__valid_always)
        output += " : " + show_range();

    return output;
}

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
