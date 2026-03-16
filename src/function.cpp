/* housepi - Epidemiological Analysis of Household Outbreaks.
 *
 * Copyright © 2026 Institut Pasteur (Paris), CNRS and INSERM.
 * Distributed under the terms of the GNU General Public License (GPLv3).
 * See the LICENSE file for details.
 */

#include <algorithm>
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
#include "registry.h"

/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////// *** FunctionRegistry : class ***  /////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

HOUSEPI_DEFINE_REGISTRY_INSTANCE_(Custom)
HOUSEPI_DEFINE_REGISTRY_INSTANCE_(Arithmetic)
HOUSEPI_DEFINE_REGISTRY_INSTANCE_(Transcendental)
HOUSEPI_DEFINE_REGISTRY_INSTANCE_(ContinuosDistrib)
HOUSEPI_DEFINE_REGISTRY_INSTANCE_(DiscreteDistrib)
HOUSEPI_DEFINE_REGISTRY_INSTANCE_(InfecProfile)

namespace hsp
{
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

        if (na::is(opt))
            return nullptr;

        std::vector<const FunctionRegistry *> registry_ptrs = {};
        for (const std::string &r : registries)
            if (r == "custom")
                registry_ptrs.push_back(&CustomRegistry::instance());

            else if (r == "transcendental")
                registry_ptrs.push_back(&TranscendentalRegistry::instance());

            else if (r == "arithmetic")
                registry_ptrs.push_back(&ArithmeticRegistry::instance());

            else if (r == "continuos")
                registry_ptrs.push_back(&ContinuosDistribRegistry::instance());

            else if (r == "discrete")
                registry_ptrs.push_back(&DiscreteDistribRegistry::instance());

            else if (r == "infecprofile")
                registry_ptrs.push_back(&InfecProfileRegistry::instance());

            else
                throw std::invalid_argument(
                    "Unknown " + err_type + " registry '" + r + "'.");

        std::vector<std::string> options = {};
        for (const FunctionRegistry *ptr : registry_ptrs)
        {
            if (ptr->count(opt))
                return ptr->create(opt, is_prim);
            const std::vector<std::string> &y = ptr->get_options();
            options.reserve(options.size() + y.size());
            options.insert(options.end(), y.begin(), y.end());
        }

        throw std::invalid_argument(
            "Unknown " + err_type + " '" + std::string(opt) + "'." +
            str::find_closest(opt, options));
    }

    /////////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////  *** Function : class ***  ////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////

    /// Initial defaults for function.
    /// NOTE: to access arguments by index, I must provide separate name/value.
    /// UPDATE: at the moment only place where to set open ranges
    /// UPDATE: not sure if setting pointer to smtg exiting the domain I still have control
    Function::Function(const std::vector<std::pair<std::string, double>> &args,
                       double domain_min, double domain_max)
    {
        // default arguments initialized as empty
        x_arg_.name("x");
        loc_.name("loc");
        lower_.name("min");
        upper_.name("max");

        // allocate other arguments
        arguments_.reserve(args.size());
        for (const auto &[key, val] : args)
            if (key == "x")
            {
                x_arg_.init_numeric(val);
            }
            else if (key == "loc")
            {
                loc_.init_numeric(val);
            }
            else
            {
                // distinguish empty defaults
                if (na::is(val))
                    arguments_.emplace_back();
                else
                    arguments_.emplace_back(val);
                arguments_.back().name(key);
            }
        arguments_.shrink_to_fit();

        // assign domain if any, otherwise leave empty
        if (!na::is(domain_min))
        {
            lower_.init_numeric(domain_min);
            below_mode_ = Mode::Interval::Closed;
        }
        if (!na::is(domain_max))
        {
            upper_.init_numeric(domain_max);
            above_mode_ = Mode::Interval::Closed;
        }
    }

    /////////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////  Function : Getter  /////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////

    /// Assert argument exists at assignment. Returns it's index.
    Attribute &Function::find_arg_(std::string_view key)
    {
        if (key == "x")
            return x_arg_;
        if (key == "loc")
            return loc_;
        if (key == "min")
            return lower_;
        if (key == "max")
            return upper_;
        for (Attribute &attr : arguments_)
            if (attr.name() == key)
                return attr;

        // raise error for unnamed argument
        if (na::is(key))
            throw std::invalid_argument(
                "Unnamed argument passed to " + err_name() + ".");

        // else format the error showing options
        std::vector<std::string> options = {"x", "loc", "min", "max"};
        for (const Attribute &attr : arguments_)
            options.emplace_back(attr.name());
        throw std::invalid_argument(
            "Unused argument '" + std::string(key) + "' in " + err_name() + "." +
            str::find_closest(key, options));
    }

    /// Wheter return value is always constant.
    /// UPDATE: at the moment also a mutable range makes the function non constant.
    bool Function::is_constant() const noexcept
    {
        if (!x_arg_.is_constant())
            return false;

        // variable arguments
        for (const Attribute &attr : arguments_)
            if (!attr.is_constant())
                return false;

        // variable range
        if (below_mode_ != Mode::Interval::None && !lower_.is_constant())
            return false;
        if (above_mode_ != Mode::Interval::None && !upper_.is_constant())
            return false;

        // otherwise constant
        return true;
    }

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
                std::string("Missing non default argument") +
                (undefined_args.size() > 1 ? "s " : " ") +
                str::format(undefined_args, str::MAX_PRINT, false) +
                " for " + err_name() + ".");
    }

    /////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////  Function : set range  //////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////

    bool Function::helper_init_range_min_(double min_val, bool is_min_open)
    {
        if (lower_.is_empty())
            throw std::runtime_error(
                "Forbidden initialization of minimum value for " + err_name() + ".");

        const double old = lower_.value();
        if (old < min_val)
            throw std::invalid_argument(
                "Invoked invalid min value " + str::format(old) + " for " + err_name() + ".");

        if (old == min_val && upper_.is_empty())
            return false;

        helper_init_standard_range_(is_min_open ? old == min_val : false, false);
        return true;
    }

    void Function::helper_init_standard_range_(bool is_min_open, bool is_max_open)
    {
        if (lower_.is_empty())
            below_mode_ = Mode::Interval::None;
        else if (is_min_open)
            below_mode_ = Mode::Interval::Open;
        else
            below_mode_ = Mode::Interval::Closed;

        if (upper_.is_empty())
            above_mode_ = Mode::Interval::None;
        else if (is_max_open)
            above_mode_ = Mode::Interval::Open;
        else
            above_mode_ = Mode::Interval::Closed;

        // check for consistency if both range limits are provided
        if (!lower_.is_callable() || !upper_.is_callable())
            return;

        bool is_wrong = lower_.value() > upper_.value();
        if (below_mode_ == Mode::Interval::Open && below_mode_ == above_mode_)
            is_wrong = is_wrong || lower_.value() == upper_.value();
        if (is_wrong)
            throw std::invalid_argument(
                "Generated an inconsistent range " + str_range() + " for " + err_name() + ".");
    }

    /////////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////  Function : show  ///////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////

    /// Return formatted function summary.
    std::string Function::str_summary(bool show_x_arg) const
    {
        std::string sms = std::string(name()) + "(";
        if (show_x_arg)
            sms += std::string(x_arg_.name()) + "=" + x_arg_.str_val();

        for (size_t c = 0; c < arguments_.size(); ++c)
        {
            if (c > 0 || show_x_arg)
                sms += ", ";
            sms += std::string(arguments_[c].name()) + "=" + arguments_[c].str_val();
        }
        return sms + ") : " + str_range();
    }

    /// Return formatted range.
    std::string Function::str_range() const
    {
        const std::string &str_min = lower_.str_val();
        const std::string &str_max = upper_.str_val();

        std::string output;

        switch (below_mode_)
        {
        case Mode::Interval::None:
            output = "(-Inf";
            break;
        case Mode::Interval::Open:
            output = "(" + str_min;
            break;
        case Mode::Interval::Closed:
            output = "[" + str_min;
            break;
        default:
            __builtin_unreachable();
        }
        switch (above_mode_)
        {
        case Mode::Interval::None:
            return output + ", Inf)";
        case Mode::Interval::Open:
            return output + ", " + str_max + ")";
        case Mode::Interval::Closed:
            return output + ", " + str_max + "]";
        }
        __builtin_unreachable();
    }

    /////////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////  Function : Arguments  ///////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////

    void Function::init_attr(std::string_view key, const Attribute &attr)
    {
        // assert that function attribute are implicitely callable
        try
        {
            helper_assert_valid_attribute_(key, attr);
        }
        catch (const std::exception &e)
        {
            throw std::invalid_argument(
                "Impossible to associate " + attr.err_name() +
                " as argument '" + std::string(key) + "' of " + err_name() + ".\n" +
                std::string(e.what()));
        }

        // transfer attribute
        find_arg_(key).transfer(attr);

        // apply location argument shift
        if (key == "loc")
            deloc_mode_ = Mode::Location::Active;
    }

    /// Assert that proposed function attribute is valid.
    void Function::helper_assert_valid_attribute_(std::string_view key, const Attribute &attr)
    {
        if (attr.is_function())
        {
            const auto &undefined_args = attr.function()->get_missing_argument();
            if (!undefined_args.empty())
                throw std::invalid_argument(
                    "Non default argument" + std::string(undefined_args.size() > 1 ? "s " : " ") +
                    str::format(undefined_args) +
                    " of " + attr.err_name() + " cannot be missing.");

            if (!attr.function()->is_callable_implicitly())
                throw std::invalid_argument(
                    "Argument 'x' of " + attr.err_name() + " cannot be missing.");
        }
        else if (attr.is_constant())
        {
            double x = attr.value();
            if (na::is(x) || (key == "x" && !is_valid(x)))
                throw std::invalid_argument(
                    "Value " + str::format(x) + " for argument '" + std::string(key) +
                    "' is outside of range " + str_range() + ".");
        }
    }

    void Function::optimize()
    {
        for (Attribute &attr : arguments_)
            attr.optimize();
        x_arg_.optimize();
        loc_.optimize();
        lower_.optimize();
        upper_.optimize();
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////  *** Function : implementations ***  /////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

namespace
{
    using namespace hsp;

#define MAKE_REAL_TO_INTEGER_FUNCTION(Name, call)                                    \
    class Name final : public Function                                               \
    {                                                                                \
    public:                                                                          \
        HOUSEPI_OVERRIDE_FUNCTION_CLONE_METHOD(Name)                                 \
        HOUSEPI_OVERRIDE_NAME_METHOD(call)                                           \
        bool is_discrete() const noexcept override final { return true; }            \
                                                                                     \
    protected:                                                                       \
        double raw_function_(double x) const override final { return std::call(x); } \
    };                                                                               \
    HOUSEPI_REGISTER_FUNCTION(Name, Arithmetic)

    MAKE_REAL_TO_INTEGER_FUNCTION(Round_, round)
    MAKE_REAL_TO_INTEGER_FUNCTION(Ceil_, ceil)
    MAKE_REAL_TO_INTEGER_FUNCTION(Floor_, floor)

    class Minus_ final : public Function
    {
    public:
        HOUSEPI_OVERRIDE_FUNCTION_CLONE_METHOD(Minus_)
        HOUSEPI_OVERRIDE_NAME_METHOD(minus)
        bool is_discrete() const noexcept override final { return x_arg_.is_discrete(); }

    protected:
        double raw_function_(double x) const override final { return -x; }
    };
    HOUSEPI_REGISTER_FUNCTION(Minus_, Arithmetic)

    class Subtraction_ final : public Function
    {
    public:
        Subtraction_() : Function({{"x2", na::DOUBLE}}) {}
        HOUSEPI_OVERRIDE_FUNCTION_CLONE_METHOD(Subtraction_)
        HOUSEPI_OVERRIDE_NAME_METHOD(sub)
        bool is_discrete() const noexcept override final { return arguments_[0].is_discrete() && x_arg_.is_discrete(); }

    protected:
        double raw_function_(double x) const override final { return x - arg(0); }
    };
    HOUSEPI_REGISTER_FUNCTION(Subtraction_, Arithmetic)

    class Division_ final : public Function
    {
    public:
        Division_() : Function({{"x2", na::DOUBLE}}) {}
        HOUSEPI_OVERRIDE_FUNCTION_CLONE_METHOD(Division_)
        HOUSEPI_OVERRIDE_NAME_METHOD(div)

    protected:
        double raw_function_(double x) const override final { return x / arg(0); }
    };
    HOUSEPI_REGISTER_FUNCTION(Division_, Arithmetic)

    class Power_ final : public Function
    {
    public:
        Power_() : Function({{"x2", na::DOUBLE}}) {}
        HOUSEPI_OVERRIDE_FUNCTION_CLONE_METHOD(Power_)
        HOUSEPI_OVERRIDE_NAME_METHOD(pow)
        bool is_discrete() const noexcept override final { return Function::is_discrete() || (arguments_[0].is_discrete() && x_arg_.is_discrete()); }

    protected:
        double raw_function_(double x) const override final { return std::pow(x, arg(0)); }
        double raw_log_func_(double x) const override final { return arg(0) * std::log(x); }
        double raw_primitive_(double x) const override final
        {
            if (std::abs(arg(0) + 1.) < 1.e-10)
                return std::log(std::abs(x));
            return std::pow(x, arg(0) + 1.) / (arg(0) + 1.);
        }
        double raw_log_ratio_(double xo, double xn) const override final { return arg(0) * std::log(xn / xo); }
    };
    HOUSEPI_REGISTER_FUNCTION(Power_, Transcendental)

    class SquareRoot_ final : public Function
    {
    public:
        SquareRoot_() : Function({}, 0.) {}
        HOUSEPI_OVERRIDE_FUNCTION_CLONE_METHOD(SquareRoot_)
        HOUSEPI_OVERRIDE_NAME_METHOD(sqrt)
        void init_truncation() override final { helper_init_range_min_(); }

    protected:
        double raw_function_(double x) const override final { return std::sqrt(x); }
        double raw_log_func_(double x) const override final { return 0.5 * std::log(x); }
        double raw_primitive_(double x) const override final { return std::pow(x, 1.5) / 1.5; }
        double raw_log_ratio_(double xo, double xn) const override final { return raw_log_func_(xn) - raw_log_func_(xo); }
    };
    HOUSEPI_REGISTER_FUNCTION(SquareRoot_, Transcendental)

    class Exponential_ final : public Function
    {
    public:
        HOUSEPI_OVERRIDE_FUNCTION_CLONE_METHOD(Exponential_)
        HOUSEPI_OVERRIDE_NAME_METHOD(exp)

    protected:
        double raw_function_(double x) const override final { return std::exp(x); }
        double raw_log_func_(double x) const override final { return x; }
        double raw_primitive_(double x) const override final { return std::exp(x); }
        double raw_log_ratio_(double xo, double xn) const override final { return xn - xo; }
    };
    HOUSEPI_REGISTER_FUNCTION(Exponential_, Transcendental)

    class Logarithm_ final : public Function
    {
    public:
        Logarithm_() : Function({}, 0.) { below_mode_ = Mode::Interval::Open; }
        HOUSEPI_OVERRIDE_FUNCTION_CLONE_METHOD(Logarithm_)
        HOUSEPI_OVERRIDE_NAME_METHOD(log)
        void init_truncation() override final { helper_init_range_min_(0., true); }

    protected:
        double raw_function_(double x) const override final { return std::log(x); }
        double raw_log_func_(double x) const override final { return std::log(std::log(x)); }
        double raw_primitive_(double x) const override final { return x * std::log(x) - x; }
        double raw_log_ratio_(double xo, double xn) const override final { return raw_log_func_(xn) - raw_log_func_(xo); }
    };
    HOUSEPI_REGISTER_FUNCTION(Logarithm_, Transcendental)

    class Sine_ final : public Function
    {
    public:
        HOUSEPI_OVERRIDE_FUNCTION_CLONE_METHOD(Sine_)
        HOUSEPI_OVERRIDE_NAME_METHOD(sin)

    protected:
        double raw_function_(double x) const override final { return std::sin(x); }
        double raw_log_func_(double x) const override final { return std::log(std::sin(x)); }
        double raw_primitive_(double x) const override final { return -std::cos(x); }
        double raw_log_ratio_(double xo, double xn) const override final { return raw_log_func_(xn) - raw_log_func_(xo); }
    };
    HOUSEPI_REGISTER_FUNCTION(Sine_, Transcendental)

    class Cosine_ final : public Function
    {
    public:
        HOUSEPI_OVERRIDE_FUNCTION_CLONE_METHOD(Cosine_)
        HOUSEPI_OVERRIDE_NAME_METHOD(cos)

    protected:
        double raw_function_(double x) const override final { return std::cos(x); }
        double raw_log_func_(double x) const override final { return std::log(std::cos(x)); }
        double raw_primitive_(double x) const override final { return std::sin(x); }
        double raw_log_ratio_(double xo, double xn) const override final { return raw_log_func_(xn) - raw_log_func_(xo); }
    };
    HOUSEPI_REGISTER_FUNCTION(Cosine_, Transcendental)

    class Tangent_ final : public Function
    {
    public:
        HOUSEPI_OVERRIDE_FUNCTION_CLONE_METHOD(Tangent_)
        HOUSEPI_OVERRIDE_NAME_METHOD(tan)

    protected:
        double raw_function_(double x) const override final { return std::tan(x); }
        double raw_log_func_(double x) const override final { return std::log(std::tan(x)); }
        double raw_log_ratio_(double xo, double xn) const override final { return raw_log_func_(xn) - raw_log_func_(xo); }
    };
    HOUSEPI_REGISTER_FUNCTION(Tangent_, Transcendental)

    /////////////////////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////  VariadicFunction : class  ///////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////

    /// Auxialiary class to deal with variadic functions.
    class VariadicFunction : public Function
    {
    public:
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
        Attribute &find_arg_(std::string_view key) override final
        {
            if (key == "x")
                return x_arg_;
            if (key == "loc")
                return loc_;
            if (key == "min")
                return lower_;
            if (key == "max")
                return upper_;

            for (Attribute &attr : arguments_)
                if (!na::is(attr.name()) && attr.name() == key)
                    return attr;

            arguments_.emplace_back();
            return arguments_.back();
        }

        // If any constant, erase all other constants arguments.
        void remove_if_not_first_const(Attribute *first_const)
        {
            if (!first_const)
                return;

            arguments_.erase(
                std::remove_if(arguments_.begin(), arguments_.end(),
                               [&first_const](Attribute &attr)
                               {
                                   return attr.is_constant() && &attr != first_const;
                               }),
                arguments_.end());
        }
    };

#define MAKE_VARIADIC_ALGEBRAIC(Class_, name_, symbol)                              \
    class Class_ final : public VariadicFunction                                    \
    {                                                                               \
    public:                                                                         \
        Class_() : VariadicFunction() {}                                            \
        HOUSEPI_OVERRIDE_FUNCTION_CLONE_METHOD(Class_)                              \
        HOUSEPI_OVERRIDE_NAME_METHOD(name_)                                         \
                                                                                    \
        void optimize() override                                                    \
        {                                                                           \
            Function::optimize();                                                   \
            Attribute *first_const = nullptr;                                       \
            if (x_arg_.is_constant())                                               \
                first_const = &x_arg_;                                              \
            for (Attribute &attr : arguments_)                                      \
                if (attr.is_constant())                                             \
                {                                                                   \
                    if (first_const)                                                \
                        first_const->set(attr.value() symbol first_const->value()); \
                    else                                                            \
                        first_const = &attr;                                        \
                }                                                                   \
            remove_if_not_first_const(first_const);                                 \
        }                                                                           \
                                                                                    \
    protected:                                                                      \
        double raw_function_(double x) const override final                         \
        {                                                                           \
            for (const Attribute &attr : arguments_)                                \
                x = x symbol attr.value();                                          \
            return x;                                                               \
        }                                                                           \
    };
    MAKE_VARIADIC_ALGEBRAIC(Summation_, sum, +)
    HOUSEPI_REGISTER_FUNCTION(Summation_, Arithmetic)
    MAKE_VARIADIC_ALGEBRAIC(Product_, prod, *)
    HOUSEPI_REGISTER_FUNCTION(Product_, Arithmetic)

#define MAKE_VARIADIC_MIN_MAX(Class_, name_, symbol)             \
    class Class_ final : public VariadicFunction                 \
    {                                                            \
    public:                                                      \
        Class_() : VariadicFunction() {}                         \
        HOUSEPI_OVERRIDE_FUNCTION_CLONE_METHOD(Class_)           \
        HOUSEPI_OVERRIDE_NAME_METHOD(name_)                      \
                                                                 \
        void optimize() override                                 \
        {                                                        \
            Function::optimize();                                \
            Attribute *first_const = nullptr;                    \
            if (x_arg_.is_constant())                            \
                first_const = &x_arg_;                           \
            for (Attribute &attr : arguments_)                   \
                if (attr.is_constant())                          \
                {                                                \
                    if (first_const)                             \
                    {                                            \
                        const double new_val = attr.value();     \
                        if (new_val symbol first_const->value()) \
                            first_const->set(new_val);           \
                    }                                            \
                    else                                         \
                        first_const = &attr;                     \
                }                                                \
            remove_if_not_first_const(first_const);              \
        }                                                        \
                                                                 \
    protected:                                                   \
        double raw_function_(double x) const override final      \
        {                                                        \
            for (const Attribute &attr : arguments_)             \
            {                                                    \
                const double new_val = attr.value();             \
                if (new_val symbol x)                            \
                    x = new_val;                                 \
            }                                                    \
            return x;                                            \
        }                                                        \
    };
    MAKE_VARIADIC_MIN_MAX(Minimum_, min, <)
    HOUSEPI_REGISTER_FUNCTION(Minimum_, Arithmetic)
    MAKE_VARIADIC_MIN_MAX(Maximum_, max, >)
    HOUSEPI_REGISTER_FUNCTION(Maximum_, Arithmetic)
}

/////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////  Attribute : init product  ///////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

/// UPDATE: this will be moved when Attribute func beomes unique ptr
void Attribute::init_product(const std::vector<const Attribute *> &other_ptrs)
{
    if (other_ptrs.size() == 1)
    {
        init_pointer(other_ptrs[0]);
        return;
    }

    auto fun = std::make_unique<Product_>();
    Attribute aux;
    aux.init_pointer(other_ptrs[0]);
    fun->init_attr("x", aux);
    for (size_t i = 1; i < other_ptrs.size(); ++i)
    {
        aux.init_pointer(other_ptrs[i]);
        fun->init_attr("", aux);
    }

    init_function(std::move(fun));
}
