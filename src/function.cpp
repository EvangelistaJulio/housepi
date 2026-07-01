/* housepi - Epidemiological Analysis of Household Outbreaks.
 *
 * Copyright © 2024-2026 Institut Pasteur (Paris), CNRS and INSERM.
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
#include "tools/parser.h"
#include "tools/string.h"
#include "tools/types.h"

#include "function.h"
#include "registry.h"

namespace housepi
{
    namespace
    {
        std::string err_type_name_(ErrType type) noexcept
        {
            switch (type)
            {
            case ErrType::Numeric:
                return "numeric";
            case ErrType::Pointer:
                return "pointer";
            case ErrType::Function:
                return "function";
            case ErrType::Parameter:
                return "parameter";
            case ErrType::Augmented:
                return "augmented";
            case ErrType::Observation:
                return "observation";
            case ErrType::HouseholdObservation:
                return "household observation";
            case ErrType::GroupRate:
                return "group rate";
            case ErrType::Data:
                return "data";
            default: // case ErrType::Empty:
                return "empty";
            }
        }

        /// Deduplicated string pool. Shared across models, thus append only.
        /// NOTE: intern should be run only during serial model construction.
        struct NamePool_
        {
        public:
            /// Intern a name and return its index.
            /// NOTE: index 0 == empty name (NA).
            NameIndex intern(std::string_view name)
            {
                auto it = lut_.find(name); // heterogeneous lookup: no std::string built
                if (it != lut_.end())
                    return it->second;

                NameIndex id = static_cast<NameIndex>(pool_.size());
                pool_.emplace_back(name);
                lut_.emplace(name, id);
                return id;
            }

            /// Resolve an index to its copied name.
            std::string at(NameIndex id) const { return pool_[id]; }

        private:
            std::vector<std::string> pool_ = {""}; /// interned names
            dict<NameIndex> lut_ = {{"", 0}};      /// reverse name-to-index
        };

        /// The single process-global name pool.
        NamePool_ &name_pool_()
        {
            static NamePool_ pool;
            return pool;
        }
    }

    /////////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////  *** Attribute : class ***  /////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////

    Attribute::Attribute() = default;
    Attribute::~Attribute() = default;

    Attribute::Attribute(const Attribute &other)
        : mode_(other.mode_),
          err_type_(other.err_type_),
          is_user_const_(other.is_user_const_),
          name_id_(other.name_id_),
          value_(other.value_),
          pAttr_(other.pAttr_),
          pFunc_(other.pFunc_ ? other.pFunc_->clone() : nullptr) {}

    Attribute::Attribute(Attribute &&) = default;

    Attribute &Attribute::operator=(Attribute other)
    {
        swap(other);
        return *this;
    }

    void Attribute::swap(Attribute &other) noexcept
    {
        std::swap(mode_, other.mode_);
        std::swap(err_type_, other.err_type_);
        std::swap(is_user_const_, other.is_user_const_);
        std::swap(name_id_, other.name_id_);
        std::swap(value_, other.value_);
        std::swap(pAttr_, other.pAttr_);
        std::swap(pFunc_, other.pFunc_);
    }

    bool Attribute::operator==(const Attribute &other) const noexcept
    {
        if (mode_ != other.mode_)
            return false;
        if (is_constant() != other.is_constant())
            return false;
        if (is_discrete() != other.is_discrete())
            return false;

        switch (mode_)
        {
        case Mode::Numeric:
            return value_ == other.value_;
        case Mode::Pointer:
            return pAttr_ == other.pAttr_;
        case Mode::Function:
            return pFunc_ == other.pFunc_;
        default:
            return true;
        }
    }

    std::string Attribute::name() const { return name_pool_().at(name_id_); }
    std::string Attribute::err_name() const { return err_type_name_(err_type_) + " '" + name() + "'"; }

    std::string Attribute::str_val() const
    {
        switch (mode_)
        {
        case Mode::Numeric:
            return str::format(value_);
        case Mode::Pointer:
            return pAttr_->str_val();
        case Mode::Function:
            if (pFunc_->is_callable_implicitly())
                return str::format(pFunc_->evaluate());
            return "<>";
        default: // case Mode::Empty:
            return str::format(na::DOUBLE);
        }
    }

#ifndef NDEBUG
    std::string Attribute::str_summary() const
    {
        std::string sms = err_name() + " ";

        sms += (is_function() && !pFunc_->is_callable_implicitly()) ? "<> | "
                                                                    : str_val() + " | ";
        sms += is_discrete() ? "discrete"
                             : "real";
        sms += is_constant() ? " constant"
                             : " variable";

        switch (mode_)
        {
        case Mode::Numeric:
            return sms + " numeric";
        case Mode::Pointer:
            return sms + " pointer -> " + pAttr_->err_name();
        case Mode::Function:
            return sms + " function " + pFunc_->str_summary();
        default: // case Mode::Empty:
            return sms + " empty";
        }
    }
#endif

    bool Attribute::is_constant() const noexcept
    {
        switch (mode_)
        {
        case Mode::Numeric:
            return is_user_const_;
        case Mode::Pointer:
            return pAttr_->is_constant();
        case Mode::Function:
            return pFunc_->is_constant();
        default: // case Mode::Empty:
            return false;
        }
    }

    bool Attribute::is_discrete() const noexcept
    {
        switch (mode_)
        {
        case Mode::Numeric:
            return value_ == std::floor(value_);
        case Mode::Pointer:
            return pAttr_->is_discrete();
        case Mode::Function:
            return pFunc_->is_discrete();
        default: // case Mode::Empty:
            return false;
        }
    }

    bool Attribute::is_callable() const noexcept
    {
        switch (mode_)
        {
        case Mode::Numeric:
            return !na::is(value_);
        case Mode::Pointer:
            return pAttr_->is_callable();
        case Mode::Function:
            return pFunc_->is_callable_implicitly();
        default: // case Mode::Empty:
            return false;
        }
    }

    /////////////////////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////  Attribute : setter //////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////

    void Attribute::name(const std::string &txt) { name_id_ = name_pool_().intern(txt); }

    void Attribute::init_numeric(double x)
    {
        if (err_type_ == ErrType::Empty)
            err_type_ = ErrType::Numeric;

        mode_ = Mode::Numeric;
        value_ = x;
        pAttr_ = nullptr;
        pFunc_ = nullptr;
        is_user_const_ = true;
    }

    void Attribute::init_pointer(const Attribute *y)
    {
        if (!y)
            throw std::invalid_argument(
                "Tried to assign empty pointer to " + err_name() + ".");

        if (err_type_ == ErrType::Empty)
            err_type_ = ErrType::Pointer;

        mode_ = Mode::Pointer;
        value_ = na::DOUBLE;
        pAttr_ = y;
        pFunc_ = nullptr;
        is_user_const_ = false;
    }

    void Attribute::init_function(std::unique_ptr<Function> f)
    {
        // check the function is different from the current one
        if (!f)
            throw std::invalid_argument(
                "Tried to assign empty function to " + err_name() + ".");
        if (f == pFunc_)
            throw std::invalid_argument(
                "Tried to assign " + err_name() + " as its own argument.");

        if (err_type_ == ErrType::Empty)
            err_type_ = ErrType::Function;
        if (name_id_ == 0)
            name(std::string(f->name()));
        mode_ = Mode::Function;
        value_ = na::DOUBLE;
        pAttr_ = nullptr;
        pFunc_ = std::move(f);
        is_user_const_ = false;
    }

    void Attribute::transfer(const Attribute &other)
    {
        switch (other.mode_)
        {
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
            init_function(other.function()->clone());
            return;
        default: // case Mode::Empty:
            // NOTE: nothing to transfer if other is pure empty
            if (other.err_type_ != ErrType::Empty)
                init_pointer(other.pointer());
            return;
        }
    }
    void Attribute::optimize()
    {
        if (is_constant())
            init_numeric(value());
        else if (is_user_const_)
            throw std::invalid_argument(
                "Impossible to set " + err_name() + " as constant. Check its dependencies.");

        // reduce function indirection
        if (is_pointer() && pointer()->is_function())
            init_function(pointer()->function()->clone());
        if (is_function())
            pFunc_->optimize();
    }

    void Attribute::init_func_attr(std::string_view key, const Attribute &attr)
    {
        switch (mode_)
        {
        case Mode::Function:
            pFunc_->init_attr(key, attr);
            return;
        default:
            throw std::runtime_error(
                "Attempted to initialize function attribute of " + err_name() + ".");
        }
    }

    void Attribute::raise_set_error_() const
    {
        std::string sms = "Tried to change non numeric ";
        if (na::is(name()))
            sms += "empty " + err_type_name_(err_type_);
        else
            sms += err_name();
        throw std::runtime_error(sms + ".");
    }

    double Attribute::value() const
    {
        switch (mode_)
        {
        case Mode::Numeric:
            return value_;
        case Mode::Pointer:
            return pAttr_->value(); // UPDATE: not inlineable
        case Mode::Function:
            return pFunc_->evaluate();
        default: // case Mode::Empty:
            return na::DOUBLE;
        }
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

        // Fill canonical parametrization at construction
        register_alt_("canonical", args);

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

    [[noreturn]] void Function::raise_undefined_call_err_(const std::string &call_name) const
    {
        throw std::invalid_argument(
            "Called undefined method '" + call_name + "' of " + err_name() + ".");
    }

    /////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////  Function : Getter  /////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////

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

        // search non standard arguments
        // NOTE: first hit picks/confirms the active alternative set
        for (size_t i = 0; i < alts_.size(); ++i)
            for (Attribute &attr : alts_[i].args)
                if (attr.name() == key)
                {
                    const std::uint8_t new_idx = static_cast<std::uint8_t>(i);
                    if (is_user_alt_ && alt_idx_ != new_idx)
                        throw std::invalid_argument(
                            "Attempt to mix parametrizations '" +
                            alts_[alt_idx_].name + "' and '" + alts_[new_idx].name + "' in " + err_name() +
                            " (argument '" + std::string(key) + "').");
                    alt_idx_ = new_idx;
                    is_user_alt_ = true;
                    return attr;
                }

        // raise error for unnamed argument
        if (na::is(key))
            throw std::invalid_argument(
                "Unnamed argument passed to " + err_name() + ".");

        // else format the error showing options from every parametrization
        std::vector<std::string> options = {"x", "loc", "min", "max"};
        for (const AltArgSet &alt : alts_)
            for (const Attribute &attr : alt.args)
                options.emplace_back(attr.name());
        throw std::invalid_argument(
            "Unused argument '" + std::string(key) + "' in " + err_name() + "." +
            str::find_closest(key, options));
    }

    Attribute &Function::find_or_append_arg_(std::string_view key)
    {
        if (key == "x")
            return x_arg_;
        if (key == "loc")
            return loc_;
        if (key == "min")
            return lower_;
        if (key == "max")
            return upper_;

        std::vector<Attribute> &args = arguments_();
        for (Attribute &attr : args)
            if (!na::is(attr.name()) && attr.name() == key)
                return attr;

        args.emplace_back();
        return args.back();
    }

    /// Wheter return value is always constant.
    /// UPDATE: at the moment also a mutable range makes the function non constant.
    bool Function::is_constant() const noexcept
    {
        if (!x_arg_.is_constant())
            return false;

        // variable arguments (only the active parametrization matters)
        for (const Attribute &attr : arguments_())
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
        for (const Attribute &attr : arguments_())
            if (attr.is_empty())
                return false;
        return true;
    }

    /// Collect names of missing arguments (from the active parametrization).
    std::vector<std::string> Function::get_missing_argument() const
    {
        std::vector<std::string> undefined_args = {};
        for (const Attribute &attr : arguments_())
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
    ///////////////////////////////////////  Function : Setter  /////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////

    void Function::register_alt_(const std::string &alt_name,
                                 const std::vector<std::pair<std::string, double>> &args)
    {
        // collect already-known names for collision detection
        std::vector<std::string> existing_names = {"x", "loc", "min", "max"};
        for (const AltArgSet &alt : alts_)
            for (const Attribute &attr : alt.args)
                existing_names.emplace_back(attr.name());

        // also check alt set names are unique
        for (const AltArgSet &alt : alts_)
            if (alt.name == alt_name)
                throw std::invalid_argument(
                    "Duplicate parametrization name '" + alt_name + "' in " + err_name() + ".");

        AltArgSet alt;
        alt.name = alt_name;
        alt.args.reserve(args.size());
        for (const auto &[key, val] : args)
        {
            for (const std::string &n : existing_names)
                if (n == key)
                    throw std::invalid_argument(
                        "Argument name '" + key + "' from parametrization '" + alt_name +
                        "' collides with an existing argument in " + err_name() + ".");
            existing_names.emplace_back(key);

            if (na::is(val))
                alt.args.emplace_back();
            else
                alt.args.emplace_back(val);
            alt.args.back().name(key);
        }
        alts_.emplace_back(std::move(alt));
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

#ifndef NDEBUG
    /// Return formatted function summary.
    std::string Function::str_summary(bool show_x_arg) const
    {
        std::string sms = std::string(name()) + "(";
        if (show_x_arg)
            sms += std::string(x_arg_.name()) + "=" + x_arg_.str_val();

        const std::vector<Attribute> &active = arguments_();
        for (size_t c = 0; c < active.size(); ++c)
        {
            if (c > 0 || show_x_arg)
                sms += ", ";
            sms += std::string(active[c].name()) + "=" + active[c].str_val();
        }
        return sms + ") : " + str_range();
    }
#endif

    /// Return formatted range.
    std::string Function::str_range() const
    {
        const std::string &str_min = lower_.str_val();
        const std::string &str_max = upper_.str_val();

        std::string output;

        switch (below_mode_)
        {
        case Mode::Interval::Open:
            output = "(" + str_min;
            break;
        case Mode::Interval::Closed:
            output = "[" + str_min;
            break;
        default: // case Mode::Interval::None:
            output = "(-Inf";
            break;
        }
        switch (above_mode_)
        {
        case Mode::Interval::Open:
            return output + ", " + str_max + ")";
        case Mode::Interval::Closed:
            return output + ", " + str_max + "]";
        default: // case Mode::Interval::None:
            return output + ", Inf)";
        }
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
        for (Attribute &attr : arguments_())
            attr.optimize();
        x_arg_.optimize();
        loc_.optimize();
        lower_.optimize();
        upper_.optimize();
    }

    /////////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////// *** FunctionRegistry : class ***  /////////////////////////////
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

    std::unique_ptr<Function> FunctionRegistry::create(std::string_view opt, bool is_prim) const
    {
        auto it = factory_map_.find(opt);
        if (it == factory_map_.end())
            return nullptr;

        auto output = it->second();
        if (is_prim)
            output->set_eval_as_prim();
        // evaluate as function by default

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

        // FIXME: this should be an error
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
}

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////  *** Function : implementations ***  /////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

HOUSEPI_DEFINE_REGISTRY_INSTANCE_(Custom)
HOUSEPI_DEFINE_REGISTRY_INSTANCE_(Arithmetic)
HOUSEPI_DEFINE_REGISTRY_INSTANCE_(Transcendental)
HOUSEPI_DEFINE_REGISTRY_INSTANCE_(ContinuosDistrib)
HOUSEPI_DEFINE_REGISTRY_INSTANCE_(DiscreteDistrib)
HOUSEPI_DEFINE_REGISTRY_INSTANCE_(InfecProfile)

namespace
{
    using namespace housepi;

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
        bool is_discrete() const noexcept override final { return arguments_(0).is_discrete() && x_arg_.is_discrete(); }

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
        bool is_discrete() const noexcept override final { return Function::is_discrete() || (arguments_(0).is_discrete() && x_arg_.is_discrete()); }

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
            for (const Attribute &attr : arguments_())
                if (!attr.is_discrete())
                    return false;
            return true;
        }

    protected:
        Attribute &find_arg_(std::string_view key) override final { return find_or_append_arg_(key); }

        /// NA-skipping fold over x and all arguments.
        /// Returns NA if no non-NA operand is found.
        template <class Op>
        double reduce_na_rm_(double x, Op op) const
        {
            bool has = !na::is(x);
            double acc = x;
            for (const Attribute &attr : arguments_())
            {
                const double v = attr.value();
                if (na::is(v))
                    continue;
                acc = has ? op(acc, v) : v;
                has = true;
            }
            return has ? acc : na::DOUBLE;
        }

        /// Reduce all single constant inside the first constant argument.
        /// Constant NA arguments are skipped.
        template <class Op>
        void reduce_const_na_rm_(Op op)
        {
            Attribute *first_const = nullptr;
            if (x_arg_.is_constant() && !na::is(x_arg_.value()))
                first_const = &x_arg_;
            for (Attribute &attr : arguments_())
                if (attr.is_constant() && !na::is(attr.value()))
                {
                    if (first_const)
                        first_const->set(op(first_const->value(), attr.value()));
                    else
                        first_const = &attr;
                }

            if (!first_const)
                return;

            std::vector<Attribute> &args = arguments_();
            args.erase(
                std::remove_if(args.begin(), args.end(),
                               [&first_const](Attribute &attr)
                               {
                                   return attr.is_constant() && &attr != first_const;
                               }),
                args.end());
        }
    };

#define MAKE_VARIADIC_ALGEBRAIC(Class_, name_, symbol)             \
    class Class_ final : public VariadicFunction                   \
    {                                                              \
    public:                                                        \
        Class_() : VariadicFunction() {}                           \
        HOUSEPI_OVERRIDE_FUNCTION_CLONE_METHOD(Class_)             \
        HOUSEPI_OVERRIDE_NAME_METHOD(name_)                        \
                                                                   \
        void optimize() override                                   \
        {                                                          \
            Function::optimize();                                  \
            reduce_const_na_rm_(                                   \
                [](double a, double b) { return a symbol b; });    \
        }                                                          \
                                                                   \
    protected:                                                     \
        double raw_function_(double x) const override final        \
        {                                                          \
            return reduce_na_rm_(                                  \
                x, [](double a, double b) { return a symbol b; }); \
        }                                                          \
    };
    MAKE_VARIADIC_ALGEBRAIC(Summation_, sum, +)
    HOUSEPI_REGISTER_FUNCTION(Summation_, Arithmetic)
    MAKE_VARIADIC_ALGEBRAIC(Product_, prod, *)
    HOUSEPI_REGISTER_FUNCTION(Product_, Arithmetic)

#define MAKE_VARIADIC_MIN_MAX(Class_, name_, symbol)                       \
    class Class_ final : public VariadicFunction                           \
    {                                                                      \
    public:                                                                \
        Class_() : VariadicFunction() {}                                   \
        HOUSEPI_OVERRIDE_FUNCTION_CLONE_METHOD(Class_)                     \
        HOUSEPI_OVERRIDE_NAME_METHOD(name_)                                \
                                                                           \
        void optimize() override                                           \
        {                                                                  \
            Function::optimize();                                          \
            reduce_const_na_rm_(                                           \
                [](double a, double b) { return b symbol a ? b : a; });    \
        }                                                                  \
                                                                           \
    protected:                                                             \
        double raw_function_(double x) const override final                \
        {                                                                  \
            return reduce_na_rm_(                                          \
                x, [](double a, double b) { return b symbol a ? b : a; }); \
        }                                                                  \
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
