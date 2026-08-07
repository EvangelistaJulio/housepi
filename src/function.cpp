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
    std::string Attribute::err_name() const { return err_type_name(err_type_) + " " + msg::val(name()); }

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
            msg::stop(
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
            msg::stop(
                "Tried to assign empty function to " + err_name() + ".");
        if (f == pFunc_)
            msg::stop(
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
    void Attribute::finalize()
    {
        if (is_constant())
            init_numeric(value());
        else if (is_user_const_)
            msg::stop(
                "Impossible to set " + err_name() + " as constant. Check its dependencies.");

        // reduce function indirection
        if (is_pointer() && pointer()->is_function())
            init_function(pointer()->function()->clone());
        if (is_function())
            pFunc_->finalize();
    }

    void Attribute::init_func_attr(std::string_view key, const Attribute &attr)
    {
        switch (mode_)
        {
        case Mode::Function:
            pFunc_->init_attr(key, attr);
            return;
        default:
            msg::stop(
                "Attempted to initialize function attribute of " + err_name() + ".");
        }
    }

    void Attribute::raise_set_error_() const
    {
        std::string sms = "Tried to change non numeric ";
        if (na::is(name()))
            sms += "empty " + err_type_name(err_type_);
        else
            sms += err_name();
        msg::stop(sms + ".");
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

    std::string Function::err_name() const
    {
        return err_type_name(err_family_) + " " + msg::val(name());
    }

    [[noreturn]] void Function::raise_undefined_call_err_(const std::string &call_name) const
    {
        msg::stop(
            "Called undefined method " + msg::val(call_name) + " of " + err_name() + ".");
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
                        msg::stop(
                            "Attempt to mix parametrizations " +
                            msg::val(alts_[alt_idx_].name) + " and " + msg::val(alts_[new_idx].name) + " in " + err_name() +
                            " (argument " + msg::arg(std::string(key)) + ").");
                    alt_idx_ = new_idx;
                    is_user_alt_ = true;
                    return attr;
                }

        // raise error for unnamed argument
        if (na::is(key))
            msg::stop(
                "Unnamed argument passed to " + err_name() + ".");

        // else format the error showing options from every parametrization
        std::vector<std::string> options = {"x", "loc", "min", "max"};
        for (const AltArgSet &alt : alts_)
            for (const Attribute &attr : alt.args)
                options.emplace_back(attr.name());
        msg::stop(
            "Unused argument " + msg::arg(std::string(key)) + " in " + err_name() + "." +
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
        args.back().name(std::string(key));
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
            msg::stop(
                std::string("Missing non default argument") +
                (undefined_args.size() > 1 ? "s " : " ") +
                str::format(undefined_args, false) +
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
                msg::stop(
                    "Duplicate parametrization name " + msg::val(alt_name) + " in " + err_name() + ".");

        AltArgSet alt;
        alt.name = alt_name;
        alt.args.reserve(args.size());
        for (const auto &[key, val] : args)
        {
            for (const std::string &n : existing_names)
                if (n == key)
                    msg::stop(
                        "Argument name " + msg::arg(key) + " from parametrization " + msg::val(alt_name) +
                        " collides with an existing argument in " + err_name() + ".");
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
            msg::stop(
                "Forbidden initialization of minimum value for " + err_name() + ".");

        const double old = lower_.value();
        if (old < min_val)
            msg::stop(
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
            msg::stop(
                "Generated an inconsistent range " + str_range() + " for " + err_name() + ".");
    }

    /////////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////  Function : show  ///////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////

    /// Return formatted range.
    std::string Function::str_range() const
    {
        std::string str_min = lower_.str_val();
        std::string str_max = upper_.str_val();
        if (deloc_mode_ == Mode::Location::Active && loc_.is_callable())
        {
            if (lower_.is_callable())
                str_min = str::format(relocate_(lower_.value()));
            if (upper_.is_callable())
                str_max = str::format(relocate_(upper_.value()));
        }

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
            msg::stop(
                "Impossible to associate " + attr.err_name() +
                " as argument " + msg::arg(std::string(key)) + " of " + err_name() + ".\n" +
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
                msg::stop(
                    "Non default argument" + std::string(undefined_args.size() > 1 ? "s " : " ") +
                    str::format(undefined_args) +
                    " of " + attr.err_name() + " cannot be missing.");

            if (!attr.function()->is_callable_implicitly())
                msg::stop(
                    "Argument " + msg::arg("x") + " of " + attr.err_name() + " cannot be missing.");
        }
        else if (attr.is_constant())
        {
            double x = attr.value();
            if (na::is(x) || (key == "x" && !is_valid(x)))
                msg::stop(
                    "Value " + str::format(x) + " for argument " + msg::arg(std::string(key)) +
                    " is outside of range " + str_range() + ".");
        }
    }

    void Function::finalize()
    {
        for (Attribute &attr : arguments_())
            attr.finalize();
        x_arg_.finalize();
        loc_.finalize();
        lower_.finalize();
        upper_.finalize();
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

    namespace
    {
        /// Resolve a registry to its singleton instance.
        const FunctionRegistry &registry_instance_(Registry r)
        {
            switch (r)
            {
            case Registry::Custom:
                return CustomRegistry::instance();
            case Registry::Arithmetic:
                return ArithmeticRegistry::instance();
            case Registry::Transcendental:
                return TranscendentalRegistry::instance();
            case Registry::Continuos:
                return ContinuosDistribRegistry::instance();
            case Registry::Discrete:
                return DiscreteDistribRegistry::instance();
            case Registry::InfecProfile:
                return InfecProfileRegistry::instance();
            }
            msg::stop("Unknown registry.");
        }

        /// Resolve registries to their singleton instances, preserving lookup order.
        std::vector<const FunctionRegistry *> collect_registries_(const std::vector<Registry> &registries,
                                                                  ErrType err_type)
        {
            if (registries.empty())
                msg::stop(
                    "No registries were provided for " + err_type_name(err_type) + ".");

            std::vector<const FunctionRegistry *> registry_ptrs = {};
            registry_ptrs.reserve(registries.size());
            for (Registry r : registries)
                registry_ptrs.push_back(&registry_instance_(r));

            return registry_ptrs;
        }
    }

    std::vector<std::string> registry_options(const std::vector<Registry> &registries,
                                              ErrType err_type)
    {
        std::vector<std::string> options = {};
        for (const FunctionRegistry *ptr : collect_registries_(registries, err_type))
        {
            const std::vector<std::string> y = ptr->get_options();
            options.reserve(options.size() + y.size());
            options.insert(options.end(), y.begin(), y.end());
        }
        return str::unique(options);
    }

    std::unique_ptr<Function> compose_func_registries(std::string_view opt,
                                                      ErrType err_type,
                                                      const std::vector<Registry> &registries,
                                                      bool is_prim)
    {
        // FIXME: this should be an error
        if (na::is(opt))
            return nullptr;

        for (const FunctionRegistry *ptr : collect_registries_(registries, err_type))
            if (ptr->count(opt))
                return ptr->create(opt, is_prim);

        msg::stop(
            "Unknown " + err_type_name(err_type) + " " + msg::val(std::string(opt)) + "." +
            str::find_closest(opt, registry_options(registries, err_type)));
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
    //////////////////////////////////  ReductionFunction : class  //////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////

    /// Auxialiary class to deal with reduction functions.
    class ReductionFunction : public Function
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
    class Class_ final : public ReductionFunction                  \
    {                                                              \
    public:                                                        \
        Class_() : ReductionFunction() {}                          \
        HOUSEPI_OVERRIDE_FUNCTION_CLONE_METHOD(Class_)             \
        HOUSEPI_OVERRIDE_NAME_METHOD(name_)                        \
                                                                   \
        void finalize() override                                   \
        {                                                          \
            Function::finalize();                                  \
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
    class Class_ final : public ReductionFunction                          \
    {                                                                      \
    public:                                                                \
        Class_() : ReductionFunction() {}                                  \
        HOUSEPI_OVERRIDE_FUNCTION_CLONE_METHOD(Class_)                     \
        HOUSEPI_OVERRIDE_NAME_METHOD(name_)                                \
                                                                           \
        void finalize() override                                           \
        {                                                                  \
            Function::finalize();                                          \
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
