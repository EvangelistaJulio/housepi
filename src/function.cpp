/* housepi - Epidemiological Analysis of Household Outbreaks.
 * Copyright 2024-2026 Institut Pasteur (Paris), France.
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#include <algorithm>
#include <cctype>
#include <cmath>
#include <memory>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

#include "tools/dictionary.h"
#include "tools/messages.h"
#include "tools/nan.h"
#include "tools/parser.h"
#include "tools/string.h"
#include "tools/vector.h"
#include "function.h"
#include "registry/functions.h"

namespace housepi
{
    namespace
    {
        /// Deduplicated string pool, shared across models, hence append-only.
        /// @note Intern only during serial model construction.
        struct NamePool_
        {
        public:
            /// Intern `name`, returning its index.
            /// * index 0 -> empty name, NA.
            NameIdx intern(std::string_view name)
            {
                // heterogeneous lookup: no std::string built
                auto it = lut_.find(name);
                if (it != lut_.end())
                    return it->second;

                NameIdx id = static_cast<NameIdx>(pool_.size());
                pool_.emplace_back(name);
                lut_.emplace(name, id);
                return id;
            }

            /// Return a copy of the name at index `id`.
            std::string at(NameIdx id) const { return pool_[id]; }

        private:
            /// Interned names.
            std::vector<std::string> pool_ = {""};
            /// Index of each name.
            dict<NameIdx> lut_ = {{"", 0}};
        };

        /// Return the process-global name pool.
        NamePool_ &name_pool_()
        {
            static NamePool_ pool;
            return pool;
        }

        /// Operator spellings, each with the function it builds.
        constexpr std::pair<std::string_view, std::string_view> operators_[] = {
            {"+", "sum"}, {"-", "minus"}, {"*", "prod"}, {"/", "div"}, {"^", "pow"}};
    }

    //////////////////////////////  Interval arithmetic  //////////////////////////////

    namespace
    {
        /// Combine openness for a sum: a bound is reached only where both operands reach theirs.
        Interval both_(Interval a, Interval b) noexcept
        {
            if (a == Interval::Unbounded || b == Interval::Unbounded)
                return Interval::Unbounded;
            return a == Interval::Closed && b == Interval::Closed ? Interval::Closed : Interval::Open;
        }

        /// Combine openness for a tied extremum: reached as soon as either operand reaches it.
        /// @note Both sides must be bounded.
        Interval either_(Interval a, Interval b) noexcept
        {
            return a == Interval::Closed || b == Interval::Closed ? Interval::Closed : Interval::Open;
        }

        /// Integers only when both lattices are.
        Lattice both_integers_(Lattice a, Lattice b) noexcept
        {
            return a == Lattice::Integers && b == Lattice::Integers ? Lattice::Integers : Lattice::Doubles;
        }

        /// Sum of codomains `a` and `b`.
        Codomain add_(const Codomain &a, const Codomain &b)
        {
            Codomain c;
            c.lattice = both_integers_(a.lattice, b.lattice);
            c.below = both_(a.below, b.below);
            if (c.below != Interval::Unbounded)
                c.inf = a.inf + b.inf;
            c.above = both_(a.above, b.above);
            if (c.above != Interval::Unbounded)
                c.sup = a.sup + b.sup;
            return c;
        }

        /// Minimum of codomains `a` and `b`.
        Codomain min_(const Codomain &a, const Codomain &b)
        {
            Codomain c;
            c.lattice = both_integers_(a.lattice, b.lattice);

            // the smaller infimum wins, a tie is reached as soon as either reaches it
            if (a.below == Interval::Unbounded || b.below == Interval::Unbounded)
                c.below = Interval::Unbounded;
            else if (a.inf == b.inf)
                c.inf = a.inf, c.below = either_(a.below, b.below);
            else
                c.inf = std::min(a.inf, b.inf), c.below = a.inf < b.inf ? a.below : b.below;

            // the smaller supremum wins, an unbounded operand never provides it
            if (a.above == Interval::Unbounded && b.above == Interval::Unbounded)
                c.above = Interval::Unbounded;
            else if (a.above == Interval::Unbounded)
                c.sup = b.sup, c.above = b.above;
            else if (b.above == Interval::Unbounded)
                c.sup = a.sup, c.above = a.above;
            else if (a.sup == b.sup)
                c.sup = a.sup, c.above = either_(a.above, b.above);
            else
                c.sup = std::min(a.sup, b.sup), c.above = a.sup < b.sup ? a.above : b.above;
            return c;
        }

        /// Maximum of codomains `a` and `b`.
        Codomain max_(const Codomain &a, const Codomain &b)
        {
            Codomain c;
            c.lattice = both_integers_(a.lattice, b.lattice);

            // the larger infimum wins, a tie is reached only where both reach it
            if (a.below == Interval::Unbounded && b.below == Interval::Unbounded)
                c.below = Interval::Unbounded;
            else if (a.below == Interval::Unbounded)
                c.inf = b.inf, c.below = b.below;
            else if (b.below == Interval::Unbounded)
                c.inf = a.inf, c.below = a.below;
            else if (a.inf == b.inf)
                c.inf = a.inf, c.below = both_(a.below, b.below);
            else
                c.inf = std::max(a.inf, b.inf), c.below = a.inf > b.inf ? a.below : b.below;

            // the larger supremum wins, an unbounded operand carries the result
            if (a.above == Interval::Unbounded || b.above == Interval::Unbounded)
                c.above = Interval::Unbounded;
            else if (a.sup == b.sup)
                c.sup = a.sup, c.above = either_(a.above, b.above);
            else
                c.sup = std::max(a.sup, b.sup), c.above = a.sup > b.sup ? a.above : b.above;
            return c;
        }

        /// Smallest interval holding codomains `a` and `b`.
        Codomain hull_(const Codomain &a, const Codomain &b)
        {
            Codomain c = min_(a, b);
            const Codomain upper = max_(a, b);
            c.sup = upper.sup, c.above = upper.above;
            return c;
        }

        /// Negation of codomain `a`.
        Codomain negate_(const Codomain &a)
        {
            Codomain c;
            c.lattice = a.lattice;
            if (a.above != Interval::Unbounded)
                c.inf = -a.sup, c.below = a.above;
            if (a.below != Interval::Unbounded)
                c.sup = -a.inf, c.above = a.below;
            return c;
        }

        /// Whether no value of `a` is negative.
        bool is_non_negative_(const Codomain &a) { return a.below != Interval::Unbounded && a.inf >= 0.; }
        /// Whether no value of `a` is positive.
        bool is_non_positive_(const Codomain &a) { return a.above != Interval::Unbounded && a.sup <= 0.; }
        /// Whether `x` is a bound reached at zero.
        bool is_reached_zero_(double x, Interval side) { return x == 0. && side == Interval::Closed; }

        /// Product of non-negative codomains `a` and `b`.
        /// A bound is reached where both factors reach theirs, or where either reaches zero.
        Codomain mul_non_negative_(const Codomain &a, const Codomain &b)
        {
            Codomain c;
            c.lattice = both_integers_(a.lattice, b.lattice);

            c.inf = a.inf * b.inf;
            c.below = is_reached_zero_(a.inf, a.below) || is_reached_zero_(b.inf, b.below)
                          ? Interval::Closed
                          : both_(a.below, b.below);

            // a factor held at zero holds the product there
            if (is_reached_zero_(a.sup, a.above) || is_reached_zero_(b.sup, b.above))
                c.sup = 0., c.above = Interval::Closed;
            else if ((c.above = both_(a.above, b.above)) != Interval::Unbounded)
                c.sup = a.sup * b.sup;
            return c;
        }

        /// Product of codomains `a` and `b`.
        /// * factor changing sign -> split at zero, products joined;
        /// * non-positive factor -> negated around the non-negative product.
        Codomain mul_(const Codomain &a, const Codomain &b)
        {
            // negative and positive parts, zero reached by both
            const auto split = [](Codomain neg)
            {
                Codomain pos = neg;
                neg.sup = pos.inf = 0.;
                neg.above = pos.below = Interval::Closed;
                return std::pair{neg, pos};
            };
            if (!is_non_negative_(a) && !is_non_positive_(a))
            {
                const auto [neg, pos] = split(a);
                return hull_(mul_(neg, b), mul_(pos, b));
            }
            if (!is_non_negative_(b) && !is_non_positive_(b))
            {
                const auto [neg, pos] = split(b);
                return hull_(mul_(a, neg), mul_(a, pos));
            }

            if (is_non_negative_(a) && is_non_negative_(b))
                return mul_non_negative_(a, b);
            if (is_non_negative_(a))
                return negate_(mul_non_negative_(a, negate_(b)));
            if (is_non_negative_(b))
                return negate_(mul_non_negative_(negate_(a), b));
            return mul_non_negative_(negate_(a), negate_(b));
        }

        /// Exponential of codomain `a`.
        Codomain exp_(const Codomain &a)
        {
            Codomain c;
            // the exponential approaches zero without ever reaching it
            c.inf = 0., c.below = Interval::Open;
            if (a.below != Interval::Unbounded)
                c.inf = std::exp(a.inf), c.below = a.below;
            if (a.above != Interval::Unbounded)
                c.sup = std::exp(a.sup), c.above = a.above;
            return c;
        }

        /// Logarithm of codomain `a`.
        Codomain log_(const Codomain &a)
        {
            Codomain c;
            // the logarithm falls away without bound as its argument approaches zero
            if (a.below != Interval::Unbounded && a.inf > 0.)
                c.inf = std::log(a.inf), c.below = a.below;
            if (a.above != Interval::Unbounded && a.sup > 0.)
                c.sup = std::log(a.sup), c.above = a.above;
            return c;
        }

        /// Reciprocal of codomain `a`.
        /// * divisor reaching zero or changing sign -> unbounded, two rays are not an interval;
        /// * divisor approaching zero -> unbounded on that side.
        Codomain recip_(const Codomain &a)
        {
            if (!a.is_positive())
                return is_non_positive_(a) && !is_reached_zero_(a.sup, a.above) ? negate_(recip_(negate_(a)))
                                                                                : Codomain();
            Codomain c;
            c.inf = 0., c.below = Interval::Open;
            if (a.above != Interval::Unbounded)
                c.inf = 1. / a.sup, c.below = a.above;
            if (a.inf > 0.)
                c.sup = 1. / a.inf, c.above = a.below;
            return c;
        }

        /// Codomain of number `v`: the interval holding only itself.
        Codomain numeric_codomain_(double v) noexcept
        {
            if (na::is(v))
                return Codomain();

            Codomain c;
            c.inf = c.sup = v;
            c.below = c.above = Interval::Closed;
            c.lattice = v == std::floor(v) ? Lattice::Integers : Lattice::Doubles;
            return c;
        }
    }

    std::string Codomain::str() const
    {
        const std::string lo = below == Interval::Unbounded ? "(-Inf" : (below == Interval::Open ? "(" : "[") + str::format(inf);
        const std::string hi = above == Interval::Unbounded ? "Inf)" : str::format(sup) + (above == Interval::Open ? ")" : "]");
        return lo + ", " + hi;
    }

    bool Codomain::is_within(const Codomain &outer) const noexcept
    {
        // a side holds where the outer side is unbounded, or bounds every value on it
        // an unbounded side reads as open at -Inf or Inf
        const double lo = below == Interval::Unbounded ? -inf::DOUBLE : inf;
        const double hi = above == Interval::Unbounded ? inf::DOUBLE : sup;
        const bool is_low_within =
            outer.below == Interval::Unbounded ||
            lo > outer.inf || (lo == outer.inf && (below != Interval::Closed || outer.below == Interval::Closed));
        const bool is_high_within =
            outer.above == Interval::Unbounded ||
            hi < outer.sup || (hi == outer.sup && (above != Interval::Closed || outer.above == Interval::Closed));
        return is_low_within && is_high_within && (!outer.on_integers() || on_integers());
    }

    //////////////////////////////  Attribute  //////////////////////////////

    /// Read the attribute through its storage mode, one branch per mode.
#define MAKE_ATTR_MODE(numeric, pointer, function, empty) \
    switch (mode_)                                        \
    {                                                     \
    case Mode::Numeric:                                   \
        return numeric;                                   \
    case Mode::Pointer:                                   \
        return pointer;                                   \
    case Mode::Function:                                  \
        return function;                                  \
    default: /* case Mode::Empty: */                      \
        return empty;                                     \
    }

    Attribute::Attribute() = default;
    Attribute::~Attribute() = default;

    Attribute::Attribute(const Attribute &other)
        : mode_(other.mode_),
          err_type_(other.err_type_),
          is_user_const_(other.is_user_const_),
          name_id_(other.name_id_),
          value_(other.value_),
          attr_(other.attr_),
          func_(other.func_ ? other.func_->clone() : nullptr) {}

    Attribute::Attribute(Attribute &&) noexcept = default;

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
        std::swap(attr_, other.attr_);
        std::swap(func_, other.func_);
    }

    bool Attribute::operator==(const Attribute &other) const noexcept
    {
        if (mode_ != other.mode_)
            return false;
        if (is_constant() != other.is_constant())
            return false;
        if (is_discrete() != other.is_discrete())
            return false;

        MAKE_ATTR_MODE(value_ == other.value_,
                       attr_ == other.attr_,
                       func_ == other.func_,
                       true);
    }

    std::string Attribute::name() const { return name_pool_().at(name_id_); }
    std::string Attribute::str() const { return str::format(err_type_) + " " + msg::val(name()); }

    std::string Attribute::str_source() const
    {
        const Attribute *source = pointer();
        if (source->is_function())
            return msg::arg(source->name()) + " definition";
        return is_pointer() ? source->str() : str::format(value_);
    }

    bool Attribute::is_constant() const noexcept
    {
        MAKE_ATTR_MODE(is_user_const_,
                       attr_->is_constant(),
                       func_->is_constant(),
                       false);
    }

    Codomain Attribute::codomain() const noexcept
    {
        MAKE_ATTR_MODE(numeric_codomain_(value_),
                       attr_->codomain(),
                       func_->codomain(),
                       Codomain());
    }
    bool Attribute::is_na() const noexcept
    {
        MAKE_ATTR_MODE(na::is(value_),
                       attr_->is_na(),
                       false,
                       true);
    }
    bool Attribute::is_callable() const noexcept
    {
        MAKE_ATTR_MODE(!na::is(value_),
                       attr_->is_callable(),
                       func_->is_callable_implicitly(),
                       false);
    }

    //////////////////////////////  Attribute : setter  //////////////////////////////

    void Attribute::name(const std::string &txt) { name_id_ = name_pool_().intern(txt); }

    void Attribute::init_numeric(double x)
    {
        if (err_type_ == ErrType::Empty)
            err_type_ = ErrType::Numeric;

        mode_ = Mode::Numeric;
        value_ = x;
        attr_ = nullptr;
        func_ = nullptr;
        is_user_const_ = true;
    }

    void Attribute::init_pointer(const Attribute *y)
    {
        if (err_type_ == ErrType::Empty)
            err_type_ = ErrType::Pointer;

        mode_ = Mode::Pointer;
        value_ = na::DOUBLE;
        attr_ = y;
        func_ = nullptr;
        is_user_const_ = false;
    }

    void Attribute::init_function(std::unique_ptr<Function> f)
    {
        if (err_type_ == ErrType::Empty)
            err_type_ = ErrType::Function;
        if (name_id_ == 0)
            name(std::string(f->name()));
        mode_ = Mode::Function;
        value_ = na::DOUBLE;
        attr_ = nullptr;
        func_ = std::move(f);
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
        // case Mode::Empty:
        default:
            // nothing to transfer from a pure empty
            if (other.err_type_ != ErrType::Empty)
                init_pointer(other.pointer());
            return;
        }
    }
    void Attribute::finalize()
    {
        if (is_constant())
            init_numeric(value());

        // reduce function indirection
        if (is_pointer() && pointer()->is_function())
            init_function(pointer()->function()->clone());
        if (is_function())
            func_->finalize();
    }

    void Attribute::init_func_attr(std::string_view key, const Attribute &attr)
    {
        func_->init_attr(key, attr);
    }

    //////////////////////////////  Function  //////////////////////////////

    /// Name the standard arguments, register the canonical parametrization, close the given bounds.
    /// * arguments come as name-value pairs, to be accessed by index.
    /// @todo Construction is the only place setting the domain bounds.
    /// @todo Check control over a pointer argument leaving the domain.
    Function::Function(const Function::default_args &args,
                       double domain_min, double domain_max)
    {
        // standard arguments initialized as empty
        x_arg_.name("x");
        domain_.loc.name("loc");
        domain_.min.name("min");
        domain_.max.name("max");

        // canonical parametrization
        register_alt_("canonical", args);

        // assign domain if any, otherwise leave empty
        if (!na::is(domain_min))
        {
            domain_.min.init_numeric(domain_min);
            domain_.below = Interval::Closed;
        }
        if (!na::is(domain_max))
        {
            domain_.max.init_numeric(domain_max);
            domain_.above = Interval::Closed;
        }
    }

    std::string Function::str() const
    {
        std::string label = name();
        if (err_family_ == ErrType::Operator)
            for (const auto &[symbol, func] : operators_)
                if (func == label)
                    label = symbol;
        return str::format(err_family_) + " " + msg::val(label);
    }

    [[noreturn]] void Function::raise_undefined_call_err_(const std::string &call_name) const
    {
        msg::stop(
            "Called undefined method " + msg::val(call_name) + " of " + str() + ".");
    }

    //////////////////////////////  Function : getter  //////////////////////////////

    Attribute *Function::find_standard_arg_(std::string_view key)
    {
        for (Attribute *attr : standard_args_())
            if (attr->name() == key)
                return attr;
        return nullptr;
    }

    Attribute &Function::find_arg_(std::string_view key)
    {
        if (Attribute *attr = find_standard_arg_(key))
            return *attr;

        // non-standard arguments: the first hit picks or confirms the active parametrization
        for (size_t i = 0; i < alts_.size(); ++i)
            for (Attribute &attr : alts_[i].args)
                if (attr.name() == key)
                {
                    const std::uint8_t new_idx = static_cast<std::uint8_t>(i);
                    if (is_user_alt_ && alt_idx_ != new_idx)
                        msg::stop(
                            "Attempt to mix parametrizations " +
                            msg::val(alts_[alt_idx_].name) + " and " + msg::val(alts_[new_idx].name) + " in " + str() +
                            " (argument " + msg::arg(std::string(key)) + ").");
                    alt_idx_ = new_idx;
                    is_user_alt_ = true;
                    return attr;
                }

        // raise error for unnamed argument
        if (na::is(key))
            msg::stop(
                "Unnamed argument passed to " + str() + ".");

        // keys a definition call can carry:
        std::vector<std::string> options;
        for (size_t i = 0; i < alts_.size(); ++i)
            // parametrization the user has picked excludes the others
            if (!is_user_alt_ || i == static_cast<size_t>(alt_idx_))
                for (const Attribute &attr : alts_[i].args)
                    options.emplace_back(attr.name());
        for (const Attribute *attr : standard_args_())
            options.emplace_back(attr->name());
        msg::stop(
            "Unused argument " + msg::arg(std::string(key)) + " in " + str() + "." +
            str::find_closest(key, options, "Show its arguments with " + registry::browse_call(name()) + "."));
    }

    Attribute &Function::append_arg_(std::string_view key)
    {
        if (Attribute *attr = find_standard_arg_(key))
            return *attr;

        std::vector<Attribute> &args = arguments_();
        args.emplace_back();
        args.back().name(std::string(key));
        return args.back();
    }

    bool Function::has_constant_args_() const noexcept
    {
        // variable arguments (only the active parametrization matters)
        for (const Attribute &attr : arguments_())
            if (!attr.is_constant())
                return false;

        // variable domain
        if (domain_.below != Interval::Unbounded && !domain_.min.is_constant())
            return false;
        if (domain_.above != Interval::Unbounded && !domain_.max.is_constant())
            return false;

        // otherwise constant
        return true;
    }

    bool Function::is_callable_explicitly() const noexcept { return get_missing_argument().empty(); }

    std::vector<std::string> Function::get_missing_argument() const
    {
        std::vector<std::string> undefined_args = {};
        for (const Attribute &attr : arguments_())
            if (attr.is_empty())
                undefined_args.emplace_back(attr.name());
        // a domain bound under a requirement is an argument of the law
        for (const auto &[attr, req] : arg_requirements_())
            if ((attr == &domain_.min || attr == &domain_.max) && attr->is_empty())
                undefined_args.emplace_back(attr->name());
        return undefined_args;
    }

    /// Throw if not callable explicitly.
    void Function::assert_is_callable_explicitly() const
    {
        const std::vector<std::string> undefined_args = get_missing_argument();
        if (!undefined_args.empty())
            msg::stop(
                std::string("Missing non default argument") +
                (undefined_args.size() > 1 ? "s " : " ") +
                str::format(undefined_args, false) +
                " for " + str() + ".");
    }

    namespace
    {
        /// Raise on an attribute that must hold an integer value but does not.
        [[noreturn]] void raise_not_discrete_(const Attribute &attr, std::string sms)
        {
            if (attr.is_numeric())
                sms += "\nChange " + msg::arg(attr.name()) + " (" + attr.str_source() + ") to an integer.";
            else
                sms += "\nChange " + attr.str_source() + " to discrete.";
            msg::stop(sms);
        }

        /// Words for the bounds of `c`, as read by the user.
        /// * side open at -Inf or Inf -> "finite", once.
        std::string str_bounds_(const Codomain &c)
        {
            std::string sms;
            bool is_finite = false;
            if (c.below != Interval::Unbounded)
            {
                if (std::isinf(c.inf))
                    is_finite = true;
                else
                    sms = (c.below == Interval::Open ? "greater than " : "at least ") + str::format(c.inf);
            }
            if (c.above != Interval::Unbounded)
            {
                if (std::isinf(c.sup))
                    is_finite = true;
                else
                    sms += (sms.empty() ? "" : " and ") +
                           std::string(c.above == Interval::Open ? "less than " : "at most ") + str::format(c.sup);
            }
            if (is_finite)
                sms += (sms.empty() ? "" : " and ") + std::string("finite");
            return sms;
        }

        /// Raise on an argument that can take a value outside its requirement `req`.
        /// * number -> its value is set outside;
        /// * constant source -> the source holds a value outside;
        /// * any other source -> it can take values outside.
        [[noreturn]] void raise_outside_(const Attribute &attr, const Codomain &req, const std::string &owner)
        {
            const std::string bounds = str_bounds_(req);
            // a number is named after its argument, a source by itself
            const std::string source = attr.is_numeric() ? msg::arg(attr.name()) : attr.str_source();
            std::string sms = msg::arg(attr.name()) + " of " + owner + " must be " + bounds;
            const bool is_finite_only = bounds == "finite";
            if (attr.is_numeric() || attr.is_constant())
                sms += ".\nChange " + source + " (" + str::format(attr.value()) + ") to " +
                       (is_finite_only ? "a finite value" : "a value " + bounds) + ".";
            else
                sms += ".\n" + source + " can take other values.\nRestrict " + source + " to " +
                       (is_finite_only ? "finite values" : "values " + bounds) + ".";
            msg::stop(sms);
        }
    }

    void Function::assert_location_matches_discrete() const
    {
        // a lattice survives an integer shift, not a fractional one
        // the check asks where the lattice lands, not whether there is one
        // an unbound location holds no value to judge
        const Attribute &loc = domain_.loc;
        if (!is_discretized() || loc.is_empty() || (loc.is_constant() && loc.is_na()) || loc.is_discrete())
            return;

        raise_not_discrete_(
            domain_.loc, msg::arg(domain_.loc.name()) + " definition is incompatible with discrete " + str() + ".");
    }

    void Function::assert_arg_requirements() const
    {
        for (const auto &[attr, req] : arg_requirements_())
        {
            // an empty or unbound argument holds no value to judge
            if (attr->is_empty() || (attr->is_constant() && attr->is_na()))
                continue;

            const Codomain values = attr->codomain();
            if (req.on_integers() && !values.on_integers())
                raise_not_discrete_(*attr, msg::arg(attr->name()) + " of " + str() + " must be discrete.");
            if (!values.is_within(req))
                raise_outside_(*attr, req, str());
        }
    }

    //////////////////////////////  Function : setter  //////////////////////////////

    void Function::register_alt_(const std::string &alt_name,
                                 const Function::default_args &args)
    {
        // known names, for collision detection
        std::vector<std::string> existing_names;
        for (const Attribute *attr : standard_args_())
            existing_names.emplace_back(attr->name());
        for (const AltArgSet &alt : alts_)
            for (const Attribute &attr : alt.args)
                existing_names.emplace_back(attr.name());

        AltArgSet alt;
        alt.name = alt_name;
        alt.args.reserve(args.size());
        for (const auto &[key, val] : args)
        {
            existing_names.emplace_back(key);

            if (na::is(val))
                alt.args.emplace_back();
            else
                alt.args.emplace_back(val);
            alt.args.back().name(key);
        }
        alts_.emplace_back(std::move(alt));
    }

    //////////////////////////////  Function : domain  //////////////////////////////

    bool Domain::init_min(double min_val,
                          bool is_min_open,
                          const std::string &owner)
    {
        if (min.is_empty())
            msg::stop(
                "Forbidden initialization of minimum value for " + owner + ".");

        const double old = min.value();
        if (old < min_val)
            msg::stop(
                "Invoked invalid min value " + str::format(old) + " for " + owner + ".");

        if (is_min_at(min_val) && max.is_empty())
            return false;

        init_standard(is_min_open ? old == min_val : false, false, owner);
        return true;
    }

    void Domain::init_standard(bool is_min_open, bool is_max_open, const std::string &owner)
    {
        if (min.is_empty() || is_min_at(-inf::DOUBLE))
            below = Interval::Unbounded;
        else if (is_min_open)
            below = Interval::Open;
        else
            below = Interval::Closed;

        if (max.is_empty() || is_max_at(inf::DOUBLE))
            above = Interval::Unbounded;
        else if (is_max_open)
            above = Interval::Open;
        else
            above = Interval::Closed;

        assert_consistent(owner);
    }

    bool Domain::is_consistent() const
    {
        // an unread bound delimits nothing to contradict
        if (!min.is_callable() || !max.is_callable())
            return true;

        if (min.value() > max.value())
            return false;
        // a collapsed interval excludes its point under any open bound
        return min.value() < max.value() ||
               (below != Interval::Open && above != Interval::Open);
    }

    void Domain::assert_consistent(const std::string &owner) const
    {
        if (!is_consistent())
            msg::stop("Inconsistent domain " + str() + " for " + owner + ".");
    }

    //////////////////////////////  Function : show  //////////////////////////////

    Codomain Domain::relocated() const
    {
        Codomain c;

        // a side is bounded only where the domain and its bound both bound it
        const Codomain c_min = min.codomain();
        c.below = both_(below, c_min.below);
        if (c.below != Interval::Unbounded)
            c.inf = c_min.inf;

        const Codomain c_max = max.codomain();
        c.above = both_(above, c_max.above);
        if (c.above != Interval::Unbounded)
            c.sup = c_max.sup;

        return loc.is_empty() ? c : add_(c, loc.codomain());
    }

    Codomain Domain::relocated(Lattice lattice) const
    {
        Codomain c = relocated();
        c.lattice = lattice;

        if (lattice == Lattice::Doubles)
        {
            if (c.below == Interval::Closed)
                c.below = Interval::Open;
            if (c.above == Interval::Closed)
                c.above = Interval::Open;
        }
        return c;
    }

    //////////////////////////////  Function : arguments  //////////////////////////////

    void Function::init_attr(std::string_view key, const Attribute &attr)
    {
        // NA constant -> error, an unnamed one named by its argument
        if (attr.is_constant() && attr.is_na())
        {
            const std::string arg = "argument " + msg::arg(std::string(key)) + " of " + str();
            msg::stop((attr.name().empty() ? "Constant " + arg + " is NA"
                                           : attr.str() + " has no value to pass as " + arg) + ".");
        }

        // the attribute must be callable implicitly
        try
        {
            helper_assert_valid_attribute_(key, attr);
        }
        catch (const std::exception &e)
        {
            msg::stop(
                "Invalid argument " + msg::arg(std::string(key)) + " for " + str() + ".\n" + msg::alert(e.what()));
        }

        // transfer attribute
        find_arg_(key).transfer(attr);
    }

    /// Assert argument `key` accepts `attr`.
    void Function::helper_assert_valid_attribute_(std::string_view key,
                                                  const Attribute &attr)
    {
        if (attr.is_function())
        {
            const auto &undefined_args = attr.function()->get_missing_argument();
            if (!undefined_args.empty())
                msg::stop(
                    "Non default argument" + std::string(undefined_args.size() > 1 ? "s " : " ") +
                    str::format(undefined_args, false) +
                    " of " + attr.function()->str() + " cannot be missing.");

            if (!attr.function()->is_callable_implicitly())
                msg::stop(
                    "Argument " + msg::arg(x_arg_.name()) + " of " + attr.function()->str() + " cannot be missing.");
        }
        else if (attr.is_constant())
        {
            double x = attr.value();
            if (key == x_arg_.name() && !is_valid(x))
                msg::stop("Value " + str::format(x) + " is outside of domain " + domain_.str() + ".");
        }
    }

    void Function::finalize()
    {
        for (Attribute &attr : arguments_())
            attr.finalize();
        x_arg_.finalize();
        domain_.loc.finalize();
        domain_.min.finalize();
        domain_.max.finalize();
        trunc_norm_mode_ = TruncNorm::Compute;
        if (trunc_mode_ != Trunc::None && has_constant_args_())
        {
            const double z = trunc_();
            trunc_norm_mode_ = TruncNorm::LookUp;
            lookup_log_trunc_ = std::log(z);
            lookup_trunc_ = z;
        }
    }

    //////////////////////////////  registry::RegFunctions  //////////////////////////////

    namespace registry
    {
        std::vector<std::string> RegFunctions::names() const
        {
            std::vector<std::string> entries;
            entries.reserve(factory_map_.size());
            for (const auto &[key, val] : factory_map_)
                entries.push_back(key);
            return entries;
        }

        std::unique_ptr<Function> RegFunctions::create(std::string_view name, bool is_prim) const
        {
            auto it = factory_map_.find(name);
            if (it == factory_map_.end())
                return nullptr;

            auto output = it->second();
            // evaluates as function unless primitive
            if (is_prim)
                output->set_eval_as_prim();

            return output;
        }

        RegFunctions &instance(Kind kind)
        {
            static RegFunctions custom, arithmetic, transcendental;
            static RegDistributions continuous, discrete, categorical, infec_profile;
            switch (kind)
            {
            case Kind::Custom:
                return custom;
            case Kind::Arithmetic:
                return arithmetic;
            case Kind::Transcendental:
                return transcendental;
            case Kind::Continuous:
                return continuous;
            case Kind::Discrete:
                return discrete;
            case Kind::Categorical:
                return categorical;
            // case Kind::InfecProfile
            default:
                return infec_profile;
            }
        }

        //////  SearchOrder  //////

        bool SearchOrder::contains(std::string_view name) const
        {
            return std::any_of(begin(), end(), [&](Kind kind)
                               { return instance(kind).contains(name); });
        }
        std::vector<std::string> SearchOrder::names() const
        {
            std::vector<std::string> out = {};
            for (Kind kind : kinds)
                vect::append(out, instance(kind).names());
            return str::unique(out);
        }

        SearchOrder distr_search_order(bool discr)
        {
            return {discr ? std::vector<Kind>{Kind::Discrete, Kind::Continuous, Kind::Categorical}
                          : std::vector<Kind>{Kind::Continuous, Kind::Discrete, Kind::Categorical}};
        }
        SearchOrder func_search_order()
        {
            return {{Kind::Custom, Kind::Transcendental, Kind::Arithmetic}};
        }
        SearchOrder func_distr_search_order(bool discr)
        {
            return {vect::concatenate(vect::concatenate(func_search_order().kinds, distr_search_order(discr).kinds),
                                      {Kind::InfecProfile})};
        }
        SearchOrder infec_profile_search_order(bool discr)
        {
            return {discr ? std::vector<Kind>{Kind::InfecProfile, Kind::Discrete, Kind::Continuous}
                          : std::vector<Kind>{Kind::InfecProfile, Kind::Continuous, Kind::Discrete}};
        }

        namespace
        {
            /// R function showing the names `kind` registers.
            std::string browser_(Kind kind)
            {
                switch (kind)
                {
                case Kind::Custom:
                case Kind::Arithmetic:
                case Kind::Transcendental:
                    return "browse_functions";
                case Kind::Continuous:
                case Kind::Discrete:
                case Kind::Categorical:
                    return "browse_distributions";
                // case Kind::InfecProfile
                default:
                    return "browse_infectivity_profiles";
                }
            }

            /// Hint naming the R functions that show the names `order` registers.
            std::string browse_hint_(const SearchOrder &order)
            {
                std::vector<std::string> fns;
                for (Kind kind : order)
                    fns.push_back(browser_(kind));
                fns = str::unique(fns);

                std::string out = "Show the available names with " + msg::run(fns[0] + "()");
                for (size_t i = 1; i < fns.size(); ++i)
                    out += (i + 1 < fns.size() ? ", " : " and ") + msg::run(fns[i] + "()");
                return out + ".";
            }
        }

        std::string browse_call(std::string_view name)
        {
            std::string browser = browser_(Kind::Custom);
            for (Kind kind : {Kind::InfecProfile, Kind::Continuous, Kind::Discrete, Kind::Categorical})
                if (instance(kind).contains(name))
                {
                    browser = browser_(kind);
                    break;
                }
            return msg::run(browser + "(\"" + std::string(name) + "\")");
        }

        std::unique_ptr<Function> create_function(std::string_view name,
                                                  ErrType err_type,
                                                  const SearchOrder &order,
                                                  bool flag)
        {
            // an operator builds the function it spells, named as spelled in errors
            const auto *op = std::find_if(std::begin(operators_), std::end(operators_),
                                          [&](const auto &o)
                                          { return o.first == name; });
            const bool is_operator = op != std::end(operators_);
            if (is_operator)
                name = op->second;

            for (Kind kind : order)
                if (auto output = instance(kind).create(name, flag))
                {
                    if (is_operator)
                        output->set_err_family(ErrType::Operator);
                    return output;
                }

            msg::stop(
                "Unknown " + str::format(err_type) + " " + msg::val(std::string(name)) + "." +
                str::find_closest(name, order.names(), browse_hint_(order)));
        }
    }
}

//////////////////////////////  Function implementations  //////////////////////////////

namespace housepi
{
    namespace
    {
#define MAKE_REAL_TO_INTEGER_FUNCTION(Name, call)                                    \
    class Name final : public Cloneable<Name, Function>                              \
    {                                                                                \
    public:                                                                          \
        static constexpr const char *NAME = #call;                                   \
                                                                                     \
        Codomain codomain() const noexcept override final                            \
        {                                                                            \
            Codomain c = x_arg_.codomain();                                          \
            c.lattice = Lattice::Integers;                                           \
            /* monotone, and an integer value is always reached */                   \
            if (c.below != Interval::Unbounded)                                      \
                c.inf = std::call(c.inf), c.below = Interval::Closed;                \
            if (c.above != Interval::Unbounded)                                      \
                c.sup = std::call(c.sup), c.above = Interval::Closed;                \
            return c;                                                                \
        }                                                                            \
                                                                                     \
    protected:                                                                       \
        double raw_function_(double x) const override final { return std::call(x); } \
    };                                                                               \
    HOUSEPI_REGISTER_FUNCTION(Name, Arithmetic)

        MAKE_REAL_TO_INTEGER_FUNCTION(Round_, round)
        MAKE_REAL_TO_INTEGER_FUNCTION(Ceil_, ceil)
        MAKE_REAL_TO_INTEGER_FUNCTION(Floor_, floor)

        class Minus_ final : public Cloneable<Minus_, Function>
        {
        public:
            static constexpr const char *NAME = "minus";

            Codomain codomain() const noexcept override final { return negate_(x_arg_.codomain()); }

        protected:
            double raw_function_(double x) const override final { return -x; }
        };
        HOUSEPI_REGISTER_FUNCTION(Minus_, Arithmetic)

        class Subtraction_ final : public Cloneable<Subtraction_, Function>
        {
        public:
            static constexpr const char *NAME = "sub";
            Subtraction_() : Cloneable({{"x2", na::DOUBLE}}) {}

            /// Subtracting is adding the negation.
            Codomain codomain() const noexcept override final
            {
                return add_(x_arg_.codomain(), negate_(arguments_(0).codomain()));
            }

        protected:
            double raw_function_(double x) const override final { return x - arg(0); }
        };
        HOUSEPI_REGISTER_FUNCTION(Subtraction_, Arithmetic)

        class Division_ final : public Cloneable<Division_, Function>
        {
        public:
            static constexpr const char *NAME = "div";
            Division_() : Cloneable({{"x2", na::DOUBLE}}) {}

            /// Dividing is multiplying by the reciprocal.
            Codomain codomain() const noexcept override final
            {
                return mul_(x_arg_.codomain(), recip_(arguments_(0).codomain()));
            }

        protected:
            double raw_function_(double x) const override final { return x / arg(0); }
        };
        HOUSEPI_REGISTER_FUNCTION(Division_, Arithmetic)

        class Power_ final : public Cloneable<Power_, Function>
        {
        public:
            static constexpr const char *NAME = "pow";
            Power_() : Cloneable({{"x2", na::DOUBLE}}) {}

            Codomain codomain() const noexcept override final
            {
                const Codomain b = x_arg_.codomain();
                const Codomain e = arguments_(0).codomain();
                // a positive base is the exponential of the exponent times its logarithm
                // a base reaching zero or below needs the sign and parity of the exponent
                Codomain c = b.is_positive() ? exp_(mul_(e, log_(b))) : Codomain();
                // an integer base raised to an integer non-negative power stays an integer
                c.lattice = b.on_integers() && e.on_integers() && e.below != Interval::Unbounded && e.inf >= 0.
                                ? Lattice::Integers
                                : Lattice::Doubles;
                return c;
            }

        protected:
            double raw_function_(double x) const override final { return std::pow(x, arg(0)); }
            double raw_log_func_(double x) const override final { return arg(0) * std::log(x); }
            double raw_primitive_(double x) const override final
            {
                if (std::abs(arg(0) + 1.) < 1.e-10)
                    return std::log(std::abs(x));
                return std::pow(x, arg(0) + 1.) / (arg(0) + 1.);
            }
        };
        HOUSEPI_REGISTER_FUNCTION(Power_, Transcendental)

        class SquareRoot_ final : public Cloneable<SquareRoot_, Function>
        {
        public:
            static constexpr const char *NAME = "sqrt";
            SquareRoot_() : Cloneable({}, 0.) {}

            Codomain codomain() const noexcept override final
            {
                const Codomain a = x_arg_.codomain();
                Codomain c;
                // the domain forbids a negative argument, so the root is never negative
                c.inf = 0., c.below = Interval::Closed;
                if (a.below != Interval::Unbounded && a.inf >= 0.)
                    c.inf = std::sqrt(a.inf), c.below = a.below;
                if (a.above != Interval::Unbounded && a.sup >= 0.)
                    c.sup = std::sqrt(a.sup), c.above = a.above;
                return c;
            }
            void truncate() override final { helper_init_domain_min_(); }

        protected:
            double raw_function_(double x) const override final { return std::sqrt(x); }
            double raw_log_func_(double x) const override final { return 0.5 * std::log(x); }
            double raw_primitive_(double x) const override final { return std::pow(x, 1.5) / 1.5; }
        };
        HOUSEPI_REGISTER_FUNCTION(SquareRoot_, Transcendental)

        class Exponential_ final : public Cloneable<Exponential_, Function>
        {
        public:
            static constexpr const char *NAME = "exp";

            Codomain codomain() const noexcept override final { return exp_(x_arg_.codomain()); }

        protected:
            double raw_function_(double x) const override final { return std::exp(x); }
            double raw_log_func_(double x) const override final { return x; }
            double raw_primitive_(double x) const override final { return std::exp(x); }
        };
        HOUSEPI_REGISTER_FUNCTION(Exponential_, Transcendental)

        class Logarithm_ final : public Cloneable<Logarithm_, Function>
        {
        public:
            static constexpr const char *NAME = "log";
            Logarithm_() : Cloneable({}, 0.) { domain_.below = Interval::Open; }

            Codomain codomain() const noexcept override final { return log_(x_arg_.codomain()); }
            void truncate() override final { helper_init_domain_min_(0., true); }

        protected:
            double raw_function_(double x) const override final { return std::log(x); }
            double raw_log_func_(double x) const override final { return std::log(std::log(x)); }
            double raw_primitive_(double x) const override final { return x * std::log(x) - x; }
        };
        HOUSEPI_REGISTER_FUNCTION(Logarithm_, Transcendental)

        class Sine_ final : public Cloneable<Sine_, Function>
        {
        public:
            static constexpr const char *NAME = "sin";

            Codomain codomain() const noexcept override final
            {
                Codomain c;
                c.inf = -1., c.sup = 1.;
                c.below = c.above = Interval::Closed;
                return c;
            }

        protected:
            double raw_function_(double x) const override final { return std::sin(x); }
            double raw_log_func_(double x) const override final { return std::log(std::sin(x)); }
            double raw_primitive_(double x) const override final { return -std::cos(x); }
        };
        HOUSEPI_REGISTER_FUNCTION(Sine_, Transcendental)

        class Cosine_ final : public Cloneable<Cosine_, Function>
        {
        public:
            static constexpr const char *NAME = "cos";

            Codomain codomain() const noexcept override final
            {
                Codomain c;
                c.inf = -1., c.sup = 1.;
                c.below = c.above = Interval::Closed;
                return c;
            }

        protected:
            double raw_function_(double x) const override final { return std::cos(x); }
            double raw_log_func_(double x) const override final { return std::log(std::cos(x)); }
            double raw_primitive_(double x) const override final { return std::sin(x); }
        };
        HOUSEPI_REGISTER_FUNCTION(Cosine_, Transcendental)

        class Tangent_ final : public Cloneable<Tangent_, Function>
        {
        public:
            static constexpr const char *NAME = "tan";

        protected:
            double raw_function_(double x) const override final { return std::tan(x); }
            double raw_log_func_(double x) const override final { return std::log(std::tan(x)); }
        };
        HOUSEPI_REGISTER_FUNCTION(Tangent_, Transcendental)

        //////////////////////////////  ReductionFunction  //////////////////////////////

        /// Base of reductions over `x` and variadic arguments.
        class ReductionFunction : public Function
        {
        public:
        protected:
            /// Whether the fold always skips operand `attr`.
            ///
            /// * structurally absent, or a literal NA -> skipped;
            /// * pointer -> never skipped, even reading NA before initialization: it answers for its target.
            static bool is_skipped_(const Attribute &attr) noexcept
            {
                return attr.is_empty() || (attr.is_numeric() && attr.is_na());
            }

            Attribute &find_arg_(std::string_view key) override final { return append_arg_(key); }

            /// Fold over `x` and every argument.
            /// * fold mode -> NA skipped;
            /// * otherwise -> any NA operand makes it NA.
            template <class Op>
            double reduce_(double x, Op op) const
            {
                switch (eval_mode_)
                {
                case Evaluate::Fold:
                    return reduce_na_rm_(x, op);
                default:
                    if (na::is(x))
                        return na::DOUBLE;
                    for (const Attribute &attr : arguments_())
                    {
                        const double v = attr.value();
                        if (na::is(v))
                            return na::DOUBLE;
                        x = op(x, v);
                    }
                    return x;
                }
            }

            /// Fold over `x` and every argument, NA skipped.
            /// * no operand left -> NA.
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

            /// Fold codomain() over `x` and every argument, NA skipped.
            /// * no operand folded -> unbounded codomain, where the runtime fold returns NA.
            template <class Op>
            Codomain reduce_codomain_(Op op) const
            {
                bool has = !is_skipped_(x_arg_);
                Codomain acc = has ? x_arg_.codomain() : Codomain();
                for (const Attribute &attr : arguments_())
                {
                    if (is_skipped_(attr))
                        continue;
                    acc = has ? op(acc, attr.codomain()) : attr.codomain();
                    has = true;
                }
                return has ? acc : Codomain();
            }

            /// Fold every constant argument into the first constant one.
            /// * constant NA argument -> skipped in fold mode, else the call stays unfolded to return NA.
            template <class Op>
            void reduce_const_(Op op)
            {
                const auto is_na_const = [](const Attribute &attr)
                { return attr.is_constant() && attr.is_na(); };
                if (eval_mode_ != Evaluate::Fold &&
                    (is_na_const(x_arg_) || std::any_of(arguments_().begin(), arguments_().end(), is_na_const)))
                    return;

                Attribute *first_const = nullptr;
                if (x_arg_.is_constant() && !x_arg_.is_na())
                    first_const = &x_arg_;
                for (Attribute &attr : arguments_())
                    if (attr.is_constant() && !attr.is_na())
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

#define MAKE_VARIADIC_ALGEBRAIC(Class_, name_, symbol, fold_)        \
    class Class_ final : public Cloneable<Class_, ReductionFunction> \
    {                                                                \
    public:                                                          \
        static constexpr const char *NAME = #name_;                  \
                                                                     \
        Codomain codomain() const noexcept override final            \
        {                                                            \
            return reduce_codomain_(fold_);                          \
        }                                                            \
                                                                     \
        void finalize() override                                     \
        {                                                            \
            Function::finalize();                                    \
            reduce_const_(                                           \
                [](double a, double b) { return a symbol b; });      \
        }                                                            \
                                                                     \
    protected:                                                       \
        double raw_function_(double x) const override final          \
        {                                                            \
            return reduce_(                                          \
                x, [](double a, double b) { return a symbol b; });   \
        }                                                            \
    };
        MAKE_VARIADIC_ALGEBRAIC(Summation_, sum, +, add_)
        HOUSEPI_REGISTER_FUNCTION(Summation_, Arithmetic)
        MAKE_VARIADIC_ALGEBRAIC(Product_, prod, *, mul_)
        HOUSEPI_REGISTER_FUNCTION(Product_, Arithmetic)

#define MAKE_VARIADIC_MIN_MAX(Class_, name_, symbol, fold_)                \
    class Class_ final : public Cloneable<Class_, ReductionFunction>       \
    {                                                                      \
    public:                                                                \
        static constexpr const char *NAME = #name_;                        \
                                                                           \
        Codomain codomain() const noexcept override final                  \
        {                                                                  \
            return reduce_codomain_(fold_);                                \
        }                                                                  \
                                                                           \
        void finalize() override                                           \
        {                                                                  \
            Function::finalize();                                          \
            reduce_const_(                                                 \
                [](double a, double b) { return b symbol a ? b : a; });    \
        }                                                                  \
                                                                           \
    protected:                                                             \
        double raw_function_(double x) const override final                \
        {                                                                  \
            return reduce_(                                                \
                x, [](double a, double b) { return b symbol a ? b : a; }); \
        }                                                                  \
    };
        MAKE_VARIADIC_MIN_MAX(Minimum_, min, <, min_)
        HOUSEPI_REGISTER_FUNCTION(Minimum_, Arithmetic)
        MAKE_VARIADIC_MIN_MAX(Maximum_, max, >, max_)
        HOUSEPI_REGISTER_FUNCTION(Maximum_, Arithmetic)
    }

    //////////////////////////////  Attribute : init product  //////////////////////////////

    /// @todo Move next to the Attribute implementation.
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
}
