/* housepi - Epidemiological Analysis of Household Outbreaks.
 *
 * Copyright © 2024-2026 Institut Pasteur (Paris), CNRS and INSERM.
 * Distributed under the terms of the GNU General Public License (GPLv3).
 * See the LICENSE file for details.
 */

#include <algorithm>
#include <stdexcept>
#include <string>
#include <vector>

#include "tools/dependency.h"
#include "tools/keys.h"
#include "tools/messages.h"
#include "tools/nan.h"
#include "tools/operators.h"
#include "tools/parser.h"
#include "tools/random.h"
#include "tools/read.h"
#include "tools/string.h"
#include "tools/types.h"

#include "data_wide.h"
#include "distribution.h"
#include "function.h"
#include "group_rate.h"
#include "individual.h"
#include "parameter.h"
#include "variable.h"

namespace housepi
{
    ///////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////  *** Proposal : class ***  /////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////

    /// Base class for Metropolis proposals.
    class Proposal
    {
    public:
        /// Default constructor.
        Proposal() = default;
        /// Virtual destructor.
        virtual ~Proposal() = default;
        /// Clone this proposal object.
        virtual std::unique_ptr<Proposal> clone() const = 0;

        /// If adaptive, reset to default proposal scale.
        virtual void reset_scale(const Distribution * /*prior*/) {}

        /// Propose a value using random number generator.
        virtual double propose(rnd_gen & /*prng*/, double /*xo*/) const = 0;
        /// Compute log-ratio of proposal densities for Metropolis-Hastings.
        virtual double log_ratio(double /*xo*/, double /*xn*/) const = 0;

        /// Adapt proposal parameters based on observed acceptance rate
        /// and current adaptation ratio.
        virtual void adapt(double /*accept_rate*/, double /*ratio*/) {};

        /// Check if proposal is adaptive.
        virtual bool is_adaptive() const noexcept = 0;

        /// Rebind proposal to a new distribution.
        virtual void rebind_distrib(const Distribution * /*new_distr*/) {}
    };

#define MAKE_PROPOSAL_CLONE_METHOD(Class_) \
    std::unique_ptr<Proposal> clone() const override { return std::make_unique<Class_>(*this); }

    ///////////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////  *** NullProposal_ : class ***  /////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////

    class NullProposal_ final : public Proposal
    {
    public:
        NullProposal_() {}
        MAKE_PROPOSAL_CLONE_METHOD(NullProposal_)

        double propose(rnd_gen & /*prng*/, double xo) const override final { return xo; }
        double log_ratio(double /*xo*/, double /*xn*/) const override final { return 0.; }

        bool is_adaptive() const noexcept override final { return false; }
    };

    ///////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////  *** SampleProposal_ : class ***  ////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////

    class SampleProposal_ final : public Proposal
    {
    public:
        SampleProposal_(const Distribution *distr) : pDistr_(distr) {}
        MAKE_PROPOSAL_CLONE_METHOD(SampleProposal_)

        void rebind_distrib(const Distribution *new_distr) override final { pDistr_ = new_distr; }

        double propose(rnd_gen &prng, double /*xo*/) const override final { return pDistr_->rand(prng); }
        double log_ratio(double xo, double xn) const override final { return pDistr_->log_ratio(xo, xn); }

        bool is_adaptive() const noexcept override final { return false; }

    private:
        const Distribution *pDistr_;
    };

    ///////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////  *** AdaptableProposal : class ***  ///////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////

    /// Adaptable version of Metropolis proposals.
    class AdaptableProposal : public Proposal
    {
    public:
        ///////////////////////////
        //////  CONSTRUCTOR  //////
        ///////////////////////////

        AdaptableProposal(
            const std::string &distr_name,
            double scale,
            bool is_discrete)
            : pDistr_(compose_distr_registries(distr_name, "proposal", {"continuos"})),
              is_adaptive_(na::is(scale))
        {
            if (is_discrete)
                pDistr_->init_discretization();
            pDistr_->init_truncation();
        }
        virtual ~AdaptableProposal() = default;

        AdaptableProposal(const AdaptableProposal &other)
            : pDistr_(other.pDistr_ ? other.pDistr_->clone_distrib() : nullptr),
              is_adaptive_(other.is_adaptive_) {}

        //////////////////////////
        //////  ADAPTATION  //////
        //////////////////////////

        void reset_scale(const Distribution *prior) override final
        {
            if (is_adaptive())
                set_scale_(get_default_scale_(prior));
        }

        void adapt(double accept_rate, double ratio) override final
        {
            set_scale_(get_scale_() * robbins_monro_(accept_rate, ratio));
        }

        /////////////////////
        //////  FLAGS  //////
        /////////////////////

        /// Wheter proposal is adaptive.
        bool is_adaptive() const noexcept override final { return is_adaptive_; }

    protected:
        std::unique_ptr<Distribution> pDistr_; /// distribution associated to proposal
        bool is_adaptive_ = true;              /// whether scale is adaptive

        /// Set default proposal scale.
        virtual double get_default_scale_(const Distribution * /*prior*/) = 0;
        /// Read the current scale from pDistr_.
        virtual double get_scale_() const = 0;
        /// Mechanism to change corresponding scale in pDistr_.
        virtual void set_scale_(double /*scale*/) = 0;

        /// Vanishing Robbins-Monro multiplier from an observed accept rate.
        /// Based on @cite andrieu_tutorial_2008
        static double robbins_monro_(double accept_rate, double ratio) noexcept
        {
            return std::exp((1. - ratio) * (accept_rate - Variable::TARGET_ACCEPT_RATE));
        }
        /// Compute InterQuartile Range (IQR) from prior distribution (Q25, Q75).
        static std::pair<double, double> get_iqr_from_prior_(const Distribution *prior)
        {
            return {prior->quantile(0.25), prior->quantile(0.75)};
        }
        /// Wheter IQR is finite and can be used to compute initial proposal scale.
        static bool is_iqr_finite_(double q25, double q75) noexcept
        {
            return std::isfinite(q25) && std::isfinite(q75) && q75 > q25;
        }
    };

    ///////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////  *** UniformProposal_ : class ***  ////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////

    class UniformProposal_ final : public AdaptableProposal
    {
    public:
        UniformProposal_(double scale, const Distribution *prior)
            : AdaptableProposal("uniform", scale, prior->is_discrete())
        {
            set_scale_(is_adaptive() ? 1. : std::fabs(scale));
        }
        MAKE_PROPOSAL_CLONE_METHOD(UniformProposal_)

        double propose(rnd_gen &prng, double xo) const override final { return xo + pDistr_->rand(prng); }
        double log_ratio(double, double) const override final { return 0.; }

    protected:
        double get_default_scale_(const Distribution *prior) override final
        {
            const auto [q25, q75] = get_iqr_from_prior_(prior);
            if (is_iqr_finite_(q25, q75))
                return q75 - q25;
            return 1.;
        }
        /// NOTE: uniform stores its scale in the bounds, not in the args.
        double get_scale_() const override final { return pDistr_->upper().value(); }
        /// Adapt standard minimum and maximum.
        void set_scale_(double scale) override final
        {
            pDistr_->set_lower(-scale);
            pDistr_->set_upper(scale);
        }
    };

    ///////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////  *** NormalProposal_ : class ***  /////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////

    class NormalProposal_ final : public AdaptableProposal
    {
    public:
        NormalProposal_(double scale, const Distribution *prior)
            : AdaptableProposal("normal", scale, prior->is_discrete())
        {
            set_scale_(is_adaptive() ? 1. : std::fabs(scale));
        }
        MAKE_PROPOSAL_CLONE_METHOD(NormalProposal_)

        double propose(rnd_gen &prng, double xo) const override final { return xo + pDistr_->rand(prng); }
        double log_ratio(double, double) const override final { return 0.; }

    protected:
        double get_default_scale_(const Distribution *prior) override final
        {
            const auto [q25, q75] = get_iqr_from_prior_(prior);
            if (is_iqr_finite_(q25, q75))
                return (q75 - q25) / 1.349;
            return 1.;
        }
        double get_scale_() const override final { return pDistr_->arg(1); }
        /// Adapt standard deviation.
        void set_scale_(double scale) override final { pDistr_->set_arg(1, scale); }
    };

    ///////////////////////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////  *** LognormalProposal_ : class ***  ///////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////

    class LognormalProposal_ final : public AdaptableProposal
    {
    public:
        LognormalProposal_(double scale, const Distribution *prior)
            : AdaptableProposal("lognormal", scale, prior->is_discrete())
        {
            // assert that value range is positive
            if (prior->lower().is_empty() || prior->lower().value() < 0)
                throw std::invalid_argument(
                    "Provided proposal 'lognormal' for non-positive");

            // associate location parameter from prior
            if (prior->loc().is_callable())
                pAttr_loc_ = prior->loc().pointer();
            mode_ = pAttr_loc_ ? Location::Active : Location::None;

            set_scale_(is_adaptive() ? 1. : std::fabs(scale));
        }
        MAKE_PROPOSAL_CLONE_METHOD(LognormalProposal_)

        double propose(rnd_gen &prng, double xo) const override final
        {
            switch (mode_)
            {
            case Location::None:
                return xo * pDistr_->rand(prng);
            default: // Location::Active:
                const double loc = pAttr_loc_->value();
                return loc + (xo - loc) * pDistr_->rand(prng);
            }
        }
        double log_ratio(double xo, double xn) const override final
        {
            switch (mode_)
            {
            case Location::None:
                return std::log(xo / xn);
            default: // Location::Active:
                const double loc = pAttr_loc_->value();
                return std::log((xo - loc) / (xn - loc));
            }
        }

    protected:
        enum class Location : std::uint8_t
        {
            None,
            Active
        } mode_ = Location::None;
        const Attribute *pAttr_loc_ = nullptr;

        double get_default_scale_(const Distribution *prior) override final
        {
            const double loc = pAttr_loc_ ? pAttr_loc_->value() : 0.;
            const auto [q25, q75] = get_iqr_from_prior_(prior);
            if (is_iqr_finite_(q25 - loc, q75 - loc))
                return std::log((q75 - loc) / (q25 - loc)) / 1.349;
            return 1.;
        }
        double get_scale_() const override final { return pDistr_->arg(1); }
        /// Adapt log of standard deviation.
        void set_scale_(double scale) override final { pDistr_->set_arg(1, scale); }
    };

    ///////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////  *** Variable : class ***  /////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////

    Variable::Variable() = default;
    Variable::~Variable() = default;

    Variable::Variable(const Variable &other)
        : Attribute(other),
          valid_mode_(other.valid_mode_),
          ratio_mode_(other.ratio_mode_),
          old_(other.old_),
          pDistr_(other.pDistr_ ? other.pDistr_->clone_distrib() : nullptr),
          pProp_(other.pProp_ ? other.pProp_->clone() : nullptr),
          accept_rate_(other.accept_rate_)
    {
        if (pProp_)
            pProp_->rebind_distrib(pDistr_.get());
    }
    Variable::Variable(Variable &&) = default;

    Variable &Variable::operator=(Variable other)
    {
        swap(other);
        return *this;
    }

    void Variable::swap(Variable &other) noexcept
    {
        Attribute::swap(other);
        std::swap(pDistr_, other.pDistr_);
        std::swap(pProp_, other.pProp_);
    }

    void Variable::init_func_attr(std::string_view key, const Attribute &attr)
    {
        pFunc_->init_attr(key, attr);
    }
    void Variable::init_func_reshaping(bool discr)
    {
        if (discr)
            pFunc_->init_discretization();
        pFunc_->init_truncation();

        // assert no missing arguments
        pFunc_->assert_is_callable_explicitly();
    }

    void Variable::init_distr(std::string_view opt)
    {
        if (na::is(opt))
            return;
        pDistr_ = compose_distr_registries(opt, "distribution", {"continuos", "discrete"});
    }

    void Variable::init_distr_attr(std::string_view key, const Attribute &attr)
    {
        pDistr_->init_attr(key, attr);
    }

    void Variable::init_distr_reshaping(bool is_discr)
    {
        if (is_discr)
            pDistr_->init_discretization();
        pDistr_->init_truncation();

        pDistr_->assert_is_callable_explicitly();

        // check that non-constant variable has non-zero dimensional range
        helper_assert_non_degenerate_range_();
    }

    void Variable::helper_assert_non_degenerate_range_()
    {
        const Attribute &lower = pDistr_->lower();
        const Attribute &upper = pDistr_->upper();

        // skip if one of the range extremes is not provided
        if (!lower.is_callable() || !upper.is_callable())
            return;

        // skip if range is not degenerate
        const double sx = lower.value();
        const double dx = upper.value();
        if (sx < dx)
            return;

        // case empty interval, lower and upper are "open"
        if (sx > dx || !pDistr_->is_valid(sx) || !pDistr_->is_valid(dx))
            throw std::invalid_argument(
                "Definition produced an inconsistent range " + pDistr_->str_range() +
                " for " + err_name() + ".");

        if (!lower.is_constant() || !upper.is_constant())
            return;

        if (!is_user_const_)
        {
            const std::string sms = "Definitions of 'min' as " + lower.err_name() +
                                    " and 'max' as " + upper.err_name() +
                                    " produced a degenerate range ";
            // parameter
            if (err_type_ == "parameter")
                throw std::invalid_argument(sms + " for " + err_name() + ".");

            // augmented
            msg::warning("Collapsed (at least) one " + err_name() + " to constant.\n" + sms + ".");
            init_numeric(sx);
            helper_init_null_proposal_();
        }
    }

    void Variable::helper_init_null_proposal_()
    {
        // assign the pointer to the proposal
        pProp_ = std::make_unique<NullProposal_>();

        // setup sample proposal functional methods
        ratio_mode_ = Ratio::None;
        valid_mode_ = Validity::Always;
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////  Variable : Setter (Proposal)  //////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////

    void Variable::init_proposal(const std::string &proposal, double scale)
    {
        if (mode_ != Mode::Empty && mode_ != Mode::Numeric)
            return;

        if (mode_ == Mode::Empty)
            init_numeric(na::DOUBLE);
        is_user_const_ = false;

        try
        {
            if (proposal == "sample")
                helper_init_sample_proposal_();
            else
                helper_init_monte_carlo_proposal_(proposal, scale);
        }
        catch (const std::exception &e)
        {
            std::string sm = e.what();
            if (!sm.empty() && sm.back() == '.')
                sm.pop_back();
            throw std::invalid_argument(sm + " for " + err_name() + ".");
        }
    }

    void Variable::helper_init_sample_proposal_()
    {
        pProp_ = std::make_unique<SampleProposal_>(pDistr_.get());
        ratio_mode_ = Ratio::None;
        valid_mode_ = Validity::Always;
    }

    void Variable::helper_init_monte_carlo_proposal_(const std::string &proposal, double scale)
    {
        const Distribution *prior = pDistr_.get();

        const std::vector<std::string> options = {"uniform", "normal", "lognormal", "sample"};
        if (proposal == options[0])
            pProp_ = std::make_unique<UniformProposal_>(scale, prior);
        else if (proposal == options[1])
            pProp_ = std::make_unique<NormalProposal_>(scale, prior);
        else if (proposal == options[2])
            pProp_ = std::make_unique<LognormalProposal_>(scale, prior);
        else
            throw std::invalid_argument(
                "Unknown proposal '" + proposal + "'." + str::find_closest(proposal, options));

        ratio_mode_ = Ratio::Compute;
        valid_mode_ = Validity::Check;
    }

    void Variable::restart(rnd_gen &prng)
    {
        if (is_constant() || !has_model())
            return;

        set(pDistr_->rand(prng));
        pProp_->reset_scale(pDistr_.get());
    }

    void Variable::optimize()
    {
        Attribute::optimize();
        if (pDistr_)
            pDistr_->optimize();
    }

    double Variable::operator()(double x) const
    {
        switch (mode_)
        {
        case Mode::Function:
            return pFunc_->evaluate(x);
        default: // cases Mode::Empty, Mode::Numeric, Mode::Pointer
            throw std::runtime_error(
                "Called undefined evaluate operator for non-function " + err_name() + ".");
        }
    }

    void Variable::clear() noexcept
    {
        mode_ = Mode::Empty;
        value_ = na::DOUBLE;
        pAttr_ = nullptr;
        pFunc_ = nullptr;
        pDistr_ = nullptr;
        pProp_ = nullptr;
    }

    bool Variable::is_discrete() const noexcept { return has_model() ? pDistr_->is_discrete()
                                                                     : Attribute::is_discrete(); }
    bool Variable::has_model() const noexcept { return pDistr_ != nullptr; }
    bool Variable::is_adaptive() const noexcept { return pProp_ && pProp_->is_adaptive(); }

    void Variable::propose(rnd_gen &prng)
    {
        set(pProp_->propose(prng, value_));
        accept_rate_ += adapt_smoothing_ * (1. - accept_rate_); // EMA for accept = TRUE
    }
    bool Variable::is_valid_move() const
    {
        switch (valid_mode_)
        {
        case Validity::Check:
            return pDistr_->is_valid(value_);
        default: // case Validity::Always:
            return true;
        }
    }
    double Variable::get_distr_log_pdf() const { return pDistr_->log_pdf(value_); }
    double Variable::get_distr_log_ratio() const { return pDistr_->log_ratio(old_, value_); }
    double Variable::get_propose_log_ratio() const { return pProp_->log_ratio(old_, value_); }

    void Variable::adapt_proposal(double ratio) { pProp_->adapt(accept_rate_, ratio); }

#ifndef NDEBUG
    std::string Variable::str_summary() const
    {
        return Attribute::str_summary() +
               (has_model() ? " | ~" + pDistr_->str_summary(false)
                            : "");
    }
#endif

    /////////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////  *** Parameter : struct ***  ////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////

    void Parameter::build(const str_table &input)
    {
        khold_.clear_optional();

        if (input.empty())
            return;

        // check entries
        const std::vector<std::string> file_headers = {"name", "value", "const", "discr",
                                                       "func", "prior", "min", "max",
                                                       "proposal", "proposal_scale"};
        const std::vector<std::string> headers = read::extract_keys(input);
        str::assert_no_missing_entries(headers, file_headers);
        str::assert_no_unused_entries(headers, file_headers);

        // build dependency graph and validate function parameters
        bool has_func_params = false;
        dict<size_t> name_to_row;
        for (size_t k = 0; k < input.at("name").size(); ++k)
        {
            const std::string &pstr = input.at("name")[k];
            name_to_row[pstr] = k;

            // validate: function parameters can only depend on other parameters
            if (!na::is(input.at("func")[k]))
            {
                has_func_params = true;
                parse::Definition parser(input.at("func")[k]);
                if (parser.is_augmented_depen())
                    throw std::invalid_argument(
                        "Parameter function '" + pstr +
                        "' cannot depend on augmented variables.");
                if (parser.is_data_wide_depen())
                    throw std::invalid_argument(
                        "Parameter function '" + pstr +
                        "' cannot depend on data variables.");
            }
        }

        // topological sort only when function parameters exist (preserves input order otherwise)
        std::vector<std::string> ordered;
        if (has_func_params)
        {
            // use input order as iteration base for stable topological sort
            str_table dep_graph;
            for (size_t k = 0; k < input.at("name").size(); ++k)
            {
                const std::string &pstr = input.at("name")[k];
                std::vector<std::string> deps;
                for (const auto col : {"func", "min", "max"})
                    deps += parse::Definition(input.at(col)[k]).list_parameter_depen();
                dep_graph[pstr] = str::unique(deps);
            }
            ordered = dependency::order_dependencies(dep_graph);
        }
        else
        {
            ordered = input.at("name");
        }

        init_param_vec_();
        clear_lists_();
        // load parameters in dependency order
        for (const std::string &pstr : ordered)
        {
            size_t k = name_to_row.at(pstr);

            if (!khold_.count(pstr))
            {
                khold_.add(pstr);
                param_vec_.emplace_back();
                param_vec_.back().name(pstr);
                param_vec_.back().set_err_type("parameter");
            }

            if (!na::is(input.at("func")[k]))
                helper_build_function_(input, k);
            else if (str::tobool(input.at("const")[k]))
                helper_build_constant_(input, k);
            else
                helper_build_to_sample_(input, k);
        }

        // initialize all parameters to standard (i.e. no cache dependencies)
        list_sample_standard = list_sample;
        for (size_t p : list_sample)
            if (at(p).is_adaptive())
                list_adaptive.push_back(p);
    }

    void Parameter::build_for_predic(const dbl_table &input)
    {
        const std::vector<std::string> headers = read::extract_keys(input);
#ifndef NDEBUG
        if (headers.size() != str::unique(headers).size())
            throw std::invalid_argument(
                "Found multiple definitions for the same parameters.");
#endif

        clear_lists_();
        simulation_values.clear();
        // update list of constant parameters to all parameters
        list_const = khold_.from_string(khold_.get_all());

        // reset definition dictionary to all constants NAs.
        for (const std::string &pstr : khold_.get_all())
            if (std::find(headers.begin(), headers.end(), pstr) == headers.end())
                simulation_values.emplace_back();
            else
                simulation_values.emplace_back(input.at(pstr));
    }

    void Parameter::init_param_vec_()
    {
        size_t n = khold_.size();
        param_vec_.clear();
        param_vec_.resize(n);
        for (size_t p = 0; p < n; ++p)
        {
            param_vec_[p].name(khold_.as_string(p));
            param_vec_[p].set_err_type("parameter");
        }
    }

    void Parameter::clear_lists_()
    {
        list_const.clear();
        list_func.clear();
        list_sample.clear();
        list_sample_standard.clear();
        list_adaptive.clear();
        list_sample_cached.clear();

        simulation_values.clear();
    }

    void Parameter::helper_build_constant_(const str_table &input, size_t k)
    {
        const std::string &pstr = input.at("name")[k]; // parameter name

#ifndef NDEBUG
        if (na::is(input.at("value")[k]))
            throw std::invalid_argument(
                "Constant parameter '" + pstr + "' had undefined value.");
#endif

        // assign constant parameter to the set of variables
        size_t p = khold_.from_string(pstr);
        param_vec_[p].transfer(read_numeric(input.at("value")[k]));
        // check if discrete as requested
        if (user_discrete_set_.count(p) && !param_vec_[p].is_discrete())
            throw std::invalid_argument(
                "Definition of " + param_vec_[p].err_name() +
                " did not return a discrete object as requested by model.");

        // update the lists of constant parameters
        list_const.push_back(p);
    }

    void Parameter::helper_build_to_sample_(const str_table &input, size_t k)
    {
        const std::string &pstr = input.at("name")[k]; // parameter name

#ifndef NDEBUG
        if (na::is(input.at("value")[k]) && na::is(input.at("prior")[k]))
            throw std::invalid_argument(
                "Non-constant parameter '" + pstr + "' had undefined prior.");
#endif

        // create a new parameter Variable
        size_t p = khold_.from_string(pstr);
        Variable &xparam = at(p);
        xparam.transfer(read_numeric(input.at("value")[k]));

        // assign a prior distribution to the parameter
        parse::Definition parser(input.at("prior")[k]);
        xparam.init_distr(parser.root.func);

        // load prior arguments
        for (const auto &child : parser.root.children)
            xparam.init_distr_attr(child.key, read_numeric(child));
        // load parameter range
        for (const auto key : {"min", "max"})
            xparam.init_distr_attr(key, read_numeric(input.at(key)[k]));
        xparam.init_distr_reshaping(str::tobool(input.at("discr")[k]));

        // assign proposal
        xparam.init_proposal(input.at("proposal")[k],
                             str::tod(input.at("proposal_scale")[k]));

        // check if discrete as requested
        if (user_discrete_set_.count(p) && !xparam.is_discrete())
            throw std::invalid_argument(
                "Definition of " + xparam.err_name() +
                " did not return a discrete object as requested by model.");

        // update the lists of parameters to sample
        list_sample.push_back(p);
    }

    void Parameter::helper_build_function_(const str_table &input, size_t k)
    {
        const std::string &pstr = input.at("name")[k];
        size_t p = khold_.from_string(pstr);
        Variable &xparam = at(p);

        parse::Definition parser(input.at("func")[k]);

        // pointer mode (terminal = alias for another parameter)
        if (parser.is_terminal())
        {
            xparam.transfer(read_attribute(parser.root));
            list_const.push_back(p);
            list_func.push_back(p);
            return;
        }

        // function mode
        xparam.init_function(compose_func_registries(
            parser.root.func, "function",
            {"custom", "transcendental", "arithmetic"}));

        // bind function arguments
        for (const auto &child : parser.root.children)
            xparam.init_func_attr(child.key, read_attribute(child));

        // range
        for (const auto key : {"min", "max"})
            xparam.init_func_attr(key, read_attribute(input.at(key)[k]));

        // discretization
        xparam.init_func_reshaping(str::tobool(input.at("discr")[k]));

        // check if discrete as requested
        if (user_discrete_set_.count(p) && !xparam.is_discrete())
            throw std::invalid_argument(
                "Definition of " + xparam.err_name() +
                " did not return a discrete object as requested by model.");

        // classify: deterministic, not sampled
        list_const.push_back(p);
        list_func.push_back(p);
    }

    void Parameter::make_cache_dependent(size_t p)
    {
        // remove a parameter from standard if it was defined as such
        if (std::find(list_sample_standard.begin(), list_sample_standard.end(), p) != list_sample_standard.end())
        {
            const std::vector<size_t> vec({p});
            list_sample_standard = list_sample_standard - vec;
        }

        // add the parameter as cache if not already present
        if (std::find(list_sample_cached.begin(), list_sample_cached.end(), p) == list_sample_cached.end())
            list_sample_cached.push_back(p);
    }

    void Parameter::restart(rnd_gen &prng)
    {
        for (size_t p : list_sample)
            at(p).restart(prng);
    }

    void Parameter::optimize()
    {
        for (size_t p : list_const)
            at(p).optimize();
        for (size_t p : list_sample)
            at(p).optimize();
    }

    void Parameter::load_sample_values(size_t r)
    {
        for (size_t p = 0; p < khold_.size(); ++p)
        {
            const std::vector<double> &vec = simulation_values[p];
            if (!vec.empty())
                at(p).set(vec[r]);
        }
    }

    void Parameter::make_cache_dependent(const parse::Entry &node)
    {
        // cache dependecy if an augmented function depends on a parameter to be inferred
        if (parse::is_parameter(node.type))
        {
            size_t p = khold_.from_string(node.value);
            if (!at(p).is_constant())
                make_cache_dependent(p);
        }
    }

    const Attribute &Parameter::read_attribute(
        const parse::Entry &node,
        const Individual *member_ptr) const
    {
        storage_ = Attribute(); // reset auxiliary storage to default empty

        if (!node.is_terminal())
        {
            helper_store_read_function_(node, member_ptr);
        }
        else if (parse::is_numeric(node.type))
        {
            storage_.init_numeric(str::tod(node.value, "argument with numeric type"));
        }
        else if (parse::is_parameter(node.type))
        {
            return at(khold_.from_string(node.value));
        }
        else if (!node.is_empty())
        {
            if (!member_ptr)
            {
                std::string sms = "Individual specific definitions as " + node.err_str() +
                                  " are not supported here.";
                // add tip if augmented (which does not require '@')
                if (parse::is_augmented(node.type) && node.is_terminal())
                    sms += "\nDid you forget rlang injection operator '!!'?";
                throw std::invalid_argument(sms);
            }

            if (parse::is_data_wide(node.type))
            {
                size_t d = pDataWideKeys_->from_string(node.value);
                storage_.init_numeric(member_ptr->data(d));
                storage_.name(node.value);
                storage_.set_err_type("data");
            }
            else if (parse::is_augmented(node.type))
            {
                return member_ptr->get_augm_variable(pAugmKeys_->from_string(node.value));
            }
            else if (parse::is_group_rate(node.type))
            {
                return member_ptr->get_group_rate_attribute(pGroupRateKeys_->from_string(node.value));
            }
            else
                throw std::invalid_argument(
                    "Unrecognized source type '" + node.type +
                    "' found while parsing " + node.err_str() + ".");
        }

        return storage_;
    }

    void Parameter::helper_store_read_function_(
        const parse::Entry &node,
        const Individual *member_ptr) const
    {
        Attribute parent;
        parent.init_function(compose_func_registries(
            node.func, "function",
            {"custom", "transcendental", "arithmetic",
             "continuos", "discrete", "infecprofile"}));

        for (const auto &child : node.children)
            parent.init_func_attr(child.key, read_attribute(child, member_ptr));

        storage_ = parent;
    }

    const Attribute &Parameter::read_numeric(
        const parse::Entry &node,
        const Individual *member_ptr) const
    {
        storage_ = read_attribute(node, member_ptr);
        if (!storage_.is_empty() && !storage_.is_constant())
            throw std::invalid_argument(
                "Definition " + node.err_str() +
                " could not be interpreted as a numeric constant.");
        return storage_;
    }
}
