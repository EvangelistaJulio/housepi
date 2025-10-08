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
 * distribution.cpp
 *
 */

#include "tools/operators.h"
#include "tools/discr.h"
#include "tools/params.h"
#include "tools/string.h"
#include "distribution.h"
#include "numerical_recipes.h"

/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////  Distribution - rand  /////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

// UPDATE: random generator classes should be redesigned to avoid update_args_...
double Distribution::rand(rnd_gen &prng)
{
    update_args_();
    return (this->*rand_)(prng);
}

/// Rejection sampling.
/// UPDATE: not efficient at all, need to use quantiles probably.
double Distribution::rejec_sample_(rnd_gen &prng, double (Distribution::*func)(rnd_gen &))
{
    double xnew = (this->*func)(prng);
    for (unsigned n = 0; !is_valid(xnew) && (n < 100); ++n)
        xnew = (this->*func)(prng);
    if (!is_valid(xnew))
        throw std::runtime_error(
            "Exceeded attempts (N > " + std::to_string(100) +
            ") in random draw of " + str_summary(false) + ".");
    return xnew;
}

/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////  Distribution - reshaping  ///////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

/// Select functional for discretization of random sampling.
void Distribution::init_discr_helper_(unsigned short discr_type)
{
    Function::init_discr_helper_(discr_type);

    switch (discr_type)
    {
    case hsp::discr::CENTER: // centered discretization
        rand_ = &Distribution::do_discr_center_rand_;
        break;
    case hsp::discr::LEFT: // left discretization
        rand_ = &Distribution::do_discr_left_rand_;
        break;
    default: // discrete defined or no discretization
        rand_ = &Distribution::raw_rand_;
    }
}

/// Select functional method in case of truncation.
void Distribution::init_trunc_helper_()
{
    unsigned short discr_type = which_discrete_case_();

    switch (which_valid_case_())
    {
    case 1: // Lower Range Truncation
        switch (discr_type)
        {
        case hsp::discr::DEFINED: // discrete defined
            trunc_ = &Distribution::do_trunc_min_discr_defined_;
            logtrunc_ = &Distribution::do_logtrunc_min_discr_defined_;
            shift_ = &Distribution::do_shift_min_discr_defined_;
            rand_ = &Distribution::do_trunc_rand_;
            break;
        case hsp::discr::CENTER: // centered discretization
            trunc_ = &Distribution::do_trunc_min_discr_center_;
            logtrunc_ = &Distribution::do_logtrunc_min_discr_center_;
            shift_ = &Distribution::do_shift_min_discr_center_;
            rand_ = &Distribution::do_trunc_rand_discr_center_;
            break;
        case hsp::discr::LEFT: // left discretization
            trunc_ = &Distribution::do_trunc_min_discr_left_;
            logtrunc_ = &Distribution::do_logtrunc_min_discr_left_;
            shift_ = &Distribution::do_shift_min_discr_left_;
            rand_ = &Distribution::do_trunc_rand_discr_left_;
            break;
        default: // no discretization
            trunc_ = &Distribution::do_trunc_min_;
            logtrunc_ = &Distribution::do_logtrunc_min_;
            shift_ = &Distribution::do_shift_min_;
            rand_ = &Distribution::do_trunc_rand_;
        }
        break;

    case 2: // Upper Range Truncation
        shift_ = &Distribution::do_shift_dummy_;
        switch (discr_type)
        {
        case hsp::discr::DEFINED: // discrete defined
            trunc_ = &Distribution::do_trunc_max_discr_defined_;
            logtrunc_ = &Distribution::do_logtrunc_max_discr_defined_;
            rand_ = &Distribution::do_trunc_rand_;
            break;
        case hsp::discr::CENTER: // centered discretization
            trunc_ = &Distribution::do_trunc_max_discr_center_;
            logtrunc_ = &Distribution::do_logtrunc_max_discr_center_;
            rand_ = &Distribution::do_trunc_rand_discr_center_;
            break;
        case hsp::discr::LEFT: // left discretization
            trunc_ = &Distribution::do_trunc_max_discr_left_;
            logtrunc_ = &Distribution::do_logtrunc_max_discr_left_;
            rand_ = &Distribution::do_trunc_rand_discr_left_;
            break;
        default: // no discretization
            trunc_ = &Distribution::do_trunc_max_;
            logtrunc_ = &Distribution::do_logtrunc_max_;
            rand_ = &Distribution::do_trunc_rand_;
        }
        break;

    case 3: // Both Sides Truncation
        switch (discr_type)
        {
        case hsp::discr::DEFINED: // discrete defined
            trunc_ = &Distribution::do_trunc_both_sides_discr_defined_;
            logtrunc_ = &Distribution::do_logtrunc_both_sides_discr_defined_;
            shift_ = &Distribution::do_shift_min_discr_defined_;
            rand_ = &Distribution::do_trunc_rand_;
            break;
        case hsp::discr::CENTER: // centered discretization
            trunc_ = &Distribution::do_trunc_both_sides_discr_center_;
            logtrunc_ = &Distribution::do_logtrunc_both_sides_discr_center_;
            shift_ = &Distribution::do_shift_min_discr_center_;
            rand_ = &Distribution::do_trunc_rand_discr_center_;
            break;
        case hsp::discr::LEFT: // left discretization
            trunc_ = &Distribution::do_trunc_both_sides_discr_left_;
            logtrunc_ = &Distribution::do_logtrunc_both_sides_discr_left_;
            shift_ = &Distribution::do_shift_min_discr_left_;
            rand_ = &Distribution::do_trunc_rand_discr_left_;
            break;
        default: // no discretization
            trunc_ = &Distribution::do_trunc_both_sides_;
            logtrunc_ = &Distribution::do_logtrunc_both_sides_;
            shift_ = &Distribution::do_shift_min_;
            rand_ = &Distribution::do_trunc_rand_;
        }
        break;

    default: // No Truncation (case 0)
        trunc_ = &Distribution::do_trunc_dummy_;
        logtrunc_ = &Distribution::do_logtrunc_dummy_;
        shift_ = &Distribution::do_shift_dummy_;
        break;
    }
}

/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////  Distribution Registry : getter  //////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

/// Override get_options method of generic function with distribution.
/// note: this differ in the choice of factory_map_ member.
str_vec DistributionRegistry::get_options() const
{
    str_vec keys(factory_map_.size());
    for (const auto &pair : factory_map_)
        keys.push_back(pair.first);
    return keys;
}

std::shared_ptr<Distribution> DistributionRegistry::get_distr(const std::string &opt) const
{
    if (factory_map_ % opt)
        return factory_map_.find(opt)->second();
    return nullptr;
}

std::shared_ptr<Function> DistributionRegistry::get(const std::string &opt, bool is_primitive) const
{
    auto output = get_distr(opt);
    if (is_primitive)
        output->set_eval_as_distr_primitive();
    else
        output->set_eval_as_distr_call();
    return output;
}

/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////  Distribution Registry : setter  //////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

std::shared_ptr<Distribution> compose_distr_registries(const std::string &opt,
                                                       const std::string &name,
                                                       const std::vector<const DistributionRegistry *> &registries)
{
    str_vec options = {};

    if (registries.empty())
        throw std::invalid_argument("No registries provided for \"" + name + "\".");
    else if (hsp::na::is(opt))
        return std::make_shared<Distribution>();
    else
        for (auto &r : registries)
            if (r->get_options() % opt)
                return r->get_distr(opt);
            else
                options += r->get_options();

    throw std::invalid_argument(
        "Unknown " + name + " \"" + opt + "\"." +
        hsp::str::find_closest(opt, options));
}

/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////  MACROS  ///////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

// Macro for left discretization when min = 0 (close).
#define INIT_DISCR_LEFT_MIN_ZERO_CLOSE()      \
    void init_discretization() override final \
    {                                         \
        init_discr_helper_(hsp::discr::LEFT);   \
    }                                         \
    void init_truncation() override final     \
    {                                         \
        if (init_range_min_(0., false))       \
            init_trunc_helper_();             \
    }

#define GET_ARGS_1 arg(0)
#define GET_ARGS_2 arg(0), arg(1)

#define MAKE_RAW_NUMR_IMPORT_FOR_DISTR(distr_, n_args)              \
    double raw_pdf_(double x) const override final                  \
    {                                                               \
        return numr::distr_::pdf(x, GET_ARGS_##n_args);             \
    }                                                               \
    double raw_cdf_(double x) const override final                  \
    {                                                               \
        return numr::distr_::cdf(x, GET_ARGS_##n_args);             \
    }                                                               \
    double raw_logpdf_(double x) const override final               \
    {                                                               \
        return numr::distr_::logpdf(x, GET_ARGS_##n_args);          \
    }                                                               \
    double raw_logratio_(double xo, double xn) const override final \
    {                                                               \
        return numr::distr_::logratio(xo, xn, GET_ARGS_##n_args);   \
    }

/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////  Uniform Distribution : class  //////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

/// Uniform Distribution derived class.
class Uniform_Distribution : public Distribution
{
public:
    Uniform_Distribution()
        : Distribution("uniform", {hsp::params::MIN, hsp::params::MAX}, {0., 1.},
                       hsp::na::DOUBLE, false, hsp::na::DOUBLE, false),
          unif_gen(arg(0), arg(1)) { init_uniform_range_(); }

    void init_discretization() override final
    {
        init_discr_helper_(hsp::discr::CENTER);
    }
    void init_truncation() override final
    {
        init_uniform_range_();
        init_range_helper_(false, false);
        // NOTE: no truncation for uniform distribution
    }

protected:
    std::uniform_real_distribution<double> unif_gen;

    double raw_pdf_(double) const override final { return 1.; }
    double raw_cdf_(double x) const override final { return x; }
    double raw_logpdf_(double) const override final { return 0.; }
    double raw_logratio_(double, double) const override final { return 0.; }
    double raw_rand_(rnd_gen &prng) override final
    {
        const double r = unif_gen(prng);
        return arg(0) * (1 - r) + arg(1) * r;
    }

private:
    /// Ensure uniform range (min, max) match argument (min, max).
    void init_uniform_range_()
    {
        // impose range minimum equal to argument minimum
        const Attribute &min_arg = arguments_[0];
        if (min_arg.is_function())
            range_min_ = min_arg.get_function();
        else
            range_min_ = min_arg.get_ptr();

        // impose range maximum equal to argument maximum
        const Attribute &max_arg = arguments_[1];
        if (max_arg.is_function())
            range_max_ = max_arg.get_function();
        else
            range_max_ = max_arg.get_ptr();
    }
};

/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////  Normal Distribution : class  //////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

/// Normal Distribution derived class.
class Normal_Distribution : public Distribution
{
public:
    Normal_Distribution()
        : Distribution("normal", {"mean", "sd"}, {0., 1.},
                       hsp::na::DOUBLE, false, hsp::na::DOUBLE, false),
          norm_gen(0., 1.) {}

    void init_discretization() override final
    {
        switch (which_valid_case_())
        {
        case 1: // lower range truncation
            init_discr_helper_(hsp::discr::LEFT);
            break;
        case 3: // both sides truncation
            init_discr_helper_(hsp::discr::LEFT);
            break;
        default:
            init_discr_helper_(hsp::discr::CENTER);
        }
    }
    void init_truncation() override final
    {
        init_range_helper_(false, false);
        init_trunc_helper_();
    }

protected:
    std::normal_distribution<double> norm_gen;

    MAKE_RAW_NUMR_IMPORT_FOR_DISTR(normal, 2)
    double raw_rand_(rnd_gen &prng) override final { return arg(0) + arg(1) * norm_gen(prng); }
};

/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////  Lognormal Distribution : class  /////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

/// Lognormal Distribution derived class.
class Lognormal_Distribution : public Distribution
{
public:
    Lognormal_Distribution()
        : Distribution("lognormal", {"meanlog", "sdlog"}, {0., 1.},
                       0., true, hsp::na::DOUBLE, false),
          norm_gen(0., 1.) {}

    INIT_DISCR_LEFT_MIN_ZERO_CLOSE()

protected:
    std::normal_distribution<double> norm_gen;

    MAKE_RAW_NUMR_IMPORT_FOR_DISTR(lognormal, 2)
    double raw_rand_(rnd_gen &prng) override final { return std::exp(arg(0) + arg(1) * norm_gen(prng)); }
};

/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////  Beta Distribution : class  //////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

/// Beta Distribution derived class.
class Beta_Distribution : public Distribution
{
public:
    Beta_Distribution()
        : Distribution("beta", {"shape1", "shape2"}, {1., 1.},
                       0., true, 1., true),
          gamma_1(arg(0), 1.), gamma_2(arg(1), 1.) {}

    void init_discretization() override final
    {
        throw std::invalid_argument(
            "Invalid option discrete for " + err_name() + ".");
    }
    void init_truncation() override final
    {
        if (get_min() < 0)
            throw std::invalid_argument(
                "Invalid min value " + hsp::str::format(range_min_.get()) + " for " + err_name() + ".");
        else if (get_max() > 1)
            throw std::invalid_argument(
                "Invalid max value " + hsp::str::format(range_max_.get()) + " for " + err_name() + ".");
        else if ((get_min() > 0) || (get_max() < 1))
        {
            init_range_helper_((get_min() == 0), (get_max() == 1));
            init_trunc_helper_();
        }
    }

protected:
    std::gamma_distribution<double> gamma_1;
    std::gamma_distribution<double> gamma_2;

    // UPDATE
    void update_args_() override final
    {
        gamma_1.param(std::gamma_distribution<double>::param_type(arg(0), 1.));
        gamma_2.param(std::gamma_distribution<double>::param_type(arg(1), 1.));
    }

    MAKE_RAW_NUMR_IMPORT_FOR_DISTR(beta, 2)
    double raw_rand_(rnd_gen &prng) override final
    {
        double x1 = gamma_1(prng);
        double x2 = gamma_2(prng);
        return x1 / (x1 + x2);
    }
};

/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////  Gamma Distribution : class  /////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

/// Gamma Distribution derived class.
class Gamma_Distribution : public Distribution
{
public:
    Gamma_Distribution()
        : Distribution("gamma", {"shape", "scale"}, {1., 1.},
                       0., false, hsp::na::DOUBLE, false),
          distr_gen(arg(0), arg(1)) {}

    INIT_DISCR_LEFT_MIN_ZERO_CLOSE()

protected:
    std::gamma_distribution<double> distr_gen;

    void update_args_() override final { distr_gen.param(std::gamma_distribution<double>::param_type(arg(0), arg(1))); } // UPDATE

    MAKE_RAW_NUMR_IMPORT_FOR_DISTR(gamma, 2)
    double raw_rand_(rnd_gen &prng) override final { return distr_gen(prng); }
};

/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////  Gamma (Mean-Varance) Distribution : class  /////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

#define GAMMA_MEANVAR_CONVERTER()     \
    const double m = arg(0);          \
    const double v = arg(1);          \
    const double shape = (m * m) / v; \
    const double scale = v / m;

/// Gamma (Mean-Varance) Distribution derived class.
class GammaMeanVar_Distribution : public Distribution
{
public:
    GammaMeanVar_Distribution()
        : Distribution("gamma_meanvar", {"mean", "var"}, {1., 1.},
                       0., false, hsp::na::DOUBLE, false),
          distr_gen(arg(0), arg(1)) {}

    INIT_DISCR_LEFT_MIN_ZERO_CLOSE()

protected:
    std::gamma_distribution<double> distr_gen;

    void update_args_() override final
    {
        GAMMA_MEANVAR_CONVERTER()
        distr_gen.param(std::gamma_distribution<double>::param_type(shape, scale));
    }

    double raw_pdf_(double x) const override final
    {
        GAMMA_MEANVAR_CONVERTER()
        return numr::gamma::pdf(x, shape, scale);
    }
    double raw_cdf_(double x) const override final
    {
        GAMMA_MEANVAR_CONVERTER()
        return numr::gamma::cdf(x, shape, scale);
    }
    double raw_logpdf_(double x) const override final
    {
        GAMMA_MEANVAR_CONVERTER()
        return numr::gamma::logpdf(x, shape, scale);
    }
    double raw_logratio_(double xo, double xn) const override final
    {
        GAMMA_MEANVAR_CONVERTER()
        return numr::gamma::logratio(xo, xn, shape, scale);
    }
    double raw_rand_(rnd_gen &prng) override final
    {
        return distr_gen(prng);
    }
};

/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////  Rayleigh Distribution : class  ////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

/// Rayleigh Distribution derived class.
class Rayleigh_Distribution : public Distribution
{
public:
    Rayleigh_Distribution()
        : Distribution("rayleigh", {"scale"}, {1.},
                       0., false, hsp::na::DOUBLE, false),
          unif_gen(0., 1.) {}

    INIT_DISCR_LEFT_MIN_ZERO_CLOSE()

protected:
    std::uniform_real_distribution<double> unif_gen;

    MAKE_RAW_NUMR_IMPORT_FOR_DISTR(rayleigh, 1)
    double raw_rand_(rnd_gen &prng) override final { return arg(0) * sqrt(-2 * std::log(unif_gen(prng))); }
};

/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////  Weibull Distribution : class  /////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

/// Weibull Distribution derived class.
class Weibull_Distribution : public Distribution
{
public:
    Weibull_Distribution()
        : Distribution("weibull", {"shape", "scale"}, {1., 1.},
                       0., false, hsp::na::DOUBLE, false),
          unif_gen(0., 1.) {}

    INIT_DISCR_LEFT_MIN_ZERO_CLOSE()

protected:  
    std::uniform_real_distribution<double> unif_gen;

    MAKE_RAW_NUMR_IMPORT_FOR_DISTR(weibull, 2)
    double raw_rand_(rnd_gen &prng) override final { return arg(1) * std::pow(-std::log(unif_gen(prng)), 1.0 / arg(0)); }
};

/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////  Exponential  Distribution : class  //////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

/// Exponential  Distribution derived class.
class Exponential_Distribution : public Distribution
{
public:
    Exponential_Distribution()
        : Distribution("exponential", {"rate"}, {1.},
                       0., false, hsp::na::DOUBLE, false),
          unif_gen(0., 1.) {}

    INIT_DISCR_LEFT_MIN_ZERO_CLOSE()

protected:
    std::uniform_real_distribution<double> unif_gen;

    MAKE_RAW_NUMR_IMPORT_FOR_DISTR(exponential, 1)
    double raw_rand_(rnd_gen &prng) override final { return -std::log(unif_gen(prng)) / arg(0); }
};

/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////  Poisson Distribution : class  ////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

#define INIT_DISCR_DEFINED_MIN_ZERO_CLOSE()      \
    void init_discretization() override final {} \
    void init_truncation() override final        \
    {                                            \
        if (init_range_min_(0., false))          \
            init_trunc_helper_();                \
    }

/// Poisson Distribution derived class.
class Poisson_Distribution : public Distribution
{
public:
    Poisson_Distribution()
        : Distribution("poisson", {"lambda"}, {1.},
                       0., false, hsp::na::DOUBLE, false),
          distr_gen(1.) { init_discr_helper_(hsp::discr::DEFINED); }

    INIT_DISCR_DEFINED_MIN_ZERO_CLOSE()

protected:
    std::poisson_distribution<unsigned> distr_gen;

    void update_args_() override final { distr_gen.param(std::poisson_distribution<unsigned>::param_type(arg(0))); } // UPDATE

    MAKE_RAW_NUMR_IMPORT_FOR_DISTR(poisson, 1)
    double raw_rand_(rnd_gen &prng) override final { return distr_gen(prng); }
};

/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////  Binomial Distribution : class  ///////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

/// Binomial Distribution derived class.
class Binomial_Distribution : public Distribution
{
public:
    Binomial_Distribution()
        : Distribution("binomial", {"size", "prob"}, {4, 0.5},
                       0., false, hsp::na::DOUBLE, false),
          distr_gen(arg(0), arg(1))
    {
        init_discr_helper_(hsp::discr::DEFINED);
        init_attr(hsp::params::MAX, get_ptr_to_arg(0));
        init_range_helper_(false, false);
    }

    void init_discretization() override final {}
    void init_truncation() override final
    {
        if (get_min() < 0.)
            throw std::invalid_argument(
                "Invalid min value " + hsp::str::format(range_min_.get()) + " for " + err_name() + ".");
        else if ((get_min() > 0.) || (get_ptr_to_arg(0) != get_ptr_to_max()))
        {
            if (get_max() > arg(0))
                throw std::invalid_argument(
                    "Invalid max value " + hsp::str::format(range_max_.get()) + " for " + err_name() + " with size " + std::to_string(arg(0)) + ".");
            init_range_helper_(false, false);
            init_trunc_helper_();
        }
    }

protected:
    std::binomial_distribution<unsigned> distr_gen;

    void update_args_() override final { distr_gen.param(std::binomial_distribution<unsigned>::param_type(arg(0), arg(1))); } // UPDATE

    MAKE_RAW_NUMR_IMPORT_FOR_DISTR(binomial, 2)
    double raw_rand_(rnd_gen &prng) override final { return distr_gen(prng); }
};

/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////  NegativeBinomial Distribution : class  ///////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

/// NegativeBinomial Distribution derived class.
class NegativeBinomial_Distribution : public Distribution
{
public:
    NegativeBinomial_Distribution()
        : Distribution("negative_binomial", {"size", "prob"}, {4, .5},
                       0., false, hsp::na::DOUBLE, false),
          distr_gen(arg(0), arg(1)) { init_discr_helper_(hsp::discr::DEFINED); }

    INIT_DISCR_DEFINED_MIN_ZERO_CLOSE()

protected:
    std::negative_binomial_distribution<unsigned> distr_gen;

    void update_args_() override final { distr_gen.param(std::negative_binomial_distribution<unsigned>::param_type(arg(0), arg(1))); } // UPDATE

    MAKE_RAW_NUMR_IMPORT_FOR_DISTR(negbinomial, 2)
    double raw_rand_(rnd_gen &prng) override final { return distr_gen(prng); }
};

/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////  Geometric Distribution : class  ///////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

/// Geometric Distribution derived class.
class Geometric_Distribution : public Distribution
{
public:
    Geometric_Distribution()
        : Distribution("geometric", {"proba"}, {0.5},
                       0., false, hsp::na::DOUBLE, false),
          distr_gen(arg(0)) { init_discr_helper_(hsp::discr::DEFINED); }

    INIT_DISCR_DEFINED_MIN_ZERO_CLOSE()

protected:
    std::geometric_distribution<unsigned> distr_gen;

    void update_args_() override final { distr_gen.param(std::geometric_distribution<unsigned>::param_type(arg(0))); } // UPDATE

    MAKE_RAW_NUMR_IMPORT_FOR_DISTR(geometric, 1)
    double raw_rand_(rnd_gen &prng) override final { return distr_gen(prng); }
};

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////  Register Distributions  ///////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

#define IMPLEMENT_HOUSEPI_DISTRIBUTION_REGISTRY(ClassName) \
    ClassName &ClassName::inst()                           \
    {                                                      \
        static ClassName _inst;                            \
        return _inst;                                      \
    }

IMPLEMENT_HOUSEPI_DISTRIBUTION_REGISTRY(ContinuosDistrRegistry)
REGISTER_HOUSEPI_DISTRIBUTION(ContinuosDistrRegistry, Uniform_Distribution)
REGISTER_HOUSEPI_DISTRIBUTION(ContinuosDistrRegistry, Normal_Distribution)
REGISTER_HOUSEPI_DISTRIBUTION(ContinuosDistrRegistry, Lognormal_Distribution)
REGISTER_HOUSEPI_DISTRIBUTION(ContinuosDistrRegistry, Beta_Distribution)
REGISTER_HOUSEPI_DISTRIBUTION(ContinuosDistrRegistry, Gamma_Distribution)
REGISTER_HOUSEPI_DISTRIBUTION(ContinuosDistrRegistry, GammaMeanVar_Distribution)
REGISTER_HOUSEPI_DISTRIBUTION(ContinuosDistrRegistry, Rayleigh_Distribution)
REGISTER_HOUSEPI_DISTRIBUTION(ContinuosDistrRegistry, Weibull_Distribution)
REGISTER_HOUSEPI_DISTRIBUTION(ContinuosDistrRegistry, Exponential_Distribution)

IMPLEMENT_HOUSEPI_DISTRIBUTION_REGISTRY(DiscreteDistrRegistry)
REGISTER_HOUSEPI_DISTRIBUTION(DiscreteDistrRegistry, Poisson_Distribution)
REGISTER_HOUSEPI_DISTRIBUTION(DiscreteDistrRegistry, Binomial_Distribution)
REGISTER_HOUSEPI_DISTRIBUTION(DiscreteDistrRegistry, NegativeBinomial_Distribution)
REGISTER_HOUSEPI_DISTRIBUTION(DiscreteDistrRegistry, Geometric_Distribution)

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////  InfecProfile_Ferretti20 : class  ///////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/// The infectivity profile is considered equal to the distribution
/// of the Time from Onset of Symptoms to Transmission 't_ost':
/// >	t_ost = <time infection infectee> - <time onset infector>
/// Its distribution is modeled as a skew-logistic distribution (generalized logistic type I),
/// where the pre-symptomatic values are rescaled as a function of the incubation time 'incub'.
/// The arguments to provide are :
///   0) "shape";
///   1) "loc";
///   2) "scale";
///   3) "mean_incub" (mean incubation period);
///   4) "incub_period" (individual incubation period);
/// Note: symtoms can occur after infection (pre-symptomatic transmission),
/// i. e. 't_ost' can be negative.
///
/// Ref:
/// ----
/// Ferretti L, Ledda A, Wymant C, Zhao L, Ledda V, Abeler-Dorner L, et al.
/// The timing of COVID-19 transmission.
/// medRxiv; 2020. p. 2020.09.04.20188516. doi:10.1101/2020.09.04.20188516
class InfecProfile_Ferretti20 : public Distribution
{
public:
    InfecProfile_Ferretti20()
        : Distribution("ferretti20",
                       {"shape", "loc", "scale", "mean_incub", "incub_period"},
                       {5.85, -4.0, 1.85, 5.8, 5.8},
                       0., false, hsp::na::DOUBLE, false) {}

    INIT_DISCR_LEFT_MIN_ZERO_CLOSE()

protected:
    double raw_pdf_(double x) const override final
    {
        const double a = arg(0); // 'shape'
        const double m = arg(1); // 'loc'
        const double s = arg(2); // 'scale'
        const double T = arg(3); // 'mean_incub'
        const double i = arg(4); // 'incub_period'
        double tOST = x - i;     // time from onset to transmission

        const double norm_pdf = (s / a) * (1. + (i / T - 1.) * std::pow(1 + std::exp(m / s), -a) - (i / T) * std::pow(1. + std::exp((m + T) / s), -a));

        if (tOST < 0)
            tOST *= (T / i);

        return std::exp((m - tOST) / s) * std::pow(1. + std::exp((m - tOST) / s), -a - 1.) / norm_pdf;
    }

    double raw_cdf_(double x) const override final
    {
        const double a = arg(0);   // 'shape'
        const double m = arg(1);   // 'loc'
        const double s = arg(2);   // 'scale'
        const double T = arg(3);   // 'mean_incub'
        const double i = arg(4);   // 'incub_period'
        const double tOST = x - i; // time from onset to transmission

        const double norm_cdf = (1. + (i / T - 1.) * std::pow(1 + exp(m / s), -a) - (i / T) * (std::pow(1. + exp((m + T) / s), -a)));

        if (tOST < 0) // Integrate from -incub_period to tOST
            return (i / T) * (std::pow(1. + exp((m - tOST * T / i) / s), -a) - std::pow(1. + exp((m + T) / s), -a)) / norm_cdf;

        else
        {
            // Integrate from -incub to 0
            const double F0 = (i / T) * (std::pow(1. + std::exp(m / s), -a) - std::pow(1. + std::exp((m + T) / s), -a));
            // Integrate from 0 to t_ost
            return (F0 + std::pow(1. + std::exp((m - tOST) / s), -a) - std::pow(1. + std::exp(m / s), -a)) / norm_cdf;
        }
    }
};

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////  Register Infectivity Profiles  ////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

IMPLEMENT_HOUSEPI_DISTRIBUTION_REGISTRY(InfecProfileRegistry)
REGISTER_HOUSEPI_DISTRIBUTION(InfecProfileRegistry, InfecProfile_Ferretti20)

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////  Attribute : init function from name  /////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/// Initialize attribute as function given the function name.
void Attribute::init_function(const std::string &s, bool is_primitive)
{
    init_function(compose_func_registries(
        s, "function",
        {&UserFuncRegistry::inst(),
         &TranscendentalRegistry::inst(),
         &ArithmeticRegistry::inst(),
         &ContinuosDistrRegistry::inst(),
         &DiscreteDistrRegistry::inst(),
         &InfecProfileRegistry::inst()},
        is_primitive));
}
