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
#include "tools/discrete.h"
#include "tools/params.h"
#include "tools/string.h"
#include "tools/math_distr.h"

#include "distribution.h"
#include "distribution_registry.h"

#include <boost/random/beta_distribution.hpp>

using std::exp;
using std::log;
using std::pow;
using std::sqrt;

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
    const unsigned short discr_type = which_discrete_case_();

    if (range_lower_.is_empty() && range_upper_.is_empty()) // No Truncation
    {
        trunc_ = &Distribution::do_trunc_dummy_;
        logtrunc_ = &Distribution::do_logtrunc_dummy_;
        shift_ = &Distribution::do_shift_dummy_;
    }
    else if (!range_lower_.is_empty() && range_upper_.is_empty()) // Lower Range Truncation
    {
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
    }
    else if (range_lower_.is_empty() && !range_upper_.is_empty()) // Upper Range Truncation
    {
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
    }
    else // Both Sides Truncation
    {
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
    }
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
        init_discr_helper_(hsp::discr::LEFT); \
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
    }                                                               \
    double raw_quantile_(double p) const override final             \
    {                                                               \
        return numr::distr_::quantile(p, GET_ARGS_##n_args);        \
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
        // UPDATE: should understand when discrete defined...
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
        return arg(0) * (1. - r) + arg(1) * r;
    }
    double raw_quantile_(double p) const override final { return arg(0) * (1. - p) + arg(1) * p; }

private:
    /// Ensure uniform range (min, max) match argument (min, max).
    void init_uniform_range_()
    {
        // impose range minimum equal to argument minimum
        const Attribute &min_arg = arguments_[0];
        if (min_arg.is_function())
            range_lower_.transfer(min_arg);
        else
            range_lower_.transfer(min_arg);

        // impose range maximum equal to argument maximum
        const Attribute &max_arg = arguments_[1];
        if (max_arg.is_function())
            range_upper_.transfer(max_arg);
        else
            range_upper_.transfer(max_arg);
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
        if (range_lower_.is_empty())
            init_discr_helper_(hsp::discr::CENTER);
        else
            init_discr_helper_(hsp::discr::LEFT);
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
    double raw_rand_(rnd_gen &prng) override final { return exp(arg(0) + arg(1) * norm_gen(prng)); }
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
                       0., true, 1., true) {}

    void init_discretization() override final
    {
        throw std::invalid_argument(
            "Invoked invalid option discrete for " + err_name() + ".");
    }
    void init_truncation() override final
    {
        if (lower() < 0)
            throw std::invalid_argument(
                "Invoked invalid min value " + hsp::str::format(lower()) + " for " + err_name() + ".");
        if (upper() > 1)
            throw std::invalid_argument(
                "Invoked invalid max value " + hsp::str::format(upper()) + " for " + err_name() + ".");
        if (lower() > 0 || upper() < 1)
        {
            init_range_helper_(lower() == 0, upper() == 1);
            init_trunc_helper_();
        }
    }

protected:
    MAKE_RAW_NUMR_IMPORT_FOR_DISTR(beta, 2)
    double raw_rand_(rnd_gen &prng) override final
    {
        return boost::random::beta_distribution<double>(arg(0), arg(1))(prng);
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
                       0., false, hsp::na::DOUBLE, false) {}

    INIT_DISCR_LEFT_MIN_ZERO_CLOSE()

protected:
    MAKE_RAW_NUMR_IMPORT_FOR_DISTR(gamma, 2)
    double raw_rand_(rnd_gen &prng) override final
    {
        return std::gamma_distribution<double>(arg(0), arg(1))(prng);
    }
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
                       0., false, hsp::na::DOUBLE, false) {}

    INIT_DISCR_LEFT_MIN_ZERO_CLOSE()

protected:
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
    double raw_quantile_(double p) const override final
    {
        GAMMA_MEANVAR_CONVERTER()
        return numr::gamma::quantile(p, shape, scale);
    }
    double raw_rand_(rnd_gen &prng) override final
    {
        GAMMA_MEANVAR_CONVERTER()
        return std::gamma_distribution<double>(shape, scale)(prng);
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
    double raw_rand_(rnd_gen &prng) override final { return arg(0) * sqrt(-2 * log(unif_gen(prng))); }
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
    double raw_rand_(rnd_gen &prng) override final
    {
        return arg(1) * pow(-log(unif_gen(prng)), 1.0 / arg(0));
    }
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
    double raw_rand_(rnd_gen &prng) override final { return -log(unif_gen(prng)) / arg(0); }
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
                       0., false, hsp::na::DOUBLE, false) { init_discr_helper_(hsp::discr::DEFINED); }

    INIT_DISCR_DEFINED_MIN_ZERO_CLOSE()

protected:
    MAKE_RAW_NUMR_IMPORT_FOR_DISTR(poisson, 1)
    double raw_rand_(rnd_gen &prng) override final
    {
        return std::poisson_distribution<unsigned>(arg(0))(prng);
    }
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
                       0., false, hsp::na::DOUBLE, false)
    {
        init_discr_helper_(hsp::discr::DEFINED);
        range_upper_.init_pointer(get_arg(0).get_ptr());
        init_range_helper_(false, false);
    }

    void init_discretization() override final {}
    void init_truncation() override final
    {
        if (lower() < 0.)
            throw std::invalid_argument(
                "Invoked invalid min value " + hsp::str::format(range_lower_.get()) +
                " for " + err_name() + ".");
        if (lower() > 0. || get_arg(0).get_ptr() != get_upper().get_ptr())
        {
            if (upper() > arg(0))
                throw std::invalid_argument(
                    "Invoked invalid max value " + hsp::str::format(range_upper_.get()) +
                    " for " + err_name() + " with size " + std::to_string(arg(0)) + ".");
            init_range_helper_(false, false);
            init_trunc_helper_();
        }
    }

protected:
    MAKE_RAW_NUMR_IMPORT_FOR_DISTR(binomial, 2)
    double raw_rand_(rnd_gen &prng) override final
    {
        return std::binomial_distribution<unsigned>(arg(0), arg(1))(prng);
    }
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
                       0., false, hsp::na::DOUBLE, false) { init_discr_helper_(hsp::discr::DEFINED); }

    INIT_DISCR_DEFINED_MIN_ZERO_CLOSE()

protected:
    MAKE_RAW_NUMR_IMPORT_FOR_DISTR(negbinomial, 2)
    double raw_rand_(rnd_gen &prng) override final
    {
        return std::negative_binomial_distribution<unsigned>(arg(0), arg(1))(prng);
    }
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
                       0., false, hsp::na::DOUBLE, false) { init_discr_helper_(hsp::discr::DEFINED); }

    INIT_DISCR_DEFINED_MIN_ZERO_CLOSE()

protected:
    MAKE_RAW_NUMR_IMPORT_FOR_DISTR(geometric, 1)
    double raw_rand_(rnd_gen &prng) override final { return std::geometric_distribution<unsigned>(arg(0))(prng); }
};

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

        const double norm_pdf = (s / a) * (1. + (i / T - 1.) * pow(1 + exp(m / s), -a) - (i / T) * pow(1. + exp((m + T) / s), -a));

        if (tOST < 0)
            tOST *= (T / i);

        return exp((m - tOST) / s) * pow(1. + exp((m - tOST) / s), -a - 1.) / norm_pdf;
    }

    double raw_cdf_(double x) const override final
    {
        const double a = arg(0);   // 'shape'
        const double m = arg(1);   // 'loc'
        const double s = arg(2);   // 'scale'
        const double T = arg(3);   // 'mean_incub'
        const double i = arg(4);   // 'incub_period'
        const double tOST = x - i; // time from onset to transmission

        const double norm_cdf = (1. + (i / T - 1.) * pow(1 + exp(m / s), -a) - (i / T) * (pow(1. + exp((m + T) / s), -a)));

        if (tOST < 0) // Integrate from -incub_period to tOST
            return (i / T) * (pow(1. + exp((m - tOST * T / i) / s), -a) - pow(1. + exp((m + T) / s), -a)) / norm_cdf;

        // Integrate from -incub to 0
        const double F0 = (i / T) * (pow(1. + exp(m / s), -a) - pow(1. + exp((m + T) / s), -a));
        // Integrate from 0 to t_ost
        return (F0 + pow(1. + exp((m - tOST) / s), -a) - pow(1. + exp(m / s), -a)) / norm_cdf;
    }
};

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

std::unique_ptr<Distribution> DistributionRegistry::get_distr(const std::string &opt) const
{
    if (factory_map_.count(opt))
        return factory_map_.find(opt)->second();
    return nullptr;
}

std::unique_ptr<Function> DistributionRegistry::get(const std::string &opt, bool is_primitive) const
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

std::unique_ptr<Distribution> compose_distr_registries(const std::string &opt,
                                                       const std::string &name,
                                                       const std::vector<const DistributionRegistry *> &registries)
{
    str_vec options = {};

    if (registries.empty())
        throw std::invalid_argument(
            "No registries were provided for '" + name + "'.");

    if (hsp::na::is(opt))
        return std::make_unique<Distribution>();

    for (auto &r : registries)
    {
        if (r->get_options() % opt)
            return r->get_distr(opt);
        options += r->get_options();
    }

    throw std::invalid_argument(
        "Unknown " + name + " '" + opt + "'." + hsp::str::find_closest(opt, options));
}

///////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////  Register Distributions  ///////////////////////////////////
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
//////////////////////////////////  Register Infectivity Profiles  ////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

IMPLEMENT_HOUSEPI_DISTRIBUTION_REGISTRY(InfecProfileRegistry)
REGISTER_HOUSEPI_DISTRIBUTION(InfecProfileRegistry, InfecProfile_Ferretti20)

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////  Attribute : init function from name  ////////////////////////////
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
