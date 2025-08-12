/*      __
 *     / /_  ____  __  __ ____ __  ____  __
 *    / __ \/ __ \/ / / / ___/ _ \/ __ \/ /
 *   / / / / /_/ / /_/ (__  )  __/ /_/ / /
 *  /_/ /_/\____/\____/____/\___/ ____/_/
 *                             /_/
 *
 * Copyright (C) 2025 Francesco Camaglia, Institut Pasteur
 *
 * distribution.cpp
 *
 */

#include "numerical_recipes.h"
#include "distribution.h"

//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//

/////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////  Uniform Distribution : class  //////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

/* Uniform Distribution derived class. */
class Uniform_Distribution : public Distribution
{
public:
    //===== CONSTRUCTOR =====//
    Uniform_Distribution()
        : Distribution("unif", {"min", "max"}, {0., 1.},
                       HSP_NA_DOUBLE, HSP_NA_DOUBLE, false, false),
          distr_gen(arg(0), arg(1))
    {
        init_uniform_range(); // ensure range matches arguments
    }
    ~Uniform_Distribution() = default;

    //===== SETTER =====//
    void init_reshaping(bool discrete) override final
    {
        unsigned short discr_case = (discrete) ? D_DISCR_TYPE_CNTR : D_DISCR_TYPE_NONE;

        init_functional_range(false, false);
        init_functional_discretization(discr_case);
        init_functional_truncation(discr_case);
    }

protected:
    std::uniform_real_distribution<double> distr_gen;

    //===== GETTER =====//
    double raw_pdf(double) const override final { return 1.; }
    double raw_cdf(double x) const override final { return x; }
    double raw_logpdf(double) const override final { return 0.; }
    double raw_logratio(double, double) const override final { return 0.; }
    double raw_rand(rnd_gen &prng) override final
    {
        const double r = distr_gen(prng);
        return arg(0) * (1 - r) + arg(1) * r;
    }
};

/////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////  Normal Distribution : class  //////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

/* Normal Distribution derived class. */
class Normal_Distribution : public Distribution
{
public:
    //===== CONSTRUCTOR =====//
    Normal_Distribution()
        : Distribution("norm", {"mean", "sd"}, {0., 1.},
                       HSP_NA_DOUBLE, HSP_NA_DOUBLE, false, false),
          distr_gen(arg(0), arg(1)) {}
    ~Normal_Distribution() = default;

    //===== SETTER =====//
    void init_reshaping(bool discrete) override final
    {
        init_functional_range(false, false);

        unsigned short discr_case = D_DISCR_TYPE_NONE;
        if (discrete)
            switch (which_valid_case())
            {
            case 1: //----- lower range truncation -----//
                discr_case = D_DISCR_TYPE_LEFT;
                break;
            case 3: //----- both sides truncation -----//
                discr_case = D_DISCR_TYPE_LEFT;
                break;
            default:
                discr_case = D_DISCR_TYPE_CNTR;
            }
        init_functional_discretization(discr_case);
        init_functional_truncation(discr_case);
    }

protected:
    std::normal_distribution<double> distr_gen;

    //===== GETTER =====//
    double raw_pdf(double x) const override final { return numr::normal::pdf(x, arg(0), arg(1)); }
    double raw_cdf(double x) const override final { return numr::normal::cdf(x, arg(0), arg(1)); }
    double raw_logpdf(double x) const override final { return numr::normal::logpdf(x, arg(0), arg(1)); }
    double raw_logratio(double xold, double xnew) const override final { return numr::normal::logratio(xold, xnew, arg(0), arg(1)); }
    double raw_rand(rnd_gen &prng) override final { return arg(0) + arg(1) * distr_gen(prng); }
};

/////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////  Lognormal Distribution : class  /////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

/* Lognormal Distribution derived class. */
class Lognormal_Distribution : public Distribution
{
public:
    //===== CONSTRUCTOR =====//
    Lognormal_Distribution()
        : Distribution("lnorm", {"meanlog", "sdlog"}, {0., 1.},
                       0., HSP_NA_DOUBLE, true, false),
          distr_gen(arg(0), arg(1)) {}
    ~Lognormal_Distribution() = default;

    //===== SETTER =====//
    void init_reshaping(bool discrete) override final { init_reshaping_min_zero(discrete, false); }

protected:
    std::lognormal_distribution<double> distr_gen;

    //===== SETTER =====//
    void update_args() override final
    {
        distr_gen.param(std::lognormal_distribution<double>::param_type(arg(0), arg(1)));
    }

    //===== GETTER =====//
    double raw_pdf(double x) const override final { return numr::lognormal::pdf(x, arg(0), arg(1)); }
    double raw_cdf(double x) const override final { return numr::lognormal::cdf(x, arg(0), arg(1)); }
    double raw_logpdf(double x) const override final { return numr::lognormal::logpdf(x, arg(0), arg(1)); }
    double raw_logratio(double xold, double xnew) const override final { return numr::lognormal::logratio(xold, xnew, arg(0), arg(1)); }
    double raw_rand(rnd_gen &prng) override final { return distr_gen(prng); }
};

/////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////  Beta Distribution : class  /////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

/* Beta Distribution derived class. */
class Beta_Distribution : public Distribution
{
public:
    //===== CONSTRUCTOR =====//
    Beta_Distribution()
        : Distribution("beta", {"shape1", "shape2"}, {1., 1.},
                       0., 1., true, true),
          gamma_1(arg(0), 1.), gamma_2(arg(1), 1.) {}
    ~Beta_Distribution() = default;

    //==== SETTER =====//
    void init_reshaping(bool discrete) override final
    {
        if (discrete)
            throw std::invalid_argument(
                "Invalid option discrete for distribution \"" + name() + "\".");
        init_functional_discretization(D_DISCR_TYPE_NONE);

        if (get_min() < 0)
            throw std::invalid_argument(
                "Invalid min value " + str_min() + " for distribution \"" + name() + "\".");
        else if (get_max() > 1)
            throw std::invalid_argument(
                "Invalid max value " + str_max() + " for distribution \"" + name() + "\".");
        else if ((get_min() > 0) || (get_max() < 1))
        {
            init_functional_range((get_min() == 0), (get_max() == 1));
            init_functional_truncation(D_DISCR_TYPE_NONE);
        }
    }

protected:
    std::gamma_distribution<double> gamma_1;
    std::gamma_distribution<double> gamma_2;

    //==== SETTER =====//
    void update_args() override final
    {
        gamma_1.param(std::gamma_distribution<double>::param_type(arg(0), 1.));
        gamma_2.param(std::gamma_distribution<double>::param_type(arg(1), 1.));
    }

    //==== GETTER =====//
    double raw_pdf(double x) const override final { return numr::beta::pdf(x, arg(0), arg(1)); }
    double raw_cdf(double x) const override final { return numr::beta::cdf(x, arg(0), arg(1)); }
    double raw_logpdf(double x) const override final
    {
        return numr::beta::logpdf(x, arg(0), arg(1));
    }
    double raw_logratio(double xold, double xnew) const override final { return numr::beta::logratio(xold, xnew, arg(0), arg(1)); }
    double raw_rand(rnd_gen &prng) override final
    {
        double x1 = gamma_1(prng);
        double x2 = gamma_2(prng);
        return x1 / (x1 + x2);
    }
};

/////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////  Gamma Distribution : class  /////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

/* Gamma Distribution derived class. */
class Gamma_Distribution : public Distribution
{
public:
    //===== CONSTRUCTOR =====//
    Gamma_Distribution()
        : Distribution("gamma", {"shape", "scale"}, {1., 1.},
                       0., HSP_NA_DOUBLE, false, false),
          distr_gen(arg(0), arg(1)) {}
    ~Gamma_Distribution() = default;

    //==== SETTER =====//
    void init_reshaping(bool discrete) override final { init_reshaping_min_zero(discrete, false); }

protected:
    std::gamma_distribution<double> distr_gen;

    //==== SETTER =====//
    void update_args() override final
    {
        distr_gen.param(std::gamma_distribution<double>::param_type(arg(0), arg(1)));
    }

    //==== GETTER =====//
    double raw_pdf(double x) const override final { return numr::gamma::pdf(x, arg(0), arg(1)); }
    double raw_cdf(double x) const override final { return numr::gamma::cdf(x, arg(0), arg(1)); }
    double raw_logpdf(double x) const override final { return numr::gamma::logpdf(x, arg(0), arg(1)); }
    double raw_logratio(double xold, double xnew) const override final { return numr::gamma::logratio(xold, xnew, arg(0), arg(1)); }
    double raw_rand(rnd_gen &prng) override final { return distr_gen(prng); }
};

/////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////  Gamma (Mean-Varance) Distribution : class  /////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

/* Gamma (Mean-Varance) Distribution derived class. */
class GammaMeanVar_Distribution : public Distribution
{
public:
    //===== CONSTRUCTOR =====//
    GammaMeanVar_Distribution()
        : Distribution("gamma_meanvar", {"mean", "var"}, {1., 1.},
                       0., HSP_NA_DOUBLE, false, false),
          distr_gen(arg(0), arg(1)) {}
    ~GammaMeanVar_Distribution() = default;

    //==== SETTER =====//
    void init_reshaping(bool discrete) override final
    {
        init_reshaping_min_zero(discrete, false);
    }

protected:
    std::gamma_distribution<double> distr_gen;

    //==== SETTER =====//
    void update_args() override final
    {
        const double m = arg(0);
        const double v = arg(1);
        distr_gen.param(std::gamma_distribution<double>::param_type((m * m) / v, v / m));
    }

    //==== GETTER =====//
    double raw_pdf(double x) const override final
    {
        const double m = arg(0);
        const double v = arg(1);
        return numr::gamma::pdf(x, (m * m) / v, v / m);
    }
    double raw_cdf(double x) const override final
    {
        const double m = arg(0);
        const double v = arg(1);
        return numr::gamma::cdf(x, (m * m) / v, v / m);
    }
    double raw_logpdf(double x) const override final
    {
        const double m = arg(0);
        const double v = arg(1);
        return numr::gamma::logpdf(x, (m * m) / v, v / m);
    }
    double raw_logratio(double xold, double xnew) const override final
    {
        const double m = arg(0);
        const double v = arg(1);
        return numr::gamma::logratio(xold, xnew, (m * m) / v, v / m);
    }
    double raw_rand(rnd_gen &prng) override final
    {
        return distr_gen(prng);
    }
};

/////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////  Rayleigh Distribution : class  ////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

/* Rayleigh Distribution derived class. */
class Rayleigh_Distribution : public Distribution
{
public:
    //===== CONSTRUCTOR =====//
    Rayleigh_Distribution()
        : Distribution("rayleigh", {"scale"}, {1.},
                       0., HSP_NA_DOUBLE, false, false),
          unif_gen(0., 1.) {}
    ~Rayleigh_Distribution() = default;

    //===== SETTER =====//
    void init_reshaping(bool discrete) override final { init_reshaping_min_zero(discrete, false); }

protected:
    std::uniform_real_distribution<double> unif_gen;

    //===== GETTER =====//
    double raw_pdf(double x) const override final { return numr::rayleigh::pdf(x, arg(0)); }
    double raw_cdf(double x) const override final { return numr::rayleigh::cdf(x, arg(0)); }
    double raw_logpdf(double x) const override final { return numr::rayleigh::logpdf(x, arg(0)); }
    double raw_logratio(double xold, double xnew) const override final { return numr::rayleigh::logratio(xold, xnew, arg(0)); }
    double raw_rand(rnd_gen &prng) override final { return arg(0) * sqrt(-2 * std::log(unif_gen(prng))); }
};

/////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////  Weibull Distribution : class  /////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

/* Weibull Distribution derived class. */
class Weibull_Distribution : public Distribution
{
public:
    //===== CONSTRUCTOR =====//
    Weibull_Distribution()
        : Distribution("weibull", {"shape", "scale"}, {1., 1.},
                       0., HSP_NA_DOUBLE, false, false),
          distr_gen(arg(0), arg(1)) {}
    ~Weibull_Distribution() = default;

    //===== SETTER =====//
    void init_reshaping(bool discrete) override final { init_reshaping_min_zero(discrete, false); }

protected:
    std::weibull_distribution<double> distr_gen;

    //===== SETTER =====//
    void update_args() override final
    {
        distr_gen.param(std::weibull_distribution<double>::param_type(arg(0), arg(1)));
    }

    //===== GETTER =====//
    double raw_pdf(double x) const override final { return numr::weibull::pdf(x, arg(0), arg(1)); }
    double raw_cdf(double x) const override final { return numr::weibull::cdf(x, arg(0), arg(1)); }
    double raw_logpdf(double x) const override final { return numr::weibull::logpdf(x, arg(0), arg(1)); }
    double raw_logratio(double xold, double xnew) const override final { return numr::weibull::logratio(xold, xnew, arg(0), arg(1)); }
    double raw_rand(rnd_gen &prng) override final { return distr_gen(prng); }
};

/////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////  Exponential  Distribution : class  //////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

/* Exponential  Distribution derived class. */
class Exponential_Distribution : public Distribution
{
public:
    //===== CONSTRUCTOR =====//
    Exponential_Distribution()
        : Distribution("exp", {"rate"}, {1.},
                       0., HSP_NA_DOUBLE, false, false),
          distr_gen(arg(0)) {}
    ~Exponential_Distribution() = default;

    //===== SETTER =====//
    void init_reshaping(bool discrete) override final { init_reshaping_min_zero(discrete, false); }

protected:
    std::exponential_distribution<double> distr_gen;

    //===== SETTER =====//
    void update_args() override final { distr_gen.param(std::exponential_distribution<double>::param_type(arg(0))); }

    //===== GETTER =====//
    double raw_pdf(double x) const override final { return numr::exponential::pdf(x, arg(0)); }
    double raw_cdf(double x) const override final { return numr::exponential::cdf(x, arg(0)); }
    double raw_logpdf(double x) const override final { return numr::exponential::logpdf(x, arg(0)); }
    double raw_logratio(double xold, double xnew) const override final { return numr::exponential::logratio(xold, xnew, arg(0)); }
    double raw_rand(rnd_gen &prng) override final { return distr_gen(prng); }
};

//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//

/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////  Poisson Distribution : class  ////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

/* Poisson Distribution derived class. */
class Poisson_Distribution : public Distribution
{
public:
    //===== CONSTRUCTOR =====//
    Poisson_Distribution()
        : Distribution("pois", {"lambda"}, {1.},
                       0., HSP_NA_DOUBLE, false, false),
          distr_gen(1.)
    {
        init_functional_discretization(D_DISCR_TYPE_DFND);
    }
    ~Poisson_Distribution() = default;

    //===== SETTER =====//
    void init_reshaping(bool) override final { init_reshaping_min_zero_discr_defined(); }

protected:
    std::poisson_distribution<unsigned> distr_gen;

    //===== SETTER =====//
    void update_args() override final { distr_gen.param(std::poisson_distribution<unsigned>::param_type(arg(0))); }

    //===== GETTER =====//
    double raw_pdf(double x) const override final { return numr::poisson::pdf(x, arg(0)); }
    double raw_cdf(double x) const override final { return numr::poisson::cdf(x, arg(0)); }
    double raw_logpdf(double x) const override final { return numr::poisson::logpdf(x, arg(0)); }
    double raw_logratio(double xold, double xnew) const override final { return numr::poisson::logratio(xold, xnew, arg(0)); }
    double raw_rand(rnd_gen &prng) override final { return distr_gen(prng); }
};

/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////  Binomial Distribution : class  ///////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

/* Binomial Distribution derived class. */
class Binomial_Distribution : public Distribution
{
public:
    //===== CONSTRUCTOR =====//
    Binomial_Distribution()
        : Distribution("binom", {"size", "prob"}, {4, 0.5},
                       0., HSP_NA_DOUBLE, false, false),
          distr_gen(arg(0), arg(1))
    {
        init_functional_discretization(D_DISCR_TYPE_DFND);
        init_attr("max", HSP_ATTR_RANGE, get_ptr_to_arg(0));
        init_functional_range(false, false);
    }
    ~Binomial_Distribution() = default;

    //===== SETTER =====//
    void init_reshaping(bool) override final
    {
        if (get_min() < 0.)
            throw std::invalid_argument(
                "Invalid min value " + str_min() + " for distribution \"" + name() + "\".");
        else if ((get_min() > 0.) || (get_ptr_to_arg(0) != get_ptr_to_max()))
        {
            if (get_max() > arg(0))
                throw std::invalid_argument(
                    "Invalid max value " + str_max() + " for distribution \"" + name() + "\" with size " + std::to_string(arg(0)) + ".");
            init_functional_range(false, false);
            init_functional_truncation(D_DISCR_TYPE_DFND);
        }
    }

protected:
    std::binomial_distribution<unsigned> distr_gen;

    //===== SETTER =====//
    void update_args() override final
    {
        distr_gen.param(std::binomial_distribution<unsigned>::param_type(arg(0), arg(1)));
    }

    //===== GETTER =====//
    double raw_pdf(double x) const override final { return numr::binomial::pdf(x, arg(0), arg(1)); }
    double raw_cdf(double x) const override final { return numr::binomial::cdf(x, arg(0), arg(1)); }
    double raw_logpdf(double x) const override final { return numr::binomial::logpdf(x, arg(0), arg(1)); }
    double raw_logratio(double xold, double xnew) const override final { return numr::binomial::logratio(xold, xnew, arg(0), arg(1)); }
    double raw_rand(rnd_gen &prng) override final { return distr_gen(prng); }
};

/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////  NegativeBinomial Distribution : class  ///////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

/* NegativeBinomial Distribution derived class. */
class NegativeBinomial_Distribution : public Distribution
{
public:
    //===== CONSTRUCTOR =====//
    NegativeBinomial_Distribution()
        : Distribution("nbinom", {"size", "prob"}, {4, .5},
                       0., HSP_NA_DOUBLE, false, false),
          distr_gen(arg(0), arg(1))
    {
        init_functional_discretization(D_DISCR_TYPE_DFND);
    }
    ~NegativeBinomial_Distribution() = default;

    //===== SETTER =====//
    void init_reshaping(bool) override final { init_reshaping_min_zero_discr_defined(); }

protected:
    std::negative_binomial_distribution<unsigned> distr_gen;

    //===== SETTER =====//
    void update_args() override final { distr_gen.param(std::negative_binomial_distribution<unsigned>::param_type(arg(0), arg(1))); }

    //===== GETTER =====//
    double raw_pdf(double x) const override final { return numr::negbinomial::pdf(x, arg(0), arg(1)); }
    double raw_cdf(double x) const override final { return numr::negbinomial::cdf(x, arg(0), arg(1)); }
    double raw_logpdf(double x) const override final { return numr::negbinomial::logpdf(x, arg(0), arg(1)); }
    double raw_logratio(double xold, double xnew) const override final { return numr::negbinomial::logratio(xold, xnew, arg(0), arg(1)); }
    double raw_rand(rnd_gen &prng) override final { return distr_gen(prng); }
};

/////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////  Geometric Distribution : class  ///////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

/* Geometric Distribution derived class. */
class Geometric_Distribution : public Distribution
{
public:
    //===== CONSTRUCTOR =====//
    Geometric_Distribution()
        : Distribution("geom", {"proba"}, {0.5},
                       0., HSP_NA_DOUBLE, false, false),
          distr_gen(arg(0))
    {
        init_functional_discretization(D_DISCR_TYPE_DFND);
    }
    ~Geometric_Distribution() = default;

    //===== SETTER =====//
    void init_reshaping(bool) override final { init_reshaping_min_zero_discr_defined(); }

protected:
    std::geometric_distribution<unsigned> distr_gen;

    //===== SETTER =====//
    void update_args() override final
    {
        distr_gen.param(std::geometric_distribution<unsigned>::param_type(arg(0)));
    }

    //===== GETTER =====//
    double raw_pdf(double x) const override final { return numr::geometric::pdf(x, arg(0)); }
    double raw_cdf(double x) const override final { return numr::geometric::cdf(x, arg(0)); }
    double raw_logpdf(double x) const override final { return numr::geometric::logpdf(x, arg(0)); }
    double raw_logratio(double xold, double xnew) const override final { return numr::geometric::logratio(xold, xnew, arg(0)); }
    double raw_rand(rnd_gen &prng) override final { return distr_gen(prng); }
};

//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//

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

//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//

///////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////// InfecProfile_Ferretti20 : class ////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
/*
The infectivity profile is considered equal to the distribution
of the Time from Onset of Symptoms to Transmission 't_ost':
>	t_ost = <time infection infectee> - <time onset infector>
Its distribution is modeled as a skew-logistic distribution (generalized logistic type I),
where the pre-symptomatic values are rescaled as a function of the incubation time 'incub'.
The arguments to provide are :
  0) "shape";
  1) "loc";
  2) "scale";
  3) "mean_incub" (mean incubation period);
  4) "incub_period" (individual incubation period);
Note: symtoms can occur after infection (pre-symptomatic transmission),
i. e. 't_ost' can be negative.
----
Ref:
----
Ferretti L, Ledda A, Wymant C, Zhao L, Ledda V, Abeler-Dorner L, et al.
The timing of COVID-19 transmission.
medRxiv; 2020. p. 2020.09.04.20188516. doi:10.1101/2020.09.04.20188516
 */
class InfecProfile_Ferretti20 : public Distribution
{
public:
    //===== CONSTRUCTOR =====//
    InfecProfile_Ferretti20()
        : Distribution("ferretti20",
                       {"shape", "loc", "scale", "mean_incub", "incub_period"},
                       {5.85, -4.0, 1.85, 5.8, 5.8},
                       0., HSP_NA_DOUBLE, false, false) {}
    ~InfecProfile_Ferretti20() = default;

    //===== SETTER =====//
    void init_reshaping(bool discrete) override final { init_reshaping_min_zero(discrete, false); }

protected:
    double raw_pdf(double x) const override final
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

    double raw_cdf(double x) const override final
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

//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//

IMPLEMENT_HOUSEPI_DISTRIBUTION_REGISTRY(InfecProfileRegistry)
REGISTER_HOUSEPI_DISTRIBUTION(InfecProfileRegistry, InfecProfile_Ferretti20)
