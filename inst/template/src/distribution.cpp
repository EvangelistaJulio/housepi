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
                       HOUSEPI_NA_DOUBLE, HOUSEPI_NA_DOUBLE, false, false),
          distr_gen(get_arg(0), get_arg(1))
    {
        init_uniform_range(); // ensure range matches arguments
    }
    ~Uniform_Distribution() = default;

    //===== SETTER =====//
    void init_reshaping(bool discrete) override final
    {
        ushort discr_case = (discrete) ? D_DISCR_TYPE_CNTR : D_DISCR_TYPE_NONE;

        init_functional_range(false, false);
        init_functional_discretization(discr_case);
        init_functional_truncation(discr_case);
    }

protected:
    std::uniform_real_distribution<double> distr_gen;

    //===== GETTER =====//
    double raw_pdf(double x) const override final
    {
        UNUSED(x);
        return 1.;
    }
    double raw_cdf(double x) const override final
    {
        return x;
    }
    double raw_logpdf(double x) const override final
    {
        UNUSED(x);
        return 0.;
    }
    double raw_logratio(double xold, double xnew) const override final
    {
        UNUSED(xold);
        UNUSED(xnew);
        return 0.;
    }
    double raw_rand(rnd_gen &prng) override final
    {
        const double r = distr_gen(prng);
        return get_arg(0) * (1 - r) + get_arg(1) * r;
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
                       HOUSEPI_NA_DOUBLE, HOUSEPI_NA_DOUBLE, false, false),
          distr_gen(get_arg(0), get_arg(1)) {}
    ~Normal_Distribution() = default;

    //===== SETTER =====//
    void init_reshaping(bool discrete) override final
    {
        init_functional_range(false, false);

        ushort discr_case = D_DISCR_TYPE_NONE;
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
    double raw_pdf(double x) const override final
    {
        return numr::normal::pdf(x, get_arg(0), get_arg(1));
    }
    double raw_cdf(double x) const override final
    {
        return numr::normal::cdf(x, get_arg(0), get_arg(1));
    }
    double raw_logpdf(double x) const override final
    {
        return numr::normal::logpdf(x, get_arg(0), get_arg(1));
    }
    double raw_logratio(double xold, double xnew) const override final
    {
        return numr::normal::logratio(xold, xnew, get_arg(0), get_arg(1));
    }
    double raw_rand(rnd_gen &prng) override final
    {
        return get_arg(0) + get_arg(1) * distr_gen(prng);
    }
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
                       0., HOUSEPI_NA_DOUBLE, false, false),
          distr_gen(get_arg(0), get_arg(1)) {}
    ~Lognormal_Distribution() = default;

    //===== SETTER =====//
    void init_reshaping(bool discrete) override final
    {
        init_min_zero_reshaping(discrete, false);
    }

protected:
    std::lognormal_distribution<double> distr_gen;

    //===== SETTER =====//
    void update_args() override final
    {
        distr_gen.param(std::lognormal_distribution<double>::param_type(get_arg(0), get_arg(1)));
    }

    //===== GETTER =====//
    double raw_pdf(double x) const override final
    {
        return numr::lognormal::pdf(x, get_arg(0), get_arg(1));
    }
    double raw_cdf(double x) const override final
    {
        return numr::lognormal::cdf(x, get_arg(0), get_arg(1));
    }
    double raw_logpdf(double x) const override final
    {
        return numr::lognormal::logpdf(x, get_arg(0), get_arg(1));
    }
    double raw_logratio(double xold, double xnew) const override final
    {
        return numr::lognormal::logratio(xold, xnew, get_arg(0), get_arg(1));
    }
    double raw_rand(rnd_gen &prng) override final
    {
        return distr_gen(prng);
    }
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
          gamma_1(get_arg(0), 1.), gamma_2(get_arg(1), 1.) {}
    ~Beta_Distribution() = default;

    //==== SETTER =====//
    void init_reshaping(bool discrete) override final
    {
        if (discrete)
            throw std::invalid_argument(
                "Invalid option discrete for distribution \"" + show_name() + "\".");
        init_functional_discretization(D_DISCR_TYPE_NONE);

        if (get_min() < 0)
            throw std::invalid_argument(
                "Invalid min value " + show_min() + " for distribution \"" + show_name() + "\".");
        else if (get_max() > 1)
            throw std::invalid_argument(
                "Invalid max value " + show_max() + " for distribution \"" + show_name() + "\".");
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
        gamma_1.param(std::gamma_distribution<double>::param_type(get_arg(0), 1.));
        gamma_2.param(std::gamma_distribution<double>::param_type(get_arg(1), 1.));
    }

    //==== GETTER =====//
    double raw_pdf(double x) const override final
    {
        return numr::beta::pdf(x, get_arg(0), get_arg(1));
    }
    double raw_cdf(double x) const override final
    {
        return numr::beta::cdf(x, get_arg(0), get_arg(1));
    }
    double raw_logpdf(double x) const override final
    {
        return numr::beta::logpdf(x, get_arg(0), get_arg(1));
    }
    double raw_logratio(double xold, double xnew) const override final
    {
        return numr::beta::logratio(xold, xnew, get_arg(0), get_arg(1));
    }
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
                       0., HOUSEPI_NA_DOUBLE, false, false),
          distr_gen(get_arg(0), get_arg(1)) {}
    ~Gamma_Distribution() = default;

    //==== SETTER =====//
    void init_reshaping(bool discrete) override final
    {
        init_min_zero_reshaping(discrete, false);
    }

protected:
    std::gamma_distribution<double> distr_gen;

    //==== SETTER =====//
    void update_args() override final
    {
        distr_gen.param(std::gamma_distribution<double>::param_type(get_arg(0), get_arg(1)));
    }

    //==== GETTER =====//
    double raw_pdf(double x) const override final
    {
        return numr::gamma::pdf(x, get_arg(0), get_arg(1));
    }
    double raw_cdf(double x) const override final
    {
        return numr::gamma::cdf(x, get_arg(0), get_arg(1));
    }
    double raw_logpdf(double x) const override final
    {
        return numr::gamma::logpdf(x, get_arg(0), get_arg(1));
    }
    double raw_logratio(double xold, double xnew) const override final
    {
        return numr::gamma::logratio(xold, xnew, get_arg(0), get_arg(1));
    }
    double raw_rand(rnd_gen &prng) override final
    {
        return distr_gen(prng);
    }
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
                       0., HOUSEPI_NA_DOUBLE, false, false),
          distr_gen(get_arg(0), get_arg(1)) {}
    ~GammaMeanVar_Distribution() = default;

    //==== SETTER =====//
    void init_reshaping(bool discrete) override final
    {
        init_min_zero_reshaping(discrete, false);
    }

protected:
    std::gamma_distribution<double> distr_gen;

    //==== SETTER =====//
    void update_args() override final
    {
        const double m = get_arg(0);
        const double v = get_arg(1);
        distr_gen.param(std::gamma_distribution<double>::param_type((m * m) / v, v / m));
    }

    //==== GETTER =====//
    double raw_pdf(double x) const override final
    {
        const double m = get_arg(0);
        const double v = get_arg(1);
        return numr::gamma::pdf(x, (m * m) / v, v / m);
    }
    double raw_cdf(double x) const override final
    {
        const double m = get_arg(0);
        const double v = get_arg(1);
        return numr::gamma::cdf(x, (m * m) / v, v / m);
    }
    double raw_logpdf(double x) const override final
    {
        const double m = get_arg(0);
        const double v = get_arg(1);
        return numr::gamma::logpdf(x, (m * m) / v, v / m);
    }
    double raw_logratio(double xold, double xnew) const override final
    {
        const double m = get_arg(0);
        const double v = get_arg(1);
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
                       0., HOUSEPI_NA_DOUBLE, false, false),
          unif_gen(0., 1.) {}
    ~Rayleigh_Distribution() = default;

    //===== SETTER =====//
    void init_reshaping(bool discrete) override final
    {
        init_min_zero_reshaping(discrete, false);
    }

protected:
    std::uniform_real_distribution<double> unif_gen;

    //===== GETTER =====//
    double raw_pdf(double x) const override final
    {
        return numr::rayleigh::pdf(x, get_arg(0));
    }
    double raw_cdf(double x) const override final
    {
        return numr::rayleigh::cdf(x, get_arg(0));
    }
    double raw_logpdf(double x) const override final
    {
        return numr::rayleigh::logpdf(x, get_arg(0));
    }
    double raw_logratio(double xold, double xnew) const override final
    {
        return numr::rayleigh::logratio(xold, xnew, get_arg(0));
    }
    double raw_rand(rnd_gen &prng) override final
    {
        return get_arg(0) * sqrt(-2 * std::log(unif_gen(prng)));
    }
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
                       0., HOUSEPI_NA_DOUBLE, false, false),
          distr_gen(get_arg(0), get_arg(1)) {}
    ~Weibull_Distribution() = default;

    //===== SETTER =====//
    void init_reshaping(bool discrete) override final
    {
        init_min_zero_reshaping(discrete, false);
    }

protected:
    std::weibull_distribution<double> distr_gen;

    //===== SETTER =====//
    void update_args() override final
    {
        distr_gen.param(std::weibull_distribution<double>::param_type(get_arg(0), get_arg(1)));
    }

    //===== GETTER =====//
    double raw_pdf(double x) const override final
    {
        return numr::weibull::pdf(x, get_arg(0), get_arg(1));
    }
    double raw_cdf(double x) const override final
    {
        return numr::weibull::cdf(x, get_arg(0), get_arg(1));
    }
    double raw_logpdf(double x) const override final
    {
        return numr::weibull::logpdf(x, get_arg(0), get_arg(1));
    }
    double raw_logratio(double xold, double xnew) const override final
    {
        return numr::weibull::logratio(xold, xnew, get_arg(0), get_arg(1));
    }
    double raw_rand(rnd_gen &prng) override final
    {
        return distr_gen(prng);
    }
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
                       0., HOUSEPI_NA_DOUBLE, false, false),
          distr_gen(get_arg(0)) {}
    ~Exponential_Distribution() = default;

    //===== SETTER =====//
    void init_reshaping(bool discrete) override final
    {
        init_min_zero_reshaping(discrete, false);
    }

protected:
    std::exponential_distribution<double> distr_gen;

    //===== SETTER =====//
    void update_args() override final
    {
        distr_gen.param(std::exponential_distribution<double>::param_type(get_arg(0)));
    }

    //===== GETTER =====//
    double raw_pdf(double x) const override final
    {
        return numr::exponential::pdf(x, get_arg(0));
    }
    double raw_cdf(double x) const override final
    {
        return numr::exponential::cdf(x, get_arg(0));
    }
    double raw_logpdf(double x) const override final
    {
        return numr::exponential::logpdf(x, get_arg(0));
    }
    double raw_logratio(double xold, double xnew) const override final
    {
        return numr::exponential::logratio(xold, xnew, get_arg(0));
    }
    double raw_rand(rnd_gen &prng) override final
    {
        return distr_gen(prng);
    }
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
                       0., HOUSEPI_NA_DOUBLE, false, false),
          distr_gen(1.)
    {
        init_functional_discretization(D_DISCR_TYPE_DFND);
    }
    ~Poisson_Distribution() = default;

    //===== SETTER =====//
    void init_reshaping(bool discrete) override final
    {
        UNUSED(discrete);
        init_min_zero_discr_defined_reshaping();
    }

protected:
    std::poisson_distribution<unsigned> distr_gen;

    //===== SETTER =====//
    void update_args() override final
    {
        distr_gen.param(std::poisson_distribution<unsigned>::param_type(get_arg(0)));
    }

    //===== GETTER =====//
    double raw_pdf(double x) const override final
    {
        return numr::poisson::pdf(x, get_arg(0));
    }
    double raw_cdf(double x) const override final
    {
        return numr::poisson::cdf(x, get_arg(0));
    }
    double raw_logpdf(double x) const override final
    {
        return numr::poisson::logpdf(x, get_arg(0));
    }
    double raw_logratio(double xold, double xnew) const override final
    {
        return numr::poisson::logratio(xold, xnew, get_arg(0));
    }
    double raw_rand(rnd_gen &prng) override final
    {
        return distr_gen(prng);
    }
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
                       0., HOUSEPI_NA_DOUBLE, false, false),
          distr_gen(get_arg(0), get_arg(1))
    {
        init_functional_discretization(D_DISCR_TYPE_DFND);
        assign_attr("max", ATTR_RANGE, get_ptr_to_arg(0));
        init_functional_range(false, false);
    }
    ~Binomial_Distribution() = default;

    //===== SETTER =====//
    void init_reshaping(bool discrete) override final
    {
        UNUSED(discrete);

        if (get_min() < 0.)
            throw std::invalid_argument(
                "Invalid min value " + show_min() + " for distribution \"" + show_name() + "\".");
        else if ((get_min() > 0.) || (get_ptr_to_arg(0) != get_ptr_to_max()))
        {
            if (get_max() > get_arg(0))
                throw std::invalid_argument(
                    "Invalid max value " + show_max() + " for distribution \"" + show_name() + "\" with size " + std::to_string(get_arg(0)) + ".");
            init_functional_range(false, false);
            init_functional_truncation(D_DISCR_TYPE_DFND);
        }
    }

protected:
    std::binomial_distribution<unsigned> distr_gen;

    //===== SETTER =====//
    void update_args() override final
    {
        distr_gen.param(std::binomial_distribution<unsigned>::param_type(get_arg(0), get_arg(1)));
    }

    //===== GETTER =====//
    double raw_pdf(double x) const override final
    {
        return numr::binomial::pdf(x, get_arg(0), get_arg(1));
    }
    double raw_cdf(double x) const override final
    {
        return numr::binomial::cdf(x, get_arg(0), get_arg(1));
    }
    double raw_logpdf(double x) const override final
    {
        return numr::binomial::logpdf(x, get_arg(0), get_arg(1));
    }
    double raw_logratio(double xold, double xnew) const override final
    {
        return numr::binomial::logratio(xold, xnew, get_arg(0), get_arg(1));
    }
    double raw_rand(rnd_gen &prng) override final
    {
        return distr_gen(prng);
    }
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
                       0., HOUSEPI_NA_DOUBLE, false, false),
          distr_gen(get_arg(0), get_arg(1))
    {
        init_functional_discretization(D_DISCR_TYPE_DFND);
    }
    ~NegativeBinomial_Distribution() = default;

    //===== SETTER =====//
    void init_reshaping(bool discrete) override final
    {
        UNUSED(discrete);
        init_min_zero_discr_defined_reshaping();
    }

protected:
    std::negative_binomial_distribution<unsigned> distr_gen;

    //===== SETTER =====//
    void update_args() override final
    {
        distr_gen.param(std::negative_binomial_distribution<unsigned>::param_type(get_arg(0), get_arg(1)));
    }

    //===== GETTER =====//
    double raw_pdf(double x) const override final
    {
        return numr::negbinomial::pdf(x, get_arg(0), get_arg(1));
    }
    double raw_cdf(double x) const override final
    {
        return numr::negbinomial::cdf(x, get_arg(0), get_arg(1));
    }
    double raw_logpdf(double x) const override final
    {
        return numr::negbinomial::logpdf(x, get_arg(0), get_arg(1));
    }
    double raw_logratio(double xold, double xnew) const override final
    {
        return numr::negbinomial::logratio(xold, xnew, get_arg(0), get_arg(1));
    }
    double raw_rand(rnd_gen &prng) override final
    {
        return distr_gen(prng);
    }
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
                       0., HOUSEPI_NA_DOUBLE, false, false),
          distr_gen(get_arg(0))
    {
        init_functional_discretization(D_DISCR_TYPE_DFND);
    }
    ~Geometric_Distribution() = default;

    //===== SETTER =====//
    void init_reshaping(bool discrete) override final
    {
        UNUSED(discrete);
        init_min_zero_discr_defined_reshaping();
    }

protected:
    std::geometric_distribution<unsigned> distr_gen;

    //===== SETTER =====//
    void update_args() override final
    {
        distr_gen.param(std::geometric_distribution<unsigned>::param_type(get_arg(0)));
    }

    //===== GETTER =====//
    double raw_pdf(double x) const override final
    {
        return numr::geometric::pdf(x, get_arg(0));
    }
    double raw_cdf(double x) const override final
    {
        return numr::geometric::cdf(x, get_arg(0));
    }
    double raw_logpdf(double x) const override final
    {
        return numr::geometric::logpdf(x, get_arg(0));
    }
    double raw_logratio(double xold, double xnew) const override final
    {
        return numr::geometric::logratio(xold, xnew, get_arg(0));
    }
    double raw_rand(rnd_gen &prng) override final
    {
        return distr_gen(prng);
    }
};

//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//

/////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////  Distribution Registry  ///////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

std::shared_ptr<Distribution> ContinuosDistrRegistry::get_distr(const std::string &opt) const
{
    if (opt == "unif")
        return std::make_shared<Uniform_Distribution>();

    else if (opt == "norm")
        return std::make_shared<Normal_Distribution>();

    else if (opt == "lnorm")
        return std::make_shared<Lognormal_Distribution>();

    else if (opt == "beta")
        return std::make_shared<Beta_Distribution>();

    else if (opt == "gamma")
        return std::make_shared<Gamma_Distribution>();

    else if (opt == "gamma_meanvar")
        return std::make_shared<GammaMeanVar_Distribution>();

    else if (opt == "rayleigh")
        return std::make_shared<Rayleigh_Distribution>();

    else if (opt == "weibull")
        return std::make_shared<Weibull_Distribution>();

    else if (opt == "exp")
        return std::make_shared<Exponential_Distribution>();

    else
        return nullptr;
}

std::shared_ptr<Distribution> DiscreteDistrRegistry::get_distr(const std::string &opt) const
{
    if (opt == "pois")
        return std::make_shared<Poisson_Distribution>();

    else if (opt == "binom")
        return std::make_shared<Binomial_Distribution>();

    else if (opt == "nbinom")
        return std::make_shared<NegativeBinomial_Distribution>();

    else if (opt == "geom")
        return std::make_shared<Geometric_Distribution>();

    else
        return nullptr;
}
