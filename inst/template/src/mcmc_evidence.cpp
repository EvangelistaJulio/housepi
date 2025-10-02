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
 * mcmc_build_param.cpp
 *
 */

#include <atomic>

#include "mcmc_impl.h"

/// Welford online's algorithm for mean and variance.
/// Ref :
/// ----
/// https://en.wikipedia.org/wiki/Algorithms_for_calculating_variance
std::pair<double, double> compute_stats(const dbl_vec &data)
{
    size_t i = 0;
    size_t N = data.size();
    double x, delta, delta2;
    double sqr_diff = 0.;

    double avg = 0.;
    while (i < N)
    {
        x = data[i++];
        delta = x - avg;
        avg += delta / i;
        delta2 = x - avg;
        sqr_diff += delta * delta2;
    }

    return {avg, sqr_diff / (N - 1.)};
}

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////  MCMC : Impl : estimate logevidence  /////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/// Trapezoidal rule to estimate logarithm of marginal likelihood (evidence) with power posteriors.
///
/// Ref
/// -----
/// Friel, N., Hurn, M. & Wyse, J.
/// Improving power posterior estimation of statistical evidence.
/// Stat Comput 24, 709–723 (2014).
std::vector<dbl_vec> MCMC::Impl::power_posteriors_(const std::string &scheme, unsigned long n_rungs)
{
    // check scheme choice
    const str_vec options = {"adaptive", "power-fractions"};
    if (!(options % scheme))
        throw std::invalid_argument("Unknown power posterior scheme \"" + scheme + "\"." +
                                    tools::str::find_closest(scheme, options, HSP_MIN_STR_DIST));

    // check number of rungs for the thermodynamic ladder
    if (n_rungs < 2)
        throw std::invalid_argument("N. of rungs must be greater than 1.");
    setup.n_rungs = n_rungs; // NOTE: set before progress bar

    // initialize mcmc sampling
    init_monte_carlo_();

    // buffer user setup
    double buffer_inv_T = setup.inv_T;

    setup.prog_declare();        // declare current progress bar
    init_model_for_inference_(); // initialize model for inference
    augm.init();                 // initialize augmented variables
    setup.prog_open((setup.c_n_sweeps * setup.n_chains) * (setup.n_rungs + 1));

    // thermodynamic ladder estimates
    dbl_vec itm = {}; // inverse temperature
    dbl_vec avg = {}; // mean loglikelihood
    dbl_vec var = {}; // variance loglikelihood

    // run power posteriors
    power_post_zero_(itm, avg, var);
    if (scheme == "adaptive")
        power_post_adaptive_(itm, avg, var);
    else if (scheme == "power-fractions")
        power_post_fractions_(itm, avg, var);

    // close progress bar
    setup.prog_close();

    // restore user setup
    setup.inv_T = buffer_inv_T;

    return {itm, avg, var};
}

/// Power posteriors at 0 inverse temperature.
void MCMC::Impl::power_post_zero_(dbl_vec &itm, dbl_vec &avg, dbl_vec &var)
{
    // UPDATE: should I use previous inference for itm = 1 if available?
    // UPDATE: always under-estimate time if starting from itm = 0

    setup.rung = 0;
    setup.inv_T = 0.0;
    run_monte_carlo_();
    std::pair<double, double> stats = compute_stats(smpl_full_loglike_);

    itm.push_back(setup.inv_T);
    avg.push_back(stats.first);
    var.push_back(stats.second);
}

/// Power posteriors with adaptive scheme for the thermodynamic ladder.
void MCMC::Impl::power_post_adaptive_(dbl_vec &itm, dbl_vec &avg, dbl_vec &var)
{
    // Inverse Temperature 1
    setup.rung = 1;
    setup.inv_T = 1.0;
    run_monte_carlo_();
    std::pair<double, double> stats = compute_stats(smpl_full_loglike_);

    itm.push_back(setup.inv_T);
    avg.push_back(stats.first);
    var.push_back(stats.second);

    // Inverse Temperatures (0, 1)
    size_t k = 0; // index to interval lower bound
    bool is_right = false;
    while (setup.rung < setup.n_rungs)
    {
        ++setup.rung;

        // adaptive choice of inverse temperature
        // intersection between tangents
        double new_t = (avg[k + 1] - itm[k + 1] * var[k + 1] - avg[k] + itm[k] * var[k]) / (var[k] - var[k + 1]);
        // otherwise weighted average
        if ((new_t <= itm[k]) || (new_t >= itm[k + 1]))
            new_t = itm[k] + var[k + 1] * (itm[k + 1] - itm[k]) / (var[k] + var[k + 1]);

        // power posterior sampling
        setup.inv_T = new_t;
        run_monte_carlo_();
        stats = compute_stats(smpl_full_loglike_);

        // find larger rectangular discretization error
        is_right = (itm[k + 1] - setup.inv_T) * (avg[k + 1] - stats.first) > (new_t - itm[k]) * (stats.first - avg[k]);

        // insert in-between for all lists
        itm.insert(itm.begin() + k + 1, setup.inv_T);
        avg.insert(avg.begin() + k + 1, stats.first);
        var.insert(var.begin() + k + 1, stats.second);

        // if right interval is chosen, index points to inserted as lower bound
        k += (is_right);
    }
}

/// Power posteriors with power fractions scheme for the thermodynamic ladder.
void MCMC::Impl::power_post_fractions_(dbl_vec &itm, dbl_vec &avg, dbl_vec &var)
{
    for (setup.rung = 1; setup.rung <= setup.n_rungs; ++setup.rung)
    {
        // inverse tempearture = (rung / #rungs)^5
        setup.inv_T = std::pow(static_cast<double>(setup.rung) / setup.n_rungs, 5);
        run_monte_carlo_();

        std::pair<double, double> stats = compute_stats(smpl_full_loglike_);
        itm.push_back(setup.inv_T);
        avg.push_back(stats.first);
        var.push_back(stats.second);
    }
}

#ifdef _HOUSEPI_RCPP
/// Estimate logevidence according to a given scheme and number of rungs.
/// Uses trapezoidal rule to estimate logarithm of marginal likelihood (evidence).
///
/// Ref.
/// -----
/// Ref: Friel, N., Hurn, M. & Wyse, J.
/// Improving power posterior estimation of statistical evidence.
/// Stat Comput 24, 709–723 (2014).
Rcpp::List MCMC::Impl::logevidence(const std::string &scheme, unsigned long n_rungs)
{
    using namespace Rcpp;

    try
    {
        double output = 0;
        const std::vector<dbl_vec> pp = power_posteriors_(scheme, n_rungs);
        const dbl_vec &inv_tem = pp[0];
        const dbl_vec &avg_vec = pp[1];
        const dbl_vec &var_vec = pp[2];

        for (unsigned long k = 0; k < setup.n_rungs; ++k)
        {
            double d_inv_tem = inv_tem[k + 1] - inv_tem[k];
            output += (avg_vec[k + 1] + avg_vec[k]) * d_inv_tem / 2.;
            output += (var_vec[k] - var_vec[k + 1]) * d_inv_tem * d_inv_tem / 12.;
        }

        return List::create(_["estimate"] = output,
                            _["scheme"] = scheme,
                            _["n_rungs"] = setup.n_rungs,
                            _["inv_temper"] = wrap(inv_tem),
                            _["mean_loglike"] = wrap(avg_vec),
                            _["var_loglike"] = wrap(var_vec));
    }
    catch (const std::exception &e)
    {
        tools::error("While estimating logevidence with power posteriors:\n" +
                     std::string(e.what()));
    }
    return List();
}

#else
/// Estimate logevidence according to a given scheme and number of rungs.
/// Uses trapezoidal rule to estimate logarithm of marginal likelihood (evidence).
///
/// Ref.
/// -----
/// Ref: Friel, N., Hurn, M. & Wyse, J.
/// Improving power posterior estimation of statistical evidence.
/// Stat Comput 24, 709–723 (2014).
void MCMC::Impl::print_logevidence(const std::string &scheme, unsigned long n_rungs)
{
    try
    {
        double output = 0;
        const std::vector<dbl_vec> pp = power_posteriors_(scheme, n_rungs);
        const dbl_vec &inv_tem = pp[0];
        const dbl_vec &avg_vec = pp[1];
        const dbl_vec &var_vec = pp[2];

        for (unsigned long k = 0; k < setup.n_rungs; ++k)
        {
            double d_inv_tem = inv_tem[k + 1] - inv_tem[k];
            output += (avg_vec[k + 1] + avg_vec[k]) * d_inv_tem / 2.;
            output += (var_vec[k] - var_vec[k + 1]) * d_inv_tem * d_inv_tem / 12.;
        }

        // print result to stdout
        std::cout << output << std::endl;
    }
    catch (const std::exception &e)
    {
        tools::error("While estimating logevidence with power posteriors:\n" +
                     std::string(e.what()));
    }
}
#endif
