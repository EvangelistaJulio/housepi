#include <cmath>
#include <random>

// FIXME : this package is a bottleneck for portability
#include <boost/math/distributions/lognormal.hpp>
#include <boost/math/distributions/normal.hpp>
#include <boost/math/distributions/gamma.hpp>
#include "utils.h"

using namespace std;

// FIXME : are these MCMC variables?
double mIncub = 1.63;
double sdIncub = 0.5;
double mu_f = -4.0;
double sigma_f = 1.85; 
double alpha_f = 5.85;
double tau = exp(mIncub + 0.5 * pow(0.5, 2));

double cf = alpha_f / (sigma_f*(1- pow(1+exp((mu_f+tau)/sigma_f), -alpha_f)));

//----------------------------------------------------------------------
// Cumulative probability density functions
//----------------------------------------------------------------------
double plnorm (double x, double mean, double sd) {
    // FIXME : comment
    boost::math::lognormal_distribution<> d(mean, sd);
    return cdf(d, x);
}

double pgamma (double x, double shape, double scale) {
    // FIXME : comment
    boost::math::gamma_distribution<> d(shape, scale);
    return cdf(d, x);
}

double ptost (double tost, double inc) {
    // FIXME : comment
    double out = 0.0;

    if (tost < 0) {
        out = inc / tau * (pow(1 + exp(-(tost * tau / inc - mu_f) / sigma_f), -alpha_f) - pow(1 + exp((mu_f + tau) / sigma_f), -alpha_f));
    } else {
        // Integrate from -inc to 0
        out = inc / tau * (pow(1 + exp(mu_f / sigma_f), -alpha_f) - pow(1 + exp((mu_f + tau) / sigma_f), -alpha_f));

        // Integrate from 0 to tost
        out += pow(1 + exp(-(tost-mu_f) / sigma_f), -alpha_f) - pow(1 + exp(mu_f / sigma_f), -alpha_f) ;
    }
    out *= cf * sigma_f / alpha_f;

    return out;
}

double dtost (double tost, double inc) {
    // FIXME : comment
    double out = 0.0;

    if (tost < 0) {
        out = exp(-(tost * tau/inc - mu_f) / sigma_f) / pow(1 + exp(-(tost * tau / inc-mu_f) / sigma_f), alpha_f + 1);
    } else {
        out = exp(-(tost - mu_f) / sigma_f) / pow(1 + exp(-(tost - mu_f) / sigma_f), alpha_f + 1);
    }
    out *= cf;

    return out;
}

//----------------------------------------------------------------------
// Probability density functions
//----------------------------------------------------------------------
double dlnorm (double x, double mean, double sd) {
    // FIXME : comment
    return exp( - pow(log(x) - mean, 2) / ( 2 * pow(sd,2) ) ) / (x * sd * 2.50662827); // FIXME : what is this constant?
}

double dgamma (double x, double shape, double scale) {
    // FIXME : comment
    double out = 0.0;
    boost::math::gamma_distribution<> dist(shape, scale);
    if (x > 0) {
        out = pdf(dist, x); // FIXME: where does pdf come from?
    } else {
        out = 0.0;
    }
    return out;
}

double logdgamma (double x, double shape, double scale) {
    // FIXME : comment // FIXME : just the logarithm of the previous 
    double out = 0.0;
    boost::math::gamma_distribution<> dist(shape, scale);
    if (x > 0) {
        out = log(pdf(dist, x)); // FIXME: where does pdf come from?
    } else {
        out = log(0.0); // FIXME: this numerical minus infinity
    }
    return out;
}

// Log probability density distributions
double logdexp (double x, double rate) {
    // FIXME : comment
    return log(rate) - rate * x;
}

double logdlnorm (double x, double mean, double sd) {
    // FIXME : comment
    return -log(x) - log(sd) - 0.918938533 - pow(log(x)-mean, 2) / (2.0 * pow(sd, 2));
}


//----------------------------------------------------------------------
// Random generator
//----------------------------------------------------------------------

// FIXME: are these disributions from boost?

double rnorm (std::mt19937_64& rnd_generator, double mean, double sd) {
    std::normal_distribution<double> dist(mean, sd);
    return dist(rnd_generator);
}

double runif (std::mt19937_64& rnd_generator, double lower_value, double upper_value) {
    std::uniform_real_distribution<double> dist(lower_value, upper_value);
    return dist(rnd_generator);
}

double rlnorm (std::mt19937_64& rnd_generator, double mean, double sd) {
    std::lognormal_distribution<double> dist(mean, sd);
    return dist(rnd_generator);
}

double rgamma (std::mt19937_64& rnd_generator, double shape, double scale) {
    std::gamma_distribution<double> dist(shape, scale);
    return dist(rnd_generator);
}

double rexp (std::mt19937_64& rnd_generator, double rate) {
    std::exponential_distribution<double> dist(rate);
    return dist(rnd_generator);
}