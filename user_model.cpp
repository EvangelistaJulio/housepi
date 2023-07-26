#include "user_model.hpp"

using namespace std;

//----------------------------------------------------------------------
// Cumulative probability density functions
//----------------------------------------------------------------------

double t_onset_sympt_to_transm_norm(double t_incub, std::vector<double> args) {
    // Normalization for t_onset_sympt_to_transm_pdf
    double norm;
    double tau = args[0];
    double mu_f = args[1];
    double sigma_f = args[2];
    double alpha_f = args[3];

    norm = 1. + (t_incub / tau - 1.) * std::pow(1 + std::exp(mu_f / sigma_f), -alpha_f) - t_incub / tau * (std::pow(1. + std::exp((mu_f + tau) / sigma_f), -alpha_f));
    norm *= sigma_f / alpha_f;
    return norm;
}

// dtost
double t_onset_sympt_to_transm_pdf (double tost, double t_incub, std::vector<double> args) {
    /* The Time from Onset of Symptoms to Transmission (TOST) is modeled as a skew-logistic distribution (type I).
    This distribution is also a function of the incubation time `t_incub`. 
    Note: symtoms can occur after infection (pre-symptomatic transmission), in other words `tost` can be negative.
    ----
    Ref:
    ----
    Ferretti, L. et al. The timing of COVID-19 transmission. 2020.09.04.20188516 
    Preprint at https://doi.org/10.1101/2020.09.04.20188516 (2020).*/
    double out = 0.0, norm;
    double tau = args[0];
    double mu_f = args[1];
    double sigma_f = args[2];
    double alpha_f = args[3];

    norm = t_onset_sympt_to_transm_norm(t_incub, args);

    if (tost < 0) {
        tost = tost * tau / t_incub;
    } 
    out = std::exp((mu_f - tost) / sigma_f) * std::pow(1. + std::exp((mu_f - tost) / sigma_f), - 1. - alpha_f);

    return out / norm;
}

// ptost
double t_onset_sympt_to_transm_cdf (double tost, double t_incub, std::vector<double> args) {
    // Cumulative distribution function for `t_onset_sympt_to_transm_pdf`.
    double out = 0.0, norm;
    double tau = args[0];
    double mu_f = args[1];
    double sigma_f = args[2];
    double alpha_f = args[3];

    norm = t_onset_sympt_to_transm_norm(t_incub, args);

    if (tost < 0) {
        // Integrate from -t_incub to tost
        out = t_incub / tau * (std::pow(1. + std::exp((mu_f - tost * tau / t_incub) / sigma_f), -alpha_f) - std::pow(1. + std::exp((mu_f + tau) / sigma_f), -alpha_f));
    } else {
        // Integrate from -t_incub to 0
        out = t_incub / tau * (std::pow(1. + std::exp(mu_f / sigma_f), -alpha_f) - std::pow(1. + std::exp((mu_f + tau) / sigma_f), -alpha_f));
        // Integrate from 0 to tost
        out += std::pow(1. + std::exp((mu_f - tost) / sigma_f), -alpha_f) - std::pow(1. + std::exp(mu_f / sigma_f), -alpha_f) ;
    }

    return out / norm;
}


//----------------------------------------------------------------------
// Probability density functions (and their logarithms)
//----------------------------------------------------------------------

double lognorm_pdf (double x, std::vector<double> args) {
    // (logarithm of the) lognormal probability density function 
    double mean = args[0];
    double sd = args[1];
    double pdf = 0.0;

    if (x < 0.0 || mean <= 0.0 || sd <= 0.0) {
        // FIXME raise a warning for wrong parameters
        pdf = 0.0; 
    } else {
        pdf = std::exp(- 0.5 * std::pow((std::log(x) - mean) / sd, 2)) * std::pow(x * sd, -1) * 0.398942280;
    }
    return pdf;
}

double log_lognorm_pdf (double x, std::vector<double> args) {
    // (logarithm of the) lognormal probability density function 
    double mean = args[0];
    double sd = args[1];
    double logpdf = 0.0;

    if (x < 0.0 || mean <= 0.0 || sd <= 0.0) {
        // FIXME raise a warning for wrong parameters
        logpdf = - INFINITY; 
    } else {
        logpdf = -std::log(x) - std::log(sd) - 0.5 * std::pow((std::log(x) - mean) / sd, 2) - 0.918938533;
    }
    return logpdf;
}

double gamma_pdf (double x, std::vector<double> args) {
    // lognormal probability density function 
    double shape = args[0];
    double scale = args[1];
    double pdf, gamma_value;

    // The pdf is defined for x >= 0 and shape, scale > 0.
    if (x < 0.0 || shape <= 0.0 || scale <= 0.0) {
        // FIXME raise a warning for wrong parameters
        pdf = 0.0; 
    } else {
        gamma_value = std::tgamma(shape); 
        pdf = std::exp(-x / scale) * std::pow(x, shape - 1.) / (gamma_value * std::pow(scale, shape));
    }
    return pdf;
}

double log_gamma_pdf (double x, std::vector<double> args) {
    // (logarithm of the) lognormal probability density function 
    // FIXME : just the logarithm of the previous ?
    double shape = args[0];
    double scale = args[1];
    double logpdf, loggamma_value;

    // The logarithm of the pdf is undefined for x < 0 and shape, scale <= 0.
    if (x < 0.0 || shape <= 0.0 || scale <= 0.0) {
        // FIXME raise a warning for wrong parameters
        logpdf = - INFINITY; 
    } else {
        loggamma_value = std::lgamma(shape); 
        logpdf = loggamma_value - shape * std::log(scale) + (shape - 1) * std::log(x) - x / scale;
    }
    return logpdf;
}


//----------------------------------------------------------------------
// Random generators
//----------------------------------------------------------------------

// FIXME: are these disributions from boost?
// are they supposed to generate a random number? or the return a probability?

double norm_rand (std::mt19937_64& rand_gen, std::vector<double> args) {
    double mean = args[0];
    double sd = args[1];
    std::normal_distribution<double> dist(mean, sd); // FIXME: BOX-MULLER
    return dist(rand_gen);
}

double unif_rand (std::mt19937_64& rand_gen, std::vector<double> args) {
    double lower_value = args[0];
    double upper_value = args[1];
    std::uniform_real_distribution<double> dist(lower_value, upper_value); 
    return dist(rand_gen);
}

double lognorm_rand (std::mt19937_64& rand_gen, std::vector<double> args) {
    double mean = args[0];
    double sd = args[1];
    std::lognormal_distribution<double> dist(mean, sd); // FIXME: boh...
    return dist(rand_gen);
}

double gamma_rand (std::mt19937_64& rand_gen, std::vector<double> args) {
    double shape = args[0];
    double scale = args[1];
    std::gamma_distribution<double> dist(shape, scale); 
    return dist(rand_gen);
}

double exp_rand (std::mt19937_64& rand_gen, std::vector<double> args) {
    double rate = args[0];
    std::exponential_distribution<double> dist(rate); 
    return dist(rand_gen);
}