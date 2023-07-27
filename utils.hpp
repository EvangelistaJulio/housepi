// user_model.hpp //

#ifndef DEF_USER_MODEL__H
#define DEF_USER_MODEL__H

#include <cmath>
#include <random>

// time from onset of symptoms to transmission
double tost_norm_const(double t_incub, std::vector<double> args);
double t_onset_sympt_to_transm_pdf (double tost, double t_incub, std::vector<double> args);
double t_onset_sympt_to_transm_cdf (double tost, double t_incub, std::vector<double> args);

// Probability distribution functions
double normal_pdf(double x, std::vector<double> args);
double log_normal_pdf(double x, std::vector<double> args);
double lognorm_pdf(double x, std::vector<double> args);
double log_lognorm_pdf(double x, std::vector<double> args);
double gamma_pdf(double x, std::vector<double> args);
double log_gamma_pdf(double x, std::vector<double> args);

// Random generators
double normal_rand(std::mt19937_64& rand_gen, std::vector<double> args = {0., 1.});
double unif_rand(std::mt19937_64& rand_gen, std::vector<double> args = {0., 1.});
double lognorm_rand(std::mt19937_64& rand_gen, std::vector<double> args = {0., 1.});
double gamma_rand(std::mt19937_64& rand_gen, std::vector<double> args = {0., 1.});
double exp_rand(std::mt19937_64& rand_gen, std::vector<double> args = {1.});

#endif // DEF_USER_MODEL__H