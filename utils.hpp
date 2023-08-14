// user_model.hpp //

#ifndef DEF_USER_MODEL__H
#define DEF_USER_MODEL__H

#include <string>
#include <cmath>
#include <random>
#include <unordered_map>
#include <iostream>

#include "setup.hpp"

typedef std::unordered_map<std::string, double> kwargs;
typedef double (*pdf_ptr)(double x, kwargs args);
typedef double (*rand_ptr)(std::mt19937_64 &rand_gen, kwargs args);

// NOTE: in C++ these objects can be improved moking rcpp pdf classes

// time from onset of symptoms to transmission
double tost_normaliz_const(double t_incub, kwargs args);
double t_onset_sympt_to_transm_pdf(double tost, double t_incub, kwargs args);
double t_onset_sympt_to_transm_cdf(double tost, double t_incub, kwargs args);

// Probability distribution functions
double unif_pdf(double x, kwargs args = {{"low", 0.}, {"up", 1.}});
double log_unif_pdf(double x, kwargs args = {{"low", 0.}, {"up", 1.}});
double normal_pdf(double x, kwargs args = {{"mean", 0.}, {"sd", 1.}});
double log_normal_pdf(double x, kwargs args = {{"mean", 0.}, {"sd", 1.}});
double lognorm_pdf(double x, kwargs args = {{"mean", 0.}, {"sd", 1.}});
double log_lognorm_pdf(double x, kwargs args = {{"mean", 0.}, {"sd", 1.}});
double gamma_pdf(double x, kwargs args = {{"shape", 0.}, {"scale", 1.}});
double log_gamma_pdf(double x, kwargs args = {{"shape", 0.}, {"scale", 1.}});
double exp_pdf(double x, kwargs args = {{"rate", 1.}});
double log_exp_pdf(double x, kwargs args = {{"rate", 1.}});

// Random generators
double unif_rand(std::mt19937_64 &rand_gen, kwargs args = {{"low", 0.}, {"up", 1.}});
double normal_rand(std::mt19937_64 &rand_gen, kwargs args = {{"mean", 0.}, {"sd", 1.}});
double lognorm_rand(std::mt19937_64 &rand_gen, kwargs args = {{"mean", 0.}, {"sd", 1.}});
double gamma_rand(std::mt19937_64 &rand_gen, kwargs args = {{"shape", 0.}, {"scale", 1.}});
double exp_rand(std::mt19937_64 &rand_gen, kwargs args = {{"rate", 1.}});

// Switchboard for distribution choice
void distribution_switchboard(std::string name, pdf_ptr* point_to_pdf, pdf_ptr* point_to_log_pdf, rand_ptr* point_to_rand);

// Linear interpolation for an indexed array
double lin_idx_interp(double x_value, const std::vector<double> &y);

#endif // DEF_USER_MODEL__H