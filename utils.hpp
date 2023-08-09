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

// NOTE: in C++ these objects can be improved moking rcpp pdf classes 

// time from onset of symptoms to transmission
double tost_normaliz_const(double t_incub, kwargs args);
double t_onset_sympt_to_transm_pdf (double tost, double t_incub, kwargs args);
double t_onset_sympt_to_transm_cdf (double tost, double t_incub, kwargs args);

// Probability distribution functions
double unif_pdf(double x, kwargs args = {{"low", 0.}, {"up", 1.}});
double log_unif_pdf(double x, kwargs args = {{"low", 0.}, {"up", 1.}});
double normal_pdf(double x, kwargs args = {{"mean", 0.}, {"sd", 1.}});
double log_normal_pdf(double x, kwargs args = {{"mean", 0.}, {"sd", 1.}});
double lognorm_pdf(double x, kwargs args = {{"mean", 0.}, {"sd", 1.}});
double log_lognorm_pdf(double x, kwargs args = {{"mean", 0.}, {"sd", 1.}});
double gamma_pdf(double x, kwargs args = {{"shape", 0.}, {"scale", 1.}});
double log_gamma_pdf(double x, kwargs args = {{"shape", 0.}, {"scale", 1.}});

// Random generators
double unif_rand(std::mt19937_64& rand_gen, kwargs args = {{"low", 0.}, {"up", 1.}});
double normal_rand(std::mt19937_64& rand_gen, kwargs args = {{"mean", 0.}, {"sd", 1.}});
double lognorm_rand(std::mt19937_64& rand_gen, kwargs args = {{"mean", 0.}, {"sd", 1.}});
double gamma_rand(std::mt19937_64& rand_gen, kwargs args = {{"shape", 0.}, {"scale", 1.}});
double exp_rand(std::mt19937_64& rand_gen, kwargs args = {{"rate", 1.}});

#endif // DEF_USER_MODEL__H