// user_model.hpp //

#ifndef DEF_USER_MODEL__H
#define DEF_USER_MODEL__H

#include <string>
#include <cmath>
#include <random>
#include <unordered_map>
#include <iostream>

#include "setup.hpp"

#ifdef CONTINUOS_TIME
typedef double timetype;
#else 
typedef int timetype;
//UPDATE : Distribution turn into discrete LUT
#endif

typedef std::unordered_map<std::string, double> kwargs;
typedef std::unordered_map<std::string, double*> kwargs_ext;

// Linear interpolation for an indexed array
double lin_idx_interp(double x_value, const std::vector<double> &y);

///////////////////////////////////////////////////////////////////////
// ////////////////////// Distribution : class ///////////////////// //
///////////////////////////////////////////////////////////////////////

class Distribution
{
protected:
    kwargs args = {};               // internal arguments (fixed at initialization)
    kwargs_ext args_ext = {};       // external arguments (pointer to variable)

public:
    // Constructors
    Distribution(kwargs new_args);
    Distribution(kwargs_ext new_args_ext);
    Distribution() = default;
    ~Distribution() {}

    // Getter (fixed)
    double get_arg(std::string argID) { return (this->*_get_arg_func_)(argID); }

    // Setter (just for fixed)
    void set_arg(std::string argID, double new_value) { args.at(argID) = new_value; }

private :
	double (Distribution::*_get_arg_func_)(std::string argID);
	double _get_arg_int_(std::string argID) { return args.at(argID); }
	double _get_arg_ext_(std::string argID) { return *(args_ext.at(argID)); }
};

///////////////////////////////////////////////////////////////////////
// /////////////////// Distribution_1D : class ///////////////////// //
///////////////////////////////////////////////////////////////////////

class Distribution_1D : public Distribution
{
public:
    using Distribution::Distribution;

    // Virtual functions (dummy)
    virtual double pdf(double x) { return x; }
    virtual double log_pdf(double x) { return std::log(pdf(x)); }
    virtual double log_ratio(double x, double y) { return std::log(pdf(x)) - std::log(pdf(y)); }
    virtual double rand(std::mt19937_64 &rng)
    {
        std::uniform_real_distribution<double> distr_gen(0.0, 1.0);
        return distr_gen(rng);
    }
};

///////////////////////////////////////////////////////////////////////
// ///////////////////// Uniform_Distr : class ///////////////////// //
///////////////////////////////////////////////////////////////////////

class Uniform_Distr : public Distribution_1D
{
private:
    std::uniform_real_distribution<double> distr_gen;

public:
    Uniform_Distr(kwargs new_args = {{"low", (double)0.0}, {"up", (double)1.0}}) : Distribution_1D(new_args), distr_gen(get_arg("low"), get_arg("up")) {}
    Uniform_Distr(kwargs_ext new_args_ext) : Distribution_1D(new_args_ext), distr_gen(get_arg("low"), get_arg("up")) {}
    ~Uniform_Distr() {}

    double pdf(double x) override;
    double log_pdf(double x) override;
    double log_ratio(double x, double y) override;
    double rand(std::mt19937_64 &rng) override { return distr_gen(rng); }
};

///////////////////////////////////////////////////////////////////////
// ///////////////////// Normal_Distr : class ////////////////////// //
///////////////////////////////////////////////////////////////////////

class Normal_Distr : public Distribution_1D
{
private:
    std::normal_distribution<double> distr_gen;

public:
    // Constructors
    Normal_Distr(kwargs new_args = {{"mean", (double)0.0}, {"sd", (double)1.0}}) : Distribution_1D(new_args), distr_gen(get_arg("mean"), get_arg("sd")) {}
    Normal_Distr(kwargs_ext new_args_ext) : Distribution_1D(new_args_ext), distr_gen(get_arg("mean"), get_arg("sd")) {}
    ~Normal_Distr() {}

    double pdf(double x) override;
    double log_pdf(double x) override;
    double log_ratio(double x, double y) override;
    double rand(std::mt19937_64 &rng) override { return distr_gen(rng); }
};

///////////////////////////////////////////////////////////////////////
// ///////////////////// Lognorm_Distr : class ///////////////////// //
///////////////////////////////////////////////////////////////////////

class Lognorm_Distr : public Distribution_1D
{
private:
    std::lognormal_distribution<double> distr_gen;

public:
    // Constructors
    Lognorm_Distr(kwargs new_args = {{"mean", (double)0.0}, {"sd", (double)1.0}}) : Distribution_1D(new_args), distr_gen(get_arg("mean"), get_arg("sd")) {}
    Lognorm_Distr(kwargs_ext new_args_ext) : Distribution_1D(new_args_ext), distr_gen(get_arg("mean"), get_arg("sd")) {}
    ~Lognorm_Distr() {}

    double pdf(double x) override;
    double log_pdf(double x) override;
    double log_ratio(double x, double y) override;
    double rand(std::mt19937_64 &rng) override { return distr_gen(rng); }
};

///////////////////////////////////////////////////////////////////////
// ////////////////////// Gamma_Distr : class ////////////////////// //
///////////////////////////////////////////////////////////////////////

class Gamma_Distr : public Distribution_1D
{
private:
    std::gamma_distribution<double> distr_gen;

public:
    // Constructors
    Gamma_Distr(kwargs new_args = {{"shape", (double)0.0}, {"scale", (double)1.0}}) : Distribution_1D(new_args), distr_gen(get_arg("shape"), get_arg("scale")) {}
    Gamma_Distr(kwargs_ext new_args_ext) : Distribution_1D(new_args_ext), distr_gen(get_arg("shape"), get_arg("scale")) {}
    ~Gamma_Distr() {}

    double pdf(double x) override;
    double log_pdf(double x) override;
    double log_ratio(double x, double y) override;
    double rand(std::mt19937_64 &rng) override { return distr_gen(rng); }
};

///////////////////////////////////////////////////////////////////////
// /////////////////// Exponential_Distr : class /////////////////// //
///////////////////////////////////////////////////////////////////////

class Exponential_Distr : public Distribution_1D
{
private:
    std::exponential_distribution<double> distr_gen;

public:
    Exponential_Distr(kwargs new_args = {{"rate", (double)1.0}}) : Distribution_1D(new_args), distr_gen(get_arg("rate")) {}
    Exponential_Distr(kwargs_ext new_args_ext) : Distribution_1D(new_args_ext), distr_gen(get_arg("rate")) {}
    ~Exponential_Distr() {}

    double pdf(double x) override;
    double log_pdf(double x) override;
    double log_ratio(double x, double y) override;
    double rand(std::mt19937_64 &rng) override { return distr_gen(rng); }
};

///////////////////////////////////////////////////////////////////////
// //////////////// Onset to Transmission : class ////////////////// //
///////////////////////////////////////////////////////////////////////

class Onset_to_Transm : public Distribution
{
private:
    kwargs args = {}; // arguments
    kwargs_ext args_ext = {}; // arguments

public:
    // Constructors
    Onset_to_Transm(kwargs new_args) : Distribution(new_args) {}
    Onset_to_Transm(kwargs_ext new_args_ext) : Distribution(new_args_ext) {}
    ~Onset_to_Transm() {}

    double pdf(timetype tost, timetype t_incub);
    double cdf(timetype tost, timetype t_incub);
};

#endif // DEF_USER_MODEL__H