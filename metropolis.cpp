// metropolis.cpp //

#include "metropolis.hpp"

using namespace std;

/////////////////////////////////////////////////////
//  Metropolis : class definition
/////////////////////////////////////////////////////

Metropolis::Metropolis(std::string strategy, kwargs init_args, kwargs init_range)
{
    // --- metropolis strategy --- //
    if (strategy == "_none_")
    {
        _propose_func_ = nullptr;
        _log_ratio_propose_func_ = nullptr;
        _adapt_func_ = nullptr;
    }
    else if (strategy == "normal")
    {
        _propose_func_ = &Metropolis::_normal_propose_;
        _log_ratio_propose_func_ = &Metropolis::_normal_log_ratio_proposal_;
        _adapt_func_ = &Metropolis::_normal_adapt_;
    }
    else if (strategy == "exp")
    {
        _propose_func_ = &Metropolis::_exp_propose_;
        _log_ratio_propose_func_ = &Metropolis::_exp_log_ratio_proposal_;
        _adapt_func_ = &Metropolis::_exp_adapt_;
    }
    else if (strategy == "gamma")
    {
        _propose_func_ = &Metropolis::_gamma_propose_;
        _log_ratio_propose_func_ = &Metropolis::_gamma_log_ratio_proposal_;
        _adapt_func_ = &Metropolis::_gamma_adapt_;
    }
    else
    {
        std::cout << "ERROR: Unknown metropolis strategy " << strategy << std::endl;
        std::cerr << "Exit...\n"
                  << std::endl;
        exit(EXIT_FAILURE);
    }

    set_args(init_args);
    set_range(init_range);
    reset_moves();
    // FIXME : a consistency valid command could be useful at the beginning
}

void Metropolis::set_args(kwargs new_args)
{
    args = new_args;
    // FIXME : a consistency valid command could be useful at the beginning
    // FIXME : how to point to global variables?
}

void Metropolis::set_range(kwargs new_range)
{
    // FIXME (FUTURE) : functional programming to go beyond map and allow for single valid in within_range
    int condition = 0;

    if (new_range.find("low") != new_range.end())
    {
        _range_low_ = new_range.at("low");
        condition += 1;
    }

    if (new_range.find("up") != new_range.end())
    {
        _range_up_ = new_range.at("up");
        condition += 2;
    }

    if (condition == 0)
    {
        _valid_func_ = &Metropolis::_valid_dummy_;
    }
    else if (condition == 1)
    {
        _valid_func_ = &Metropolis::_valid_range_low_;
    }
    else if (condition == 2)
    {
        _valid_func_ = &Metropolis::_valid_range_up_;
    }
    else
    {
        if (_range_low_ >= _range_up_)
        {
            std::cerr << "ERROR: Provided an inconsistent range: {" << new_range.at("low") << ", " << new_range.at("up") << std::endl;
            std::cerr << "Exit...\n"
                      << std::endl;
            exit(EXIT_FAILURE);
        }
        else
        {
            _valid_func_ = &Metropolis::_valid_range_both_;
        }
    }
}

void Metropolis::reset_moves()
{
    num_propos = 0;
    num_accept = 0;
    accept_rate = 0.0;
}

void Metropolis::update_accept_rate(bool move_accepted, int step)
{
    // update accept rate
    accept_rate = ((move_accepted + accept_rate * (step - 1)) / step);
}

//----------------------------------------------------------------------
// Strategy : functional methods
//----------------------------------------------------------------------

// --- normal --- //

double Metropolis::_normal_propose_(std::mt19937_64& rand_gen, double old_value)
{
    return old_value + normal_rand(rand_gen, args);
}

double Metropolis::_normal_log_ratio_proposal_(double old_value, double new_value)
{
    return log_normal_pdf(new_value, args) - log_normal_pdf(old_value, args);
}

void Metropolis::_normal_adapt_(int step, int num_adapt_iter)
{
    double delta, diff;

    delta = 1.0 - step / num_adapt_iter;
    diff = DEFAULT_ACCEPT_RATE - accept_rate;
    set_args(
        {{"mean", args.at("mean")},
         {"sd", (1.0 + delta * diff) * args.at("sd")}});
}

// --- exp --- //

double Metropolis::_exp_propose_(std::mt19937_64& rand_gen, double old_value)
{
    return old_value * exp_rand(rand_gen, args);
}

double Metropolis::_exp_log_ratio_proposal_(double old_value, double new_value)
{
    return log_exp_pdf(new_value, args) - log_exp_pdf(old_value, args);
}

void Metropolis::_exp_adapt_(int step, int num_adapt_iter)
{
    double delta, diff;

    delta = 1.0 - step / num_adapt_iter;
    diff = DEFAULT_ACCEPT_RATE - accept_rate;
    set_args(
        {{"rate", (1.0 + delta * diff) * args.at("rate")}});
}

// --- gamma --- //

double Metropolis::_gamma_propose_(std::mt19937_64& rand_gen, double old_value)
{
    return old_value + gamma_rand(rand_gen, args);
}

double Metropolis::_gamma_log_ratio_proposal_(double old_value, double new_value)
{
    return log_gamma_pdf(new_value, args) - log_gamma_pdf(old_value, args);
}

void Metropolis::_gamma_adapt_(int step, int num_adapt_iter)
{
    double delta, diff;

    delta = 1.0 - step / num_adapt_iter;
    diff = DEFAULT_ACCEPT_RATE - accept_rate;
    set_args(
        {{"shape", args.at("shape")},
         {"scale", (1.0 + delta * diff) * args.at("scale")}});
}

//----------------------------------------------------------------------
// Validation : functional methods
//----------------------------------------------------------------------

bool Metropolis::_valid_dummy_(double new_value)
{
    return (true);
}

bool Metropolis::_valid_range_low_(double new_value)
{
    return (new_value > _range_low_);
}

bool Metropolis::_valid_range_up_(double new_value)
{
    return (new_value < _range_up_);
}

bool Metropolis::_valid_range_both_(double new_value)
{
    return _valid_range_low_(new_value) || _valid_range_up_(new_value);
}