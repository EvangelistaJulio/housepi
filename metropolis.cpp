// metropolis.cpp //

#include "metropolis.hpp"

using namespace std;

/////////////////////////////////////////////////////
//  Metropolis : class definition
/////////////////////////////////////////////////////

Metropolis::Metropolis(
    double new_value,
    std::string strat_name,
    kwargs strat_args,
    kwargs value_range,
    std::string prior_name,
    kwargs prior_args)
{
    value = new_value;
    value_0 = new_value;

    // FIXME : a consistency valid command could be useful at the beginning
    // FIXME: implement consistency

    // --- metropolis strategy --- //
    if (strat_name == "_none_")
    {
        infer = false;
        _propose_func_ = nullptr;
        _propose_log_ratio_func_ = nullptr;
        _adapt_func_ = nullptr;
    }
    else if (strat_name == "prior")
    {
        infer = true;
        value_0 = 0.0;
        _propose_func_ = &Metropolis::_prior_propose_;                                 
        _propose_log_ratio_func_ = &Metropolis::_normal_propose_log_ratio_; // FIXME: 0.! also log_ratio_prior
        _adapt_func_ = &Metropolis::_normal_adapt_;                          // FIXME : skip
    }
    else if (strat_name == "normal")
    {
        infer = true;
        _propose_func_ = &Metropolis::_normal_propose_;
        _propose_log_ratio_func_ = &Metropolis::_normal_propose_log_ratio_;
        _adapt_func_ = &Metropolis::_normal_adapt_;
    }
    else if (strat_name == "lognorm")
    {
        infer = true;
        _propose_func_ = &Metropolis::_lognorm_propose_;
        _propose_log_ratio_func_ = &Metropolis::_lognorm_propose_log_ratio_;
        _adapt_func_ = &Metropolis::_lognorm_adapt_;
    }
    else if (strat_name == "exp")
    {
        infer = true;
        _propose_func_ = &Metropolis::_exp_propose_;
        _propose_log_ratio_func_ = &Metropolis::_exp_propose_log_ratio_;
        _adapt_func_ = &Metropolis::_exp_adapt_;
    }
    else if (strat_name == "gamma")
    {
        infer = true;
        _propose_func_ = &Metropolis::_gamma_propose_;
        _propose_log_ratio_func_ = &Metropolis::_gamma_propose_log_ratio_;
        _adapt_func_ = &Metropolis::_gamma_adapt_;
    }
    else
    {
        std::cerr << "ERROR: Unknown metropolis strategy " << strat_name << std::endl;
        std::cerr << "Exit...\n"
                  << std::endl;
        exit(EXIT_FAILURE);
    }

    set_propose_args(strat_args);
    set_range(value_range);
    reset_moves();

    if (infer == false)
    {
        prior_name = "_none_";
    }
    distribution_switchboard(prior_name, &_prior_pdf_, &_prior_log_pdf_, &_prior_rand_);

}

void Metropolis::set_propose_args(kwargs new_args)
{
    propose_args = new_args;
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

double Metropolis::_prior_propose_log_ratio_(double old_value, double new_value)
{
    return _prior_log_pdf_(new_value, get_prior_args()) - _prior_log_pdf_(old_value, get_prior_args());
}

double Metropolis::_prior_propose_(std::mt19937_64 &rand_gen, double _)
{
    return _prior_rand_(rand_gen, get_prior_args());
}

// --- normal --- //

double Metropolis::_normal_propose_(std::mt19937_64 &rand_gen, double old_value)
{
    return old_value + normal_rand(rand_gen, propose_args);
}

double Metropolis::_normal_propose_log_ratio_(double old_value, double new_value)
{
    double epsilon = new_value - old_value;

    return log_normal_pdf(-epsilon, propose_args) - log_normal_pdf(epsilon, propose_args);
}

void Metropolis::_normal_adapt_(int step, int num_adapt_iter)
{
    _adapter_("sd", step, num_adapt_iter);
}

// --- lognorm --- //

double Metropolis::_lognorm_propose_(std::mt19937_64 &rand_gen, double old_value)
{
    return old_value * lognorm_rand(rand_gen, propose_args);
}

double Metropolis::_lognorm_propose_log_ratio_(double old_value, double new_value)
{
    double epsilon = new_value / old_value;
    double i_epsilon = old_value / new_value;

    return log_lognorm_pdf(i_epsilon, propose_args) - log_lognorm_pdf(epsilon, propose_args);
}

void Metropolis::_lognorm_adapt_(int step, int num_adapt_iter)
{
    _adapter_("sd", step, num_adapt_iter);
}

// --- gamma --- //

double Metropolis::_gamma_propose_(std::mt19937_64 &rand_gen, double old_value)
{
    return gamma_rand(rand_gen, propose_args); // if drawing from prior, now old_value!!!
}

double Metropolis::_gamma_propose_log_ratio_(double old_value, double new_value)
{
    double epsilon = new_value - old_value;

    return 0.; // if drawing from prior, now old_value!!!
}

void Metropolis::_gamma_adapt_(int step, int num_adapt_iter)
{
    _adapter_("scale", step, num_adapt_iter);
}

// --- exp --- //

double Metropolis::_exp_propose_(std::mt19937_64 &rand_gen, double old_value)
{
    return old_value * exp_rand(rand_gen, propose_args);
}

double Metropolis::_exp_propose_log_ratio_(double old_value, double new_value)
{
    double epsilon = new_value / old_value;
    double i_epsilon = old_value / new_value;

    return log_exp_pdf(i_epsilon, propose_args) - log_exp_pdf(epsilon, propose_args);
}

void Metropolis::_exp_adapt_(int step, int num_adapt_iter)
{
    _adapter_("rate", step, num_adapt_iter);
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