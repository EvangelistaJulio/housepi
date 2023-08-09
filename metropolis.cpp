// metropolis.cpp //

#include "metropolis.hpp"

using namespace std;

// UPDATE: too baroque maybe...

/////////////////////////////////////////////////////
//  Metropolis : Class
/////////////////////////////////////////////////////

Metropolis::Metropolis(std::map<std::string, double> init_args, std::vector<double> init_range)
{
    set_args(init_args);
    set_range(init_range);
    reset_moves();
}

void Metropolis::set_args(std::map<std::string, double> new_args)
{
    // UPDATE : a consistency check command could be useful
    /*
    // check metropolis_rate values
    if (new_rate <= 0. || new_rate >= 1.)
    {
        std::cout << "ERROR: Provided a rate not in [0,1] : " << new_rate << std::endl;
        std::cerr << "Exit...\n"
                  << std::endl;
        exit(EXIT_FAILURE);
    }
    else
    {
        rate = new_rate;
    }
    */
    args = new_args;
}

void Metropolis::set_range(std::vector<double> new_range)
{
    if (new_range.size() != 2)
    {
        std::cout << "ERROR: Provided a range of size: " << new_range.size() << std::endl;
        std::cerr << "Exit...\n"
                  << std::endl;
        exit(EXIT_FAILURE);
    }
    else if (new_range[0] >= new_range[1])
    {
        std::cout << "ERROR: Provided an inconsistent range: {" << new_range[0] << ", " << new_range[1] << std::endl;
        std::cerr << "Exit...\n"
                  << std::endl;
        exit(EXIT_FAILURE);
    }
    else
    {
        range = new_range;
    }
}

void Metropolis::reset_moves()
{
    num_propos = 0;
    num_accept = 0;
    accept = 0.0;
}

bool Metropolis::within_range(double new_value)
{
    return (new_value > range[0]) && (new_value < range[1]);
}

//----------------------------------------------------------------------
// Metropolis : Parameter Strategy auxiliary classes
//----------------------------------------------------------------------

void Strategy_Param::adapt_rate(bool move_accepted, int step, int num_adapt_iter)
{
    double delta, diff;

    accept = ((move_accepted + (step - 1) * accept) / step);

    if (num_adapt_iter > 0 && step <= num_adapt_iter)
    {
        delta = 1.0 - step / num_adapt_iter;
        diff = accept - DEFAULT_ACCEPT_RATE;
        set_args({{"rate", 1.0 + delta * diff}});
    }
}

// --- Strategy_Param_exp --- //

double Strategy_Param_exp::propose(std::mt19937_64 &rand_gen, double old_value)
{
    return old_value * std::exp(args["rate"] * normal_rand(rand_gen));
}

double Strategy_Param_exp::log_ratio_proposal() { return args["rate"]; }

// --- Strategy_Param_normal --- //

double Strategy_Param_normal::propose(std::mt19937_64 &rand_gen, double old_value)
{
    return old_value + normal_rand(rand_gen, {{"mean", 0.0}, {"sd", args["rate"]}});
}

double Strategy_Param_normal::log_ratio_proposal() { return 0.; }

//----------------------------------------------------------------------
// Metropolis : Augmented_Variable Strategy auxiliary classes
//----------------------------------------------------------------------
