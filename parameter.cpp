#include "parameter.hpp"

using namespace std;

/////////////////////////////////////////////////////
//  Metropolis_Strategy : Class
/////////////////////////////////////////////////////

Metropolis_Strategy::Metropolis_Strategy(double new_rate, std::vector<double> new_range)
{
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

    // check range values
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

    reset_moves();
}

void Metropolis_Strategy::reset_moves()
{
    num_propos = 0;
    num_accept = 0;
    accept = 0.0;
}

bool Metropolis_Strategy::within_range(double new_value)
{
    return (new_value > range[0]) && (new_value < range[1]);
}

void Metropolis_Strategy::adapt_rate(bool move_accepted, int step, int num_adapt_iter)
{
    double delta, diff;

    accept = ((move_accepted + (step - 1) * accept) / step);
    
    if (num_adapt_iter > 0 && step <= num_adapt_iter)
    {
        delta = 1.0 - step / num_adapt_iter;
        diff = accept - MCMC_OPTIMAL_ACCEPTANCE;
        set_rate(1.0 + delta * diff);
    }
}

// --- Strategy_exp ---
double Strategy_exp::propose_move(std::mt19937_64 &rand_gen, double old_value)
{
    double new_value = 0.0;
    new_value = old_value * std::exp(rate * normal_rand(rand_gen));

    return new_value;
}

double Strategy_exp::log_ratio_proposal()
{
    return rate;
}

// --- Strategy_normal ---
double Strategy_normal::propose_move(std::mt19937_64 &rand_gen, double old_value)
{
    double new_value = 0.0;
    std::vector<double> args = {0.0, rate};
    new_value = old_value + normal_rand(rand_gen, args);

    return new_value;
}

double Strategy_normal::log_ratio_proposal()
{
    return 0.;
}

/////////////////////////////////////////////////////
//  Parameter : Class
/////////////////////////////////////////////////////

// ----------------------------
//  Constructors
// ----------------------------

Parameter::Parameter(
    std::string name,
    double value,
    bool infer,
    std::string prior,
    std::vector<double> args,
    std::string metropolis_strategy,
    double metropolis_rate,
    std::vector<double> range)
{
    p_name = name;
    p_value = value;

    // --- inference mode setup --- //
    p_infer = infer;

    if (infer == true)
    {
        // --- prior function --- //
        if (prior == "unif")
        {
            p_prior_rand = unif_rand;
            p_prior_log_pdf = log_unif_pdf;
        }
        else if (prior == "normal")
        {
            p_prior_rand = normal_rand;
            p_prior_log_pdf = log_normal_pdf;
        }
        else if (prior == "lognorm")
        {
            p_prior_rand = lognorm_rand;
            p_prior_log_pdf = log_lognorm_pdf;
        }
        else if (prior == "gamma")
        {
            p_prior_rand = gamma_rand;
            p_prior_log_pdf = log_gamma_pdf;
        }
        else
        {
            std::cout << "Unknown prior " << prior << std::endl;
            exit(EXIT_FAILURE);
        }

        p_prior = prior;
        p_prior_args = args;

        // --- metropolis strategy --- //
        if (metropolis_strategy == "exp")
        {
            _metropolis = new Strategy_exp(metropolis_rate, range);
        }
        else if (metropolis_strategy == "normal")
        {
            _metropolis = new Strategy_normal(metropolis_rate, range);
        }
        else
        {
            std::cout << "ERROR: Unknown metropolis strategy " << metropolis_strategy << std::endl;
            std::cerr << "Exit...\n"
                      << std::endl;
            exit(EXIT_FAILURE);
        }

        reset_moves();
    }
}

// ----------------------------
//  Methods : Getter
// ----------------------------

double Parameter::get_value() const
{
    return p_value;
}

// ----------------------------
//  Methods : Setter
// ----------------------------

void Parameter::set_value(double value)
{
    p_value = value;
}

void Parameter::set_rate_metropolis(double new_rate)
{
    _metropolis->set_rate(new_rate);
}

// ----------------------------
//  Methods : Computation
// ----------------------------

double Parameter::generate_from_prior(std::mt19937_64 &rand_gen)
{
    return p_prior_rand(rand_gen, p_prior_args);
}

double Parameter::log_ratio_prior(double new_value)
{
    return p_prior_log_pdf(new_value, p_prior_args) - p_prior_log_pdf(p_value, p_prior_args);
}

///////////////////////////////////////////////////////////
//////// Par_Dict : Class
///////////////////////////////////////////////////////////

// Getter
double Par_Dict::value(std::string parID) const
{
    return parameters.at(parID).get_value();
}

bool Par_Dict::infer(std::string parID) const
{
    return parameters.at(parID).p_infer;
}

int Par_Dict::num_propos(std::string parID) const
{
    return parameters.at(parID).get_num_propos();
}

int Par_Dict::num_accept(std::string parID) const
{
    return parameters.at(parID).get_num_accept();
}

// Setter
void Par_Dict::add_param(std::string parID, Parameter new_param)
{
    parameters.emplace(parID, new_param);
    if (new_param.p_infer == true)
    {
        to_be_inferred.push_back(parID);
    }
}

void Par_Dict::set_value(std::string parID, double new_value)
{
    parameters.find(parID)->second.set_value(new_value);
}

void Par_Dict::init_value(std::string parID, std::mt19937_64 &m_rand_gen)
{
    double new_value = parameters.find(parID)->second.generate_from_prior(m_rand_gen);
    set_value(parID, new_value);
}

void Par_Dict::reset_moves(std::string parID)
{
    parameters.find(parID)->second.reset_moves();
}

void Par_Dict::set_num_propos_pp(std::string parID)
{
    parameters.find(parID)->second.set_num_propos_pp();
}

void Par_Dict::set_num_accept_pp(std::string parID)
{
    parameters.find(parID)->second.set_num_accept_pp();
}

void Par_Dict::adapt_rate(std::string parID, bool move_accepted, int step, int m_num_adapt_iter)
{
    parameters.find(parID)->second.adapt_rate(move_accepted, step, m_num_adapt_iter);
}

// Computations
bool Par_Dict::within_range(std::string parID, double new_value)
{
    return parameters.find(parID)->second.within_range(new_value);
}

double Par_Dict::propose_move(std::string parID, std::mt19937_64 &m_rand_gen)
{
    return parameters.find(parID)->second.propose_move(m_rand_gen);
}

double Par_Dict::log_ratio_proposal(std::string parID)
{
    return parameters.find(parID)->second.log_ratio_proposal();
}

double Par_Dict::log_ratio_prior(std::string parID, double new_value)
{
    return parameters.find(parID)->second.log_ratio_prior(new_value);
}

// Other
void Par_Dict::print()
{
    for (const auto &entry : parameters)
    {
        Parameter param = entry.second;
        std::cout << param.p_name << "\t" << param.get_value() << "\t" << param.p_infer << endl;
    }
}
