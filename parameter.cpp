// parameter.cpp //

#include "parameter.hpp"

using namespace std;

// UPDATE: too baroque maybe...

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
    kwargs args,
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
        // FIXME : this should be internal to Strategy_Param
        if (metropolis_strategy == "exp")
        {
            _tmp_ = new Strategy_Param_exp(metropolis_rate, range);
        }
        else if (metropolis_strategy == "normal")
        {
            _tmp_ = new Strategy_Param_normal(metropolis_rate, range);
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

//  NOTE: aliases... but ugly...

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

double Par_Dict::propose(std::string parID, std::mt19937_64 &m_rand_gen)
{
    return parameters.find(parID)->second.propose(m_rand_gen);
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
    std::cout << "_NAME_" << TAB << "_VALUE_" << TAB << "_INFER_" << endl;
    for (const auto &entry : parameters)
    {
        Parameter param = entry.second;
        std::cout << param.p_name << TAB << param.get_value();
        std::cout << TAB << (param.p_infer > 0 ? "true" : "false")  << endl;
    }
}
