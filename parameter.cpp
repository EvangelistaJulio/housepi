// parameter.cpp //

#include "parameter.hpp"

using namespace std;

// UPDATE: too baroque maybe...

///////////////////////////////////////////////////////////
//////// Parameter : class definition /////////////////////
///////////////////////////////////////////////////////////

// ----------------------------
//  Constructors
// ----------------------------

Parameter::Parameter(
    double value,
    bool infer,
    std::string metropolis_strategy,
    kwargs metropolis_args,
    kwargs range,
    std::string prior_name,
    kwargs prior_args)
{
    p_value = value;
    p_infer = infer;

    // --- inference mode setup --- //
    if (infer == true)
    {
        metropolis = new Metropolis(metropolis_strategy, metropolis_args, range);

        // --- prior function --- //
        if (prior_name == "_none_")
        {
            p_prior_rand = nullptr;
            p_prior_log_pdf = nullptr;
            prior_args = {};
        }
        else if (prior_name == "unif")
        {
            p_prior_rand = unif_rand;
            p_prior_log_pdf = log_unif_pdf;
        }
        else if (prior_name == "normal")
        {
            p_prior_rand = normal_rand;
            p_prior_log_pdf = log_normal_pdf;
        }
        else if (prior_name == "lognorm")
        {
            p_prior_rand = lognorm_rand;
            p_prior_log_pdf = log_lognorm_pdf;
        }
        else if (prior_name == "gamma")
        {
            p_prior_rand = gamma_rand;
            p_prior_log_pdf = log_gamma_pdf;
        }
        else
        {
            std::cout << "Unknown prior " << prior_name << std::endl;
            exit(EXIT_FAILURE);
        }
        p_prior_args = prior_args;
    }
}

// ----------------------------- //
// --- Methods : Computation --- //
// ----------------------------- //

/* generate_from_prior */
double Parameter::generate_from_prior(std::mt19937_64& rand_gen)
{
    return p_prior_rand(rand_gen, p_prior_args);
}

/* log_ratio_prior */
double Parameter::log_ratio_prior(double old_value, double new_value)
{
    return p_prior_log_pdf(new_value, p_prior_args) - p_prior_log_pdf(old_value, p_prior_args);
}

///////////////////////////////////////////////////////////
// Parameter_Dictionary : Class 
///////////////////////////////////////////////////////////

// --- Getter --- //

/* get_value */
double Parameter_Dictionary::get_value(std::string parID) const
{
    return parameters.at(parID).get_value();
}

/* get_infer */
bool Parameter_Dictionary::get_infer(std::string parID) const
{
    return parameters.at(parID).p_infer;
}

/* get_accept_rate */
double Parameter_Dictionary::get_accept_rate(std::string parID) const
{
    return parameters.at(parID).get_accept_rate();
}

// --- Setter --- //

/* add_param */
void Parameter_Dictionary::add_param(std::string parID, Parameter new_param)
{
    parameters.emplace(parID, new_param);
    if (new_param.p_infer == true)
    {
        to_be_inferred.push_back(parID);
    }
}

/* set_value */
void Parameter_Dictionary::set_value(std::string parID, double new_value)
{
    parameters.find(parID)->second.set_value(new_value);
}

/* init_value */
void Parameter_Dictionary::init_value(std::string parID, std::mt19937_64 &m_rand_gen)
{
    double new_value = parameters.find(parID)->second.generate_from_prior(m_rand_gen);
    set_value(parID, new_value);
}

/* reset_moves */
void Parameter_Dictionary::reset_moves(std::string parID)
{
    parameters.find(parID)->second.reset_moves();
}

/* set_num_propos_pp */
void Parameter_Dictionary::set_num_propos_pp(std::string parID)
{
    parameters.find(parID)->second.set_num_propos_pp();
}

/* set_num_accept_pp */
void Parameter_Dictionary::set_num_accept_pp(std::string parID)
{
    parameters.find(parID)->second.set_num_accept_pp();
}

/* update_accept_rate */
void Parameter_Dictionary::update_accept_rate(std::string parID, bool move_accepted, int step)
{
    parameters.find(parID)->second.update_accept_rate(move_accepted, step);
}

/* adapt */
void Parameter_Dictionary::adapt(std::string parID, int step, int m_num_adapt_iter)
{
    parameters.find(parID)->second.adapt(step, m_num_adapt_iter);
}

// --- Computations --- //

/* within_range */
bool Parameter_Dictionary::valid_param(std::string parID, double new_value)
{
    return parameters.find(parID)->second.valid(new_value);
}

/* propose */
double Parameter_Dictionary::propose(std::string parID, std::mt19937_64 &m_rand_gen, double old_value)
{
    return parameters.find(parID)->second.propose(m_rand_gen, old_value);
}

/* log_ratio_proposal */
double Parameter_Dictionary::log_ratio_proposal(std::string parID, double old_value, double new_value)
{
    return parameters.find(parID)->second.log_ratio_proposal(old_value, new_value);
}

/* log_ratio_prior */
double Parameter_Dictionary::log_ratio_prior(std::string parID, double old_value, double new_value)
{
    return parameters.find(parID)->second.log_ratio_prior(old_value, new_value);
}

// --- Dispaly --- //

/* print */
void Parameter_Dictionary::print()
{
    std::cout << "_NAME_" << TAB << "_VALUE_" << TAB << "_INFER_" << endl;
    for (const auto &entry : parameters)
    {
        Parameter param = entry.second;
        std::cout << entry.first << TAB << param.get_value();
        std::cout << TAB << (param.p_infer > 0 ? "true" : "false") << endl;
    }
}
