/*
 * varia.cpp
 *
 * Copyright (C) 2024 Francesco Camaglia, Institut Pasteur
 */

#include "varia.h"

//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//

///////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////  Base_Proposer : constructor  ////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

Base_Proposer::Base_Proposer(Base_Distribution *distr,
                             double rule_start,
                             std::vector<int> args_to_adapt) : distr_ptr(distr)
{
    _args_to_adapt_ = args_to_adapt;

    if (rule_start == INF_DOUBLE)
    {
        _is_start_from_zero_ = false;
    }
    else
    {
        _is_start_from_zero_ = true;
        set_rule_start_value(rule_start);
    }
};

///////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////  Base_Proposer : methods  //////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/* Increase spread if acceptation rate is high and viceversa. */
void Base_Proposer::adapt(double ratio)
{
    double factor = 1. + (1. - ratio) * (get_accept_rate() - TARGET_ACCEPT_RATE);

    for (const int &argID : _args_to_adapt_)
    {
        (distr_ptr->set_arg)(argID, ((distr_ptr->get_arg)(argID)) * factor);
    }
}

/* Reset acceptation rate. */
void Base_Proposer::reset_accept_rate()
{
    accept_rate = 0.0;
    num_proposed_moves = 0;
};

/* Update the parameter accept rate. */
void Base_Proposer::update_accept_rate(bool is_accepted)
{
    accept_rate = (double)(accept_rate * num_proposed_moves + is_accepted) / (num_proposed_moves + 1);
    num_proposed_moves += 1;
}

/* Switchboard to select the custom metropolis proposer rule. */
Base_Proposer *proposer_switchboard(std::string rule_name, dictargs rule_args, double rule_from)
{
    Base_Proposer *output = new Base_Proposer(new Base_Distribution());

    if (rule_name == "unif")
    {
        output = new Uniform_Proposer(rule_args, rule_from);
    }
    else if (rule_name == "norm")
    {
        output = new Normal_Proposer(rule_args, rule_from);
    }
    else if (rule_name == "lnorm")
    {
        output = new Lognorm_Proposer(rule_args, rule_from);
    }
    else if (!tools::is_na(rule_name))
    {
        throw std::invalid_argument("Unknown proposer rule '" + rule_name + "' for Metropolis-Hastings.");
    }

    return output;
}

//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//

///////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////  Parameter : methods  ///////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/* Init the values for the Bayesian diagnostics. */
void Parameter::init_diagnostics(ulong sample_size)
{
    dgn_count = 0;
    dgn_value.clear();
    dgn_value = std::vector<double>(sample_size, std::numeric_limits<double>::quiet_NaN());
}

/* Update the values for the Bayesian diagnostics. */
void Parameter::update_diagnostics()
{
    dgn_value[dgn_count++] = get_value();
}

/* Message to be used at the moment of display diagnostics. (TEMPORARY) */
std::string Parameter::show_diagnostics() const
{
    using namespace tools;
    std::ostringstream oss;

    // Welford online algorithm to compute statistics
    double avg = 0, var = 0;
    size_t N = dgn_value.size();
    double x, delta, delta2;
    double sqr_diff = 0.0;

    avg = 0.0;
    for (size_t i = 0; i < N; ++i)
    {
        x = dgn_value[i];
        delta = x - avg;
        avg += delta / (i + 1);
        delta2 = x - avg;
        sqr_diff += delta * delta2;
    }
    var = sqr_diff / (N - 1);

    oss << avg << C_DELIM << std::sqrt(var) << C_DELIM << truth(infer)
        << C_DELIM << show_distr_name() << C_DELIM << show_distr_arg_names() << C_DELIM << show_distr_arg_values()
        << C_DELIM << show_prop_name() << C_DELIM << show_prop_arg_names() << C_DELIM << show_prop_arg_values()
        << C_DELIM << show_range_min() << C_DELIM << show_range_max();
    return oss.str();
}
