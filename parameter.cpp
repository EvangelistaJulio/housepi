// parameter.cpp //

#include "parameter.hpp"

using namespace std;

//----------------------------------------------------------------------
// Strategy : functional methods
//----------------------------------------------------------------------

// --- prior --- //
// FIXME : missing

/* adapt : increase spread argument accept_rate is small and viceversa */
double Base_Proposer::_adapter_(int step, int num_adapt_iter)
{
    return 1. + (1. - step / num_adapt_iter) * (DEFAULT_ACCEPT_RATE - get_accept_rate());
}

/* update accept rate */
void Base_Proposer::update_accept_rate(bool is_accepted, int step)
{
    accept_rate = ((is_accepted + accept_rate * (step - 1)) / step);
}

/* reset moves */
void Base_Proposer::reset_moves()
{
    num_propos = 0;
    num_accept = 0;
    accept_rate = 0.0;
}

/* set range */
void Base_Proposer::set_range(kwargs new_range)
{
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
        _valid_func_ = &Base_Proposer::_valid_dummy_;
    }
    else if (condition == 1)
    {
        _valid_func_ = &Base_Proposer::_valid_range_low_;
    }
    else if (condition == 2)
    {
        _valid_func_ = &Base_Proposer::_valid_range_up_;
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
            _valid_func_ = &Base_Proposer::_valid_range_both_;
        }
    }
}

///////////////////////////////////////////////////////////
/////////////// Parameter : constructor ///////////////////
///////////////////////////////////////////////////////////

Parameter::Parameter(
    double *new_value_ptr,
    std::string prior_name,
    kwargs prior_args,
    std::string strat_name,
    kwargs strat_args,
    kwargs value_range) : _base_(nullptr, nullptr)
{
    p_value_ptr = new_value_ptr;

    // --- metropolis strategy --- //
    Base_Proposer *proposer;

    if (strat_name == UNDEF)
    {
        p_infer = false;
        proposer = nullptr;
    }
    else if (strat_name == "normal")
    {
        p_infer = true;
        proposer = new Normal_Proposer(strat_args, value_range);
    }
    else if (strat_name == "lognorm")
    {
        p_infer = true;
        proposer = new Lognorm_Proposer(strat_args, value_range);
    }
    else
    {
        std::cerr << "ERROR: Unknown metropolis strategy " << strat_name << std::endl;
        std::cerr << "Exit...\n"
                  << std::endl;
        exit(EXIT_FAILURE);
    }

    // --- prior function --- //
    Distribution_1D *distribution;

    if (p_infer == true)
    {
        if (prior_name == UNDEF)
        {
            distribution = nullptr;
        }
        else if (prior_name == "unif")
        {
            distribution = new Uniform_Distr(prior_args);
        }
        else if (prior_name == "normal")
        {
            distribution = new Normal_Distr(prior_args);
        }
        else if (prior_name == "lognorm")
        {
            distribution = new Lognorm_Distr(prior_args);
        }
        else if (prior_name == "gamma")
        {
            distribution = new Gamma_Distr(prior_args);
        }
        else if (prior_name == "exp")
        {
            distribution = new Exponential_Distr(prior_args);
        }
        else
        {
            std::cerr << "\nERROR: Unknown distribution " << prior_name << std::endl;
            std::cerr << "Exit...\n"
                      << std::endl;
            exit(EXIT_FAILURE);
        }

        _base_ = Base_Parameter(distribution, proposer);
    }
}

///////////////////////////////////////////////////////////
/////////////// Augmented : constructor ///////////////////
///////////////////////////////////////////////////////////

// FIXME : can I compact this through templates of a common wrapper of Base_Parameter

Augmented::Augmented(
    double new_value,
    std::string model_name,
    kwargs_ext model_args_ext,
    std::string strat_name,
    kwargs strat_args,
    kwargs value_range) : _base_(nullptr, nullptr)
{
    a_value = new_value;
    a_value_0 = new_value;

    // --- metropolis strategy --- //
    Base_Proposer *proposer;

    if (strat_name == UNDEF)
    {
        a_infer = false;
        proposer = nullptr;
    }
    else if (strat_name == "normal")
    {
        a_infer = true;
        proposer = new Normal_Proposer(strat_args, value_range);
    }
    else if (strat_name == "lognorm")
    {
        a_infer = true;
        proposer = new Lognorm_Proposer(strat_args, value_range);
    }
    else
    {
        std::cerr << "ERROR: Unknown metropolis strategy " << strat_name << std::endl;
        std::cerr << "Exit...\n"
                  << std::endl;
        exit(EXIT_FAILURE);
    }

    // --- prior function --- //
    Distribution_1D *distribution;

    if (a_infer == true)
    {
        if (model_name == UNDEF)
        {
            distribution = nullptr;
        }
        else if (model_name == "unif")
        {
            distribution = new Uniform_Distr(model_args_ext);
        }
        else if (model_name == "normal")
        {
            distribution = new Normal_Distr(model_args_ext);
        }
        else if (model_name == "lognorm")
        {
            distribution = new Lognorm_Distr(model_args_ext);
        }
        else if (model_name == "gamma")
        {
            distribution = new Gamma_Distr(model_args_ext);
        }
        else if (model_name == "exp")
        {
            distribution = new Exponential_Distr(model_args_ext);
        }
        else
        {
            std::cerr << "\nERROR: Unknown distribution " << model_name << std::endl;
            std::cerr << "Exit...\n"
                      << std::endl;
            exit(EXIT_FAILURE);
        }

        _base_ = Base_Parameter(distribution, proposer);
    }
}

///////////////////////////////////////////////////////////
//////// Parameter_Dictionary : methods ///////////////////
///////////////////////////////////////////////////////////

/* add_param */
void Parameter_Dictionary::add_param(
    const std::string k,
    double new_value,
    std::string prior_name,
    kwargs prior_args,
    std::string strat_name,
    kwargs strat_args,
    kwargs value_range)
{
    param_values[_counter_] = new_value;
    Parameter new_param(&param_values[_counter_], prior_name, prior_args, strat_name, strat_args, value_range);

    if (params.find(k) == params.end())
    {
        params.emplace(k, new_param);
    }
    else
    {
        std::cerr << "ERROR: Parameter `" << k << "` has multiple definitions." << std::endl;
        std::cerr << "Exit...\n"
                  << std::endl;
    }

    if (new_param.get_infer() == true)
    {
        to_be_inferred.push_back(k);
    }
    else
    {
        not_to_be_inferred.push_back(k);
    }
    _counter_++;
}

/* print */
void Parameter_Dictionary::print()
{
    std::cout << "_NAME_" << TAB << "_VALUE_" << TAB << "_INFER_" << endl;
    for (std::string k : to_be_inferred)
    {
        std::cout << k << TAB << get_value(k);
        std::cout << TAB << "true" << endl;
    }
    for (std::string k : not_to_be_inferred)
    {
        std::cout << k << TAB << get_value(k);
        std::cout << TAB << "false" << endl;
    }
}