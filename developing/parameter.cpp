// parameter.cpp //

#include "parameter.hpp"

using namespace std;

//----------------------------------------------------------------------
// Strategy : functional methods
//----------------------------------------------------------------------

// --- prior --- //
// FIXME : missing

/* adapt : increase spread argument accept_rate is small and viceversa */
double Metropolis_Hastings::_adapter_(int step, int num_adapt_iter)
{
    return 1. + (1. - step / num_adapt_iter) * (DEFAULT_ACCEPT_RATE - get_accept_rate());
}

/* update accept rate */
void Metropolis_Hastings::update_accept_rate(bool is_accepted, int step)
{
    accept_rate = ((is_accepted + accept_rate * (step - 1)) / step);
}

/* reset moves */
void Metropolis_Hastings::reset_moves()
{
    num_propos = 0;
    num_accept = 0;
    accept_rate = 0.0;
}

/* set range */
void Metropolis_Hastings::set_range(kwargs new_range)
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
        _valid_func_ = &Metropolis_Hastings::_valid_dummy_;
    }
    else if (condition == 1)
    {
        _valid_func_ = &Metropolis_Hastings::_valid_range_low_;
    }
    else if (condition == 2)
    {
        _valid_func_ = &Metropolis_Hastings::_valid_range_up_;
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
            _valid_func_ = &Metropolis_Hastings::_valid_range_both_;
        }
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